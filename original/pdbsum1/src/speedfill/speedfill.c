/***********************************************************************

  speedfill.c - Program to gap-fill given chain or pair of chains using
                new fast cavity-finding algorithm

************************************************************************

Date:-         2 Aug 2000
Dir update:-   2 Dec 2012

Written by:-   Roman Laskowski
               Department of Crystallography       
               Birkbeck College
               University of London

               Program uses kd-tree routines written by Jonathan Barker.

@@@ Write out DNA chains if within range

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      File name         Description
-----------      ---------         -----------
CATHPARAM            -             Input parameter file giving all the
                                   directory paths for the files used by
                                   the program
<filename>       pdb_name          Input PDB file holding coords of
                                   protein to be processed
rasmol.dfn       dfn_name          Input file giving SECSURF's
                                   interpretation of the contants of the
                                   PDB file (proteins, DNA, ligands, etc)
speedfill.pdb    out_name          Output file of gap-sphere coords and
                                   radii
speedfill.arr        -             Output file of saved grid array giving
                                    masked points and points occupied by
                                   each cluster
speedfill.depth  depth_name        Output file listing all atoms and
                                   showing which cleft, if any, they are
                                   close to and their "depth" within that
                                   cleft. Optionally also contains all
                                   the surface vertices
speedfill.rasmol pdb_out_name      Output RasMol script file of the
                                   structure and selected gap-maps
speedfill.html   html_name         Output HTML file showing analysis of
                                   clefts
surfnet.fil      surfnet_name      Output file giving header records for
                                   running SURFNET on the speedfill.pdb
                                   file
speedfill.r3d    r3d_name          Output Raster3D file of the protein and
                                   clefts
chains.dat       chains_dat        Input file giving protein chain 
                                   equivalences
conserv.dat      conserv_dat       Input file giving residue conservation
                                   scores
XXXXc.gradesPE   conserv_dat       Input file giving ConSurf residue
                                   conservation scores

Compilation
-----------

cc -c speedfill.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -c kdtree.c
cc -c ras3d.c
cc -o speedfill speedfill.o common.o reapar.o readpdb.o kdtree.o ras3d.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c speedfill.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -O2 -funroll-loops -c kdtree.c
cc -O2 -funroll-loops -c ras3d.c
cc -o speedfill speedfill.o common.o reapar.o readpdb.o kdtree.o ras3d.o -lm -lz

Notes
-----

Speedfill for docking tests:

     speedfill -f file.new -centroid -col 3 -bv -profunc -log

To generate HTML page, run as follows:

     speedfill pdb_code -c X -html -ntop 10 -surfnet


Profiling:-

speedfill m002 -c H -c I -c A -c B -col 3 -html -ntop 10 -d -save

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> interpret_code
         -> open_output_file (in common.c)
   -> prepend_dir
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
   -> form_file_names
         -> check_file (in common.c)
   -> get_dfn_data
         -> string_truncate (in common.c)
         -> create_rasdef_entry
         -> get_ligand_residues
               -> get_token
               -> format_res_number
               -> create_rasdef_entry
   -> read_pdb
         -> string_truncate (in common.c)
         -> r_get_atom_details (in readpdb.c)
         -> get_element (in readpdb.c)
               -> is_hydrogen_atom (in readpdb.c)
                      -> check_for_metal (in readpdb.c)
         -> check_for_duplicate
         -> assign_to_object
         -> create_residue_record
         -> get_atom_radius
         -> create_atom_record
               -> get_atom_type
   -> calc_atom_limits
         -> get_max_min_coords
   -> mark_inrange_ligands
   -> write_raster3d_headers
   -> calculate_rascale
   -> write_raster3d_protein
         -> get_atom_colour
         -> get_colour_rgb
   -> read_chain_equivs
         -> create_chain_record
   -> get_conservation_scores
         -> check_chain_match
   -> get_consurf_scores
         -> get_tokens (in common.c)
         -> remove_leading_blanks (in common.c)
         -> format_res_num
         -> check_chain_match
   -> calc_grid_size
   -> read_delphi_file
   -> create_grid_array
   -> create_atom_pointers_array
   -> create_index_arrays
   -> sort_index_arrays
         -> shell_sort
   -> mask_grid_points
   -> mark_potential_gap_points
         -> get_atom_neighbour_list
               -> get_atom_start_end
                     -> binary_atom_search
               -> transfer_atom_list
               -> filter_atom_coord_range
         -> mark_sphere_locations
               -> get_nearest_grid_point
   -> place_spheres
         -> calculate_ijk_index
               -> shell_sort
         -> grow_sphere
         -> check_sphere
   -> generate_gap_spheres
         -> create_sphere_record
   -> cluster_spheresx
         -> join_clusters
   -> transfer_spheresx
   -> create_sphere_kdtree
   -> sort_clusters
         -> create_arrays
         -> count_spheres
   -> write_spheres
   -> open_rasmol_file
   -> write_script_headers
   -> open_output_file (in common.c)
   -> write_pymol_headers
   -> write_surfnet_headers
   -> write_out_coords
         -> write_script_coords
         -> write_pymol_script
         -> write_surfnet_line
   -> read_vertices
         -> string_truncate
         -> r_get_atom_details (in readpdb.c)
         -> create_vertex_record
   -> write_all_vertices
       -> write_vertices
       -> write_script_clusters
               -> write_restype_colours
                     -> get_rasmol_colour
                          -> get_colour_rgb
               -> write_restype_pymol
               -> write_conserv_colours
                     -> get_rasmol_colour
                          -> get_colour_rgb
               -> write_conserv_pymol
   -> calc_exp_table
   -> contour_clusters
         -> sort_by_volume
         -> get_protein_name
         -> write_html_headers
               -> profunc_html_header
         -> reset_grid_array
               -> get_grid_value
                     -> extract_grid_value
               -> get_grid_value (in ras3d.c)
         -> update_map
               -> get_nearest_grid_point
               -> get_grid_value (as above)
               -> encode_grid_value
         -> calculate_gap_volume
               -> get_grid_value (as above)
         -> get_contour_points
               -> get_grid_value (as above)
               -> create_vertex_record
               -> place_vertex_in_sort_list
               -> splice_bands_together
         -> calc_vertex_energies
               -> calculate_energy
         -> create_vertex_array
         -> calc_vertex_burial
         -> save_cluster_points
         -> generate_polygons
               -> binary_edge_search
               -> get_grid_value (as above)
               -> make_line
               -> make_polygons
                     -> generate_triangles
                           -> create_polygon_record
                                 -> create_polist_record
         -> order_polygon_vertices
               -> orien3
                     -> vecprd
                     -> dot3
               -> calculate_normal
         -> calculate_vertex_normals
         -> write_raster3d_polygons
               -> get_colour_cpk
                     -> get_colour_rgb
               -> get_colour_restype
                     -> get_colour_rgb
               -> get_colour_conserv
                     -> get_colour_rgb
               -> get_chain_colour
               -> get_colour_rgb
         -> filterx_vertices
               -> get_neighbour_list
                     -> get_start_end
                           -> binary_search
                     -> transfer_vertex_list
                     -> filter_coord_range
         -> filter_vertices    ... (new method)
               -> create_vertex_kdtree
                     -> KdTree_create
                           -> KdTreeNode_create
               -> Annulus_create
               -> KdTree_query
               -> KdTreeQuery_next
               -> Annulus_free
               -> KdTreeQuery_free
               -> KdTree_free
                     -> KdTreeNode_free
         -> calc_vertex_burial
               -> get_nearest_grid_point
               -> get_grid_value (as above)
         -> get_ligands_in_gap
               -> get_nearest_grid_point
               -> get_grid_value (as above)
         -> calculate_gap_depth
               -> mark_exposed_vertices
                     -> get_neighbour_list
                           -> get_start_end
                                 -> binary_search
                           -> transfer_atom_list
                           -> filter_coord_range
               -> get_vertex_joins
                     -> create_vjoin_record
               -> get_depth
         -> principal_components
               -> calc_eigen_values
               -> eigen_sort
         -> sort_vertices
               -> perform_vertex_sort
         -> write_vertices
         -> analyse_vertices
               -> get_restype
         -> write_cleft_atoms
         -> relink_vertices
         -> write_script_clusters
               -> get_rasmol_colour
                    -> get_colour_rgb
               -> write_restype_colours
                     -> get_rasmol_colour
                          -> get_colour_rgb
               -> write_conserv_colours
                     -> get_rasmol_colour
                          -> get_colour_rgb
         -> write_surfnet_gap
         -> order_results
         -> get_closest_sites
         -> write_html_gap
         -> write_db_gap
         -> close_html_page
         -> write_speedfill_arr
   -> write_script_tail
   -> write_pymol_tail
   -> rewrite_tmp_pml
         -> insert_pdb_file
               -> string_chomp (in common.c)
   -> close_surfnet_par
   -> close_and_reopen
   -> transfer_rasmol_file
   -> copy_file (in common.c)
   -> delete_temp_files
         -> call_system (in common.c)
   -> write_atom_depths
   -> save_vertices


*/


#include "common.h"
#include "speedfill.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
void call_system(char *command,int dos);
int check_file(char *file_name);
void copy_file(char *infile,char *outfile);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
char *form_filename(char *pdb_code,char *filename_template);
void free_tokens(char **token,int ntokens);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void remove_leading_blanks(char *string);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* Function prototypes from readpdb.c */
int r_find_pdb_file(char *pdb_code,char *pdb_template[NPDB_DIR],
                    char *file_name,int pdb_type);
void r_get_atom_details(char *input_line,int *atom_number,char *atom_name,
                        char *chain,char *res_name,char *res_num,
                        double *xval,double *yval,double *zval,
                        double *bval,double *occup,char *element,
                        int *rec_type,char wanted_chain);
FILE *r_open_pdb_file(char *file_name,int *gz);

/* In ras3d.c */
int get_grid_value(int grid_value,int *mskd,int contour_level);


/***********************************************************************

interpret_code  -  Extract the PDB code and chain-id or ligand-id
                   from the entered string

***********************************************************************/

void interpret_code(char *string,int len,char pdb_code[5],char *chain_id,
                    char ligand_id[6])
{
  char number_string[10];

  int domain_no, ligid_pos;

  /* If string long enough to contain the PDB code, take it to be the
     first 4 characters */
  if (len >= 4)
    {
      /* Extract the PDB code */
      strncpy(pdb_code,string,4);
      pdb_code[4] = '\0';

      /* Check if there is a chain-id or ligand identifier appended to
         the PDB code */
      ligid_pos = 5;
      if (len > 4)
        {
          /* If first character after the code is not an underscore, then
             have a chain-id */
          if (string[4] != '_')
            {
              /* Extract the chain identifier */
              *chain_id = string[4];
              if (*chain_id == '0')
                *chain_id = ' ';

              /* Ligand identifier starts one position beyond
                 chain-id */
              ligid_pos = ligid_pos + 1;

              /* Check if there's a domain number attached */
              if (len > 5)
                {
                  if (string[5] != '_')
                    {
                      /* Extract the domain number */
                      number_string[0] = string[5];
                      number_string[1] = '\0';
                      domain_no = atoi(number_string);

                      /* Ligand identifier starts one position beyond
                         domain number */
                      ligid_pos = ligid_pos + 1;
                    }
                }
            }
        }

      /* If input string contains the ligand identifier, then store
         that */
      if (len > 9)
        {
          /* Extract the ligand identifier */
          strncpy(ligand_id,string+ligid_pos,5);
          ligand_id[5] = '\0';
        }
    }
}
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           char *file_name,char *output_file_name,
                           struct parameters *params,char pdb_code[5],
                           char ch[MAX_CHAINS],char ligand_id[6],
                           int *nchain,char *log_name,char *rasmol_out_name,
                           char *out_dir,FILE **logfile_ptr,
                           double *grid_sep)
{
  char chain, out_name[FILENAME_LEN], number_string[12];

  char *string_ptr;

  int i, igap, ipos, itoken, len;
  int wanted;

  double g_sep;

  FILE *lgfile_ptr;

  /* Initialise variables */
  file_name[0] = '\0';
  g_sep = *grid_sep;
  pdb_code[0] = '\0';
  chain = ' ';
  ligand_id[0] = '\0';
  *nchain = 0;
  for (i = 0; i < MAX_CHAINS; i++)
    ch[i] = ' ';

  /* Initialise parameters */
  params->min_gap_radius = 1.4;
  params->max_gap_radius = 1.8;
  params->max_gap_span = 10.0;
  params->ntop_clusters = 4;
  params->indiv_atom_radii = FALSE;
  params->include_H_atoms = FALSE;
  params->include_waters = FALSE;
  params->select = FALSE;
  params->from_www = FALSE;
  params->generate_html = FALSE;
  params->buried_vertices_only = FALSE;
  params->vertex_colouring = STANDARD;
  params->write_surfnet_input = FALSE;
  params->read_vertices = FALSE;
  params->save_vertices = FALSE;
  params->generate_r3d = FALSE;
  strcpy(params->user_id,"0");
  params->profunc = FALSE;
  params->centroid = FALSE;
  params->store_atom_pointers = FALSE;
  params->col_by_conservation = FALSE;
  params->conservation_colours[0] = '\0';
  params->use_consurf_data = FALSE;
  params->have_rescons = FALSE;
  params->check_all_spheres = FALSE;
  params->no_rasmol = FALSE;
  params->no_protein = FALSE;
  params->no_gaps = FALSE;
  params->no_gif = FALSE;
  params->named_file = FALSE;
  params->depth_analysis = FALSE;
  params->save_grid = FALSE;
  params->filter_vertices = TRUE;
  for (igap = 0; igap < MAX_TOP_GAPS; igap++)
    params->show_gap[igap] = FALSE;
  params->delphi_file[0] = '\0';
  params->delete_temp_files = TRUE;
  params->jmol_resnames = FALSE;
  params->metatom_radius =  1.5;
  params->ligatom_radius =  -1.0;
  params->ligbond_radius =  0.0;
  strcpy(params->ligand_colour,"CPK");
  params->raster3d_size = 24;
  params->flip = FALSE;
  params->pdb_type = STANDARD_PDB;
  params->relinks = FALSE;
  params->sph_resno = FALSE;
  params->pymol = FALSE;
  params->pymol_www = FALSE;
  params->run_64 = FALSE;

  /* If nothing entered on the command line, show options available */
  if (num < 1)
    {
      printf("\n");
      printf("*** ERROR. Correct usage is:\n");
      printf("\n");
      printf("    speedfill { pdb_code[chain|ligand_no] | -f filename");
      printf(" [-c chain] }\n");
      printf("              [-min value]  [-max value]  [-span value]\n");
      printf("              [-r value]  [-ntop number]  [-b]  [-h]  [-w]\n");
      printf("              [-log]  [-www]  [-html]  [[-show number] ...]\n");
      printf("              [-outdir name]  [-bv]  [-col type]  [-surfnet]\n");
      printf("              [-noras]  [-noprot]  [-conser]  [-consurf]  ");
      printf("[-nogaps]  \n");
      printf("              [-nogif]  [-d]  [-centroid]\n");
      printf("              [-save]  [-read]  [-profunc]  [-r3d] ");
      printf("[-g grid_spacing]\n");
      printf("              [-e delphi_file]  [-nofilt]  [-r3d_size ntiles]");
      printf("  [UPLOAD | MCSG | AF]\n");
      printf("              [-havecons]  [-nodelete]  [-jmol]  [-flip]  "
             "[-wgrid]  [-sphres]  [-64]\n");
      printf("              [-pdbsum1]  [-relinks]  -[pymol]  "
             "[-pymol_www]  [-id user_id]\n ");
      printf("\n");
      printf("where\n");
      printf("   * pdb_code     = 4-character PDB identifier\n");
      printf("   * chain        = chain-id to be gap-filled\n");
      printf("   * ligand_no    = Number _MM_NN corresponding to");
      printf(" ligplotMM_NN.res file holding\n");
      printf("                    list of residues to be output\n");
      printf("                    (Default is 01_01)\n");
      printf("   * filename     = full name (and path) of PDB file\n");
      printf("   * -min value   = minimum gap-sphere radius ");
      printf("(default = %4.2f)\n",params->min_gap_radius);
      printf("   * -max value   = maximum gap-sphere radius ");
      printf("(default = %4.2f)\n",params->max_gap_radius);
      printf("   * -span value  = maximum gap-span ");
      printf("(default = %4.2f)\n",params->max_gap_span);
      printf("   * -ntop number = number of top gap regions to be output\n");
      printf("   * -b           = use radii given in B-value column\n");
      printf("   * -bv          = show buried vertices only\n");
      printf("   * -col type    = colouring of residues by\n");
      printf("                       0 = standard\n");      
      printf("                       1 = cpk\n");      
      printf("                       2 = residue type\n");      
      printf("                       3 = residue conservation\n");      
      printf("                       4 = buried/exposed\n");
      printf("                       5 = by chain\n");
      printf("                       6 = electrostatic potential (name of ");
      printf("Delphi input file must\n");
      printf("                           specified with -e option)\n");
      printf("                       7 = B-value\n");
      printf("   * -h           = don't exclude H-atoms\n");
      printf("   * -w           = don't exclude waters\n");
      printf("   * -show number = gap-region to be output to RasMol file\n");
      printf("   * -www         = called from Web page, so output RasMol");
      printf(" file to standard output\n");
      printf("   * -log         = write screen output to logfile,");
      printf(" speedfill.out\n");
      printf("   * -html        = generate HTML page showing gap analysis\n");
      printf("   * -surfnet     = generate lines for surfnet.par file\n");
      printf("   * -outdir name = directory to which files are to be ");
      printf("written\n");
      printf("   * -conser      = show protein residues coloured by");
      printf(" conservation score\n");
      printf("   * -consurf     = show protein residues coloured by ConSurf");
      printf(" conservation score\n");
      printf("   * -nogaps      = show protein residues and ligands only\n");
      printf("   * -nogif       = don't output gif file to speedfill.html\n");
      printf("   * -noras       = generate html file only (ie no RasMol)\n");
      printf("   * -noprot      = no protein coords to be written to output");
      printf(" RasMol file)\n");
      printf("   * -d           = perform depth analysis of all atoms\n");
      printf("   * -save        = save vertices to the speedfill.depth ");
      printf("file\n");
      printf("   * -read        = read in the vertices from the ");
      printf("speedfill.depth file\n");
      printf("   * -wgrid       = write out the grid to speedfill.arr ");
      printf("file\n");
      printf("   * -r3d         = write out a .r3d file for rendering");
      printf(" by Raster3D\n");
      printf("   * -g grid_spacing = default is 0.75\n");
      printf("   * -profunc     = Flag indicating that called by ProFunc\n");
      printf("   * -e delphi_file = Colour vertices by electrostatic ");
      printf("potential read in\n");
      printf("                    from the Delphi file given\n");
      printf("   * -nofilt      = No filtering of vertices prior to output ");
      printf("to RasMol file\n");
      printf("   * -r3d_size ntiles = Size out output Raster3D image (no. ");
      printf("of tiles\n");
      printf("   * UPLOAD       = Uploaded PDB file for PDBsum generation\n");
      printf("   * MCSG         = MCSG structure\n");
      printf(      " AF       = Alpha Fold structure\n");
      printf("   * -centroid    = Write out the centroid of each gap region "
             "to the RasMol script file");
      printf("   * -havecons    = Flag indicating that -read file holds "
             "residue conservation\n");
      printf("   * -nodelete    = Don't delete temporary files\n");
      printf("   * -jmol        = Write residue names for gap regions as "
             "VTA, VTB, etc for residue\n"
             "                    conservation colouring\n");
      printf("   * -flip        = Rotate starting coords through 180 "
             "degrees about x-axis\n");
      printf("   * -sphres      = list spheres as separate residues\n");
      printf("   * -pdbsum1     = PDBsum1 structure\n");
      printf("   * -relinks     = Relative links in PDBsum1 pages\n");
      printf("   * -pymol       = Write out PyMOL scripts\n");
      printf("   * -pymol_www   = Write out PyMOL scripts for web output\n");
      printf("   * -id user_id  = User id (tmp directory holding files)\n");
      printf("   * -64          = run on 64-bit processor\n");
      printf("\n");
      printf("\n");
      printf("Examples:-\n");
      printf("\n");
      printf("   speedfill -f pdb8tim.ent -c A\n");
      printf("   speedfill 1aaw_01_01\n");
      printf("\n");
      printf("\n");
      exit(1);
    }

  /* Loop through the command arguments */
  for (itoken = 1; itoken < num + 1; itoken++)
    {
      /* Get the length of this token */
      len = strlen(strptrs[itoken]);

      /* Check if this is the number of clusters to be written out parameter */
      if (!strcmp(strptrs[itoken],"-ntop"))
        {
          /* Get the number of top clusters from the next command-line
             parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->ntop_clusters = atoi(number_string);
              itoken++;

              /* Check that maximum number of top clusters not exceeded */
              if (params->ntop_clusters > MAX_TOP_GAPS)
                {
                  printf("*** Warning. Maximum number of top gap");
                  printf(" regions exceeded. Reset to %d\n",
                         MAX_TOP_GAPS);
                  params->ntop_clusters = MAX_TOP_GAPS;
                }
            }
        }

      /* Check if this is the minimum gap_sphere radius */
      else if (!strcmp(strptrs[itoken],"-min"))
        {
          /* Get the radius from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->min_gap_radius = atof(number_string);
              itoken++;
            }
        }

      /* Check if this is the maximum gap_sphere radius */
      else if (!strcmp(strptrs[itoken],"-max"))
        {
          /* Get the radius from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->max_gap_radius = atof(number_string);
              itoken++;
            }
        }

      /* Check if this is the maximum gap_span */
      else if (!strcmp(strptrs[itoken],"-span"))
        {
          /* Get the value from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->max_gap_span = atof(number_string);
              itoken++;
            }
        }

      /* Check for vertex colouring scheme */
      else if (!strcmp(strptrs[itoken],"-col"))
        {
          /* Get the value from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->vertex_colouring = atoi(number_string);
              if (params->vertex_colouring < 0 ||
                  params->vertex_colouring > NCOLOUR_SCHEMES - 1)
                params->vertex_colouring = STANDARD;
              itoken++;
            }
        }

      /* Check if this identifies a gap-region to be output */
      else if (!strcmp(strptrs[itoken],"-show"))
        {
          /* Get the value from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              /* Check if all gaps to be shown */
              if (!strcmp(strptrs[itoken + 1],"all"))
                {
                  /* Loop over all gaps to set to true */
                  for (igap = 0; igap < NTOP_BS; igap++)
                    params->show_gap[igap] = TRUE;

                  /* Set selection flag */
                  params->select = TRUE;
                }
              
              /* Otherwise, determine gap region to be shown */
              else
                {
                  /* Determine gap region to be shown */
                  strcpy(number_string,strptrs[itoken + 1]);
                  igap = atoi(number_string);

                  /* If within limits, store the gap region to be shown */
                  if (igap > 0 && igap <= MAX_TOP_GAPS)
                    {
                      params->show_gap[igap - 1] = TRUE;
                      params->select = TRUE;
                    }
                }
              
              /* Increment token-counter */
              itoken++;
            }
        }

      /* Check if the -www flag set indicating that called from WWW page */
      else if (!strcmp(strptrs[itoken],"-www"))
        {
          /* Set RasMol output file to go to standard output */
          rasmol_out_name[0] = '\0';

          /* Set screen output to go to a logfile */
          strcpy(log_name,"speedfill.out");

          /* Set flag */
          params->from_www = TRUE;
        }

      /* Check if the -noras flag set indicating that called from WWW page
         to generate the HTML analysis of binding sites only (ie no
         RasMol output required) */
      else if (!strcmp(strptrs[itoken],"-noras"))
        {
          /* Set screen output to go to a logfile */
          strcpy(log_name,"speedfill.out");

          /* Set flag */
          params->generate_html = TRUE;
          params->no_rasmol = TRUE;
        }

      /* Check if the -noprot flag set indicating that protein and ligand
         coords aren't to be written out to the RasMol output file */
      else if (!strcmp(strptrs[itoken],"-noprot"))
        {
          params->no_protein = TRUE;
        }

      /* Check if the -log flag set indicating that screen output is to
         go to a log file, speedfill.out */
      else if (!strcmp(strptrs[itoken],"-log"))
        {
          /* Set screen output to go to a logfile */
          strcpy(log_name,"speedfill.out");
        }

      /* Check if the -html flag set indicating that HTML summary of gap
         region analysis is required */
      else if (!strcmp(strptrs[itoken],"-html"))
        {
          /* Set flag */
          params->generate_html = TRUE;
        }

      /* Check if the -conser flag set indicating that protein residues
         are to be coloured by residue conservation */
      else if (!strcmp(strptrs[itoken],"-conser"))
        {
          /* Set flag */
          params->col_by_conservation = TRUE;
        }

      /* Check if the -connsurf flag set indicating that protein residues
         are to be coloured by residue conservation obtained from the
         ConSurf data files */
      else if (!strcmp(strptrs[itoken],"-consurf"))
        {
          /* Set flag */
          params->col_by_conservation = TRUE;
          params->use_consurf_data = TRUE;
        }

      /* Check if the -nogaps flag set indicating that no gap regions
         are to be generated and only protein residues and ligands
         are to be displayed */
      else if (!strcmp(strptrs[itoken],"-nogaps"))
        {
          /* Set flag */
          params->no_gaps = TRUE;
        }

      /* Check if the -nogif flag set indicating that no gif
         files are to be written out to the speedfill.html file */
      else if (!strcmp(strptrs[itoken],"-nogif"))
        {
          /* Set flag */
          params->no_gif = TRUE;
        }

      /* Check if the -bv flag set indicating that only buried vertices
         are to be output for the gap regions */
      else if (!strcmp(strptrs[itoken],"-bv"))
        {
          /* Set flag */
          params->buried_vertices_only = TRUE;
        }

      /* Check if the -surfnet flag set indicating that surfnet.par file
         to be generated */
      else if (!strcmp(strptrs[itoken],"-surfnet"))
        {
          /* Set flag */
          params->write_surfnet_input = TRUE;
        }

      /* Check whether vertices are to be saved to the speedfill.depth
         file */
      else if (!strcmp(strptrs[itoken],"-save"))
        {
          /* Set flag */
          params->save_vertices = TRUE;
        }

      /* Check whether vertices are to be read in from the speedfill.depth
         file */
      else if (!strcmp(strptrs[itoken],"-read"))
        {
          /* Set flag */
          params->read_vertices = TRUE;
        }

      /* Check whether grid array is to be saved to the speedfill.arr
         file */
      else if (!strcmp(strptrs[itoken],"-wgrid"))
        {
          /* Set flag */
          params->save_grid = TRUE;
        }

      /* Check whether a .r3d file is required */
      else if (!strcmp(strptrs[itoken],"-r3d"))
        {
          /* Set flag */
          params->generate_r3d = TRUE;
        }

      /* Check whether we are in the ProFunc system */
      else if (!strcmp(strptrs[itoken],"-profunc"))
        {
          /* Set flag */
          params->profunc = TRUE;
        }

      /* Check whether centroid to be written out */
      else if (!strcmp(strptrs[itoken],"-centroid"))
        {
          /* Set flag */
          params->centroid = TRUE;
        }

      /* Check if the -read file contains residue conservation */
      else if (!strcmp(strptrs[itoken],"-havecons"))
        {
          /* Set flag */
          params->have_rescons = TRUE;
        }

      /* Check if temporary files to be retained */
      else if (!strcmp(strptrs[itoken],"-nodelete"))
        {
          /* Set flag */
          params->delete_temp_files = FALSE;
        }

      /* Check if special residue names for JMol are required */
      else if (!strcmp(strptrs[itoken],"-jmol"))
        {
          /* Set flag */
          params->jmol_resnames = TRUE;
        }

      /* Check whether coords are to be flipped at start */
      else if (!strcmp(strptrs[itoken],"-flip"))
        {
          /* Set flag */
          params->flip = TRUE;
        }

      /* Check for -outdir option preceding the name of output directory */
      else if (!strcmp(strptrs[itoken],"-outdir"))
        {
          /* Get the directory name from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(out_dir,strptrs[itoken + 1]);
              itoken++;
            }
        }

      /* Check if the -b flag set defining individual atomic radii to be
         taken from the B-value column of the PDB file */
      else if (!strcmp(strptrs[itoken],"-b"))
        {
          /* Set flag */
          params->indiv_atom_radii = TRUE;
        }

      /* Check for -c option preceding the chain-identifier */
      else if (!strcmp(strptrs[itoken],"-c") ||
               !strcmp(strptrs[itoken],"-chain"))
        {
          /* Get the chain-id from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              chain = strptrs[itoken + 1][0];
              itoken++;

              /* If can, then store in list of chains */
              if (*nchain < MAX_CHAINS - 1)
                {
                  ch[*nchain] = chain;
                  (*nchain)++;
                }
            }
        }

      /* Check for -e option preceding the name of the Delphi file
         to be read to colour surface by electrostatic potential */
      else if (!strcmp(strptrs[itoken],"-e"))
        {
          /* Get the filename from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(params->delphi_file,strptrs[itoken + 1]);
              itoken++;

              /* Set vertex colouring mode to electrostatic */
              params->vertex_colouring = ELECTROSTATIC;
            }
        }

      /* Check for -f option preceding the full filename */
      else if (!strcmp(strptrs[itoken],"-f"))
        {
          /* Get the filename from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(file_name,strptrs[itoken + 1]);
              itoken++;

              /* Set flag that this is a prespecified file, and not
                 one in the standard PDB location */
              params->named_file = TRUE;

              /* Extract the PDB code from the file name */
              string_ptr = strstr(file_name,".new");
              if (string_ptr == NULL)
                string_ptr = strstr(file_name,".pdb");

              /* If have the right type of extension, extract PDB,
                 if possible */
              if (string_ptr != NULL)
                {
                  /* Get the extension's start position */
                  ipos = string_ptr - file_name;

                  /* If it looks OK, extract the PDB code */
                  if (ipos - 4 > -1)
                    {
                      /* Save the PDB code */
                      strncpy(pdb_code,file_name + ipos - 4,4);
                      pdb_code[4] = '\0';

                    }
                }
            }
        }

      /* Check for -g option preceding the grid-spacing */
      else if (!strcmp(strptrs[itoken],"-g"))
        {
          /* Get the spacing from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              *grid_sep = atof(number_string);
              if (*grid_sep == 0.0)
                *grid_sep = g_sep;
              itoken++;
            }
        }

      /* Check if the -h flag set requiring inclusion of H-atoms */
      else if (!strcmp(strptrs[itoken],"-h"))
        {
          /* Set flag */
          params->include_H_atoms = TRUE;
        }

      /* Check if the -w flag set requiring inclusion of waters */
      else if (!strcmp(strptrs[itoken],"-w"))
        {
          /* Set flag */
          params->include_waters = TRUE;
        }

      /* Check if the -d flag set requiring depth analysis */
      else if (!strcmp(strptrs[itoken],"-d"))
        {
          /* Set flag */
          params->depth_analysis = TRUE;
        }

      /* Check if the -nofilt flag set requiring that vertices not be
         filtered */
      else if (!strcmp(strptrs[itoken],"-nofilt"))
        {
          /* Set flag */
          params->filter_vertices = FALSE;
        }

      /* Check for -r3d_size option preceding the Raster3D image size */
      else if (!strcmp(strptrs[itoken],"-r3d_size"))
        {
          /* Get the size from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(number_string,strptrs[itoken + 1]);
              params->raster3d_size = atoi(number_string);
              itoken++;
            }
        }

      /* Check for the user id */
      else if (!strcmp(strptrs[itoken],"-id"))
        {
          /* Get the id from the next command-line parameter */
          if (itoken + 1 < num + 1)
            {
              strcpy(params->user_id,strptrs[itoken + 1]);
              itoken++;
            }
        }

      /* Check for the UPLOAD option */
      else if (!strcmp(strptrs[itoken],"UPLOAD"))
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

      /* Check for the -relinks option */
      else if (!strcmp(strptrs[itoken],"-relinks"))
        params->relinks = TRUE;

      /* Check for the PyMOL option */
      else if (!strcmp(strptrs[itoken],"-pymol"))
        params->pymol = TRUE;

      /* Check for the PyMOL web option */
      else if (!strcmp(strptrs[itoken],"-pymol_www"))
        params->pymol_www = TRUE;

      /* Check for the sphres option */
      else if (!strcmp(strptrs[itoken],"-sphres"))
        params->sph_resno = TRUE;

      /* Otherwise, take token to be PDB code plus possible chain
         identifier or ligand-id appended */
      else if (len > 3 && len < 12)
        {
          /* Set flag */
          wanted = FALSE;

          /* Check that first character is number or letter */
          if (strptrs[itoken][0] >= '0' && strptrs[itoken][0] <= '9')
            wanted = TRUE;
          else if (strptrs[itoken][0] >= 'a' && strptrs[itoken][0] <= 'z')
            wanted = TRUE;
          else if (len == 4 && itoken == 1 &&
                   strptrs[itoken][0] >= 'A' && strptrs[itoken][0] <= 'Z')
            wanted = TRUE;

          /* If could be PDB code, then interpret */
          if (wanted == TRUE)
            {
              /* Interpret the code */
              interpret_code(strptrs[itoken],len,pdb_code,&chain,ligand_id);

              /* If can, then store in list of chains */
              if (chain != ' ' && *nchain < MAX_CHAINS - 1)
                {
                  ch[*nchain] = chain;
                  (*nchain)++;
                }
            }
        }
    }

  /* Open the output log file */
  if (log_name[0] != '\0')
    {
      if (out_dir[0] != '\0')
        {
          strcpy(out_name,out_dir);
          strcat(out_name,DIR_SEP);
          strcat(out_name,log_name);
        }
      else
        strcpy(out_name,log_name);
    }
  else
    out_name[0] = '\0';
  lgfile_ptr = open_output_file(out_name,TRUE);
  *logfile_ptr = lgfile_ptr;

  /* Show the file to be processed */
  fprintf(lgfile_ptr,"Run parameters:-\n");
  fprintf(lgfile_ptr,"\n");
  fprintf(lgfile_ptr,"   PDB file:-  %s\n",file_name);

  /* If no chains identified, have just chain blank */
  if (*nchain == 0)
    {
      ch[*nchain] = ' ';
      (*nchain)++;
    }
      
  /* Show chain(s) */
  if (*nchain > 1 || ch[0] != ' ')
    {
      if (*nchain == 1)
        fprintf(lgfile_ptr,"   Chain-id:-  %c\n",ch[0]);
      else
        {
          fprintf(lgfile_ptr,"   Chain-ids:-  ");
          for (i = 0; i < *nchain; i++)
            {
              fprintf(lgfile_ptr,"%c",ch[i]);
              if (i < *nchain - 1)
                fprintf(lgfile_ptr,", ");
            }
          fprintf(lgfile_ptr,"\n");
        }
    }

  /* If gap regions not individually selected, then show required number
     of top ones */
  if (params->select == FALSE)
    {
      for (igap = 0; igap < params->ntop_clusters; igap++)
        params->show_gap[igap] = TRUE;
    }

  /* Otherwise, work out which is the top cluster to be shown */
  else
    {
      /* Initialise */
      params->ntop_clusters = 1;

      /* Loop through all possible gap regions and save whichever is the
         highest selected one */
      for (igap = 0; igap < MAX_TOP_GAPS; igap++)
        if (params->show_gap[igap] == TRUE)
          params->ntop_clusters = igap + 1;
    }

  /* Determine whether grid of atom pointers required for obtaining the
     atom associated with each grid-point */
  if (params->vertex_colouring != STANDARD)
    {
      /* Only have non-standard colouring for the buried-vertices-only
         option */
      if (params->buried_vertices_only == FALSE &&
          (params->vertex_colouring != BURIED_EXPOSED &&
           params->save_vertices == FALSE))
        params->vertex_colouring = STANDARD;

      /* Otherwise, set flag to ensure that array of atom pointers
         corresponding to the grid-points is created */
      else
        params->store_atom_pointers = TRUE;
    }

  /* If generating the HTML pages, or performing the depth analysis,
     then need to store the atom pointers */
  if (params->generate_html == TRUE || params->depth_analysis == TRUE ||
      params->pdb_type == PDBSUM1)
    params->store_atom_pointers = TRUE;

  /* If checking all spheres in gap generation, then need to store atom
     pointers */
  if (params->check_all_spheres == TRUE)
    params->store_atom_pointers = TRUE;

  /* If reading in the vertices, then don't want to also save them */
  if (params->read_vertices == TRUE)
    params->save_vertices = FALSE;

  /* If saving vertices, then need to store atom pointers */
  if (params->save_vertices == TRUE)
    params->store_atom_pointers = TRUE;

  /* If colouring by electrostatic potential selected, check that
     name of Delphi file supplied */
  if (params->vertex_colouring == ELECTROSTATIC &&
      params->read_vertices == FALSE && params->delphi_file[0] == '\0' )
    params->vertex_colouring = STANDARD;

  /* Check that JMol residue names only generated for residue
     conservation */
  if (params->jmol_resnames == TRUE &&
      params->vertex_colouring != CONSERVATION)
    params->jmol_resnames = FALSE;

  /* Show the run-time parameters */
  fprintf(lgfile_ptr,"\n");
  fprintf(lgfile_ptr,"   Minimum gap-sphere radius   %8.3f\n",
          params->min_gap_radius);
  fprintf(lgfile_ptr,"   Maximum gap-sphere radius   %8.3f\n",
          params->max_gap_radius);
  fprintf(lgfile_ptr,"   Maximum gap-span            %8.3f\n",
          params->max_gap_span);
  fprintf(lgfile_ptr,"   Top clusters to process     %4d\n",
          params->ntop_clusters);
  fprintf(lgfile_ptr,"\n");
}
/***********************************************************************

prepend_dir  -  Prepend output directory name to the names of the input
                and output files

***********************************************************************/

void prepend_dir(char *out_dir,char *arr_name,char *html_name,
                 char *depth_name,char *out_name,char *pdb_out_name,
                 char *r3d_name,char *script_out_name,char *surfnet_name,
                 char *rasmol_out_name,char *clefts_pml,char *tmp_pml,
                 char *dummy_name,struct parameters *params)
{
  char tmp_name[FILENAME_LEN];
 
  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,html_name);
  strcpy(html_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,arr_name);
  strcpy(arr_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,depth_name);
  strcpy(depth_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,dummy_name);
  strcpy(dummy_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,out_name);
  strcpy(out_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,pdb_out_name);
  strcpy(pdb_out_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,r3d_name);
  strcpy(r3d_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,script_out_name);
  strcpy(script_out_name,tmp_name);

  strcpy(tmp_name,out_dir);
  strcat(tmp_name,DIR_SEP);
  strcat(tmp_name,surfnet_name);
  strcpy(surfnet_name,tmp_name);

  if (rasmol_out_name[0] != '\0')
    {
      strcpy(tmp_name,out_dir);
      strcat(tmp_name,DIR_SEP);
      strcat(tmp_name,rasmol_out_name);
      strcpy(rasmol_out_name,tmp_name);
    }

  if (params->delphi_file[0] != '\0')
    {
      strcpy(tmp_name,out_dir);
      strcat(tmp_name,DIR_SEP);
      strcat(tmp_name,params->delphi_file);
      strcpy(params->delphi_file,tmp_name);
    }

  if (params->pymol_www == TRUE)
    {
      strcpy(tmp_name,out_dir);
      strcat(tmp_name,DIR_SEP);
      strcat(tmp_name,clefts_pml);
      strcpy(clefts_pml,tmp_name);

      strcpy(tmp_name,out_dir);
      strcat(tmp_name,DIR_SEP);
      strcat(tmp_name,tmp_pml);
      strcpy(tmp_pml,tmp_name);
    }
}
/***********************************************************************

form_file_names  -  Form the input filenames of the PDB file and ligplot
                    output file giving list of ligand residues and
                    interacting protein residues

***********************************************************************/

void form_file_names(char *pdbsum_dir,char *pdb_code,char chain,
                     char *ligand_id,char *pdb_name,char *reslist_name,
                     char *res_cons_name,char *dfn_name,char *chains_dat,
                     char *out_dir,int pdb_type,int profunc,int consurf,
                     char *consurf_dir,char *profunc_dir)
{
  char chain2[2], middle[3];
  int ipdb;
  int found;

  FILE *file_ptr;

  /* Initialise variables */
  found = FALSE;
  if (pdb_name[0] != '\0')
    found = TRUE;

  /* Link to the appropriate LIGPLOT .res file */
  strcpy(reslist_name,pdbsum_dir);
  strcat(reslist_name,"ligplot");
  strcat(reslist_name,ligand_id);
  strcat(reslist_name,".res");

  /* Form name of residue conservation file */
  if (profunc == TRUE)
    {
      /* Residue conservation file */
      if (out_dir[0] != '\0')
        {
          strcpy(res_cons_name,out_dir);
          strcat(res_cons_name,DIR_SEP);
          strcat(res_cons_name,"conserv.dat");

          /* See if file exists */
          found = check_file(res_cons_name);

          /* If not found, look in ProFunc directory */
          if (found == FALSE)
            {
              /* Form name of the file */
              strcpy(res_cons_name,profunc_dir);
              strcat(res_cons_name,DIR_SEP);
              strcat(res_cons_name,"profunc");
              strcat(res_cons_name,pdb_code);
              strcat(res_cons_name,DIR_SEP);
              strcat(res_cons_name,"conserv.dat");
            }
        }
      else
        strcpy(res_cons_name,"conserv.dat");

      /* RasMol definition file */
      if (out_dir[0] != '\0')
        {
          strcpy(dfn_name,out_dir);
          strcat(dfn_name,DIR_SEP);
          strcat(dfn_name,"rasmol.dfn");
        }
      else
        strcpy(dfn_name,"rasmol.dfn");

      /* Chain definitions file */
      if (out_dir[0] != '\0')
        {
          strcpy(chains_dat,out_dir);
          strcat(chains_dat,DIR_SEP);
          strcat(chains_dat,"chains.dat");
        }
      else
        strcpy(chains_dat,"chains.dat");
    }

  else
    {
      /* If using ConSurf residue conservation file, point to 
         appropariate directory */
      if (consurf == TRUE)
        {
          /* Get middle part of PDB code */
          strncpy(middle,pdb_code + 1,2);
          middle[2] = '\0';

          /* Form name of the corresponding ConSurf file */
          strcpy(res_cons_name,consurf_dir);
          strcat(res_cons_name,DIR_SEP);
          strcat(res_cons_name,middle);
          strcat(res_cons_name,DIR_SEP);
          strcat(res_cons_name,pdb_code);
          if (chain == ' ')
            strcat(res_cons_name,"_");
          else
            {
              chain2[0] = chain;
              chain2[1] = '\0';
              strcat(res_cons_name,chain2);
            }
          strcat(res_cons_name,".gradesPE");
        }

      /* Otherwise, use the conserv.dat file in the PDBsum directory */
      else
        {
          strcpy(res_cons_name,pdbsum_dir);
          strcat(res_cons_name,"conserv.dat");
        }

      /* RasMol definition file */
      strcpy(dfn_name,pdbsum_dir);
      if (pdb_type == PDBSUM1)
        strcat(dfn_name,"../newsec/");
      strcat(dfn_name,"rasmol.dfn");

      /* Chain definitions file */
      strcpy(chains_dat,pdbsum_dir);
      if (pdb_type == PDBSUM1)
        strcat(chains_dat,"../newsec/");
      strcat(chains_dat,"chains.dat");
    }
}
/***********************************************************************

create_rasdef_entry  -  Create a rasdef entry

***********************************************************************/

struct rasdef *create_rasdef_entry(struct rasdef **fst_rasdef_ptr,
                                   struct rasdef **lst_rasdef_ptr,
                                   int type,int number,char *res_name,
                                   char *res_num,char chain,char *colour)
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
  rasdef_ptr->type = type;
  rasdef_ptr->number = number;
  strcpy(rasdef_ptr->res_name,res_name);
  strcpy(rasdef_ptr->res_num,res_num);
  rasdef_ptr->chain = chain;
  strcpy(rasdef_ptr->colour,colour);

  /* Initialise other fields */
  rasdef_ptr->first_atom_ptr = NULL;
  rasdef_ptr->last_atom_ptr = NULL;

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

void get_ligand_residues(char *input_line,struct rasdef **fst_rasdef_ptr,
                         struct rasdef **lst_rasdef_ptr,int number)
{
  char chain, res_name[4], res_num[6];
  char colour[COLNAM_LEN], number_string[10];
  char token_string[TOKEN_LEN + 1], tmp_string[TOKEN_LEN + 1];
  char *string_ptr;

  int len, lpos, token_len, type;
  int name_end;
  int done;

  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  /* Initialise variables */
  chain = ' ';
  done = FALSE;
  type = LIGAND;
  colour[0] = '\0';

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

          /* Remove any trailing blanks from the name */
          len = string_truncate(res_name,4);

          /* Right justify, if not 3 characters long */
          if (len != 3)
            {
              strcpy(tmp_string,res_name);
              if (len == 2)
                {
                  strcpy(res_name," ");
                  strcat(res_name,tmp_string);
                }
              else if (len == 1)
                {
                  strcpy(res_name,"  ");
                  strcat(res_name,tmp_string);
                }
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

          /* Create a rasdef entry for this ligand residue */
          rasdef_ptr
            = create_rasdef_entry(&first_rasdef_ptr,
                                  &last_rasdef_ptr,type,number,
                                  res_name,res_num,chain,colour);
        }
    }

  /* Return current pointers */
  *fst_rasdef_ptr = first_rasdef_ptr;
  *lst_rasdef_ptr = last_rasdef_ptr;
}
/***********************************************************************

get_dfn_data  -  Read in and store the RasMol definitions from the
                 rasmol.dfn file

***********************************************************************/

void get_dfn_data(char *dfn_name,struct rasdef **fst_rasdef_ptr,
                  struct rasdef **lst_rasdef_ptr)
{
  char input_line[LINELEN + 1];
  char chain, res_name[4], res_num[6];
  char colour[COLNAM_LEN];

  int len, line, nchain, nligand, nmetal, type;
  int first;

  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;

  FILE *file_ptr;

  /* Initialise */
  first = TRUE;
  line = 0;
  nchain = 0;
  nligand = 0;
  nmetal = 0;

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = NULL;
  last_rasdef_ptr = NULL;

  /* Open the rasmol.dfn input file */
  if ((file_ptr = fopen(dfn_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

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
                                      res_name,res_num,chain,colour);
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
              strncpy(res_num,input_line+1,4);
              res_num[4] = '\0';

              /* Get the metal colour */
              strncpy(colour,input_line+5,COLNAM_LEN - 1);
              colour[COLNAM_LEN - 1] = '\0';

              /* Create a rasdef entry */
              rasdef_ptr
                = create_rasdef_entry(&first_rasdef_ptr,
                                      &last_rasdef_ptr,type,nmetal,
                                      res_name,res_num,chain,colour);
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
                                      res_name,res_num,chain,colour);
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
                                      res_name,res_num,chain,colour);
            }
        }
    }

  /* Return pointer to start of RasMol definitions */
  *fst_rasdef_ptr = first_rasdef_ptr;
  *lst_rasdef_ptr = last_rasdef_ptr;
}
/***********************************************************************

check_if_want_residue  -  Check whether the residue is wanted, or is a
                          ligand residue

***********************************************************************/

int check_if_want_residue(char in_chain,char res_name[4],
                          char chain[MAX_CHAINS],int nchain,
                          int *in_ligand,int chain_end,char athet)
{
  int ichain;
  int wanted;

  /* Initialise variables */
  *in_ligand = FALSE;
  wanted = FALSE;

  /* If all chains to be includes, then residue is wanted */
  if (nchain ==  1 && chain[0] == ' ')
    wanted = TRUE;

  /* Otherwise, check chains in list */
  else
    {
      /* Check whether it belongs to the required chain */
      for (ichain = 0; ichain < nchain; ichain++)
        {
          if (in_chain == chain[ichain])
            wanted = TRUE;
        }
    }

  /* If this is a HETATM, then only consider it as a possible ligand
     if we already have our required chain(s) */
  if (athet == 'H')
    {
      /* If required chain has ended, then this may be a ligand */
      if (chain_end == TRUE)
        {
          *in_ligand = TRUE;
          wanted = TRUE;
        }
    }

  /* Exclude waters */
  if (!strcmp(res_name,"HOH"))
    wanted = FALSE;

  /* Return whether given atom is wanted */
  return(wanted);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      char res_name[4],char res_num[6],
                                      char chain)
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
  residue_ptr->natoms = 0;
  residue_ptr->type = PROTEIN;
  residue_ptr->in_range = TRUE;
  residue_ptr->conservation_score = -1.0;
  residue_ptr->checked = FALSE;
  residue_ptr->first_atom_ptr = NULL;
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return new residue pointer */
  return(residue_ptr);
}
/***********************************************************************

get_atom_type  -  Get atom type from its name

***********************************************************************/

int get_atom_type(char atom_name[5])
{
  int type;

  /* Set type to default value */
  type = UNKNOWN;

  /* Nitrogen */
  if (atom_name[1] == 'N' && atom_name[0] != 'Z' && atom_name[0] != 'M')
    type = NITROGEN;

  /* Carbon */
  else if (atom_name[1] == 'C')
    type = CARBON;

  /* Oxygen */
  else if (atom_name[1] == 'O' && atom_name[0] != 'C' && atom_name[0] != 'M')
    type = OXYGEN;

  /* Sulphur */
  else if (atom_name[1] == 'S')
    type = SULPHUR;

  /* Phosphorus */
  else if (atom_name[1] == 'P')
    type = PHOSPHORUS;

  /* Iron */
  else if (!strncmp(atom_name,"FE",2))
    type = IRON;

  /* Return the atom type */
  return(type);
}
/***********************************************************************

get_atom_radius  -  Get atom's van der Waals radius

***********************************************************************/

double get_atom_radius(char *elmnt,int ca_only)
{
  char element[3];

  int atom_name, itype;
  int found;

  static char *elem[] = {
    "XX", " N", " C", " O", " S", " P", "FE", "ca"
  };

  static double vdw_radius[NATOM_TYPES] = {
    2.00,  /*  0. Unknown - default */
    1.65,  /*  1. Nitrogen */
    1.87,  /*  2. Carbon */
    1.40,  /*  3. Oxygen */
    1.85,  /*  4. Sulphur */
    1.90,  /*  5. Phosphorus */
    1.10,  /*  6. Iron */
    6.00   /*  7. CA-only */
  };

  /* Store the element type */
  strcpy(element,elmnt);

  /* For CA-only chains, want CAs represented by a single large radius */
  if (ca_only == TRUE)
    strcpy(element,"ca");

  /* Determine the atom-type from its element identifier */
  atom_name = 0;
  found = FALSE;
  for (itype = 1; itype < NATOM_TYPES && found == FALSE; itype++)
    {
      /* Check for match */
      if (!strcmp(element,elem[itype]))
        {
          /* Store type and set flag */
          atom_name = itype;
          found = TRUE;
        }
    }

  /* Return the radius */
  return(vdw_radius[atom_name]);
}
/***********************************************************************

create_atom_record  -  Create a new atom record

***********************************************************************/

struct atom *create_atom_record(struct atom **fst_atom_ptr,
                                struct atom **lst_atom_ptr,
                                struct residue *residue_ptr,
                                int atom_number,char *atom_name,
                                char *element,double x,double y,double z,
                                double bvalue,double radius)
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
      printf("***        Program rasurf terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  atom_ptr->atom_number = atom_number;
  strcpy(atom_ptr->atom_name,atom_name);
  strcpy(atom_ptr->element,element);
  atom_ptr->coords[0] = x;
  atom_ptr->coords[1] = y;
  atom_ptr->coords[2] = z;
  atom_ptr->bvalue = bvalue;
  atom_ptr->radius = radius;
  atom_ptr->residue_ptr = residue_ptr;

  /* Initialise remaining fields */
  atom_ptr->depth = 0.0;
  atom_ptr->icluster = -1;
  atom_ptr->ndepth = 0;

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

check_for_duplicate  -  Check for an earlier copy of the current atom

***********************************************************************/

int check_for_duplicate(struct atom *first_atom_ptr,
                        struct atom *first_lig_atom_ptr,
                        char *atom_name,char *res_name,char *res_num,
                        char chain)
{
  int have_atom;

  struct atom *atom_ptr;

  /* Initialise variables */
  have_atom = FALSE;

  /* Get the first of the stored protein atom records */
  atom_ptr = first_atom_ptr;

  /* Loop through all atoms */
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

  /* Get the first of the stored ligand atom records */
  atom_ptr = first_lig_atom_ptr;

  /* Loop through all atoms */
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
      else if (rasdef_ptr->chain == chain && rasdef_ptr->res_name[0] == '\0')
        assigned_rasdef_ptr = rasdef_ptr;

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }

  /* Return the object to which this residue is to be assigned */
  return(assigned_rasdef_ptr);
}
/***********************************************************************

read_pdb  -  Scan through the PDB file to get the x-, y- z-coord limits
             of the atoms required

***********************************************************************/

int read_pdb(char *pdb_name,struct atom **fst_atom_ptr,
             struct atom **fst_lig_atom_ptr,int *nlig,
             struct residue **fst_residue_ptr,int *nres,
             struct parameters params,char *chain,int *nchn,
             struct rasdef **fst_rasdef_ptr,
             struct rasdef *last_rasdef_ptr,FILE *logfile_ptr,
             int all_chains,int flip)
{
  char input_line[LINELEN + 1];
  char atom_name[5], in_chain, res_name[4], res_num[6], occup[4];
  char element[3], wanted_chain;
  char last_chain, last_res_name[4], last_res_num[6];
  char rs_name[4], rs_num[6];
  char colour[COLNAM_LEN];

  int atom_number, natoms, nlig_atoms, nlines, nread, nresid;
  int ichain, len, nchain, rec_type, type;
  int chain_end, have_atom, have_chain, have_rasdefs, keep_reading;
  int ca_only, eof, gzipped, hydrogen, in_ligand, wanted;

  double bvalue, occupancy, radius, x, y, z;

  struct atom *atom_ptr, *first_atom_ptr, *last_atom_ptr;
  struct atom *first_lig_atom_ptr, *last_lig_atom_ptr;
  struct rasdef *rasdef_ptr;
  struct rasdef *first_rasdef_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise counts */
  fprintf(logfile_ptr,"Reading PDB file: %s ...\n",pdb_name);
  chain_end = FALSE;
  eof = FALSE;
  have_chain = FALSE;
  keep_reading = TRUE;
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  last_chain = '\0';
  natoms = 0;
  nchain = *nchn;
  nlig_atoms = 0;
  nread = 0;
  nlines = 0;
  nresid = 0;
  wanted_chain = '\0';

  /* Initialise pointers */
  first_lig_atom_ptr = NULL;
  last_lig_atom_ptr = NULL;
  first_atom_ptr = NULL;
  last_atom_ptr = NULL;
  residue_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  rasdef_ptr = NULL;
  first_rasdef_ptr = *fst_rasdef_ptr;

  /* If have rasdef details, then set flag */
  if (first_rasdef_ptr != NULL)
    have_rasdefs = TRUE;
  else
    have_rasdefs = FALSE;

  /* Open the PDB file, if found */
  if (pdb_name[0] != '\0')
    {
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
          if (gzfile_ptr == NULL)
            eof = TRUE;
        }
      else
        {
          file_ptr = fopen(pdb_name,"r");
          if (file_ptr == NULL)
            eof = TRUE;
        }
#endif
    }
  else
    eof = TRUE;

  /* If file not opened, then abort */
  if (eof == TRUE)
    {
      printf("echo \"*** ERROR. Unable to open PDB file ");
      printf("[%s] \"\n",pdb_name);
      printf("echo \" \"\n");
      exit(-1);
    }

  /* Initialise flag */
  eof = FALSE;

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
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Increment count of lines read in */
          nlines++;

          /* Check for the end of an NMR structure model */
          if (!strncmp("ENDMDL",input_line,6))
            keep_reading = FALSE;

          /* Check for the end of the current chain */
          else if (!strncmp("TER   ",input_line,6))
            {
              /* If this is the chain we're after, then mark it as having
                 finished now */
              if (have_chain == TRUE)
                chain_end = TRUE;

              /* Set last chain to NULL */
              last_chain = '\0';
            }

          /* Check for ATOM or HETATM record */
          else if (keep_reading == TRUE &&
                   (!strncmp(input_line,"ATOM",4) ||
                    !strncmp(input_line,"HETATM",6)))
            {
              /* Increment count of coord records read in */
              nread++;

              /* Get all the atom details */
              r_get_atom_details(input_line,&atom_number,atom_name,
                                 &in_chain,res_name,res_num,&x,&y,&z,
                                 &bvalue,&occupancy,element,&rec_type,
                                 wanted_chain);

              /* Rotate through 180 degrees about the x-axis, if required */
              if (flip == TRUE)
                {
                  y = - y;
                  z = - z;
                }

              /* Determine whether this is a wanted ATOM or HETATM record */
              wanted = TRUE;

              /* NMR structures in the pdb have hydrogens and pseudo-atoms, 
                 so we need to remove these */
              if (!strcmp(element," H"))
                hydrogen = TRUE;
              else
                hydrogen = FALSE;
              if (hydrogen == TRUE)
                wanted = FALSE;

              /* If this is a water, then don't want it */
              if (!strcmp(res_name,"HOH"))
                wanted = FALSE;

              /* If this is an alternate-occupancy atom, check whether
                 we already have it */
              strncpy(occup,input_line+56,3);
              occup[3] = '\0';
              if (wanted == TRUE &&
                  (input_line[16] != ' ' || strncmp(occup,"1.0",3)))
                {
                  /* Check for earlier copy of this atom */
                  have_atom
                    = check_for_duplicate(first_atom_ptr,first_lig_atom_ptr,
                                          atom_name,res_name,res_num,
                                          in_chain);

                  /* If have this atom already, then don't want this
                     second copy */
                  if (have_atom == TRUE)
                    wanted = FALSE;
                }

              /* If atom still wanted, then prepare to store */
              if (wanted == TRUE)
                {
                  /* Check whether atom belongs to a new residue */
                  if (in_chain != last_chain ||
                      strncmp(res_num,last_res_num,5) ||
                      strncmp(res_name,last_res_name,3))
                    {
                      /* Initialise C-alpha-only flag */
                      ca_only = FALSE;
                      wanted = FALSE;
                      in_ligand = FALSE;
                      have_chain = FALSE;

                      /* If reading from a given file and this is a new
                         chain, then create a new rasdef record */
                      if ((all_chains == TRUE || have_rasdefs == FALSE) &&
                          in_chain != last_chain)
                        {
                          /* Define type */
                          type = PROTEIN;
                          chain[nchain] = in_chain;
                          nchain++;

                          /* Define residue details */
                          rs_name[0] = '\0';
                          rs_num[0] = '\0';
                          colour[0] = '\0';

                          /* Create a rasdef entry */
                          rasdef_ptr
                            = create_rasdef_entry(&first_rasdef_ptr,
                                                  &last_rasdef_ptr,type,
                                                  nchain,rs_name,rs_num,
                                                  in_chain,colour);
                        }

                      /* Determine which of the objects, if any, in the
                         rasmol.dfn file this residue belongs to */
                      rasdef_ptr
                        = assign_to_object(first_rasdef_ptr,atom_name,
                                           res_name,res_num,in_chain);

                      /* If this residue is wanted, then save */
                      if (rasdef_ptr != NULL)
                        {
                          /* If protein or DNA chain, check if it belongs to
                             the requested chain */
                          if (rasdef_ptr->type == PROTEIN ||
                              rasdef_ptr->type == NUCLEIC_ACID ||
                              rasdef_ptr->type == CA_ONLY)
                            {
                              /* If all chains to be included, then residue
                                 is wanted */
                              if (all_chains == TRUE)
                                wanted = TRUE;

                              /* Otherwise, check chains in list */
                              else
                                {
                                  /* Check whether it belongs to the 
                                     required chain */
                                  for (ichain = 0; ichain < nchain; ichain++)
                                    {
                                      if (in_chain == chain[ichain])
                                        wanted = TRUE;
                                    }

                                  /* If not a wanted chain, but is a DNA
                                     chain, store in case close to the
                                     wanted chain */
                                  /* @@@ */
                                }

                              /* Set flag that we have a chain */
                              if (wanted == TRUE)
                                have_chain = TRUE;
                            }

                          /* Otherwise, is wanted as a ligand */
                          else
                            {
                              /* Set flags */
                              wanted = TRUE;
                              in_ligand = TRUE;
                            }
                        }

                      /* If this residue is wanted, then save */
                      if (wanted == TRUE)
                        {
                          /* Create a new residue record */
                          residue_ptr
                            = create_residue_record(&first_residue_ptr,
                                                    &last_residue_ptr,
                                                    res_name,res_num,in_chain);

                          /* If C-alpha only chain, then set flag */
                          if (rasdef_ptr->type == CA_ONLY)
                            ca_only = TRUE;

                          /* Increment count of residues read in */
                          nresid++;
                        }

                      /* Otherwise, residue and its atoms aren't wanted */
                      else
                        residue_ptr = NULL;

                      /* Save the current residue details */
                      strcpy(last_res_name,res_name);
                      strcpy(last_res_num,res_num);
                      last_chain = in_chain;
                    }

                  /* If have an object to which this atom's residue has
                     been assigned, then can store the atom */
                  if (residue_ptr != NULL && rasdef_ptr != NULL)
                    {
                      /* Get the atom's radius */
                      radius = get_atom_radius(element,ca_only);

                      /* Create non-ligand atom record */
                      if (in_ligand == FALSE)
                        {
                          atom_ptr
                            = create_atom_record(&first_atom_ptr,&last_atom_ptr,
                                                 residue_ptr,natoms,atom_name,
                                                 element,x,y,z,bvalue,radius);

                          /* Increment count of atoms retained */
                          if (atom_ptr != NULL)
                            natoms++;
                        }

                      /* Create ligand atom record */
                      else
                        {
                          atom_ptr
                            = create_atom_record(&first_lig_atom_ptr,
                                                 &last_lig_atom_ptr,
                                                 residue_ptr,nlig_atoms,atom_name,
                                                 element,x,y,z,bvalue,radius);

                          /* Increment count of atoms retained */
                          if (atom_ptr != NULL)
                            nlig_atoms++;

                          /* Mark residue as a ligand residue */
                          residue_ptr->type = rasdef_ptr->type;
                        }

                      /* If individual atomic radii to be used, take the value
                         read in from the B-value column */
                      if (params.indiv_atom_radii == TRUE && atom_ptr != NULL)
                        atom_ptr->radius = bvalue;
                    }
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

  /* If no atoms match the criteria, then abort */
  if (natoms == 0)
    {
      printf("*** ERROR. No atoms found which satisfy the selection ");
      printf("criteria\n");
      exit(1);
    }

  /* Show number of atoms read in */
  fprintf(logfile_ptr,"  Number of atoms read in:       %8d\n",nread);
  fprintf(logfile_ptr,"  Number stored:                 %8d\n",natoms);
  fprintf(logfile_ptr,"  Number of residues stored:     %8d\n",nresid);
  fprintf(logfile_ptr,"\n");
  fprintf(logfile_ptr,"  Number of ligand atoms:        %8d\n",nlig_atoms);
  fprintf(logfile_ptr,"\n");
  if (all_chains == TRUE)
    {
      fprintf(logfile_ptr,"  Chains: ");
      for (ichain = 0; ichain < nchain; ichain++)
        fprintf(logfile_ptr," [%c]",chain[ichain]);
      fprintf(logfile_ptr,"\n");
      fprintf(logfile_ptr,"\n");
    }

  /* Store the first atom- and residue-record pointers */
  *fst_atom_ptr = first_atom_ptr;
  *fst_lig_atom_ptr = first_lig_atom_ptr;
  *fst_residue_ptr = first_residue_ptr;
  *fst_rasdef_ptr = first_rasdef_ptr;

  /* Return number of ligand atoms and residues stored */
  *nlig = nlig_atoms;
  *nres = nresid;
  *nchn = nchain;

  /* Return the number of atoms retained */
  return(natoms);
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

calc_atom_limits  -  Calculate the atom limits

***********************************************************************/

void calc_atom_limits(struct atom *first_atom_ptr,struct limits *limit,
                      FILE *logfile_ptr,double max_gap_radius)
{
  int icoord;

  double half_box;

  struct atom *atom_ptr;

  /* Initialise variables */
  limit->first = TRUE;
  limit->max_radius = 0.0;
  limit->min_radius = 9999.9;
  limit->npoints = 0;
  for (icoord = 0; icoord < 3; icoord++)
    {
      limit->valmin[icoord] = 0.0;
      limit->valmax[icoord] = 0.0;
      limit->c_of_g[icoord] = 0.0;
    }

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop through all atoms to get the maximum and minimum coords */
  while (atom_ptr != NULL)
    {
      /* Update the maximum and minimum atom coordinates */
      get_max_min_coords(limit,atom_ptr->coords[0],atom_ptr->coords[1],
                         atom_ptr->coords[2]);

      /* If this is the largest atom radius encountered so far, then
         store */
      if (atom_ptr->radius > limit->max_radius)
        limit->max_radius = atom_ptr->radius;

      /* Repeat for smallest atom radius */
      if (atom_ptr->radius < limit->min_radius)
        limit->min_radius = atom_ptr->radius;

      /* Get pointer to the next atom in the linked list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Use maximum radius to compute extra border to be added round
     atom limits */
  half_box = limit->max_radius + 2.0 * max_gap_radius;

  /* Adjust all the limits */
  for (icoord = 0; icoord < 3; icoord++)
    {
      /* Extend the limits equally in all directions */
      limit->valmin[icoord] = limit->valmin[icoord] - half_box;
      limit->valmax[icoord] = limit->valmax[icoord] + half_box;

      /* Calculate centre of gravity */
      limit->c_of_g[icoord]
        = (limit->valmin[icoord] + limit->valmax[icoord]) / 2.0;
    }

  /* Show the coordinate limits */
  fprintf(logfile_ptr,"  Coordinate limits:\n");
  fprintf(logfile_ptr,"      x:  %8.3f %8.3f \n",limit->valmin[0],
          limit->valmax[0]);
  fprintf(logfile_ptr,"      y:  %8.3f %8.3f \n",limit->valmin[1],
          limit->valmax[1]);
  fprintf(logfile_ptr,"      z:  %8.3f %8.3f \n",limit->valmin[2],
          limit->valmax[2]);
  fprintf(logfile_ptr,"\n");
  fprintf(logfile_ptr,"  Min radius:     %8.3f\n",limit->min_radius);
  fprintf(logfile_ptr,"  Max radius:     %8.3f\n",limit->max_radius);
  fprintf(logfile_ptr,"\n");
}
/***********************************************************************

get_atom_colour  -  Determine colour of atom

***********************************************************************/

void get_atom_colour(char *element,char *colour_name)
{
  int colour, ielem;

  static char *elements[] = { " X", " C", " N", " O", " P", " S" };

  static char *colour_names[] = { "green", "dark_grey", "blue", "red",
                                    "purple", "yellow" };

  static int nelem = 6;

  /* Initialise variables */
  colour = 0;

  /* Loop over the elements to locate given atom type */
  for (ielem = 1; ielem < nelem && colour == 0; ielem++)
    {
      /* If have a match, store the colour */
      if (!strcmp(element,elements[ielem]))
        colour = ielem;
    }

  /* Save the colour */
  strcpy(colour_name,colour_names[colour]);
}
/***********************************************************************

write_raster3d_protein  -  Write out the protein's covalent bonds to the
                           output Raster3D file

***********************************************************************/

void write_raster3d_protein(FILE *outfile_ptr,
                            struct atom *first_atom_ptr,
                            double c_of_g[3],double rascale,
                            struct parameters params)
{
  char colour_name[COLNAM_LEN], other_colour_name[COLNAM_LEN];

  double radius, rgb[3], x1, x2, y1, y2, z1, z2, xm, ym, zm;
  double adj_x1, adj_x2, adj_y1, adj_y2, adj_z1, adj_z2;
  double distance_real, distance_x, distance_y, distance_z, dist_sqrd;

  struct atom *atom_ptr, *other_atom_ptr;
  struct residue *residue_ptr;

  /* Set bond radius and apply scaling factor */
  radius = rascale * BOND_RADIUS;

  /* Get pointer to first ligand atom */
  atom_ptr = first_atom_ptr;

  /* Loop through to pick out all the atoms */
  while (atom_ptr != NULL)
    {
      /* Get atom's residue */
      residue_ptr = atom_ptr->residue_ptr;

      /* Retrieve this atom's coordinates */
      x1 = atom_ptr->coords[0];
      y1 = atom_ptr->coords[1];
      z1 = atom_ptr->coords[2];

      /* Calculate atom's scaled coords */
      adj_x1 = rascale * (x1 - c_of_g[0]);
      adj_y1 = rascale * (y1 - c_of_g[1]);
      adj_z1 = rascale * (z1 - c_of_g[2]);

      /* Determine this atom's colour */
      get_atom_colour(atom_ptr->element,colour_name);

      /* Get pointer to next atom record */
      other_atom_ptr = atom_ptr->next_atom_ptr;

      /* Loop through all other atoms to find those bonded
         to the current atom */
      while (other_atom_ptr != NULL)
        {
          /* Retrieve this atom's coordinates */
          x2 = other_atom_ptr->coords[0];
          y2 = other_atom_ptr->coords[1];
          z2 = other_atom_ptr->coords[2];

          /* Calculate distance between atoms 1 & 2 */
          distance_x = (x1 - x2) * (x1 - x2);
          distance_y = (y1 - y2) * (y1 - y2);
          distance_z = (z1 - z2) * (z1 - z2);
          distance_real = distance_x + distance_y + distance_z;

          /* Calculate cut-off distance for covalent connectivity,
             based on atom types involved */
          dist_sqrd =  3.6;
          if (!strncmp(atom_ptr->element," P",2) ||
              !strncmp(other_atom_ptr->element," P",2))
            dist_sqrd =  3.8;
          if (!strncmp(atom_ptr->element,"FE",2) ||
              !strncmp(other_atom_ptr->element,"FE",2))
            dist_sqrd =  5.0;
          if (atom_ptr->element[1] == 'S' ||
              other_atom_ptr->element[1] == 'S')
            dist_sqrd =  5.0;

          /* If distance satisfies the distance criteria for covalent
             bonding, then atoms are bonded */
          if (distance_real < dist_sqrd)
            {
              /* Calculate scaled coords of this atom */
              adj_x2 = rascale * (x2 - c_of_g[0]);
              adj_y2 = rascale * (y2 - c_of_g[1]);
              adj_z2 = rascale * (z2 - c_of_g[2]);

              /* Determine this atom's colour */
              get_atom_colour(other_atom_ptr->element,other_colour_name);

              /* If colours are the same, then draw as a single
                 bond */
              if (!strcmp(colour_name,other_colour_name))
                {
                  /* Line is Raster3D object of type 3 */
                  fprintf(outfile_ptr,"3\n");

                  /* Write out first end's coords and radius */
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f %8.4f ",
                          adj_x1,adj_y1,adj_z1,radius);

                  /* Write out second end's coords and radius */
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f %8.4f\n",
                          adj_x2,adj_y2,adj_z2,radius);

                  /* Convert colour name to RGB triplet */
                  get_colour_rgb(colour_name,rgb);

                  /* Write out bond colour */
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f\n",rgb[0],
                          rgb[1],rgb[2]);
                }

              /* Otherwise, show as a split bond */
              else
                {
                  /* Calculate coords of mid-point */
                  xm = (adj_x1 + adj_x2) / 2;
                  ym = (adj_y1 + adj_y2) / 2;
                  zm = (adj_z1 + adj_z2) / 2;

                  /* Print first half of the bond-line */
                  fprintf(outfile_ptr,"3\n");
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f %8.4f ",
                          adj_x1,adj_y1,adj_z1,radius);
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f %8.4f\n",
                          xm,ym,zm,radius);

                  /* Write out its colour */
                  get_colour_rgb(colour_name,rgb);
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f\n",rgb[0],
                          rgb[1],rgb[2]);

                  /* Print second half of the bond-line */
                  fprintf(outfile_ptr,"3\n");
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f %8.4f ",
                          adj_x2,adj_y2,adj_z2,radius);
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f %8.4f\n",
                          xm,ym,zm,radius);

                  /* Write out its colour */
                  get_colour_rgb(other_colour_name,rgb);
                  fprintf(outfile_ptr,"%8.4f %8.4f %8.4f\n",rgb[0],
                          rgb[1],rgb[2]);
                }
            }

          /* Get pointer to next atom record */
          other_atom_ptr = other_atom_ptr->next_atom_ptr;
        }

      /* Get pointer to next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

mark_inrange_ligands  -  Determine which ligands lie within the atom
                         limits

***********************************************************************/

void mark_inrange_ligands(struct atom *first_lig_atom_ptr,
                          struct limits *limit)
{
  double x, y, z;

  struct atom *atom_ptr;

  /* Get pointer to the first ligand atom */
  atom_ptr = first_lig_atom_ptr;

  /* Loop through all atoms to initialise in-range flags */
  while (atom_ptr != NULL)
    {
      /* Unset the in-range flag for this atom's residue */
      atom_ptr->residue_ptr->in_range = FALSE;

      /* Get pointer to the next atom in the linked list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Get pointer to the first ligand atom again */
  atom_ptr = first_lig_atom_ptr;

  /* Loop through all atoms to set in-range flags if ligand atom is
     within the coordinate limits */
  while (atom_ptr != NULL)
    {
      /* Get the atom's coords */
      x = atom_ptr->coords[0];
      y = atom_ptr->coords[1];
      z = atom_ptr->coords[2];

      /* Check whether atom lies within the bounds */
      if (x >= limit->valmin[0] && x <= limit->valmax[0] &&
          y >= limit->valmin[1] && y <= limit->valmax[1] &&
          z >= limit->valmin[2] && z <= limit->valmax[2])
        {
          /* Set the in-range flag for this atom's residue */
          atom_ptr->residue_ptr->in_range = TRUE;
        }

      /* Get pointer to the next atom in the linked list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

compact_number  -  Convert number into shortest string for output

***********************************************************************/

int compact_number(double number,char *number_string,char *format)
{
  char temp_string[20];
  char *string_ptr;

  int ipos, len, nblanks, ppos;
  int done;

  /* Initialise variables */
  len = 0;
  for (ipos = 0; ipos < 20; ipos++)
    temp_string[ipos] = '\0';

  /* Write out the number using the given format string */
  sprintf(temp_string,format,number);
  strcpy(number_string,temp_string);

  /* Get the length of the string */
  len = strlen(number_string);

  /* If number contains a decimal point, remove any redundant zeros
     from the end */
  if ((string_ptr = strstr(number_string,".")) != NULL)
    {
      /* Get the position of the decimal point */
      ppos = string_ptr - number_string;

      /* Remove any trailing zeroes */
      done = FALSE;
      for (ipos = len - 1; ipos > ppos - 1 && done == FALSE; ipos--)
        {
          /* If this is a zero, then reduce length of string */
          if (number_string[ipos] == '0' || number_string[ipos] == '.')
            number_string[ipos] = '\0';

          /* Otherwise, set flag that we're done */
          else
            done = TRUE;
        }
    }

  /* Remove any redundant blank spaces from the front */
  done = FALSE;
  nblanks = 0;
  for (ipos = 0; ipos < ppos && done == FALSE; ipos++)
    {
      /* If this is a blank space, increment count */
      if (number_string[ipos] == ' ')
        nblanks++;

      /* Otherwise, set flag that we're done */
      else
        done = TRUE;
    }

  /* Transfer to temporary string */
  strcpy(temp_string,number_string + nblanks);

  /* Copy back */
  strcpy(number_string,temp_string);

  /* Get the string length */
  len = strlen(number_string);

  /* Return the string length */
  return(len);
}
/***********************************************************************

write_raster3d_ligatoms  -  Write out the metals and ligand atoms to the
                            output Raster3D file

***********************************************************************/

void write_raster3d_ligatoms(FILE *outfile_ptr,
                             struct atom *first_atom_ptr,
                             double c_of_g[3],double rascale,
                             struct parameters params)
{
  char colour_name[COLNAM_LEN], number_string[20];

  int icol;
  int write_sphere;

  double radius, rgb[3], x, y, z;

  struct atom *atom_ptr;
  struct rasdef *rasdef_ptr;
  struct residue *residue_ptr;

  /* Get pointer to first ligand atom */
  atom_ptr = first_atom_ptr;

  /* Loop through to pick out all the metals and plot as spheres */
  while (atom_ptr != NULL)
    {
      /* Initialise flag */
      write_sphere = FALSE;

      /* Get atom's residue */
      residue_ptr = atom_ptr->residue_ptr;

      /* Check for metal */
      if (residue_ptr != NULL && residue_ptr->type == METAL)
        {
          /* Set flag to write out as sphere */
          write_sphere = TRUE;

          /* Get colour from corresponding rasdef pointer */
          rasdef_ptr = residue_ptr->rasdef_ptr;

          /* If have a rasdef record, get colour */
          if (rasdef_ptr != NULL)
            strcpy(colour_name,rasdef_ptr->colour);

          /* Otherwise, set default colour as green */
          else
            strcpy(colour_name,"green");

          /* Set sphere radius */
          radius = params.metatom_radius;
        }

      /* Otherwise, if this is a ligand atom and we are showing them
         in spacefill or as ball-and-stick, write as sphere also */
      else if (residue_ptr != NULL && residue_ptr->type == LIGAND)
        {
          /* Set flag to write out as sphere */
          write_sphere = TRUE;

          /* If have a user-defined colour, use that */
          if (strcmp(params.ligand_colour,"CPK") &&
              strcmp(params.ligand_colour,"NONE"))
            strcpy(colour_name,params.ligand_colour);

          /* Otherwise, define sphere colour according to the atom type */
          else
            {
              get_atom_colour(atom_ptr->element,colour_name);
              if (!strcmp(colour_name,"dark_grey") &&
                  !strcmp(params.background_colour,"black"))
                strcpy(colour_name,"light_grey");
            }

          /* Set sphere radius */
          if (params.ligatom_radius < 0.0)
            radius = atom_ptr->radius;
          else if (params.ligatom_radius == 0.0)
            radius = 0.0;
          else
            radius = params.ligatom_radius;
        }

      /* If radius is zero, then don't want atom output */
      if (radius <= 0.0)
        write_sphere = FALSE;

      /* If writing as sphere, output appropriate Raster3D commands
         for colour and size */
      if (write_sphere == TRUE)
        {
          /* Sphere is object type 2 */
          fprintf(outfile_ptr,"2\n");

          /* Get atom's centre and apply scaling factor */
          x = rascale * (atom_ptr->coords[0] - c_of_g[0]);
          y = rascale * (atom_ptr->coords[1] - c_of_g[1]);
          z = rascale * (atom_ptr->coords[2] - c_of_g[2]);

          /* Apply scaling factor to radius */
          radius = rascale * radius;

          /* Write out sphere centre, radius and colour */
          compact_number(x,number_string,"%8.4f");
          fprintf(outfile_ptr,"%s",number_string);
          compact_number(y,number_string,"%8.4f");
          fprintf(outfile_ptr," %s",number_string);
          compact_number(z,number_string,"%8.4f");
          fprintf(outfile_ptr," %s",number_string);
          compact_number(radius,number_string,"%8.4f");
          fprintf(outfile_ptr," %s",number_string);

          /* Convert colour name to RGB triplet */
          get_colour_rgb(colour_name,rgb);

          /* Write out polygon colour */
          for (icol = 0; icol < 3; icol++)
            {
              /* Compact the RGB value */
              compact_number(rgb[icol],number_string,"%8.4f");

              /* Write out */
              fprintf(outfile_ptr," %s",number_string);
            }
          fprintf(outfile_ptr,"\n");
        }

      /* Get pointer to next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  char chain_id,char copy_of)
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
  chain_ptr->chain_id = chain_id;
  chain_ptr->copy_of = copy_of;
  chain_ptr->next_chain_ptr = NULL;

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return new chain pointer */
  return(chain_ptr);
}
/***********************************************************************

read_chain_equivs  -  Read in the chain equivalences from the chains.dat
                      file

***********************************************************************/

int read_chain_equivs(char *file_name,struct chain **fst_chain_ptr,
                      char *chain,int all_chains,int nchain)
{
  char input_line[LINELEN + 1];
  char chain_id, copy_of;

  int len, nchains, nlines;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  if (all_chains == TRUE)
    nchain = 0;
  nlines = 0;

  /* Initialise chain pointers */
  chain_ptr = NULL;
  first_chain_ptr = NULL;
  last_chain_ptr = NULL;

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open input file [%s]\n",file_name);
    }

  /* If OK, read through the file v*/
  else
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment count of lines read in */
          nlines++;

          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINELEN);

          /* Get the line length */
          len = strlen(input_line);

          /* Get the chain details */
          if (len > 0)
            {
              /* Get the chain */
              chain_id = input_line[0];

              /* If this is a copy of another chain, stroe that */
              if (len > 2)
                copy_of = input_line[2];
              else
                copy_of = '\0';

              /* Create a chain record */
              chain_ptr
                = create_chain_record(&first_chain_ptr,&last_chain_ptr,
                                      chain_id,copy_of);

              /* If need all chains, store this one in the array and
                 increment chain count */
              if (all_chains == TRUE && nchain < MAX_CHAINS)
                {
                  /* Store this chain */
                  chain[nchain] = chain_id;

                  /* Increment chain count */
                  nchain++;
                }
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Return number of chains stored in chain array */
  return(nchain);
}
/***********************************************************************

check_chain_match  -  Check whether the two given chains are equivalent

***********************************************************************/

int check_chain_match(struct chain *first_chain_ptr,char chain1,
                      char chain2)
{
  int chain_match;

  struct chain *chain_ptr;

  /* Initialise variables */
  chain_match = FALSE;

  /* If have an exact match, then need go no further */
  if (chain1 == chain2)
    chain_match = TRUE;

  /* Check for blank chain id */
  else if (chain1 == ' ' && chain2 == 'A')
    chain_match = TRUE;

  /* Otherwise, search through the chains list to see if these are
     equivalent chains */
  else
    {
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
    }

  /* Return answer */
  return(chain_match);
}
/***********************************************************************

get_conservation_scores  -  Pick up all the residue conservation scores
                            from the conserv.dat file in this PDB file's
                            PDBsum directory

***********************************************************************/

void get_conservation_scores(char *filename,struct residue *first_residue_ptr,
                             struct chain *first_chain_ptr,
                             FILE *logfile_ptr)
{
  char input_line[LINELEN + 1], number_string[20];
  char chain, last_chain, res_name[4], res_num[6];

  int nlines, nmatched;
  int chain_match, found, have_file;

  double conserv;

  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Reading conserv.dat file ...\n");
  have_file = TRUE;
  nlines = 0;
  nmatched = 0;

  /* Open input residue conservation scores file */
  if ((file_ptr = fopen(filename,"r")) ==  NULL)
    {
      fprintf(logfile_ptr,"\n*** Unable to open file [%s]\n",filename);
      have_file = FALSE;
    }

  /* Loop while reading in records from the PDB file */
  while (have_file == TRUE && fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment count of lines read in */
      nlines++;

      /* Get the residue details and conservation score */
      strncpy(res_name,input_line,3);
      res_name[3] = '\0';
      chain = input_line[4];
      strncpy(res_num,input_line+5,5);
      res_num[5] = '\0';
      strncpy(number_string,input_line+10,7);
      number_string[7] = '\0';
      conserv = atof(number_string);

      /* Initialise flags */
      chain_match = FALSE;
      found = FALSE;
      last_chain = '\0';

      /* Get pointer to the first of the stored residues */
      residue_ptr = first_residue_ptr;

      /* Loop through all the residues until current one is found */
      while (residue_ptr != NULL)
        {
          /* Initialise flags */
          found = FALSE;

          /* If this is a new chain, check whether it matches the one
             just read in from the file */
          if (residue_ptr->chain != last_chain)
            {
              /* Search through chain-equivalences list to see if these
                 two chains match */
              chain_match
                = check_chain_match(first_chain_ptr,chain,
                                    residue_ptr->chain);

              /* Store the chain */
              last_chain = residue_ptr->chain;
            }

          /* If residue details match, then set flag */
          if (!strcmp(residue_ptr->res_name,res_name) &&
              !strcmp(residue_ptr->res_num,res_num) &&
              chain_match == TRUE)
            found = TRUE;

          /* If residue number and chain match, check whether we have
             a special residue */
          else if (!strcmp(residue_ptr->res_num,res_num) &&
                   chain_match == TRUE)
            {
              /* Check for MSE */
              if (!strcmp(residue_ptr->res_name,"MSE") &&
                  !strcmp(res_name,"MET"))
                found = TRUE;
            }

          /* If residue found, then store conservation score */
          if (found == TRUE)
            {
              /* Store conversation score */
              residue_ptr->conservation_score = conserv;

              /* Increment count of matched residues */
              nmatched++;
            }

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
        }
    }

  /* Show number of matched residues */
  fprintf(logfile_ptr,"    Number of matched residues = %d\n",nmatched);
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

get_consurf_scores  -  Pick up all the residue conservation scores
                       from this PDB file's ConSurf data file

***********************************************************************/

void get_consurf_scores(char *pdb_code,int all_chains,char *chain_list,
                        int nchain,char *consurf_dir,
                        struct chain *first_chain_ptr,
                        struct residue *first_residue_ptr,
                        FILE *logfile_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, last_chain, number_string[20], res_name[4], res_num[6];
  char chain2[2], middle[3], tmp[LINELEN + 1];

  char *token[MAX_TOKEN], *string_ptr;

  int ichain, len, nlines, nmatched, ntokens;
  int chain_match, found, have_file, reading;

  double conserv;

  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Loop over the required chains */
  for (ichain = 0; ichain < nchain; ichain++)
    {
      /* Initialise variables */
      fprintf(logfile_ptr,"Reading ConSurf data for chain [%c] ...\n",
              chain_list[ichain]);
      have_file = TRUE;
      nlines = 0;
      nmatched = 0;
      reading = FALSE;

      /* Get middle part of PDB code */
      strncpy(middle,pdb_code + 1,2);
      middle[2] = '\0';

      /* Form name of the corresponding ConSurf file */
      strcpy(file_name,consurf_dir);
      strcat(file_name,DIR_SEP);
      strcat(file_name,middle);
      strcat(file_name,DIR_SEP);
      strcat(file_name,pdb_code);
      if (chain_list[ichain] == ' ')
        strcat(file_name,"_");
      else
        {
          chain2[0] = chain_list[ichain];
          chain2[1] = '\0';
          strcat(file_name,chain2);
        }
      strcat(file_name,".gradesPE");

      /* Open input residue conservation scores file */
      if ((file_ptr = fopen(file_name,"r")) ==  NULL)
        {
          fprintf(logfile_ptr,"\n*** Unable to open file [%s]\n",
                  file_name);
          have_file = FALSE;
        }
      else
        fprintf(logfile_ptr,"Reading %s ...\n",file_name);

      /* Loop while reading in records from the PDB file */
      while (have_file == TRUE &&
             fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment count of lines read in */
          nlines++;

          /* If reading in the residue data, process this line */
          if (reading == TRUE)
            {
              /* Unpack the line of tab-separated fields */
              ntokens = get_tokens(input_line,token,MAX_TOKEN,'\t');
              if (ntokens > MAX_TOKEN)
                ntokens = MAX_TOKEN;

              /* If have sufficient data, then store details */
              if (ntokens > 6)
                {
                  /* Pre-process the residue field, removing blanks
                     from start */
                  strcpy(tmp,token[2]);
                  remove_leading_blanks(tmp);

                  /* Get the residue details */
                  strncpy(res_name,tmp,3);
                  res_name[3] = '\0';
                  len = strlen(tmp);
                  chain = tmp[len - 1];
                  if (chain == ':')
                    chain = ' ';
                  strcpy(number_string,tmp + 3);
                  string_ptr = strstr(number_string,":");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  strcpy(res_num,number_string);

                  /* Format the residue number, adding the relevant
                     numbers of blanks either side of the number */
                  format_res_num(res_num);

                  /* Get the conservation score */
                  conserv = atof(token[5]) / 10.0 + 0.01;

                  /* Initialise flag */
                  chain_match = FALSE;
                  last_chain = '\0';
                  found = FALSE;

                  /* Get pointer to the first of the stored residues */
                  residue_ptr = first_residue_ptr;

                  /* Loop through all the residues until current one
                     is found */
                  while (residue_ptr != NULL)
                    {
                      /* Initialise flag */
                      found = FALSE;

                      /* If this is a new chain, check whether it matches
                         the one just read in from the file */
                      if (residue_ptr->chain != last_chain)
                        {
                          /* Search through chain-equivalences list to see
                             if these two chains match */
                          chain_match
                            = check_chain_match(first_chain_ptr,chain,
                                                residue_ptr->chain);

                          /* Store the chain */
                          last_chain = residue_ptr->chain;
                        }

                      /* If residue details match, then set flag */
                      if (chain_match == TRUE &&
                          !strcmp(residue_ptr->res_name,res_name) &&
                          !strcmp(residue_ptr->res_num,res_num))
                        found = TRUE;

                      /* If residue number and chain match, check whether
                         we have a special residue */
                      else if (chain_match == TRUE &&
                               !strcmp(residue_ptr->res_num,res_num))
                        {
                          /* Check for MSE */
                          if (!strcmp(residue_ptr->res_name,"MSE") &&
                              !strcmp(res_name,"MET"))
                            found = TRUE;
                        }

                      /* If residue found, then store conservation score */
                      if (found == TRUE)
                        {
                          /* Store conversation score */
                          residue_ptr->conservation_score = conserv;

                          /* Increment count of matched residues */
                          nmatched++;
                        }

                      /* Get pointer to the next residue */
                      residue_ptr = residue_ptr->next_residue_ptr;
                    }
                }

              /* Free up memory taken up by the tokens */
              free_tokens(token,ntokens);
            }

          /* Otherwise check for column titles line signifying the start
             of the data */
          else if (!strncmp(input_line," POS",4))
            reading = TRUE;
        }

      /* If found file, close it */
      if (have_file == TRUE)
        fclose(file_ptr);
    }

  /* Show number of matched residues */
  fprintf(logfile_ptr,"    Number of matched residues = %d\n",nmatched);
}
/***********************************************************************

calc_grid_size  -  Calculate the size of the grid-array onto which to
                   map the sphere clusters

***********************************************************************/

double calc_grid_size(double valmin[3],double valmax[3],int npoints[3],
                      double max_radius,double grid_sep,FILE *logfile_ptr)
{
  int icoord, loop, lower[3], upper[3];
  int done;

  long array_size, narray;

  double diff, grid_spacing;

  /* Define the grid spacing */
  /*  grid_spacing = sqrt(3.0) * max_radius; */
  grid_spacing = grid_sep;

  /* Add an extra border to the max and min values */
  for (icoord = 0; icoord < 3; icoord++)
    {
      valmin[icoord] = valmin[icoord] - grid_spacing;
      valmax[icoord] = valmax[icoord] + grid_spacing;
    }

  /* Initialise flag */
  done = FALSE;

  /* Loop until have a suitable grid */
  for (loop = 0; loop < MAX_LOOPS && done == FALSE; loop++)
    {
      /* Calculate lowest grid value for each coordinate */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Get nearest integer equivalent of min coord */
          lower[icoord] = (int) (valmin[icoord] / grid_spacing);
          if (valmin[icoord] < 0.0)
            lower[icoord] = lower[icoord] - 1;
          diff = fabs((double) lower[icoord] - valmin[icoord] / grid_spacing);
          if (diff < 0.0001)
            lower[icoord] = lower[icoord] - 1;

          /* Repeat for maximum coordinate */
          upper[icoord] = (int) (valmax[icoord] / grid_spacing);
          if (valmax[icoord] < 0.0)
            upper[icoord] = upper[icoord] - 1;
          diff = fabs((double) upper[icoord] - valmax[icoord] / grid_spacing);
          if (diff < 0.0001)
            upper[icoord] = upper[icoord] + 1;
          upper[icoord] = upper[icoord] + 1;
          
          /* Calculate number of grid points required for this coordinate */
          npoints[icoord] = upper[icoord] - lower[icoord] + 1;

          /* Adjust the minimum and maximum values to correspond to the
             grid limits */
          valmin[icoord] = lower[icoord] * grid_spacing;
          valmax[icoord] = upper[icoord] * grid_spacing;
        }

      /* Calculate array size required and check that it isn't too large */
      array_size = npoints[0] * npoints[1] * npoints[2];
      fprintf(logfile_ptr,"Loop %d. Array size = %d x %d x %d = %ld\n",
             loop,npoints[0],npoints[1],npoints[2],array_size);
      fflush(logfile_ptr);
      if (array_size > MAX_ARRAY_SIZE)
        {
          /* Adjust the grid spacing to reduce the array size */
          grid_spacing = grid_spacing * (double) array_size
            / (double) (MAX_ARRAY_SIZE + loop + 1);
        }

      /* Otherwise, set flag that we're done */
      else
        done = TRUE;
    }

  /* Show size of grid array created */
  narray = npoints[0] * npoints[1] * npoints[2];
  fprintf(logfile_ptr,"Grid-size = %d x %d x %d = %ld\n",npoints[0],
          npoints[1],npoints[2],narray);
  fprintf(logfile_ptr,"\n");

  /* Show ajusted coordinate limits */
  fprintf(logfile_ptr,"  Adjusted coordinate limits:\n");
  fprintf(logfile_ptr,"      x:  %8.3f %8.3f \n",valmin[0],valmax[0]);
  fprintf(logfile_ptr,"      y:  %8.3f %8.3f \n",valmin[1],valmax[1]);
  fprintf(logfile_ptr,"      z:  %8.3f %8.3f \n",valmin[2],valmax[2]);
  fprintf(logfile_ptr,"\n");
  fflush(logfile_ptr);

  /* Calculate array size required and check that it isn't too large */
  array_size = npoints[0] * npoints[1] * npoints[2];
  if (array_size > MAX_ARRAY_SIZE)
    {
      fprintf(logfile_ptr,
              "*** ERROR. Maximum array limit exceeded. Program aborted\n");
      exit(1);
    }

  /* Return the calculated grid spacing */
  return(grid_spacing);
}
/***********************************************************************

read_delphi_file  -  Read in the grid holding the electrostatic potential

***********************************************************************/

int read_delphi_file(char *file_name,float **potential,double *grid_sep,
                     double *valmin,double *min_energy,double *max_energy,
                     struct limits limit,FILE *logfile_ptr)
{
#define HEADER_SIZE      106

  char header[HEADER_SIZE];

  int i, ipoint, nbytes, ngrid, nhalf_grid, npoints, npts;

  float origin[3], scale;
  float *grid;

  double c_of_g[3], grid_separation;

  FILE *file_ptr;

  /* Initialise variables */
  *min_energy = 0.0;
  *max_energy = 0.0;

  /* Calculate centre of gravity */
  for (i = 0; i < 3; i++)
    c_of_g[i] = (limit.valmin[i] + limit.valmax[i]) / 2;

  /* Show name of file */
  fprintf(logfile_ptr,"Reading Delphi file: %s ...\n",file_name);

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"rb")) !=  NULL)
    {
      /* Read the title */
      fread(header,sizeof(char)*HEADER_SIZE,1,file_ptr);

      /* Read in the total number of grid-points */
      fread(&nbytes,sizeof(int),1,file_ptr);
      npoints = nbytes / sizeof(float);

      /* Create the array to hold the Delphi electrostatic potential */
      grid = (float *) malloc(npoints * sizeof(float));
      if (grid == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for Delphi ");
          printf("potential\n");
          exit(1);
        }

      /* Read in the grid */
      npts = fread(grid,sizeof(float),npoints,file_ptr);
      if (npts != npoints)
        {
          printf("*** ERROR. Fewer Delphi grid points read in than");
          printf("expected\n");
          printf("***        Read: %d.  Expected: %d\n",npts,npoints);
          exit(0);
        }

      /* Read in number of bytes again as a check */
      fread(&nbytes,sizeof(int),1,file_ptr);
      if (nbytes != (int) (npoints * sizeof(float)))
        {
          printf("*** ERROR. Error in reading Delphi array\n");
          printf("***        Bytes: %d.  Expected: %d\n",nbytes,
                 npoints * (int) sizeof(float));
          exit(0);
        }

      /* Read in header text */
      fread(header,sizeof(char),28,file_ptr);

      /* Read in scale factor */
      fread(&scale,sizeof(float),1,file_ptr);

      /* Read in coords of grid origin */
      fread(origin,sizeof(float),3,file_ptr);

      /* Get Delphi grid size */
      fread(&ngrid,sizeof(int),1,file_ptr);
      if (ngrid * ngrid * ngrid != npoints)
        {
          printf("*** ERROR. Mismatch in Delphi grid size\n");
          printf("***        Grid: %d x %d x %d = %d.  Have: %d\n",
                 ngrid,ngrid,ngrid,ngrid * ngrid * ngrid,npoints);
          exit(0);
        }

      /* Close the file */
      fclose(file_ptr);

      /* Loop through the potential energy grid to get maximum and
         minimum energy values */
      for (ipoint = 0; ipoint < npoints; ipoint++)
        {
          /* Check for highest and lowest values so far */
          if (grid[ipoint] < *min_energy)
            *min_energy = grid[ipoint];
          if (grid[ipoint] > *max_energy)
            *max_energy = grid[ipoint];
        }
    }

  /* Can't find file, so show error message */
  else
    {
      printf("*** ERROR. File not found: %s\n",file_name);
      exit(0);
    }

  /* Calculate grid separation in Angstroms */
  grid_separation = (float) (1.0 / scale);

  /* Calculate number of grid-spacings from centre of the grid to
     its corner */
  nhalf_grid = (ngrid - 1) / 2;

  /* Calculate the adjusted coordinates of the corner of the grid */
  for (i = 0; i < 3; i++)
    {
      /* Calculate raw coords */
      valmin[i] = origin[i] - nhalf_grid * grid_separation;

      /* Adjust to C of G of whole molecule */
      valmin[i] = valmin[i] - c_of_g[i];
    }

  /* Return the grid separation */
  *grid_sep = grid_separation;

  /* Show statistics for data read in */
  fprintf(logfile_ptr,"Grid size: %d x %d x %d\n",ngrid,ngrid,ngrid);
  fprintf(logfile_ptr,"Origin:      %8.3f %8.3f %8.3f\n",origin[0],
          origin[1],origin[2]);
  fprintf(logfile_ptr,"Min. energy: %8.3f\n",*min_energy);
  fprintf(logfile_ptr,"Max. energy: %8.3f\n",*max_energy);
  fprintf(logfile_ptr,"Array range:\n");
  fprintf(logfile_ptr,"  x:   %8.3f %8.3f\n",valmin[0],
          valmin[0] + (ngrid -1) * grid_separation);
  fprintf(logfile_ptr,"  y:   %8.3f %8.3f\n",valmin[1],
          valmin[1] + (ngrid -1) * grid_separation);
  fprintf(logfile_ptr,"  z:   %8.3f %8.3f\n",valmin[2],
          valmin[2] + (ngrid -1) * grid_separation);

  /* Return pointer to electrostatic potential array */
  *potential = grid;

  /* Return grid size */
  return(ngrid);
}
/***********************************************************************

create_grid_array  -  Create the grid array

***********************************************************************/

void create_grid_array(int **grid,int npoints[3])
{
  int ipoint, ngrid_points;
  int *grid_map;

  /* Calculate overall size of grid map */
  ngrid_points = npoints[0] * npoints[1] * npoints[2];

  /* Create the 3D grid map as a 1D array */
  grid_map = (int *) malloc(ngrid_points * sizeof(int));
  if (grid_map == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for grid map\n");
      exit(1);
    }

  /* Initialise grid array */
  for (ipoint = 0; ipoint < ngrid_points; ipoint++)
    grid_map[ipoint] = EMPTY;

  /* Return the grid-map pointer */
  *grid = grid_map;
}
/***********************************************************************

create_atom_pointers_array  -  Create a grid of atom pointers, each
                               corresponding to the atom that overlaps
                               that grid pointer

***********************************************************************/

void create_atom_pointers_array(int **grd_atom_map,int npoints[3])
{
  int ipoint, ngrid_points;
  int *grid_atom_map;

  /* Calculate overall size of grid map */
  ngrid_points = npoints[0] * npoints[1] * npoints[2];

  /* Create the 3D grid map as a 1D array */
  grid_atom_map = (int *) malloc(ngrid_points * sizeof(int));
  if (grid_atom_map == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for atom pointers ");
      printf("array\n");
      exit(1);
    }

  /* Initialise grid array */
  for (ipoint = 0; ipoint < ngrid_points; ipoint++)
    grid_atom_map[ipoint] = -1;

  /* Return the pointer to the start of the array */
  *grd_atom_map = grid_atom_map;
}
/***********************************************************************

create_index_arrays  -  Create the arrays holding the x-, y- and z-
                        sort-indices

***********************************************************************/

void create_index_arrays(struct atom *first_atom_ptr,
                         struct atom ***atom_index_ptr,int natoms)
{
  int i, iatom, ipos;
  struct atom **atom_ind_ptr;
  struct atom *atom_ptr;

  /* Create the three index arrays as one large array */
  atom_ind_ptr = (struct atom **) malloc(4 * natoms * sizeof(struct atom *));

  /* Check that sufficient memory was allocated */
  if (atom_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sort-index array\n");
      exit(1);
    }

  /* Loop over the three coordinates */
  for (i = 0; i < 3; i++)
    {
      /* Get pointer to the first atom */
      iatom = 0;
      atom_ptr = first_atom_ptr;

      /* Loop through all atoms to get current index array to point to
         each atom in turn */
      while (atom_ptr != NULL)
        {
          /* Get array position */
          ipos = i * natoms + iatom;

          /* Initialise the index array */
          atom_ind_ptr[ipos] = atom_ptr;

          /* Get pointer to the next atom in the linked list */
          atom_ptr = atom_ptr->next_atom_ptr;
          iatom++;
        }
    }

  /* Return the array pointer */
  *atom_index_ptr = atom_ind_ptr;
}
/***********************************************************************

shell_sort  -  Shell-sort which sorts the double values into ascending
               order

***********************************************************************/

void shell_sort(int *index,double *value,int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  double aln2i = 1.4426950, tiny = 1.e-5;
  double calc;
  double dble, dble1;

  double val, other_val;

  /*--Initialise values */
  n = nindex;
  dble = n;
  dble1 = log(dble);
  calc = dble1;
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

              /* Get the two vals */
              val = value[index[indi]];
              other_val = value[index[indl]];

              /* If in wrong order, the swap sort-pointers */
              if (val > other_val)
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

sort_index_arrays  -  Sort atoms according to their x-, y- and z-
                       coordinates

***********************************************************************/

void sort_index_arrays(struct atom **atom_index_ptr,int natoms)
{
  int i, iatom, iorder, ipos, jorder, store_pos;
  int done;
  int *index;

  double *coordinate;

  struct atom *atom_ptr, *store_atom_ptr;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(natoms * sizeof(int));
  coordinate = (double *) malloc(natoms * sizeof(double));

  /* Loop over the three coordinates */
  for (i = 0; i < 3; i++)
    {

      /* Loop over all the atoms to transfer the current set of
         coordinates into the temporary sort arrays */
      for (iatom = 0; iatom < natoms; iatom++)
        {
          /* Get array position */
          ipos = i * natoms + iatom;

          /* Get atom pointer */
          atom_ptr = atom_index_ptr[ipos];

          /* Transfer the appropriate coordinate value */
          coordinate[iatom] = atom_ptr->coords[i];
          index[iatom] = iatom;
        }

      /* Sort the coordinate values */
      shell_sort(index,coordinate,natoms);

      /* Loop through the sorted list to place the atom-pointers into
         the sorted pointers array for this coordinate */
      for (iorder = 0; iorder < natoms; iorder++)
        {
          /* Save the atom pointer at this order position as it will
             be replaced by the atom pointer we wish to place here */
          store_pos = i * natoms + iorder;
          store_atom_ptr = atom_index_ptr[store_pos];

          /* Get the atom number of the atom that we want to place at
             this position */
          iatom = index[iorder];

          /* Get this atom's position in the original array */
          ipos = i * natoms + iatom;

          /* If atom not in its required position, then need to shift it */
          if (ipos != store_pos)
            {
              /* Get the atom's pointer */
              atom_ptr = atom_index_ptr[ipos];

              /* Transfer it to the required position */
              atom_index_ptr[store_pos] = atom_ptr;

              /* Put the saved pointer into the position from which the
                 other pointer had been taken from */
              atom_index_ptr[ipos] = store_atom_ptr;
              done = FALSE;

              /* Find the index-pointer to the displaced atom and adjust
                 it to point to the atom's new position */
              for (jorder = iorder + 1; jorder < natoms && done == FALSE;
                   jorder++)
                {
                  /* If this is the one we're after, perform the
                     amendment */
                  if (index[jorder] == iorder)
                    {
                      /* Amend the index entry and set flag */
                      index[jorder] = iatom;
                      done = TRUE;
                    }
                }
            }
        }
    }
}
/***********************************************************************

mask_grid_points  -  Use all the atom radii to mask the grid-points in
                     the grid-array that are within the protein

***********************************************************************/

int mask_grid_points(int *grid,int *grid_atom_map,int npoints[3],
                     struct atom **atom_index_ptr,int natoms,double valmin[3],
                     double valmax[3],double grid_sep,
                     FILE *logfile_ptr,struct parameters params)
{
  int i, iatom, icoord, ipoint, j, jatom, k;
  int end_point[3], start_point[3];
  int ncalc, nmasked;

  double calc, dist_sqrd, dx_sqrd, dy_sqrd, dz_sqrd, radius_sqrd, x, y, z;
  double cell_volume, curr_dist, dist, volume;

  struct atom *atom_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Masking grid points for %d atoms ...\n",natoms);
  ncalc = 0;
  nmasked = 0;

  /* Loop over all the atoms, using the x-sorted atom list */
  for (iatom = 0; iatom < natoms; iatom++)
    {
      /* Get atom pointer */
      atom_ptr = atom_index_ptr[iatom];

      /* Get the square of this atom's radius */
      radius_sqrd = atom_ptr->radius * atom_ptr->radius;

      /* Loop over the coordinates to determine start and end grid-points
         in each direction to cover current atom */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Get start-point */
          calc  = atom_ptr->coords[icoord] - (atom_ptr->radius + 2 * grid_sep);
          start_point[icoord] = (int) ((calc - valmin[icoord]) / grid_sep);
          if (calc - valmin[icoord] < 0.0)
            start_point[icoord] = start_point[icoord] - 1;

          /* Adjust if out of bounds */
          if (start_point[icoord] < 0)
            start_point[icoord] = 0;
          if (start_point[icoord] > npoints[icoord] - 1)
            start_point[icoord] = npoints[icoord] - 1;

          /* Repeat for end-point */
          calc  = atom_ptr->coords[icoord] + (atom_ptr->radius + 2 * grid_sep);
          end_point[icoord] = (int) ((calc - valmin[icoord]) / grid_sep) + 1;
          if (calc - valmin[icoord] < 0.0)
            end_point[icoord] = end_point[icoord] - 1;

          /* Adjust if out of bounds */
          if (end_point[icoord] < 0)
            end_point[icoord] = 0;
          if (end_point[icoord] > npoints[icoord] - 1)
            end_point[icoord] = npoints[icoord] - 1;
        }

      /* Loop over all the x grid-points */
      for (i = start_point[0]; i < end_point[0] + 1; i++)
        {
          /* Get this grid-point's x-coordinate */
          x = valmin[0] + i * grid_sep;
          dx_sqrd = (x - atom_ptr->coords[0])
            * (x - atom_ptr->coords[0]);
          ncalc++;

          /* Loop over all the y grid-points */
          for (j = start_point[1]; j < end_point[1] + 1; j++)
            {
              /* Get this grid-point's y-coordinate */
              y = valmin[1] + j * grid_sep;
              dy_sqrd = (y - atom_ptr->coords[1])
                * (y - atom_ptr->coords[1]);
              ncalc++;

              /* Loop over all the z grid-points */
              for (k = start_point[2]; k < end_point[2] + 1; k++)
                {
                  /* Get this grid-point's z-coordinate */
                  z = valmin[2] + k * grid_sep;
                  dz_sqrd = (z - atom_ptr->coords[2])
                    * (z - atom_ptr->coords[2]);
                  ncalc++;

                  /* Get the distance between the grid-point and the
                     atom centre */
                  dist_sqrd = dx_sqrd + dy_sqrd + dz_sqrd;

                  /* Get the array location of this grid-point */
                  ipoint = i * npoints[1] * npoints[2]
                    + j * npoints[2] + k;

                  /* If grid-point is within the radius, then set */
                  if (dist_sqrd <= radius_sqrd)
                    {
                      /* If not already masked, increment count of
                         masked points */
                      if (grid[ipoint] != MASKED)
                        nmasked++;

                      /* Update the grid-point */
                      grid[ipoint] = MASKED;

                      /* If need to store the atom pointer, then do so */
                      if (params.store_atom_pointers == TRUE)
                        {
                          /* Store reference to this atom at this
                             grid-point */
                          grid_atom_map[ipoint] = 1000 * iatom;
                        }
                    }

                  /* Otherwise, only interested if storing atom pointers */
                  else if (params.store_atom_pointers == TRUE)
                    {
                      /* If grid-point has no atom associated with it, then
                         set current distance to a large value */
                      if (grid_atom_map[ipoint] == -1)
                        curr_dist = 10.0;

                      /* Otherwise, unpack atom number and its distance */
                      else
                        {
                          /* Get the distance from the previous atom
                             associated with this grid-point */
                          jatom = (int) (grid_atom_map[ipoint] / 1000);
                          curr_dist = (grid_atom_map[ipoint]
                                               - 1000 * jatom) / 100;
                        }

                      /* If previous distance was not zero, then determine
                         distance from current atom's surface to this
                         grid-point */
                      if (curr_dist > 0.001)
                        {
                          /* Get distance to centre */
                          dist = sqrt(dist_sqrd);

                          /* Get distance to atom's surface */
                          dist = fabs(dist - atom_ptr->radius);

                          /* If this distance is lower than previous one
                             store reference to current atom */
                          if (dist < curr_dist && dist < 9.99)
                            {
                              /* Convert atom reference and distance to a
                                 single integer value */
                              grid_atom_map[ipoint] = 1000 * iatom
                                + (int) (dist * 100 + 0.5);
                            }
                        }
                    }
                }
            }
        }
    }

  /* Calculate cell and protein volumes */
  cell_volume = grid_sep * grid_sep * grid_sep;
  volume = nmasked * cell_volume;

  /* Show the protein volume */
  fprintf(logfile_ptr,"Number of masked grid-points:  %8d\n",nmasked);
  fprintf(logfile_ptr,"Corresponding protein volume:  %8.3f\n",volume);

  /* Return number of masked points */
  return(nmasked);
}
/***********************************************************************

calculate_ijk_index  -  Calculate index of relative grid-points
                        corresponding to larger and larger sphere sizes

***********************************************************************/

int calculate_ijk_index(int relative_ijk[MAX_SPHERE_POINTS][3],
                        int index_ijk[MAX_SPHERE_POINTS],
                        double distance_ijk[MAX_SPHERE_POINTS],
                        double grid_sep,double max_radius,FILE *logfile_ptr)
{
  int nijk;
  int i, j, k, n, nall, napprox, npts;

  double dist_sqrd, volume, x, x_sqrd, y, y_sqrd, z, z_sqrd;
  double max_rad_sqrd;

  /* Initialise variables */
  nijk = 0;
  max_rad_sqrd = max_radius * max_radius;

  /* Calculate extent of array required to encompass the maximum sphere */
  n = (int) ((max_radius / grid_sep + 1) + 0.5);
  npts = 2 * n + 1;
  nall = npts * npts * npts;
  fprintf(logfile_ptr,"Maximum sphere size allowed = %.3f\n",max_radius);
  fprintf(logfile_ptr,"Extent in grid-points %3d x %3d x %3d = %7d\n",
          npts, npts, npts,nall);

  /* Loop over all the grid-points, calculating their distance from the
     origin */
  for (i = -n; i < n + 1; i++)
    {
      /* Calculate relative x-coordinate */
      x = i * grid_sep;
      x_sqrd = x * x;

      for (j = -n; j < n + 1; j++)
        {
          /* Calculate relative y-coordinate */
          y = j * grid_sep;
          y_sqrd = y * y;

          for (k = -n; k < n + 1; k++)
            {
              /* Calculate relative z-coordinate */
              z = k * grid_sep;
              z_sqrd = z * z;

              /* Calculate the distance of this point from the origin */
              dist_sqrd = x_sqrd + y_sqrd + z_sqrd;

              /* If this is within the maximum sphere, then store this
                 point */
              if (dist_sqrd <= max_rad_sqrd)
                {
                  /* Check that the maximum number of points haven't been
                     exceeded */
                  if (nijk >= MAX_SPHERE_POINTS - 1)
                    {
                      /* Too few points in array */
                      fprintf(logfile_ptr,
                              "*** ERROR. Number of points, MAXPTS,");
                      fprintf(logfile_ptr,
                              " insufficient to encompass maximum ");
                      fprintf(logfile_ptr,
                              "sphere, MAXRAD: %8.3f\n",max_radius);

                      /* Calculate approximate number required */
                      volume = (4.0 / 3.0) * 4.0 * PI
                        * max_radius * max_radius * max_radius;
                      napprox
                        = (int) (volume / (grid_sep * grid_sep * grid_sep));
                      fprintf(logfile_ptr,"***        Need approx. ");
                      fprintf(logfile_ptr,"%d points\n",napprox);
                      exit(1);
                    }

                  /* Store this point */
                  relative_ijk[nijk][0] = i;
                  relative_ijk[nijk][1] = j;
                  relative_ijk[nijk][2] = k;
                  index_ijk[nijk] = nijk;
                  distance_ijk[nijk] = sqrt(dist_sqrd);

                  /* Increment number of points stored */
                  nijk++;
                }
            }
        }
    }

  /* Write out number of points generated */
  fprintf(logfile_ptr,"Number of reference grid-points generated = %d\n",
          nijk);

  /* Sort the grid-points in ascending distance from the origin */
  shell_sort(index_ijk,distance_ijk,nijk);

  /* Return number of grid-points in relative arrays */
  return(nijk);
}
/***********************************************************************

binary_atom_search  -  Search for the given coordinate value in the given
                       x-, y- or z-coordinate index and return the first
                       array position whose value is higher than that
                       supplied. If the value supplied is larger than any
                       value in the array, then -1 is returned.

***********************************************************************/

int binary_atom_search(struct atom **atom_index_ptr,int natoms,
                       double search_value,int icoord)
{
  int bottom_pos, ipoint, ipos, new_pos, top_pos;
  int found;

  struct atom *atom_ptr;

  /* Make starting point half-way into the index */
  ipos = natoms / 2;
  bottom_pos = 0;
  top_pos = natoms - 1;

  /* Initialise variables */
  found = FALSE;

  /* Loop until required array position found */
  while (found == FALSE)
    {
      /* Get array position */
      ipoint = icoord * natoms + ipos;

      /* Get the associated atom pointer */
      atom_ptr = atom_index_ptr[ipoint];

      /* If array value is equal to or higher than the search value,
         check whether it is the value we want */
      if (atom_ptr->coords[icoord] >= search_value)
        {
          /* If this is the first atom in the sorted list, then
             this is the one we want */
          if (ipos == 0)
            found = TRUE;

          /* Otherwise, check whether the previous atom's coordinate
             is below the search value */
          else
            {
              /* Get the pointer to the previous atom */
              atom_ptr = atom_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after */
              if (atom_ptr->coords[icoord] < search_value)
                found = TRUE;

              /* Otherwise, we need to search earlier in the list */
              else
                {
                  /* Save the current position as the uppermost position
                     to search on */
                  top_pos = ipos;

                  /* Set pointer to halfway between the last search
                     low-search and here */
                  ipos = (ipos + bottom_pos) / 2;
                }
            }

        }

      /* If current atom's coords are below the search value, have to
         look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_pos = ipos;

          /* If have reached the top of the list, then our value is
             higher than any value in the list, so return -1 */
          if (ipos >= natoms - 1)
            {
              /* Set ipos to -1 and flag as done */
              ipos = -1;
              found = TRUE;
            }

          /* Otherwise, can search further up the list */
          else
            {
              /* Set pointer to halfway between the last search high-search
                 and here */
              new_pos = (ipos + top_pos) / 2;

              /* Adjust if going back to the same point */
              if (new_pos == ipos)
                new_pos++;
              ipos = new_pos;
            }
        }
    }

  /* Return the required position in the index array */
  return(ipos);
}
/***********************************************************************

get_atom_start_end  -  Get start- and end-atoms given the start and
                       end-coords to span

***********************************************************************/

int get_atom_start_end(struct atom **atom_index_ptr,int natoms,
                       int icoord,double value,double max_radius,
                       int start_point[3],int end_point[3])
{
  int nsearch;

  double from_coord, to_coord;

  /* Initialise variables */
  start_point[icoord] = 0;
  end_point[icoord] = -1;

  /* Determine range of coords required */
  from_coord = value - max_radius;
  to_coord = value + max_radius;

  /* Perform a binary search to locate the starting atom for the
     lower limit of the range */
  start_point[icoord] = binary_atom_search(atom_index_ptr,natoms,from_coord,
                                           icoord);

  /* Repeat for upper limit of the range */
  end_point[icoord] = binary_atom_search(atom_index_ptr,natoms,to_coord,
                                         icoord);

  /* If both ends of the range are off either end, then there are no
     atoms within the range */
  if (start_point[icoord] == -1 || end_point[icoord] == 0)
    {
      /* Set start and end points such that no atoms are processed */
      start_point[icoord] = 0;
      end_point[icoord] = -1;
    }

  /* Check for upper-range value being off the end */
  else if (end_point[icoord] == -1)
    end_point[icoord] = natoms - 1;

  /* Otherwise, reduce end position by 1 */
  else
    end_point[icoord]--;

  /* Calculate number of atoms in list */
  if (start_point[icoord] > end_point[icoord])
    nsearch = 0;
  else
    nsearch = end_point[icoord] - start_point[icoord] + 1;

  /* Return number of atoms in list */
  return(nsearch);
}
/***********************************************************************

transfer_atom_list  -  Transfer atoms from the x-sorted list to the
                       all-purpose sort-stack

***********************************************************************/

int transfer_atom_list(struct atom **atom_index_ptr,int natoms,
                       int start_point,int end_point,
                       int ref_atom_number)
{
  int ipoint, jpoint, nstored;

  struct atom *atom_ptr;

  /* Initialise variables */
  nstored = 0;

  /* Loop over the atoms that have been selected on their x-coords */
  for (ipoint = start_point; ipoint < end_point + 1; ipoint++)
    {
      /* Get this atom's pointer */
      atom_ptr = atom_index_ptr[ipoint];

      /* Transfer only if atom-number is higher than the atom being
         compared against */
      if (atom_ptr->atom_number > ref_atom_number)
        {
          /* Get the location of where to store the pointer */
          jpoint = 3 * natoms + nstored;

          /* Store the pointer */
          atom_index_ptr[jpoint] = atom_ptr;

          /* Initialise the atom's checked field */
          atom_ptr->checked = 0;

          /* Increment count of atoms stored */
          nstored++;
        }
    }

  /* Return number of atoms stored in filtered list */
  return(nstored);
}
/***********************************************************************

filter_atom_coord_range  -  Filter out the unwanted atoms from the current
                            list

***********************************************************************/

int filter_atom_coord_range(int icoord,struct atom **atom_index_ptr,
                            int natoms,int nlist,double from_coord,
                            double to_coord,double middle_from_coord,
                            double middle_to_coord)
{
  int ipoint, jpoint, nstored;
  int wanted;

  struct atom *atom_ptr;

  /* Initialise variables */
  nstored = 0;

  /* Loop over tha atoms that have been selected so far */
  for (ipoint = 3 * natoms; ipoint < 3 * natoms + nlist; ipoint++)
    {
      /* Get this atom's pointer */
      atom_ptr = atom_index_ptr[ipoint];

      /* Determine whether the atom lies in range */
      wanted = FALSE;
      if (atom_ptr->coords[icoord] >= from_coord &&
          atom_ptr->coords[icoord] <= to_coord)
        wanted = TRUE;

      /* If atom too close on this coordinate, increment counter */
      if (atom_ptr->coords[icoord] > middle_from_coord &&
          atom_ptr->coords[icoord] < middle_to_coord)
        {
          /* Increment counter */
          atom_ptr->checked++;

          /* If counter has reached 3, then atom too close on all
             3 coordinates, and so can comfortably be rejected */
          if (atom_ptr->checked == 3)
            wanted = FALSE;
        }

      /* If atom's coordinates lie within the required range, copy
         the atom's pointer to the filtered list */
      if (wanted == TRUE)
        {
          /* Get the location of where to store the pointer */
          jpoint = 3 * natoms + nstored;

          /* Store the pointer */
          atom_index_ptr[jpoint] = atom_ptr;

          /* Increment count of atoms stored */
          nstored++;
        }
    }

  /* Return number of atoms stored in filtered list */
  return(nstored);
}
/***********************************************************************

get_atom_neighbour_list  -  Get list of atoms within the gap-span
                            distance of the current atom

***********************************************************************/

int get_atom_neighbour_list(struct atom *atom_ptr,
                            struct atom **atom_index_ptr,int natoms,
                            double reach_dist,double too_close_dist)
{
  int icoord, nlist;
  int end_point[3], start_point[3];

  double from_coord, to_coord;
  double middle_from_coord, middle_to_coord;

  /* Initialise variables */
  nlist = 0;

  /* Get start- and end-atoms in the sorted x-coord index */
  nlist = get_atom_start_end(atom_index_ptr,natoms,0,atom_ptr->coords[0],
                             reach_dist,start_point,end_point);

  /* If have any atoms, continue */
  if (nlist > 0)
    {
      /* Transfer atoms from the x-sorted list to the all-purpose
         sort-stack */
      nlist = transfer_atom_list(atom_index_ptr,natoms,start_point[0],
                                 end_point[0],atom_ptr->atom_number);

      /* Loop over the coordinates to further filter the sort-list */
      for (icoord = 0; icoord < 3 && nlist > 0; icoord++)
        {
          /* Determine range of coordinates required */
          from_coord = atom_ptr->coords[icoord] - reach_dist;
          to_coord = atom_ptr->coords[icoord] + reach_dist;

          /* Get the range in the middle that isn't wanted (ie atoms too
             close to the current one) */
          middle_from_coord = atom_ptr->coords[icoord] - too_close_dist;
          middle_to_coord = atom_ptr->coords[icoord] + too_close_dist;

          /* Filter out the unwanted atoms from the current list */
          nlist = filter_atom_coord_range(icoord,atom_index_ptr,natoms,
                                          nlist,from_coord,to_coord,
                                          middle_from_coord,
                                          middle_to_coord);
        }
    }

  /* Return number of atom pointers in the list */
  return(nlist);
}
/***********************************************************************

get_nearest_grid_point  -  Calculate the nearest grid-point to the given
                           x-, y- and z-coords

***********************************************************************/

int get_nearest_grid_point(int npoints[3],double valmin[3],double valmax[3],
                           double grid_sep,double coords[3],
                           int nearest_point[3])
{
  int icoord, ipoint;

  /* Find the nearest grid-point to the given coords */
  for (icoord = 0; icoord < 3; icoord++)
    {
      /* Calculate the nearest grid-point in this direction */
      nearest_point[icoord]
        = (int) ((coords[icoord] - valmin[icoord]) / grid_sep + 0.5);

      /* Adjust if out of bounds */
      if (nearest_point[icoord] < 0)
        nearest_point[icoord] = 0;
      if (nearest_point[icoord] > npoints[icoord] - 1)
        nearest_point[icoord] = npoints[icoord] - 1;
    }

  /* Get the array location of this grid-point */
  ipoint = nearest_point[0] * npoints[1] * npoints[2]
    + nearest_point[1] * npoints[2] + nearest_point[2];

  /* Return the array location */
  return(ipoint);
}
/***********************************************************************

get_sphere_positions  -  Determine how many gap-spheres, and their
                         relative spacing, can be placed within the
                         given distance

***********************************************************************/

int get_sphere_positions(double distance,struct parameters params,
                         double min_gap_diameter,double sphere_sep,
                         double *step_size)
{
  int ndiam, nspheres;

  /* Initialise variables */
  nspheres = 0;
  *step_size = 0.0;

  /* If can fit gap spheres between these two atoms, then proceed */
  if (distance >= min_gap_diameter && distance <= params.max_gap_span)
    {
      /* If separation is less than the maximum gap diameter,
         then can only fit a single sphere between the two atoms */
      if (distance <= min_gap_diameter)
        {
          /* Define number of spheres as 1 */
          nspheres = 1;
        }

      /* Otherwise, determine how many spheres of maximum size we can
         place between these two atoms */
      else
        {
          /* Calculate number of sphere diameters between the two
             atoms */
          ndiam = (int) ((distance - min_gap_diameter) / sphere_sep);

          /* Calculate number of loops to go through */
          nspheres = ndiam + 2;

          /* Calculate step-size */
          *step_size = (distance - min_gap_diameter) / (nspheres - 1);
        }
    }

  /* Return the number of spheres to be placed */
  return(nspheres);
}
/***********************************************************************

mark_sphere_locations  -  Update the grid-points along the line joining
                          the current pair of atoms

***********************************************************************/

void mark_sphere_locations(int *grid,int npoints[3],double valmin[3],
                           double valmax[3],double grid_sep,
                           double coords[3],double other_coords[3],
                           double distance,int nspheres,double step_size,
                           double radius,double min_gap_radius)
{
  int icoord, ipoint, isphere, nearest_point[3];

  double fraction, pcoords[3], sphere_pos, start_pos;

  /* Initialise variables */

  /* Calculate the start position for the first sphere */
  start_pos = radius + min_gap_radius;
  if (nspheres == 1)
    start_pos = distance / 2.0;

  /* Loop over all the possible sphere locations */
  for (isphere = 0; isphere < nspheres; isphere++)
    {
      /* Determine position of the centre of this sphere */
      sphere_pos = start_pos + isphere * step_size;

      /* Calculate how far along the line between the two atoms the
         sphere-centre lies */
      fraction = sphere_pos / distance;

      /* Calculate coords of the centre of the sphere */
      for (icoord = 0; icoord < 3; icoord++)
        {
          pcoords[icoord] = coords[icoord]
            + fraction * (other_coords[icoord] - coords[icoord]);
        }

      /* Find the closest grid-point to this coordinate */
      ipoint = get_nearest_grid_point(npoints,valmin,valmax,grid_sep,
                                      pcoords,nearest_point);

      /* Update the grid-point if not already masked */
      if (grid[ipoint] == EMPTY)
        grid[ipoint] = SPHERE_LOC;
    }
}
/***********************************************************************

mark_potential_gap_points  -  Mark all grid-points that lie in potential
                              gap-sphere positions (ie between a pair of
                              atoms the right distance apart)

***********************************************************************/

void mark_potential_gap_points(int *grid,int npoints[3],double valmin[3],
                               double valmax[3],double grid_sep,
                               struct atom **atom_index_ptr,int natoms,
                               struct atom *first_atom_ptr,double min_radius,
                               double max_radius,struct parameters params,
                               int relative_ijk[MAX_SPHERE_POINTS][3],
                               int index_ijk[MAX_SPHERE_POINTS],
                               double distance_ijk[MAX_SPHERE_POINTS],
                               FILE *logfile_ptr)
{
  int ilist, ipoint, nlist, nspheres;

  double distance, dist_sqrd, reach_dist, too_close_dist;
  double min_gap_diameter, sphere_sep, step_size;
  double surf_sep;

  struct atom *atom_ptr, *other_atom_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Marking potential gap points ...\n");
  min_gap_diameter = 2.0 * params.min_gap_radius;
  sphere_sep = min_gap_diameter;
  if (sphere_sep > grid_sep)
    sphere_sep = grid_sep;

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop over all the atoms */
  while (atom_ptr != NULL)
    {
      /* Calculate maximum distance any other atom can be from current
         atom such that maximum spanning distance lies between them */
      reach_dist = params.max_gap_span + atom_ptr->radius + max_radius;

      /* Calculate minimum distance between the two atom which will still
         allow a minimum-sized gap sphere to be placed between them */
      too_close_dist = min_gap_diameter + atom_ptr->radius
        + min_radius;

      /* Convert the minimum distance to a minimum separation in all
         three coords that will ensure minimum distance */
      too_close_dist = too_close_dist / SQRT_THREE;

      /* Get list of atoms within the gap-span distance of the current
         atom */
      nlist = get_atom_neighbour_list(atom_ptr,atom_index_ptr,natoms,
                                      reach_dist,too_close_dist);

      /* Loop over all the sufficiently-close atoms to find locations
         where a gap-sphere might be placed */
      for (ilist = 0; ilist < nlist; ilist++)
        {
          /* Get this atom's pointer */
          ipoint = 3 * natoms + ilist;
          other_atom_ptr = atom_index_ptr[ipoint];

          /* Calculate the distance between these two atoms */
          dist_sqrd
            = (atom_ptr->coords[0] - other_atom_ptr->coords[0])
            * (atom_ptr->coords[0] - other_atom_ptr->coords[0])
            + (atom_ptr->coords[1] - other_atom_ptr->coords[1])
            * (atom_ptr->coords[1] - other_atom_ptr->coords[1])
            + (atom_ptr->coords[2] - other_atom_ptr->coords[2])
            * (atom_ptr->coords[2] - other_atom_ptr->coords[2]);
          distance = sqrt(dist_sqrd);

          /* Calculate surface-surface distance */
          surf_sep = distance - atom_ptr->radius - other_atom_ptr->radius;

          /* Calculate fractional points along line joining the two
             atoms between which we can start and end the sphere-
             filling process */
          nspheres = get_sphere_positions(surf_sep,params,min_gap_diameter,
                                          sphere_sep,&step_size);

          /* If have spheres that can be placed between these two atoms,
             update the grid-points along the line joining them */
          if (nspheres > 0)
            {
              mark_sphere_locations(grid,npoints,valmin,valmax,grid_sep,
                                    atom_ptr->coords,other_atom_ptr->coords,
                                    distance,nspheres,step_size,
                                    atom_ptr->radius,params.min_gap_radius);
            }
        }

      /* Get pointer to the next atom in the linked list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

grow_sphere  -  Generate largest sphere possible at the given grid-point

***********************************************************************/

double grow_sphere(int *grid,int npoints[3],int current_i,int current_j,
                  int current_k,
                  int relative_ijk[MAX_SPHERE_POINTS][3],
                  int index_ijk[MAX_SPHERE_POINTS],
                  double distance_ijk[MAX_SPHERE_POINTS],int nijk,
                  double grid_sep,int *fitted_max,
                  int grid_ref[MAX_SPHERE_POINTS],int *nrf)
{
  int i, igrid, ipoint, ipos, j, jgrid, k, nref;
  int done, fitted, end_sphere;

  double radius, sphere_radius;

  /* Initialise variables */
  done = FALSE;
  fitted = FALSE;
  *fitted_max = FALSE;
  nref = 0;
  sphere_radius = 0.0;
  radius = 0.0;

  /* Loop until largest sphere found, or maximum sphere reached */
  for (ipos = 0; ipos < nijk && done == FALSE; ipos++)
    {
      /* Get this grid-point's corresponding entry number */
      igrid = index_ijk[ipos];

      /* Use the relative grid-position to determine the actual
         grid position */
      i = current_i + relative_ijk[igrid][0];
      j = current_j + relative_ijk[igrid][1];
      k = current_k + relative_ijk[igrid][2];

      /* Save the radius at this grid-position */
      radius = distance_ijk[igrid];

      /* Check that none of these is off the end of the grid */
      if (i >= 0 && i < npoints[0] &&
          j >= 0 && j < npoints[1] &&
          k >= 0 && k < npoints[2])
        {
          /* Get the array-location of this grid-point */
          ipoint = i * npoints[1] * npoints[2]
            + j * npoints[2] + k;

          /* If have hit an occupied grid-point, need to stop */
          if (grid[ipoint] == MASKED)
            {
              /* Set flag to indicate that have hit an atom */
              done = TRUE;
            }

          /* Otherwise, check whether this point represents the last
             point for a given sphere radius */
          else
            {
              /* Set flag that a sphere has been fitted */
              fitted = TRUE;

              /* Initialise flag */
              end_sphere = FALSE;

              /* If this is the last point, then have come to the end of
                 the maximum sphere */
              if (ipos == nijk - 1)
                {
                  /* Set flags  */
                  end_sphere = TRUE;
                  *fitted_max = TRUE;
                }

              /* Otherwise, check sphere-radius of the following point */
              else
                {
                  /* Get the entry for the next point */
                  jgrid = index_ijk[ipos + 1];

                  /* If it corresponds to a greater sphere, then have closed
                     another shell */
                  if (distance_ijk[jgrid] > distance_ijk[igrid] + 0.000001)
                    end_sphere = TRUE;
                }

              /* If have closed off another shell, store shell radius */
              if (end_sphere == TRUE)
                {
                  /* Store current shell radius as the sphere radius */
                  sphere_radius = distance_ijk[igrid];
                }

              /* Store reference to this unmasked grid-point */
              grid_ref[nref] = ipoint;
              nref++;
            }
        }

      /* If off end, then don't want to go any further and regard any
         spheres fitted so far as not truly fitted */
      else
        {
          done = TRUE;
          fitted = FALSE;
        }
    }

  /* If didn't hit the maximum sphere, then return an average radius
     between that of the completed shell and the next uncompleted one */
  if (fitted == TRUE && *fitted_max == FALSE)
    sphere_radius = (sphere_radius + radius) / 2.0;

  /* Return the number of grid reference points stored */
  *nrf = nref;

  /* Return radius of grown sphere */
  return(sphere_radius);
}
/***********************************************************************

check_sphere  -  Check whether the current gap-sphere interpenetrates
                 any neighbouring atoms

***********************************************************************/

double check_sphere(int *grid_atom_map,int npoints[3],
                   double x,double y,double z,double sphere_radius,
                   int *fitted_max,
                   struct atom **atom_index_ptr,int natoms,
                   int grid_ref[MAX_SPHERE_POINTS],int nref,
                   double min_radius)
{
  int iatom, ipoint, iref, jref;
  int cut, found;

  double dist, dist2, limit, xa, ya, za;

  struct atom *atom_ptr;

  /* Initialise variables */
  cut = FALSE;

  /* Loop through all the stored grid-points to pick up the atoms
     that lie close to them */
  for (iref = 0; iref < nref; iref++)
    {
      /* Get the current grid-point */
      ipoint = grid_ref[iref];

      /* If there is no atom associated with this gird-point, then
         ignore */
      if (grid_atom_map[ipoint] == -1)
        grid_ref[iref] = -1;

      /* Otherwise, retrieve the atom reference */
      else
        {
          /* Get atom reference */
          iatom = (int) (grid_atom_map[ipoint] / 1000);
          grid_ref[iref] = iatom;

          /* Initialise flag */
          found = FALSE;

          /* Loop through all the references stored so far to check
             whether we already have this one */
          for (jref = 0; jref < iref && found == FALSE; jref++)
            {
              /* If already have this one, then ignore */
              if (grid_ref[jref] == iatom)
                {
                  /* Mark as found */
                  found = TRUE;
                  grid_ref[iref] = -1;
                }
            }
        }
    }

  /* Loop a second time to calculate the distance between each atom
     and the current gap sphere */
  for (iref = 0; iref < nref && sphere_radius > 0.0; iref++)
    {
      /* If have an atom reference, get the corresponding atom */
      if (grid_ref[iref] != -1)
        {
          /* Get atom reference */
          iatom = grid_ref[iref];

          /* Get the pointer to this atom */
          atom_ptr = atom_index_ptr[iatom];

          /* Get this atom's coordinates */
          xa = atom_ptr->coords[0];
          ya = atom_ptr->coords[1];
          za = atom_ptr->coords[2];

          /* Calculate the centre-centre distance from the atom
             to sphere */
          dist2 = (x - xa) * (x - xa)
            + (y - ya) * (y - ya)
            + (z - za) * (z - za);

          /* Get length of combined radii */
          limit = sphere_radius + atom_ptr->radius;

          /* If distance is smaller than this limit, then atom
             must interpenetrate the gap-sphere */
          if (dist2 < limit * limit)
            {
              /* Set flag that sphere has been cut */
              cut = TRUE;

              /* Calculate the surface-to-surface distance */
              dist = sqrt(dist2);

              /* If sphere within the atom, then abort altogether */
              if (dist < atom_ptr->radius)
                sphere_radius = 0.0;

              /* Otherwise, adjust sphere's radius until it just
                 touches this atom */
              else
                {
                  /* Adjust radius */
                  sphere_radius = dist - atom_ptr->radius;

                  /* If radius has gone below the minimum, then
                     abandon this sphere */
                  if (sphere_radius < min_radius)
                    sphere_radius = 0.0;
                }
            }

        }      
    }

  /* If had a maximum-sized sphere here and it was cut, then mark
     it as no longer being of maximum size */
  if (*fitted_max == TRUE && cut == TRUE)
    *fitted_max = FALSE;
  
  /* Return radius of grown sphere */
  return(sphere_radius);
}
/***********************************************************************

place_spheres  -  Place spheres at the gap-sphere grid-points and
                  determine whether they match the acceptance criteria

***********************************************************************/

void place_spheres(int *grid,int *grid_atom_map,int npoints[3],
                   double valmin[3],double valmax[3],double grid_sep,
                   struct atom **atom_index_ptr,int natoms,
                   double max_radius,struct parameters params,
                   int relative_ijk[MAX_SPHERE_POINTS][3],
                   int index_ijk[MAX_SPHERE_POINTS],
                   double distance_ijk[MAX_SPHERE_POINTS],
                   FILE *logfile_ptr)
{
  int i, ipoint, j, k;
  int ncut, ngrid, nijk, nlost, nmax, nmax_lost, nref;
  int grid_ref[MAX_SPHERE_POINTS];
  int last_x_store, last_y_store, store_radius;
  int fitted_max, save_max;

  double reach_dist, x, y, z;
  double save_radius, sphere_radius;

  /* Initialise variables */
  fprintf(logfile_ptr,"Placing spheres ...\n");
  last_x_store = -1;
  last_y_store = -1;
  ncut = 0;
  ngrid = 0;
  nlost = 0;
  nmax = 0;
  nmax_lost = 0;
  reach_dist = params.max_gap_span + max_radius;

  /* Calculate index of relative grid-points corresponding to larger
     and larger sphere sizes */
  nijk = calculate_ijk_index(relative_ijk,index_ijk,distance_ijk,
                             grid_sep,params.max_gap_radius,logfile_ptr);

  /* Loop over all the x grid-points */
  for (i = 0; i < npoints[0]; i++)
    {
      /* Get this grid-point's x-coordinate */
      x = valmin[0] + i * grid_sep;

      /* Loop over all the y grid-points */
      for (j = 0; j < npoints[1]; j++)
        {
          /* Get this grid-point's y-coordinate */
          y = valmin[1] + j * grid_sep;

          /* Loop over all the z grid-points */
          for (k = 0; k < npoints[2]; k++)
            {
              /* Get this grid-point's z-coordinate */
              z = valmin[2] + k * grid_sep;

              /* Get the array location of this grid-point */
              ipoint = i * npoints[1] * npoints[2]
                + j * npoints[2] + k;

              /* If grid-point is a potential gap-sphere location, then
                 check what size of gap-sphere can be placed here
                 without overlapping any atoms */
              if (grid[ipoint] == SPHERE_LOC)
                {
                  /* Generate largest sphere possible at this point */
                  sphere_radius
                    = grow_sphere(grid,npoints,i,j,k,relative_ijk,
                                  index_ijk,distance_ijk,nijk,grid_sep,
                                  &fitted_max,grid_ref,&nref);

                  /* Check whether any of the neighbouring atoms
                     interpenetrate current gap-sphere */
                  if (sphere_radius > params.min_gap_radius &&
                      params.check_all_spheres == TRUE)
                    {
                      /* Save current settings */
                      save_radius = sphere_radius;
                      save_max = fitted_max;

                      /* Perform check on sphere */
                      sphere_radius
                        = check_sphere(grid_atom_map,npoints,x,y,z,
                                       sphere_radius,
                                       &fitted_max,atom_index_ptr,natoms,
                                       grid_ref,nref,
                                       params.min_gap_radius);

                      /* See if sphere was reduced at all */
                      if (save_max == TRUE && fitted_max == FALSE)
                        nmax_lost++;
                      else if (sphere_radius < save_radius)
                        {
                          /* Increment count of spheres cut */
                          ncut++;

                          /* If sphere dropped below minimum, increment
                             count of spheres lost */
                          if (sphere_radius <= params.min_gap_radius)
                            nlost++;
                        }
                    }

                  /* If maximum-sized sphere can be fitted here, mark
                     the grid-point accordingly */
                  if (fitted_max == TRUE)
                    {
                      grid[ipoint] = MAX_SPHERE;
                      nmax++;
                    }

                  /* Otherwise, mark to note that a smaller sphere
                     can be placed here */
                  else if (sphere_radius > params.min_gap_radius)
                    {
                      /* Convert sphere radius into a scaled integer
                         value for storing in the array */
                      store_radius = (int) (GRID_SCALE_FACTOR * sphere_radius);
                      if (store_radius > SMALLER_SPHERE)
                        {
                          grid[ipoint] = store_radius;
                          ngrid++;
                        }
                    }
                }
            }
        }
    }

  /* Show numbers of spheres fitted */
  fprintf(logfile_ptr,"\n");
  fprintf(logfile_ptr,"Number of max-size spheres fitted   = %d\n",nmax);
  fprintf(logfile_ptr,"Number of under-sized spheres       = %d\n",ngrid);

  /* If any were cut, show how many */
  if (params.check_all_spheres == TRUE)
    {
      fprintf(logfile_ptr,"\n");
      fprintf(logfile_ptr,"Number of max-size spheres cut      = %d\n",
              nmax_lost);
      fprintf(logfile_ptr,"Number of spheres cut               = %d\n",
              ncut);
      fprintf(logfile_ptr,"Number of spheres lost              = %d\n",
              nlost);
    }
}
/***********************************************************************

create_sphere_record  -  Create and initialise a new sphere record

***********************************************************************/

struct sphere *create_sphere_record(struct sphere **fst_sphere_ptr,
                                    struct sphere **lst_sphere_ptr,
                                    double radius,double x,double y,
                                    double z)
{
  struct sphere *sphere_ptr, *first_sphere_ptr, *last_sphere_ptr;

  /* Initialise sphere pointers */
  sphere_ptr = NULL;
  first_sphere_ptr = *fst_sphere_ptr;
  last_sphere_ptr = *lst_sphere_ptr;

  /* Allocate memory for structure to hold sphere info */
  sphere_ptr = (struct sphere *) malloc(sizeof(struct sphere));
  if (sphere_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct sphere\n");
      exit (1);
    }

  /* If this is the very first sphere, then store pointer */
  if (first_sphere_ptr == NULL)
    first_sphere_ptr = sphere_ptr;

  /* Add link from previous sphere to the current one */
  if (last_sphere_ptr != NULL)
    last_sphere_ptr->next_sphere_ptr = sphere_ptr;
  last_sphere_ptr = sphere_ptr;

  /* Store the sphere's radius and type */
  sphere_ptr->radius = radius;

  /* Store sphere's coordinates */
  sphere_ptr->coords[0] = x;
  sphere_ptr->coords[1] = y;
  sphere_ptr->coords[2] = z;

  /* Initialise the other elements of the structure */
  sphere_ptr->cluster = 0;
  sphere_ptr->next_sphere_ptr = NULL;

  /* Return current pointers */
  *fst_sphere_ptr = first_sphere_ptr;
  *lst_sphere_ptr = last_sphere_ptr;

  /* Return new sphere pointer */
  return(sphere_ptr);
}
/***********************************************************************

generate_gap_spheres  -  Generate the gap spheres from the valid
                         grid-points

***********************************************************************/

int generate_gap_spheres(int *grid,int npoints[3],double valmin[3],
                         double valmax[3],double grid_sep,
                         struct sphere **fst_over_sphere_ptr,
                         struct sphere **lst_over_sphere_ptr,
                         struct sphere **fst_under_sphere_ptr,
                         double min_gap_radius,double max_gap_radius,
                         int *nu,int *no,FILE *logfile_ptr)
{
  int nover, nunder;
  int i, ipoint, j, k, nspheres;

  double sphere_radius;
  double x, y, z;

  struct sphere *first_over_sphere_ptr, *last_over_sphere_ptr;
  struct sphere *first_under_sphere_ptr, *last_under_sphere_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Generating gap spheres ...\n");
  nspheres = 0;
  nover = 0;
  nunder = 0;

  /* Initialise pointers */
  first_over_sphere_ptr = NULL;
  first_under_sphere_ptr = NULL;
  last_over_sphere_ptr = NULL;
  last_under_sphere_ptr = NULL;

  /* Loop over all the x grid-points */
  for (i = 0; i < npoints[0]; i++)
    {
      /* Get this grid-point's x-coordinate */
      x = valmin[0] + i * grid_sep;

      /* Loop over all the y grid-points */
      for (j = 0; j < npoints[1]; j++)
        {
          /* Get this grid-point's y-coordinate */
          y = valmin[1] + j * grid_sep;

          /* Loop over all the z grid-points */
          for (k = 0; k < npoints[2]; k++)
            {
              /* Get this grid-point's z-coordinate */
              z = valmin[2] + k * grid_sep;

              /* Get the array location of this grid-point */
              ipoint = i * npoints[1] * npoints[2]
                + j * npoints[2] + k;

              /* If grid-point holds the maximum-sized gap-sphere,
                 create a sphere record */
              if (grid[ipoint] == MAX_SPHERE)
                {
                  /* Get the sphere's radius */
                  sphere_radius = max_gap_radius;

                  /* Increment count */
                  nover++;

                  /* Create an over-sized gap-sphere */
                  create_sphere_record(&first_over_sphere_ptr,
                                       &last_over_sphere_ptr,
                                       sphere_radius,x,y,z);

                  /* Increment count of spheres generated */
                  nspheres++;
                }

              /* Otherwise, if grid-point holds an under-sized sphere,
                 retrieve its radius and store */
              else if (grid[ipoint] > SMALLER_SPHERE)
                {
                  /* Extract the sphere's radius */
                  sphere_radius = grid[ipoint] / GRID_SCALE_FACTOR;

                  /* Increment count */
                  nunder++;


                  /* Create an under-sized gap-sphere */
                  create_sphere_record(&first_under_sphere_ptr,
                                       &last_under_sphere_ptr,
                                       sphere_radius,x,y,z);

                  /* Increment count of spheres generated */
                  nspheres++;
                }
            }
        }
    }

  /* Return current pointers */
  *fst_over_sphere_ptr = first_over_sphere_ptr;
  *lst_over_sphere_ptr = last_over_sphere_ptr;
  *fst_under_sphere_ptr = first_under_sphere_ptr;

  /* Return counts of spheres */
  *nu = nunder;
  *no = nover;

  /* Return number of spheres generated */
  fprintf(logfile_ptr,"  Number of spheres over size limit:  %8d\n",
          nover);
  fprintf(logfile_ptr,"  Number of spheres under size limit: %8d\n",
          nunder);
  fprintf(logfile_ptr,"  Number of gap-spheres generated:    %8d\n",
          nspheres);
  fprintf(logfile_ptr,"\n");

  /* If have no spheres, then print warning message */
  if (nspheres == 0 || nover == 0)
    {
      fprintf(logfile_ptr,"*** Warning. No over-sized spheres generated!\n");
      fprintf(logfile_ptr,"\n");
    }

  return(nspheres);
}
/***********************************************************************

join_clusters  -  Assign higher-numbered cluster to the overlapping
                  lower-numbered cluster

***********************************************************************/

void join_clusters(struct sphere *first_sphere_ptr,int cluster,
                   int other_cluster)
{
  struct sphere *sphere_ptr;

  /* Get the first of the spheres */
  sphere_ptr = first_sphere_ptr;

  /* Loop through all the spheres */
  while (sphere_ptr != NULL)
    {
      /* If this sphere's cluster-number needs to be altered, then do so */
      if (sphere_ptr->cluster == other_cluster)
        sphere_ptr->cluster = cluster;

      /* Get pointer to the next sphere */
      sphere_ptr = sphere_ptr->next_sphere_ptr;
    }
}
/***********************************************************************

cluster_spheresx  -  Cluster all the over-size spheres into separate
                     non-interacting clusters and assign any associated
                     under-spheres to those clusters ... Old method

***********************************************************************/

int cluster_spheresx(struct sphere *first_over_sphere_ptr,
                     FILE *logfile_ptr,double grid_sep)
{
  int higher_cluster, ncluster;
  int new_cluster;

  double d_sqrd, x, x1, y, y1, z, z1;
  double overlap;
  double radius, other_radius;

  struct sphere *sphere_ptr, *other_sphere_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Clustering spheres ...\n");
  ncluster = 0;

  /* Get the first of the over-sized spheres */
  sphere_ptr = first_over_sphere_ptr;

  /* Loop through all the over-spheres to get spheres they overlap */
  while (sphere_ptr != NULL)
    {
      /* Initialise flag */
      new_cluster = FALSE;

      /* Get sphere's radius */
      radius = sphere_ptr->radius;

      /* Get its coordinates */
      x = sphere_ptr->coords[0];
      y = sphere_ptr->coords[1];
      z = sphere_ptr->coords[2];

      /* If sphere has no cluster number, then assign next one to it */
      if (sphere_ptr->cluster == 0)
        {
          /* Increment cluster number and assign to this sphere */
          ncluster++;
          sphere_ptr->cluster = ncluster;

          /* Set flag to indicate this is currently an orphan sphere */
          new_cluster = TRUE;
        }

      /* Get pointer to the next sphere */
      other_sphere_ptr = sphere_ptr->next_sphere_ptr;

      /* Loop through all the other spheres to find which penetrate
         this one */
      while (other_sphere_ptr != NULL)
        {
          /* Get other sphere's radius */
          other_radius = other_sphere_ptr->radius;

          /* Define the overlap distance as the sum of the two radii */
          /*      overlap = radius + other_radius; */
          /*      overlap = 1.01 * grid_sep; */
          overlap = (radius + other_radius) / 2;

          /* Get its coordinates */
          x1 = other_sphere_ptr->coords[0];
          y1 = other_sphere_ptr->coords[1];
          z1 = other_sphere_ptr->coords[2];

          /* Calculate distance squared */
          d_sqrd = (x - x1) * (x - x1) + (y - y1) * (y - y1)
            + (z - z1) * (z - z1);

          /* If the two spheres overlap then both belong to the same
             cluster */
          if (d_sqrd < overlap * overlap)
            {
              /* Easy case is where second sphere hasn't yet been assigned
                 to a cluster */
              if (other_sphere_ptr->cluster == 0)
                {
                  /* Assign it to our current cluster number */
                  other_sphere_ptr->cluster = sphere_ptr->cluster;

                  /* Reset flag */
                  new_cluster = FALSE;
                }

              /* Otherwise, sphere has already been assigned to a cluster */
              else if (sphere_ptr->cluster != other_sphere_ptr->cluster)
                {
                  /* If our sphere was a new sphere, then reassign its
                     cluster number */
                  if (new_cluster == TRUE)
                    {
                      /* Assign current sphere to other sphere's cluster */
                      sphere_ptr->cluster = other_sphere_ptr->cluster;

                      /* Reset flag */
                      new_cluster = FALSE;

                      /* Decrement cluster number */
                      ncluster--;
                    }

                  /* Otherwise, need to join two sphere-clusters into one */
                  else
                    {
                      /* Join the two clusters, using the smaller of the
                         two numbers */
                      if (sphere_ptr->cluster > other_sphere_ptr->cluster)
                        {
                          /* Store the higher cluster-number */
                          higher_cluster = sphere_ptr->cluster;

                          /* Join the clusters */
                          join_clusters(first_over_sphere_ptr,
                                        other_sphere_ptr->cluster,
                                        sphere_ptr->cluster);
                        }
                      else
                        {
                          /* Store the higher cluster-number */
                          higher_cluster = other_sphere_ptr->cluster;

                          /* Join the clusters */
                          join_clusters(first_over_sphere_ptr,
                                        sphere_ptr->cluster,
                                        other_sphere_ptr->cluster);
                        }

                      /* If we have renamed the current cluster number,
                         then reset */
                      if (higher_cluster == ncluster)
                        ncluster--;
                    }
                }
            }

          /* Get pointer to the next sphere */
          other_sphere_ptr = other_sphere_ptr->next_sphere_ptr;
        }

      /* Get pointer to the next sphere */
      sphere_ptr = sphere_ptr->next_sphere_ptr;
    }

  /* Show number of sphere clusters */
  fprintf(logfile_ptr,"  Number of clusters:                 %8d\n",
          ncluster);
  fprintf(logfile_ptr,"\n");

  /* Return number of clusters */
  return(ncluster);
}
/***********************************************************************

transfer_spheresx  -  Assign to clusters any of the under-size spheres
                      that penetrate any of the over-size spheres and
                      transfer to the over-size sphere list

***********************************************************************/

int transfer_spheresx(struct sphere *first_over_sphere_ptr,
                      struct sphere *last_over_sphere_ptr,
                      struct sphere *first_under_sphere_ptr,
                      int nover,int nunder,FILE *logfile_ptr)
{
  int ndeleted, nsaved, nspheres;

  double d_sqrd, x, x1, y, y1, z, z1;
  double closest_dist2, overlap;
  double under_radius, over_radius;

  struct sphere *under_sphere_ptr, *over_sphere_ptr;
  struct sphere *closest_sphere_ptr;
  struct sphere *last_under_sphere_ptr;
  struct sphere *next_under_sphere_ptr, *save_under_sphere_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Assigning under-sized spheres to clusters ...\n");
  ndeleted = 0;
  nsaved = 0;
  nspheres = nover;

  /* Initialise pointers */
  last_under_sphere_ptr = NULL;

  /* Get pointer to the first of the under-spheres */
  under_sphere_ptr = first_under_sphere_ptr;

  /* Loop through all the under-spheres to determine which need to be
     transferred to the over-size list */
  while (under_sphere_ptr != NULL)
    {
      /* Get pointer to next under-sphere */
      next_under_sphere_ptr = under_sphere_ptr->next_sphere_ptr;

      /* Save the current under-sphere pointer */
      save_under_sphere_ptr = under_sphere_ptr;

      /* Initialise this under-spherte's closest distance to an 
         over-sphere and pointer to over-sphere */
      closest_sphere_ptr = NULL;
      closest_dist2 = 0.0;

      /* Get sphere's radius */
      under_radius = under_sphere_ptr->radius;

      /* Get its coordinates */
      x = under_sphere_ptr->coords[0];
      y = under_sphere_ptr->coords[1];
      z = under_sphere_ptr->coords[2];

      /* Get pointer to the first of the over-spheres */
      over_sphere_ptr = first_over_sphere_ptr;

      /* Loop through all the most-recently added over-spheres to
         see if any overlap our current under-sphere */
      while (over_sphere_ptr != NULL)
        {
          /* Get over sphere's radius */
          over_radius = over_sphere_ptr->radius;

          /* Define the overlap distance as the radius of the 
             larger of the two spheres */
          overlap = over_radius;

          /* Get its coordinates */
          x1 = over_sphere_ptr->coords[0];
          y1 = over_sphere_ptr->coords[1];
          z1 = over_sphere_ptr->coords[2];

          /* Calculate distance squared */
          d_sqrd = (x - x1) * (x - x1) + (y - y1) * (y - y1)
            + (z - z1) * (z - z1);

          /* If the under-sphere overlaps the over-sphere then want
             to retain it */
          if (d_sqrd < overlap * overlap)
            {
              /* Check whether the sphere was previously found to
                 overlap an over-sphere */
              if (closest_sphere_ptr != NULL)
                {
                  /* Check if distance is smaller than what we had before */
                  if (d_sqrd < closest_dist2)
                    {
                      /* Store pointer and distance */
                      closest_sphere_ptr = over_sphere_ptr;
                      closest_dist2 = d_sqrd;
                    }
                }

              /* Otherwise, store this over-sphere */
              else
                {
                  /* Store pointer and distance */
                  closest_sphere_ptr = over_sphere_ptr;
                  closest_dist2 = d_sqrd;
                }
            }

          /* Get pointer to the next over-sphere */
          over_sphere_ptr = over_sphere_ptr->next_sphere_ptr;
        }

      /* If sphere to be retained, then transfer it to the over-size
         list and delete it from the under-size list by adjusting all 
         the pointers to by-pass it */
      if (closest_sphere_ptr != NULL)
        {
          /* Assign to the same cluster as the closest over-size sphere */
          under_sphere_ptr->cluster = closest_sphere_ptr->cluster;

          /* Add link from last over-sized sphere to the current
             under-sized one */
          if (last_over_sphere_ptr != NULL)
            last_over_sphere_ptr->next_sphere_ptr = under_sphere_ptr;
          last_over_sphere_ptr = under_sphere_ptr;
          under_sphere_ptr->next_sphere_ptr = NULL;

          /* Remove under-sphere from under-sized list by by-passing
             it */

          /* If this is the very first under-sphere, then set pointer to
             next one to bypass it */
          if (under_sphere_ptr == first_under_sphere_ptr)
            first_under_sphere_ptr = next_under_sphere_ptr;

          /* Otherwise, set last pointer around current sphere */
          else
            last_under_sphere_ptr->next_sphere_ptr 
              = next_under_sphere_ptr;

          /* Increment count of deleted spheres */
          nsaved++;
        }

      /* Otherwise, store the current under-sphere pointer */
      else
        last_under_sphere_ptr = save_under_sphere_ptr;

      /* Get pointer to the next under-sphere */
      under_sphere_ptr = next_under_sphere_ptr;
    }

  /* Update sphere counts */
  nspheres = nover + nsaved;
  ndeleted = (nover + nunder) - nspheres;

  /* Show numbers of spheres removed */
  fprintf(logfile_ptr,"  Number of spheres removed:          %8d\n",
          ndeleted);
  fprintf(logfile_ptr,"\n");
  fprintf(logfile_ptr,"  Number of spheres remaining:        %8d\n",
          nspheres);

  /* Return new number of spheres */
  return(nspheres);
}
/***********************************************************************

create_sphere_kdtree -  Create the kd-tree for the over-spheres

***********************************************************************/

void create_sphere_kdtree(struct sphere *first_over_sphere_ptr,
                          int nover,KdTree **kdt_ptr,
                          struct sphere ***sphere_index_ptr)
{
  int isphere;

  double **coords_array;

  struct sphere *sphere_ptr;
  struct sphere **sphere_ind_ptr;

  KdTree *kdtree_ptr;

  /* Create the index array of sphere pointers */
  sphere_ind_ptr
    = (struct sphere **) malloc(nover * sizeof(struct sphere *));

  /* Check that sufficient memory was allocated */
  if (sphere_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sphere ");
      printf("sort-index array\n");
      exit(1);
    }

  /* Create the array holding the x-values */
  coords_array = (double**) calloc(nover,sizeof(double));

  /* Check that sufficient memory was allocated */
  if (coords_array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for coords_array\n");
      exit(1);
    }

  /* Get pointer to the first sphere */
  isphere = 0;
  sphere_ptr = first_over_sphere_ptr;

  /* Loop through all spheres to place each in the kd-tree */
  while (sphere_ptr != NULL)
    {
      /* Store the sphere pointer in the sphere index array */
      sphere_ind_ptr[isphere] = sphere_ptr;

      /* Store the coords of this point in the array */
      coords_array[isphere] = (double *) sphere_ptr->coords;

      /* Increment sphere count */
      isphere++;

      /* Get pointer to the next sphere in the linked list */
      sphere_ptr = sphere_ptr->next_sphere_ptr;
    }

  /* Update number of spheres */
  nover = isphere;

  /* Create the kd-tree */
  kdtree_ptr = KdTree_create(coords_array,nover,3);

  /* Return pointer to start of kd-tree */
  *kdt_ptr = kdtree_ptr;

  /* Return the index array pointer */
  *sphere_index_ptr = sphere_ind_ptr;
}
/***********************************************************************

cluster_spheres  -  Cluster all the over-size spheres into separate
                    non-interacting clusters and assign any associated
                    under-spheres to those clusters

***********************************************************************/

int cluster_spheres(struct sphere *first_over_sphere_ptr,int nover,
                    FILE *logfile_ptr,double grid_sep,
                    KdTree *kdtree_ptr,struct sphere **sphere_index_ptr,
                    double max_gap_radius)
{
  int jsphere;
  int higher_cluster, ncluster, ndone;
  int new_cluster;

  double d_sqrd, x, x1, y, y1, z, z1;
  double overlap;
  double percen, plast, radius, other_radius;
  double inner_radius, outer_radius;

  struct sphere *sphere_ptr, *other_sphere_ptr;

  KdTreeQuery *kdquery_ptr;
  Region *region_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Clustering spheres ...\n");
  inner_radius = 0.0;
  ncluster = 0;
  ndone = 0;
  plast = 0.0;

  /* Get the first of the over-sized spheres */
  sphere_ptr = first_over_sphere_ptr;

  /* Loop through all the over-spheres to get spheres they overlap */
  while (sphere_ptr != NULL)
    {
      /* Calculate percentage progress */
      ndone++;
      percen = 100.0 * (double) ndone / (double) nover;
      if (percen - plast >= PCDIFF)
        {
          printf("%5.1f%% done\n",percen);
          plast = percen;
        }

      /* Initialise flag */
      new_cluster = FALSE;

      /* Get sphere's radius */
      radius = sphere_ptr->radius;

      /* Get its coordinates */
      x = sphere_ptr->coords[0];
      y = sphere_ptr->coords[1];
      z = sphere_ptr->coords[2];

      /* If sphere has no cluster number, then assign next one to it */
      if (sphere_ptr->cluster == 0)
        {
          /* Increment cluster number and assign to this sphere */
          ncluster++;
          sphere_ptr->cluster = ncluster;

          /* Set flag to indicate this is currently an orphan sphere */
          new_cluster = TRUE;
        }

      /* Calculate a search distance for locating all neighbouring
         spheres to find which penetrate this one */
      outer_radius = sphere_ptr->radius + max_gap_radius;

      /* Create the region of space to be searched for nearby
         vertices */
      region_ptr
        = Annulus_create(sphere_ptr->coords,inner_radius,outer_radius,3);

      /* Perform the search to retrieve all neighbouring spheres */
      kdquery_ptr = KdTree_query(kdtree_ptr,region_ptr);

      /* Loop through the neighbouring spheres to find which penetrate
         this one */
      while ((jsphere = KdTreeQuery_next(kdquery_ptr)) >= 0)
        {
          /* Get the sphere pointer returned by the query */
          other_sphere_ptr = sphere_index_ptr[jsphere];

          /* If this is not the same sphere as we're currently
             considering, and it corresponds to an over-sphere,
             calculate the distance between them */
          if (other_sphere_ptr != sphere_ptr)
            {
              /* Get other sphere's radius */
              other_radius = other_sphere_ptr->radius;

              /* Define the overlap distance as the sum of the two radii */
              /*      overlap = radius + other_radius; */
              /*      overlap = 1.01 * grid_sep; */
              overlap = (radius + other_radius) / 2;

              /* Get its coordinates */
              x1 = other_sphere_ptr->coords[0];
              y1 = other_sphere_ptr->coords[1];
              z1 = other_sphere_ptr->coords[2];

              /* Calculate distance squared */
              d_sqrd = (x - x1) * (x - x1) + (y - y1) * (y - y1)
                + (z - z1) * (z - z1);

              /* If the two spheres overlap then both belong to the same
                 cluster */
              if (d_sqrd < overlap * overlap)
                {
                  /* Easy case is where second sphere hasn't yet been assigned
                     to a cluster */
                  if (other_sphere_ptr->cluster == 0)
                    {
                      /* Assign it to our current cluster number */
                      other_sphere_ptr->cluster = sphere_ptr->cluster;

                      /* Reset flag */
                      new_cluster = FALSE;
                    }

                  /* Otherwise, sphere has already been assigned to a cluster */
                  else if (sphere_ptr->cluster != other_sphere_ptr->cluster)
                    {
                      /* If our sphere was a new sphere, then reassign its
                         cluster number */
                      if (new_cluster == TRUE)
                        {
                          /* Assign current sphere to other sphere's cluster */
                          sphere_ptr->cluster = other_sphere_ptr->cluster;

                          /* Reset flag */
                          new_cluster = FALSE;

                          /* Decrement cluster number */
                          ncluster--;
                        }

                      /* Otherwise, need to join two sphere-clusters into one */
                      else
                        {
                          /* Join the two clusters, using the smaller of the
                             two numbers */
                          if (sphere_ptr->cluster > other_sphere_ptr->cluster)
                            {
                              /* Store the higher cluster-number */
                              higher_cluster = sphere_ptr->cluster;

                              /* Join the clusters */
                              join_clusters(first_over_sphere_ptr,
                                            other_sphere_ptr->cluster,
                                            sphere_ptr->cluster);
                            }
                          else
                            {
                              /* Store the higher cluster-number */
                              higher_cluster = other_sphere_ptr->cluster;

                              /* Join the clusters */
                              join_clusters(first_over_sphere_ptr,
                                            sphere_ptr->cluster,
                                            other_sphere_ptr->cluster);
                            }

                          /* If we have renamed the current cluster number,
                             then reset */
                          if (higher_cluster == ncluster)
                            ncluster--;
                        }
                    }
                }
            }
        }

      /* Get pointer to the next sphere */
      sphere_ptr = sphere_ptr->next_sphere_ptr;
    }

  /* Show number of sphere clusters */
  fprintf(logfile_ptr,"  Number of clusters:                 %8d\n",
          ncluster);
  fprintf(logfile_ptr,"\n");

  /* Return number of clusters */
  return(ncluster);
}
/***********************************************************************

transfer_spheres  -  Assign to clusters any of the under-size spheres
                     that penetrate any of the over-size spheres and
                     transfer to the over-size sphere list

***********************************************************************/

int transfer_spheres(struct sphere *first_over_sphere_ptr,
                     struct sphere *last_over_sphere_ptr,
                     struct sphere *first_under_sphere_ptr,
                     int nover,int nunder,FILE *logfile_ptr,
                     KdTree *kdtree_ptr,struct sphere **sphere_index_ptr,
                     double max_gap_radius)
{
  int jsphere, ndeleted, nsaved, nspheres;

  double d_sqrd, x, x1, y, y1, z, z1;
  double closest_dist2, overlap;
  double under_radius, over_radius;
  double inner_radius, outer_radius;

  struct sphere *under_sphere_ptr, *over_sphere_ptr;
  struct sphere *closest_sphere_ptr;
  struct sphere *last_under_sphere_ptr;
  struct sphere *next_under_sphere_ptr, *save_under_sphere_ptr;

  KdTreeQuery *kdquery_ptr;
  Region *region_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Assigning under-sized spheres to clusters ...\n");
  inner_radius = 0.0;
  ndeleted = 0;
  nsaved = 0;
  nspheres = nover;

  /* Initialise pointers */
  last_under_sphere_ptr = NULL;

  /* Get pointer to the first of the under-spheres */
  under_sphere_ptr = first_under_sphere_ptr;

  /* Loop through all the under-spheres to determine which need to be
     transferred to the over-size list */
  while (under_sphere_ptr != NULL)
    {
      /* Get pointer to next under-sphere */
      next_under_sphere_ptr = under_sphere_ptr->next_sphere_ptr;

      /* Save the current under-sphere pointer */
      save_under_sphere_ptr = under_sphere_ptr;

      /* Initialise this under-sphere's closest distance to an 
         over-sphere and pointer to over-sphere */
      closest_sphere_ptr = NULL;
      closest_dist2 = 0.0;

      /* Get sphere's radius */
      under_radius = under_sphere_ptr->radius;

      /* Get its coordinates */
      x = under_sphere_ptr->coords[0];
      y = under_sphere_ptr->coords[1];
      z = under_sphere_ptr->coords[2];

      /* Calculate a search distance for locating all neighbouring
         spheres to find which penetrate this one */
      outer_radius = under_sphere_ptr->radius + max_gap_radius;

      /* Create the region of space to be searched for nearby
         vertices */
      region_ptr
        = Annulus_create(under_sphere_ptr->coords,inner_radius,
                         outer_radius,3);

      /* Perform the search to retrieve all neighbouring spheres */
      kdquery_ptr = KdTree_query(kdtree_ptr,region_ptr);

      /* Loop through the neighbouring spheres to find which penetrate
         this one */
      while ((jsphere = KdTreeQuery_next(kdquery_ptr)) >= 0)
        {
          /* Get the sphere pointer returned by the query */
          over_sphere_ptr = sphere_index_ptr[jsphere];

          /* If this is not the same sphere as we're currently
             considering, calculate the distance between them */
          if (over_sphere_ptr != under_sphere_ptr)
            {
              /* Get over sphere's radius */
              over_radius = over_sphere_ptr->radius;

              /* Define the overlap distance as the radius of the 
                 larger of the two spheres */
              overlap = over_radius;

              /* Get its coordinates */
              x1 = over_sphere_ptr->coords[0];
              y1 = over_sphere_ptr->coords[1];
              z1 = over_sphere_ptr->coords[2];

              /* Calculate distance squared */
              d_sqrd = (x - x1) * (x - x1) + (y - y1) * (y - y1)
                + (z - z1) * (z - z1);

              /* If the under-sphere overlaps the over-sphere then want
                 to retain it */
              if (d_sqrd < overlap * overlap)
                {
                  /* Check whether the sphere was previously found to
                     overlap an over-sphere */
                  if (closest_sphere_ptr != NULL)
                    {
                      /* Check if distance is smaller than what we had before */
                      if (d_sqrd < closest_dist2)
                        {
                          /* Store pointer and distance */
                          closest_sphere_ptr = over_sphere_ptr;
                          closest_dist2 = d_sqrd;
                        }
                    }

                  /* Otherwise, store this over-sphere */
                  else
                    {
                      /* Store pointer and distance */
                      closest_sphere_ptr = over_sphere_ptr;
                      closest_dist2 = d_sqrd;
                    }
                }
            }
        }

      /* If sphere to be retained, then transfer it to the over-size
         list and delete it from the under-size list by adjusting all 
         the pointers to by-pass it */
      if (closest_sphere_ptr != NULL)
        {
          /* Assign to the same cluster as the closest over-size sphere */
          under_sphere_ptr->cluster = closest_sphere_ptr->cluster;

          /* Add link from last over-sized sphere to the current
             under-sized one */
          if (last_over_sphere_ptr != NULL)
            last_over_sphere_ptr->next_sphere_ptr = under_sphere_ptr;
          last_over_sphere_ptr = under_sphere_ptr;
          under_sphere_ptr->next_sphere_ptr = NULL;

          /* Remove under-sphere from under-sized list by by-passing
             it */

          /* If this is the very first under-sphere, then set pointer to
             next one to bypass it */
          if (under_sphere_ptr == first_under_sphere_ptr)
            first_under_sphere_ptr = next_under_sphere_ptr;

          /* Otherwise, set last pointer around current sphere */
          else
            last_under_sphere_ptr->next_sphere_ptr 
              = next_under_sphere_ptr;

          /* Increment count of deleted spheres */
          nsaved++;
        }

      /* Otherwise, store the current under-sphere pointer */
      else
        last_under_sphere_ptr = save_under_sphere_ptr;

      /* Free up the memory taken for defining the region and
         making the query */
      Annulus_free(region_ptr);
      KdTreeQuery_free(kdquery_ptr);

      /* Get pointer to the next under-sphere */
      under_sphere_ptr = next_under_sphere_ptr;
    }

  /* Update sphere counts */
  nspheres = nover + nsaved;
  ndeleted = (nover + nunder) - nspheres;

  /* Show numbers of spheres removed */
  fprintf(logfile_ptr,"  Number of spheres removed:          %8d\n",
          ndeleted);
  fprintf(logfile_ptr,"\n");
  fprintf(logfile_ptr,"  Number of spheres remaining:        %8d\n",
          nspheres);

  /* Return new number of spheres */
  return(nspheres);
}
/***********************************************************************

create_arrays  -  Create the arrays for sphere-counts and sort-index

***********************************************************************/

void create_arrays(int **clust_sph,int **ind_clust,double **vol,
                   int nclusters)
{
  int icluster;
  int *index_clust, *sph_count;

  double *volume;

  /* Create the array of clustered sphere-counts */
  sph_count = (int *) malloc(nclusters * sizeof(int));
  if (sph_count == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sphere-count array\n");
      exit(1);
    }

  /* Repeat for their sort-index */
  index_clust = (int *) malloc(nclusters * sizeof(int));
  if (index_clust == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sort-index array\n");
      exit(1);
    }

  /* Create array for cluster volumes */
  volume = (double *) malloc(nclusters * sizeof(double));
  if (volume == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for volumes array\n");
      exit(1);
    }

  /* Initialise all arrays */
  for (icluster = 0; icluster < nclusters; icluster++)
    {
      index_clust[icluster] = icluster;
      sph_count[icluster] = 0;
      volume[icluster] = 0.0;
    }

  /* Return the array pointers */
  *ind_clust = index_clust;
  *clust_sph = sph_count;
  *vol = volume;
}
/***********************************************************************

count_spheres  -  Calculate the volume occupied by each cluster
                  of spheres by filling a grid of points

***********************************************************************/

void count_spheres(struct sphere *first_sphere_ptr,int nover,
                   int nclusters,int *sph_count,double *volume)
{
  int icluster, isphere;

  double radius, vol;

  struct sphere *sphere_ptr;

  /* Initialise sphere-count */
  isphere = 0;

  /* Get the first of the spheres */
  sphere_ptr = first_sphere_ptr;

  /* Loop through all the spheres to find which cluster they belong to */
  while (sphere_ptr != NULL)
    {
      /* Get this sphere's cluster */
      icluster = sphere_ptr->cluster;

      /* If a valid cluster-number, then increment appropriate count */
      if (icluster > 0 && icluster <= nclusters)
        {
          /* Increment count */
          sph_count[icluster - 1]++;

          /* Get this sphere's radius and calculate its volume */
          radius = sphere_ptr->radius;
          vol = radius * radius * radius;
          volume[icluster - 1] = volume[icluster - 1] + vol;
        }

      /* Increment gap-sphere count */
      isphere++;

      /* Get pointer to the next sphere */
      sphere_ptr = sphere_ptr->next_sphere_ptr;
    }
}
/***********************************************************************

sort_clusters  -  Sort clusters according to their numbers of spheres

***********************************************************************/

void sort_clusters(struct sphere *first_sphere_ptr,int nclusters,
                   int nover,struct parameters params)
{
  int icluster, iorder, jcluster;
  int *index_clust, *sph_count;
  int swap_cluster;
  int swapped;

  double swap_radius, swap_x, swap_y, swap_z;

  double *volume;

  struct sphere *sphere_ptr, *next_free_sphere_ptr, *next_sphere_ptr;
  struct sphere *start_sphere_ptr;

  /* Create the arrays for sphere-counts and sort-index */
  create_arrays(&sph_count,&index_clust,&volume,nclusters);

  /* Count the numbers of spheres belonging to each cluster */
  count_spheres(first_sphere_ptr,nover,nclusters,sph_count,volume);

  /* Sort the index-array using a simple bubble-sort such that the 
     clusters are ordered in descending order of their volume */

  /* Initialise swap flag */
  swapped = TRUE;

  /* Keep looping while any elements are still being swapped */
  while (swapped == TRUE)
    {
      /* Unset the swap flag */
      swapped = FALSE;

      /* Loop through all the clusters */
      for (iorder = 0; iorder < nclusters - 1; iorder++)
        {
          /* Get the clusters pointed to by this index element and the
             next one */
          icluster = index_clust[iorder];
          jcluster = index_clust[iorder + 1];

          /* If their volumes are the wrong way round, then need to
             swap them */
          // if (sph_count[icluster] < sph_count[jcluster])
          if (volume[icluster] < volume[jcluster])
            {
              /* Swap the index elements */
              index_clust[iorder] = jcluster;
              index_clust[iorder + 1] = icluster;

              /* Set flag to indicate that a swap took place in this loop */
              swapped = TRUE;
            }
        }
    }

  /* For first loop start at the very beginning of the spheres */
  start_sphere_ptr = first_sphere_ptr;

  /* Loop over all the clusters */
  for (iorder = 0; iorder < nclusters; iorder++)
    {
      /* Get the cluster-number corresponding to this sort position */
      icluster = index_clust[iorder] + 1;

      /* Initialise pointers */
      next_free_sphere_ptr = NULL;

      /* Get the first of the spheres */
      sphere_ptr = start_sphere_ptr;

      /* Loop through all the spheres to find those belonging to the
         current cluster */
      while (sphere_ptr != NULL)
        {
          /* If sphere belongs to the current cluster move it to the
             next free slot, swapping it which the sphere currently
             there */
          if (sphere_ptr->cluster == icluster)
            {
              /* If we have a valid free slot, then perform the swap */
              if (next_free_sphere_ptr != NULL)
                {
                  /* Store details of unwanted sphere */
                  swap_radius = next_free_sphere_ptr->radius;
                  swap_x = next_free_sphere_ptr->coords[0];
                  swap_y = next_free_sphere_ptr->coords[1];
                  swap_z = next_free_sphere_ptr->coords[2];
                  swap_cluster = next_free_sphere_ptr->cluster;

                  /* Move details across */
                  next_free_sphere_ptr->radius = sphere_ptr->radius;
                  next_free_sphere_ptr->coords[0] = sphere_ptr->coords[0];
                  next_free_sphere_ptr->coords[1] = sphere_ptr->coords[1];
                  next_free_sphere_ptr->coords[2] = sphere_ptr->coords[2];
                  next_free_sphere_ptr->cluster = sphere_ptr->cluster;

                  /* Transfer unwanted sphere into old slot */
                  sphere_ptr->radius = swap_radius;
                  sphere_ptr->coords[0] = swap_x;
                  sphere_ptr->coords[1] = swap_y;
                  sphere_ptr->coords[2] = swap_z;
                  sphere_ptr->cluster = swap_cluster;

                  /* Get the next sphere after the one in the previous
                     free slot */
                  next_sphere_ptr = next_free_sphere_ptr->next_sphere_ptr;

                  /* Check next sphere to see if it can be the next free
                     slot */
                  if (next_sphere_ptr != NULL)
                    {
                      /* If next sphere belongs to a different cluster
                         then can be the next free slot */
                      if (next_sphere_ptr->cluster != icluster)
                        next_free_sphere_ptr = next_sphere_ptr;

                      /* Otherwise, blank it out */
                      else
                        next_free_sphere_ptr = NULL;
                    }
                }
            }

          /* If sphere belongs to a different cluster, and don't yet have
             a free slot, then this is it */
          else if (next_free_sphere_ptr == NULL)
            next_free_sphere_ptr = sphere_ptr;

          /* Get pointer to the next sphere */
          sphere_ptr = sphere_ptr->next_sphere_ptr;
        }

      /* If have a valid free slot, this will be where the next loop
         will start */
      if (next_free_sphere_ptr != NULL)
        start_sphere_ptr = next_free_sphere_ptr;

      /* Otherwise, end here */
      else
        iorder = nclusters + 1;
    }

  /* Free memory taken by arrays */
  free(sph_count);
  free(index_clust);
}
/***********************************************************************

write_spheres  -  Write out the spheres to a PDB file for input into
                  SURFNET

***********************************************************************/

int write_spheres(char out_name[FILENAME_LEN],
                  struct sphere *first_sphere_ptr,
                  struct parameters params,FILE *logfile_ptr)
{
  char atom_name[5], chain, res_name[4], res_num[6];

  static char CLUSTER_CHAIN[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

  int atom_number, cluster, nout, last_cluster, ncluster, resno;
  int done;

  double occupancy;

  struct sphere *sphere_ptr;

  FILE *file_out;

  /* Initialise variables */
  fprintf(logfile_ptr,"Writing spheres to %s ...\n",out_name);
  strcpy(atom_name," C  ");
  strcpy(res_name,"SPH");
  chain = ' ';
  done = FALSE;
  last_cluster = -1;
  ncluster = 0;
  nout = 0;
  occupancy = 1.0;
  resno = 0;

  /* Open the output file, if required */
  if (params.from_www == FALSE && params.no_gif == FALSE)
    {
      if ((file_out = fopen(out_name,"w")) ==  NULL)
        {
          printf("\n*** Unable to open output PDB file\n");
          exit(1);
        }
    }

  /* Get pointer to the first sphere */
  sphere_ptr = first_sphere_ptr;

  /* Loop through all spheres to write out */
  while (sphere_ptr != NULL && done == FALSE)
    {
      /* Increment count of spheres written out */
      nout++;

      /* Get the cluster number */
      cluster = sphere_ptr->cluster;

      /* If the cluster number has changed, assign the next ordered
         cluster number to this sphere */
      if (cluster != last_cluster)
        {
          /* Get the next ordered cluster number */
          ncluster++;
        }

      /* If outputting only the top gap, then end if this is a smaller one */
      if (ncluster > params.ntop_clusters)
        done = TRUE;

      /* Otherwise, write out this sphere */
      else
        {
          /* Replace the cluster-number with the ordered cluster number */
          sphere_ptr->cluster = ncluster;

          /* Form the residue number */
          if (ncluster > 9999)
            sprintf(res_num,"%5d",ncluster);
          else
            sprintf(res_num,"%4d ",ncluster);

          /* Get the atom-number, checking for overflow */
          atom_number = nout;
          if (atom_number > 99999)
            atom_number = 99999;

          /* If assigning each sphere to a different residue, get atom
             and residue numbers */
          if (params.sph_resno == TRUE)
            {
              /* If this is a new cluster, reset residue number and set
                 chain */
              if (cluster != last_cluster)
                {
                  /* Reset residue number */
                  resno = 0;

                  /* Store the cluster number in the chain and occupancy
                     fields */
                  if (ncluster < strlen(CLUSTER_CHAIN))
                    chain = CLUSTER_CHAIN[ncluster - 1];
                  else
                    chain = ' ';
                  if (ncluster < 100)
                    occupancy = ncluster;
                  else
                    occupancy = 0;
                }

              /* Increment residue count and residue number */
              resno++;
              if (resno > 9999)
                sprintf(res_num,"%5d",resno);
              else
                sprintf(res_num,"%4d ",resno);
            }

          /* Write out the sphere record */
          if (params.from_www == FALSE && params.no_gif == FALSE)
            {
              fprintf(file_out,"ATOM  %5d ",atom_number);
              fprintf(file_out,"%4s ",atom_name);
              fprintf(file_out,"%3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                      res_name,chain,res_num,sphere_ptr->coords[0],
                      sphere_ptr->coords[1],sphere_ptr->coords[2],occupancy,
                      sphere_ptr->radius);
            }
        }

      /* Save the current cluster number */
      last_cluster = cluster;

      /* Get pointer to the next sphere */
      sphere_ptr = sphere_ptr->next_sphere_ptr;
    }

  /* Close the output file */
  if (params.from_www == FALSE && params.no_gif == FALSE)
    fclose(file_out);

  /* Show number of spheres written out */
  fprintf(logfile_ptr,"  Number of spheres written out:      %8d in %d",
          nout,ncluster);
  fprintf(logfile_ptr," clusters\n");

  /* Return the number of sphere clusters */
  return(ncluster);
}
/***********************************************************************

open_rasmol_file  -  Open the output RasMol coordinates output file

***********************************************************************/

void open_rasmol_file(FILE **out_ptr,char out_name[FILENAME_LEN])
{
  FILE *outfile_ptr;

  /* Open the output file, or write to stdout */
  if (out_name[0] == '\0')
    outfile_ptr = stdout;
  else
    {
      if ((outfile_ptr = fopen(out_name,"w")) ==  NULL)
        {
          printf("\n*** Unable to open output RasMol output file [%s]\n",
                 out_name);
          exit(1);
        }
    }

  /* Return the pointer to the output file */
  *out_ptr = outfile_ptr;
}
/***********************************************************************

write_surfnet_headers  -  Write out the header lines to the surfnet.par
                          file

***********************************************************************/

void write_surfnet_headers(FILE *file_ptr)
{
  /* Write header records */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"Parameters for SURFNET              (surfnet.par)\n");
  fprintf(file_ptr,"----------------------\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"TITLE\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"INPUT FILE\n");
  fprintf(file_ptr,"speedfill.raspdb\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"OUTPUT FILES\n");
}
/***********************************************************************

write_script_headers  -  Write all the header records to the output
                         RasMol script file

***********************************************************************/

void write_script_headers(FILE *outfile_ptr,char *pdb_name,
                          char *pdbsum_dir,struct parameters params)
{
  /* Write out the starting lines */
  fprintf(outfile_ptr,"#!rasmol -script\n");
  fprintf(outfile_ptr,"# Generated by SPEEDFILL - R A Laskowski, Aug 2000\n");
  fprintf(outfile_ptr,"\n");

  /* Show PDB code and PDB file name */
  fprintf(outfile_ptr,"# PDB code: %s\n",pdb_name);
  fprintf(outfile_ptr,"\n");

  /* Initilise and indicate that ATOM records to follow */
  fprintf(outfile_ptr,"zap\n");
  if (params.pdb_type == PDBSUM1 && params.relinks == FALSE)
    fprintf(outfile_ptr,"load %s/clefts/clefts.pdb\n",pdbsum_dir);
  else
    fprintf(outfile_ptr,"load inline\n");

  /* Set background colour */
  fprintf(outfile_ptr,"background black\n");
  fprintf(outfile_ptr,"set ambient 60\n");
  fprintf(outfile_ptr,"set specular off\n");
  fprintf(outfile_ptr,"slab off\n");

  /* Rotate by 180 about x-axis to be in same orientation as in PDB file */
  fprintf(outfile_ptr,"rotate x 180\n");
  fprintf(outfile_ptr,"\n");

  /* Switch all objects off */
  fprintf(outfile_ptr,"set bonds off\n");
  fprintf(outfile_ptr,"set axes off\n");
  fprintf(outfile_ptr,"set boundingbox off\n");
  fprintf(outfile_ptr,"set unitcell off\n");
  fprintf(outfile_ptr,"set bondmode and\n");
  fprintf(outfile_ptr,"dots off\n");
  fprintf(outfile_ptr,"\n");

  /* Initialise all object colours */
  fprintf(outfile_ptr,"# Initialise colours\n");
  fprintf(outfile_ptr,"select all\n");
  fprintf(outfile_ptr,"colour bonds none\n");
  fprintf(outfile_ptr,"colour backbone none\n");
  fprintf(outfile_ptr,"colour hbonds none\n");
  fprintf(outfile_ptr,"colour ssbonds none\n");
  fprintf(outfile_ptr,"colour ribbons none\n");

  /* Switch everything off */
  fprintf(outfile_ptr,"wireframe off\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"\n");

  /* Initialise RasMol's x-term output commands */
  fprintf(outfile_ptr,"echo \" \"\n");
  if (params.no_gaps == FALSE)
    {
      fprintf(outfile_ptr,"echo \"   +---------------------+\"\n");
      if (params.buried_vertices_only == FALSE)
        fprintf(outfile_ptr,"echo \"   |    Binding sites    |\"\n");
      else
        fprintf(outfile_ptr,"echo \"   |  Binding surfaces   |\"\n");
      fprintf(outfile_ptr,"echo \"   +---------------------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
    }
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_pymol_headers  -  Write the header records for the PyMOL script

***********************************************************************/

void write_pymol_headers(FILE *outfile_ptr,char *pdbsum_dir,
                         char *pdb_code,int relinks,char *clefts_pdb,
                         struct c_file_data file_name_ptr,
                         struct parameters params)
{
  char file_name[FILENAME_LEN], number_string[20];

  /* Write out the starting lines */
  fprintf(outfile_ptr,"#\n");
  if (params.pdb_type == PDBSUM1)
    fprintf(outfile_ptr,"# PDBsum1 ");
  else
    fprintf(outfile_ptr,"# PDBsum ");
  fprintf(outfile_ptr,"PyMOL script for speedfill output\n");
  fprintf(outfile_ptr,"#\n");

  /* Form the name of the coordinate file */
  if (params.pymol_www == FALSE)
    {
      /* File is in the PDBsum directory */
      strcpy(clefts_pdb,pdbsum_dir);
      if (params.pdb_type == PDBSUM1)
        strcat(clefts_pdb,"/clefts");
      strcat(clefts_pdb,"/clefts.pdb");

      /* Load the structure */
      fprintf(outfile_ptr,"load %s\n",clefts_pdb);
    }

  /* For the web version, need to call the catcoords script */
  else
    {
      strcpy(file_name,file_name_ptr.other_dir[URL_EBI]);
      strcat(file_name,file_name_ptr.other_dir[PDBSUM_CGI_DIR]);
      strcat(file_name,"/catcoords.pl?dir=");
      sprintf(number_string,"%s",params.user_id);
      strcat(file_name,number_string);
      strcat(file_name,"&prog=speedfill, PDBfile");

      /* Load the structure */
      fprintf(outfile_ptr,"load %s\n",file_name);
    }

  /* Write out the header records */
  fprintf(outfile_ptr,"set depth_cue=0\n");
  fprintf(outfile_ptr,"set ray_trace_fog=0\n");
  fprintf(outfile_ptr,"set cartoon_ring_mode,1\n");
  fprintf(outfile_ptr,"hide all\n");
  fprintf(outfile_ptr,"set auto_color, 0\n");

  /* Define the colours */
  fprintf(outfile_ptr,"set_color black,    [  0.000,  0.000,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color salmon,   [  0.941,  0.502,  0.502 ]\n");
  fprintf(outfile_ptr,"set_color red,      [  1.000,  0.000,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color green,    [  0.565,  1.000,  0.565 ]\n");
  fprintf(outfile_ptr,"set_color blue,     [  0.000,  0.000,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color yellow,   [  1.000,  1.000,  0.392 ]\n");
  fprintf(outfile_ptr,"set_color purple,   [  0.886,  0.490,  0.984 ]\n");
  fprintf(outfile_ptr,"set_color brown,    [  0.863,  0.824,  0.471 ]\n");
  fprintf(outfile_ptr,"set_color sky_blue, [  0.608,  0.886,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color skyblue,  [  0.608,  0.886,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color orange,   [  1.000,  0.400,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color cyan,     [  0.200,  0.800,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color grey,     [  0.584,  0.651,  0.718 ]\n");

  /* Parameters for surface rendering */
  fprintf(outfile_ptr,"set solvent_radius, 0.2\n");
  fprintf(outfile_ptr,"set transparency, 0.2\n");

  /* Write out background colour */
  if (params.pymol_www == FALSE)
    fprintf(outfile_ptr,"bg_color white\n");
  else
    fprintf(outfile_ptr,"bg_color black\n");
}
/***********************************************************************

write_script_coords  -  Write out the script commands to define the
                        atom coordinates just written out

***********************************************************************/

void write_script_coords(FILE *outfile_ptr,int loop,
                         int first_atom_number,int last_atom_number,
                         struct parameters params)
{
  char object_name[20];

  /* Print comment */
  fprintf(outfile_ptr,"\n");
  if (loop == 0)
    {
      fprintf(outfile_ptr,"# Protein\n");
      strcpy(object_name,"protall");
    }
  else
    {
      fprintf(outfile_ptr,"# Ligand\n");
      strcpy(object_name,"ligall");
    }

  /* Define the atom-number range for the current object */
  fprintf(outfile_ptr,"define %s atomno >= %d & atomno <= %d\n",
          object_name,first_atom_number,last_atom_number);

  /* Select and render */
  fprintf(outfile_ptr,"select %s\n",object_name);
  if (loop == 0)
    fprintf(outfile_ptr,"wireframe on\n");
  else
    {
      fprintf(outfile_ptr,"spacefill 1.2\n");
    }
  fprintf(outfile_ptr,"colour cpk\n");

  /* Show details in RasMol command window for user's info */
  fprintf(outfile_ptr,"echo \"   %s:   ",object_name);
  if (loop == 0)
    fprintf(outfile_ptr,"shown as wireframe\"\n");
  else
    fprintf(outfile_ptr,"shown as spacefill\"\n");

  /* If colouring by conservation, colour by B-value */
  if (params.col_by_conservation == TRUE && loop == 0)
    {
      fprintf(outfile_ptr,
              "select protall & temperature >= 0 & temperature <= 69\n");
      fprintf(outfile_ptr,"colour [255,255,255]\n");
      fprintf(outfile_ptr,
              "select protall & temperature > 69 & temperature <= 80\n");
      fprintf(outfile_ptr,"colour [255,102,0]\n");
      fprintf(outfile_ptr,
              "select protall & temperature > 80 & temperature <= 90\n");
      fprintf(outfile_ptr,"colour [255,105,180]\n");
      fprintf(outfile_ptr,
              "select protall & temperature > 90 & temperature <= 100\n");
      fprintf(outfile_ptr,"colour [255,0,0]\n");
    }
}
/***********************************************************************

write_pymol_script  -  Write out the PyMOL script commands to define the
                       atom coordinates just written out

***********************************************************************/

void write_pymol_script(FILE *outfile_ptr,int loop,
                        int first_atom_number,int last_atom_number,
                        int nclusters,struct parameters params)
{
  char object[20], object_type[20];

  /* If nothing to show, return */
  if (last_atom_number < first_atom_number)
    return;

  /* Print comment */
  fprintf(outfile_ptr,"\n");
  if (loop == 0)
    {
      fprintf(outfile_ptr,"# Protein\n");
      strcpy(object,"protall");
    }
  else
    {
      fprintf(outfile_ptr,"# Ligand\n");
      strcpy(object,"ligall");
    }

  /* Define the atom-number range for the current object */
  fprintf(outfile_ptr,"create %s, id %d-%d\n",
          object,first_atom_number,last_atom_number);

  /* First loop is protein */
  if (loop == 0)
    {
      /* Show as wireframe */
      fprintf(outfile_ptr,"show sticks, %s\n",object);

      /* Colour by atom type */
      fprintf(outfile_ptr,"util.cbag (%s and not name n+ca+c+o)\n",
              object);
      fprintf(outfile_ptr,
              "color grey, %s and elem C and not name n+ca+c+o\n",
              object);
    }

  /* Show ligands in spacefill */
  else
    {
      /* Colour by atom type */
      fprintf(outfile_ptr,"util.cbag (%s and not elem C)\n",
              object);
      if (params.pymol_www == TRUE)
        fprintf(outfile_ptr,"color grey, %s and elem C\n",object);
      else
        fprintf(outfile_ptr,"color black, %s and elem C\n",object);

      /* Show as spheres */
      fprintf(outfile_ptr,"show spheres, %s\n",object);
    }
}
/***********************************************************************

write_surfnet_line  -  Write out the line to the surfnet.par file for
                       protein or ligand

***********************************************************************/

void write_surfnet_line(FILE *file_ptr,int loop,
                        int first_atom_number,int last_atom_number)
{
  char object_name[20], object_type;

  double object_size;

  /* Print comment */
  if (loop == 0)
    {
      object_type = 'B';
      object_size = 0.1;
      strcpy(object_name,"protein");
    }
  else
    {
      object_type = 'A';
      object_size = 1.6;
      strcpy(object_name,"ligand ");
    }

  /* Check that have at least one atom in this object */
  if (first_atom_number <= last_atom_number)
    {
      /* Define the atom-number range for the current object */
      fprintf(file_ptr,"%s      %c %6d   %6d  100.0",object_name,
              object_type,first_atom_number,last_atom_number);

      /* Define colours */
      fprintf(file_ptr,"       1        10      %4.2f\n",object_size);
    }
}
/***********************************************************************

write_out_coords  -  Write out the atom and ligand coordinates to the
                     RasMol script and coordinate output files

***********************************************************************/

void write_out_coords(FILE *outfile_ptr,FILE *scriptfile_ptr,
                      FILE *pmlfile_ptr,struct atom *first_atom_ptr,
                      struct atom *first_lig_atom_ptr,int *atom_no,
                      FILE *logfile_ptr,int nclusters,
                      struct parameters params,FILE *surfnet_file_ptr)
{
  static char *keyword[2] = { "ATOM  ", "HETATM" };

  int atom_number, first_atom_number, last_atom_number, nout, loop;
  int done;

  double bvalue, occupancy;

  struct atom *atom_ptr;

  /* Initialise variables */
  fprintf(logfile_ptr,"Writing atoms to RasMol output files ...\n");
  atom_number = *atom_no;
  done = FALSE;
  nout = 0;
  occupancy = 1.0;

  /* Loop twice, once to write out the atom coordinates, and then
     to write out the ligand coordinates */
  for (loop = 0; loop < 2; loop++)
    {
      /* Get pointer to the first atom or ligand-atom, depending on which
         loop this is */
      if (loop == 0)
        atom_ptr = first_atom_ptr;
      else
        atom_ptr = first_lig_atom_ptr;

      /* Initialise start and end atoms numbers */
      first_atom_number = atom_number + 1;
      last_atom_number = atom_number;

      /* Loop through all atoms to be written out */
      while (atom_ptr != NULL)
        {
          /* Write out only if atom is within the grid-array bounds */
          if (loop == 0 || atom_ptr->residue_ptr->in_range == TRUE)
            {
              /* Increment count of atoms written out */
              nout++;

              /* If colouring by residue conservation, get the 
                conservation score */
              if (params.col_by_conservation == TRUE)
                bvalue = atom_ptr->residue_ptr->conservation_score;
              else
                bvalue = atom_ptr->bvalue;

              /* Get the atom-number, checking for overflow */
              atom_number++;
              if (atom_number > 99999)
                atom_number = 99999;

              /* Update the last atom number */
              last_atom_number = atom_number;

              /* Write out the atom record */
              fprintf(outfile_ptr,"%s%5d ",keyword[loop],atom_number);
              fprintf(outfile_ptr,"%4s ",atom_ptr->atom_name);
              fprintf(outfile_ptr,"%3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                      atom_ptr->residue_ptr->res_name,
                      atom_ptr->residue_ptr->chain,
                      atom_ptr->residue_ptr->res_num,
                      atom_ptr->coords[0],atom_ptr->coords[1],
                      atom_ptr->coords[2],occupancy,bvalue);
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next_atom_ptr;
        }

      /* Write a TER record out */
      fprintf(outfile_ptr,"TER   \n");

      /* Write out the RasMol script commands for the atoms just
         written out */
      if (scriptfile_ptr != NULL)
        write_script_coords(scriptfile_ptr,loop,first_atom_number,
                            last_atom_number,params);

      /* For PyMOL script, write to the .pml file */
      if (pmlfile_ptr != NULL)
        write_pymol_script(pmlfile_ptr,loop,first_atom_number,
                           last_atom_number,nclusters,params);

      /* If writing to a surfnet.par file, write out the line for this
         object */
      if (params.write_surfnet_input == TRUE)
        write_surfnet_line(surfnet_file_ptr,loop,first_atom_number,
                           last_atom_number);
    }

  /* Show number of atoms written out */
  fprintf(logfile_ptr,"  Total number of atoms written out:  %8d\n",nout);

  /* Return the new atom number */
  *atom_no = atom_number;
}
/***********************************************************************

create_vertex_record  -  Create and initialise a new vertex record

***********************************************************************/

struct vertex *create_vertex_record(struct vertex **fst_vertex_ptr,
                                    struct vertex **lst_vertex_ptr,
                                    double coords[3],int vertex_no,
                                    int edge_no,struct atom *atom_ptr)
{
  int i;

  struct vertex *vertex_ptr, *first_vertex_ptr, *last_vertex_ptr;

  /* Initialise vertex pointers */
  vertex_ptr = NULL;
  first_vertex_ptr = *fst_vertex_ptr;
  last_vertex_ptr = *lst_vertex_ptr;

  /* Allocate memory for structure to hold vertex info */
  vertex_ptr = (struct vertex *) malloc(sizeof(struct vertex));
  if (vertex_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct vertex\n");
      exit (1);
    }

  /* If this is the very first vertex, then store pointer */
  if (first_vertex_ptr == NULL)
    first_vertex_ptr = vertex_ptr;

  /* Add link from previous vertex to the current one */
  if (last_vertex_ptr != NULL)
    last_vertex_ptr->next_vertex_ptr = vertex_ptr;
  last_vertex_ptr = vertex_ptr;

  /* Store vertex's coordinates */
  for (i = 0; i < 3; i++)
    {
      vertex_ptr->coords[i] = coords[i];
      vertex_ptr->normal[i] = 0.0;
      vertex_ptr->next_sorted_vertex_ptr[i] = NULL;
      vertex_ptr->sort_pos[i] = -1;
    }

  /* Initialise the other elements of the structure */
  vertex_ptr->vertex_no = vertex_no;
  vertex_ptr->ratio = 0.0;
  vertex_ptr->bvalue = 0.0;
  vertex_ptr->deleted = FALSE;
  vertex_ptr->cluster = 0;
  vertex_ptr->edge_no = edge_no;
  vertex_ptr->npolygons = 0;
  vertex_ptr->exposed = FALSE;
  vertex_ptr->atom_ptr = atom_ptr;
  vertex_ptr->first_polist_ptr = NULL;
  vertex_ptr->last_polist_ptr = NULL;
  vertex_ptr->next_vertex_list_ptr = NULL;
  vertex_ptr->next_vertex_ptr = NULL;

  /* Return current pointers */
  *fst_vertex_ptr = first_vertex_ptr;
  *lst_vertex_ptr = last_vertex_ptr;

  /* Return new vertex pointer */
  return(vertex_ptr);
}
/***********************************************************************

read_vertices  -  Read through the binding-site analysis file to
                  find the top few binding sites that most closely
                  resemble this one

***********************************************************************/

int read_vertices(char *file_name,struct vertex **fst_vertex_ptr,
                  struct atom *first_atom_ptr,int natoms,int *nclust,
                  double *volume)
{
  char input_line[LINELEN + 1], number_string[20];
  char burexp, atom_name[5], chain, res_name[4], res_num[6];
  char element[3], wanted_chain;

  char *string_ptr;

  int atom_number, iatom, icluster, icoord, len, nclusters, nvertices;
  int edge_no, rec_type;
  int found, have_file, started;

  double bvalue, depth, cluster, coords[3], ratio, x, y, z;

  struct atom *atom_ptr;
  struct residue *residue_ptr;
  struct vertex *vertex_ptr;
  struct vertex *first_vertex_ptr, *last_vertex_ptr;

  struct atom **atom_index_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  edge_no = 0;
  have_file = TRUE;
  icluster = 0;
  nclusters = 0;
  nvertices = 0;
  started = FALSE;
  wanted_chain = '\0';

  /* Initialise pointers */
  first_vertex_ptr = NULL;
  last_vertex_ptr = NULL;

  /* Create an index array for the atoms to be read in from the
     speedfill.depth file */
  atom_index_ptr
    = (struct atom **) malloc(natoms * sizeof(struct atom *));

  /* Check that sufficient memory was allocated */
  if (atom_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sort-index ");
      printf("in read_vertices\n");
      exit(1);
    }

  /* Initialise the array of atom pointers */ 
  for (iatom = 0; iatom < natoms; iatom++)
    atom_index_ptr[iatom] = NULL;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If this is an atom record, store its cluster
             number and depth details */
          if (!strncmp(input_line,"ATOM  ",6))
            {
              /* Get the atom details from this line */
              r_get_atom_details(input_line,&atom_number,atom_name,
                                 &chain,res_name,res_num,&x,&y,&z,
                                 &depth,&cluster,element,&rec_type,
                                 wanted_chain);

              /* Get pointer to first atom record */
              atom_ptr = first_atom_ptr;
              found = FALSE;

              /* Loop until find the matching atom */
              while (atom_ptr != NULL && found == FALSE)
                {
                  /* Get this atom's residue */
                  residue_ptr = atom_ptr->residue_ptr;

                  /* Check whether we have found the correct atom
                     record */
                  if (!strcmp(residue_ptr->res_num,res_num) &&
                      !strcmp(atom_ptr->atom_name,atom_name) &&
                      !strcmp(residue_ptr->res_name,res_name) &&
                      chain == residue_ptr->chain)
                    {
                      /* Set flag that we have a match */
                      found = TRUE;

                      /* Store pointer to atom in the index array */
                      if (atom_number -1 < natoms)
                        atom_index_ptr[atom_number - 1] = atom_ptr;

                      /* Store the cluster number and depth */
                      atom_ptr->icluster = (int) (cluster + 0.5);
                      atom_ptr->depth = depth;
                    }

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next_atom_ptr;
                }
            }

          /* Otherwise, check for the start of the vertex list */
          else if (!strncmp(input_line,"REMARK  Vertices:",17))
            {
              /* Set flag to say we are at the start of the list of
                 vertices */
              started = TRUE;
            }

          /* Check for start of next cluster */
          else if (!strncmp(input_line,"Cluster:",8))
            {
              /* Read in the cluster number */
              strncpy(number_string,input_line+9,4);
              number_string[4] = '\0';
              icluster = atoi(number_string);

              /* Retrieve the cluster's volume */
              string_ptr = strstr(input_line,"Volume:");
              if (string_ptr != NULL)
                {
                  strcpy(number_string,string_ptr+8);
                  volume[icluster] = atof(number_string);
                }

              /* Increment cluster count */
              nclusters++;
            }

          /* If have started reading in the vertices, read in coords */
          else if (started == TRUE && len > 0)
            {
              /* Get whether the vertex is buried or exposed */
              burexp = input_line[0];

              /* Read in the atom number associated with this vertex */
              strncpy(number_string,input_line+1,5);
              number_string[5] = '\0';
              iatom = atoi(number_string) - 1;

              /* Read in the vertex coordinates */
              for (icoord = 0; icoord < 3; icoord++)
                {
                  strncpy(number_string,input_line+6+icoord*8,8);
                  number_string[8] = '\0';
                  coords[icoord] = atof(number_string);
                }

              /* Read in the vertex's burial ratio */
              strncpy(number_string,input_line+30,8);
              number_string[8] = '\0';
              ratio = atof(number_string);

              /* Read in the vertex's B-value */
              strncpy(number_string,input_line+38,8);
              number_string[8] = '\0';
              bvalue = atof(number_string);

              /* If vertex linked to an atom, retrieve the
                 corresponding atom pointer */
              if (iatom >= 0 && iatom < natoms)
                atom_ptr = atom_index_ptr[iatom];
              else
                atom_ptr = NULL;

              /* Create a new vertex record */
              vertex_ptr
                = create_vertex_record(&first_vertex_ptr,
                                       &last_vertex_ptr,coords,
                                       nvertices,edge_no,atom_ptr);

              /* If vertex is exposed, set flag */
              if (burexp == 'e')
                vertex_ptr->exposed = TRUE;

              /* Store the cluster number and burial ratio */
              vertex_ptr->cluster = icluster;
              vertex_ptr->ratio = ratio;
              vertex_ptr->bvalue = bvalue;

              /* Increment count of vertices */
              nvertices++;
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return vertex pointer */
  *fst_vertex_ptr = first_vertex_ptr;

  /* Return number of clusters read in */
  *nclust = nclusters;

  /* Return number of vertices stored */
  return(nvertices);
}
/***********************************************************************

create_vertex_index_arrays  -  Create the arrays holding the x-, y- and
                               z- sort-indices for the vertices

***********************************************************************/

void create_vertex_index_arrays(struct vertex *first_sorted_vertex_ptr[3],
                                struct vertex ***vertex_index_ptr,
                                int nvertices)
{
  int icoord, ivertex, ipos;
  struct vertex **vertex_ind_ptr;
  struct vertex *vertex_ptr;

  /* Create the three index arrays as one large array */
  vertex_ind_ptr
    = (struct vertex **) malloc(3 * nvertices * sizeof(struct vertex *));

  /* Check that sufficient memory was allocated */
  if (vertex_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for sort-index array\n");
      exit(1);
    }

  /* Loop over the three coordinates */
  for (icoord = 0; icoord < 3; icoord++)
    {
      /* Get pointer to the first vertex */
      ivertex = 0;
      vertex_ptr = first_sorted_vertex_ptr[icoord];

      /* Loop through all vertices to get current index array to point to
         each vertex in turn */
      while (vertex_ptr != NULL)
        {
          /* Get array position */
          ipos = icoord * nvertices + ivertex;

          /* Initialise the index array */
          vertex_ind_ptr[ipos] = vertex_ptr;

          /* Store this vertex's position in the list */
          vertex_ptr->sort_pos[icoord] = ivertex;

          /* Get pointer to the next vertex in the linked list */
          vertex_ptr = vertex_ptr->next_sorted_vertex_ptr[icoord];
          ivertex++;
        }
    }

  /* Return the array pointer */
  *vertex_index_ptr = vertex_ind_ptr;
}
/***********************************************************************

binary_search  -  Search for the given coordinate value in the given
                  x-, y- or z-coordinate index and return the first array
                  position whose value is higher than that supplied.
                  If the value supplied is larger than any value in the
                  array, then -1 is returned.

***********************************************************************/

int binary_search(struct vertex **vertex_index_ptr,int nvertices,
                  double search_value,int icoord)
{
  int bottom_pos, ipoint, ipos, new_pos, top_pos;
  int found;

  struct vertex *vertex_ptr;

  /* Make starting point half-way into the index */
  ipos = nvertices / 2;
  bottom_pos = 0;
  top_pos = nvertices - 1;

  /* Initialise variables */
  found = FALSE;

  /* Loop until required array position found */
  while (found == FALSE)
    {
      /* Get array position */
      ipoint = icoord * nvertices + ipos;

      /* Get the associated vertex pointer */
      vertex_ptr = vertex_index_ptr[ipoint];

      /* If array value is equal to or higher than the search value,
         check whether it is the value we want */
      if (vertex_ptr->coords[icoord] >= search_value)
        {
          /* If this is the first vertex in the sorted list, then
             this is the one we want */
          if (ipos == 0)
            found = TRUE;

          /* Otherwise, check whether the previous vertex's coordinate
             is below the search value */
          else
            {
              /* Get the pointer to the previous vertex */
              vertex_ptr = vertex_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after */
              if (vertex_ptr->coords[icoord] < search_value)
                found = TRUE;

              /* Otherwise, we need to search earlier in the list */
              else
                {
                  /* Save the current position as the uppermost position
                     to search on */
                  top_pos = ipos;

                  /* Set pointer to halfway between the last search
                     low-search and here */
                  ipos = (ipos + bottom_pos) / 2;
                }
            }

        }

      /* If current vertex's coords are below the search value, have to
         look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_pos = ipos;

          /* If have reached the top of the list, then our value is
             higher than any value in the list, so return -1 */
          if (ipos >= nvertices - 1)
            {
              /* Set ipos to -1 and flag as done */
              ipos = -1;
              found = TRUE;
            }

          /* Otherwise, can search further up the list */
          else
            {
              /* Set pointer to halfway between the last search high-search
                 and here */
              new_pos = (ipos + top_pos) / 2;

              /* Adjust if going back to the same point */
              if (new_pos == ipos)
                new_pos++;
              ipos = new_pos;
            }
        }
    }

  /* Return the required position in the index array */
  return(ipos);
}
/***********************************************************************

get_start_end  -  Get start- and end-vertices given the start and end-coords
                  to span

***********************************************************************/

int get_start_end(struct vertex **vertex_index_ptr,int nvertices,
                  int icoord,double value,double max_radius,
                  int start_point[3],int end_point[3])
{
  int nsearch;

  double from_coord, to_coord;

  /* Initialise variables */
  start_point[icoord] = 0;
  end_point[icoord] = -1;

  /* Determine range of coords required */
  from_coord = value - max_radius;
  to_coord = value + max_radius;

  /* Perform a binary search to locate the starting vertex for the
     lower limit of the range */
  start_point[icoord] = binary_search(vertex_index_ptr,nvertices,from_coord,
                                      icoord);

  /* Repeat for upper limit of the range */
  end_point[icoord] = binary_search(vertex_index_ptr,nvertices,to_coord,
                                    icoord);

  /* If both ends of the range are off either end, then there are no
     vertices within the range */
  if (start_point[icoord] == -1 || end_point[icoord] == 0)
    {
      /* Set start and end points such that no vertices are processed */
      start_point[icoord] = 0;
      end_point[icoord] = -1;
    }

  /* Check for upper-range value being off the end */
  else if (end_point[icoord] == -1)
    end_point[icoord] = nvertices - 1;

  /* Otherwise, reduce end position by 1 */
  else
    end_point[icoord]--;

  /* Calculate number of vertices in list */
  if (start_point[icoord] > end_point[icoord])
    nsearch = 0;
  else
    nsearch = end_point[icoord] - start_point[icoord] + 1;

  /* Return number of vertices in list */
  return(nsearch);
}
/***********************************************************************

transfer_vertex_list  -  Transfer vertices from the x-sorted list to the
                         all-purpose sort-stack

***********************************************************************/

int transfer_vertex_list(struct vertex **vertex_index_ptr,int nvertices,
                         int start_point[3],int end_point[3],
                         int ref_vertex_number,
                         struct vertex **fst_vertex_list_ptr,int filter)
{
  int ipoint, nstored;
  int wanted;

  struct vertex *vertex_ptr;
  struct vertex *first_vertex_list_ptr, *last_vertex_list_ptr;

  /* Initialise variables */
  nstored = 0;
  first_vertex_list_ptr = NULL;
  last_vertex_list_ptr = NULL;

  /* Loop over the vertices that have been selected on their x-coords */
  for (ipoint = start_point[0]; ipoint < end_point[0] + 1; ipoint++)
    {
      /* Get this vertex's pointer */
      vertex_ptr = vertex_index_ptr[ipoint];

      /* Set wanted flag */
      wanted = TRUE;

      /* If filtering out lower-numbered vertices, check that vertex
         number is higher than that of the one being compared against */
      if (vertex_ptr->deleted == TRUE ||
          (filter == TRUE && vertex_ptr->vertex_no <= ref_vertex_number))
        wanted = FALSE;
          
      /* Filter out vertex if it doesn't lie within the acceptable
         range of y-values */
      if (vertex_ptr->sort_pos[1] < start_point[1] ||
          vertex_ptr->sort_pos[1] > end_point[1])
        wanted = FALSE;

      /* Ditto for z-values */
      else if (vertex_ptr->sort_pos[2] < start_point[2] ||
               vertex_ptr->sort_pos[2] > end_point[2])
        wanted = FALSE;

      /* If vertex is wanted, then store */
      if (wanted == TRUE)
        {
          /* Get the location of where to store the pointer */
          /*          jpoint = 3 * nvertices + nstored; */

          /* Store the pointer */
          /*          vertex_index_ptr[jpoint] = vertex_ptr; */
          
          /* If this is the very first vertex, then store pointer */
          if (first_vertex_list_ptr == NULL)
            first_vertex_list_ptr = vertex_ptr;
  
          /* Add link from previous vertex to the current one */
          if (last_vertex_list_ptr != NULL)
            last_vertex_list_ptr->next_vertex_list_ptr = vertex_ptr;
          last_vertex_list_ptr = vertex_ptr;

          /* Vertex is current end of list */
          vertex_ptr->next_vertex_list_ptr = NULL;
  
          /* Increment count of vertices stored */
          nstored++;
        }
    }

  /* Return pointer to first vertex in the selected list */
  *fst_vertex_list_ptr = first_vertex_list_ptr;

  /* Return number of vertices stored in filtered list */
  return(nstored);
}
/***********************************************************************

get_neighbour_list  -  Get list of vertices within the gap-span distance
                       of the current vertex

***********************************************************************/

int get_neighbour_list(struct vertex *vertex_ptr,
                       struct vertex **vertex_index_ptr,int nvertices,
                       double reach_dist,
                       struct vertex **fst_vertex_list_ptr,int filter)
{
  int icoord, n, nlist;
  int end_point[3], start_point[3];

  struct vertex *first_vertex_list_ptr;

  /* Initialise variables */
  nlist = 0;
  first_vertex_list_ptr = NULL;

  /* Get start- and end-vertices in the sorted x-, y- and z-coord
     indices */
  for (icoord = 0; icoord < 3; icoord++)
    {
      n = get_start_end(vertex_index_ptr,nvertices,icoord,
                        vertex_ptr->coords[icoord],reach_dist,start_point,
                        end_point);
      nlist = nlist + n;
    }

  /* If have any vertices, continue */
  if (nlist > 0)
    {
      /* Transfer vertices from the x-sorted list to the all-purpose
         sort-stack */
      nlist
        = transfer_vertex_list(vertex_index_ptr,nvertices,start_point,
                               end_point,vertex_ptr->vertex_no,
                               &first_vertex_list_ptr,filter);
    }

  /* Return pointer to first vertex in the selected list */
  *fst_vertex_list_ptr = first_vertex_list_ptr;

  /* Return number of vertex pointers in the list */
  return(nlist);
}
/***********************************************************************

filterx_vertices  -  Filter out vertices that are too close to one
                     another for display in RasMol (Old method)

***********************************************************************/

int filterx_vertices(struct vertex *first_vertex_ptr,int nvertices,
                     struct vertex **vertex_index_ptr,
                     int buried_vertices_only)
{
  int ndeleted, nlist;
  int ndel, nstop;
  int done, filter, pair_combined;

  double dist_sqrd, x, y, z, other_x, other_y, other_z;
  double too_close_dist, too_close_sqrd;
  double reach_dist;

  struct vertex *vertex_ptr, *other_vertex_ptr;
  struct vertex *first_vertex_list_ptr;


  /* Initialise variables */
  done = FALSE;
  filter = TRUE;
  ndeleted = 0;
  too_close_dist = TOO_CLOSE_DIST;
  if (buried_vertices_only == TRUE)
    too_close_dist = TOO_CLOSE_BV;
  too_close_sqrd = too_close_dist * too_close_dist;
  reach_dist = 1.5 * too_close_dist;

  /* Calculate a reasonable time to stop the filtering as being after
     any loop during which fewer than 1% of the vertices were deleted */
  nstop = (int) (0.01 * nvertices);

  /* Loop until all very short vertex-vertex distances have been
     eliminated */
  while (done == FALSE)
    {
      /* Initialise count of vertices deleted in this loop */
      ndel = 0;

      /* Initialise flag saying whether any vertices combined in
         current loop */
      pair_combined = FALSE;

      /* Get pointer to the first vertex */
      vertex_ptr = first_vertex_ptr;

      /* Loop through all vertices */
      while (vertex_ptr != NULL)
        {
          /* Get the coordinates of this vertex */
          x = vertex_ptr->coords[0];
          y = vertex_ptr->coords[1];
          z = vertex_ptr->coords[2];

          /* Check whether vertex has already been deleted */
          if (vertex_ptr->deleted == FALSE)
            {
              /* Get list of neighbouring vertices that might be within
                 joining distance */
              nlist
                = get_neighbour_list(vertex_ptr,vertex_index_ptr,
                                     nvertices,reach_dist,
                                     &first_vertex_list_ptr,filter);
            }

          /* Otherwise, there is no list to consider */
          else
            {
              nlist = 0;
              first_vertex_list_ptr = NULL;
            }

          /* Get pointer to the first of the neighbouring vertices */
          other_vertex_ptr = first_vertex_list_ptr;

          /* Loop through the neighbour-list of nearby vertices to
             knock out those that are too close */
          while (other_vertex_ptr != NULL)
            {
              /* Get the coordinates of this vertex */
              other_x = other_vertex_ptr->coords[0];
              other_y = other_vertex_ptr->coords[1];
              other_z = other_vertex_ptr->coords[2];

              /* Calculate the distance between the two vertices */
              dist_sqrd = (x - other_x) * (x - other_x)
                + (y - other_y) * (y - other_y)
                + (z - other_z) * (z - other_z);

              /* If the distance is less than the cut-off, then
                 want to replace the two vertices with a single one
                 half-way between them */
              if (dist_sqrd < too_close_sqrd)
                {
                  /* Combine the vertex coords to shift the current
                     vertex to halfway */
                  vertex_ptr->coords[0] = (x + other_x) / 2.0;
                  vertex_ptr->coords[1] = (y + other_y) / 2.0;
                  vertex_ptr->coords[2] = (z + other_z) / 2.0;

                  /* If current vertex has no associated atom, then transfer
                     atom from vertex that is about to be deleted */
                  if (vertex_ptr->atom_ptr == NULL)
                    vertex_ptr->atom_ptr = other_vertex_ptr->atom_ptr;

                  /* Delete the other vertex */
                  other_vertex_ptr->deleted = TRUE;

                  /* Reduce count of vertices */
                  ndel++;
                  ndeleted++;

                  /* Set flag that vertices combined in this pass */
                  pair_combined = TRUE;
                }

              /* Get pointer to next vertex in the neightbour list */
              other_vertex_ptr = other_vertex_ptr->next_vertex_list_ptr;
            }

          /* Get pointer to the next vertex */
          vertex_ptr = vertex_ptr->next_vertex_ptr;
        }

      /* If no pairs of vertices combined this time round, then can
         finish here */
      if (pair_combined == FALSE ||  ndel < nstop)
        done = TRUE;
    }

  /* Adjust the final number of vertices by the number that were
     deleted */
  nvertices = nvertices - ndeleted;

  /* Return number of vertices retained */
  return(nvertices);
}
/***********************************************************************

create_vertex_array -  Create the arrays holding the x-, y- and z-
                       sort-indices

***********************************************************************/

void create_vertex_array(struct vertex *first_vertex_ptr,
                         struct vertex ***vertex_index_ptr,int nvertices)
{
  int ivertex;
  struct vertex **vertex_ind_ptr;
  struct vertex *vertex_ptr;

  /* Create the index array */
  vertex_ind_ptr
          = (struct vertex **) malloc(nvertices * sizeof(struct vertex *));

  /* Check that sufficient memory was allocated */
  if (vertex_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for vertex ");
      printf("sort-index array\n");
      exit(1);
    }

  /* Get pointer to the first vertex */
  ivertex = 0;
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to get current index array to point to
     each vertex in turn */
  while (vertex_ptr != NULL)
    {
      /* Include only non-deleted vertices */
      if (vertex_ptr->deleted == FALSE)
        {
          /* Initialise the index array */
          vertex_ind_ptr[ivertex] = vertex_ptr;

          /* Increment count of stiored vertices */
          ivertex++;
        }

      /* Get pointer to the next vertex in the linked list */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Return the array pointer */
  *vertex_index_ptr = vertex_ind_ptr;
}
/***********************************************************************

create_vertex_kdtree -  Create the kd-tree for the vertices

***********************************************************************/

int create_vertex_kdtree(struct vertex *first_vertex_ptr,int nvertices,
                         KdTree **kdt_ptr,
                         struct vertex **vertex_index_ptr)
{
  int ivertex;

  double **coords_array;

  struct vertex *vertex_ptr;
  KdTree *kdtree_ptr;

  /* Create the array holding the x-, y- and z-values */
  coords_array = (double**) calloc(nvertices,sizeof(double));

  /* Check that sufficient memory was allocated */
  if (coords_array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for coords_array\n");
      exit(1);
    }

  /* Get pointer to the first vertex */
  ivertex = 0;
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to place each in the kd-tree */
  while (vertex_ptr != NULL)
    {
      /* Only include vertices that haven't already been deleted */
      if (vertex_ptr->deleted == FALSE)
        {
          /* Store the vertex pointer in the vertex index array */
          vertex_index_ptr[ivertex] = vertex_ptr;

          /* Store the coords of this point in the array */
          coords_array[ivertex] = (double *) vertex_ptr->coords;

          /* Increment vertex count */
          ivertex++;
        }

      /* Get pointer to the next vertex in the linked list */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Update number of vertices */
  nvertices = ivertex;

  /* Create the kd-tree */
  kdtree_ptr = KdTree_create(coords_array,nvertices,3);

  /* Return pointer to start of kd-tree */
  *kdt_ptr = kdtree_ptr;

  /* Return current number of vertices */
  return(nvertices);
}
/***********************************************************************

filter_vertices  -  Filter out vertices that are too close to one
                    another (replacing each pair of too-close points
                    by another, midway between them)

***********************************************************************/

int filter_vertices(struct vertex *first_vertex_ptr,int nvertices,
                    struct vertex **vertex_index_ptr,
                    int buried_vertices_only,int save_vertices)
{
  int ivertex, jvert, jvertex, ndeleted;
  int ndel, nloops, nstop;
  int done, pair_combined;

  double inner_radius, outer_radius;
  double dist_sqrd, x, y, z, other_x, other_y, other_z, too_close_sqrd;

  struct vertex *vertex_ptr, *other_vertex_ptr;

  KdTree *kdtree_ptr;
  KdTreeQuery *kdquery_ptr;
  Region *region_ptr;

  /* Initialise variables */
  done = FALSE;
  inner_radius = 0.0;
  outer_radius = TOO_CLOSE_DIST;
  if (buried_vertices_only == TRUE || save_vertices == TRUE)
    outer_radius = TOO_CLOSE_BV;
  too_close_sqrd = outer_radius * outer_radius;
  ndeleted = 0;
  nloops = 0;

  /* Calculate a reasonable time to stop the filtering as being after
     any loop during which fewer than 1% of the vertices were deleted */
  nstop = (int) ((float) 0.1 * nvertices);

  /* Loop until vertices have been sufficiently filtered */
  while (done == FALSE)
    {
      /* Increment loop count */
      nloops++;

      /* Initialise count of vertices deleted in this loop */
      ndel = 0;

      /* Initialise flag saying whether any vertices combined in
         current loop */
      pair_combined = FALSE;

      /* Generate the kd-tree for the vertices */
      nvertices
        = create_vertex_kdtree(first_vertex_ptr,nvertices,&kdtree_ptr,
                               vertex_index_ptr);

      /* Loop through all vertices */
      for (ivertex = 0; ivertex < nvertices; ivertex++)
        {
          /* Get the current vertex */
          vertex_ptr = vertex_index_ptr[ivertex];

          /* Process only if this vertex hasn't already been deleted */
          if (vertex_ptr->deleted == FALSE)
            {
              /* Get the coordinates of this vertex */
              x = vertex_ptr->coords[0];
              y = vertex_ptr->coords[1];
              z = vertex_ptr->coords[2];

              /* Create the region of space to be searched for nearby
                 vertices */
              region_ptr
                = Annulus_create(vertex_ptr->coords,inner_radius,
                                 outer_radius,3);

              /* Perform the search */
              kdquery_ptr = KdTree_query(kdtree_ptr,region_ptr);

              /* Loop through all the vertices matched by the search */
              jvert = 0;
              while ((jvertex = KdTreeQuery_next(kdquery_ptr)) >= 0)
                {
                  /* Get the vertex pointer returned by the query */
                  other_vertex_ptr = vertex_index_ptr[jvertex];

                  /* If this is not the same vertex as we're currently
                     considering, calculate the distance between them */
                  if (jvertex != ivertex &&
                      other_vertex_ptr->deleted == FALSE)
                    {
                      /* Get the coordinates of this vertex */
                      other_x = other_vertex_ptr->coords[0];
                      other_y = other_vertex_ptr->coords[1];
                      other_z = other_vertex_ptr->coords[2];

                      /* Calculate the distance between the two vertices */
                      dist_sqrd = (x - other_x) * (x - other_x)
                        + (y - other_y) * (y - other_y)
                        + (z - other_z) * (z - other_z);

                      /* If the distance is less than the cut-off, then
                         want to replace the two vertices with a single one
                         half-way between them */
                      if (dist_sqrd < too_close_sqrd)
                        {
                          /* Combine the vertex coords to shift the current
                             vertex to halfway */
                          vertex_ptr->coords[0] = (x + other_x) / (float) 2.0;
                          vertex_ptr->coords[1] = (y + other_y) / (float) 2.0;
                          vertex_ptr->coords[2] = (z + other_z) / (float) 2.0;

                          /* If current vertex has no associated atom,
                             then transfer atom from vertex that is about
                             to be deleted */
                          if (vertex_ptr->atom_ptr == NULL)
                            vertex_ptr->atom_ptr = other_vertex_ptr->atom_ptr;

                          /* Delete the other vertex */
                          other_vertex_ptr->deleted = TRUE;

                          /* Reduce count of vertices */
                          ndel++;
                          ndeleted++;

                          /* Set flag that vertices combined in this pass */
                          pair_combined = TRUE;
                        }
                    }

                  /* Increment count of vertices matching the query */
                  jvert++;
                }

              /* Free up the memory taken for defining the region and
                 making the query */
              Annulus_free(region_ptr);
              KdTreeQuery_free(kdquery_ptr);
            }
        }

      /* Adjust total number of vertices */
      nvertices = nvertices - ndel;

      /* Free up the memory taken for the kd-tree and the vertex index
         array */
      KdTree_free(kdtree_ptr);

      /* If no pairs of vertices combined this time round, then can
         finish here */
      if (pair_combined == FALSE || ndel < nstop)
        done = TRUE;
    }

  /* Return the final number of filtered vertices */
  return(nvertices);
}
/***********************************************************************

write_vertices  -  Write out the vertices to the PDB file for display
                   by RasMol

***********************************************************************/

int write_vertices(FILE *outfile_ptr,int cluster,int *atom_no,
                   struct vertex **vert_atom_index_ptr,int *index,
                   int nvertices,int *nres,struct parameters params)
{
  char atom_name[5], chain, res_name[4], res_num[6];

  int i, ipoint, ivertex, j, ncentroid, npts;
  int atom_number, ires, nout;
  int iband, in_band;


  double bvalue, centroid[2][3], occupancy, ratio, totweight, weight;
  double bvalue_out;

  struct vertex *vertex_ptr;

  /* Initialise variables */
  atom_number = *atom_no;
  ncentroid = 0;
  ires = cluster;
  nout = 0;
  bvalue = 10.0;
  occupancy = 1.0;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 3; j++)
      centroid[i][j] = 0.0;
  totweight = 0.0;
  npts = 0;
  for (iband = 0; iband < NCONSERV; iband++)
    nres[iband] = 0;

  /* Write out a REMARK record to the PDB file showing which cleft
     region this is */
  fprintf(outfile_ptr,"REMARK  Cleft number %5d\n",cluster + 1);

  /* Loop through all the vertices */
  for (ivertex = 0; ivertex < nvertices; ivertex++)
    {
      /* Get pointer to this vertex */
      ipoint = index[ivertex];
      vertex_ptr = vert_atom_index_ptr[ipoint];

      /* Get the buried ratio, marking non-exposed points as non-zero */
      ratio = vertex_ptr->ratio;
      if (ratio == 0.0 && vertex_ptr->exposed == FALSE)
        ratio = 0.01;

      /* Write out vertex, checking for buried-only option */
      if (params.buried_vertices_only == FALSE || 
          vertex_ptr->ratio > 0.01)
        {
          /* Increment count of vertices written out */
          nout++;

          /* If have an associated atom, then use its name */
          if (params.buried_vertices_only == TRUE &&
              vertex_ptr->atom_ptr != NULL)
            {
              /* Transfer atom details across */
              strcpy(atom_name,vertex_ptr->atom_ptr->atom_name);
              strcpy(res_name,vertex_ptr->atom_ptr->residue_ptr->res_name);
              strcpy(res_num,vertex_ptr->atom_ptr->residue_ptr->res_num);
              chain = vertex_ptr->atom_ptr->residue_ptr->chain;
              if (!strcmp(atom_name," DUM"))
                {
                  strcpy(atom_name," C  ");
                  strcpy(res_name,"VTX");
                  strcpy(res_num,"   1 ");
                  chain = 'V';
                }
              occupancy = 1.0;
            }

          /* Otherwise, use defaults */
          else
            {
              /* Define dummy atom of residue VTX */
              strcpy(atom_name," C  ");
              strcpy(res_name,"VTX");
              chain = 'V';
              occupancy = 1.0;

              /* Form the residue number */
              sprintf(res_num,"%4d ",ires);
            }

          /* Get the atom-number, checking for overflow */
          atom_number++;
          if (atom_number > 99999)
            atom_number = 99999;

          /* If standard colouring, store the exposure ratio in the
             B-value column */
          bvalue = 0.0;
          if (params.vertex_colouring == STANDARD)
            bvalue = 0.8 * ratio;
          else if (params.vertex_colouring == BURIED_EXPOSED)
            {
              if (vertex_ptr->exposed == TRUE)
                bvalue = 10.0;
              else
                bvalue = 0.0;
            }

          /* Otherwise, get the conservation score from the corresponding
             residue, if there is one */
          else
            {
              /* If have read in the residue conservation from the -read
                 file, then use this */
              if (params.have_rescons == TRUE)
                bvalue = vertex_ptr->bvalue;

              /* Otherwise, get it directly from the residue */
              else if (vertex_ptr->atom_ptr != NULL)
                bvalue
                  = vertex_ptr->atom_ptr->residue_ptr->conservation_score;
            }

          /* Store the B-value used */
          vertex_ptr->bvalue = bvalue;

          /* If colouring by conservation, and result to be viewed by
             JMol, need to modify residue name according to conservation */
          if (params.vertex_colouring == CONSERVATION &&
              params.jmol_resnames == TRUE)
            {
              /* Initialise conservation range */
              in_band = -1;

              /* Determine which range this B-value lies in */
              for (iband = 1; iband < NCONSERV + 1 && in_band == -1; iband++)
                if (bvalue < RANGE[iband])
                  in_band = iband;

              /* Form residue name */
              strcpy(res_name,"VTA");
              if (in_band > 1)
                res_name[2] = (char) ((int) 'A' + in_band - 1);

              /* Update counts of residues in each band */
              if (in_band > 0 && in_band < NCONSERV + 1)
                nres[in_band - 1]++;
            }

          /* Determine the B-value to write out */
          bvalue_out = bvalue;
          if (params.jmol_resnames == FALSE)
            bvalue_out = 100 * bvalue;

          /* Write out the vertex record */
          fprintf(outfile_ptr,"ATOM  %5d ",atom_number);
          fprintf(outfile_ptr,"%4s ",atom_name);
          fprintf(outfile_ptr,"%3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  res_name,chain,res_num,vertex_ptr->coords[0],
                  vertex_ptr->coords[1],vertex_ptr->coords[2],occupancy,
                  bvalue_out);

          /* Add to computation of centroid(s) */
          if (params.vertex_colouring == CONSERVATION && bvalue > 0.0)
            {
              /* Use the residue conservation to weight the current
                 vertex coord */
              weight = bvalue * bvalue * bvalue;

              /* Add to centroid coords */
              for (i = 0; i < 3; i++)
                centroid[1][i]
                  = centroid[1][i] + weight * vertex_ptr->coords[i];
              totweight = totweight + weight;
            }

          /* If computing position of cluster centroid, add to coords */
          if (params.centroid == TRUE)
            {
              /* Add to centroid coords */
              for (i = 0; i < 3; i++)
                centroid[0][i] = centroid[0][i] + vertex_ptr->coords[i];

              /* Increment count of points used */
              npts++;
            }
        }
    }

  /* If showing centroid, write it out */
  if (params.centroid == TRUE)
    {
      /* Loop twice for each centroid */
      for (i = 0; i < 2; i++)
        {
          /* If can write out this centroid, do so */
          if ((i == 0 && npts > 0) || (i == 1 && totweight > 0.0))
            {
              /* Calculate the averaged coords of the centroid position */
              for (j = 0; j < 3; j++)
                {
                  if (i == 0)
                    centroid[i][j] = centroid[i][j] / npts;
                  else
                    centroid[i][j] = centroid[i][j] / totweight;
                }

              /* Get the atom-number, checking for overflow */
              atom_number++;
              if (atom_number > 99999)
                atom_number = 99999;

              /* Define dummy atom of residue CNT */
              strcpy(atom_name," C  ");
              strcpy(res_name,"CNT");
              chain = 'V';
              occupancy = 1.0;
              bvalue = 1.0;

              /* Form the residue number */
              sprintf(res_num,"%4d ",ires);

              /* Write out the centroid record */
              fprintf(outfile_ptr,"ATOM  %5d ",atom_number);
              fprintf(outfile_ptr,"%4s ",atom_name);
              fprintf(outfile_ptr,
                      "%3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                      res_name,chain,res_num,centroid[i][0],centroid[i][1],
                      centroid[i][2],occupancy,bvalue);

              /* Increment count of centroids */
              ncentroid++;
            }
        }
    }

  /* Write a TER record out */
  fprintf(outfile_ptr,"TER   \n");

  /* Return current atom-number */
  *atom_no = atom_number;

  /* Return number of cluster centroids written out */
  return(ncentroid);
}
/***********************************************************************

get_rasmol_colour  -  Retrieve the corresponding RasMol colour triplet
                      from the supplied colour name

***********************************************************************/

void get_rasmol_colour(char *colour_name,char *colour_numbers)
{
  char number_string[20];
  
  int i, icol;
  
  double rgb[3];
  
  /* Initialise variables */
  strcpy(colour_numbers,"[");
  
  /* Get the RGB equivalent of the colour */
  get_colour_rgb(colour_name,rgb);
  
  /* Convert RGB values to RasMol colour string */
  for (i = 0; i < 3; i++)
    {
      /* Convert double to integer in range 0-255 */
      icol = (int) (rgb[i] * 255);

      /* Add to colour numbers string */
      sprintf(number_string,"%d",icol);
      strcat(colour_numbers,number_string);

      /* Add comma or close-brackets */
      if (i == 2)
        strcat(colour_numbers,"]");
      else
        strcat(colour_numbers,",");
    }
}
/***********************************************************************

write_restype_colours  -  Colour the binding site vertices according to
                          residue types

***********************************************************************/

void write_restype_colours(FILE *file_ptr,char *object_type)
{
#define NTYPES      7

  char colour_numbers[COLNO_LEN];

  int iamino, itype, nout;

  static char amino[]={"ACDEFGHIKLMNPQRSTVWY"};
  
  static char *sas_colour_name[] = { 
    "grey", "yellow", "red", "red", "purple", "brown", "skyblue", "grey",
    "skyblue", "grey", "grey", "green", "brown", "green", "skyblue",
    "green", "green", "grey", "purple", "purple", };

  static char *aa_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
    "SER", "THR", "VAL", "TRP", "TYR" };

  static char *colour_name[] = { 
    "skyblue", "red", "green", "grey", "purple", "brown", "yellow" };

  static char *colour_desc[] = { 
    "positive", "negative", "neutral", "aliphatic", "aromatic",
    "Pro&Gly", "cysteine" };

  static int res_type[20]
    = { 3, 6, 1, 1, 4, 5, 0, 3, 0, 3, 3, 2, 5, 2, 0, 2, 2, 3, 4, 4 };

  
  /* Loop over the different residue types */
  for (itype = 0; itype < NTYPES; itype++)
    {
      /* Start output line for this residue type */
      fprintf(file_ptr,"select %s & (",object_type);

      /* Initialise count of amino acids of this type */
      nout = 0;

      /* Loop through all the amino acids to pick out those belonging
         to this residue type */
      for (iamino = 0; iamino < 20; iamino++)
        {
          /* If this residue belongs to current type, then write out
             name */
          if (res_type[iamino] == itype)
            {
              /* If not the first one, add an OR */
              if (nout > 0)
                fprintf(file_ptr," | ");

              /* Write out this amino acid */
              fprintf(file_ptr,"%s ",aa_name[iamino]);

              /* Increment count of amino acids of this type */
              nout++;
            }
        }

      /* Close off this residue type's line */
      fprintf(file_ptr,")\n");

      /* Get the RasMol colour-number string for this colour */
      get_rasmol_colour(colour_name[itype],colour_numbers);

      /* Write out its colour */
      fprintf(file_ptr,"colour %s\n",colour_numbers);
    }
}
/***********************************************************************

write_restype_pymol  -  Colour the binding site vertices according to
                        residue types

***********************************************************************/

void write_restype_pymol(FILE *file_ptr,char *object_type)
{
#define NTYPES      7

  char colour_numbers[COLNO_LEN];

  int iamino, itype, nout;

  static char amino[]={"ACDEFGHIKLMNPQRSTVWY"};
  
  static char *aa_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
    "SER", "THR", "VAL", "TRP", "TYR" };

  static char *colour_name[] = { 
    "skyblue", "salmon", "green", "grey", "purple", "brown", "yellow" };

  static int res_type[20]
    = { 3, 6, 1, 1, 4, 5, 0, 3, 0, 3, 3, 2, 5, 2, 0, 2, 2, 3, 4, 4 };

  
  /* Loop over the different residue types */
  for (itype = 0; itype < NTYPES; itype++)
    {
      /* Start output line for this residue type */
      fprintf(file_ptr,"color %s, %s and (",colour_name[itype],
              object_type);

      /* Initialise count of amino acids of this type */
      nout = 0;

      /* Loop through all the amino acids to pick out those belonging
         to this residue type */
      for (iamino = 0; iamino < 20; iamino++)
        {
          /* If this residue belongs to current type, then write out
             name */
          if (res_type[iamino] == itype)
            {
              /* If not the first one, add an OR */
              if (nout > 0)
                fprintf(file_ptr," or ");

              /* Write out this amino acid */
              fprintf(file_ptr,"resn %s",aa_name[iamino]);

              /* Increment count of amino acids of this type */
              nout++;
            }
        }

      /* Close off this residue type's line */
      fprintf(file_ptr,")\n");
    }
}
/***********************************************************************

write_conserv_pymol  -  Colour the binding site vertices according to
                        residue conservation

***********************************************************************/

void write_conserv_pymol(FILE *file_ptr,char *object_type)
{
  /* Colour according to conservation */
  fprintf(file_ptr,"spectrum B, rainbow, %s, minimum=0, maximum=100\n",
          object_type);
}
/***********************************************************************

write_conserv_colours  -  Colour the binding site vertices according to
                          residue conservation

***********************************************************************/

void write_conserv_colours(FILE *file_ptr,char *object_type,
                           int jmol_resnames,int *nres)
{
  char op[3], colour_numbers[COLNO_LEN], res_name[4];

  int iband, temp1, temp2;

  static char *colour_name[] = { 
    "blue", "purple", "skyblue", "cyan", "grey", "green", "yellow",
    "orange", "pink", "red" };
  
  /* Loop over the different conservation bands */
  for (iband = 0; iband < NCONSERV; iband++)
    {
      /* Write out this band only if have any residues in it -
         otherwise JMol just stops the script(!!) */
      if (jmol_resnames == FALSE || nres[iband] > 0)
        {
          /* Get the range of temperature values */
          temp1 = (int) (10000 * RANGE[iband]);
          temp2 = (int) (10000 * RANGE[iband + 1]);

          /* Form equivalent JMol residue name */
          strcpy(res_name,"VTA");
          res_name[2] = (char) ((int) 'A' + iband);

          /* Select the residues by temperature-factor range */
          if (iband == 0)
            strcpy(op,">=");
          else
            strcpy(op,">");
          if (jmol_resnames == TRUE)
            fprintf(file_ptr,"select %s & %s\n",object_type,res_name);
          else
            fprintf(file_ptr,
                    "select %s & temperature %s %d & temperature <= %d\n",
                    object_type,op,temp1,temp2);

          /* Get the RasMol colour-number string for this colour */
          get_rasmol_colour(colour_name[iband],colour_numbers);

          /* Write out the colour for this conservation band */
          fprintf(file_ptr,"colour %s\n",colour_numbers);
        }
    }
}
/***********************************************************************

write_script_clusters  -  Write out the script commands to define the
                          cluster surface vertices

***********************************************************************/

void write_script_clusters(FILE *outfile_ptr,FILE *pmlfile_ptr,
                           int icluster,int first_atom_number,
                           int last_atom_number,double volume,
                           int ncentroids,int *nres,
                           struct parameters params)
{
#define NRANGES       7

  char object_type[10];
  char colour_numbers[COLNO_LEN];

  static char *RANGE_COLOUR[NRANGES] = {
    "purple", "blue", "cyan", "green", "yellow", "orange", "red" };

  int cluster_num, colour, i, diff, irange, value;
  int have_pymol;

  static char *col_name[] = {
    "red", "purple", "yellow", "blue", "green",
    "brown", "pink", "olivegreen", "magenta", "cyan",
    "limegreen", "violet", "orange", "white"
  };

  /* Initialise variables */
  colour = icluster % 13;
  cluster_num = icluster + 1;
  sprintf(object_type,"gap%d",cluster_num);
  value = 0;
  diff = (int) (8000.0 / NRANGES);
  have_pymol = FALSE;
  if (pmlfile_ptr != NULL)
    have_pymol = TRUE;

  /* Print comment */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Cluster %d\n",cluster_num);
  if (have_pymol == TRUE)
    {
      fprintf(pmlfile_ptr,"\n");
      fprintf(pmlfile_ptr,"# Cluster %d\n",cluster_num);
    }

  /* Define the atom-number range for the cluster */
  fprintf(outfile_ptr,"define %s atomno >= %d & atomno <= %d\n",
          object_type,first_atom_number,last_atom_number - ncentroids);
  if (have_pymol == TRUE)
    fprintf(pmlfile_ptr,"create %s, id %d-%d\n",
            object_type,first_atom_number,last_atom_number - ncentroids);

  /* Select and render */
  fprintf(outfile_ptr,"select %s\n",object_type);
  if (params.jmol_resnames == TRUE)
    fprintf(outfile_ptr,"wireframe 0.1\n");
  else
    fprintf(outfile_ptr,"wireframe on\n");
  if (have_pymol == TRUE)
    {
      if (params.pdb_type == PDBSUM1 || params.pymol_www == TRUE)
        fprintf(pmlfile_ptr,"show line, %s\n",object_type);
      else
        fprintf(pmlfile_ptr,"show surface, %s\n",object_type);
    }
  
  /* Colour according to user-defined scheme */
  if (params.vertex_colouring == STANDARD)
    {
      /* Get the RasMol colour-number string for this colour */
      get_rasmol_colour(col_name[colour],colour_numbers);

      /* Write out its colour */
      fprintf(outfile_ptr,"colour %s\n",colour_numbers);
      if (have_pymol == TRUE)
        fprintf(pmlfile_ptr,"color %s, %s\n",col_name[colour],
                object_type);
    }
  else if (params.vertex_colouring == CPK)
    {
      /* Colour by atom type */
      fprintf(outfile_ptr,"colour cpk\n");

      /* Make carbons dark grey */
      /*      fprintf(outfile_ptr,"select %s & carbon\n",object_type);
              fprintf(outfile_ptr,"colour darkgrey\n"); */
    }
  else if (params.vertex_colouring == RESIDUE_TYPE)
    {
      if (have_pymol == TRUE)
        write_restype_pymol(pmlfile_ptr,object_type);
      else
        write_restype_colours(outfile_ptr,object_type);
    }
  else if (params.vertex_colouring == CONSERVATION)
    {
      if (have_pymol == TRUE)
        write_conserv_pymol(pmlfile_ptr,object_type);
      else
        write_conserv_colours(outfile_ptr,object_type,
                              params.jmol_resnames,nres);
    }
  else if (params.vertex_colouring == BURIED_EXPOSED)
    fprintf(outfile_ptr,"colour temperature\n");

  /* If colouring by B-value, write the commands */
  else if (params.vertex_colouring == BY_BVALUE)
    {
      /* Loop over the residue conservation ranges to colour accordingly */
      for (irange = 0; irange < NRANGES; irange++)
        {
          /* Write out the colour for this range */
          fprintf(outfile_ptr,"select %s & temperature > %d\n",
                  object_type,value);
          fprintf(outfile_ptr,"colour %s\n",RANGE_COLOUR[irange]);

          /* Increment the value */
          value = value + diff;
        }
    }

  /* Show details in RasMol command window for user's info */
  fprintf(outfile_ptr,"echo \"   %s:   Volume =%8.2f A^3",
          object_type,volume);
  if (params.vertex_colouring == STANDARD)
    fprintf(outfile_ptr," (%s)",col_name[colour]);
  fprintf(outfile_ptr,"\"\n");

  /* If have any cluster centroids, show these */
  for (i = 0; i < ncentroids; i++)
    {
      /* Show first centroid */
      fprintf(outfile_ptr,"# Centroid %d\n",i + 1);

      /* Define atom corresponding to centroid */
      fprintf(outfile_ptr,"define centroid%d_%d atomno = %d\n",
              cluster_num,i + 1,last_atom_number - ncentroids + i + 1);

      /* Select and render */
      fprintf(outfile_ptr,"select centroid%d_%d\n",cluster_num,i + 1);
      fprintf(outfile_ptr,"spacefill 0.3\n");
  
      /* Colour according to which centroid this is */
      if (i == 0)
        fprintf(outfile_ptr,"colour white\n");
      else
        fprintf(outfile_ptr,"colour red\n");
    }
}
/***********************************************************************

write_all_vertices  -  Process the list of vertices just read in from
                       file to write them out cluster by cluster

***********************************************************************/

void write_all_vertices(FILE *outfile_ptr,FILE *scriptfile_ptr,
                        FILE *pmlfile_ptr,
                        struct vertex *first_vertex_ptr,int nvertices,
                        int nclusters,int atom_number,double *volume,
                        struct parameters params)
{
  int icluster, ivertex, nres[NCONSERV], nvert;
  int first_atom_number, last_atom_number, ncentroids;
  int done;

  int *index;

  struct vertex *vertex_ptr, *first_of_cluster_vertex_ptr;
  struct vertex *next_vertex_ptr, *last_of_cluster_vertex_ptr;

  struct vertex **vertex_index_ptr;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(nvertices * sizeof(int));
  vertex_index_ptr
    = (struct vertex **) malloc(nvertices * sizeof(struct vertex *));

  /* Loop over all the sphere clusters to be contoured */
  for (icluster = 0; icluster < nclusters; icluster++)
    {
      /* Initialise variables */
      ivertex = 0;
      done = FALSE;

      /* If this cluster is not wanted, set flag as done */
      if (params.show_gap[icluster] == FALSE)
        done = TRUE;

      /* Initialise first and last atom numbers for this cluster */
      first_atom_number = atom_number + 1;
      last_atom_number = atom_number;

      /* Initialise vertex pointers */
      first_of_cluster_vertex_ptr = NULL;
      last_of_cluster_vertex_ptr = NULL;
      next_vertex_ptr = NULL;

      /* Get pointer to the first vertex */
      vertex_ptr = first_vertex_ptr;

      /* Loop through all vertices */
      while (vertex_ptr != NULL && done == FALSE)
        {
          /* If vertex belongs to the current cluster the store in
             the array */
          if (vertex_ptr->cluster == icluster)
            {
              /* If total number of vertices not exceeded, then store
                 vertex pointer */
              if (ivertex < nvertices)
                {
                  /* Store all the details in the arrays */
                  index[ivertex] = ivertex;
                  vertex_index_ptr[ivertex] = vertex_ptr;

                  /* If this is the first of this cluster's vertices,
                     store pointer to it */
                  if (first_of_cluster_vertex_ptr == NULL)
                    first_of_cluster_vertex_ptr = vertex_ptr;

                  /* Store as cluster's last vertex */
                  last_of_cluster_vertex_ptr = vertex_ptr;
                }

              /* Increment vertex count */
              ivertex++;
            }

          /* Check if we have done the current cluster */
          else if (vertex_ptr->cluster > icluster)
            done = TRUE;

          /* Get pointer to the next vertex */
          vertex_ptr = vertex_ptr->next_vertex_ptr;
        }

      /* Store number of vertices */
      nvert = ivertex;
      if (nvert > nvertices)
        nvert = nvertices;

      /* If showing the whole gap region, then need to filter the 
         vertices */
      if (nvert > 1 && params.buried_vertices_only == FALSE &&
          params.filter_vertices == TRUE)
        {
          /* Temporarilty terminate cluster's list of vertices for
             filtering */
          if (last_of_cluster_vertex_ptr != NULL)
            {
              /* Save the vertex this points to */
              next_vertex_ptr
                = last_of_cluster_vertex_ptr->next_vertex_ptr;

              /* Terminate the linked list */
              last_of_cluster_vertex_ptr->next_vertex_ptr = NULL;
            }

          /* Filter the vertices */
          nvert = filter_vertices(first_of_cluster_vertex_ptr,nvert,
                                  vertex_index_ptr,
                                  params.buried_vertices_only,
                                  params.save_vertices);

          /* Reinstate the linked list terminated prior to the filter */
          if (last_of_cluster_vertex_ptr != NULL)
            last_of_cluster_vertex_ptr->next_vertex_ptr
              = next_vertex_ptr;
        }

      /* Write out all the vertices making up the current cluster */
      ncentroids
        = write_vertices(outfile_ptr,icluster,&atom_number,
                         vertex_index_ptr,index,nvert,nres,params);
      last_atom_number = atom_number;

      /* Write out the RasMol script commands to display the cluster
         jusr written out */
      if (last_atom_number > first_atom_number)
        write_script_clusters(scriptfile_ptr,pmlfile_ptr,icluster,
                              first_atom_number,last_atom_number,
                              volume[icluster],ncentroids,nres,params);
    }
}
/***********************************************************************

get_protein_name  -  Pick up the protein names from the PDBsum 
                     main.html files

***********************************************************************/

void get_protein_name(char *pdb_code,char *pdbsum_dir,char *protein_name)
{
  char file_name[FILENAME_LEN];
  char input_line[LINELEN + 1];
  char *string1_ptr, *string2_ptr;

  int got_title, have_file, next_line, wanted;
  int len, line;

  FILE *file_ptr;

  /* Initialise various flags  */
  got_title = FALSE;
  have_file = TRUE;
  line = 0;
  next_line = FALSE;
  protein_name[0] = '\0';

  /* Form the current filename */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"main.html");

  /* Open the PDBsum file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in the protein name from it */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             got_title == FALSE)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Check for end of title */
          if (next_line == TRUE && line > 0 && input_line[0] == '<')
            {
              /* Set flags to end here */
              next_line = FALSE;
              got_title = TRUE;
            }

          /* If this is the title, then write out */
          else if (next_line == TRUE && strncmp(input_line,"   ",3))
            {
              /* Check that line isn't one of the unwanted ones */
              wanted = TRUE;
              if (!strncmp(input_line,"Chain:",6) ||
                  !strncmp(input_line," Chain:",7) ||
                  !strncmp(input_line,"Engineered:",11) ||
                  !strncmp(input_line," Engineered:",12))
                wanted = FALSE;

              /* Store this line */
              if (wanted == TRUE)
                {
                  /* Get line length */
                  if (len > 0 && strlen(protein_name) + len < NAMLEN)
                    {
                      /* Append to name */
                      strcat(protein_name,input_line);

                      /* Increment line-count */
                      line++;
                    }
                }
            }

          /* If this is the line before the protein name, set flag */
          if (got_title == FALSE)
            {
              /* Check for the keywords */
              string1_ptr = strstr(input_line,"<b>Structure:</b>");
              string2_ptr = strstr(input_line,"<b>Title:</b>");

              /* If this is the line before the protein name, set flag */
              if (string1_ptr != NULL || string2_ptr != NULL)
                next_line = TRUE;
            }                  
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show warning message */
  else
    {
      printf("*** Warning. PDBsum main.html file not found:\n");
      printf("***          %s\n",file_name);
    }
}
/***********************************************************************

profunc_html_header  -  Write the ProFunc header lines to the HTML file

***********************************************************************/

void profunc_html_header(FILE *file_ptr)
{
  /* Write out the header lines */
  fprintf(file_ptr,"<HTML>\n");
  fprintf(file_ptr,"<HEAD>\n");
  fprintf(file_ptr,"<TITLE>Binding site analysis</TITLE>\n");
  fprintf(file_ptr,"</HEAD>\n");
  // fprintf(file_ptr,"<BODY BGCOLOR=\"#008ED5\" text=\"BLACK\"");
  fprintf(file_ptr,"<BODY BGCOLOR=\"WHITE\" text=\"BLACK\"");
  fprintf(file_ptr," link=\"YELLOW\" vlink=\"ORANGE\" alink=\"RED\">\n");
  fprintf(file_ptr,"\n");

  /* Table containing top bar */
  fprintf(file_ptr,"<TABLE WIDTH=780>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH ALIGN=LEFT>\n");
  fprintf(file_ptr,"   <A HREF=\"*WWW home*\"\n");
  fprintf(file_ptr,"onMouseOver=\"window.status='Go to ProFunc");
  fprintf(file_ptr," home page';");
  fprintf(file_ptr," return true;\">\n");
  fprintf(file_ptr,"   <IMG BORDER=0 SRC=\"gif/logo.gif\"></A>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TD ALIGN=\"RIGHT\">\n");
  fprintf(file_ptr,"   <TABLE BORDER=1>\n");
  fprintf(file_ptr,"    <TR>\n");
  //  fprintf(file_ptr,"     <TH BGCOLOR=\"#008ED5\">\n");
  fprintf(file_ptr,"     <TH BGCOLOR=\"WHITE\">\n");
  // fprintf(file_ptr,"      <FONT COLOR=WHITE SIZE=+2>&nbsp;Binding sites");
  fprintf(file_ptr,"      <FONT COLOR=#888888 SIZE=+2>&nbsp;Binding sites");
  fprintf(file_ptr,"</FONT>\n");
  fprintf(file_ptr,"     </TH>\n");
  fprintf(file_ptr,"    </TR>\n");
  fprintf(file_ptr,"   </TABLE>\n");
  fprintf(file_ptr,"  </TD>\n");
  fprintf(file_ptr," </TR>\n");

  /* Back button */
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TD COLSPAN=2>\n");
  fprintf(file_ptr,"   <A HREF=\"/cgi-bin/profunc/GetResults.pl?");
  fprintf(file_ptr,"*ID*\"\n");
  fprintf(file_ptr,"     onMouseOver=\"window.status='Back to menu ");
  fprintf(file_ptr,"of analyses'; return true;\">\n");
  fprintf(file_ptr,"   <IMG BORDER=0 ALT=\"Back\" SRC=\"gif/back.gif\">");
  fprintf(file_ptr,"</A>\n");
  fprintf(file_ptr,"  </TD>\n");
  fprintf(file_ptr," </TR>\n");

  /* Page title */
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH COLSPAN=2>\n");
  // fprintf(file_ptr,"   <H1><FONT COLOR=WHITE>Binding site analysis");
  fprintf(file_ptr,"   <H1><FONT COLOR=#888888>Binding site analysis");
  fprintf(file_ptr,"</FONT></H1>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"</TABLE>\n");
}
/***********************************************************************

write_html_headers  -  Write the header records to the HTML file

***********************************************************************/

void write_html_headers(FILE *file_ptr,char *pdb_code,
                        char chain[MAX_CHAINS],int nchain,
                        char *protein_name,char *www_pdbsum_dir,
                        struct parameters params)
{
  char add_chain[10], chain_id[2], chain_string[10 * MAX_CHAINS];
  char file_name[FILENAME_LEN];

  int ichain;

  static char f_start[] = "<FONT COLOR=\"BLACK\">";
  static char f_end[] = "</FONT>";

  /* If called from ProFunc, display ProFunc banner */
  if (params.profunc == TRUE)
    profunc_html_header(file_ptr);

  /* Otherwise, show standard page */
  else
    {
      /* Write out header records */
      fprintf(file_ptr,"<HTML>\n");
      fprintf(file_ptr,"\n");
      fprintf(file_ptr,"<HEAD>\n");
      fprintf(file_ptr,"<TITLE>Binding site analysis");
      if (pdb_code[0] != '\0')
        fprintf(file_ptr," for %s",pdb_code);

      /* If have several chains, then list these */
      if (nchain == 1 && chain[0] != ' ')
        fprintf(file_ptr," chain %c",chain[0]);
      else if (nchain > 1)
        {
          /* Prepare to list chains */
          fprintf(file_ptr," chains");

          /* Loop through all the chains to write out */
          for (ichain = 0; ichain < nchain; ichain++)
            fprintf(file_ptr," %c",chain[ichain]);
        }
      fprintf(file_ptr,"</TITLE>\n");
      fprintf(file_ptr,"</HEAD>\n");
      fprintf(file_ptr,"\n");

      /* Set page's background colour */
      fprintf(file_ptr,"<BODY BGCOLOR=\"#9AC0CD\">\n");
      fprintf(file_ptr,"\n");

      /* Show headings */
      fprintf(file_ptr,"<CENTER>\n");
      fprintf(file_ptr,"<H1>Binding site analysis");
      if (pdb_code[0] != '\0')
        fprintf(file_ptr," for %s",pdb_code);

      /* If have several chains, then list these */
      if (nchain == 1 && chain[0] != ' ')
        {
          chain_id[0] = chain[0];
          chain_id[1] = '\0';
          fprintf(file_ptr," chain %c",chain[0]);
        }
      else if (nchain > 1)
        {
          /* Prepare to list chains */
          chain_id[0] = '\0';
          fprintf(file_ptr," chains");

          /* Loop through all the chains to write out */
          for (ichain = 0; ichain < nchain; ichain++)
            fprintf(file_ptr," %c",chain[ichain]);
        }
      else
        chain_id[0] = '\0';
      fprintf(file_ptr,"</H1>\n");
      fprintf(file_ptr,"\n");
    }

  /* Show figure */
  if (params.no_gif == FALSE)
    {
      if (params.profunc == TRUE)
        {
          fprintf(file_ptr,"<CENTER>\n");
          fprintf(file_ptr,"<TABLE WIDTH=780>\n");
        }
      else
        fprintf(file_ptr,"<TABLE>\n");
      fprintf(file_ptr,"<TR>\n");
      fprintf(file_ptr,"<TH>\n");
      if (params.profunc == TRUE)
        fprintf(file_ptr,"<IMG SRC=\"tmpgif/speedfill.jpg\">\n");
      // fprintf(file_ptr,"<IMG SRC=\"tmpgif/speedfill.gif\">\n");
      else
        fprintf(file_ptr,"<IMG SRC=\"bindsite_%s.gif\">\n",chain_id);
      fprintf(file_ptr,"</TH>\n");
      fprintf(file_ptr,"</TR>\n");
      fprintf(file_ptr,"<TD ALIGN=CENTER>\n");
      if (params.profunc == TRUE)
        fprintf(file_ptr,"<FONT COLOR=#888888>\n");
      // fprintf(file_ptr,"<FONT COLOR=WHITE>\n");
      fprintf(file_ptr,"<EM>Clefts and cavities in protein structure");
      if (pdb_code[0] != '\0')
        fprintf(file_ptr," <B>%s</B>",pdb_code);
      fprintf(file_ptr,",<BR>\n");
      fprintf(file_ptr,"the colours corresponding to the entries in the ");
      fprintf(file_ptr,"table below.\n");
      fprintf(file_ptr,"</EM>\n");
      if (params.profunc == TRUE)
        fprintf(file_ptr,"</FONT>\n");
      fprintf(file_ptr,"</TD>\n");
      fprintf(file_ptr,"</TR>\n");
      fprintf(file_ptr,"</TABLE>\n");
      if (params.profunc == TRUE)
        fprintf(file_ptr,"</CENTER>\n");
      fprintf(file_ptr,"\n");
    }
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"</CENTER>\n");

  /* Form link to PDBsum page */
  strcpy(file_name,www_pdbsum_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,pdb_code);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"main.html");

  /* Write out the link */
  fprintf(file_ptr,"<H3>");
  fprintf(file_ptr,"<A HREF=\"%s\">%s</A>&nbsp;",file_name,pdb_code);
  fprintf(file_ptr,"<EM><FONT COLOR=BLUE>");

  /* Write out the protein name */
  fprintf(file_ptr,"%s ",protein_name);
  fprintf(file_ptr,"</FONT></EM></H3>");

  /* Print explanatory text */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"<P>\n");
  if (params.profunc == TRUE)
    fprintf(file_ptr,"<FONT COLOR=#888888>\n");
  // fprintf(file_ptr,"<FONT COLOR=WHITE>\n");
  fprintf(file_ptr,"The following table shows the \"gap\" regions (clefts\n");
  fprintf(file_ptr,"and cavities) in the protein's surface.\n");
  fprintf(file_ptr,"The gaps are ordered by decreasing volume (in\n");
  fprintf(file_ptr,"&Aring;<SUP>3</SUP>), and the table includes various\n");
  fprintf(file_ptr,"parameters that are indicative of binding sites.\n");
  fprintf(file_ptr,"(Ignore data shown in grey).\n");
  if (params.profunc == TRUE)
    fprintf(file_ptr,"</FONT>\n");
  fprintf(file_ptr,"<P>\n");

  /* Start the form */
  fprintf(file_ptr,"<FORM METHOD=\"POST\" ");
  fprintf(file_ptr,"ACTION=\"/cgi-bin/cath/RunSpeedfill.pl\">");
  fprintf(file_ptr,"\n");

  /* Store the PDB code and chain-id */
  fprintf(file_ptr,"<INPUT TYPE=\"hidden\" NAME=\"source\" ");
  fprintf(file_ptr,"VALUE=\"bsite\">\n");
  fprintf(file_ptr,"<INPUT TYPE=\"hidden\" NAME=\"pdb\" ");
  fprintf(file_ptr,"VALUE=\"%s\">",pdb_code);
  if (nchain > 0)
    {
      /* Initialise the string holding all the chain identifiers */
      chain_string[0] = '\0';

      /* Loop through all the chains to add each one in */
      for (ichain = 0; ichain < nchain; ichain++)
        {
          if (chain[ichain] != ' ')
            {
              sprintf(add_chain,"-chain %c ",chain[ichain]);
              strcat(chain_string,add_chain);
            }
        }

      /* Write out the chain string */
      if (chain_string[0] != '\0')
        {
          fprintf(file_ptr,"<INPUT TYPE=\"hidden\" NAME=\"chain\" ");
          fprintf(file_ptr,"VALUE=\"%s\">",chain_string);
        }
    }
  fprintf(file_ptr,"\n");

  /* Store the RasMol script flag */
  fprintf(file_ptr,"<INPUT TYPE=\"hidden\" NAME=\"script\" ");
  fprintf(file_ptr,"VALUE=\"YES\">");
  fprintf(file_ptr,"\n");

  /* Start the table */
  fprintf(file_ptr,"<CENTER>\n");
  fprintf(file_ptr,"<TABLE BORDER=1 CELLPADDING=5 ");
  fprintf(file_ptr,"BGCOLOR=\"#FDF5E6\">\n");

  /* Write out the column headings */
  fprintf(file_ptr,"<TR>\n");
  fprintf(file_ptr,"<TH>%sGap<BR>region%s</TH>\n",f_start,f_end);
  fprintf(file_ptr,"<TH>%sVolume%s</TH>\n",f_start,f_end);
  fprintf(file_ptr,"<TH>%sR1<BR>ratio%s</TH>\n",f_start,f_end);
  fprintf(file_ptr,"<TH BGCOLOR=\"#CCCCCC\" COLSPAN=2>%s",f_start);
  fprintf(file_ptr,"Accessible<BR>vertices%s</TH>\n",f_end);
  fprintf(file_ptr,"<TH BGCOLOR=\"#CCCCCC\" COLSPAN=2>%s",f_start);
  fprintf(file_ptr,"Buried<BR>vertices%s</TH>\n",f_end);
  fprintf(file_ptr,"<TH BGCOLOR=\"#CCCCCC\" COLSPAN=2>%s",f_start);
  fprintf(file_ptr,"Ave.<BR>depth%s</TH>\n",f_end);
  fprintf(file_ptr,"<TH>%sResidue<BR>type%s<SUP>",f_start,f_end);
  fprintf(file_ptr,"<FONT COLOR=BLUE>*</SUP></FONT></TH>\n");
  fprintf(file_ptr,"<TH>%sResidue<BR>conservation%s<SUP>",f_start,f_end);
  fprintf(file_ptr,"<FONT COLOR=BLUE>**</SUP></FONT></TH>\n");
  fprintf(file_ptr,"<TH ALIGN=LEFT>%sLigands%s</TH>\n",f_start,f_end);
  fprintf(file_ptr,"</TR>\n");
  fprintf(file_ptr,"\n");
}
/***********************************************************************

encode_grid_value  -  Encode the current grid value and mask-flag into
                      a single value

***********************************************************************/

int encode_grid_value(int grid_value,int masked)
{
  int encoded_value;

  /* Encode value and flag */
  encoded_value = 10 * grid_value + masked;

  /* Return the encoded value */
  return(encoded_value);
}
/***********************************************************************

reset_grid_array  -  Re-initialise the grid array so as to retain a
                     memory of the grid-points masked by the atoms in
                     the structure

***********************************************************************/

void reset_grid_array(int *grid,int ngrid_points,int contour_level)
{
  int ipoint;
  int current_value, masked;

  /* Loop over all the points in the grid */
  for (ipoint = 0; ipoint < ngrid_points; ipoint++)
    {
      /* Get the grid value and masked flag */
      current_value = get_grid_value(grid[ipoint],&masked,contour_level);

      /* Reset grid-value to zero while retaining the mask flag */
      grid[ipoint] = masked;
    }
}
/***********************************************************************

update_map  -  Use given sphere to update points on the grid-map via
               the Gaussian distribution

***********************************************************************/

void update_map(int *grid,int npoints[3],double valmin[3],double valmax[3],
                double grid_sep,double coords[3],double radius,
                int *exp_lookup_table,int nradii,int nvalues,
                int contour_level,struct parameters params)
{
  int from_point[3], icoord, ipoint, nearest_point[3], nsteps, to_point[3];
  int i, j, k;
  int current_value, masked, update_value;
  int ipos, iradius, ivalue;

  double dist_sqrd, dx_sqrd, dy_sqrd, dz_sqrd, x, y, z;
  double constant;
  /*  double grid_value; */

  /* Calculate the constant in the Gaussian function given the radius
     of the current sphere */
  constant = (- log(2.0) / (radius * radius));

  /* Determine the closest radius in the lookup-table */
  iradius = (int) ((radius - params.min_gap_radius) / RADSTEP + 0.5);
  if (iradius < 0)
    iradius = 0;
  if (iradius > nradii - 1)
    iradius = nradii - 1;

  /* Find the nearest grid-point to the given position, and calculate
     its coordinates */
  ipoint = get_nearest_grid_point(npoints,valmin,valmax,grid_sep,
                                  coords,nearest_point);

  /* Calculate number of grid-points either side of nearest that need to
     be updated */
  nsteps = (int) (radius / grid_sep + 1);

  /* Calculate start- and end-points in each direction */
  for (icoord = 0; icoord < 3; icoord++)
    {
      /* Start-point */
      from_point[icoord] = nearest_point[icoord] - nsteps;
      if (from_point[icoord] < 0)
        from_point[icoord] = 0;

      /* End-point */
      to_point[icoord] = nearest_point[icoord] + nsteps;
      if (to_point[icoord] > npoints[icoord] - 1)
        to_point[icoord] = npoints[icoord] - 1;
    }

  /* Loop through the points either side of nearest one in x-direction */
  for (i = from_point[0]; i < to_point[0]; i++)
    {
      /* Get the x-coordinate */
      x = valmin[0] + i * grid_sep;
      dx_sqrd = (x - coords[0]) * (x - coords[0]);

      /* Loop through the points either side of nearest one in y-direction */
      for (j = from_point[1]; j < to_point[1]; j++)
        {
          /* Get the y-coordinate */
          y = valmin[1] + j * grid_sep;
          dy_sqrd = (y - coords[1]) * (y - coords[1]);

          /* Loop through the points either side of nearest one in
             z-direction */
          for (k = from_point[2]; k < to_point[2]; k++)
            {
              /* Get the z-coordinate */
              z = valmin[2] + k * grid_sep;
              dz_sqrd = (z - coords[2]) * (z - coords[2]);

              /* Calculate the distance between the grid-point and the
                 sphere centre */
              dist_sqrd = dx_sqrd + dy_sqrd + dz_sqrd;

              /* Calculate position in table */
              ivalue = (int) (dist_sqrd / DISTEP + 0.5);
              if (ivalue > nvalues - 1)
                update_value = 0;
              else
                {
                  /* Calculate location in the table */
                  ipos = iradius * nvalues + ivalue;

                  /* Retrieve the exponential value from the lookup
                     table */
                  update_value = exp_lookup_table[ipos];
                }

              /* Old method 

              Calculate update value for current grid-point using
              a spherically-symmetric Gaussian function
              grid_value = 2 * CONTOUR_LEVEL * exp(constant * dist_sqrd);
              update_value = (int) (GRID_SCALE_FACTOR * grid_value);
              */

              /* Get the array location of this grid-point */
              ipoint = i * npoints[1] * npoints[2]
                + j * npoints[2] + k;

              /* Extract this grid-point's current value */
              current_value
                = get_grid_value(grid[ipoint],&masked,contour_level);

              /* If new value is higher than grid-point's current value,
                 then update */
              if (update_value > current_value)
                grid[ipoint] = encode_grid_value(update_value,masked);
            }
        }
    }
}
/***********************************************************************

calc_grid_limits  -  Calculate the limits of the non-zero region of the
                     grid

***********************************************************************/

void calc_grid_limits(int *grid,int npoints[3],
                      struct grid_limits *glimits,int contour_level)
{
  int grid_value, i, icoord, ipoint, j, k;
  int masked;

  /* Initialise variables */
  for (icoord = 0; icoord < 3; icoord++)
    {
      glimits->minpos[icoord] = npoints[icoord] - 1;
      glimits->maxpos[icoord] = 0;
    }

  /* Loop through the x points */
  for (i = 0; i < npoints[0]; i++)
    {
      /* Loop through the y points */
      for (j = 0; j < npoints[1]; j++)
        {
          /* Loop through the z points */
          for (k = 0; k < npoints[2]; k++)
            {
              /* Get the array location of this grid-point */
              ipoint = i * npoints[1] * npoints[2]
                + j * npoints[2] + k;

              /* Get the grid-value at the current point */
              grid_value
                = get_grid_value(grid[ipoint],&masked,contour_level);

              /* If grid value is non-zero, then update boundaries */
              if (grid_value > 0)
                {
                  if (i < glimits->minpos[0])
                    glimits->minpos[0] = i;
                  if (i > glimits->maxpos[0])
                    glimits->maxpos[0] = i;
                  if (j < glimits->minpos[1])
                    glimits->minpos[1] = j;
                  if (j > glimits->maxpos[1])
                    glimits->maxpos[1] = j;
                  if (k < glimits->minpos[2])
                    glimits->minpos[2] = k;
                  if (k > glimits->maxpos[2])
                    glimits->maxpos[2] = k;
                }
            }
        }
    }
  
  /* Extend limits by one grid point either side */
  for (icoord = 0; icoord < 3; icoord++)
    {
      if (glimits->minpos[icoord] > 0)
        glimits->minpos[icoord]--;
      if (glimits->maxpos[icoord] < npoints[icoord] - 1)
        glimits->maxpos[icoord]++;
    }
}
/***********************************************************************

place_vertex_in_sort_list  -  Find the location of the current vertex
                              in the list sorted by the given coordinate

***********************************************************************/

void place_vertex_in_sort_list(int icoord,int iband,
                               struct vertex *vertex_ptr,
                               struct vertex **first_sorted_vertex_ptr,
                               struct vertex **last_sorted_vertex_ptr,
                               struct vertex *sort_vertex_ptr[3])
{
  int found;

  struct vertex *prev_vertex_ptr, *other_vertex_ptr;

  /* Initialise flag */
  found = FALSE;

  /* Initialise pointers */
  prev_vertex_ptr = NULL;

  /* Get pointer to current start-point in the sorted list */
  other_vertex_ptr = sort_vertex_ptr[icoord];

  /* Loop through vertices in sorted list until find one with a
     higher coordinate value than the current vertex, or we hit
     the end of the list */
  while (other_vertex_ptr != NULL && found == FALSE)
    {
      /* If coordinate value is higher, then we have found the
         place where to splice in our vertex */
      if (other_vertex_ptr->coords[icoord] >= vertex_ptr->coords[icoord])
        {
          /* Set flag that we have found a suitable location to
             place our vertex */
          found = TRUE;

          /* Get our vertex to point to the next vertex in the
             sorted list */
          vertex_ptr->next_sorted_vertex_ptr[icoord]
            = other_vertex_ptr;
        }

      /* Otherwise, store this pointer */
      else
        prev_vertex_ptr = other_vertex_ptr;

      /* Get pointer to next vertex in this sorted list */
      other_vertex_ptr
        = other_vertex_ptr->next_sorted_vertex_ptr[icoord];
    }

  /* If there is an earlier vertex in the list, get it to
     point to our vertex */
  if (prev_vertex_ptr != NULL)
    {
      /* Store the link */
      prev_vertex_ptr->next_sorted_vertex_ptr[icoord]
        = vertex_ptr;
    }

  /* Otherwise, this must be the first entry in the list
     so save its pointer as first in this list */
  else
    {
      /* If we have no previous first, then this will also be the
         current last */
      if (first_sorted_vertex_ptr[iband] == NULL)
        last_sorted_vertex_ptr[iband] = vertex_ptr;

      first_sorted_vertex_ptr[iband] = vertex_ptr;

      /* Set the search-start for current coordinate to this vertex */
      sort_vertex_ptr[icoord] = vertex_ptr;
    }

  /* If not found, then must have been added to the end */
  if (found == FALSE)
    last_sorted_vertex_ptr[iband] = vertex_ptr;
      
}
/***********************************************************************

splice_bands_together  -  Splice the separate bands of sorted vertices
                          together into one long, sorted list

***********************************************************************/

void splice_bands_together(int icoord,int nbands,
                           struct vertex **first_vertex_sort_ptr,
                           struct vertex **last_vertex_sort_ptr,
                           struct vertex *first_sorted_vertex_ptr[3])
{
  int iband;

  struct vertex *last_vertex_ptr;

  /* Initialise variables */
  last_vertex_ptr = NULL;

  /* Loop through all the bands */
  for (iband = 0; iband < nbands; iband++)
    {
      /* If have a pointer in this band, then process */
      if (first_vertex_sort_ptr[iband] != NULL)
        {
          /* If this is the very first non-null pointer, then store
             as such */
          if (first_sorted_vertex_ptr[icoord] == NULL)
            first_sorted_vertex_ptr[icoord] = first_vertex_sort_ptr[iband];

          /* If have a pointer from a previous band, get it to point
             to here */
          if (last_vertex_ptr != NULL)
            last_vertex_ptr->next_sorted_vertex_ptr[icoord]
              = first_vertex_sort_ptr[iband];

          /* Store the last pointer of this band */
          last_vertex_ptr = last_vertex_sort_ptr[iband];
        }
    }
}
/***********************************************************************

calculate_gap_volume  -  Calculate the volume of the current gap region

***********************************************************************/

double calculate_gap_volume(int *grid,int npoints[3],double valmin[3],
                           double valmax[3],double grid_sep,
                           int contour_level,int ngrid_points)
{
  int ipoint;
  int grid_value, masked;

  double cell_volume, volume;

  /* Initialise variables */
  volume = 0.0;

  /* Calculate cell volume */
  cell_volume = grid_sep * grid_sep * grid_sep;

  /* Loop through the grid array */
  for (ipoint = 0; ipoint < ngrid_points; ipoint++)
    {
      /* Extract the current grid-value */
      grid_value = get_grid_value(grid[ipoint],&masked,contour_level);

      /* If grid-point is inside the gap region, determine what its
         contribution to the volume is */
      if (grid_value > contour_level)
        {
          /* Add cell-volume to overall volume */
          volume = volume + cell_volume;
        }
    }

  /* Return gap volume */
  return(volume);
}
/***********************************************************************

get_contour_points  -  Get all contouring points in the current map

***********************************************************************/

int get_contour_points(int *grid,int *grid_atom_map,int npoints[3],
                       struct atom **atom_index_ptr,int natoms,
                       double valmin[3],double valmax[3],double grid_sep,
                       struct vertex **fst_vertex_ptr,
                       struct vertex *first_sorted_vertex_ptr[3],
                       int contour_level,struct grid_limits glimits,
                       int icluster,struct parameters params)
{
  int i, iatom, icoord, idir, ip, ipoint, j, jpoint, k, nvertices;
  int istart, iend, jstart, jend, kstart, kend;
  int ipos, jpos, kpos, ni, nj, nk;
  int edge_no, grid_value, masked, other_value;
  int add_point[3][3];

  static int add_on[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  double x, y, z, other_x, other_y, other_z;
  double add_sep[3][3], cross[3];
  double diff1, diff2, fraction;

  struct atom *atom_ptr;
  struct vertex *first_vertex_ptr, *last_vertex_ptr;
  struct vertex *vertex_ptr, *sorted_vertex_ptr[3];
  struct vertex **first_vertex_xsort_ptr, **first_vertex_ysort_ptr;
  struct vertex **first_vertex_zsort_ptr;
  struct vertex **last_vertex_xsort_ptr, **last_vertex_ysort_ptr;
  struct vertex **last_vertex_zsort_ptr;

  /* Initialise variables */
  edge_no = 0;
  nvertices = 0;

  /* Precalculate commonly used values for the loops below */
  for (idir = 0; idir < 3; idir++)
    {
      /* Calculate the three relative grid-point locations */
      add_point[idir][0] = add_on[idir][0] * npoints[1] * npoints[2];
      add_point[idir][1] = add_on[idir][1] * npoints[2];
      add_point[idir][2] = add_on[idir][2];

      /* Calculate relative grid-separations */
      for (i = 0; i < 3; i++)
        add_sep[idir][i] = add_on[idir][i] * grid_sep;
    }

  /* Initialise pointers */
  first_vertex_ptr = NULL;
  last_vertex_ptr = NULL;
  for (icoord = 0; icoord < 3; icoord++)
    first_sorted_vertex_ptr[icoord] = NULL;

  /* Define range of non-zero array */
  istart = glimits.minpos[0];
  iend = glimits.maxpos[0] + 1;
  if (iend > npoints[0] - 1)
    iend = npoints[0] - 1;
  jstart = glimits.minpos[1];
  jend = glimits.maxpos[1] + 1;
  if (jend > npoints[1] - 1)
    jend = npoints[1] - 1;
  kstart = glimits.minpos[2];
  kend = glimits.maxpos[2] + 1;
  if (kend > npoints[2] - 1)
    kend = npoints[2] - 1;

  /* Allocate memory for start- and end pointers for vertex sort-lists
     in each subscript band */

  /* Bands in x-direction */
  ni = iend - istart;
  first_vertex_xsort_ptr
    = (struct vertex **) malloc(ni * sizeof(struct vertex *));
  last_vertex_xsort_ptr
    = (struct vertex **) malloc(ni * sizeof(struct vertex *));
  for (i = 0; i < ni; i++)
    {
      first_vertex_xsort_ptr[i] = NULL;
      last_vertex_xsort_ptr[i] = NULL;
    }

  /* Bands in y-direction */
  nj = jend - jstart;
  first_vertex_ysort_ptr
    = (struct vertex **) malloc(nj * sizeof(struct vertex *));
  last_vertex_ysort_ptr
    = (struct vertex **) malloc(nj * sizeof(struct vertex *));
  for (j = 0; j < nj; j++)
    {
      first_vertex_ysort_ptr[j] = NULL;
      last_vertex_ysort_ptr[j] = NULL;
    }

  /* Bands in z-direction */
  nk = kend - kstart;
  first_vertex_zsort_ptr
    = (struct vertex **) malloc(nk * sizeof(struct vertex *));
  last_vertex_zsort_ptr
    = (struct vertex **) malloc(nk * sizeof(struct vertex *));
  for (k = 0; k < nk; k++)
    {
      first_vertex_zsort_ptr[k] = NULL;
      last_vertex_zsort_ptr[k] = NULL;
    }

  /* Initialise i-position */
  ipos = 0;

  /* Loop through the x points */
  for (i = istart; i < iend; i++)
    {
      /* Get the x-coordinate */
      x = valmin[0] + i * grid_sep;

      /* Reset x-pointer for sorting vertices in x */
      sorted_vertex_ptr[0] = first_vertex_xsort_ptr[ipos];

      /* Initialise j-position */
      jpos = 0;

      /* Loop through the y points */
      for (j = jstart; j < jend; j++)
        {
          /* Get the y-coordinate */
          y = valmin[1] + j * grid_sep;

          /* Reset y-pointer for sorting vertices in y */
          sorted_vertex_ptr[1] = first_vertex_ysort_ptr[jpos];

          /* Initialise k-position */
          kpos = 0;

          /* Loop through the z points */
          for (k = kstart; k < kend; k++)
            {
              /* Get the z-coordinate */
              z = valmin[2] + k * grid_sep;

              /* Get the array location of this grid-point */
              ipoint = i * npoints[1] * npoints[2]
                + j * npoints[2] + k;

              /* Get the grid-value at the current point */
              grid_value
                = get_grid_value(grid[ipoint],&masked,contour_level);

              /* Reset z-pointer for sorting vertices in z */
              sorted_vertex_ptr[2] = first_vertex_zsort_ptr[kpos];

              /* Loop over the three possible directions */
              for (idir = 0; idir < 3; idir++)
                {
                  /* Get the array location of the adjacent point */
                  jpoint = ipoint + add_point[idir][0]
                    + add_point[idir][1] + add_point[idir][2];

                  /* Get the grid-value of the adjacent point */
                  other_value
                    = get_grid_value(grid[jpoint],&masked,contour_level);

                  /* Check for a cross-over in the contour levels */
                  if ((grid_value <= contour_level &&
                       other_value > contour_level) ||
                      (grid_value > contour_level &&
                       other_value <= contour_level))
                    {
                      /* Get the coordinates of the second point */
                      other_x = x + add_sep[idir][0];
                      other_y = y + add_sep[idir][1];
                      other_z = z + add_sep[idir][2];

                      /* Calculate fractional distance along line
                         joining the two grid-points at which the
                         cross-over occurs */
                      diff1 = abs(grid_value - contour_level);
                      diff2 = abs(other_value - contour_level);
                      fraction = (diff1 / (diff1 + diff2));
                      if (fraction > 1.0)
                        fraction = 1.0;

                      /* Calculate the coordinates of the cross-over
                         point */
                      cross[0] = x + fraction * add_sep[idir][0];
                      cross[1] = y + fraction * add_sep[idir][1];
                      cross[2] = z + fraction * add_sep[idir][2];

                      /* Initialise atom pointer associated with vertex */
                      atom_ptr = NULL;

                      /* If have atom pointers stored in the grid-atom-map
                         then retrieve closest atom */
                      if (params.store_atom_pointers == TRUE)
                        {
                          /* Get the closer of the two grid-points */
                          if (fraction < 0.5)
                            ip = ipoint;
                          else
                            ip = jpoint;

                          /* Get the atom associated with this grid-point */
                          if (grid_atom_map[ip] > -1)
                            {
                              /* Extract the atom number */
                              iatom = grid_atom_map[ip] / 1000;

                              /* If a valid atom number, then retrieve the
                                 corresponding atom pointer */
                              if (iatom > -1 && iatom < natoms)
                                atom_ptr = atom_index_ptr[iatom];
                            }
                        }

                      /* Create a vertex record for this cross-over point */
                      vertex_ptr
                        = create_vertex_record(&first_vertex_ptr,
                                               &last_vertex_ptr,cross,
                                               nvertices,edge_no,atom_ptr);

                      /* Place the vertex into its rightful position
                         in each of the three sorted lists */
                      place_vertex_in_sort_list(0,ipos,vertex_ptr,
                                                first_vertex_xsort_ptr,
                                                last_vertex_xsort_ptr,
                                                sorted_vertex_ptr);
                      place_vertex_in_sort_list(1,jpos,vertex_ptr,
                                                first_vertex_ysort_ptr,
                                                last_vertex_ysort_ptr,
                                                sorted_vertex_ptr);
                      place_vertex_in_sort_list(2,kpos,vertex_ptr,
                                                first_vertex_zsort_ptr,
                                                last_vertex_zsort_ptr,
                                                sorted_vertex_ptr);

                      /* As next z-value is bound to be the same
                         or higher than this one, make current
                         vertex the start vertex for the search */
                      sorted_vertex_ptr[2] = vertex_ptr;

                      /* Increment count of vertices */
                      nvertices++;

                      /* If have an associated atom, store cluster number */
                      if (atom_ptr != NULL)
                        atom_ptr->icluster = icluster;
                    }

                  /* Calculate this edge's reference number */
                  edge_no++;
                }

              /* Increment k-position */
              kpos++;
            }

          /* Increment j-position */
          jpos++;
        }

      /* Increment i-position */
      ipos++;
    }

  /* Loop through the x-bands to splice together the x-sorted vertex
     list */
  splice_bands_together(0,ni,first_vertex_xsort_ptr,last_vertex_xsort_ptr,
                        first_sorted_vertex_ptr);

  /* Loop through the y-bands to splice together the y-sorted vertex
     list */
  splice_bands_together(1,nj,first_vertex_ysort_ptr,last_vertex_ysort_ptr,
                        first_sorted_vertex_ptr);

  /* Loop through the z-bands to splice together the z-sorted vertex
     list */
  splice_bands_together(2,nk,first_vertex_zsort_ptr,last_vertex_zsort_ptr,
                        first_sorted_vertex_ptr);

  /* Return vertex pointer */
  *fst_vertex_ptr = first_vertex_ptr;

  /* Return number of vertices stored */
  return(nvertices);
}
/***********************************************************************

calculate_energy  -  Determine the colour based on the
                             electrostatic potential

***********************************************************************/

double calculate_energy(float *potential,int ngrid,double grid_sep,
                        double valmin[3],double coord[3])
{
  int i, icoord, ip, ipoint[3], ix, iy, iz, j, k;

  int outside;

  double calc, energy, factor[2][3], pcoord;

  /* Initialise variables */
  energy = 0.0;
  outside = FALSE;

  /* Loop over the 3 coordinates */
  for (icoord = 0; icoord < 3 && outside == FALSE; icoord++)
    {
      /* Find the reference corner grid-point defining the 3D cell
         in which the coordinates are located */
      calc = (coord[icoord] - valmin[icoord]) / grid_sep;
      if (calc < 0.0)
        ip = (int) calc - 1;
      else
        ip = (int) calc;

      /* Store the grid-point reference */
      ipoint[icoord] = ip;

      /* Check whether point lies outside the grid */
      if (ip < 0 || ip > ngrid - 2)
        outside = TRUE;

      /* Calculate the interpolation factor */
      else
        {
          /* Calculate the coordinate of this the nearest grid-point */
          pcoord = ip * grid_sep + valmin[icoord];

          /* Calculate factor depending on distance from point */
          factor[1][icoord] = (coord[icoord] - pcoord) / grid_sep;
          factor[0][icoord] = 1.0 - factor[1][icoord];
        }
    }

  /* If point is inside grid, calculate interpolated energy */
  if (outside == FALSE)
    {
      /* Loop over the 2 z-coords */
      for (k = 0; k < 2; k++)
        {
          /* Calculate array subscript */
          iz = ipoint[2] + k;

          /* Loop over the 2 y-coords */
          for (j = 0; j < 2; j++)
            {
              /* Calculate array subscript */
              iy = ipoint[1] + j;

              /* Loop over the 2 x-coords */
              for (i = 0; i < 2; i++)
                {
                  /* Calculate array subscript */
                  ix = ipoint[0] + i;

                  /* Calculate location in array */
                  ip = iz * ngrid * ngrid + iy * ngrid + ix;

                  /* Calculate interpolated value */
                  energy = energy
                    + factor[i][0] * factor[j][1] * factor[k][2]
                    * potential[ip];
                }
            }
        }
    }

  /* Return calculated energy */
  return(energy);
}
/***********************************************************************

calc_vertex_energies  -  Calculate the potential energy at each vertex
                         and convert to a colour

***********************************************************************/

void calc_vertex_energies(struct vertex *first_vertex_ptr,
                          float *potential,int ngrid,double grid_sep,
                          double *valmin)
{
  double energy, value;

  struct vertex *vertex_ptr;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to calculate energies */
  while (vertex_ptr != NULL)
    {
      /* Calculate the electrostatic energy at this point */
      energy = calculate_energy(potential,ngrid,grid_sep,valmin,
                                vertex_ptr->coords);

      /* Cut off at extreme values */
      if (energy < NEG_UPPER)
        energy = NEG_UPPER;
      else if (energy > POS_UPPER)
        energy = POS_UPPER;

      /* If energy neither very positive nor negative, set B-value
         to 40.0 */
      if (energy > NEG_LOWER && energy < POS_LOWER)
        vertex_ptr->bvalue = 40.0;

      /* Otherwise convert energy to a B-value either side of 40.0 */
      else
        {
          /* Check for positive energy */
          if (energy > 0.0)
            {
              /* Scale the energy to lie in the range 0.0 -> 1.0 */
              value = (energy - POS_LOWER) / (POS_UPPER - POS_LOWER);

              /* Bump up the lower values */
              if (value > 0.0)
                value = sqrt(value);
              else
                value = 0.0;

              /* Convert to a B-value */
              vertex_ptr->bvalue = 40.0 + 40.0 * value;
            }

          /* Check for negative energy */
          else if (energy < 0.0)
            {
              /* Scale the energy to lie in the range 0.0 -> 1.0 */
              value = (energy - NEG_LOWER) / (NEG_UPPER - NEG_LOWER);

              /* Bump up the lower values */
              if (value > 0.0)
                value = sqrt(value);
              else
                value = 0.0;

              /* Convert to a B-value */
              vertex_ptr->bvalue = 40.0 - 40.0 * value;
            }
        }

      /* Get pointer to the next vertex in the linked list */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

}
/***********************************************************************

get_ligands_in_gap  -  Get all contouring points in the current map

***********************************************************************/

int get_ligands_in_gap(int *grid,int npoints[3],double valmin[3],
                       double valmax[3],double grid_sep,
                       struct atom *first_lig_atom_ptr,
                       int contour_level,int ngrid_points,
                       char ligand_list[LISTLEN],
                       char ligand_list2[LISTLEN])
{
  char chain_str[4], add_list[20], res_num[6];

  int ipoint, istart, len, natoms, nearest_point[3];
  int grid_value, masked;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  residue_ptr = NULL;
  natoms = 0;

  /* Get pointer to the first ligand-atom */
  atom_ptr = first_lig_atom_ptr;

  /* Loop through all the ligand atoms */
  while (atom_ptr != NULL)
    {
      /* Get the nearest grid-point to this atom */
      ipoint = get_nearest_grid_point(npoints,valmin,valmax,grid_sep,
                                      atom_ptr->coords,nearest_point);

      /* If this is a valid grid-point, determine whether it lies
         within the current gap region */
      if (ipoint >= 0 && ipoint < ngrid_points)
        {
          /* Get the grid value at this point */
          grid_value = get_grid_value(grid[ipoint],&masked,contour_level);

          /* If value is above the contour level, then this ligand
             atom is inside this gap region */
          if (grid_value > contour_level)
            {
              /* If this is a new residue, write out */
              if (atom_ptr->residue_ptr != residue_ptr)
                {
                  /* Store the residue pointer */
                  residue_ptr = atom_ptr->residue_ptr;

                  /* Add residue details to ligands list */
                  sprintf(add_list,"[%s %s %c] ",residue_ptr->res_name,
                          residue_ptr->res_num,residue_ptr->chain);

                  /* Get current list length */
                  len = strlen(ligand_list);

                  /* If not already too long, add current ligand to list */
                  if (len < LISTLEN - 33)
                    {
                      /* Add on to ligand list */
                      strcat(ligand_list,add_list);
                    }

                  /* Get first non-blank char of the residue number */
                  istart = 0;
                  strcpy(res_num,residue_ptr->res_num);
                  if (res_num[0] == ' ')
                    {
                      istart = 1;
                      if (res_num[1] == ' ')
                        {
                          istart = 2;
                          if (res_num[2] == ' ')
                            istart = 3;
                        }
                    }

                  /* If no insertion code, then crop */
                  if (res_num[4] == ' ')
                    res_num[4] ='\0';

                  /* Repeat for second list */
                  sprintf(add_list,"%s&nbsp;%s",residue_ptr->res_name,
                          res_num+istart);

                  /* If have a chain, then add */
                  if (residue_ptr->chain != ' ')
                    {
                      /* Make a string of the chain identifier */
                      sprintf(chain_str,"[%c]",residue_ptr->chain);

                      /* Append to residue name and number */
                      strcat(add_list,chain_str);
                    }

                  /* Get current list length */
                  len = strlen(ligand_list2);

                  /* If not already too long, add current ligand to list */
                  if (len < LISTLEN - 33)
                    {
                      /* If not first entry, append a comma */
                      if (ligand_list2[0] != '\0')
                        strcat(ligand_list2,", ");

                      /* Add on to ligand list */
                      strcat(ligand_list2,add_list);
                    }
                }

              /* Increment count of atoms in this gap region */
              natoms++;
            }
        }

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* If any atoms found in this gap region, add the number to the ligand
     list */
  if (natoms > 0)
    {
      /* Get current list length */
      len = strlen(ligand_list);

      /* If have space, add the number */
      if (len < LISTLEN - 13)
        {
          sprintf(add_list,"%5d atoms",natoms);
          strcat(ligand_list,add_list);
        }

      /* Repeat for second list */
      len = strlen(ligand_list2);

      /* If have space, add the number */
      if (len < LISTLEN - 11)
        {
          if (natoms == 1)
            sprintf(add_list,"(1 atom)");
          else
            sprintf(add_list," (%d&nbsp;atoms)",natoms);
          strcat(ligand_list2,add_list);
        }
    }

  /* Return the number of atoms within this gap region */
  return(natoms);
}
/***********************************************************************

get_restype  -  Return the residue type for the given residue name

***********************************************************************/

int get_restype(char *res_name,char *col_desc,char *col_name)
{  
  int aa_type, iamino, itype;

  static char *aa_name[] = {
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
    "SER", "THR", "VAL", "TRP", "TYR" };

  static char *colour_name[] = { 
    "skyblue", "red", "green", "grey", "purple", "brown", "yellow" };

  static char *colour_desc[] = { 
    "positive", "negative", "neutral", "aliphatic", "aromatic",
    "Pro&Gly", "cysteine" };

  static int res_type[20]
    = { 3, 6, 1, 1, 4, 5, 0, 3, 0, 3, 3, 2, 5, 2, 0, 2, 2, 3, 4, 4 };

  /* Set defaults */
  itype = -1;
  strcpy(col_name,"black");
  strcpy(col_desc,"unknown");

  /* Initialise for search */
  aa_type = -1;

  /* Loop until residue name found */
  for (iamino = 0; iamino < NAMINO && aa_type == -1; iamino++)
    {
      if (!strcmp(res_name,aa_name[iamino]))
        {
          /* Match found, so retrieve details */
          aa_type = iamino;
          itype = res_type[aa_type];
          strcpy(col_name,colour_name[itype]);
          strcpy(col_desc,colour_desc[itype]);
        }
    }

  /* Return the residue type */
  return(itype);
}
/***********************************************************************

analyse_vertices  -  Compile residue-type and residue conservation stats
                     for the current gap region

***********************************************************************/

void analyse_vertices(struct vertex **vert_atom_index_ptr,
                      int *index,int nvertices,
                      struct atom *first_atom_ptr,
                      struct residue *first_residue_ptr,
                      int *count_type,int *count_cons)
{
  char col_desc[30], col_name[20];
    
  int iband, ipoint, irange, itype, ivertex;

  struct atom *atom_ptr;
  struct residue *residue_ptr;
  struct vertex *vertex_ptr;

  /* Initialise counts */
  for (itype = 0; itype < NRESTYPES; itype++)
    count_type[itype] = 0;
  for (iband = 0; iband < NCONSERV; iband++)
    count_cons[iband] = 0;

  /* Get pointer to the first of the stored residues */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to unset the checked flag */
  while (residue_ptr != NULL)
    {
      /* Initalise the residue flag */
      residue_ptr->checked = FALSE;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Repeat for all atoms */
  atom_ptr = first_atom_ptr;

  /* Loop over all the atoms */
  while (atom_ptr != NULL)
    {
      /* Unset flag */
      atom_ptr->checked = FALSE;

      /* Get pointer to the next atom in the linked list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Loop through all the vertices */
  for (ivertex = 0; ivertex < nvertices; ivertex++)
    {
      /* Get pointer to this vertex */
      ipoint = index[ivertex];
      vertex_ptr = vert_atom_index_ptr[ipoint];

      /* If have an associated atom, then access it residue */
      if (vertex_ptr->atom_ptr != NULL)
        {
          /* Get the atom and flag it as adjacent to the current cleft */
          atom_ptr = vertex_ptr->atom_ptr;
          atom_ptr->checked = TRUE;

          /* Get the residue */
          residue_ptr = atom_ptr->residue_ptr;

          /* If residue not already flagged as processed, then process
             now */
          if (residue_ptr->checked == FALSE)
            {
              /* Get the residue type */
              itype = get_restype(residue_ptr->res_name,col_desc,
                                  col_name);

              /* Update the residue-type stats */
              if (itype != -1)
                count_type[itype]++;

              /* Get the conservation score's band */
              if (residue_ptr->conservation_score >= 0.0)
                {
                  /* Initialise band */
                  iband = -1;

                  /* Loop until hit the correct band */
                  for (irange = 0; irange < NCONSERV && iband == -1; irange++)
                    {
                      /* If conservation value less than the upper bound
                         of this range, then have the band */
                      if (residue_ptr->conservation_score < RANGE[irange + 1])
                        {
                          /* Update the counts */
                          iband = irange;
                          count_cons[iband]++;
                        }
                    }
                }

              /* Flag as done */
              residue_ptr->checked = TRUE;
            }
        }
    }
}
/***********************************************************************

write_cleft_atoms  -  Write out atoms adjacent to the current cleft to
                      the .dbt file

***********************************************************************/

int write_cleft_atoms(int cluster,struct atom *first_atom_ptr,
                      FILE *file_ptr,int atom_number)
{
  float bvalue, occupancy;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  bvalue = (float) cluster;
  occupancy = 1.0;

  /* Write out cluster number */
  fprintf(file_ptr,"REMARK  CLEFT NUMBER %d\n",cluster);

  /* Get the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop over all the atoms */
  while (atom_ptr != NULL)
    {
      /* If atom adjacent to current cluster, then write out */
      if (atom_ptr->checked == TRUE)
        {
          /* Get the corresponding residue */
          residue_ptr = atom_ptr->residue_ptr;

          /* Increment atom count*/
          atom_number++;
          if (atom_number > 99999)
            atom_number = 99999;

          /* Write out the atom record */
          fprintf(file_ptr,"ATOM  %5d ",atom_number);
          fprintf(file_ptr,"%4s ",atom_ptr->atom_name);
          fprintf(file_ptr,"%3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  residue_ptr->res_name,residue_ptr->chain,
                  residue_ptr->res_num,atom_ptr->coords[0],
                  atom_ptr->coords[1],atom_ptr->coords[2],occupancy,
                  bvalue);
        }

      /* Get pointer to the next atom in the linked list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Return the current atom number */
  return(atom_number);
}
/***********************************************************************

relink_vertices  -  Relink vertices for later writing out to
                    speedfill.depth file

***********************************************************************/

int relink_vertices(struct vertex **fst_vertex_ptr,
                    struct vertex **lst_vertex_ptr,
                    struct vertex **vert_atom_index_ptr,int *index,
                    int nvertices,int icluster)
{
  int ipoint, ivertex, nvert;

  struct vertex *vertex_ptr;
  struct vertex *first_vertex_ptr, *last_vertex_ptr;

  /* Initialise variables */
  nvert = 0;

  /* Initialise pointers */
  first_vertex_ptr = NULL;
  last_vertex_ptr = NULL;

  /* Loop through all the vertices */
  for (ivertex = 0; ivertex < nvertices; ivertex++)
    {
      /* Get pointer to this vertex */
      ipoint = index[ivertex];
      vertex_ptr = vert_atom_index_ptr[ipoint];

      /* Skip if vertex is a deleted one */
      if (vertex_ptr->deleted == FALSE)
        {
          /* Store cluster to which vertex belongs */
          vertex_ptr->cluster = icluster;

          /* If this is the first vertex so far, then store as first */
          if (first_vertex_ptr == NULL)
            first_vertex_ptr = vertex_ptr;

          /* Add link from previous vertex to the current one */
          if (last_vertex_ptr != NULL)
            last_vertex_ptr->next_vertex_ptr = vertex_ptr;
          last_vertex_ptr = vertex_ptr;

          /* Make current vertex point to nothing */
          vertex_ptr->next_vertex_ptr = NULL;

          /* Increment count of vertices */
          nvert++;
        }
    }

  /* Return start and end vertex pointers */
  *fst_vertex_ptr = first_vertex_ptr;
  *lst_vertex_ptr = last_vertex_ptr;

  /* Return number of vertices */
  return(nvert);
}
/***********************************************************************

write_surfnet_gap  -  Write out the line to the surfnet.par line that
                      will define the current gap region

***********************************************************************/

void write_surfnet_gap(FILE *file_ptr,int icluster,
                       int first_atom_number,int last_atom_number)
{
  char gap_name[6];

  int cluster_num, colour;

  /* Initialise variables */
  colour = icluster % 13;
  cluster_num = icluster + 1;

  /* Form filename for current gap region */
  if (cluster_num < 10)
    sprintf(gap_name,"gap0%d",cluster_num);
  else
    sprintf(gap_name,"gap%d",cluster_num);

  /* Write out the surfnet.par line */
  fprintf(file_ptr,"%s        B %6d   %6d  100.0      -%d",gap_name,
          first_atom_number,last_atom_number,colour + 1);
  fprintf(file_ptr,"        10      0.01\n");
}
/***********************************************************************

calc_vertex_burial  -  Determine how many of the vertices are close to
                       the protein, and how many are on the solvent side

***********************************************************************/

double calc_vertex_burial(int *grid,
                         int npoints[3],double valmin[3],double valmax[3],
                         double grid_sep,int ngrid_points,
                         int contour_level,
                         struct vertex *first_vertex_ptr,double *percen2,
                         struct parameters params)
{
  int i, ipoint, irel, j, jpoint, k, nearest_point[3], rel_loc[27];
  int nburied, nfree, nvert;
  int grid_value, masked;
  int nmasked, nempty, nover, nunder;

  double percen;

  struct vertex *vertex_ptr;

  /* Initialise variables */
  nburied = 0;
  nfree = 0;
  nvert = 0;
  nover = 0;
  nunder = 0;

  /* Initialise the relative locations of the 26 surrounding grid-points */
  irel = 0;
  for (i = -1; i < 2; i++)
    for (j = -1; j < 2; j++)
      for (k = -1; k < 2; k++)
        {
          /* Calculate relative location */
          rel_loc[irel] = i * npoints[1] * npoints[2] + j * npoints[2] + k;

          /* Increment subscript */
          irel++;
        }

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices */
  while (vertex_ptr != NULL)
    {
      /* Process only if vertex hasn't been deleted */
      if (vertex_ptr->deleted == FALSE)
        {
          /* Find the closest grid-point to this vertex */
          ipoint = get_nearest_grid_point(npoints,valmin,valmax,grid_sep,
                                          vertex_ptr->coords,nearest_point);

          /* Initialise counts of masked and unmasked neighbouring vertices */
          nmasked = 0;
          nempty = 0;

          /* Loop over the surrounding grid-points */
          for (irel = 0; irel < 27; irel++)
            {
              /* Get relative location of this grid-point */
              jpoint = ipoint + rel_loc[irel];

              /* If this is a valid location, then extract its grid-value */
              if (jpoint >= 0 && jpoint < ngrid_points)
                {
                  /* Get the grid value */
                  grid_value
                    = get_grid_value(grid[jpoint],&masked,contour_level);

                  /* If grid-point is outside the gap region, then consider */
                  if (grid_value < contour_level)
                    {
                      if (masked == MASKED)
                        {
                          /* Increment count of vertices masked by protein
                             atoms */
                          nmasked++;
                        }
                      else
                        nempty++;
                    }
                }
            }

          /* Store ratio */
          if (nmasked + nempty > 0)
            vertex_ptr->ratio = 100.0 * nmasked
              / (nmasked + nempty);
          else
            vertex_ptr->ratio = 0.0;

          /* If ratio is zero, increment free count */
          if (nmasked == 0)
            nfree++;
          else
            nburied++;

          /* Increment counts for second percentage */
          if (vertex_ptr->ratio > 90.0)
            nover++;
          else
            nunder++;
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Calculate the percentage of buried to total vertex points */
  percen = 100.0 * nburied / (nburied + nfree);

  /* Calculate second percentage of numbers of high-ratio points */
  *percen2 = 100.0 * nover / (nover + nunder);

  /* Return the percentage */
  return(percen);
}
/***********************************************************************

save_cluster_points  -  Save the current cluster's grid points

***********************************************************************/

void save_cluster_points(int *grid,int *save_grid,int ngrid_points,
                         int icluster,int contour_level)
{
  int ipoint;
  int grid_value, masked;

  /* Loop through the grid array */
  for (ipoint = 0; ipoint < ngrid_points; ipoint++)
    {
      /* If grid point is empty, check whether it is occupied by the 
         current cluster */
      if (save_grid[ipoint] == EMPTY - 1)
        {
          /* Extract the current grid-value */
          grid_value
            = get_grid_value(grid[ipoint],&masked,contour_level);

          /* If grid-point is inside the gap region, save the cluster
             number in the save array*/
          if (grid_value > contour_level)
            save_grid[ipoint] = icluster + 1;
        }
    }
}
/***********************************************************************

markx_exposed_vertices  -  Determine which vertices are "exposed"
                           (Old method)

***********************************************************************/

 int markx_exposed_vertices(struct vertex *first_vertex_ptr,int nvertices,
                           struct vertex **vertex_index_ptr)
{
  int nexposed;
  int nlist;
  int filter;

  double dist_sq, other_x, other_y, other_z, reach_dist, x, y, z;

  struct vertex *vertex_ptr, *other_vertex_ptr;
  struct vertex *first_vertex_list_ptr;

  /* Initialise variables */
  filter = FALSE;
  nexposed = 0;

  /* Calculate distance away to search when looking for vertex neighbours */
  reach_dist = SQRT_THREE * VERTEX_DIST;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to determine which are "exposed" */
  while (vertex_ptr != NULL)
    {
      /* If this is an "exposed" vertex, then consider its neighbours */
      if (vertex_ptr->deleted == FALSE && vertex_ptr->ratio == 0.0)
        {
          /* Get this vertex's coords */
          x = vertex_ptr->coords[0];
          y = vertex_ptr->coords[1];
          z = vertex_ptr->coords[2];

          /* Mark this vertex as "exposed" */
          vertex_ptr->exposed = TRUE;

          /* Get list of neighbouring vertices that might be within
             joining distance */
          nlist
            = get_neighbour_list(vertex_ptr,vertex_index_ptr,
                                 nvertices,reach_dist,
                                 &first_vertex_list_ptr,filter);

          /* Get pointer to the the first of the neighbouring vertices */
          other_vertex_ptr = first_vertex_list_ptr;

          /* Loop through all other vertices to find those that are
             close enough to the current one */
          while (other_vertex_ptr != NULL && vertex_ptr->exposed == TRUE)
            {
              /* Process only if vertex hasn't been deleted */
              if (other_vertex_ptr->deleted == FALSE)
                {
                  /* Get this vertex's coords */
                  other_x = other_vertex_ptr->coords[0];
                  other_y = other_vertex_ptr->coords[1];
                  other_z = other_vertex_ptr->coords[2];

                  /* Calculate the distance between the two vertices */
                  dist_sq = (x - other_x) * (x - other_x)
                    + (y - other_y) * (y - other_y)
                    + (z - other_z) * (z - other_z);

                  /* If vertices close enough to be neighbours, then check
                     whether second vertex is also exposed */
                  if (dist_sq < VERTEX_DIST_SQRD)
                    {
                      /* If second vertex does not have a zero ratio, then
                         can't consider current vertex as exposed */
                      if (other_vertex_ptr->ratio != 0.0)
                        {
                          /* Flag current vertex as not exposed */
                          vertex_ptr->exposed = FALSE;
                        }
                    }
                }

              /* Get pointer to the next other vertex */
              other_vertex_ptr = other_vertex_ptr->next_vertex_list_ptr;
            }

          /* If vertex has got through and is still marked as exposed,
             increment count */
          if (vertex_ptr->exposed == TRUE)
            nexposed++;
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Return number of exposed vertices */
  return(nexposed);
}
/***********************************************************************

mark_exposed_vertices  -  Determine which vertices are "exposed"

***********************************************************************/

int mark_exposed_vertices(struct vertex *first_vertex_ptr,int nvertices,
                          struct vertex **vertex_index_ptr,
                          KdTree *kdtree_ptr)
{
  int jvert, jvertex;
  int nexposed;

  double inner_radius, outer_radius;

  struct vertex *vertex_ptr, *other_vertex_ptr;

  KdTreeQuery *kdquery_ptr;
  Region *region_ptr;

  /* Initialise variables */
  inner_radius = 0.0;
  nexposed = 0;

  /* Calculate distance away to search when looking for vertex neighbours */
  outer_radius = VERTEX_DIST;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to determine which are "exposed" */
  while (vertex_ptr != NULL)
    {
      /* If this is an "exposed" vertex, then consider its neighbours */
      if (vertex_ptr->deleted == FALSE && vertex_ptr->ratio == 0.0)
        {
          /* Mark this vertex as "exposed" */
          vertex_ptr->exposed = TRUE;

          /* Create the region of space to be searched for nearby
             vertices */
          region_ptr
            = Annulus_create(vertex_ptr->coords,inner_radius,
                             outer_radius,3);

          /* Perform the search */
          kdquery_ptr = KdTree_query(kdtree_ptr,region_ptr);

          /* Loop through all the vertices matched by the search */
          jvert = 0;
          while ((jvertex = KdTreeQuery_next(kdquery_ptr)) >= 0 &&
                 vertex_ptr->exposed == TRUE)
            {
              /* Get the vertex pointer returned by the query */
              other_vertex_ptr = vertex_index_ptr[jvertex];

              /* If second vertex does not have a zero ratio, then
                 can't consider current vertex as exposed */
              if (other_vertex_ptr->ratio != 0.0)
                {
                  /* Flag current vertex as not exposed */
                  vertex_ptr->exposed = FALSE;
                }
            }

          /* If vertex has got through and is still marked as exposed,
             increment count */
          if (vertex_ptr->exposed == TRUE)
            nexposed++;
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Return number of exposed vertices */
  return(nexposed);
}
/***********************************************************************

create_vjoin_record  -  Create and initialise a new vjoin record

***********************************************************************/

struct vjoin *create_vjoin_record(struct vjoin **fst_vjoin_ptr,
                                  struct vjoin **lst_vjoin_ptr,
                                  struct vertex *vertex1_ptr,
                                  struct vertex *vertex2_ptr)
{
  struct vjoin *vjoin_ptr, *first_vjoin_ptr, *last_vjoin_ptr;

  /* Initialise pointers */
  vjoin_ptr = NULL;
  first_vjoin_ptr = *fst_vjoin_ptr;
  last_vjoin_ptr = *lst_vjoin_ptr;

  /* Allocate memory for structure to hold vjoin info */
  vjoin_ptr = (struct vjoin *) malloc(sizeof(struct vjoin));
  if (vjoin_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct vjoin\n");
      exit (1);
    }

  /* If this is the very first vjoin, then store pointer */
  if (first_vjoin_ptr == NULL)
    first_vjoin_ptr = vjoin_ptr;

  /* Add link from previous vjoin to the current one */
  if (last_vjoin_ptr != NULL)
    last_vjoin_ptr->next_vjoin_ptr = vjoin_ptr;
  last_vjoin_ptr = vjoin_ptr;

  /* Store pointers to the paired vertices */
  vjoin_ptr->vertex1_ptr = vertex1_ptr;
  vjoin_ptr->vertex2_ptr = vertex2_ptr;

  /* Initialise the other elements of the structure */
  vjoin_ptr->next_vjoin_ptr = NULL;

  /* Return current pointers */
  *fst_vjoin_ptr = first_vjoin_ptr;
  *lst_vjoin_ptr = last_vjoin_ptr;

  /* Return new vjoin pointer */
  return(vjoin_ptr);
}
/***********************************************************************

getx_vertex_joins  -  Determine which pairs of "exposed" vertices are
                      neighbours (Old method)

***********************************************************************/

int getx_vertex_joins(struct vertex *first_vertex_ptr,
                      struct vjoin **fst_vjoin_ptr,int nvertices,
                      struct vertex **vertex_index_ptr)
{
  int njoins, nlist;
  int filter;

  double dist_sq, other_x, other_y, other_z, reach_dist, x, y, z;

  struct vertex *vertex_ptr, *other_vertex_ptr;
  struct vjoin *first_vjoin_ptr, *last_vjoin_ptr;
  struct vertex *first_vertex_list_ptr;

  /* Initialise variables */
  filter = TRUE;
  first_vjoin_ptr = NULL;
  last_vjoin_ptr = NULL;
  njoins = 0;

  /* Calculate distance away to search when looking for vertex neighbours */
  reach_dist = SQRT_THREE * VERTEX_DIST;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to find pairs of "exposed" ones that are
     within the cut-off distance */
  while (vertex_ptr != NULL)
    {
      /* If this is an "exposed" vertex, then consider its neighbours */
      if (vertex_ptr->deleted == FALSE && vertex_ptr->exposed == TRUE)
        {
          /* Get this vertex's coords */
          x = vertex_ptr->coords[0];
          y = vertex_ptr->coords[1];
          z = vertex_ptr->coords[2];

          /* Get list of neighbouring vertices that might be within
             joining distance */
          nlist
            = get_neighbour_list(vertex_ptr,vertex_index_ptr,
                                 nvertices,reach_dist,
                                 &first_vertex_list_ptr,filter);

          /* Get pointer to the next vertex */
          other_vertex_ptr = first_vertex_list_ptr;

          /* Loop through all other vertices to find those that are
             close enough to the current one */
          while (other_vertex_ptr != NULL)
            {
              /* If this vertex is also exposed, calculate the distance
                 between them */
              if (other_vertex_ptr->deleted == FALSE &&
                  other_vertex_ptr->exposed == TRUE)
                {
                  /* Get this vertex's coords */
                  other_x = other_vertex_ptr->coords[0];
                  other_y = other_vertex_ptr->coords[1];
                  other_z = other_vertex_ptr->coords[2];

                  /* Calculate the squared distance */
                  dist_sq = (x - other_x) * (x - other_x)
                    + (y - other_y) * (y - other_y)
                    + (z - other_z) * (z - other_z);

                  /* If vertices close enough, create a virtual bond
                     between them */
                  if (dist_sq < VERTEX_DIST_SQRD)
                    {
                      /* Create a vjoin record */
                      create_vjoin_record(&first_vjoin_ptr,&last_vjoin_ptr,
                                          vertex_ptr,other_vertex_ptr);

                      /* Increment count of vertex pairs joined */
                      njoins++;
                    }
                }

              /* Get pointer to the next other vertex */
              other_vertex_ptr = other_vertex_ptr->next_vertex_list_ptr;
            }
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Return pointer to the first of the vertex pairs */
  *fst_vjoin_ptr = first_vjoin_ptr;

  /* Return number of vertex joins found */
  return(njoins);
}
/***********************************************************************

get_vertex_joins  -  Determine which pairs of "exposed" vertices are
                     neighbours

***********************************************************************/

int get_vertex_joins(struct vertex *first_vertex_ptr,
                     struct vjoin **fst_vjoin_ptr,int nvertices,
                     struct vertex **vertex_index_ptr,
                     KdTree *kdtree_ptr)
{
  int jvert, jvertex;
  int njoins;

  double inner_radius, outer_radius;

  struct vertex *vertex_ptr, *other_vertex_ptr;
  struct vjoin *first_vjoin_ptr, *last_vjoin_ptr;

  KdTreeQuery *kdquery_ptr;
  Region *region_ptr;

  /* Initialise variables */
  inner_radius = 0.0;
  first_vjoin_ptr = NULL;
  last_vjoin_ptr = NULL;
  njoins = 0;

  /* Calculate distance away to search when looking for vertex neighbours */
  outer_radius = VERTEX_DIST;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to find pairs of "exposed" ones that are
     within the cut-off distance */
  while (vertex_ptr != NULL)
    {
      /* If this is an "exposed" vertex, then consider its neighbours */
      if (vertex_ptr->deleted == FALSE && vertex_ptr->exposed == TRUE)
        {
          /* Create the region of space to be searched for nearby
             vertices */
          region_ptr
            = Annulus_create(vertex_ptr->coords,inner_radius,
                             outer_radius,3);

          /* Perform the search */
          kdquery_ptr = KdTree_query(kdtree_ptr,region_ptr);

          /* Loop through all the vertices matched by the search */
          jvert = 0;
          while ((jvertex = KdTreeQuery_next(kdquery_ptr)) >= 0)
            {
              /* Get the vertex pointer returned by the query */
              other_vertex_ptr = vertex_index_ptr[jvertex];

              /* If this vertex is also exposed, calculate the distance
                 between them */
              if (other_vertex_ptr->exposed == TRUE &&
                  other_vertex_ptr->vertex_no > vertex_ptr->vertex_no)
                {
                  /* Create a vjoin record */
                  create_vjoin_record(&first_vjoin_ptr,&last_vjoin_ptr,
                                      vertex_ptr,other_vertex_ptr);

                  /* Increment count of vertex pairs joined */
                  njoins++;
                }
            }
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Return pointer to the first of the vertex pairs */
  *fst_vjoin_ptr = first_vjoin_ptr;

  /* Return number of vertex joins found */
  return(njoins);
}
/***********************************************************************

get_depth  -  Calculate the maximum depth from the given triangle

***********************************************************************/

double get_depth(double x1,double y1,double z1,double x2,double y2,
                 double z2,double x3,double y3,double z3,
                 struct vertex *first_vertex_ptr)
{
  double ax, ay, az, bx, by, bz, cx, cy, cz, nx, ny, nz;
  double midx, midy, midz, x, y, z, vx, vy, vz;
  double depth, dist_sqrd, max_depth;
  double magnitude;

  struct atom *atom_ptr;
  struct vertex *vertex_ptr;

  /* Initialise variables */
  max_depth = 0.0;

  /* Calculate the mid-point of this triangle */
  midx = (x1 + x2 + x3) / 3.0;
  midy = (y1 + y2 + y3) / 3.0;
  midz = (z1 + z2 + z3) / 3.0;

  /* Define two vectors in the plane of the triangle */
  ax = x2 - x1;
  ay = y2 - y1;
  az = z2 - z1;
  bx = x3 - x1;
  by = y3 - y1;
  bz = z3 - z1;

  /* Calculate the normal to this triangle */
  nx = ay * bz - az * by;
  ny = az * bx - ax * bz;
  nz = ax * by - ay * bx;

  /* Normalise it */
  magnitude = sqrt(nx * nx + ny * ny + nz * nz);
  nx = nx / magnitude;
  ny = ny / magnitude;
  nz = nz / magnitude;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices to find pairs of "exposed" ones that are
     within the cut-off distance */
  while (vertex_ptr != NULL)
    {
      /* If this is a buried vertex, then see how far it lies from the
         normal to the triangle */
      if (vertex_ptr->deleted == FALSE && vertex_ptr->exposed == FALSE)
        {
          /* Get this vertex's coords */
          x = vertex_ptr->coords[0];
          y = vertex_ptr->coords[1];
          z = vertex_ptr->coords[2];

          /* Form the vector from vertex to the midpoint of the triangle */
          vx = x - midx;
          vy = y - midy;
          vz = z - midz;

          /* Use the cross product with the normal to calculate its distance
             from the normal */
          cx = vy * nz - vz * ny;
          cy = vz * nx - vx * nz;
          cz = vx * ny - vy * nx;
          dist_sqrd = cx * cx + cy * cy + cz * cz;

          /* If this distance is within the cut-off, then calculate the
             depth */
          if (dist_sqrd < VERTEX_DIST_SQRD)
            {
              /* Calculate the distance between the vertex and the mid-point
                 of the triangle */
              depth = sqrt((x - midx) * (x - midx)
                                   + (y - midy) * (y - midy)
                                   + (z - midz) * (z - midz));

              /* If this is the largest depth encountered so far, then
                 store */
              if (depth > max_depth)
                max_depth = depth;

              /* If have an atom associated with this vertex, store
                 depth for atom */
              atom_ptr = vertex_ptr->atom_ptr;
              if (atom_ptr != NULL)
                {
                  /* Add to depth accumulator and increment count of
                     depths for averaging */
                  atom_ptr->depth = atom_ptr->depth + depth;
                  atom_ptr->ndepth++;
                }
            }
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Return maximum depth here */
  return(max_depth);
}
/***********************************************************************

calculate_gap_depth  -  Calculate the "depth" of the current gap region
                        by considering triangles of fully "exposed"
                        vertices and calculating the distance along the
                        normal to the nearest "buried" vertex

***********************************************************************/

double calculate_gap_depth(struct vertex *first_vertex_ptr,int nvertices,
                           struct vertex **vertex3_index_ptr,
                           struct vertex **vertex_index_ptr,
                           int filter_method)
{
  int nexposed, njoins, ntriang;
  int done;

  double average_depth, depth, max_depth;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;

  struct vertex *vertex1_ptr, *vertex2_ptr, *vertex3_ptr;
  struct vjoin *vjoin_ptr, *first_vjoin_ptr;
  struct vjoin *first_vjoin2_ptr, *vertex1_vjoin_ptr, *vertex2_vjoin_ptr;
  struct vjoin *next_vjoin_ptr;

  KdTree *kdtree_ptr;

  /* Initialise variables */
  average_depth = 0.0;
  depth = 0.0;
  max_depth = 0.0;
  ntriang = 0;

  /* Use old or new method, the latter employing kd-trees for quickly
     finding vertex neighbour lists */
  if (filter_method == OLD)
    {
      /* Determine which vertices are "exposed" - ie totally surrounded by
         other vertices that have a ratio of 0.0 */
      nexposed
        = markx_exposed_vertices(first_vertex_ptr,nvertices,
                                 vertex3_index_ptr);
      /* Get all pairs of "exposed" vertices that are neighbours */
      njoins
        = getx_vertex_joins(first_vertex_ptr,&first_vjoin_ptr,nvertices,
                            vertex3_index_ptr);
    }

  /* New method uses kd-trees for the vertices */
  else
    {
      /* Generate the kd-tree for the vertices */
      nvertices
        = create_vertex_kdtree(first_vertex_ptr,nvertices,&kdtree_ptr,
                               vertex_index_ptr);

      /* Determine which vertices are "exposed" - ie totally surrounded by
         other vertices that have a ratio of 0.0 */
      nexposed
        = mark_exposed_vertices(first_vertex_ptr,nvertices,
                                vertex_index_ptr,kdtree_ptr);

      /* Get all pairs of "exposed" vertices that are neighbours */
      njoins
        = get_vertex_joins(first_vertex_ptr,&first_vjoin_ptr,nvertices,
                           vertex_index_ptr,kdtree_ptr);
    }

  /* Get pointer to the first vertex pair */
  vjoin_ptr = first_vjoin_ptr;

  /* Loop over all the stored vertex pairs */
  while (vjoin_ptr != NULL)
    {
      /* Get the two vertices making up this vertex pair */
      vertex1_ptr = vjoin_ptr->vertex1_ptr;
      vertex2_ptr = vjoin_ptr->vertex2_ptr;

      /* Get their coords */
      x1 = vertex1_ptr->coords[0];
      y1 = vertex1_ptr->coords[1];
      z1 = vertex1_ptr->coords[2];
      x2 = vertex2_ptr->coords[0];
      y2 = vertex2_ptr->coords[1];
      z2 = vertex2_ptr->coords[2];

      /* Get the next join record */
      next_vjoin_ptr = vjoin_ptr->next_vjoin_ptr;
      first_vjoin2_ptr = NULL;
      done = FALSE;

      /* Loop through all subsequent vertex-join records to find the
         first of the joins involving our current vertex 2 */
      while (next_vjoin_ptr != NULL && done == FALSE)
        {
          /* If it is our vertex 2, then have the start-point we are
             looking for */
          if (next_vjoin_ptr->vertex1_ptr == vertex2_ptr)
            {
              /* Save the starting record */
              first_vjoin2_ptr = next_vjoin_ptr;
              done = TRUE;
            }

          /* Otherwise, check whether the vertex number is higher and
             we have overshot (ie vertex 2 has no start-point of its own) */
          else if (next_vjoin_ptr->vertex1_ptr->vertex_no
                   > vertex2_ptr->vertex_no)
            done = TRUE;

          /* Get pointer to the next vertex pair */
          next_vjoin_ptr = next_vjoin_ptr->next_vjoin_ptr;
        }

      /* If have a valid start-point, start looking for triplets involving
         our current vertices 1 and 2 */
      if (first_vjoin2_ptr != NULL)
        {
          /* Get the next join record */
          vertex1_vjoin_ptr = vjoin_ptr->next_vjoin_ptr;

          /* Loop until all neighbours of current vertex 1 have been
             looked at */
          while (vertex1_vjoin_ptr != NULL &&
                 vertex1_vjoin_ptr->vertex1_ptr == vertex1_ptr)
            {
              /* Get the other vertex from this join */
              vertex3_ptr = vertex1_vjoin_ptr->vertex2_ptr;

              /* Only consider vertex 3 if its number is higher than that
                 of vertex 2 */
              if (vertex3_ptr->vertex_no > vertex2_ptr->vertex_no)
                {
                  /* Get this vertex's coords */
                  x3 = vertex3_ptr->coords[0];
                  y3 = vertex3_ptr->coords[1];
                  z3 = vertex3_ptr->coords[2];

                  /* Get first of vertex 2's joins */
                  vertex2_vjoin_ptr = first_vjoin2_ptr;

                  /* Loop until all neighbours of current vertex 2 have been
                     looked at */
                  while (vertex2_vjoin_ptr != NULL &&
                         vertex2_vjoin_ptr->vertex1_ptr == vertex2_ptr)
                    {
                      /* If the second vertex in this join matches our
                         current vertex 3, then have a triangle */
                      if (vertex2_vjoin_ptr->vertex2_ptr == vertex3_ptr)
                        {
                          /* Calculate the longest depth from this triangle
                             to a "buried" vertex */
                          depth = get_depth(x1,y1,z1,x2,y2,z2,x3,y3,z3,
                                            first_vertex_ptr);

                          /* If this is the deepest encountered so far, then
                             store */
                          if (depth > max_depth)
                            max_depth = depth;
                          average_depth = average_depth + depth;

                          /* Increment count of triangles considered */
                          ntriang++;
                        }

                      /* Get pointer to the next vertex pair */
                      vertex2_vjoin_ptr = vertex2_vjoin_ptr->next_vjoin_ptr;
                    }
                }

              /* Get pointer to the next vertex pair */
              vertex1_vjoin_ptr = vertex1_vjoin_ptr->next_vjoin_ptr;
            }
        }

      /* Get pointer to the next vertex pair */
      vjoin_ptr = vjoin_ptr->next_vjoin_ptr;
    }

  /* Calculate the average depth */
  if (ntriang > 0)
    average_depth = average_depth / ntriang;
  printf("Number of triangles %d, average gap depth = %8.3f\n",
         ntriang,average_depth);

  /* Return average gap depth */
  return(average_depth);
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

principal_components  -  Calculate the transformation to align the
                          supplied atom coords with their principal axes
                          along the x-, y-, and z-directions

***********************************************************************/

void principal_components(struct vertex *first_vertex_ptr,
                          double eigen_value[3])
{
  int nzeroy, nzeroz;
  int array_size, i, ipos, j, npoints, nsize;
  double amatrix[6], e1, matrix[3][3], rmatrix[9], x, y, z;
  double xxsum, xysum, xzsum, yysum, yzsum, zzsum;
  struct vertex *vertex_ptr;

  /* Initialise counts */
  npoints = 0;
  nzeroy = 0;
  nzeroz = 0;
  xxsum = 0.0;
  xysum = 0.0;
  xzsum = 0.0;
  yysum = 0.0;
  yzsum = 0.0;
  zzsum = 0.0;

  /* Initialise the eigen values */
  for (i = 0; i < 3; i++)
    eigen_value[i] = 1.0;

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices */
  while (vertex_ptr != NULL)
    {
      /* Process only if vertex hasn't been deleted */
      if (vertex_ptr->deleted == FALSE)
        {
          /* Retrieve the coordinates for this vertex */
          x = vertex_ptr->coords[0];
          y = vertex_ptr->coords[1];
          z = vertex_ptr->coords[2];

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
        }

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
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

      /* Normalise the eigen values to make largest one 1.0 */
      e1 = eigen_value[0];
      for (i = 0; i < 3; i++)
        eigen_value[i] = eigen_value[i] / e1;
    }
}
/***********************************************************************

perform_vertex_sort  -  Shell-sort which sorts the vertices into
                        ascending order of associated atom name

***********************************************************************/

void perform_vertex_sort(int *index,struct atom **atom_ind_ptr,int nindex)
{
  char key[50], other_key[50];

  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  int ktest;
  double aln2i = 1.4426950, tiny = 1.e-5;
  double calc;
  double dble, dble1;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /*--Initialise values */
  n = nindex;
  dble = n;
  dble1 = log(dble);
  calc = dble1;
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

              /* Get the atom pointer for the first vertex */
              atom_ptr = atom_ind_ptr[index[indi]];

              /* If valid, then form key */
              if (atom_ptr != NULL)
                {
                  /* Get the corresponding residue pointer */
                  residue_ptr = atom_ptr->residue_ptr;

                  /* String the atom and residue details together */
                  sprintf(key,"%c%s%s%s",residue_ptr->chain,
                          residue_ptr->res_num,residue_ptr->res_name,
                          atom_ptr->atom_name);
                }
              else
                key[0] = '\0';

              /* Repeat for second vertex */
              atom_ptr = atom_ind_ptr[index[indl]];

              /* If valid, then form key */
              if (atom_ptr != NULL)
                {
                  /* Get the corresponding residue pointer */
                  residue_ptr = atom_ptr->residue_ptr;

                  /* String the atom and residue details together */
                  sprintf(other_key,"%c%s%s%s",residue_ptr->chain,
                          residue_ptr->res_num,residue_ptr->res_name,
                          atom_ptr->atom_name);
                }
              else
                key[0] = '\0';

              /* If in wrong order, the swap sort-pointers */
              ktest = strcmp(key,other_key);
              if (ktest < 0)
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

sort_vertices  -  Sort the vertices by their associated atoms

***********************************************************************/

void sort_vertices(struct vertex *first_vertex_ptr,
                   struct vertex ***vert_atom_index_ptr,
                   int **ind,int nvertices,struct parameters params)
{
  int ivertex;
  int *index;

  struct atom *atom_ptr;
  struct atom **atom_ind_ptr;
  struct vertex *vertex_ptr;
  struct vertex **vert_ind_ptr;

  /* Initialise variables */
  ivertex = 0;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(nvertices * sizeof(int));
  atom_ind_ptr = (struct atom **) malloc(nvertices * sizeof(struct atom *));
  vert_ind_ptr
    = (struct vertex **) malloc(nvertices * sizeof(struct vertex *));

  /* Get pointer to the first of the vertex records */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all the vertices to store atom pointers */
  while (vertex_ptr != NULL)
    {
      /* Process only if vertex hasn't been deleted */
      if (vertex_ptr->deleted == FALSE)
        {
          /* Get the atom associated with this vertex */
          atom_ptr = vertex_ptr->atom_ptr;

          /* If total number of vertices not exceeded, then store
             atom pointer */
          if (ivertex < nvertices)
            {
              /* Store all the details in the arrays */
              index[ivertex] = ivertex;
              atom_ind_ptr[ivertex] = atom_ptr;
              vert_ind_ptr[ivertex] = vertex_ptr;
            }

          /* Increment vertex count */
          ivertex++;
        }

      /* Get the next vertex record in the sort-list */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Sort the vertices */
  if (params.vertex_colouring != STANDARD)
    perform_vertex_sort(index,atom_ind_ptr,nvertices);

  /* Return the array pointers */
  *vert_atom_index_ptr = vert_ind_ptr;
  *ind = index;
}
/***********************************************************************

order_results  -  Order the given set of values

***********************************************************************/

void order_results(double value[MAX_TOP_GAPS],int nvalues,
                   int order[MAX_TOP_GAPS])
{
  int i, index[MAX_TOP_GAPS], iorder, j;
  int swapped;

  /* Initialise flag */
  swapped = TRUE;

  /* Initialise the order indices */
  for (iorder = 0; iorder < nvalues; iorder++)
    index[iorder] = iorder;

  /* Keep looping while any elements are still being ordered */
  while (swapped == TRUE)
    {
      /* Unset the swap flag */
      swapped = FALSE;

      /* Loop through all the values */
      for (iorder = 0; iorder < nvalues - 1; iorder++)
        {
          /* Get the values pointed to by this index element and the
             next one */
          i = index[iorder];
          j = index[iorder + 1];

          /* If their sizes are the wrong way round, then need to
             swap them */
          if (value[i] < value[j])
            {
              /* Swap the index elements */
              index[iorder] = j;
              index[iorder + 1] = i;

              /* Set flag to indicate that a swap took place in this loop */
              swapped = TRUE;
            }
        }
    }

  /* Use the index to determine the ranking of each value */
  for (iorder = 0; iorder < nvalues; iorder++)
    {
      /* Get which value is at this rank */
      i = index[iorder];

      /* Store that value's ranking position */
      order[i] = iorder;
    }
}
/***********************************************************************

write_html_gap  -  Write out the gap statistics for current gap to the
                   HTML file

***********************************************************************/

void write_html_gap(FILE *file_ptr,int icluster,double volume,double ratio1,
                    double percen,int order1,double percen2,int order2,
                    double depth,int order3,char *ligand_list,
                    int *count_type,int *count_cons)
{
  char checked[8], fcourier[30], fs[30], fe[8];
  char tmp_string[50];
  char res_cons[2][MAX_PATTN_LEN], res_types[2][MAX_PATTN_LEN]; 

  int colour, icons, itype, len_string;
  int tens, units;
  int have_tens;

  static char count_ch[]={".123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
  
  static char *col_name[] = {
    "#FF0000", "#9900FF", "#FFD700", "#0000FF", "#00FF00",
    "#8B4513", "#FF69B4", "#339966", "#FF00FF", "#00FFFF",
    "#FF6600", "#EE82EE", "#FF9900", "#FFFFFF"
  };

  static char *res_type_colour[] = { 
    "skyblue", "red", "green", "gray", "purple", "brown", "yellow" };

  static char *res_cons_colour[] = { 
    "blue", "purple", "skyblue", "cyan", "gray", "green", "yellow",
    "orange", "pink", "red" };

  /* Initialise variables */
  len_string = strlen(count_ch);

  /* Get this cluster's colour */
  colour = (icluster - 1) % 13;

  /* Create the FONT start- and end-commands */
  sprintf(fcourier,"<FONT FACE=\"COURIER\">");
  sprintf(fs,"<FONT COLOR=\"%s\">",col_name[colour]);
  strcpy(fe,"</FONT>");

  /* If one of the top 4 clusters, then show checked */
  if (icluster < 5)
    strcpy(checked,"CHECKED");
  else
    checked[0] = '\0';

  /* Write out the details */
  fprintf(file_ptr,"<TR>\n");

  /* First column contains the check-box and the gap number */
  fprintf(file_ptr,"<TD ALIGN=CENTER>");
  fprintf(file_ptr," <INPUT TYPE=\"checkbox\" NAME=\"gap\" ");
  fprintf(file_ptr,"VALUE=\"-show %d\" %s> ",icluster,checked);
  fprintf(file_ptr,"%s%d%s</TD>\n",fs,icluster,fe);

  /* Write out remaining columns */
  fprintf(file_ptr,"<TD ALIGN=CENTER>%s%12.2f%s</TD>\n",fs,volume,fe);
  if (icluster == 1)
    fprintf(file_ptr,"<TD ALIGN=CENTER>%s%8.2f%s</TD>\n",fs,ratio1,fe);
  else
    fprintf(file_ptr,"<TD ALIGN=CENTER>%s-%s</TD>\n",fs,fe);
  fprintf(file_ptr,"<TD BGCOLOR=\"#CCCCCC\" ALIGN=CENTER>%s%8.2f%%%s</TD>\n",
          fs,percen,fe);
  fprintf(file_ptr,"<TD BGCOLOR=\"#CCCCCC\" ALIGN=CENTER>%s%2d%s</TD>\n",
          fs,order1,fe);
  fprintf(file_ptr,"<TD BGCOLOR=\"#CCCCCC\" ALIGN=CENTER>%s%8.2f%%%s</TD>\n",
          fs,percen2,fe);
  fprintf(file_ptr,"<TD BGCOLOR=\"#CCCCCC\" ALIGN=CENTER>%s%2d%s</TD>\n",
          fs,order2,fe);
  fprintf(file_ptr,"<TD BGCOLOR=\"#CCCCCC\" ALIGN=CENTER>%s%8.2f%s</TD>\n",
          fs,depth,fe);
  fprintf(file_ptr,"<TD BGCOLOR=\"#CCCCCC\" ALIGN=CENTER>%s%2d%s</TD>\n",
          fs,order3,fe);

  /* Initialise patterns */
  res_types[0][0] = '\0';
  res_cons[0][0] = '\0';
  res_types[1][0] = '\0';
  res_cons[1][0] = '\0';

  /* Form the pattern of residue-type counts for this gap */
  have_tens = FALSE;
  for (itype = 0; itype < NRESTYPES; itype++)
    {
      /* Check whether count is a 2-digit number */
      if (count_type[itype] > 9)
        {
          /* Set flag */
          have_tens = TRUE;

          /* Get the tens number */
          tens = count_type[itype] / 10;  

          /* Get the colour for this residue type */
          sprintf(tmp_string,"<FONT COLOR=%s>%d</FONT>",
                  res_type_colour[itype],tens);
        }

      /* Otherwise, take tens position to be blank */
      else
        strcpy(tmp_string,"&nbsp;");

      /* Add to tens line */
      strcat(res_types[0],tmp_string);

      /* Get the units digit */
      units = count_type[itype] % 10;

      /* Get the colour for this residue type */
      if (count_type[itype] == 0)
        sprintf(tmp_string,"<FONT COLOR=%s>.</FONT>",
                res_type_colour[itype]);
      else
        sprintf(tmp_string,"<FONT COLOR=%s>%d</FONT>",
                res_type_colour[itype],units);

      /* Add onto end of current string */
      strcat(res_types[1],tmp_string);
    }

  /* Write out the residue type stats */
  fprintf(file_ptr,"<TH>");
  if (have_tens == TRUE)
    fprintf(file_ptr,"%s%s%s<BR>\n",fcourier,res_types[0],fe);
  fprintf(file_ptr,"%s%s%s</TH>\n",fcourier,res_types[1],fe);

  /* Repeat for the residue conservation counts */
  have_tens = FALSE;
  for (icons = 0; icons < NCONSERV; icons++)
    {
      /* Check whether count is a 2-digit number */
      if (count_cons[icons] > 9)
        {
          /* Set flag */
          have_tens = TRUE;

          /* Get the tens number */
          tens = count_cons[icons] / 10;  

          /* Get the colour for this residue type */
          sprintf(tmp_string,"<FONT COLOR=%s>%d</FONT>",
                  res_cons_colour[icons],tens);
        }

      /* Otherwise, take tens position to be blank */
      else
        strcpy(tmp_string,"&nbsp;");
      
      /* Add to tens line */
      strcat(res_cons[0],tmp_string);

      /* Get the units digit */
      units = count_cons[icons] % 10;

      /* Get the colour for this residue conservation range */
      if (count_cons[icons] == 0)
        sprintf(tmp_string,"<FONT COLOR=%s>.</FONT>",
                res_cons_colour[icons]);
      else
        sprintf(tmp_string,"<FONT COLOR=%s>%d</FONT>",
                res_cons_colour[icons],units);

      /* Add onto end of current string */
      strcat(res_cons[1],tmp_string);
    }

  /* Write out the residue conservation stats */
  fprintf(file_ptr,"<TH>");
  if (have_tens == TRUE)
    fprintf(file_ptr,"%s%s%s<BR>\n",fcourier,res_cons[0],fe);
  fprintf(file_ptr,"%s%s%s</TH>\n",fcourier,res_cons[1],fe);

  /* Write out the ligands in this gap region, if any */
  if (ligand_list[0] != '\0')
    fprintf(file_ptr,"<TD><B>%s</B></TD>\n",ligand_list);
  else
    fprintf(file_ptr,"<TD>&nbsp;</TD>\n");
  fprintf(file_ptr,"</TR>\n");
  fprintf(file_ptr,"\n");
}
/***********************************************************************

write_db_gap  -  Write out the gap statistics for current gap to the
                 .dbt file

***********************************************************************/

void write_db_gap(FILE *file_ptr,int icluster,double volume,double ratio1,
                  double percen,int order1,double percen2,int order2,
                  double depth,int order3,char *ligand_list,
                  int *count_type,int *count_cons)
{
  char checked[8], fcourier[30], fs[30], fe[8];
  char tmp_string[50];
  char res_cons[2][MAX_PATTN_LEN], res_types[2][MAX_PATTN_LEN]; 

  int colour, icons, itype, len_string;
  int tens, units;
  int have_tens;

  static char count_ch[]={".123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
  
  static char *col_name[] = {
    "#FF0000", "#9900FF", "#FFD700", "#0000FF", "#00FF00",
    "#8B4513", "#FF69B4", "#339966", "#FF00FF", "#00FFFF",
    "#FF6600", "#EE82EE", "#FF9900", "#FFFFFF"
  };

  static char *res_type_colour[] = { 
    "skyblue", "red", "green", "gray", "purple", "brown", "yellow" };

  static char *res_cons_colour[] = { 
    "blue", "purple", "skyblue", "cyan", "gray", "green", "yellow",
    "orange", "pink", "red" };

  /* Initialise variables */
  len_string = strlen(count_ch);

  /* Get this cluster's colour */
  colour = (icluster - 1) % 13;

  /* Write out gap colour */
  sprintf(fs,"%s",col_name[colour]);
  fprintf(file_ptr,"%s::",fs);

  /* If one of the top 4 clusters, then show checked */
  if (icluster < 5)
    strcpy(checked,"CHECKED");
  else
    checked[0] = '\0';

  /* Write out remaining columns */
  fprintf(file_ptr,"%f::",volume);
  if (icluster == 1)
    fprintf(file_ptr,"%f::",ratio1);
  else
    fprintf(file_ptr,"0::");
  fprintf(file_ptr,"%f::",percen);
  fprintf(file_ptr,"%d::",order1);
  fprintf(file_ptr,"%f::",percen2);
  fprintf(file_ptr,"%d::",order2);
  fprintf(file_ptr,"%f::",depth);
  fprintf(file_ptr,"%d::",order3);

  /* Initialise patterns */
  res_types[0][0] = '\0';
  res_cons[0][0] = '\0';
  res_types[1][0] = '\0';
  res_cons[1][0] = '\0';

  /* Form the pattern of residue-type counts for this gap */
  have_tens = FALSE;
  for (itype = 0; itype < NRESTYPES; itype++)
    {
      sprintf(tmp_string,"%d;",count_type[itype]); 

      /* Add to tens line */
      strcat(res_types[0],tmp_string);

      /* Get the units digit */
      units = count_type[itype] % 10;

      /* Add onto end of current string */
      strcat(res_types[1],tmp_string);
    }

  /* Write out the residue type stats */
  fprintf(file_ptr,"%s",res_types[1]);

  /* Repeat for the residue conservation counts */
  have_tens = FALSE;
  for (icons = 0; icons < NCONSERV; icons++)
    {
      /* Add to tens line */
      strcat(res_cons[0],tmp_string);

      /* Add onto end of current string */
      strcat(res_cons[1],tmp_string);
    }

  /* Write out the ligands in this gap region, if any */
  if (ligand_list[0] != '\0')
    fprintf(file_ptr,"::%s",ligand_list);
  else
    fprintf(file_ptr,"::&nbsp;");
  fprintf(file_ptr,"\n");
}
/***********************************************************************

close_html_page  -  Close off the given HTML page

***********************************************************************/

void close_html_page(FILE *file_ptr,int nclusters,int profunc)
{
  /* Close off the table of gap regions */
  fprintf(file_ptr,"</TABLE>\n");

  /* Show the number of entries found */
  if (nclusters > 1)
    {
      fprintf(file_ptr,"<P ALIGN=RIGHT>\n");
      fprintf(file_ptr,"<B><FONT COLOR=BLUE>Total number of gap regions");
      fprintf(file_ptr," analysed: %d</FONT></B>\n",
              nclusters);
    }

  /* Outer table */
  fprintf(file_ptr,"<P>\n");
  fprintf(file_ptr,"<CENTER>\n");
  fprintf(file_ptr,"<TABLE WIDTH=80%% BORDER=1 CELLPADDING=5 ");
  fprintf(file_ptr,"BGCOLOR=\"#FFFFB2\">\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TD>");

  /* Show instructions for viewing in RasMol */
  fprintf(file_ptr,"<H3>Viewing in RasMol</H3>\n");
  fprintf(file_ptr,"<P>\n");
  fprintf(file_ptr,"To view the protein and its gap regions in RasMol\n");
  fprintf(file_ptr,"select the regions to be displayed from the \n");
  fprintf(file_ptr,"above table, choose the display-options\n");
  fprintf(file_ptr,"below, and then click on the\n");
  fprintf(file_ptr,"\"<FONT COLOR=RED>DISPLAY</FONT>\" button.\n");

  /* Show table of options */
  fprintf(file_ptr,"<P>\n");
  fprintf(file_ptr,"<CENTER>\n");
  fprintf(file_ptr,"<TABLE>\n");
  fprintf(file_ptr,"<TR>\n");
  fprintf(file_ptr,"<TH><FONT COLOR=BLUE>Show:</FONT></TH>\n");
  fprintf(file_ptr,"<TH ALIGN=LEFT>\n");
  fprintf(file_ptr,"<INPUT TYPE=\"radio\" NAME=\"surf\" VALUE=\" \" CHECKED>");
  fprintf(file_ptr," Binding-site(s)\n");
  fprintf(file_ptr,"</TH>\n");
  fprintf(file_ptr,"<TD ALIGN=CENTER><FONT COLOR=BLUE>&nbsp;or&nbsp;");
  fprintf(file_ptr,"</FONT></TD>\n");
  fprintf(file_ptr,"<TH ALIGN=LEFT>\n");
  fprintf(file_ptr,"<INPUT TYPE=\"radio\" NAME=\"surf\" VALUE=\"-bv\">");
  fprintf(file_ptr," Binding-surface(s)\n");
  fprintf(file_ptr,"</TH>\n");
  fprintf(file_ptr,"</TR>\n");
  fprintf(file_ptr,"\n");

  fprintf(file_ptr,"<TR>\n");
  fprintf(file_ptr,"<TH>&nbsp;</TH>\n");
  fprintf(file_ptr,"<TH><IMG WIDTH=100 ");
  if (profunc == TRUE)
    fprintf(file_ptr,"SRC=\"gif/bsite.gif\"></TH>\n");
  else
    fprintf(file_ptr,"SRC=\"/bsm/pdbsum/gif/bsite.gif\"></TH>\n");
  fprintf(file_ptr,"<TH>&nbsp;</TH>\n");
  fprintf(file_ptr,"<TH><IMG WIDTH=100 ");
  if (profunc == TRUE)
    fprintf(file_ptr,"SRC=\"gif/bsurf.gif\"></TH>\n");
  else
    fprintf(file_ptr,"SRC=\"/bsm/pdbsum/gif/bsurf.gif\"></TH>\n");
  fprintf(file_ptr,"</TR>\n");
  fprintf(file_ptr,"\n");

  fprintf(file_ptr,"<TR>\n");
  fprintf(file_ptr,"<TH>&nbsp;</TH>\n");
  fprintf(file_ptr,"<TH>\n");
  fprintf(file_ptr," <TABLE>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr," <TD>Coloured by:</TD>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr," <TD>\n");
  fprintf(file_ptr," <SPACER TYPE=\"block\" WIDTH=\"10\">\n");
  fprintf(file_ptr," <INPUT TYPE=\"radio\" NAME=\"xcol\" VALUE=\"-col 0\" "
          "CHECKED>");
  fprintf(file_ptr," gap region (as in table above)<BR>\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr," </TD>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr," </TABLE>\n");
  fprintf(file_ptr,"</TH>\n");
  fprintf(file_ptr,"<TH>&nbsp;</TH>\n");

  fprintf(file_ptr,"<TH>\n");
  fprintf(file_ptr," <TABLE>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr," <TD>Coloured by:</TD>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"\n");

  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr," <TD>\n");
  fprintf(file_ptr," <SPACER TYPE=\"block\" WIDTH=\"10\">\n");
  fprintf(file_ptr," <INPUT TYPE=\"radio\" NAME=\"col\" VALUE=\"-col 0\" "
          "CHECKED>");
  fprintf(file_ptr," gap region (as in table above)<BR>\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr," <SPACER TYPE=\"block\" WIDTH=\"10\">\n");
  fprintf(file_ptr," <INPUT TYPE=\"radio\" NAME=\"col\" VALUE=\"-col 1\">");
  fprintf(file_ptr," closest atom type<BR>\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr," <SPACER TYPE=\"block\" WIDTH=\"10\">\n");
  fprintf(file_ptr," <INPUT TYPE=\"radio\" NAME=\"col\" VALUE=\"-col 2\">");
  fprintf(file_ptr," residue type\n");
  fprintf(file_ptr," <FONT COLOR=\"BLUE\"><SUP>*</SUP></FONT><BR>\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr," <SPACER TYPE=\"block\" WIDTH=\"10\">\n");
  fprintf(file_ptr," <INPUT TYPE=\"radio\" NAME=\"col\" VALUE=\"-col 3\">");
  fprintf(file_ptr," residue conservation\n");
  fprintf(file_ptr," <FONT COLOR=\"BLUE\"><SUP>**</SUP></FONT><BR>\n");
  fprintf(file_ptr," </TD>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"</TABLE>\n");
  fprintf(file_ptr,"</TH>\n");
  fprintf(file_ptr,"</TR>\n");
  fprintf(file_ptr,"\n");

  fprintf(file_ptr,"<TR>\n");
  fprintf(file_ptr,"<TD COLSPAN=4>&nbsp;</TD>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"\n");

  fprintf(file_ptr,"<TR>\n");
  fprintf(file_ptr,"<TH COLSPAN=4><FONT COLOR=RED><INPUT TYPE=\"submit\"");
  fprintf(file_ptr,"VALUE=\"DISPLAY\"></FONT>\n");
  /* No longer have viewers.html file, so can't link to it
  if (profunc == TRUE)
    {
      fprintf(file_ptr,"     <A HREF=\"/cgi-bin/runprog.pl?prog=funchtml");
      fprintf(file_ptr,"&html_name=viewers.html&*ID*&temp=NO\">");
    }
  else
    fprintf(file_ptr,"<A HREF=\"/bsm/pdbsum/gif/viewers.html\">");
  fprintf(file_ptr,"<IMG ALIGN=BOTTOM BORDER=0\n");
  if (profunc == TRUE)
    fprintf(file_ptr,"SRC=\"gif/bubble_help.gif\"\n");
  else
    fprintf(file_ptr,"SRC=\"/bsm/pdbsum/gif/bubble_help.gif\"\n");
  fprintf(file_ptr,"ALT=\"Help!\"></A></TH>\n");
  fprintf(file_ptr,"</TR>\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"</TABLE>\n");
  fprintf(file_ptr,"</CENTER>\n");
  */
  
  /* Close off RasMol options table */
  fprintf(file_ptr,"  </TD>");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"</TABLE>\n");
  fprintf(file_ptr,"</CENTER>\n");

  /* Open table showing key to colour codes */
  fprintf(file_ptr,"<P>\n");
  fprintf(file_ptr,"<TABLE WIDTH=90%% CELLPADDING=5>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TD>");
  // fprintf(file_ptr,"<H3><FONT COLOR=WHITE>Key to residue colour ");
  fprintf(file_ptr,"<H3><FONT COLOR=#888888>Key to residue colour ");
  fprintf(file_ptr,"codes</FONT></H3>\n");
  fprintf(file_ptr,"<P>\n");

  /* Show key for residue type colourings */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"<P>\n");
  fprintf(file_ptr,"<CENTER>\n");
  fprintf(file_ptr,"<TABLE BORDER=1 WIDTH=90%% BGCOLOR=\"#FDF5E6\">\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH COLSPAN=7 ALIGN=LEFT>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"BLUE\"><SUP>*</SUP></FONT>\n");
  fprintf(file_ptr,"   <B>Residue-type colouring:</B>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"SKYBLUE\">Positive</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"RED\">Negative</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"GREEN\">Neutral</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"GRAY\">Aliphatic</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"PURPLE\">Aromatic</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"BROWN\">Pro &amp; Gly</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"YELLOW\">Cysteine</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"SKYBLUE\">H,K,R</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"RED\">D,E</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"GREEN\">S,T,N,Q</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"GRAY\">A,V,L,I,M</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"PURPLE\">F,Y,W</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"BROWN\">P,G</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"YELLOW\">C</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"</TABLE>\n");

  /* Show key for residue conservation colouring */
  fprintf(file_ptr,"<P>\n");
  fprintf(file_ptr,"<TABLE BORDER=1 WIDTH=90%% BGCOLOR=\"#FDF5E6\">\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH COLSPAN=10 ALIGN=LEFT>\n");
  fprintf(file_ptr,"   <FONT COLOR=\"BLUE\"><SUP>**</SUP></FONT>\n");
  fprintf(file_ptr,"   <B>Colouring by residue-conservation score:</B>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=BLUE>blue</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=PURPLE>purple</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=SKYBLUE>skyblue</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=CYAN>cyan</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=GRAY>grey</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=GREEN>green</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=YELLOW>yellow</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=ORANGE>orange</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=PINK>pink</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=RED>Red</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr," <TR>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=BLUE>0.0-0.1</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=PURPLE>0.1-0.2</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=SKYBLUE>0.2-0.3</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=CYAN>0.3-0.4</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=GRAY>0.4-0.5</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=GREEN>0.5-0.6</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=YELLOW>0.6-0.7</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=ORANGE>0.7-0.9</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=PINK>0.8-0.9</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr,"  <TH>\n");
  fprintf(file_ptr,"   <FONT COLOR=RED>0.9-1.0</FONT>\n");
  fprintf(file_ptr,"  </TH>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"</TABLE>\n");
  fprintf(file_ptr,"</CENTER>\n");
  fprintf(file_ptr,"\n");

  /* Close off key */
  fprintf(file_ptr,"  </TD>\n");
  fprintf(file_ptr," </TR>\n");
  fprintf(file_ptr,"</TABLE>\n");
  fprintf(file_ptr,"</CENTER>\n");

  /* Close off the form */
  fprintf(file_ptr,"</FORM>\n");

  /* Write out the closing records */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"</BODY>\n");
  fprintf(file_ptr,"</HTML>\n");
}
/***********************************************************************

write_speedfill_arr  -  Write out the grid points to the speedfill.arr
                        file

***********************************************************************/

void write_speedfill_arr(FILE *file_ptr,int *grid,int npoints[3],
                         double valmin[3],double valmax[3],
                         double grid_sep,int ngrid_points)
{
  /* Write out the number of points, grid separation and maximum and
     minimum coords */
  fwrite(npoints,3 * sizeof(int),1,file_ptr);

  /* Write out the grid spacing */
  fwrite(&grid_sep,sizeof(double),1,file_ptr);

  /* Write out the coords of the minimum corner point */
  fwrite(valmin,3 * sizeof(double),1,file_ptr);

  /* Repeat for the maximum corner point */
  fwrite(valmax,3 * sizeof(double),1,file_ptr);

  /* Write out the array itself */
  fwrite(grid,ngrid_points * sizeof(int),1,file_ptr);

  /* Close the array file */
  fclose(file_ptr);
}
/***********************************************************************

get_closest_sites  -  Read through the binding-site analysis file to
                      find the top few binding sites that most closely
                      resemble this one

***********************************************************************/

void get_closest_sites(int cluster,double volume,double percen,
                       double percen2,double depth,double eigen_value2,
                       double eigen_value3,int *count_type,
                       char *pdbsum_data_dir)
{
  char file_name[FILENAME_LEN];
  char input_line[LINELEN + 1];

  int nsites;
  int have_file;

  FILE *file_ptr;

  /* Initialise variables */
  have_file = TRUE;
  nsites = 0;

  /* Form name of the binding sites file */
  strcpy(file_name,pdbsum_data_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"binding_sites.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* If this is a binding-site record, then process */
          if (input_line[0] == '+' && (input_line[1] == '+' ||
                                       input_line[1] == '-'))
            {
              /* Extract the data from this line */
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }
}
/***********************************************************************

calc_exp_table -  Calculate a lookup table of exponential distribution
                  values for a range of possible sphere radii

***********************************************************************/

void calc_exp_table(int **exp_lookup_table,int *nr,int *nv,
                    double grid_sep,struct parameters params)
{
  int ipos, iradius, nradii, ivalue, nvalues;
  int *exp_table;

  double constant, dist_sqrd, grid_value, radius;

  /* Calculate the number of radii to be represented by the look-up
     tables */
  nradii
    = (int) ((params.max_gap_radius - params.min_gap_radius)
             / RADSTEP) + 2;

  /* Calculate distance increments for each table */
  dist_sqrd
    = (params.max_gap_radius + grid_sep)
    * (params.max_gap_radius + grid_sep);
  nvalues = (int) (3.0 * dist_sqrd / DISTEP) + 2;

  /* Create the exponential look-up tables as one large array */
  exp_table = (int *) malloc(nradii * nvalues * sizeof(int));
  if (exp_table == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for exp table\n");
      exit(1);
    }

  /* Initialise */
  ipos = 0;

  /* Loop over all the radii */
  for (iradius = 0; iradius < nradii; iradius++)
    {
      /* Calculate the corresponding radius */
      radius = params.min_gap_radius + iradius * RADSTEP;

      /* Calculate the scaling factor for this radius */
      constant =  (- log( 2.0) / (radius * radius));

      /* Loop over all the table values */
      for (ivalue = 0; ivalue < nvalues; ivalue++)
        {
          /* Calculate corresponding distance squared */
          dist_sqrd = ivalue * DISTEP;

          /* Calculate the current exponential value and store in
             the table */
          grid_value
            = 2 * CONTOUR_LEVEL *  exp(constant * dist_sqrd);
          exp_table[ipos] = (int) (GRID_SCALE_FACTOR * grid_value);

          /* Increment table position */
          ipos++;
        }
    }

  /* Return table parameters */
  *nr = nradii;
  *nv = nvalues;

  /* Return the pointer to the start of the lookup table */
  *exp_lookup_table = exp_table;
}
/***********************************************************************

sort_by_volume  -  Sort clusters according to the volumes they occupy

***********************************************************************/

void sort_by_volume(double *volume,int *index,int nclusters)
{
  int icluster, ipos, jcluster;
  int swapped;

  /* Sort the index-array using a simple bubble-sort such that the 
     clusters are ordered in descending order of their volume */

  /* Initialise swap flag */
  swapped = TRUE;

  /* Keep looping while any elements are still being swapped */
  while (swapped == TRUE)
    {
      /* Unset the swap flag */
      swapped = FALSE;

      /* Loop through all the clusters */
      for (ipos = 0; ipos < nclusters - 1; ipos++)
        {
          /* Get the clusters pointed to by this index element and the
             next one */
          icluster = index[ipos];
          jcluster = index[ipos + 1];

          /* If their volumes are the wrong way round, then need to
             swap them */
          if (volume[icluster] < volume[jcluster])
            {
              /* Swap the index elements */
              index[ipos] = jcluster;
              index[ipos + 1] = icluster;

              /* Set flag to indicate that a swap took place in this loop */
              swapped = TRUE;
            }
        }
    }
}
/***********************************************************************

contour_clusters  -  Write out the surface of each cluster as a
                     contour-map (in PDB format) for reading into RasMol

***********************************************************************/

int contour_clusters(FILE *outfile_ptr,FILE *scriptfile_ptr,
                     FILE *pmlfile_ptr,
                     FILE *r3d_file_ptr,int *grid,int *save_grid,
                     int npoints[3],double valmin[3],double valmax[3],
                     double grid_sep,struct sphere *first_sphere_ptr,
                     int nclusters,struct parameters params,
                     struct vertex **fst_vertex_ptr,int *atom_no,
                     struct atom *first_lig_atom_ptr,char pdb_code[5],
                     char chain[MAX_CHAINS],int nchain,
                     char *pdb_name,char *pdbsum_data_dir,FILE *logfile_ptr,
                     FILE *html_file_ptr,FILE *db_file_ptr,
                     FILE *at_file_ptr,FILE *surfnet_file_ptr,
                     FILE *arr_file_ptr,
                     int *grid_atom_map,struct atom **atom_index_ptr,
                     int natoms,char *pdbsum_dir,char *www_pdbsum_dir,
                     struct atom *first_atom_ptr,
                     struct residue *first_residue_ptr,
                     int *exp_lookup_table,int nradii,int nvalues,
                     double *volume,struct limits atom_limits,
                     double rascale,int filter_method,
                     float *delphi_potential,int ndelphi_grid,
                     double delphi_grid_sep,double *delphi_valmin)
 {
  char ch, code[5], inchn, star;
  char ligand_list[MAX_TOP_GAPS][LISTLEN];
  char ligand_list2[MAX_TOP_GAPS][LISTLEN];
  char protein_name[NAMLEN];

  int nres[NCONSERV], ntop_clusters;
  int ipoint, ngrid_points, nligatoms[MAX_TOP_GAPS], nout;
  int narsize, npolygons, nvertices, tot_vertices;
  int atom_count, i, ipos, itype;
  int order1[MAX_TOP_GAPS], order2[MAX_TOP_GAPS], order3[MAX_TOP_GAPS];
  int order4[MAX_TOP_GAPS];
  int cluster, colour, icluster, wclusters;
  int atom_number, first_atom_number, last_atom_number;
  int contour_level, ncentroids;
  int done, wanted;
  int count_type[MAX_TOP_GAPS][NRESTYPES];
  int count_cons[MAX_TOP_GAPS][NCONSERV];

  static char *colour_name[] = {
    "red", "purple", "yellow", "blue", "green",
    "brown", "pink", "olivegreen", "magenta", "cyan",
    "limegreen", "violet", "orange", "white"
  };

  int *index;

  double depth[MAX_TOP_GAPS];
  double pca3[MAX_TOP_GAPS];
  double percen[MAX_TOP_GAPS], percen2[MAX_TOP_GAPS], ratio1;
  double eigen_value[MAX_TOP_GAPS][3];

  struct atom *dummy_atom_ptr, *lst_atom_ptr, *fst_atom_ptr;
  struct grid_limits glimits;
  struct limits vertex_limits;
  struct polygon *first_polygon_ptr;
  struct residue *dummy_residue_ptr, *lst_residue_ptr, *fst_residue_ptr;
  struct sphere *sphere_ptr, *start_sphere_ptr;
  struct vertex *first_vertex_ptr, *last_vertex_ptr;
  struct vertex *very_first_vertex_ptr, *very_last_vertex_ptr;
  struct vertex *first_sorted_vertex_ptr[3];

  struct vertex **vert_atom_index_ptr;
  struct vertex **vertex_index_ptr, **vertex_index3_ptr;

  /* Initialise variables */
  atom_count = 0;
  inchn = ' ';
  protein_name[0] = '\0';
  ntop_clusters = params.ntop_clusters;
  fprintf(logfile_ptr,"Contouring gap regions ...\n");
  fprintf(logfile_ptr,"+\n");
  fprintf(logfile_ptr,"+ Gap-region summary for ");
  fflush(logfile_ptr);
  strcpy(code,pdb_code);
  if (pdb_code[0] != '\0')
    fprintf(logfile_ptr,"%s\n",pdb_code);
  else
    {
      fprintf(logfile_ptr,"%s\n",pdb_name);

      /* Find location of the dot */
      ipos = 0;
      for (i = 0; i < (int) strlen(pdb_name); i++)
        {
          if (pdb_name[i] == '.')
            ipos = i;
        }

      /* Strip out the PDB code (assumed to be last 4 chars before the dot) */
      if (ipos > 3)
        {
          strncpy(code,pdb_name+(ipos - 4),4);
          code[4] = '\0';
        }
      else
        strcpy(code,"  - ");
    }
  fflush(logfile_ptr);

  /* If writing Raster3D file, output the material properties of the
     cluster surfaces */
  if (params.generate_r3d == TRUE)
    {
      fprintf(r3d_file_ptr,"8\n");
      fprintf(r3d_file_ptr,"25 0.40 1.0 1.0 1.0 0.1 0 0 0 0\n");
    }

  /* Create a dummy residue record for the dummy atom */
  fst_residue_ptr = NULL;
  lst_residue_ptr = NULL;
  dummy_residue_ptr
    = create_residue_record(&fst_residue_ptr,&lst_residue_ptr,"XXX",
                            "   1 ",'X');

  /* Create a dummy atom for vertices to link to */
  fst_atom_ptr = NULL;
  lst_atom_ptr = NULL;
  dummy_atom_ptr
    = create_atom_record(&fst_atom_ptr,&lst_atom_ptr,dummy_residue_ptr,
                         0," DUM","X",0.0,0.0,0.0,0.0,0.0);

  /* If generating an HTML page showing the binding-site analysis, write
     header details */
  if (params.generate_html == TRUE)
    {
      /* Pick up the protein's name from the PDBsum directory */
      get_protein_name(pdb_code,pdbsum_dir,protein_name);

      /* Write header details to HTML file */
      write_html_headers(html_file_ptr,pdb_code,chain,nchain,protein_name,
                         www_pdbsum_dir,params);
    }

  /* If called from a WWW page, then don't need the detailed analysis */
  if (params.from_www == FALSE)
    {
      /* Print column headings */
      fprintf(logfile_ptr,"+\n");
      fprintf(logfile_ptr,"+ PDB                    R1   Accessible   Bur");
      fprintf(logfile_ptr,"ied      Ave.");
      fprintf(logfile_ptr,"    Principal \n");
      fprintf(logfile_ptr,"+ code  Gap     Volume  ratio  vertices   vert");
      fprintf(logfile_ptr,"ices    depth");
      fprintf(logfile_ptr,"    components             +ve -ve Neu Ali Aro");
      fprintf(logfile_ptr," P/G Cys  ");
      fprintf(logfile_ptr,"Ligands\n");
      fprintf(logfile_ptr,"+ ----  ---    -------  -----  --------   ----");
      fprintf(logfile_ptr,"----    -----");
      fprintf(logfile_ptr,"    -------------------    --- --- --- --- ---");
      fprintf(logfile_ptr," --- ---  ");
      fprintf(logfile_ptr,"-------\n");
      fprintf(logfile_ptr,"+\n");
      fflush(logfile_ptr);
    }

  /* Initialise counts */
  atom_number = *atom_no;
  tot_vertices = 0;

  /* Determine how many clusters are to be written out */
  wclusters = ntop_clusters;
  if (ntop_clusters > nclusters)
    wclusters = nclusters;

  /* Initialise pointers */
  very_first_vertex_ptr = NULL;
  very_last_vertex_ptr = NULL;
  first_vertex_ptr = NULL;

  /* Calculate overall size of grid map */
  ngrid_points = npoints[0] * npoints[1] * npoints[2];

  /* Calculate integer version of contour level */
  contour_level = (int) (CONTOUR_LEVEL * GRID_SCALE_FACTOR);

  /* Loop over all the points in the grid to initialise such that only
     the masked points are retained */
  for (ipoint = 0; ipoint < ngrid_points; ipoint++)
    {
      /* If grid-point not masked, then set to empty */
      if (grid[ipoint] != MASKED)
        grid[ipoint] = EMPTY;

      /* If array to be saved, copy value to save-array */
      if (params.save_grid == TRUE)
        save_grid[ipoint] = grid[ipoint] - 1;
    }

  /* Get pointer to the first sphere */
  start_sphere_ptr = first_sphere_ptr;

  /* Loop over all the sphere clusters to be contoured */
  for (icluster = 0; icluster < wclusters; icluster++)
    {
      /* Show cluster being processed */
      printf("\n");
      printf("Cluster %d\n",icluster);

      /* Determine whether this cluster is wanted */
      if (params.show_gap[icluster] == TRUE)
        wanted = TRUE;
      else
        wanted = FALSE;

      /* Reset grid array */
      if (wanted == TRUE)
        reset_grid_array(grid,ngrid_points,contour_level);

      /* Initialise first and last atom numbers for this cluster */
      first_atom_number = atom_number + 1;
      last_atom_number = atom_number;

      /* Initialise lists of ligands in this gap region */
      ligand_list[icluster][0] = '\0';
      ligand_list2[icluster][0] = '\0';

      /* Get pointer to the first sphere */
      sphere_ptr = start_sphere_ptr;

      /* Initialise flag */
      done = FALSE;
      nout = 0;

      /* Loop through all spheres for the current cluster */
      while (sphere_ptr != NULL && done == FALSE)
        {
          /* Get the cluster number */
          cluster = sphere_ptr->cluster;

          /* If sphere belongs to the required cluster, then use it to
             update the grid */
          if (cluster == icluster + 1)
            {
              /* Update the grid map with current sphere, if want
                 this cluster */
              if (wanted == TRUE)
                {
                  update_map(grid,npoints,valmin,valmax,grid_sep,
                             sphere_ptr->coords,sphere_ptr->radius,
                             exp_lookup_table,nradii,nvalues,contour_level,
                             params);
                  nout++;
                }
            }

          /* Otherwise, save current sphere-pointer as the start of the
             next cluster */
          else
            {
              /* Save current sphere-pointer as the start of the
                 next cluster */
              start_sphere_ptr = sphere_ptr;

              /* Set flag to indicate that this cluster is now done */
              done = TRUE;
            }

          /* Get pointer to the next sphere */
          sphere_ptr = sphere_ptr->next_sphere_ptr;
        }

      /* Process grid-map generated for this cluster only if cluster is
         wanted */
      if (wanted == TRUE)
        {
          /* Calculate the volume of the gap region */
          //fprintf(logfile_ptr," Process grid-map generated for this "
          //        "cluster ...\n");
          volume[icluster]
            = calculate_gap_volume(grid,npoints,valmin,valmax,grid_sep,
                                   contour_level,ngrid_points);

          /* Calculate the non-zero limits of the grid array */
          calc_grid_limits(grid,npoints,&glimits,contour_level);
  
          /* Get all contouring points in the current map */
          //fprintf(logfile_ptr,"  get_contour_points\n");
          nvertices
            = get_contour_points(grid,grid_atom_map,npoints,
                                 atom_index_ptr,natoms,valmin,valmax,
                                 grid_sep,&first_vertex_ptr,
                                 first_sorted_vertex_ptr,contour_level,
                                 glimits,icluster,params);
          tot_vertices = tot_vertices + nvertices;
          narsize = nvertices;

          /* If colouring by electrostatic potential, get colour for
             each vertex */
          if (params.vertex_colouring == ELECTROSTATIC)
            calc_vertex_energies(first_vertex_ptr,delphi_potential,
                                 ndelphi_grid,delphi_grid_sep,
                                 delphi_valmin);

          /* Create index array for vertex pointers */
          create_vertex_array(first_vertex_ptr,&vertex_index_ptr,
                              nvertices);

          /* Calculate degree of each vertex's burial, if required */
          if (params.from_www == FALSE ||
              params.buried_vertices_only == TRUE ||
              params.vertex_colouring != STANDARD)
            {
              /* Determine how many of the vertices are close to the
                 protein, and how many are on the solvent side */
              //fprintf(logfile_ptr," calc_vertex_burial\n");
              percen[icluster]
                = calc_vertex_burial(grid,npoints,valmin,valmax,
                                     grid_sep,ngrid_points,contour_level,
                                     first_vertex_ptr,&percen2[icluster],
                                     params);
            }

          /* If array to be saved, update the save array with the points
             from the current cluster */
          if (params.save_grid == TRUE)
            save_cluster_points(grid,save_grid,ngrid_points,icluster,
                                contour_level);

          /* If Raster3D output required, generate polygons and write
             out .r3d file */
          if (params.generate_r3d == TRUE)
            {
              /* Generate all the polygons to tile the given contour
                 surface */
              fprintf(logfile_ptr," generate_polygons...\n");
              npolygons
                = generate_polygons(grid,npoints,atom_limits.valmin,
                                    atom_limits.valmax,grid_sep,
                                    contour_level,glimits,
                                    vertex_index_ptr,nvertices,
                                    &first_polygon_ptr,dummy_atom_ptr,
                                    params.buried_vertices_only);

              /* Order each polygon's vertices so that they are in an
                 anti-clockwise sense as viewed from outside */
              // order_polygon_vertices(first_polygon_ptr,npoints,
              //                       atom_limits.valmin,atom_limits.valmax,
              //                       grid_sep);

              /* Calculate averaged normals at each of the vertices */
              fprintf(logfile_ptr," calculate_vertex_normals ...\n");
              calculate_vertex_normals(first_vertex_ptr,&vertex_limits);

              /* Get the colour for this cluster */
              colour = icluster % 14;

              /* Write out the polygons to the Raster3D output file */
              fprintf(logfile_ptr," writing %d Raster3d polygons ...\n",
                      npolygons);
              write_raster3d_polygons(r3d_file_ptr,first_polygon_ptr,
                                      atom_limits.c_of_g,rascale,inchn,
                                      params.vertex_colouring,
                                      params.conservation_colours,
                                      colour_name[colour]);
            }

          /* Apply the required vertex-filtering method to remove vertices
             to reduce the number of bonds RasMol draws when connecting
             them */
          if (params.filter_vertices == TRUE && filter_method == OLD)
            {
              /* Create the 3 index arrays that will hold vertex addresses
                 sorted by their x-, y- and z-coordinates */
              create_vertex_index_arrays(first_sorted_vertex_ptr,
                                         &vertex_index3_ptr,nvertices);

              /* If showing binding regions, filter out vertices that are
                 too close to one another */
              nvertices
                = filterx_vertices(first_vertex_ptr,nvertices,
                                   vertex_index3_ptr,
                                   params.buried_vertices_only);
            }

          /* New method uses kd-trees */
          else if (params.filter_vertices == TRUE)
            {
              /* Process all the vertices, knocking out any that are
                 too close to one another */
              // fprintf(logfile_ptr," filter_vertices ... %d\n",nvertices);
              nvertices
                = filter_vertices(first_vertex_ptr,nvertices,
                                  vertex_index_ptr,
                                  params.buried_vertices_only,
                                  params.save_vertices);
            }

          /* Create index array for vertex pointers */
          create_vertex_array(first_vertex_ptr,&vertex_index_ptr,
                              nvertices);

          /* If called from a WWW page, then don't need the detailed
             analysis */
          if (params.from_www == FALSE)
            {
              /* Determine which ligands, if any, lie within the current gap
                 region */
              nligatoms[icluster]
                = get_ligands_in_gap(grid,npoints,valmin,valmax,grid_sep,
                                     first_lig_atom_ptr,contour_level,
                                     ngrid_points,ligand_list[icluster],
                                     ligand_list2[icluster]);

              /* Use the buried/exposed vertex information to calculate the
                 average "depth" of the current gap region */
              depth[icluster]
                = calculate_gap_depth(first_vertex_ptr,narsize,
                                      vertex_index3_ptr,vertex_index_ptr,
                                      filter_method);

              /* Calculate the principal components */
              principal_components(first_vertex_ptr,eigen_value[icluster]);
            }

          /* If colouring by atom details, then sort vertices by closest
             atom */
          sort_vertices(first_vertex_ptr,&vert_atom_index_ptr,
                        &index,nvertices,params);

          /* If have a previous cluster, link its list of vertices to
             the current list */

          /* Skip unless outputting for RasMol */
          if (params.no_rasmol == FALSE)
            {
             /* Write out the vertices to the RasMol script file */
              ncentroids
                = write_vertices(outfile_ptr,icluster,&atom_number,
                                 vert_atom_index_ptr,index,nvertices,
                                 nres,params);
              last_atom_number = atom_number;
            }

          /* Analyse residue types and conservation scores for the current
             gap volume */
          analyse_vertices(vert_atom_index_ptr,index,nvertices,
                           first_atom_ptr,first_residue_ptr,
                           count_type[icluster],count_cons[icluster]);

          /* If generating HTML data, write out atoms adjacent to the
             current cleft to the .dbt file */
          if (params.generate_html == TRUE || params.pdb_type == PDBSUM1)
            atom_count
              = write_cleft_atoms(icluster + 1,first_atom_ptr,at_file_ptr,
                                  atom_count);

          /* Relink vertices for later writing out to speedfill.depth
             file */
          nvertices
            = relink_vertices(&first_vertex_ptr,&last_vertex_ptr,
                              vert_atom_index_ptr,index,nvertices,icluster);

          /* Incorporate verices into full list from previous
             clusters */
          if (nvertices > 0)
            {
              /* If this is the very first cluster, then store pointer */
              if (very_first_vertex_ptr == NULL)
                very_first_vertex_ptr = first_vertex_ptr;

              /* Add link from previous vertex to the current one */
              if (very_last_vertex_ptr != NULL)
                very_last_vertex_ptr->next_vertex_ptr = first_vertex_ptr;
              very_last_vertex_ptr = last_vertex_ptr;
            }

          /* Write out RasMol script commands to display current cluster */
          if (params.no_rasmol == FALSE &&
              last_atom_number > first_atom_number)
            {
              /* Write RasMol script commands */
              write_script_clusters(scriptfile_ptr,pmlfile_ptr,icluster,
                                    first_atom_number,last_atom_number,
                                    volume[icluster],ncentroids,nres,
                                    params);

              /* If writing out surfnet.par file, then write out a line
                 for current cluster */
              if (params.write_surfnet_input == TRUE)
                write_surfnet_gap(surfnet_file_ptr,icluster,first_atom_number,
                                  last_atom_number);
            }
        }
    }

  /* If called from a WWW page, then don't need the detailed analysis */
  if (params.from_www == FALSE)
    {
      /* Calculate the volume ratio (ie largest gap region to
         second-largest) */
      if (volume[1] > 0.0)
        ratio1 = volume[0] / volume[1];
      else
        ratio1 = 0.0;

      /* If there's a ligand in the largest gap region, mark all entries */
      if (nligatoms[0] > 4)
        star = '*';
      else
        star = ' ';

      /* Sort the various values to order from highest to lowest */
      order_results(percen,wclusters,order1);
      order_results(percen2,wclusters,order2);
      order_results(depth,wclusters,order3);

      /* Extract the third component from the PCA analysis, and order these */
      for (icluster = 0; icluster < wclusters; icluster++)
        pca3[icluster] = eigen_value[icluster][2];
      order_results(pca3,wclusters,order4);

      /* Loop over all the clusters again to write out the details */
      for (icluster = 0; icluster < wclusters; icluster++)
        {
          /* Read through the binding-site analysis file to find the
             top few binding sites that most closely resemble this one */
          get_closest_sites(icluster,volume[icluster],percen[icluster],
                            percen2[icluster],depth[icluster],
                            eigen_value[icluster][1],eigen_value[icluster][2],
                            count_type[icluster],pdbsum_data_dir);

          /* If ligand in this gap, set marker to "+" */
          if (nligatoms[icluster] > 4)
            ch = '+';
          else
            ch = '-';

          /* Summarise the data on this gap region */
          fprintf(logfile_ptr,
                  "+%c%s%c %2d. %10.2f %6.2f %7.2f%% %2d%7.2f%% %2d%7.2f ",
                  ch,code,star,icluster + 1,volume[icluster],ratio1,
                  percen[icluster],order1[icluster] + 1,
                  percen2[icluster],order2[icluster] + 1,
                  depth[icluster]);
          fprintf(logfile_ptr,
                  "%2d %3.1f %6.4f %6.4f %2d  ",order3[icluster] + 1,
                  eigen_value[icluster][0],eigen_value[icluster][1],
                  eigen_value[icluster][2],order4[icluster] + 1);
          for (itype = 0; itype < NRESTYPES; itype++)
            fprintf(logfile_ptr,"%3d ",count_type[icluster][itype]);

          if (ligand_list[icluster][0] != '\0')
            fprintf(logfile_ptr,"  %s\n",ligand_list[icluster]);
          else
            fprintf(logfile_ptr,"\n");

          /* If writing out HTML file, then add to that */
          if (params.generate_html == TRUE)
            write_html_gap(html_file_ptr,icluster + 1,volume[icluster],ratio1,
                           percen[icluster],order1[icluster] + 1,
                           percen2[icluster],order2[icluster] + 1,
                           depth[icluster],order3[icluster] + 1,
                           ligand_list2[icluster],count_type[icluster],
                           count_cons[icluster]);

          /* Write out clefts data to the .dbt file */
          if (params.generate_html == TRUE || params.pdb_type == PDBSUM1)
            write_db_gap(db_file_ptr,icluster + 1,volume[icluster],ratio1,
                         percen[icluster],order1[icluster] + 1,
                         percen2[icluster],order2[icluster] + 1,
                         depth[icluster],order3[icluster] + 1,
                         ligand_list2[icluster],count_type[icluster],
                         count_cons[icluster]);
        }

      /* Close off summary print */
      fprintf(logfile_ptr,"+\n");
      fprintf(logfile_ptr,
              "+ ----------------------------------------------------");
      fprintf(logfile_ptr,"-------------------------\n");
      fprintf(logfile_ptr,
              "+  * signifies that the largest gap region contains a");
      fprintf(logfile_ptr," ligand\n");

      /* Close of cluster data in .dbt file */
      if (params.generate_html == TRUE || params.pdb_type == PDBSUM1)
        fprintf(db_file_ptr,"##END_DATA\n");
    }

  /* If generating an HTML page, close it off */
  if (params.generate_html == TRUE)
    close_html_page(html_file_ptr,wclusters,params.profunc);

  /* If saving the grid array, write out to speedfill.arr file */
  if (params.save_grid == TRUE)
    write_speedfill_arr(arr_file_ptr,save_grid,npoints,valmin,valmax,
                        grid_sep,ngrid_points);

  /* Return vertex pointer and current atom-number */
  *fst_vertex_ptr = very_first_vertex_ptr;
  *atom_no = atom_number;

  /* Return number of vertices stored */
  return(tot_vertices);
}
/***********************************************************************

write_script_tail  -  Write all the final records to the output RasMol
                      script file

***********************************************************************/

void write_script_tail(FILE *outfile_ptr)
{
  /* Close off script */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"select all\n");

  fprintf(outfile_ptr,"echo \" \"\n");
  fprintf(outfile_ptr,"exit\n");
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_pymol_tail  -  Write all the final lines to the PyMOL script

***********************************************************************/

void write_pymol_tail(FILE *outfile_ptr,struct parameters params)
{
  char png_name[FILENAME_LEN];

  /* Close off script */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"rebuild\n");
  fprintf(outfile_ptr,"zoom all, complete=1, buffer=2\n");
  fprintf(outfile_ptr,"\n");

  /* If not PDBsum1, add the rendering commands */
  if (params.pdb_type != PDBSUM1 && params.pymol_www == FALSE)
    {
      /* Form name of output .png file */
      strcpy(png_name,"clefts.png");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"# Render %s\n",png_name);

      /* Write out ray-trace command */
      fprintf(outfile_ptr,"ray 240, 240\n");
      fprintf(outfile_ptr,"viewport 240, 240\n");

      /* Write out command to generate png file */
      fprintf(outfile_ptr,"png %s\n",png_name);
    }
}
/***********************************************************************

insert_pdb_file  -  Insert the given PDB file

***********************************************************************/

int insert_pdb_file(FILE *outfile_ptr,char *file_name,char *pdb_code)
{
  char end_char, input_line[LINELEN + 1];

  int len, nout;

  FILE *file_ptr;

  /* Initialise variables */
  nout = 0;
  end_char = '\\';

  /* Open the input file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr == NULL)
    {
      printf("*** Warning. Unable to open file [%s] \"\n",file_name);
      return(nout);
    }

  /* Start file insertion */
  fprintf(outfile_ptr,"cmd.read_pdbstr(\"\"\"");

  /* Loop through the file, copying to output file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* If a PDB file, ony want ATOM, HETATM, TER and CONECT records */
      if (!strncmp(input_line,"HEADER",6) ||
          !strncmp(input_line,"ATOM  ",6) ||
          !strncmp(input_line,"HETATM",6) ||
          !strncmp(input_line,"TER   ",6) ||
          !strncmp(input_line,"CONECT",6))
        {
          /* If line too long, then truncate */
          if (len > 66)
            input_line[66] = '\0';
      
          /* Write the line out */
          fprintf(outfile_ptr,"%s%c\n",input_line,end_char);

          /* Increment count of lines output */
          nout++;
        }
    }

  /* Close off file insertion */
  fprintf(outfile_ptr,"MASTER%c\n",end_char);
  fprintf(outfile_ptr,"END%c\n",end_char);
  fprintf(outfile_ptr,"\"\"\",\"%s\")\n",pdb_code);

  /* Close the input file */
  fclose(file_ptr);

  /* Return number of lines written out */
  return(nout);
}
/***********************************************************************

rewrite_tmp_pml  -  Rewrite the tmp.pml file as clefts.pml

***********************************************************************/

void rewrite_tmp_pml(char *tmp_pml,char *clefts_pml,int relinks,
                     char *clefts_pdb,char *pdb_code)
{
  char input_line[LINELEN + 1];

  FILE *file_ptr, *outfile_ptr;

  /* Open the clefts.pml output file */
  outfile_ptr = open_output_file(clefts_pml,TRUE);

  /* Open the input file */
  if ((file_ptr = fopen(tmp_pml,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the line at the last non-blank character */
          string_chomp(input_line);

          /* If this is the load line, keep or replace by inline PDB file */
          if (!strncmp(input_line,"load ",5))
            {
              /* If have relative links, insert the clefts.pdb file at
                 this point */
              if (relinks == TRUE)
                {
                  /* Insert the PDB file */
                  insert_pdb_file(outfile_ptr,clefts_pdb,pdb_code);
                  
                  /* Blank out the input line */
                  input_line[0] = '\0';
                }
            }

          /* Write the line out */
          fprintf(outfile_ptr,"%s\n",input_line);
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show warning message */
  else
    {
      printf("*** Warning. Unable to open input file [%s]\n",
             tmp_pml);
    }

  /* Close the output file */
  fclose(outfile_ptr);
}
/***********************************************************************

close_surfnet_par  -  Write the final records to the surfnet.par file
                      and close

***********************************************************************/

void close_surfnet_par(FILE *file_ptr)
{
  /* Write the SURFNET parameters */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"SURFNET\n");
  fprintf(file_ptr,"  0.0   0.0   0.0\n");
  fprintf(file_ptr,"-1000 10000   \n");
  fprintf(file_ptr,"6.0 6.0 6.0 \n");
  fprintf(file_ptr,"0.80 \n");
  fprintf(file_ptr,"1.0  \n");
  fprintf(file_ptr,"4.0  \n");
  fprintf(file_ptr,"Q  \n");
  fprintf(file_ptr,"Y  \n");
  fprintf(file_ptr,"1.0\n");
  fprintf(file_ptr,"N  \n");
  fprintf(file_ptr,"Y   \n");
  fprintf(file_ptr,"10.0 \n");
  fprintf(file_ptr,"Y   \n");
  fprintf(file_ptr,"1.4 \n");
  fprintf(file_ptr,"1.8 \n");
  fprintf(file_ptr,"4  \n");
  fprintf(file_ptr,"\n");

  /* Write the SURFACE keyword */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"SURFACE\n");

  /* Write the SURFPLOT parameters */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"SURFPLOT\n");
  fprintf(file_ptr,"B      \n");
  fprintf(file_ptr,"0          \n");
  fprintf(file_ptr,"0            \n");
  fprintf(file_ptr,"N             \n");
  fprintf(file_ptr,"Y             \n");
  fprintf(file_ptr,"Y             \n");
  fprintf(file_ptr,"0.20          \n");
  fprintf(file_ptr,"N             \n");
  fprintf(file_ptr,"N             \n");
  fprintf(file_ptr," 0.0          \n");
  fprintf(file_ptr," 0.0          \n");
  fprintf(file_ptr," 0.0          \n");
  fprintf(file_ptr," 1.0000   0.0000    0.0000    \n"); 
  fprintf(file_ptr," 0.0000   1.0000    0.0000    \n");
  fprintf(file_ptr," 0.0000   0.0000    1.0000     \n");
  fprintf(file_ptr,"  50.0 100.0 100.0    \n");
  fprintf(file_ptr," 40 32         \n");
  fprintf(file_ptr,"N              \n");
  fprintf(file_ptr,"\n");

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

close_and_reopen  -  Close given RasMol output files and re-open in
                     read-mode

***********************************************************************/

FILE *close_and_reopen(FILE *file_ptr,char out_name[FILENAME_LEN])
{
  FILE *newfile_ptr;

  /* Close the currently-open file */
  fclose(file_ptr);

  /* Re-open in read-mode */
  if ((newfile_ptr = fopen(out_name,"r")) ==  NULL)
    {
      printf("\n*** Unable to re-open output RasMol output file [%s]\n",
             out_name);
      exit(1);
    }

  /* Return the file pointer */
  return(newfile_ptr);
}
/***********************************************************************

transfer_rasmol_file  -  Transfer all the records from the given RasMol
                         file into the final file

***********************************************************************/

void transfer_rasmol_file(FILE *file_ptr,FILE *output_file_ptr)
{
  char input_line[LINELEN + 1];

  int nlines;

  /* Initialise counts */
  nlines = 0;

  /* Loop while reading in records from the input file */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment count of lines read in */
      nlines++;

      /* Write out */
      fprintf(output_file_ptr,"%s",input_line);
    }

  /* Close and delete the input file */
  fclose(file_ptr);
}
/***********************************************************************

delete_temp_files  -  Delete the temporary RasMol coordinate and script
                      files

***********************************************************************/

void delete_temp_files(char *file1,char *file2,int dos)
{
  char command[FILENAME_LEN + 1];

  /* Form the command to delete the first file */
  sprintf(command,"%s %s",RM_COMMAND,file1);
  call_system(command,dos);

  /* Repeat for the second file */
  sprintf(command,"%s %s",RM_COMMAND,file2);
  call_system(command,dos);
}
/***********************************************************************

write_atom_depths  -  Write out the atom depths and the clusters (if any)
                      to which they belong

***********************************************************************/

void write_atom_depths(FILE *outfile_ptr,struct atom *first_atom_ptr,
                       struct atom *first_lig_atom_ptr)
{
  static char *keyword[2] = { "ATOM  ", "HETATM" };

  int atom_number, loop;

  double bvalue, occupancy;

  struct atom *atom_ptr;

  /* Initialise variables */
  atom_number = 0;
  occupancy = 1.0;

  /* Loop twice, once to write out the atom coordinates, and then
     to write out the ligand coordinates */
  for (loop = 0; loop < 2; loop++)
    {
      /* Get pointer to the first atom or ligand-atom, depending on which
         loop this is */
      if (loop == 0)
        atom_ptr = first_atom_ptr;
      else
        atom_ptr = first_lig_atom_ptr;

      /* Loop through all atoms to be written out */
      while (atom_ptr != NULL)
        {
          /* Write out only if atom is within the grid-array bounds */
          if (loop == 0 || atom_ptr->residue_ptr->in_range == TRUE)
            {
              /* Calculate atom's average depth */
              if (atom_ptr->ndepth != 0)
                atom_ptr->depth = atom_ptr->depth / atom_ptr->ndepth;

              /* Store depth in B-value field */
              bvalue = atom_ptr->depth;

              /* Store atom's cluster in occupancy field */
              occupancy = atom_ptr->icluster + 1;

              /* Get the atom-number, checking for overflow */
              atom_number++;
              if (atom_number > 99999)
                atom_number = 99999;

              /* Store the atom number */
              atom_ptr->atom_number = atom_number;

              /* Write out the atom record */
              fprintf(outfile_ptr,"%s%5d ",keyword[loop],atom_number);
              fprintf(outfile_ptr,"%4s ",atom_ptr->atom_name);
              fprintf(outfile_ptr,"%3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                      atom_ptr->residue_ptr->res_name,
                      atom_ptr->residue_ptr->chain,
                      atom_ptr->residue_ptr->res_num,
                      atom_ptr->coords[0],atom_ptr->coords[1],
                      atom_ptr->coords[2],occupancy,bvalue);
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next_atom_ptr;
        }

      /* Write a TER record out */
      fprintf(outfile_ptr,"TER   \n");
    }
}
/***********************************************************************

save_vertices  -  Write out all the vertices to the speedfill.depth file

***********************************************************************/

void save_vertices(FILE *file_ptr,struct vertex *first_vertex_ptr,
                   double *volume)
{
  char burexp;

  int atom_number, icluster, last_cluster, nvertices;

  double x, y, z;

  struct vertex *vertex_ptr;

  /* Initialise variables */
  last_cluster = -1;
  nvertices = 0;

  /* Write out separator record */
  fprintf(file_ptr,"REMARK  Vertices:\n");

  /* Get pointer to the first vertex */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all vertices */
  while (vertex_ptr != NULL)
    {
      /* Get the cluster this vertex belongs to */
      icluster = vertex_ptr->cluster;

      /* If this vertex belongs to a new cluster, write out a cluster
         record */
      if (icluster != last_cluster)
        {
          /* Write out the cluster record */
          fprintf(file_ptr,"Cluster: %4d   Volume: %8.2f\n",icluster,
                  volume[icluster]);
          
          /* Save the cluster number */
          last_cluster = icluster;
        }

      /* Get the coordinates of this vertex */
      x = vertex_ptr->coords[0];
      y = vertex_ptr->coords[1];
      z = vertex_ptr->coords[2];

      /* If current vertex has an associated atom, get the atom
         details */
      if (vertex_ptr->atom_ptr != NULL)
        {
          /* Get the corresponding atom number */
          atom_number = vertex_ptr->atom_ptr->atom_number;
        }
      else
        atom_number = 0;

      /* Determine whether vertex is buried or exposed */
      if (vertex_ptr->exposed == TRUE)
        burexp = 'e';
      else
        burexp = 'b';

      /* Write out this vertex */
      fprintf(file_ptr,"%c%5d%8.3f%8.3f%8.3f%8.2f%8.2f\n",burexp,
              atom_number,x,y,z,vertex_ptr->ratio,vertex_ptr->bvalue);

      /* Increment count of vertices written out */
      nvertices++;

      /* Get pointer to the next vertex */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char chain[MAX_CHAINS], ligand_id[6], pdb_code[5];
  char pdb_name[FILENAME_LEN], reslist_name[FILENAME_LEN];
  char dfn_name[FILENAME_LEN], pdbsum_dir[FILENAME_LEN];
  char depth_name[FILENAME_LEN], dummy_name[FILENAME_LEN];
  char log_name[FILENAME_LEN], out_name[FILENAME_LEN];
  char chains_dat[FILENAME_LEN], out_dir[FILENAME_LEN];
  char pdb_out_name[FILENAME_LEN], script_out_name[FILENAME_LEN];
  char clefts_pdb[FILENAME_LEN], clefts_pml[FILENAME_LEN];
  char clefts_rasmol[FILENAME_LEN], tmp[FILENAME_LEN];
  char rasmol_out_name[FILENAME_LEN], res_cons_name[FILENAME_LEN];
  char html_name[FILENAME_LEN], surfnet_name[FILENAME_LEN];
  char arr_name[FILENAME_LEN], r3d_name[FILENAME_LEN];
  char at_name[FILENAME_LEN], db_name[FILENAME_LEN];
  char tmp_pml[FILENAME_LEN];

  char *tmp_dir;

  int natoms, nchain, nclusters, nlig_atoms, npoints[3], nresid;
  int ndelphi_grid, nspheres, nvertices;
  int len, nmasked, nover, nunder;
  int relative_ijk[MAX_SPHERE_POINTS][3], index_ijk[MAX_SPHERE_POINTS];
  int atom_number;
  int get_cons;
  int cluster_method = OLD;
/* Use old vertix filtering method as new one doesn't seem to work
   any more(?) */
  int filter_method = NEW;
//  int filter_method = OLD;
  int nradii, nvalues;
  int all_chains, dos, dir_type, file_open, os_type, write_coords;

  int *exp_lookup_table;
  int *grid_atom_map;
  int *grid_map, *save_grid_map;

  float *delphi_potential;

  double grid_sep = 0.75;
  double grid_spacing, rascale;
  double delphi_grid_sep, delphi_valmin[3], delphi_min, delphi_max;
  double distance_ijk[MAX_SPHERE_POINTS];
  double volume[MAX_TOP_GAPS];

  struct c_file_data file_name_ptr;
  struct limits atom_limits;
  struct parameters params;

  struct atom *first_atom_ptr;
  struct atom *first_lig_atom_ptr;
  struct atom **atom_index_ptr;
  struct chain *first_chain_ptr;
  struct rasdef *first_rasdef_ptr, *last_rasdef_ptr;
  struct residue *first_residue_ptr;
  struct sphere *first_sphere_ptr;
  struct sphere *first_under_sphere_ptr;
  struct sphere *first_over_sphere_ptr, *last_over_sphere_ptr;
  struct sphere **sphere_index_ptr;
  struct vertex *first_vertex_ptr;

  KdTree *kdtree_ptr;

  FILE *arr_file_ptr;
  FILE *at_file_ptr;
  FILE *depth_file_ptr;
  FILE *dummy_file_ptr;
  FILE *html_file_ptr;
  FILE *db_file_ptr;
  FILE *logfile_ptr;
  FILE *surfnet_file_ptr;
  FILE *pmlfile_ptr;
  FILE *raspdb_file_ptr;
  FILE *rascrp_file_ptr;
  FILE *rasmol_file_ptr;
  FILE *r3d_file_ptr;

  /* Initialise pointers */
  first_atom_ptr = NULL;
  first_lig_atom_ptr = NULL;
  first_rasdef_ptr = NULL;
  last_rasdef_ptr = NULL;
  first_residue_ptr = NULL;
  first_sphere_ptr = NULL;
  first_vertex_ptr = NULL;

  /* Initialise variables */
  all_chains = FALSE;
  atom_number = 0;
  dir_type = NOT_FOUND;
  reslist_name[0] = '\0';
  arr_file_ptr = NULL;
  file_open = FALSE;
  html_file_ptr = NULL;
  at_file_ptr = NULL;
  db_file_ptr = NULL;
  rascale = 1.0;
  write_coords = FALSE;

  /* Initialise file pointers */
  rascrp_file_ptr = NULL;
  pmlfile_ptr = NULL;

  /* Initialise default file names */
  strcpy(chains_dat,"chains.dat");
  strcpy(clefts_pdb,"clefts.pdb");
  strcpy(clefts_pml,"clefts.pml");
  strcpy(clefts_rasmol,"clefts.rasmol");
  strcpy(dfn_name,"rasmol.dfn");
  strcpy(arr_name,"speedfill.arr");
  strcpy(html_name,"speedfill.html");
  strcpy(at_name,"speedfill.atoms");
  strcpy(db_name,"speedfill.dbt");
  strcpy(out_name,"speedfill.pdb");
  strcpy(pdb_out_name,"speedfill.raspdb");
  strcpy(script_out_name,"speedfill.rascrp");
  strcpy(rasmol_out_name,"speedfill.rasmol");
  strcpy(res_cons_name,"conserv.dat");
  strcpy(r3d_name,"speedfill.r3d");
  strcpy(surfnet_name,"surfnet.fil");
  strcpy(depth_name,"speedfill.depth");
  strcpy(dummy_name,"speedfill.done");
  strcpy(tmp_pml,"tmp.pml");
  log_name[0] = '\0';
  nclusters = 0;
  out_dir[0] = '\0';

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,pdb_name,out_name,&params,
                        pdb_code,chain,ligand_id,&nchain,log_name,
                        rasmol_out_name,out_dir,&logfile_ptr,&grid_sep);
  if (nchain == 1 && chain[0] == ' ')
    {
      all_chains = TRUE;
      nchain = 0;
    }

  /* Get required paths and file names */
  c_get_paths(&file_name_ptr,params.run_64);

  /* For PDBsum1, use temporary PyMOL .pml file */
  if (params.pdb_type != PDBSUM1)
    strcpy(tmp_pml,"clefts.pml");

  /* Get whether operating system is DOS */
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

  /* If have an output directory, prepend it to the names of the
     input and output files */
  if (out_dir[0] != '\0')
    prepend_dir(out_dir,arr_name,html_name,depth_name,out_name,
                pdb_out_name,r3d_name,script_out_name,surfnet_name,
                rasmol_out_name,clefts_pml,tmp_pml,dummy_name,&params);

  /* Determine whether residue conservation scores will be
     required */
  if (params.col_by_conservation == TRUE ||
      params.vertex_colouring == CONSERVATION ||
      params.generate_html == TRUE)
    get_cons = TRUE;
  else
    get_cons = FALSE;

  /* If we don't have a PDB file name, look for it in the PDB
     appropriate directories */
  if (pdb_name[0] == '\0')
    {
      /* Form the name of the PDBsum directory */
      dir_type = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                                 params.pdb_type,pdbsum_dir);
      if (dir_type == NOT_FOUND)
        {
          printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
                 pdb_code);
          exit(-1);
        }

      /* Locate the PDB file */
      dir_type = r_find_pdb_file(pdb_code,file_name_ptr.pdb_template,
                                 pdb_name,params.pdb_type);

      /* If file not found, use the file.pdb in the PDBsum directory */
      if (dir_type == NOT_FOUND)
        {
          strcpy(pdb_name,pdbsum_dir);
          strcat(pdb_name,"/file.pdb");
        }
    }
  else
    pdbsum_dir[0] = '\0';

  /* Form the required input file names */
  form_file_names(pdbsum_dir,pdb_code,chain[0],ligand_id,
                  pdb_name,reslist_name,res_cons_name,
                  dfn_name,chains_dat,out_dir,params.pdb_type,
                  params.profunc,params.use_consurf_data,
                  file_name_ptr.other_dir[NEW_CONSURF_DIR],
                  file_name_ptr.other_dir[PDBSUM_PROFUNC_DIR]);

  /* Read in the RasMol definitions */
  if (params.named_file == FALSE || params.profunc == TRUE ||
      params.pdb_type == PDBSUM1)
    get_dfn_data(dfn_name,&first_rasdef_ptr,&last_rasdef_ptr);

  /* Read through the PDB file to get the max and min coords of the
     part required */
  natoms = read_pdb(pdb_name,&first_atom_ptr,&first_lig_atom_ptr,
                    &nlig_atoms,&first_residue_ptr,&nresid,params,chain,
                    &nchain,&first_rasdef_ptr,last_rasdef_ptr,
                    logfile_ptr,all_chains,params.flip);

  /* Calculate the x-, y- and z- atom limits */
  calc_atom_limits(first_atom_ptr,&atom_limits,logfile_ptr,
                   params.max_gap_radius);

  /* Determine which ligands lie within the atom limits */
  mark_inrange_ligands(first_lig_atom_ptr,&atom_limits);

  /* If writing the Raster3D file, open file, write headers, and write
     out the protein in wireframe */
  if (params.generate_r3d == TRUE)
    {
      /* Open the output Raster3D file */
      r3d_file_ptr = open_output_file(r3d_name,TRUE);

      /* Write out the header records to the Raster3D output file */
      write_raster3d_headers(r3d_file_ptr,"speedfill","white",
                             params.raster3d_size,
                             params.raster3d_size);

      /* Calculate scaling factors for Raster3D output */
      rascale = calculate_rascale(atom_limits,BORDER);

      /* Write out the protein in wireframe */
      write_raster3d_protein(r3d_file_ptr,first_atom_ptr,
                             atom_limits.c_of_g,rascale,params);

      /* Write out any metal and ligand atoms */
      write_raster3d_ligatoms(r3d_file_ptr,first_lig_atom_ptr,
                              atom_limits.c_of_g,rascale,params);
    }

  /* If require the conservation scores, read these in */
  if (get_cons == TRUE &&
      (params.read_vertices == FALSE || params.have_rescons == FALSE))
    {
      /* Read in the chain-equivalences for correctly assigning the
         conservation scores */
      nchain = read_chain_equivs(chains_dat,&first_chain_ptr,chain,
                                 all_chains,nchain);

      /* Read in the appropriate conservation scores */
      if (params.use_consurf_data)
        get_consurf_scores(pdb_code,all_chains,chain,nchain,
                           file_name_ptr.other_dir[NEW_CONSURF_DIR],
                           first_chain_ptr,first_residue_ptr,
                           logfile_ptr);
      else
        get_conservation_scores(res_cons_name,first_residue_ptr,
                                first_chain_ptr,logfile_ptr);
    }

  /* Open the output file containing surfnet.par lines, if required */
  if (params.write_surfnet_input == TRUE)
    {
      /* Open the file and write out the file header records */
      surfnet_file_ptr = open_output_file(surfnet_name,TRUE);
      write_surfnet_headers(surfnet_file_ptr);
    }

  /* For PDBsum1, use real PDBsum directory */
  if (params.pdb_type == PDBSUM1)
    {
      /* Form the name of the PDBsum directory */
      tmp_dir
        = form_filename(pdb_code,
                        file_name_ptr.other_dir[PDBSUM_TEMPLATE]);
      strcpy(pdbsum_dir,tmp_dir);
    }

  /* RasMol output, if required */
  if (params.no_rasmol == FALSE)
    {
      /* Open output RasMol script and coordinate output files */
      open_rasmol_file(&rascrp_file_ptr,script_out_name);

      /* Write out the script header details */
      write_script_headers(rascrp_file_ptr,pdb_name,pdbsum_dir,
                           params);

      /* Set flag to write out protein and ligand, if required */
      if (params.no_protein == FALSE)
        write_coords = TRUE;
    }

  /* PyMOL output, if required */
  if (params.pymol == TRUE || params.pymol_www == TRUE)
    {
      /* Open output PyMOL script and file */
      pmlfile_ptr = open_output_file(tmp_pml,TRUE);

      /* Write header records */
      write_pymol_headers(pmlfile_ptr,pdbsum_dir,pdb_code,
                          params.relinks,clefts_pdb,file_name_ptr,
                          params);

      /* Set flag to write out protein and ligand coords */
      write_coords = TRUE;
    }

  /* Open the rasmol.pdb file */
  if (write_coords == TRUE || params.no_rasmol == FALSE)
    {
      open_rasmol_file(&raspdb_file_ptr,pdb_out_name);
      file_open = TRUE;
    }

  /* Write out the atom and ligand coordinates to the output
     PDB file  output files */
  if (write_coords == TRUE)
    {
      /* Open the coordinates file */
      if (file_open == FALSE)
        open_rasmol_file(&raspdb_file_ptr,pdb_out_name);

      /* Write out the coordinates */
      write_out_coords(raspdb_file_ptr,rascrp_file_ptr,pmlfile_ptr,
                       first_atom_ptr,first_lig_atom_ptr,
                       &atom_number,logfile_ptr,nclusters,params,
                       surfnet_file_ptr);
    }

  /* For PDBsum1, reset PDBsum directory */
  if (params.pdb_type == PDBSUM1)
    pdbsum_dir[0] = '\0';

  /* If reading in vertices from the speedfill.depth file, pick them
     up now */
  if (params.read_vertices == TRUE)
    {
      /* Read in the vertices and match to atoms */
      nvertices
        = read_vertices(depth_name,&first_vertex_ptr,first_atom_ptr,
                        natoms,&nclusters,volume);
      fprintf(logfile_ptr,"Number of vertices read in from file:");
      fprintf(logfile_ptr,"%d\n",nvertices);

      /* Write the vertices out to the RasMol script and PDB files */
      if (params.no_rasmol == FALSE)
        write_all_vertices(raspdb_file_ptr,rascrp_file_ptr,pmlfile_ptr,
                           first_vertex_ptr,nvertices,nclusters,
                           atom_number,volume,params);
    }

  /* Generate gap regions, if required */
  else if (params.no_gaps == FALSE)
    {
      /* Calculate the size of the grid-array onto which to map the sphere
         clusters */
      grid_spacing = calc_grid_size(atom_limits.valmin,atom_limits.valmax,
                                    npoints,params.max_gap_radius,grid_sep,
                                    logfile_ptr);

      /* If colouring by electrostatic potential, read in the Delphi grid */
      if (params.vertex_colouring == ELECTROSTATIC ||
          (params.save_vertices == TRUE && params.delphi_file[0] != '\0'))
        ndelphi_grid
          = read_delphi_file(params.delphi_file,&delphi_potential,
                             &delphi_grid_sep,delphi_valmin,&delphi_min,
                             &delphi_max,atom_limits,logfile_ptr);

      /* Create the grid array */
      create_grid_array(&grid_map,npoints);

      /* Create a copy of the grid array, if required, to hold grid-points
         for protein atoms and gap spheres for writing out to speedfill.arr
         file */
      if (params.save_grid == TRUE)
        create_grid_array(&save_grid_map,npoints);

      /* If need to link vertices to atoms, then create an array of atom
         pointers */
      if (params.store_atom_pointers == TRUE)
        create_atom_pointers_array(&grid_atom_map,npoints);

      /* Create the 3 index arrays that will hold atom addresses sorted
         by their x-, y- and z-coordinates */
      create_index_arrays(first_atom_ptr,&atom_index_ptr,natoms);

      /* Sort the atoms by their x-, y- and z-coordinates, storing the 
         sorted pointers in the 3 arrays of atom_index_ptr */
      sort_index_arrays(atom_index_ptr,natoms);

      /* Use all the atom radii to mask the grid-points in the grid-array
         that are within the protein */
      nmasked
        = mask_grid_points(grid_map,grid_atom_map,npoints,atom_index_ptr,
                           natoms,atom_limits.valmin,atom_limits.valmax,
                           grid_sep,logfile_ptr,params);

      /* Mark all grid-points that lie in potential gap-sphere positions
         (ie between a pair of atoms the right distance apart) */
      mark_potential_gap_points(grid_map,npoints,atom_limits.valmin,
                                atom_limits.valmax,grid_sep,atom_index_ptr,
                                natoms,first_atom_ptr,atom_limits.min_radius,
                                atom_limits.max_radius,params,relative_ijk,
                                index_ijk,distance_ijk,logfile_ptr);

      /* Place spheres at the free grid-points and determine whether they
         match the acceptance criteria */
      place_spheres(grid_map,grid_atom_map,npoints,atom_limits.valmin,
                    atom_limits.valmax,grid_sep,atom_index_ptr,natoms,
                    atom_limits.max_radius,params,relative_ijk,index_ijk,
                    distance_ijk,logfile_ptr);

      /* Generate the gap-filling spheres */
      nspheres = generate_gap_spheres(grid_map,npoints,atom_limits.valmin,
                                      atom_limits.valmax,grid_sep,
                                      &first_over_sphere_ptr,
                                      &last_over_sphere_ptr,
                                      &first_under_sphere_ptr,
                                      params.min_gap_radius,
                                      params.max_gap_radius,
                                      &nunder,&nover,logfile_ptr);
      fprintf(logfile_ptr,"Number of spheres generated:  %d\n",nspheres);

      /* If have generated any spheres, then continue */
      if (nspheres > 0 && nover > 0)
        {
          /* Use old or new clustering method, as required */
          if (cluster_method == OLD)
            {
              /* Cluster the spheres into non-interacting clusters */
              nclusters = cluster_spheresx(first_over_sphere_ptr,logfile_ptr,
                                           grid_sep);

              /* Transfer all spheres below the size threshold that
                 overlap with above-threshold spheres to the over-sized
                 spheres list */
              nspheres = transfer_spheresx(first_over_sphere_ptr,
                                           last_over_sphere_ptr,
                                           first_under_sphere_ptr,nover,
                                           nunder,logfile_ptr);
            }

          /* New clustering method employs kd-trees for quickly finding
             neighbouring spheres to cluster */
          else
            {
              /* Generate the kd-tree for the spheres */
              create_sphere_kdtree(first_over_sphere_ptr,nover,
                                   &kdtree_ptr,&sphere_index_ptr);

              /* Cluster the spheres into non-interacting clusters */
              nclusters = cluster_spheres(first_over_sphere_ptr,nover,
                                          logfile_ptr,grid_sep,kdtree_ptr,
                                          sphere_index_ptr,
                                          params.max_gap_radius);

              /* Transfer all spheres below the size threshold that
                 overlap with above-threshold spheres to the
                 over-sized spheres list */
              nspheres = transfer_spheres(first_over_sphere_ptr,
                                          last_over_sphere_ptr,
                                          first_under_sphere_ptr,nover,
                                          nunder,logfile_ptr,kdtree_ptr,
                                          sphere_index_ptr,
                                          params.max_gap_radius);

              /* Free up the memory taken for the kd-tree and the sphere
                 index array */
              KdTree_free(kdtree_ptr);
              free(sphere_index_ptr);
            }

          /* Get pointer to first sphere */
          first_sphere_ptr = first_over_sphere_ptr;

          /* Sort clusters according to their number of spheres */
          sort_clusters(first_sphere_ptr,nclusters,nover,params);

          /* Write out the spheres to a PDB file for input into SURFNET */
          nclusters = write_spheres(out_name,first_sphere_ptr,params,
                                    logfile_ptr);
        }

      /* Precalculate a table of exponentials for use in updating the
         map with spheres */
      calc_exp_table(&exp_lookup_table,&nradii,&nvalues,grid_sep,params);

      /* Open the output speedfill.arr file, if required */
      if (params.save_grid == TRUE)
        {
          if ((arr_file_ptr = fopen(arr_name,"wb")) ==  NULL)
            {
              printf("\n*** Warning. Unable to open output file [%s]\n",
                     arr_name);
              params.save_grid = FALSE;
            }
        }

      /* Open the output HTML file, if required */
      if (params.generate_html == TRUE)
        html_file_ptr = open_output_file(html_name,TRUE);

      /* Open the output database files, if required */
      if (params.generate_html == TRUE || params.pdb_type == PDBSUM1)
        {
          at_file_ptr = open_output_file(at_name,TRUE);
          db_file_ptr = open_output_file(db_name,TRUE);
        }

      /* Write out the surface of each cluster as a contour-map (in
         PDB format) for reading into RasMol */
      nvertices = contour_clusters(raspdb_file_ptr,rascrp_file_ptr,pmlfile_ptr,
                                   r3d_file_ptr,grid_map,save_grid_map,
                                   npoints,atom_limits.valmin,
                                   atom_limits.valmax,grid_sep,
                                   first_sphere_ptr,nclusters,
                                   params,&first_vertex_ptr,
                                   &atom_number,first_lig_atom_ptr,
                                   pdb_code,chain,nchain,pdb_name,
                                   file_name_ptr.pdbsum_data_dir[dir_type],
                                   logfile_ptr,html_file_ptr,
                                   db_file_ptr,at_file_ptr,surfnet_file_ptr,
                                   arr_file_ptr,grid_atom_map,atom_index_ptr,
                                   natoms,pdbsum_dir,
                                   file_name_ptr.other_dir[PDBSUM_DIR],
                                   first_atom_ptr,
                                   first_residue_ptr,exp_lookup_table,
                                   nradii,nvalues,volume,atom_limits,
                                   rascale,filter_method,delphi_potential,
                                   ndelphi_grid,delphi_grid_sep,
                                   delphi_valmin);
      fprintf(logfile_ptr,"\n");
      fprintf(logfile_ptr,"Number of vertices created = %d\n",nvertices);
      fflush(logfile_ptr);
    }

  /* Write all the final records to the output RasMol script file */
  if (params.no_rasmol == FALSE)
    write_script_tail(rascrp_file_ptr);

  /* Ditto for PyMOL script */
  if (params.pymol == TRUE || params.pymol_www == TRUE)
    {
      /* CLose off PyMOL script */
      write_pymol_tail(pmlfile_ptr,params);

      /* Close the file */
      fclose(pmlfile_ptr);
    }

  /* Close off the surfnet.par file, if we have one */
  if (params.write_surfnet_input == TRUE)
    close_surfnet_par(surfnet_file_ptr);

  /* Finish off RasMol processing, if required */
  if (params.no_rasmol == FALSE)
    {
      /* Open the output RasMol file, if required */
      rasmol_file_ptr = open_output_file(rasmol_out_name,TRUE);

      /* Close the two RasMol output files and re-open in read-mode */
      raspdb_file_ptr = close_and_reopen(raspdb_file_ptr,pdb_out_name);
      rascrp_file_ptr = close_and_reopen(rascrp_file_ptr,script_out_name);

      /* Combine the RasMol coordinates and script files into a single
         file for viewing in RasMol */
      transfer_rasmol_file(rascrp_file_ptr,rasmol_file_ptr);
      transfer_rasmol_file(raspdb_file_ptr,rasmol_file_ptr);

      /* If in PDBsum1, rename speedfill.raspdb to clefts.pdb, and
         speedfill.rascrp to clefts.rasmol */
      if (params.pdb_type == PDBSUM1 && params.relinks == FALSE)
        {
          /* Rename the files */
          copy_file(pdb_out_name,clefts_pdb);
          copy_file(script_out_name,clefts_rasmol);
        }

      /* For the relative links, use speedfill.rasmol as clefts.rasmol */
      else if (params.pdb_type == PDBSUM1 && params.relinks == TRUE)
        {
          copy_file(pdb_out_name,clefts_pdb);
          copy_file(rasmol_out_name,clefts_rasmol);

          /* For PyMOL, rewrite the tmp.pml script */
          if (params.pymol == TRUE || params.pymol_www == TRUE)
            rewrite_tmp_pml(tmp_pml,clefts_pml,params.relinks,clefts_pdb,
                            pdb_code);
        }
      
      /* Delete the coordinate and script files */
      if (params.delete_temp_files == TRUE)
        delete_temp_files(script_out_name,pdb_out_name,dos);
    }

  /* If atom depth analysis required, then write out */
  if (params.depth_analysis == TRUE || params.save_vertices == TRUE)
    {
      /* Open the output file */
      depth_file_ptr = open_output_file(depth_name,TRUE);

      /* Write out the depth analysis */
      write_atom_depths(depth_file_ptr,first_atom_ptr,first_lig_atom_ptr);

      /* Write out the vertices, if required */
      if (params.save_vertices == TRUE)
        save_vertices(depth_file_ptr,first_vertex_ptr,volume);

      /* Close the depth file */
      fclose(depth_file_ptr);
    }

  /* Open and close the dummy output file */
  if (params.save_vertices == TRUE)
    {
      dummy_file_ptr = open_output_file(dummy_name,TRUE);
      fprintf(dummy_file_ptr,"Done\n");
      fclose(dummy_file_ptr);
    }

  return(0);
}
