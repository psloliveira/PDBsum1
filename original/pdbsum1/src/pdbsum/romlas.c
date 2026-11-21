/***********************************************************************

  romlas.c - Program to generate a RasMol script (if required) and then
             output all the required ATOM, HETATM and CONECT records
             from the given PDB file for display by RasMol.

             The program is used for sending the ligands (and interacting
             protein residues, where required) to RasMol from the ligand
             and LIGPLOT pages in PDBsum. The ligand details come from
             the relevant ligplot_mm_nn.res file in the PDBsum directory.
             (File created by running LIGPLOT with the -r option when
             generating the PDBsum data).

             It is also used for fitting the selected structure from a
             SAS alignment and sending the relevant fitted coords to
             RasMol for display.

             Program is called from the RunRomLas.pl script in the CATH
             cgi directory.

************************************************************************

Date:-         16 May 1997
Amended:-      23 Feb 1998  -  Addition of RasMol script option and
                               display of multiple ligands from families
                               of structures
               14 Jul 1999  -  Tidy-up of PROSITE-pattern routines
Dir update:-    2 Oct 2012

Written by:-   Roman Laskowski

               (Fitting routines qikfit and matfit written by Andrew Martin,
                based on FORTRAN versions written by Mike Sutcliffe.)

Test runs (add &verbose=TRUE to get parameters sent to romlas):-

PROSITE patterns:-

   http://www.biochem.ucl.ac.uk/cgi-bin/cath/RunRomLas.pl?
      pdb=1argB&origin=prosite&script=YES&prosite_id=PS00105

   romlas 1argB -mc -chain -prosite PS00105 -f
   romlas 1arg -mc -allchains -prosite PS00105 -f

Ligands:-

   http://www.biochem.ucl.ac.uk/cgi-bin/cath/RunRomLas.pl?
      origin=pdbsum&pdb=1arg_01_01&residues=-l&label_atoms=YES&script=YES

   romlas 1arg_01_01 -l -la -f

   http://www.biochem.ucl.ac.uk/cgi-bin/cath/RunRomLas.pl?
      origin=pdbsum&pdb=1arg_01_01&residues=-lr&label_ligand=YES&
      label_protein=YES&script=YES

   romlas 1arg_01_01 -lr -ll -lp -f

Het groups:-

   http://www.biochem.ucl.ac.uk/cgi-bin/cath/RunRomLas.pl?
      origin=pdbsum&pdb=1auf_hetADG&residues=-l&script=YES

   romlas 1auf_hetADG -l -f

SAS:-

   http://www.biochem.ucl.ac.uk/cgi-bin/cath/RunRomLas.pl?
      aln_file=/tmp/sas13892.aln&origin=SAS&pdb=1sro&script=YES&
      residues=-bb&fit=-fbb&colour=-col

   romlas 1sro -bb -fbb -cath /tmp/sas13892.aln -col -f

Binding sites:-

   http://www.biochem.ucl.ac.uk/cgi-bin/cath/RunRomLas.pl?
   origin=bsite&family=ice&mask=1iceA0&pdb=1iceA_01_01&script=YES&
   colours=YES

   romlas 1iceA_01_01 -bs -fam ice -col -bsm 1iceA0 -f

InterProScan motifs:-

   romlas 1argB -mc -chain -rlist ".100.A:.104.A+.167.A:.174.A" -f

Colouring by exon

   romlas abf9 -chain -rlist "...1.:..19." -f -clist exons -bcart abf900

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      File name         Description
-----------      ---------         -----------
CATHPARAM            -             Input parameter file giving all the
                                   directory paths for the files used by
                                   the program
<filename>.pdb   pdb_name          Input PDB file holding coords of
                                   protein to be processed
ligplot.res      ligplot_res_name  Input file listing all the residue
                                   names and numbers to be read in from
                                   PDB file and output to standard output.
romlas.out       out_name          Output file (created when -f option
                                   given - otherwise output goes to stdout)
<sasname>.aln    cath_name         CATH-format alignment file used as the
                                   basis for fitting the structure onto
                                   one another, when called from SAS.
sel.library      prosite_datafile  File of matches to PROSITE patterns,
                                   created by Kasuya's programs.
cath.domains     cath_domains_file CATH domains file giving residue-ranges
                                   for each domain in each PDB file
<family-name><chain><_domain>.mat  Transformation relative to representative
                                   sequence of the given family
<family-name><chain><_domain>.res  Residues making up the binding site of
                                   the sequence, where the binding site is
                                   as defined relative to representative
                                   sequence of the given family
<family-name><chain><_domain>.vtx  Binding-site mask vertices
<family-name><chain><_domain>.pNN  X-SITE probe region vertices, where NN
                                   corresponds to the probe-number
promotif.dat                       Input file in PDBsum directory holding
                                   all motif info from PROMOTIF

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> create_sequence_record
         -> store_string (in common.c)
         -> get_tokens (in common.c)
         -> create_reslist_record
   -> assign_display_options
   -> c_get_paths (in common.c)
         -> get_parameter (in reapar.c)
   -> reverse_path_order
   -> find_fixed_dir
   -> open_output_file (in common.c)
   -> write_script_headers
         -> form_filename (in common.c)
   -> write_pymol_headers
         -> form_filename (in common.c)
   -> get_cora_alignment
         -> get_tokens (in common.c)
         -> create_sequence_record
         -> free_tokens (in common.c)
         -> get_aa_name
         -> create_list_record
         -> create_fitlist_record
   -> form_file_names
         -> r_find_pdb_file (in readpdb.c)
         -> find_pdbsum_dir (in common.c)
   -> get_dfn_data
         -> string_truncate (in common.c)
         -> create_rasdef_entry
         -> get_ligand_residues
               -> get_token (in common.c)
               -> format_rasdef_number
               -> create_rasdef_entry
   -> get_prosite_patterns
         -> get_prosite_strings
         -> check_for_duplicate_pattern
   -> get_residue_list
         -> create_list_record
   -> get_motif
         -> get_helices
               -> string_truncate (in common.c)
         -> string_truncate (in common.c)
         -> get_beta_turn
               -> create_range_record
         -> get_gamma_turn
               -> create_range_record
         -> get_bulge
               -> create_range_record
         -> get_hairpin
               -> create_range_record
         -> get_helix
               -> create_range_record
               -> format_resname
                     -> to_lower (in common.c)
                     -> remove_leading_blanks (in common.c)
         -> get_strand
               -> create_range_record
               -> format_resname
                     -> to_lower (in common.c)
                     -> remove_leading_blanks (in common.c)
   -> get_site
         -> r_open_pdb_file (in readpdb.c)
         -> string_truncate
         -> create_list_record
   -> get_atom_records
         -> read_in_domain_data
               -> extract_domain_info (in common.c)
                     -> get_tokens (in common.c)
                     -> free_tokens (in common.c)
               -> extract_domain_data -- No longer used
                     -> get_token (in common.c)
                     -> format_residue_number (in common.c)
         -> string_truncate
         -> get_atom_details
         -> assign_to_object
         -> check_if_want_residue
               -> get_aa_code (in common.c)
               -> check_if_in_pattern
               -> check_if_in_list
               -> get_domain
               -> create_residue_record
         -> check_ranges
         -> check_reslist
         -> check_if_want_atom
               -> check_for_hydrogen
         -> process_conect_record
   -> get_transformation
   -> calculate_fit
         -> matfit
               -> qikfit
   -> apply_transformation
   -> get_hbonds
         -> create_hbond_entry
   -> write_pymol_residues
         -> add_residue_range
   -> write_script_residues
         -> add_residue_range
   -> write_script_reslist
         -> remove_blanks
         -> get_exon_colour
         -> get_rasmol_colour
   -> write_script_sascolour
         -> form_residue_def
   -> write_script_sidechain
         -> add_residue_range
         -> get_rasmol_colour
   -> write_script_metals
         -> check_for_metal
   -> write_script_hbonds
   -> write_script_promotif
         -> colour_chains
         -> remove_leading_blanks
   -> write_script_site
         -> colour_chains
         -> string_truncate (in common.c)
         -> show_ligands
   -> write_script_csa
         -> render_objects
               -> get_chn_colour
               -> remove_blanks (in common.c)
               -> to_lower (in common.c)
   -> check_distance
   -> get_max_min_coords
   -> read_surface_vertices
         -> find_pdbsum_dir (in common.c)
   -> write_script_mask
   -> write_script_file_labels
   -> write_script_tail
   -> write_atom_records
   -> write_surface_vertices
         -> get_transformation
   -> write_dummy_atoms
   -> write_conect_records


Compile as follows:-
------------------

cc -c romlas.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -o romlas romlas.o common.o reapar.o readpdb.o -lm -lz

cc -O2 -funroll-loops -c romlas.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o romlas romlas.o common.o reapar.o readpdb.o -lm -lz


*/

#include "common.h"
#include "readpdb.h"
#include "romlas.h"

/* Prototypes
   ---------- */
 
/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int extract_domain_info(char *input_line,char *wolf_code,char *cath_code,
                        char domain_details[C_MAXDOMCHAIN][16]);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
void format_residue_number(char residue_number[6],char *token);
char *form_filename(char *pdb_code,char *filename_template);
void free_tokens(char **token,int ntokens);
char get_aa_code(char *res_name);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void remove_blanks(char *string);
void remove_leading_blanks(char *string);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);
void to_lower(char *search_string,int length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* Function prototypes from readpdb.c */
int r_find_pdb_file(char *pdb_code,char *pdb_template[NPDB_DIR],
                    char *file_name,int pdb_type);
#ifdef TEST
FILE *r_open_pdb_file(char *file_name,int *gz);
#else
gzFile r_open_pdb_file(char *file_name,int *gz);
#endif


/***********************************************************************

create_sequence_record  -  Create and initialise a new sequence record
                           for the alignment details of the current PDB
                           sequence

***********************************************************************/

void create_sequence_record(struct sequence **last_sequence_ptr,
                            int ifile,char pdb_code[5],char chain,
                            int domain_no,char ligand_id[10],
                            int entry_type)
{
  char chain_ch[2], domain_string[6], label[LABEL_LEN];

  struct sequence *sequence_ptr;

  /* Allocate memory for structure to hold sequence info */
  sequence_ptr = (struct sequence *) malloc(sizeof(struct sequence));
  if (sequence_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct");
      printf(" sequence\n");
      exit (1);
    }

  /* Initialise the elements of the structure */
  strcpy(sequence_ptr->pdb_code,pdb_code);
  sequence_ptr->chain = chain;
  sequence_ptr->domain_no = domain_no;
  strcpy(sequence_ptr->ligand_id,ligand_id);

  /* Initialise display option for this entry */
  sequence_ptr->display_option = NONE;

  /* Identify display options according to entry type */
  if (entry_type == SITE_RES)
    sequence_ptr->display_option = BINDING_SITE;
  else if (entry_type == M_C)
    sequence_ptr->display_option = MAIN_CHAIN;
  else if (entry_type == C_A)
    sequence_ptr->display_option = C_ALPHA;
  else if (entry_type == CART)
    sequence_ptr->display_option = CARTOON;
  else if (entry_type == LIG_BX)
    sequence_ptr->display_option = LIGAND_AND_BOX;
  else if (entry_type == ALL_AT)
    sequence_ptr->display_option = ALL_ATOMS;
  /* NOTE: Not yet implemented */
  else if (entry_type == DOTS)
    sequence_ptr->display_option = DOT_SURFACE;

  /* Initialise other elements of the structure */
  sequence_ptr->reslist_name[0] = '\0';
  sequence_ptr->alignment_column = -1;
  sequence_ptr->grey = TRUE;
  sequence_ptr->nconect = 0;
  strcpy(sequence_ptr->colour_name,"white");
  strcpy(sequence_ptr->colour_numbers,"[255,255,255]");

  /* Form the file-label that will identify the structure in RasMol */
  strcpy(label,pdb_code);
  if (!strcmp(pdb_code,"uplo"))
    strcpy(label,"query");
  if (chain != ' ')
    {
      chain_ch[0] = chain;
      chain_ch[1] = '\0';
      strcat(label,"(");
      strcat(label,chain_ch);
      strcat(label,")");
    }
  if (domain_no != 0)
    {
      sprintf(domain_string,"/d%d",domain_no);
      strcat(label,domain_string);
    }
  strcpy(sequence_ptr->label,label);

  /* Initialise pointers */
  sequence_ptr->first_list_ptr = NULL;
  sequence_ptr->last_list_ptr = NULL;
  sequence_ptr->first_atom_ptr = NULL;
  sequence_ptr->first_conect_ptr = NULL;
  sequence_ptr->first_hbond_ptr = NULL;
  sequence_ptr->first_residue_ptr = NULL;
  sequence_ptr->last_residue_ptr = NULL;
  sequence_ptr->first_prosite_pattern_ptr = NULL;
  sequence_ptr->fitlist_ptr = NULL;
  sequence_ptr->first_fitlist_ptr = NULL;
  sequence_ptr->last_fitlist_ptr = NULL;

  /* Save current pointer */
  *last_sequence_ptr = sequence_ptr;
 
  /* Store pointer in array of sequence pointers for ease of access */
  stored_sequence_ptr[ifile] = sequence_ptr;
}
/***********************************************************************

create_reslist_record  -  Create and initialise a new reslist record

***********************************************************************/

struct reslist *create_reslist_record(struct reslist **fst_reslist_ptr,
                                      struct reslist **lst_reslist_ptr,
                                      char *res_num1,char *res_num2,
                                      char chain,int nrange)
{
  struct reslist *reslist_ptr;
  struct reslist *first_reslist_ptr, *last_reslist_ptr;

  /* Initialise reslist pointers */
  reslist_ptr = NULL;
  first_reslist_ptr = *fst_reslist_ptr;
  last_reslist_ptr = *lst_reslist_ptr;

  /* Create reslist entry */
  reslist_ptr =
    (struct reslist*) malloc(sizeof(struct reslist));
  if (reslist_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct reslist\n");
      printf("***        Program romlas terminated with error.\n");
      exit (1);
    }

  /* Store the details */
  strcpy(reslist_ptr->res_num1,res_num1);
  strcpy(reslist_ptr->res_num2,res_num2);
  reslist_ptr->chain = chain;
  reslist_ptr->number = nrange;
  reslist_ptr->residue1_ptr = NULL;
  reslist_ptr->residue2_ptr = NULL;

  /* Initialise pointer to next entry */
  reslist_ptr->next_reslist_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_reslist_ptr == NULL)
    first_reslist_ptr = reslist_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_reslist_ptr->next_reslist_ptr = reslist_ptr;
                      
  /* Save the current residue's pointer */
  last_reslist_ptr = reslist_ptr;

  /* Return current pointers */
  *fst_reslist_ptr = first_reslist_ptr;
  *lst_reslist_ptr = last_reslist_ptr;
  return(reslist_ptr);
}
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

int get_command_arguments(char *argv[],int num,
                          char protein_family[FILENAME_LEN],
                          char cath_name[FILENAME_LEN],
                          struct parameters *params)
{
  char chain_id, ligand_id[10], pdb_code[5];
  char number_string[20], prosite_id[PROSITE_PATLEN];
  char chain, res_num1[6], res_num2[6], res_range[20];
  char tmp[FILENAME_LEN];

  char *string1_ptr, *string2_ptr, *token[MAX_TOKEN];

  int iarg, occurrence, i, ipos, itoken, len, len_id, nfiles, nmasks;
  int domain_no, entry_type, nrange, ntokens;
  int all_off, found, wolf;

  struct match_pattern *match_pattern_ptr, *last_match_pattern_ptr;
  struct reslist *last_reslist_ptr, *reslist_ptr;
  struct surface *last_surface_ptr, *surface_ptr;
  struct sequence *sequence_ptr;

  /* Initialise variables */
  entry_type = PDB;
  first_match_pattern_ptr = NULL;
  last_match_pattern_ptr = NULL;
  first_surface_ptr = NULL;
  nfiles = 0;
  nmasks = 0;
  nrange = 0;
  all_off = FALSE;
  wolf = FALSE;

  /* Initialise sequence pointer */
  last_surface_ptr = NULL;
  sequence_ptr = NULL;

  /* Initialise parameters */
  params->binding_sites = FALSE;
  params->write_to_file = FALSE;
  params->cath_alignment = FALSE;
  params->csa = FALSE;
  params->jmol = FALSE;
  params->pymol = FALSE;
  params->prosite_pattern = FALSE;
  params->promotif = FALSE;
  params->pdb_site = FALSE;
  params->site_id[0] = '\0';
  params->motif_type = 0;
  params->motif_no = 0;
  params->colour_by_pattern = FALSE;
  params->display_ligand_only = FALSE;
  params->display_ligand_and_residues = FALSE;
  params->display_backbone_only = FALSE;
  params->display_cartoon = FALSE;
  params->display_mainchain_only = FALSE;
  params->display_all_atoms = FALSE;
  params->display_binding_site_waters = FALSE;
  params->whole_protein = FALSE;
  params->add_sidechains = FALSE;
  params->all_chains = FALSE;
  params->whole = FALSE;
  params->metal = FALSE;
  params->single_chain = FALSE;
  params->single_domain = FALSE;
  params->include_hets = FALSE;
  params->fit_to_backbone = FALSE;
  params->fit_to_mainchain = FALSE;
  params->rasmol_script = TRUE;
  params->residue_list = &null;
  strcpy(params->rlist_colour,"cpk");
  params->colour_by_structure = FALSE;
  params->profunc = FALSE;
  params->sas_colour = FALSE;
  params->side_sas = FALSE;
  params->label_prot_residues = FALSE;
  params->label_ligand = FALSE;
  params->label_atoms = FALSE;
  params->nsite_surface_files = 0;
  params->work_dir[0] = '\0';
  params->xsite_spheres = FALSE;
  params->reverse_search = FALSE;
  params->upload = FALSE;
  params->use_wolf = FALSE;
  params->split_pdb = FALSE;
  params->box_cutoff = 10.0;
  params->all_entries = FALSE;
  params->fit_out = FALSE;
  params->atnum_reset = FALSE;
  for (i = 0; i < MAX_CHAINS; i++)
    params->chains[i] = '\0';
  params->nchains = 0;
  params->pdb_file = &null;
  params->run_64 = FALSE;

  /* Initialise residue list pointers */
  params->first_reslist_ptr = NULL;
  last_reslist_ptr = NULL;
  reslist_ptr = NULL;
  params->pdb_type = STANDARD_PDB;

  /* Initialise protein family name */
  cath_name[0] = '\0';
  protein_family[0] = '\0';

  /* If nothing entered on the command line, show options available */
  if (num < 2)
    {
      printf("\n");
      printf("*** ERROR. Correct usage is:\n");
      printf("\n");
      printf("    romlas [-wolf] pdb_code[_ligand_no/hetgp]\n");
      printf("           [pdb_code[_ligand_no] ...]");
      printf(" [-la] [-ll] [-lp] [-l]  [-m]\n");
      printf("           [-rlist residue_list]  [-clist colour]  ");
      printf("[-csa residue_list]\n");
      printf("           [-f] [-bs] [-NS] [-col | -colsas] [-sidesas]\n");
      printf("           [-fam fam_name] "
             "[-cath aln_name|-cathall aln_name]\n");
      printf("           [-bb | -cart | -mc | -all] [-fbb | -fmc] [-xs] [-fasta]");
      printf(" [-chain] [-allchains]\n");
      printf("           [-whole] [-domain] [-hets] [-dir work_dir]\n");
      printf("           [-prosite pattern_id/occurrence]\n");
      printf("           [-promotif motif-type motif-no\n");
      printf("           [-site site-id\n");
      printf("           [-bsm wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-bsr wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-ball wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-bmc wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-bca wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-bcart wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-bdots wolf_code1 [wolf_code2 ...]]\n");
      printf("           [-bx wolf_code1_NN [wolf_code2_NN ...]]\n");
      printf("           [-rev]  [-box cutoff]  [-v]  [-split]\n");
      printf("           [UPLOAD | MCSG]  [-pdbsum1]  [-fit_out]  "
             "[-atnum_reset]\n");
      printf("           -profunc  [-pdb filename]  [-64]\n");
      printf("\n");
      printf("where\n");
      printf("     * [-wolf]     = Indicates that Wolf codes rather than");
      printf(" PDB codes will be supplied\n");
      printf("     * [-wolf]     = Indicates that Wolf codes rather than");
      printf(" PDB codes will be supplied\n");
      printf("     * pdb_code(s) = 4-character PDB code or list");
      printf(" of codes (with\n");
      printf("                     chain appended, where relevant)\n");
      printf("     * ligand_no   = Number MM_NN corresponding to");
      printf(" ligplotMM_NN.res file holding\n");
      printf("                     list of residues to be output\n");
      printf("                     (Default is 01_01)\n");
      printf("     * hetgp       = hetXXX corresponding to chetXXX.res");
      printf(" file holding residue to be output\n");
      printf("     * [-la]       = Label all ligand atoms\n");
      printf("     * [-ll]       = Label ligand residues\n");
      printf("     * [-lp]       = Label protein residues\n");
      printf("     * [-l]        = Ligand residue(s) only\n");
      printf("     * [-lr]       = Ligand and interacting protein");
      printf(" residues only\n");
      printf("     * [-m]        = Ligand is a metal\n");
      printf("     * [-f]        = Output to be written to romlas.out");
      printf(" file\n");
      printf("     * [-dir work_dir] = Location of input PDB file and ");
      printf("rasmol.dfn file\n");
      printf("     * [-rlist residue_list] = List of residues to be ");
      printf("highlighted\n");
      printf("     * [-clist colour] = Colour for list of residues ");
      printf("defined by -rlist command\n");
      printf("     * [-csa residue_list] = List of CSA catalytic residues\n");
      printf("     * [-bb]       = Backbone representation only\n");
      printf("     * [-mc]       = Display mainchain atoms only\n");
      printf("     * [-all]      = Display all atoms\n");
      printf("     * [-bs]       = Prog called from Binding Sites\n");
      printf("     * [-fbb]      = Fit to backbone atoms\n");
      printf("     * [-fmc]      = Fit to mainchain atoms\n");
      printf("     * [-xs]       = X-SITE sphere distribs\n");
      printf("     * [-NS]       = No RasMol script to be output\n");
      printf("     * [-col]      = Each structure/pattern to be of a ");
      printf("different colour\n");
      printf("     * [-colsas]   = Colour residues according to the SAS");
      printf(" colour\n");
      printf("     * [-sidesas]  = Show whole sidechain for residues");
      printf(" coloured by SAS\n");
      printf("     * [-fam fam_name] = Transform structures using the");
      printf(" fam_name.mat file\n");
      printf("                     in its PDBsum directory, as generated");
      printf(" by genfam.scr, where\n");
      printf("                     fam_name is the name of the family (e");
      printf("g hivprot\n");
      printf("     * [-cath aln_name] = Use the CATH alignment file ");
      printf("aln_name to fit\n");
      printf("                     structures to one another on the fly\n");
      printf("     * [-cathall aln_name] = As -cath, but all PDB codes to ");
      printf("be included\n");
      printf("     * [-prosite pattern_id/occurrence] = Show structure ");
      printf("with residues\n");
      printf("                     matching given PROSITE pattern ");
      printf("highlighted.\n");
      printf("                     The occurrence indicates which hit of ");
      printf("the pattern is required.\n");
      printf("     * [-promotif type number] = Show structure with residues");
      printf("belonging to given\n");
      printf("                     motif highlighted\n");
      printf("     * [-site site-id] = Show SITE residues defined in PDB");
      printf(" file for given site\n");
      printf("     * [-fasta]    = Alignment is from FASTA (via SAS)\n");
      printf("     * [-chain]    = Display only the chain given for each ");
      printf("PDB entry\n");
      printf("     * [-allchains]= Display whole PDB entry\n");
      printf("     * [-whole]    = REad in whole chain from PDB entry, not "
             "just the aligned residues\n");
      printf("     * [-domain]   = Display only the domain given for each ");
      printf("PDB entry\n");
      printf("     * [-hets]     = Include nearby HET groups in SAS ");
      printf("alignment\n");
      printf("     * [-bsm wolf_code(s)]   = Binding-site mask for given ");
      printf("Wolf code(s)\n");
      printf("     * [-bsr wolf_code(s)]   = Binding-site residues for ");
      printf("given Wolf code(s)\n");
      printf("     * [-ball wolf_code(s)]  = All atoms for given Wolf ");
      printf("code(s)\n");
      printf("     * [-bmc wolf_code(s)]   = Main-chain atoms only for ");
      printf("given Wolf code(s)\n");
      printf("     * [-bca wolf_code(s)]   = C-alpha atoms only for ");
      printf("given Wolf code(s)\n");
      printf("     * [-bcart]    = Protein cartoon representation\n");
      printf("     * [-bdots wolf_code(s)] = Dot-surface representation ");
      printf("for given Wolf code(s)\n");
      printf("     * [-bx wolf_code_NN(s)] = X-SITE distribution for ");
      printf("probe-type NN\n");
      printf("     * [-rev]      = Paths to files to be searched in reverse");
      printf(" order - not used\n");
      printf("     * [-box cutoff] = Cutoff distance for extracting atoms");
      printf(" around ligand\n");
      printf("     * [-fit_out]  = Only residues used in superposition to "
             "be written out\n");
      printf("     * [-atnum_reset] = Reset atom-numbers for each structure "
             "written out\n");
      printf("     * [-v]        = Verbose mode\n");
      printf("     * [-split]    = Write out each PDB file as a separate "
             "output file\n");
      printf("     * UPLOAD      = Uploaded PDB file for PDBsum generation\n");
      printf("     * MCSG        = MCSG structure\n");
      printf("     * -pdbsum1    = Runnin in PDBsum1\n");
      printf("     * -profunc    = PDB file is in ProFunc directory\n");
      printf("     * -jmol       = Output is for JSmol\n");
      printf("     * -pymol      = Output PyMOL commands (only for PDBsum1)\n");
      printf("     * -pdb file   = PDB file to be loaded by script\n");
      printf("     * -64         = Run on 64-bit machine\n");
      printf("\n");
      printf("\n");
      printf("Examples:-\n");
      printf("\n");
      printf("   romlas  1eca_01_01  -l\n");
      printf("\n");
      printf("   romlas  7hvp  1cpi  4hvp_02_01  -l  -fam hivprot\n");
      printf("\n");
      printf("   romlas  7hvpA  1cpiA  4hvpA  -l  -cath");
      printf(" alnfile.003.040.050.010\n");
      printf("\n");
      exit(1);
    }

  /* Loop through the command arguments */
  for (iarg = 1; iarg < num + 1; iarg++)
    {
      /* Get the length of this token */
      len = strlen(argv[iarg]);

      /* Check for the -fam option giving a protein family name */
      if (!strncmp(argv[iarg],"-fam",4))
        {
          /* Get the family name from the next token, if there is one */
          if (iarg < num)
            {
              strcpy(protein_family,argv[iarg + 1]);
              iarg++;
            }          
        }

      /* Check for the -cath option giving a CATH alignment file name */
      else if (!strcmp(argv[iarg],"-cath"))
        {
          /* Get the file name for the alignment from the next token,
             if there is one */
          if (iarg < num)
            {
              strcpy(cath_name,argv[iarg + 1]);
              iarg++;
              params->cath_alignment = TRUE;
            }          
        }

      /* Check for the -cathall option giving a CATH alignment file name */
      else if (!strcmp(argv[iarg],"-cathall"))
        {
          /* Get the file name for the alignment from the next token,
             if there is one */
          if (iarg < num)
            {
              strcpy(cath_name,argv[iarg + 1]);
              iarg++;
              params->cath_alignment = TRUE;

              /* Set flag that all PDB entries required */
              params->all_entries = TRUE;
            }
        }

      /* Check for the -prosite option for showing residues matching
         given PROSITE pattern */
      else if (!strncmp(argv[iarg],"-prosite",8))
        {
          /* Get the file name for the alignment from the next token,
             if there is one */
          if (iarg < num)
            {
              /* Get the PROSITE pattern-id from the next token */
              strcpy(prosite_id,argv[iarg + 1]);
              iarg++;

              /* Separate out the pattern-id and occurrence number */
              len_id = len;
              for (ipos = 0; ipos < len; ipos++)
                if (prosite_id[ipos] == '/')
                  len_id = ipos;

              /* If have an occurrence number appended to the PROSITE
                 pattern, then extract */
              occurrence = 0;
              if (len > len_id)
                {
                  strcpy(number_string,prosite_id+len_id+1);
                  occurrence = atoi(number_string);
                }

              /* Allocate memory for this PROSITE pattern */
              match_pattern_ptr 
                = (struct match_pattern *) 
                malloc(sizeof(struct match_pattern));
              if (match_pattern_ptr == NULL)
                {
                  printf("*** Can't allocate memory for struct ");
                  printf("match_pattern\n");
                  exit (1);
                }

              /* Add link from previous pattern record to the current one */
              if (last_match_pattern_ptr != NULL)
                last_match_pattern_ptr->next_match_pattern_ptr
                  = match_pattern_ptr;
              else
                first_match_pattern_ptr = match_pattern_ptr;
              last_match_pattern_ptr = match_pattern_ptr;

              /* Store the pattern-id and occurrence number */
              strncpy(match_pattern_ptr->prosite_id,prosite_id,len_id);
              match_pattern_ptr->prosite_id[len_id] = '\0';
              match_pattern_ptr->occurrence = occurrence;
              match_pattern_ptr->chain = ' ';

              /* Initialise pointer to next pattern */
              match_pattern_ptr->next_match_pattern_ptr = NULL;
            }          

          /* Set flag to indicate that showing PROSITE patterns */
          params->prosite_pattern = TRUE;
        }

      /* Check for the -promotif option identifying which secondary 
         structure motif to be displayed */
      else if (!strcmp(argv[iarg],"-promotif"))
        {
          /* Get the motif type from the next token, if there is one */
          if (iarg < num)
            {
              strcpy(number_string,argv[iarg + 1]);
              params->motif_type = atoi(number_string);
              iarg++;
            }          

          /* Get the motif number from the next token, if there is one */
          if (iarg < num)
            {
              strcpy(number_string,argv[iarg + 1]);
              params->motif_no = atoi(number_string);
              iarg++;
              if (params->motif_type > -1 && params->motif_type < NMOTIFS)
                params->promotif = TRUE;
            }          
        }

      /* Check for the -site option giving which SITE residues to be
         displayed */
      else if (!strcmp(argv[iarg],"-site"))
        {
          /* Get the site-id from the next token, if there is one */
          if (iarg < num)
            {
              len = strlen(argv[iarg + 1]);
              if (len < 5)
                {
                  strcpy(params->site_id,argv[iarg + 1]);
                  params->pdb_site = TRUE;
                }
              iarg++;
            }          
        }

      /* Check for the -rlist option giving list of residues to be
         highlighted */
      else if (!strcmp(argv[iarg],"-rlist") ||
               !strcmp(argv[iarg],"-csa"))
        {
          /* Get the residue list from the next token */
          if (iarg < num)
            {
              /* If this is a CSA residue list, set flag and set colour
                 of highlighted residues to CPK */
              if (!strcmp(argv[iarg],"-csa"))
                {
                  params->csa = TRUE;
                  strcpy(params->rlist_colour,"cpk");
                }

              /* Store the residue list */
              store_string(&(params->residue_list),argv[iarg + 1]);

              /* Get the length of the residue list */
              len = strlen(params->residue_list);

              /* Replace all the dots by spaces */
              for (ipos = 0; ipos < len; ipos++)
                if (params->residue_list[ipos] == '.')
                  params->residue_list[ipos] = ' ';

              /* Initialise flag */
              string1_ptr = params->residue_list;

              /* Extract the tokens on this line */
              ntokens = get_tokens(string1_ptr,token,MAX_TOKEN,'+');

              /* Loop over the tokens to extract the individual residue
                 ranges */
              for (itoken = 0; itoken < ntokens; itoken++)
                {
                  /* Get length of the current token */
                  len = strlen(token[itoken]);
                  res_num1[0] = '\0';
                  res_num2[0] = '\0';

                  /* Check whether this is a single residue */
                  if (len == 5 || len == 6)
                    {
                      /* If no chain defined, set to space */
                      if (len == 5)
                        chain = ' ';
                      else
                        chain = token[itoken][5];

                      /* Get the residue number */
                      strncpy(res_num1,token[itoken],5);
                      res_num1[5] = '\0';

                      /* Save as both residues */
                      strcpy(res_num2,res_num1);
                    }

                  /* Otherwise, have two residues defining the start
                     and end of a residue range */
                  else if (len < 14)
                    {
                      /* Extract the start- and end-residues from this
                         range */
                      strncpy(res_num1,token[itoken],5);
                      res_num1[5] = '\0';
                      chain = ' ';
                      ipos = 5;
                      if (token[itoken][ipos] != ':')
                        {
                          chain = token[itoken][ipos];
                          ipos++;
                        }
                      strncpy(res_num2,token[itoken]+ipos+1,5);
                      res_num2[5] = '\0';
                    }

                  /* If have a valid residue range, store as a reslist
                     record */
                  if (res_num1[0] != '\0' && res_num2[0] != '\0')
                    {
                      /* Create a reslist record to store this residue
                         range */
                      reslist_ptr
                        = create_reslist_record(&(params->first_reslist_ptr),
                                                &last_reslist_ptr,
                                                res_num1,res_num2,chain,
                                                nrange);

                      /* If this is the first chain, store as wanted */
                      if (params->chains[0] == '\0')
                        {
                          /* Store chain */
                          params->chains[0] = chain;
                          params->nchains++;
                          params->chains[params->nchains] = '\0';
                        }

                      /* Otherwise, see if we already have this chain
                         in our list */
                      else
                        {
                          /* Initialise flag */
                          found = FALSE;

                          /* Loop over the chains stored so far */
                          for (i = 0; i < params->nchains; i++)
                            {
                              /* If have this one, set flag */
                              if (params->chains[i] == chain)
                                found = TRUE;
                            }

                          /* If chain is a new one, then store */
                          if (found == FALSE &&
                              params->nchains < MAX_CHAINS)
                            {
                              /* Store chain */
                              params->chains[params->nchains] = chain;
                              params->nchains++;
                              params->chains[params->nchains] = '\0';
                            }
                        }

                      /* Increment count of residue ranges created */
                      nrange++;
                    }
                }

              /* Increment token count */
              iarg++;

              /* Free up the memory taken by the tokens */
              free_tokens(token,ntokens);
            }          
        }

      /* Check for the -clist option giving the colour for the list of
         residues defined by the -rlist command */
      else if (!strcmp(argv[iarg],"-clist"))
        {
          /* Get the colour from the next token */
          if (iarg < num)
            {
              /* Store the colour */
              strcpy(params->rlist_colour,argv[iarg + 1]);
              iarg++;
            }          
        }

      /* Check for the -pdb option giving the name of the PDB file
         to be loaded by the script */
      else if (!strcmp(argv[iarg],"-pdb"))
        {
          /* Get the file name from the next token */
          if (iarg < num)
            {
              /* Store the colour */
              store_string(&(params->pdb_file),argv[iarg + 1]);
              iarg++;
            }          
        }

      /* Check for the -box option requiring atoms within cutoff distance
         of ligand */
      else if (!strncmp(argv[iarg],"-box",4))
        {
          /* Get the cutoff distance */
          if (iarg < num)
            {
              strcpy(number_string,argv[iarg + 1]);
              params->box_cutoff = atof(number_string);
              iarg++;
            }          
        }

      /* Check for the -fit_out option indicating that only residues
         used in the 3D superposition are to be written out */
      else if (!strcmp(argv[iarg],"-fit_out"))
        {
          /* Set flag */
          params->fit_out = TRUE;
        }

      /* Check for the -atnum_reset option indicating that atom
         numbering to be reset for each structure written out (and
         CONECT records to be written out in the same block) */
      else if (!strcmp(argv[iarg],"-atnum_reset"))
        {
          /* Set flag */
          params->atnum_reset = TRUE;
        }

      /* Check for the -ball option indicating that all protein atoms
         to be displayed */
      else if (!strncmp(argv[iarg],"-ball",5))
        {
          /* Set flag to interpret all subsequent PDB codes for all atoms */
          entry_type = ALL_AT;
        }

      /* Check for the -bmc option indicating that only main-chain atoms
         to be displayed */
      else if (!strncmp(argv[iarg],"-bmc",4))
        {
          /* Set flag to interpret all subsequent PDB codes for m/c only */
          entry_type = M_C;
        }

      /* Check for the -bcart option indicating that protein cartoon
         representation to be displayed */
      else if (!strncmp(argv[iarg],"-bcart",6))
        {
          /* Set flag to interpret all subsequent PDB codes for cartoon */
          entry_type = CART;
        }

      /* Check for the -bsm option indicating that dot-surface
         representation to be displayed */
      else if (!strncmp(argv[iarg],"-bdots",6))
        {
          /* Set flag to interpret all subsequent PDB codes for dot-surface */
          entry_type = DOTS;
        }

      /* Check for the -bca option indicating that only C-alpha atoms
         to be displayed */
      else if (!strncmp(argv[iarg],"-bca",4))
        {
          /* Set flag to interpret all subsequent PDB codes for CA only */
          entry_type = C_A;
        }

      /* Check for the -bsm option indicating that binding-site mask
         to be displayed */
      else if (!strncmp(argv[iarg],"-bsm",4))
        {
          /* Set flag to interpret all subsequent PDB codes as mask-ids */
          entry_type = MASK;
        }

      /* Check for the -bx option indicating X-SITE probe distribution
         to be displayed */
      else if (!strncmp(argv[iarg],"-bx",3))
        {
          /* Set flag to interpret all subsequent PDB codes as X-SITE
             distributions */
          entry_type = XSITE;
        }

      /* Check for the -bsr option indicating that binding-site residue
         sets to be displayed */
      else if (!strncmp(argv[iarg],"-bsr",4))
        {
          /* Set flag to interpret all subsequent PDB codes as
             binding-site residue sets */
          entry_type = SITE_RES;
          all_off = TRUE;
        }

      /* Check for the -dir option giving the name of the working
         directory */
      else if (!strncmp(argv[iarg],"-dir",4))
        {
          /* Get the directory name */
          if (iarg < num)
            {
              strcpy(params->work_dir,argv[iarg + 1]);
              iarg++;
            }          
        }

      /* Check for the -fasta flag indicating alignment came from a FASTA
         run in SAS (rather than from a CORA alignment) */
      else if (!strncmp(argv[iarg],"-fasta",6))
        params->from_fasta = TRUE;

      /* Check for the -allchains option indicating that whole structure
         to be displayed */
      else if (!strncmp(argv[iarg],"-allchains",10))
        params->all_chains = TRUE;

      /* Check for the -whole option indicating that whole chain to be
         read in from the PDB file, not just the aligned residues */
      else if (!strcmp(argv[iarg],"-whole"))
        params->whole = TRUE;

      /* Check for the -chain option for displaying only a single chain */
      else if (!strncmp(argv[iarg],"-chain",6))
        params->single_chain = TRUE;

      /* Check for the -domain option for displaying only a single domain */
      else if (!strncmp(argv[iarg],"-domain",7))
        params->single_domain = TRUE;

      /* Check for the -hets option for displaying HET groups in SAS */
      else if (!strncmp(argv[iarg],"-hets",5))
        params->include_hets = TRUE;

      /* Check for the -fbb option, indicating that structures to be
         fitted to backbone atoms */
      else if (!strncmp(argv[iarg],"-fbb",4))
        params->fit_to_backbone = TRUE;

      /* Check for the -fmc option, indicating that structures to be
         fitted to mainchain atoms */
      else if (!strncmp(argv[iarg],"-fmc",4))
        params->fit_to_mainchain = TRUE;

      /* Check for the -f option */
      else if (!strncmp(argv[iarg],"-f",2))
        params->write_to_file = TRUE;

      /* Check for the -bs option, indicating that from the binding-site
         Web pages */
      else if (!strncmp(argv[iarg],"-bs",3))
        params->binding_sites = TRUE;

      /* Check for the -bb option, indicating that only backbone to be
         displayed */
      else if (!strncmp(argv[iarg],"-bb",3))
        params->display_backbone_only = TRUE;

      /* Check for the -mc option, indicating that only mainchain atoms
         to be displayed */
      else if (!strncmp(argv[iarg],"-mc",3))
        params->display_mainchain_only = TRUE;

      /* Check for the -cart option, indicating that cartoon representartion
         to be displayed */
      else if (!strncmp(argv[iarg],"-cart",5))
        params->display_cartoon = TRUE;
      
      /* Check for the -all option, indicating that all atoms to be
         displayed */
      else if (!strncmp(argv[iarg],"-all",4))
        params->display_all_atoms = TRUE;

      /* Check for the -xs option, indicating that sphere-approximations
         to be shown for X-SITE distribs */
      else if (!strncmp(argv[iarg],"-xs",3))
        params->xsite_spheres = TRUE;

      /* Check for the -rev option, indicating that paths are to be
         searched in reverse order */
      else if (!strncmp(argv[iarg],"-rev",4))
        params->reverse_search = TRUE;

      /* Check for the -NS option, indicating that only coords to be
         output (ie no RasMol script) */
      else if (!strncmp(argv[iarg],"-NS",3))
        params->rasmol_script = FALSE;

      /* Check for the -colsas option, indicating that residues are to be
         coloured according to the SAS colouring scheme */
      else if (!strncmp(argv[iarg],"-colsas",7))
        params->sas_colour = TRUE;

      /* Check for the -sidesas option, indicating that sideshains to be
         shown of residues coloured according to SAS colouring scheme */
      else if (!strncmp(argv[iarg],"-sidesas",8))
        params->side_sas = TRUE;

      /* Check for the -col option, indicating that each structure to be
         of a different colour */
      else if (!strncmp(argv[iarg],"-col",4))
        params->colour_by_structure = TRUE;

      /* Check for -la option: label all atoms */
      else if (!strncmp(argv[iarg],"-la",3))
        params->label_atoms = TRUE;

      /* Check for -ll option: label ligand residues */
      else if (!strncmp(argv[iarg],"-ll",3))
        params->label_ligand = TRUE;

      /* Check for -lp option: label ligand residues */
      else if (!strncmp(argv[iarg],"-lp",3))
        params->label_prot_residues = TRUE;

      /* Check for -lr option: display ligand residues and the protein 
         residues they interact with */
      else if (!strncmp(argv[iarg],"-lr",3))
        params->display_ligand_and_residues = TRUE;

      /* Check for -m option: ligand is a metal ion */
      else if (!strcmp(argv[iarg],"-m"))
        params->metal = TRUE;

      /* Check for -v option: verbose mode */
      else if (!strcmp(argv[iarg],"-v"))
        verbose = TRUE;

      /* Check for -split option: output PDB files separately */
      else if (!strcmp(argv[iarg],"-split"))
        params->split_pdb = TRUE;

      /* Check for -wolf flag indicating that following codes will be Wolf
         codes rather than PDB codes */
      else if (!strncmp(argv[iarg],"-wolf",5))
        {
          wolf = TRUE;
          params->use_wolf = TRUE;
        }

      /* Check for ligand-only option */
      else if (!strncmp(argv[iarg],"-l",2))
        params->display_ligand_only = TRUE;

      /* Check for the UPLOAD option */
      else if (!strcmp(argv[iarg],"UPLOAD"))
        params->pdb_type = UPLOAD_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(argv[iarg],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the PDBsum1 option */
      else if (!strcmp(argv[iarg],"-pdbsum1"))
        params->pdb_type = PDBSUM1;

      /* Check for the JSmol option */
      else if (!strcmp(argv[iarg],"-jmol"))
        params->jmol = TRUE;

      /* Check for the PyMOL option */
      else if (!strcmp(argv[iarg],"-pymol"))
        params->pymol = TRUE;

      /* Check for the -profunc option */
      else if (!strcmp(argv[iarg],"-profunc"))
        params->profunc = TRUE;

      /* Check for the 64-bit option */
      else if (!strcmp(argv[iarg],"-64"))
        params->run_64 = TRUE;

      /* Skip the -relinks option */
      else if (!strcmp(argv[iarg],"-relinks"))
	{
	}

      /* If long enough to be a PDB code, then store */
      else if (len > 3)
        {
          /* If this is a PDB code giving a mask-id, then store */
          if (entry_type == MASK || entry_type == XSITE)
            {
              /* Allocate memory for structure to hold mask info */
              surface_ptr = (struct surface *) malloc(sizeof(struct surface));
              if (surface_ptr == NULL)
                {
                  printf("*** Can't allocate memory for struct");
                  printf(" mask\n");
                  exit (1);
                }

              /* If this is an X-SITE surface, then extract the probe
                 number */
              surface_ptr->xsite_probe[0] = '\0';
              if (entry_type == XSITE && len > 8)
                {
                  /* Get the probe number */
                  surface_ptr->xsite_probe[0] = argv[iarg][7];
                  surface_ptr->xsite_probe[1] = argv[iarg][8];
                  surface_ptr->xsite_probe[2] = '\0';
                }

              /* Update the sequence links to this mask */
              if (first_surface_ptr == NULL)
                first_surface_ptr = surface_ptr;

              /* Add link from previous mask to the current one */
              if (last_surface_ptr != NULL)
                last_surface_ptr->next_surface_ptr = surface_ptr;
              last_surface_ptr = surface_ptr;

              /* Initialise the elements of the structure */
              strcpy(surface_ptr->site_surface_code,argv[iarg]);
              strcpy(surface_ptr->colour_name,"white");
              strcpy(surface_ptr->colour_numbers,"[255,255,255]");
              surface_ptr->first_vertex_ptr = NULL;
              surface_ptr->next_surface_ptr = NULL;
              nmasks++;

              /* Initialise grid-spacing to default value */
              surface_ptr->grid_spacing = 1.6;

              /* Store the mask number */
              surface_ptr->surface_number = nmasks;
            }

          /* Otherwise, take it to be the code of a PDB file for its
             structure */
          else
            {
              /* Check that max. number of codes not exceeded */
              if (nfiles < MAXPDB - 1)
                {
                  /* Extract the PDB code */
                  strncpy(pdb_code,argv[iarg],4);
                  pdb_code[4] = '\0';
                  chain_id = ' ';
                  domain_no = 0;

                  /* Check if there is a chain-id straight after it */
                  ipos = 5;
                  if (len > 4)
                    {
                      if (argv[iarg][4] != '_')
                        {
                          /* Ligand identifier starts one position beyond
                             chain-id */
                          ipos = ipos + 1;
                          chain_id = argv[iarg][4];
                          if (chain_id == '0')
                            chain_id = ' ';

                          /* Check if there's a domain number attached */
                          if (len > 5)
                            if (argv[iarg][5] != '_')
                              {
                                /* Ligand identifier starts one position beyond
                                   domain number */
                                ipos = ipos + 1;
                                number_string[0] = argv[iarg][5];
                                number_string[1] = '\0';

                                /* Check for 2-digit domain identifier */
                                if (len == 7)
                                  {
                                    number_string[1] = argv[iarg][6];
                                    number_string[2] = '\0';
                                  }

                                /* Get the domain number */
                                domain_no = atoi(number_string);
                              }
                        }
                    }

                  /* If input string contains the ligand or Het Group
                     identifier, then store that */
                  if (len > 9)
                    {
                      /* Extract the ligand identifier */
                      strncpy(ligand_id,argv[iarg]+ipos,5);
                      ligand_id[5] = '\0';

                      /* Check whether this is a Het group identifier */
                      if (!strncmp(ligand_id,"het",3))
                        {
                          /* Copy het group name across */
                          strncpy(ligand_id,argv[iarg]+ipos,6);
                          ligand_id[6] = '\0';
                        }

                      /* Set flag to indicate we have a LIGPLOT id */
                      entry_type = LIGPLOT;

                      /* If have a "b" on the end, then want all protein
                         atoms within a given box around the ligand */
                      len = strlen(argv[iarg]+ipos);
                      if (len > 5 && argv[iarg][ipos+5] == 'b')
                        entry_type = LIG_BX;
                    }

                  /* Otherwise use default identifier (ie first instance
                     of first ligand) */
                  else if (wolf == FALSE)
                    strcpy(ligand_id,"01_01");

                  /* If this is a bind-site residue identifier, then
                     blank out ligand-id */
                  if (entry_type == SITE_RES)
                    ligand_id[0] = '\0';

                  /* Create a sequence record corresponding to this PDB code */
                  create_sequence_record(&sequence_ptr,nfiles,pdb_code,
                                         chain_id,domain_no,ligand_id,
                                         entry_type);

                  /* Increment count of PDB files */
                  nfiles++;
                }
            }
        }
    }

  /* Check that sensible/consistent parameters chosen */

  /* For CATH alignment file for fitting, from SAS */
  if (params->cath_alignment == TRUE)
    {
      /* If no display option, then set default */
      if (params->display_ligand_only == FALSE &&
          params->display_ligand_and_residues == FALSE &&
          params->display_backbone_only == FALSE &&
          params->display_cartoon == FALSE &&
          params->display_mainchain_only == FALSE &&
          params->display_all_atoms == FALSE)
        params->display_backbone_only = TRUE;

      /* If no fitting option, then set default */
      if (params->fit_to_backbone == FALSE &&
          params->fit_to_mainchain == FALSE)
        params->fit_to_backbone = TRUE;
    }

  /* Set default for PROSITE patterns */
  if (params->prosite_pattern == TRUE &&
      params->display_backbone_only == FALSE &&
      params->display_cartoon == FALSE &&
      params->display_mainchain_only == FALSE &&
      params->display_all_atoms == FALSE)
    params->display_mainchain_only = TRUE;

  /* For PROSITE patterns, if protein backbone or mainchain to be
     displayed, then set appropriate flag */
  if ((params->prosite_pattern == TRUE || params->side_sas == TRUE ||
       params->residue_list[0] != '\0') &&
      (params->display_backbone_only == TRUE ||
       params->display_cartoon == TRUE ||
       params->display_mainchain_only == TRUE ||
       params->display_all_atoms == TRUE))
    {
      /* Set flag to display whole protein */
      if (params->single_chain == FALSE && params->single_domain == FALSE)
        params->whole_protein = TRUE;

      /* Set flag to show key sidechains */
      params->add_sidechains = TRUE;
    }

  /* For binding sites database, set appropriate flag to show single
     chain if mainchain or CA-only, etc option selected */
  if (params->binding_sites == TRUE)
    params->whole_protein = TRUE;

  /* If showing binding-site residues, then switch all-chains option off */
  if (all_off == TRUE)
    params->all_chains = FALSE;

  /* Save the number of binding-site masks */
  params->nsite_surface_files = nmasks;

  /* If have a working directory name, prefix any CATH name */
  if (params->work_dir[0] != '\0')
    {
      strcpy(tmp,params->work_dir);
      strcat(tmp,DIR_SEP);
      strcat(tmp,cath_name);
      strcpy(cath_name,tmp);
    }

  /* STOP if PyMOL option selected as not yet implemented */
  if (params->pymol == TRUE)
    {
      printf("*** ERROR. PyMOL option net yet implemented!\n");
      exit(1);
    }

  /* Return the number of PDB codes entered */
  return(nfiles);
}
/***********************************************************************

assign_display_options  -  Assign each sequence's display options
                           according to the overall parameters

***********************************************************************/

void assign_display_options(int nfiles,struct parameters params,
                            int *nligands)
{
  int ifile;

  struct sequence *sequence_ptr;

  /* Loop over all the files */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* If no display option defined for this sequence, then use
         global parameters to decide how it is to be displayed */
      if (sequence_ptr->display_option == NONE)
        {
          /* Apply overall options */
          if (params.display_all_atoms == TRUE)
            sequence_ptr->display_option = ALL_ATOMS;
          else if (params.display_mainchain_only == TRUE)
            sequence_ptr->display_option = MAIN_CHAIN;
          else if (params.display_backbone_only == TRUE)
            sequence_ptr->display_option = C_ALPHA;
          else if (params.display_cartoon == TRUE)
            sequence_ptr->display_option = CARTOON;
          else if (params.display_ligand_and_residues == TRUE)
            sequence_ptr->display_option = LIGAND_AND_RES;
          else if (params.display_ligand_only == TRUE)
            sequence_ptr->display_option = LIGAND_ONLY;

          /* If display option is still unassigned, and we are coming
             from the binding site database, then this must be a ligand */
          if (sequence_ptr->display_option == NONE &&
              params.binding_sites == TRUE)
            sequence_ptr->display_option = LIGAND_ONLY;
        }

      /* Increment count of ligands according to display option */
      if (sequence_ptr->display_option == LIGAND_AND_RES ||
          sequence_ptr->display_option == LIGAND_ONLY)
        (*nligands)++;
    }
}
/***********************************************************************

create_list_record  -  Create and initialise a new list record

***********************************************************************/

void create_list_record(struct list **new_list_ptr,
                        struct sequence *sequence_ptr,char res_name[4],
                        char res_num[6],char chain)
{
  struct list *list_ptr;

  /* Allocate memory for structure to hold list info */
  list_ptr = (struct list *) malloc(sizeof(struct list));
  if (list_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct list\n");
      exit (1);
    }

  /* Update the sequence links to this list */
  if (sequence_ptr->first_list_ptr == NULL)
    sequence_ptr->first_list_ptr = list_ptr;

  /* Add link from previous list to the current one */
  if (sequence_ptr->last_list_ptr != NULL)
    sequence_ptr->last_list_ptr->next_list_ptr = list_ptr;
  sequence_ptr->last_list_ptr = list_ptr;

  /* Initialise the elements of the structure */
  list_ptr->blank = FALSE;
  list_ptr->aa_code = ' ';
  strcpy(list_ptr->res_name,res_name);
  strcpy(list_ptr->res_num,res_num);
  list_ptr->chain = chain;
  list_ptr->is_ligand = FALSE; 
  list_ptr->interacting_residue = FALSE; 
  list_ptr->fit = FALSE; 

  /* Initialise the pointers */
  list_ptr->reference_list_ptr = NULL;
  list_ptr->residue_ptr = NULL;
  list_ptr->next_list_ptr = NULL;

  /* Return current pointer */
  *new_list_ptr = list_ptr;
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

create_fitlist_record  -  Create and initialise a new fitlist record

***********************************************************************/

struct fitlist *create_fitlist_record(struct fitlist **fst_fitlist_ptr,
                                      struct fitlist **lst_fitlist_ptr)
{
  struct fitlist *fitlist_ptr;
  struct fitlist *first_fitlist_ptr, *last_fitlist_ptr;

  /* Initialise fitlist pointers */
  fitlist_ptr = NULL;
  first_fitlist_ptr = *fst_fitlist_ptr;
  last_fitlist_ptr = *lst_fitlist_ptr;

  /* Create fitlist entry */
  fitlist_ptr =
    (struct fitlist*) malloc(sizeof(struct fitlist));
  if (fitlist_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct fitlist\n");
      printf("***        Program romlas terminated with error.\n");
      exit (1);
    }

  /* Initialise details */
  fitlist_ptr->range[0] = '\0';
  fitlist_ptr->done = FALSE;

  /* Initialise pointer to next entry */
  fitlist_ptr->next_fitlist_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_fitlist_ptr == NULL)
    first_fitlist_ptr = fitlist_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_fitlist_ptr->next_fitlist_ptr = fitlist_ptr;
                      
  /* Save the current residue's pointer */
  last_fitlist_ptr = fitlist_ptr;

  /* Return current pointers */
  *fst_fitlist_ptr = first_fitlist_ptr;
  *lst_fitlist_ptr = last_fitlist_ptr;
  return(fitlist_ptr);
}
/***********************************************************************

get_cora_alignment  -  Read in the CORA-format alignment file

***********************************************************************/

int get_cora_alignment(char file_name[FILENAME_LEN],int *nseq,
                       int *nresid,int *nf,struct parameters params)
{
  char ligand_id[10], input_line[LINELEN + 1], wolf_code[8];
  char number_string[10], pdb_code[5], residue_details[12];
  char aa_code, chain_id, res_name[4], res_num[6], res_number[2][6];

  char *token[MAX_TOKEN];

  int column, domain_no, i, ifile, ipos, irange, iresid, iseq, line;
  int ntokens;
  int len, nchains, nfiles, nmatch, nresidues, nseqs;
  int column_width, start_pos, version;
  int a, aa_num, entry_type, z;
  int done, fit, fit_flagged, found, gap_ref;

  struct fitlist *fitlist_ptr;
  struct list *list_ptr, *reference_list_ptr;
  struct sequence *sequence_ptr;

  FILE *file_ptr;

  static char chain_list[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";

  /* Initialise variables */
  a = 'a';
  z = 'z';
  chain_id = 'A';
  domain_no = 0;
  entry_type = ALL_AT;
  fit_flagged = FALSE;
  iresid = 0;
  ligand_id[0] = '\0'; 
  line = 0;
  nchains = strlen(chain_list);
  nfiles = *nf;
  nseqs = 0;
  nresidues = 0;
  number_string[0] = '\0';

  /* Initialise CORA-format parameters */
  column_width = 10;
  start_pos = 13;
  version = 0;

  /* Open the alignments file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Unable to open alignments file [%s]\n",file_name);
      exit(1);
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
          version = 1;
        }

      /* Process only if this is not a comment line */
      if (input_line[0] != '#')
        {
          /* If this is the first record, read in number of sequences */
          if (line == 0)
            nseqs = atoi(input_line);

          /* If this is the second record, read in the sequence names and
             store only those corresponding to the PDB sequences to be
             displayed */
          else if (line == 1)
            {
              /* Extract the PDB codes on this line */
              ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

              /* If fewer sequences read in than expected, show warning
                 message */
              if (ntokens < nseqs)
                {
                  /* Show warning message */
                  printf("*** Warning. Some sequences not read in from "
                         "CORA file. Increase MAX_TOKEN to %d\n",nseqs);

                  /* Save the number of sequences to be used */
                  nseqs = ntokens;
                }

              /* Loop over all the sequences */
              for (iseq = 0; iseq < nseqs; iseq++)
                {
                  /* Extract the wolf code containing PDB code, chain ID and
                     domain */
                  strcpy(wolf_code,token[iseq]);

                  /* Initialise flag */
                  found = FALSE;
                  if (!strcmp(wolf_code,"noPDB!"))
                    found = TRUE;

                  /* If reading all entries, then create a sequence record
                     for this one */
                  if (params.all_entries == TRUE &&
                      strcmp(wolf_code,"noPDB!") && nfiles < MAXPDB)
                    {
                      /* Get the length of the code */
                      len = strlen(wolf_code);

                      /* Extract the PDB code */
                      strncpy(pdb_code,wolf_code,4);
                      pdb_code[4] = '\0';
                      chain_id = ' ';
                      domain_no = 0;

                      /* Check if there is a chain-id straight after it */
                      i = 5;
                      if (len > 4)
                        {
                          if (wolf_code[4] != '_')
                            {
                              /* Ligand identifier starts one position beyond
                                 chain-id */
                              i = i + 1;
                              chain_id = wolf_code[4];
                              if (chain_id == '0')
                                chain_id = ' ';

                              /* Check if there's a domain number attached */
                              if (len > 5)
                                if (wolf_code[5] != '_')
                                  {
                                    /* Ligand identifier starts one position
                                       beyond domain number */
                                    i = i + 1;
                                    number_string[0] = wolf_code[5];
                                    number_string[1] = '\0';

                                    /* Check for 2-digit domain identifier */
                                    if (len == 7)
                                      {
                                        number_string[1] = wolf_code[6];
                                        number_string[2] = '\0';
                                      }

                                    /* Get the domain number */
                                    domain_no = atoi(number_string);
                                  }
                            }
                        }

                      /* Create a sequence record */
                      create_sequence_record(&sequence_ptr,nfiles,pdb_code,
                                             chain_id,domain_no,ligand_id,
                                             entry_type);

                      /* Store the alignment position */
                      sequence_ptr->alignment_column = iseq;

                      /* Increment count of PDB files */
                      nfiles++;

                      /* Set flag */
                      found = TRUE;
                    }

                  /* Loop through all the PDB codes and determine if any of
                     them match this one */
                  for (ifile = 0; ifile < nfiles && found == FALSE; ifile++)
                    {
                      /* Get the corresponding sequence number for
                         this entry */
                      sequence_ptr = stored_sequence_ptr[ifile];

                      /* Check that PDB code agrees */
                      if (!strncmp(wolf_code,sequence_ptr->pdb_code,4) ||
                          (ifile == 0 && !strncmp(wolf_code,"uplo",4)))
                        {
                          /* Check that the chain-id agrees, too */
                          if (wolf_code[4] == sequence_ptr->chain ||
                              (wolf_code[4] == '0' && 
                               sequence_ptr->chain == ' '))
                            {
                              /* Get the length of the wolf code */
                              len = strlen(wolf_code);
                              domain_no = 0;

                              /* Extract the domain-number, if there
                                 is one */
                              if (len > 5)
                                {
                                  /* Get the first character of the domain
                                     id */
                                  number_string[0] = wolf_code[5];
                                  number_string[1] = '\0';

                                    /* Check for 2-digit domain identifier */
                                    if (len == 7)
                                      {
                                        number_string[1] = wolf_code[6];
                                        number_string[2] = '\0';
                                      }

                                    /* Get the domain number */
                                    domain_no = atoi(number_string);
                                }

                              /* Compare domain number */
                              if (domain_no == 0 ||
                                  domain_no == sequence_ptr->domain_no)
                                    found = TRUE;
                            }
                        }

                      /* If have a match, store column-position in the
                         alignment */
                      if (found == TRUE)
                        sequence_ptr->alignment_column = iseq;
                    }
                }

              /* Free up the tokens */
              free_tokens(token,ntokens);
            }

          /* If this is the third record, read in the total number of 
             residues in the alignment */
          else if (line == 2)
            nresidues = atoi(input_line);

          /* Otherwise, get the alignment at this residue position */
          else if (line > 2 && iresid < nresidues)
            {
              /* Initialise variables for this alignment line */
              reference_list_ptr = NULL;
              nmatch = 0;
              fit = FALSE;

              /* Check to see if the reference sequence has a blank at
                 this position */
              column = 0;

              /* Calculate start of this residue's details */
              ipos = start_pos + column * column_width;

              /* Check whether blank */
              strncpy(residue_details,input_line+ipos,column_width);
              residue_details[column_width] = '\0';

              /* Get the residue name, number and chain-id from the
                 residue details */
              if (version == 0)
                aa_code = residue_details[6];
              else
                aa_code = residue_details[7];

              /* If reference sequence not blank, set flag */
              if (aa_code == '0')
                gap_ref = TRUE;
              else
                gap_ref = FALSE;

              /* Loop through all the stored sequences to see if 
                 any of them have a fit to the reference file at this
                 residue position */
              for (ifile = 1; ifile < nfiles && fit == FALSE; ifile++)
                {
                  /* Get the corresponding sequence record pointer */
                  sequence_ptr = stored_sequence_ptr[ifile];

                  /* Retrieve the alignment column */
                  column = sequence_ptr->alignment_column;

                  /* Calculate start of this residue's details */
                  ipos = start_pos + column * column_width;
                  if (version == 0)
                    ipos = ipos + 8;
                  else
                    ipos = ipos + 9;

                  /* Get the fit-flag */
                  if (input_line[ipos] == 'f' || input_line[ipos] == 'F')
                    fit = TRUE;
                }

              /* Loop through all the stored sequences */
              for (ifile = 0; ifile < nfiles &&
                     (gap_ref == FALSE || params.whole == TRUE); ifile++)
                {
                  /* Get the corresponding sequence record pointer */
                  sequence_ptr = stored_sequence_ptr[ifile];

                  /* Get the current fit list */
                  fitlist_ptr = sequence_ptr->fitlist_ptr;

                  /* Retrieve the alignment column */
                  column = sequence_ptr->alignment_column;

                  /* Calculate start of this residue's details */
                  ipos = start_pos + column * column_width;

                  /* Extract the residue details for this residue */
                  strncpy(residue_details,input_line+ipos,column_width);
                  residue_details[column_width] = '\0';

                  /* Get the residue name, number and chain-id from the
                     residue details */
                  if (version == 0)
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
                  if (version == 0)
                    strncpy(res_num,residue_details,5);
                  else
                    strncpy(res_num,residue_details+1,5);
                  res_num[5] = '\0';
                  if (ifile == 0)
                    strcpy(res_number[0],res_num);
                  else
                    strcpy(res_number[1],res_num);
                  chain_id = sequence_ptr->chain;

                  /* Create a new residue record and attach it to this
                     sequence's list of residues */
                  create_list_record(&list_ptr,sequence_ptr,res_name,
                                     res_num,chain_id);

                  /* If this is a gap, identify as a blank record */
                  if (aa_code == '0')
                    {
                      /* Mark as a gap */
                      list_ptr->blank = TRUE;
                    }

                  /* Update this residue's details */
                  else
                    {
                      /* Store the one-letter residue code */
                      list_ptr->aa_code = aa_code;
                      list_ptr->blank = FALSE;

                      /* Get the residue colour */
                      if (version == 0)
                        strncpy(number_string,residue_details+7,2);
                      else
                        strncpy(number_string,residue_details+8,2);
                      number_string[2] = '\0';

                      /* Check whether this field contains a fit-residue
                         marker */
                      if ((ifile == 0 && fit == TRUE) ||
                          (number_string[1] == 'F' ||
                           number_string[1] == 'f'))
                        {
                          /* Blank out the field */
                          strcpy(number_string,"  ");

                          /* Set flag to fit only on flagged residues */
                          fit_flagged = TRUE;

                          /* Flag this residue as a fitting residue */
                          list_ptr->fit = TRUE;

                          /* If don't yet have a fitlist record, or
                             the current one is complete, then create
                             a new one */
                          if (fitlist_ptr == NULL ||
                              fitlist_ptr->done == TRUE)
                            fitlist_ptr = sequence_ptr->fitlist_ptr
                              = create_fitlist_record(&(sequence_ptr->
                                                        first_fitlist_ptr),
                                                      &(sequence_ptr->
                                                        last_fitlist_ptr));

                          /* If first of the range, then store
                             residue number in start and end positions */
                          if (fitlist_ptr->range[0] == '\0')
                            {
                              /* Store in first position */
                              strncpy(fitlist_ptr->range,res_num,4);
                              fitlist_ptr->range[4] = '\0';
                              strcat(fitlist_ptr->range,"-");

                              /* Store in second position */
                              strncat(fitlist_ptr->range,res_num,4);
                              fitlist_ptr->range[9] = '\0';
                            }

                          /* Otherwise, store in second position */
                          else
                            {
                              strncpy(fitlist_ptr->range+5,res_num,4);
                              fitlist_ptr->range[9] = '\0';
                            }
                        }

                      /* If not a fitting residue, close off any
                         fitting list, if there is one */
                      else if (fitlist_ptr != NULL)
                        fitlist_ptr->done = TRUE;
                      list_ptr->colour = atoi(number_string);

                      /* Number in .aln file encodes for FASTA similarity
                         marker and residue colour, so need to extract
                         just the colour */
                      done = FALSE;
                      while (done == FALSE)
                        {
                          if (list_ptr->colour >= 20)
                            list_ptr->colour = list_ptr->colour - 20;
                          else
                            done = TRUE;
                        }

                      /* If residue-colour is not blank, then unset grey
                         flag for this structure */
                      if (strncmp(number_string,"  ",2))
                        sequence_ptr->grey = FALSE;

                      /* If reference sequence has a blank at this position,
                         then there's no reference residue */
                      if (gap_ref == TRUE)
                        {
                          reference_list_ptr = NULL;
                          if (ifile == 0)
                            list_ptr = NULL;
                        }

                      /* If this is the first structure, then store the
                         residue pointer as the reference residue */
                      if (ifile == 0)
                        reference_list_ptr = list_ptr;
                      else
                        {
                          list_ptr->reference_list_ptr
                            = reference_list_ptr;
                          nmatch++;
                        }
                    }
                }

              /* If no residues other than the reference one in the alignment
                 set its reference pointer to NULL so it is not included in
                 the C of G calculation when performing the fitting later on */
              if (nmatch == 0)
                {
                  if (reference_list_ptr != NULL)
                    reference_list_ptr->reference_list_ptr = NULL;
                }

              /* Otherwise, set the reference pointer to itself */
              else
                {
                  if (reference_list_ptr != NULL)
                    reference_list_ptr->reference_list_ptr
                      = reference_list_ptr;
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

  /* Return number of sequences and total length of the alignment */
  *nf = nfiles;
  *nseq = nseqs;
  *nresid = nresidues;

  /* Return whether fitting to be performed only on flagged residues */
  return(fit_flagged);
}
/***********************************************************************

form_file_names  -  Form the input filenames of the PDB file and ligplot
                    output file giving list of ligand residues and
                    interacting protein residues

***********************************************************************/

void form_file_names(struct sequence *sequence_ptr,char *pdb_name,
                     char *dfn_name,char *reslist_name,
                     char *protein_family,char *matfile_name,
                     char *grow_out_name,char *promotif_dat,int ifile,
                     struct c_file_data file_name_ptr,
                     struct parameters *params)
{
  char ch, chain_id[2], domain_id[10], het_name[6], wolf_code[7];
  char first_char[2], id[FILENAME_LEN], pdbsum_dir[FILENAME_LEN];

  char *fname, *string_ptr;

  int ipdb, len, type, upload;
  int dir, dir_type;

  FILE *file_ptr;

  /* Initialise variables */
  dfn_name[0] = '\0';
  dir_type = NOT_FOUND;
  strcpy(pdbsum_dir,"./");
  type = params->pdb_type;
  upload = FALSE;

  /* If using Wolf code, rather than PDB code, form it */
  if (params->use_wolf == TRUE)
    {
      /* Copy the PDB code */
      strcpy(wolf_code,sequence_ptr->pdb_code);

      /* Add the chain-identifier */
      if (sequence_ptr->chain == ' ')
        chain_id[0] = '0';
      else
        chain_id[0] = sequence_ptr->chain;
      chain_id[1] = '\0';
      strcat(wolf_code,chain_id);

      /* Add domain number  */
      sprintf(domain_id,"%d",sequence_ptr->domain_no);
      strcat(wolf_code,domain_id);
    }

  /* For an uploaded file, use the file.new file in the given
     working directory */
  if (!strcmp(sequence_ptr->pdb_code,"uplo") ||
      params->profunc == TRUE)
    {
      /* If have a working directory, use that */
      if (params->work_dir[0] != '\0')
        {
          /* Define name of PDB input file */
          strcpy(pdb_name,params->work_dir);
          strcat(pdb_name,DIR_SEP);

          /* Define name of rasmol.dfn file */
          strcpy(dfn_name,params->work_dir);
          strcat(dfn_name,DIR_SEP);
        }
      else
        {
          pdb_name[0] = '\0';
          dfn_name[0] = '\0';
        }

      /* Add file name */
      strcat(pdb_name,"file.new");
      strcat(dfn_name,"rasmol.dfn");

      /* Set flag that first file is an uploaded file */
      params->upload = TRUE;
    }

  /* For PDBsum1, file is in directory above */
  else if (params->pdb_type == PDBSUM1)
    {
      /* PDB file */
      sprintf(pdb_name,"../%s.pdb",sequence_ptr->pdb_code);
      if (params->pdb_file != &null)
        sprintf(pdb_name,"%s",params->pdb_file);

      /* rasmol.dfn file */
      strcpy(dfn_name,"../newsec/rasmol.dfn");
    }

  /* Otherwise, use given PDB code */
  else
    {
      /* See if this is an Alpha Fold model */
      ch = sequence_ptr->pdb_code[0];
      if (ch >= 'A' && ch <= 'Z')
        type = ALPHA_FOLD_PDB;
      else
        type = STANDARD_PDB;

      /* Assume standard PDB file */
      dir = r_find_pdb_file(sequence_ptr->pdb_code,
                            file_name_ptr.pdb_template,pdb_name,
                            type);

      /* If not found, try PDBsum upload directory */
      if (dir < 0)
        {
          type = UPLOAD_PDB;
          dir = r_find_pdb_file(sequence_ptr->pdb_code,
                                file_name_ptr.pdb_template,pdb_name,
                                type);
        }
    }

  /* Form the name of this PDB entry's PDBsum directory */
  dir_type = find_pdbsum_dir(sequence_ptr->pdb_code,
                             file_name_ptr.pdbsum_template,
                             type,pdbsum_dir);

  /* If this entry is for the binding-site residues, then link to
     appropriate family-name entry */
  if (sequence_ptr->display_option == BINDING_SITE)
    {
      /* Append the protein family name */
      strcpy(reslist_name,pdbsum_dir);
      strcat(reslist_name,protein_family);

      /* Add the chain-identifier, if there is one */
      if (sequence_ptr->chain != ' ')
        {
          chain_id[0] = sequence_ptr->chain;
          chain_id[1] = '\0';
          strcat(reslist_name,chain_id);
        }

      /* Add domain number to file name, if not zero */
      if (sequence_ptr->domain_no != 0)
        {
          sprintf(domain_id,"_%d",sequence_ptr->domain_no);
          strcat(reslist_name,domain_id);
        }

      /* Add .res extension */
      strcat(reslist_name,".res");
    }

  /* Otherwise, link to the appropriate LIGPLOT .res file */
  else
    {
      /* Form name of .res file */
      strcpy(reslist_name,pdbsum_dir);

      /* Check whether a Het Group name, ligand or metal */
      if (!strncmp(sequence_ptr->ligand_id,"het",3))
        {
          /* Get the name of the Het Group */
          strcpy(het_name,sequence_ptr->ligand_id+3);

          /* Get the first character */
          first_char[0] = het_name[0];
          first_char[1] = '\0';

          /* Use .res file in the Het Groups gif directory */
          strcpy(reslist_name,file_name_ptr.other_dir[HETGIF_REAL_DIR]);
          strcat(reslist_name,DIR_SEP);
          strcat(reslist_name,first_char);
          strcat(reslist_name,DIR_SEP);
          strcat(reslist_name,het_name);
          strcat(reslist_name,DIR_SEP);
        }
      else if (params->metal == TRUE)
        {
          strcat(reslist_name,DIR_SEP);
          strcat(reslist_name,"metplot");
        }
      else
        {
          strcat(reslist_name,DIR_SEP);
          strcat(reslist_name,"ligplot");
        }
      strcat(reslist_name,sequence_ptr->ligand_id);
      strcat(reslist_name,".res");

      /* Form name of grow.out file */
      strcpy(grow_out_name,pdbsum_dir);
      strcat(grow_out_name,"grow.out");

      /* Form name of promotif.dat file */
      strcpy(promotif_dat,pdbsum_dir);
      strcat(promotif_dat,"promotif.dat");

      /* For PDBsum1, point to the relevant files */
      if (params->pdb_type == PDBSUM1)
        {
          strcpy(promotif_dat,"../promotif/promotif.dat");
          strcpy(grow_out_name,"../newsec/grow.out");

          /* If have a PDB file name, form id */
          if (params->pdb_file != &null)
            {
              strcpy(id,params->pdb_file);
              string_ptr = strstr(id,".pdb");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              strcat(id,".res");
              strcpy(reslist_name,id);
            }
        }
    }

  /* Form name of rasmol.dfn file if we don't already have it*/
  if (dfn_name[0] == '\0')
    {
      strcpy(dfn_name,pdbsum_dir);
      strcat(dfn_name,"rasmol.dfn");
    }

  /* Form name of matrix file for the given family that holds the
     transformation of the current structure to the reference structure
     for this family */
  if (protein_family[0] == '\0')
    matfile_name[0] = '\0';
  else
    {
      strcpy(matfile_name,pdbsum_dir);
      strcat(matfile_name,protein_family);

      /* If have a chain-identifier, then add that to the name */
      if (sequence_ptr->chain != ' ')
        {
          chain_id[0] = sequence_ptr->chain;
          chain_id[1] = '\0';
          strcat(matfile_name,chain_id);
        }

      /* If have a domain number, add it to the name */
      if (sequence_ptr->domain_no != 0)
        {
          sprintf(domain_id,"%d",sequence_ptr->domain_no);
          strcat(matfile_name,"_");
          strcat(matfile_name,domain_id);
        }

      /* Add the .mat extension */
      strcat(matfile_name,".mat");
    }
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

get_dfn_data  -  Read in and store the RasMol definitions from the
                 rasmol.dfn file

***********************************************************************/

int get_dfn_data(char *pdb_code,char *dfn_name,
                 struct rasdef **fst_rasdef_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, last_chain, res_name[4], res_num[6];
  char colour[COLNAM_LEN];
  char ligand_from[7], ligand_to[7];

  int count, line, ndef, nligand, nmetal, type;
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
  nligand = 0;
  nmetal = 0;

  /* Initialise rasdef pointers */
  rasdef_ptr = NULL;
  first_rasdef_ptr = NULL;
  last_rasdef_ptr = NULL;

  /* Form filename of the rasmol.dfn file */
  strcpy(file_name,dfn_name);

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

  /* Return whether definitions read in */
  return(have_dfn);
}
/***********************************************************************

get_prosite_strings  -  Get the residue- and colour-strings from the
                        PROSITE pattern line

***********************************************************************/

void get_prosite_strings(char input_line[LINELEN + 1],
                         char residue_string[STRING_LEN],
                         char colours_string[STRING_LEN])
{
  char ch;

  int ipos, ncols, nres, spos;
  int colours, done;

  /* Initialise variables */
  colours = FALSE;
  done = FALSE;
  ipos = 63;
  ncols = 0;
  nres = 0;
  spos = 0;
  while (done == FALSE)
    {
      /* Get the next character from the input line */
      ch = input_line[ipos];

      /* If this is the end of the line, then end here */
      if (ipos > LINELEN || ch == '\0' || ch == '\n')
        done = TRUE;

      /* Otherwise, if non-blank, then store */
      else if (ch != ' ' && ch != '\t' && ch != '?' && ch != '>'
               && ch != '<' && ch != '!' && ch != '|')
        {
          /* Check that string not already too long */
          if (spos < STRING_LEN - 1)
            {
              /* Store one-letter residue code, or colour */
              if (colours == FALSE)
                {
                  residue_string[spos] = ch;
                  nres++;
                }
              else
                {
                  colours_string[spos] = ch;
                  ncols++;
                }

              /* Increment position in output string */
              spos++;
            }
        }

      /* If it is an exclamation mark, need to delete previous
         character as well */
      else if (ch == '!')
        {
          /* If have a previous character, then delete */
          if (spos > 0)
            spos--;
        }

      /* If it is a blank or a TAB, determine what to do about it */
      else if (ch == ' ' || ch == '\t')
        {
          /* If have been reading in the amino acid sequence,
             then flip flag to expect the colours sequence next */
          if (colours == FALSE && nres > 0)
            {
              colours = TRUE;
              spos = 0;
            }
        }

      /* Go to the next input character position */
      ipos++;
    }

  /* Terminate both strings */
  residue_string[nres] = '\0';
  colours_string[ncols] = '\0';
}
/***********************************************************************

check_for_duplicate_pattern  -  Check that we don't already have this
                                PROSITE pattern

***********************************************************************/

int check_for_duplicate_pattern(struct prosite_pattern 
                                *first_prosite_pattern_ptr,
                                char prosite_id[PROSITE_PATLEN],
                                char chain,char first_res_num[6],
                                char last_res_num[6],
                                char residue_string[STRING_LEN],
                                char colours_string[STRING_LEN])
{
  int wanted;

  struct prosite_pattern *prosite_pattern_ptr;

  /* Initialise variables */
  wanted = TRUE;

  /* Start at first pattern */
  prosite_pattern_ptr = first_prosite_pattern_ptr;

  /* Loop through all the patterns */
  while (prosite_pattern_ptr != NULL && wanted == TRUE)
    {
      /* Check if have a duplicate */
      if (!strcmp(prosite_pattern_ptr->prosite_id,prosite_id) &&
          prosite_pattern_ptr->chain == chain &&
          !strcmp(prosite_pattern_ptr->first_res_num,first_res_num) &&
          !strcmp(prosite_pattern_ptr->last_res_num,last_res_num) &&
          !strcmp(prosite_pattern_ptr->residue_string,residue_string) &&
          !strcmp(prosite_pattern_ptr->colours_string,colours_string))
        wanted = FALSE;

      /* Go to next pattern */
      prosite_pattern_ptr = prosite_pattern_ptr->next_prosite_pattern_ptr;
    }

  /* Return whether pattern is wanted */
  return(wanted);
}
/***********************************************************************

get_prosite_patterns  -  Pick up all this sequence's PROSITE required
                         patterns

***********************************************************************/

void get_prosite_patterns(struct sequence *sequence_ptr,
                          struct c_file_data file_name_ptr,
                          struct parameters params)
{
  char input_line[LINELEN + 1];
  char chain, last_chain;
  char first_res_num[6], last_res_num[6];
  char last_prosite_id[PROSITE_PATLEN], prosite_id[PROSITE_PATLEN];
  char colours_string[STRING_LEN], residue_string[STRING_LEN];
  char file_name[FILENAME_LEN];

  int occurrence, position;
  int wanted;

  struct match_pattern *match_pattern_ptr;
  struct prosite_pattern *prosite_pattern_ptr, *last_prosite_pattern_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  last_chain = '\0';
  occurrence = 0;
  last_prosite_id[0] = '\0';
  last_prosite_pattern_ptr = NULL;

  /* Get appropriate PROSITE file */
  if (params.pdb_type == STANDARD_PDB)
    strcpy(file_name,file_name_ptr.other_dir[PROSITE_FILE]);
  else
    strcpy(file_name,file_name_ptr.other_dir[UPLOAD_PROSITE_FILE]);

  /* Open the PROSITE matches file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If generating a RasMol script, then don't display error message */
      if (params.rasmol_script == TRUE)
        return;
      else
        {
          printf("*** ERROR. File not found: %s\n",file_name);
          exit(1);
        }
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Get the chain-id */
      chain = input_line[6];

      /* Check for PDB code and chain-id */
      if (!strncmp(input_line,sequence_ptr->pdb_code,4) &&
          (sequence_ptr->chain == ' ' || 
           chain == sequence_ptr->chain))
        {
          /* Get this pattern-id */
          strncpy(prosite_id,input_line+10,7);
          prosite_id[7] = '\0';

          /* If this is not the same chain, or the same pattern, as the
             previous pattern for this structure, then re-initialise
             the occurrence count */
          if (chain != last_chain || strncmp(prosite_id,last_prosite_id,7))
            occurrence = 0;

          /* Otherwise, increment the occurrence count */
          occurrence++;

          /* Save the pattern-id and chain */
          strncpy(last_prosite_id,prosite_id,7);
          last_prosite_id[7] = '\0';
          last_chain = chain;

          /* Initialise flag */
          wanted = FALSE;

          /* Set pointer to the first of the PROSITE patterns to be
             matched */
          match_pattern_ptr = first_match_pattern_ptr;

          /* Loop through all the PROSITE patterns */
          while (match_pattern_ptr != NULL && wanted == FALSE)
            {
              /* Check if the pattern-id matches */
              if (!strncmp(prosite_id,match_pattern_ptr->prosite_id,7))
                {
                  /* Check the occurrence count, if required */
                  if (match_pattern_ptr->occurrence == 0)
                    wanted = TRUE;
                  else if (occurrence == match_pattern_ptr->occurrence)
                    wanted = TRUE;
                }

              /* Get pointer to the next pattern */
              match_pattern_ptr = match_pattern_ptr->next_match_pattern_ptr;
            }

          /* If this pattern is wanted then check that we don't already
             have it */
          if (wanted == TRUE)
            {
              /* Get start- and end-residue numbers */
              strncpy(first_res_num,input_line+51,5);
              first_res_num[5] = '\0';
              strncpy(last_res_num,input_line+57,5);
              last_res_num[5] = '\0';

              /* Get sequential residue position */
              sscanf(input_line+46,"%d",&position);

              /* Get the strings giving the one-letter residue codes and
                 their colours */
              get_prosite_strings(input_line,residue_string,colours_string);

              /* Check that we don't already have this pattern */
              wanted = 
                check_for_duplicate_pattern(sequence_ptr->
                                            first_prosite_pattern_ptr,
                                            prosite_id,chain,first_res_num,
                                            last_res_num,residue_string,
                                            colours_string);
            }

          /* If still wanted, then store details */
          if (wanted == TRUE)
            {
              /* Allocate memory for this PROSITE pattern */
              prosite_pattern_ptr = (struct prosite_pattern *) 
                malloc(sizeof(struct prosite_pattern));
              if (prosite_pattern_ptr == NULL)
                {
                  printf("*** Can't allocate memory for struct ");
                  printf("prosite_pattern\n");
                  exit (1);
                }

              /* Add link from previous pattern record to the current one */
              if (last_prosite_pattern_ptr != NULL)
                last_prosite_pattern_ptr->next_prosite_pattern_ptr
                  = prosite_pattern_ptr;
              else
                sequence_ptr->first_prosite_pattern_ptr = prosite_pattern_ptr;
              last_prosite_pattern_ptr = prosite_pattern_ptr;

              /* Store the pattern-id and occurrence number */
              strncpy(prosite_pattern_ptr->prosite_id,prosite_id,7);
              prosite_pattern_ptr->prosite_id[7] = '\0';
              prosite_pattern_ptr->occurrence = occurrence;
              prosite_pattern_ptr->chain = chain;

              /* Get start- and end-residue numbers */
              strcpy(prosite_pattern_ptr->first_res_num,first_res_num);
              strcpy(prosite_pattern_ptr->last_res_num,last_res_num);

              /* Get sequential residue position */
              prosite_pattern_ptr->position = position;

              /* Get the strings giving the one-letter residue codes and
                 their colours */
              strcpy(prosite_pattern_ptr->residue_string,residue_string);
              strcpy(prosite_pattern_ptr->colours_string,colours_string);

              /* Initialise pointer to next pattern */
              prosite_pattern_ptr->next_prosite_pattern_ptr = NULL;
            }          
        }
    }

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

write_script_headers  -  Write all the header records to the output
                         RasMol script file

***********************************************************************/

void write_script_headers(FILE *outfile_ptr,int nfiles,
                          struct c_file_data file_name_ptr,
                          struct parameters params)
{
  char *pdbsum_dir;

  int dir_type, type;

  struct match_pattern *match_pattern_ptr;

  /* Initialise variables */
  first_sequence_ptr = stored_sequence_ptr[0];

  /* Write out the starting lines */
  fprintf(outfile_ptr,"#!rasmol -script\n");
  fprintf(outfile_ptr,"# Generated by PDBsum - R A Laskowski, Feb 1998\n");
  fprintf(outfile_ptr,"\n");

  /* Show PDB code and PDB file name */
  if (nfiles == 1)
    fprintf(outfile_ptr,"# PDB code: %s\n",first_sequence_ptr->pdb_code);
  else
    fprintf(outfile_ptr,"# Multiple PDB files\n");
  fprintf(outfile_ptr,"\n");

  /* Initilise and indicate that ATOM records to follow */
  fprintf(outfile_ptr,"zap\n");

  /* For PDBsum1, read PDB file on disk */
  if (params.pdb_type == PDBSUM1)
    {
      /* Get the name of the PDBsum directory */
      pdbsum_dir
        = form_filename(first_sequence_ptr->pdb_code,
                        file_name_ptr.other_dir[PDBSUM_TEMPLATE]);

      /* Write command to load the PDB file */
      if (params.pdb_file != &null)
        fprintf(outfile_ptr,"load inline\n");
      else
        fprintf(outfile_ptr,"load %s/%s.pdb\n",pdbsum_dir,
                first_sequence_ptr->pdb_code);
    }
  else
    fprintf(outfile_ptr,"load inline\n");

  /* Set background colour according to where prog has been called from */
  if (params.binding_sites == TRUE)
    fprintf(outfile_ptr,"background [255,204,102]\n");
  else if (params.csa == TRUE)
    fprintf(outfile_ptr,"background white\n");
  else if (params.residue_list[0] != '\0')
    if (!strcmp(params.rlist_colour,"exons"))
      fprintf(outfile_ptr,"background white\n");
    else
      fprintf(outfile_ptr,"background black\n");
  else if (params.cath_alignment == TRUE)
    fprintf(outfile_ptr,"background black\n");
  else if (params.prosite_pattern == TRUE)
    fprintf(outfile_ptr,"background black\n");
  else if (params.promotif == TRUE)
    fprintf(outfile_ptr,"background black\n");
  else if (params.pdb_site == TRUE)
    fprintf(outfile_ptr,"background black\n");
  else
    fprintf(outfile_ptr,"background [0,0,136]\n");
  /*    fprintf(outfile_ptr,"background [153,153,255]\n"); */
  /*    fprintf(outfile_ptr,"background [255,255,179]\n"); */
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
  if (params.prosite_pattern == TRUE)
    {
      /* Initialise pointer to the very first PROSITE pattern */
      match_pattern_ptr = first_match_pattern_ptr;

      /* Show start of display panel */
      fprintf(outfile_ptr,"echo \"   +-----------------------------+\"\n");

      /* If no PROSITE patterns selected, then say so */
      if (match_pattern_ptr == NULL)
        fprintf(outfile_ptr,"echo \"   | No PROSITE pattern selected |\"\n");

      /* Loop through all the PROSITE patterns */
      while (match_pattern_ptr != NULL)
        {
          /* Show the pattern */
          if (match_pattern_ptr->occurrence == 0)
            fprintf(outfile_ptr,"echo \"   | PROSITE pattern %s     |\"\n",
                    match_pattern_ptr->prosite_id);
          else
            fprintf(outfile_ptr,"echo \"   | PROSITE pattern %s /%2d |\"\n",
                    match_pattern_ptr->prosite_id,
                    match_pattern_ptr->occurrence);

          /* Get pointer to the next pattern */
          match_pattern_ptr = match_pattern_ptr->next_match_pattern_ptr;
        }

      /* Close off the panel */
      fprintf(outfile_ptr,"echo \"   +-----------------------------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"echo \"* Residues corresponding to the ");
      fprintf(outfile_ptr,"PROSITE pattern are shown in colour\"\n");

      /* Write out the common lines */
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"\n");
    }
  else if (params.binding_sites == TRUE)
    {
      fprintf(outfile_ptr,"echo \"   +---------------+\"\n");
      fprintf(outfile_ptr,"echo \"   | Binding sites |\"\n");
      fprintf(outfile_ptr,"echo \"   +---------------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"\n");
    }
  else if (params.pdb_site == TRUE)
    {
      fprintf(outfile_ptr,"echo \"   +-------------+\"\n");
      fprintf(outfile_ptr,"echo \"   | Active site |\"\n");
      fprintf(outfile_ptr,"echo \"   +-------------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"\n");
    }
  else if (params.csa == TRUE)
    {
      fprintf(outfile_ptr,"echo \"   +------------------------+\"\n");
      fprintf(outfile_ptr,"echo \"   | CSA catalytic residues |\"\n");
      fprintf(outfile_ptr,"echo \"   +------------------------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"\n");
    }
  else if (params.promotif == TRUE)
    {
      fprintf(outfile_ptr,"echo \"   +----------+\"\n");
      fprintf(outfile_ptr,"echo \"   | PROMOTIF |\"\n");
      fprintf(outfile_ptr,"echo \"   +----------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"\n");
    }
  else if (params.cath_alignment == FALSE)
    {
      fprintf(outfile_ptr,"echo \"   +-----------------+\"\n");
      fprintf(outfile_ptr,"echo \"   | Ligand in %s |\"\n",
              first_sequence_ptr->pdb_code);
      fprintf(outfile_ptr,"echo \"   +-----------------+\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"\n");
    }

  /* Output start of list of user-defined objects in this RasMol script */
  fprintf(outfile_ptr,"echo \"* User-defined sets present in structure:-\"\n");
  fprintf(outfile_ptr,"echo \" \"\n");
}
/***********************************************************************

write_pymol_headers  -  Write the header records for the PyMOL script

***********************************************************************/

void write_pymol_headers(FILE *outfile_ptr,int nfiles,
                         struct c_file_data file_name_ptr,
                         struct parameters params)
{
  char *pdbsum_dir;

  int dir_type, type;

  struct match_pattern *match_pattern_ptr;

  /* Initialise variables */
  first_sequence_ptr = stored_sequence_ptr[0];

  /* Write out the starting lines */
  fprintf(outfile_ptr,"#\n");
  fprintf(outfile_ptr,"# PDBsum1 PyMol script for PDB code %s\n",
          first_sequence_ptr->pdb_code);
  fprintf(outfile_ptr,"#\n");

  /* Get the name of the PDBsum directory */
  pdbsum_dir
    = form_filename(first_sequence_ptr->pdb_code,
                    file_name_ptr.other_dir[PDBSUM_TEMPLATE]);

  /* Load the structure */
  fprintf(outfile_ptr,"load %s/%s.pdb\n",pdbsum_dir,
          first_sequence_ptr->pdb_code);

  /* Write out the header records */
  fprintf(outfile_ptr,"bg_color white\n");
  fprintf(outfile_ptr,"set depth_cue=0\n");
  fprintf(outfile_ptr,"set ray_trace_fog=0\n");
  fprintf(outfile_ptr,"set cartoon_ring_mode,1\n");
  fprintf(outfile_ptr,"hide all\n");
  fprintf(outfile_ptr,"set auto_color, 0\n");

  /* Define the colours */
  fprintf(outfile_ptr,"set_color black, [  0.000,  0.000,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color red, [  1.000,  0.000,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color green, [  0.000,  1.000,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color blue, [  0.000,  0.000,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color yellow, [  1.000,  1.000,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color purple, [  0.800,  0.400,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color brown, [  0.800,  0.600,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color sky_blue, [  0.000,  0.600,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color skyblue, [  0.000,  0.600,  1.000 ]\n");
  fprintf(outfile_ptr,"set_color orange, [  1.000,  0.400,  0.000 ]\n");
  fprintf(outfile_ptr,"set_color cyan, [  0.200,  0.800,  1.000 ]\n");
}
/***********************************************************************

get_transformation  -  Read in this structure's transformation to fit it
                       onto the reference structure of the given family

***********************************************************************/

int get_transformation(char *matfile_name,double matrix[4][4])
{
  char input_line[LINELEN + 1];

  int i, j, row, status;
  int havr, havm, havrot, transform;

  double mobile_cog[3], ref_cog[3], rot_matrix[3][3];

  FILE *file_ptr;

  /* Initialise */
  havr = FALSE;
  havm = FALSE;
  havrot = FALSE;
  row = 0;
  status = 0;
  transform = FALSE;
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
        {
          if (i == j)
            matrix[i][j] = 1.0;
          else
            matrix[i][j] = 0.0;
        }
    }

  /* Open the input file */
  if ((file_ptr = fopen(matfile_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* If looking for the first key-word line, check if we have it */
          if (status == 0)
            {
              /* Reference C of G */
              if (!strncmp(input_line+3,"Reference CofG",14))
                status = 1;

              /* Mobile C of G */
              if (!strncmp(input_line+3,"Mobile CofG",11))
                status = 2;

              /* Rotation matrix */
              if (!strncmp(input_line+3,"Rotation matrix",15))
                {
                  status = 3;
                  row = 0;
                }
            }

          /* If current line is Reference C of G, then store */
          else if (status == 1)
            {
              /* Store the 3 CofG values */
              sscanf(input_line,"%lf %lf %lf",&ref_cog[0],&ref_cog[1],
                     &ref_cog[2]);

              /* Set flag and reset status */
              havr = TRUE;
              status = 0;
            }

          /* If current line is Mobile C of G, then store */
          else if (status == 2)
            {
              /* Store the 3 CofG values */
              sscanf(input_line,"%lf %lf %lf",&mobile_cog[0],&mobile_cog[1],
                     &mobile_cog[2]);

              /* Set flag and reset status */
              havm = TRUE;
              status = 0;
            }

          /* If current line is Rotation matrix, then store */
          else if (status == 3)
            {
              /* Store the 3 CofG values */
              sscanf(input_line,"%lf %lf %lf",&rot_matrix[0][row],
                     &rot_matrix[1][row],&rot_matrix[2][row]);

              /* Increment row */
              row++;

              /* If have done last row of the matrix, then reset status */
              if (row > 2)
                {
                  havrot = TRUE;
                  status = 0;
                }
            }
        }

      /* If have read in all the required data, then can calculate
         transformation matrix */
      if (havr == TRUE && havm == TRUE && havrot == TRUE)
        {
          /* Set flag that we have a valid transformation */
          transform = TRUE;

          /* Transfer the rotation matrix across into the final matrix */
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
              matrix[i][j] = rot_matrix[i][j];
          
          /* Fill in the final column */
          for (i = 0; i < 3; i++)
            matrix[i][3] = 0.0;
          matrix[3][3] = 1.0;

          /* Calculate bottom line of final matrix */
          for (j = 0; j < 3; j++)
            matrix[3][j] = - ref_cog[0] * matrix[0][j]
              - ref_cog[1] * matrix[1][j]
              - ref_cog[2] * matrix[2][j] + mobile_cog[j];
        }
    }

  /* Return whether have a valid transformation matrix */
  return(transform);
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
  if (verbose == TRUE)
    {
      printf("Rotation matrix:-\n");
      for(i=0;i<3;i++)
        {
          printf("     ");
          for(j=0;j<3;j++)
            printf("%10.4f ",rm[j][i]);
          printf("\n");
        }
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
  11.02.13 Return of most significant outliers (RAL)
*/
double matfit(COOR    *x1,        /* First coord array    */
              COOR    *x2,        /* Second coord array   */
              REAL    rm[3][3],   /* Rotation matrix      */
              int     n,          /* Number of points     */
              REAL    *wt1,       /* Weight array         */
              BOOL    column,     /* Column-wise output   */
              double  *atom_rmsd) /* rmsd of each atom   */
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

      /* Save this atom's rmsd */
      atom_rmsd[ipoint] = dist;
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

calculate_fit  -  Calculate the best fit between the current structure
                  and the reference structure (ie the first in the list)

***********************************************************************/

int calculate_fit(struct sequence *sequence_ptr,double matrix[4][4],
                  double *fit_cog,double *rmsdev,int *ns,int fit_flagged,
                  struct parameters params)
{
  int i, ifit, ipos, istored, j, loop, nfit, nfree, nover, nstored;
  int done, transform, used, wanted;

  double mean_x, mean_y, mean_z, mean_fit_x, mean_fit_y, mean_fit_z;
  double atom_rmsd[MAX_FIT_ATOMS], rmsd, rmsd_cutoff;

  BOOL column;

  REAL  rot_matrix[3][3];

  COOR coords1[MAX_FIT_ATOMS], coords2[MAX_FIT_ATOMS];

  struct atom *atom_ptr, *reference_atom_ptr;
  struct residue *residue_ptr, *reference_residue_ptr;
  struct residue *residue_index_ptr[MAX_FIT_ATOMS];

  /* Initialise variables */
  rmsd_cutoff = 0.0;
  if (verbose == TRUE)
    {
      printf("Calculating fit for %s%c%d\n",sequence_ptr->pdb_code,
             sequence_ptr->chain,sequence_ptr->domain_no);
      if (fit_flagged == TRUE)
        printf("  ... only fitting flagged residues\n");
    }

  /* Initialise the matrix to the identity matrix */
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
        {
          if (i == j)
            matrix[i][j] = 1.0;
          else
            matrix[i][j] = 0.0;
        }
    }

  /* If fitting to backbone only, then need to store CA-atoms only */
  if (params.fit_to_backbone == TRUE)
    nfit = 1;
  else
    nfit = 4;

  /* Initialise flag */
  done = FALSE;
  loop = 0;

  /* Loop until we have a reasonable structural alignment */
  while (done == FALSE && loop < MAX_LOOPS)
    {
      /* Initialise values */
      mean_x = 0.0;
      mean_y = 0.0;
      mean_z = 0.0;
      mean_fit_x = 0.0;
      mean_fit_y = 0.0;
      mean_fit_z = 0.0;
      nstored = 0;
      rmsd = 0.0;
      column = FALSE;
      transform = FALSE;

      /* Get pointer to this structure's first residue */
      residue_ptr = sequence_ptr->first_residue_ptr;

      /* Loop through all the residues to pick up the aligned residues */
      while (residue_ptr != NULL)
        {
          /* Unset the fit-used flag */
          residue_ptr->fit_used = FALSE;

          /* Determine whether this residue is wanted for fitting */
          wanted = TRUE;
          if (fit_flagged == TRUE && residue_ptr->fit == FALSE)
            wanted = FALSE;

          /* If the first loop, or the rmsd from the previous loop was
             below the cut-off, zero this residue's rmsd value */
          if (loop == 0 || residue_ptr->rmsd < rmsd_cutoff)
            residue_ptr->rmsd = 0;

          /* Otherwise, residue is not wanted */
          else
            wanted = FALSE;

          /* If residue is to be included in the fitting, then get its C-alpha
             atom to fit on */
          if (wanted == TRUE && residue_ptr->reference_residue_ptr != NULL)
            {
              /* Get the pointer to the reference residue */
              reference_residue_ptr = residue_ptr->reference_residue_ptr;

              /* Initialise flag */
              used = FALSE;

              /* Loop through the atoms to be fitted between these two
                 residues */
              for (ifit = 0; ifit < nfit; ifit++)
                {
                  /* Get the pointer to the residue's relevant atom, and 
                     to its equivalent on the reference residue */
                  if (ifit == 0)
                    {
                      atom_ptr = residue_ptr->CA_atom_ptr;
                      reference_atom_ptr = reference_residue_ptr->CA_atom_ptr;
                    }
                  else if (ifit == 1)
                    {
                      atom_ptr = residue_ptr->N_atom_ptr;
                      reference_atom_ptr = reference_residue_ptr->N_atom_ptr;
                    }
                  else if (ifit == 2)
                    {
                      atom_ptr = residue_ptr->C_atom_ptr;
                      reference_atom_ptr = reference_residue_ptr->C_atom_ptr;
                    }
                  else if (ifit == 3)
                    {
                      atom_ptr = residue_ptr->O_atom_ptr;
                      reference_atom_ptr = reference_residue_ptr->O_atom_ptr;
                    }

                  /* If both atoms present, then extract their coords */
                  if (atom_ptr != NULL && reference_atom_ptr != NULL &&
                      nstored < MAX_FIT_ATOMS)
                    {
                      /* Store the residue pointer */
                      residue_index_ptr[nstored] = residue_ptr;
                      atom_rmsd[nstored] = 0;

                      /* Get the atom's coordinates */
                      coords1[nstored].x = atom_ptr->x;
                      coords1[nstored].y = atom_ptr->y;
                      coords1[nstored].z = atom_ptr->z;
                      coords2[nstored].x = reference_atom_ptr->x;
                      coords2[nstored].y = reference_atom_ptr->y;
                      coords2[nstored].z = reference_atom_ptr->z;

                      /* Accumulate mean values */
                      mean_x = mean_x + atom_ptr->x;
                      mean_y = mean_y + atom_ptr->y;
                      mean_z = mean_z + atom_ptr->z;
                      mean_fit_x = mean_fit_x + reference_atom_ptr->x;
                      mean_fit_y = mean_fit_y + reference_atom_ptr->y;
                      mean_fit_z = mean_fit_z + reference_atom_ptr->z;
        
                      /* Increment count of stored coords */
                      nstored++;
                      used = TRUE;
                    }
                }

              /* Mark both residues as having been used in the fitting */
              if (used == TRUE)
                {
                  residue_ptr->fit_used = TRUE;
                  reference_residue_ptr->fit_used = TRUE;
                }
            }

          /* Get pointer to the next residue in the linked list */
          residue_ptr = residue_ptr->next_residue_ptr;
        }

      /* Calculate centres of mass of the two sets of stored C-alpha atom
         coordinates */
      if (nstored >= MIN_STORED)
        {
          /* Set flag to indicate that we have a transformation */
          transform = TRUE;

          /* Calculate mean values */
          mean_x = mean_x / nstored;
          mean_y = mean_y / nstored;
          mean_z = mean_z / nstored;
          mean_fit_x = mean_fit_x / nstored;
          mean_fit_y = mean_fit_y / nstored;
          mean_fit_z = mean_fit_z / nstored;

          /* Store the mean values of the structure being fitted to give the
             translation to be applied to it */
          matrix[3][0] = mean_x;
          matrix[3][1] = mean_y;
          matrix[3][2] = mean_z;

          /* Loop through all the stored coordinates to centre them at
             the origin */
          for (istored = 0; istored < nstored; istored++)
            {
              /* Adjust both sets of coords */
              coords1[istored].x = coords1[istored].x - mean_x;
              coords1[istored].y = coords1[istored].y - mean_y;
              coords1[istored].z = coords1[istored].z - mean_z;
              coords2[istored].x = coords2[istored].x - mean_fit_x;
              coords2[istored].y = coords2[istored].y - mean_fit_y;
              coords2[istored].z = coords2[istored].z - mean_fit_z;
            }
          if (verbose == TRUE)
            {
              printf("   Ref CofG:  %10.4f %10.4f %10.4f \n",
                     mean_x,mean_y,mean_z);
              printf("   Mob CofG:  %10.4f %10.4f %10.4f \n",
                     mean_fit_x,mean_fit_y,mean_fit_z);
            }

          /* Perform the least-squares fit */
          rmsd = matfit(coords1,coords2,rot_matrix,nstored,NULL,column,
                        atom_rmsd);

          /* Define the new rmsd cut off */
          rmsd_cutoff = (RMSD_FACTOR * rmsd) * (RMSD_FACTOR * rmsd);

          /* Initialise count of atoms over the rmsd cut-off */
          nover = 0;

          /* Transfer the atom-rmsd values to each residue */
          for (istored = 0; istored < nstored; istored++)
            {
              /* Get this residue */
              residue_ptr = residue_index_ptr[istored];

              /* If atom's rmsd is greater than the stored rmsd,
                 then overwrite */
              if (atom_rmsd[istored] > residue_ptr->rmsd)
                residue_ptr->rmsd = atom_rmsd[istored];

              /* If value is over the cut-off, increment count */
              if (residue_ptr->rmsd > rmsd_cutoff)
                nover++;
            }

          /* If no atoms over the rmsd cut-off, then can end here */
          if (nover == 0)
            done = TRUE;

          /* Store the rotation matrix obtained */
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
              matrix[i][j] = rot_matrix[i][j];
        }
          
      /* Increment loop counter */
      loop++;
    }

/* Store the C of G of the fitted portion of the reference structure */
  fit_cog[0] = mean_fit_x;
  fit_cog[1] = mean_fit_y;
  fit_cog[2] = mean_fit_z;

  /* Return the RMSD value and whether have a valid transformation
     matrix */
  *rmsdev = rmsd;
  *ns = nstored;
  return(transform);
}
/***********************************************************************

apply_transformation  -  Apply the given transformation to all the
                         atoms of the current PDB file

***********************************************************************/

void apply_transformation(struct sequence *sequence_ptr,double matrix[4][4],
                          double *fit_cog)
{
  double x, y, z;

  struct atom *atom_ptr;

  /* Initialise pointer to first atom of this structure */
  atom_ptr = sequence_ptr->first_atom_ptr;

  /* Loop through all the atoms, writing each one out */
  while (atom_ptr != NULL)
    {
      /* Get this atom's coordinates */
      x = atom_ptr->x - matrix[3][0];
      y = atom_ptr->y - matrix[3][1];
      z = atom_ptr->z - matrix[3][2];

      /* Apply the transformation */
      atom_ptr->x = x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2];
      atom_ptr->y = x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2];
      atom_ptr->z = x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2];

      /* Shift back to the reference frame of the reference structure */
      atom_ptr->x = atom_ptr->x + fit_cog[0];
      atom_ptr->y = atom_ptr->y + fit_cog[1];
      atom_ptr->z = atom_ptr->z + fit_cog[2];

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

get_residue_list  -  Read in the list of residues from the appropriate
                     .res file

***********************************************************************/

void get_residue_list(struct sequence *sequence_ptr,
                      char reslist_name[FILENAME_LEN],
                      int *nresid,struct parameters params)
{
  char input_line[LINELEN + 1];
  char chain, res_name[4], res_num[6];

  int is_ligand, wanted;

  struct list *list_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  *nresid = 0;

  /* Open the ligplot.res file */
  if ((file_ptr = fopen(reslist_name,"r")) ==  NULL)
    {
      /* If generating a RasMol script, then don't display error message */
      if (params.rasmol_script == TRUE)
        return;
      else
        {
          printf("*** ERROR. File not found: %s\n",reslist_name);
          exit(1);
        }
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Determine whether this residue is wanted */
      is_ligand = FALSE;
      wanted = TRUE;

      /* Check whether this is a ligand residue */
      if (input_line[10] == '*')
        is_ligand = TRUE;

      /* Check for comment line */
      if (input_line[0] == '#')
        wanted = FALSE;

      /* If only the ligand residues are required, check if this is one
         of those */
      else if (sequence_ptr->display_option == LIGAND_ONLY &&
               is_ligand == FALSE)
        wanted = FALSE;

      /* If this residue is wanted then store details */
      if (wanted == TRUE)
        {
          /* Get the residue name, number and chain */
          strncpy(res_name,input_line,3);
          res_name[3] = '\0';
          strncpy(res_num,input_line+3,5);
          res_num[5] = '\0';
          chain = input_line[8];

          /* Create a new residue record and attach it to this
             sequence's list of residues */
          create_list_record(&list_ptr,sequence_ptr,res_name,
                             res_num,chain);

          /* Store whether this is a ligand residue */
          list_ptr->is_ligand = is_ligand;
          if (is_ligand == FALSE)
            list_ptr->interacting_residue = TRUE;

          /* Increment count of residues read in */
          (*nresid)++;
        }
    }

  /* Close the file */
  fclose(file_ptr);
}
/***********************************************************************

get_helices  -  Pick up the ids of two interacting helices from the
                promotif.dat file

***********************************************************************/

void get_helices(char chain,int motif_no,char *file_name,char *chains,
                 int *helix_no,FILE *outfile_ptr)
{
  char input_line[LINELEN + 1], number_string[20];

  int nmotif;
  int done, started, wanted;

  struct range *first_range_ptr, *last_range_ptr;
   
  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  nmotif = 0;
  started = FALSE;
  wanted = FALSE;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      fprintf(outfile_ptr,"echo \"*** ERROR. Unable to open file ");
      fprintf(outfile_ptr,"[%s] \"\n",file_name);
      fprintf(outfile_ptr,"echo \" \"\n");
      return;
    }

  /* Loop through the file until we pick up the motif we're after */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Check whether this is a START record */
      if (!strncmp(input_line,"START:",6))
        {
          /* If have already started, then this mark the end of our search */
          if (started == TRUE)
            done = TRUE;

          /* Otherwise, check whether this is the right motif */
          else if (!strcmp(input_line+6,"HELIX_INTERACTIONS"))
            {
              /* Set flag that we have found the right motif section */
              started = TRUE;
            }
        }

      /* Otherwise, if this is a motif record and we have started,
         see if this is the one we're after */
      else if (started == TRUE && !strncmp(input_line+1,"::",2))
        {
          /* Initialise flag */
          wanted = FALSE;

          /* If this is the right chain-id, increment motif count */
          if (input_line[0] == chain)
            {
              /* Increment motif count */
              nmotif++;

              /* Determine whether this record is wanted */
              if (nmotif == motif_no)
                wanted = TRUE;
            }

          /* If this is the correct motif, extract the helix ids
             and set flag */
          if (wanted == TRUE)
            {
              /* Get the two chain identifiers */
              chains[0] = input_line[3];
              chains[1] = input_line[14];

              /* Get the two helix numbers */
              strncpy(number_string,input_line+6,3);
              number_string[3] = '\0';
              helix_no[0] = atoi(number_string);
              strncpy(number_string,input_line+17,3);
              number_string[3] = '\0';
              helix_no[1] = atoi(number_string);

              /* Set flag that motif found */
              done = TRUE;
            }
        }
    }

  /* Close input file */
  fclose(file_ptr);
}
/***********************************************************************

create_range_record  -  Create a range record

***********************************************************************/

struct range *create_range_record(struct range **fst_range_ptr,
                                  struct range **lst_range_ptr,
                                  char chain,char *res_num1,
                                  char *res_name1,char *res_num2,
                                  char *res_name2,char *colour,
                                  char *label,int sidechain,
                                  float wireframe,int hbonds,
                                  char *range_name,char *name)
{
  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Create range entry */
  range_ptr = (struct range*)malloc(sizeof(struct range));
  if (range_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct range\n");
      printf("***        Program romlas terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the data fields */
  range_ptr->chain = chain;
  strcpy(range_ptr->res_name[0],res_name1);
  strcpy(range_ptr->res_name[1],res_name2);
  strcpy(range_ptr->res_num[0],res_num1);
  strcpy(range_ptr->res_num[1],res_num2);
  strcpy(range_ptr->colour,colour);
  strcpy(range_ptr->label,label);
  strcpy(range_ptr->range_name,range_name);
  strcpy(range_ptr->name,name);
  range_ptr->sidechain = sidechain;
  range_ptr->wireframe = wireframe;
  range_ptr->hbonds = hbonds;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_range_ptr == NULL)
    first_range_ptr = range_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_range_ptr->next_range_ptr = range_ptr;
                      
  /* Save the current residue's pointer */
  last_range_ptr = range_ptr;

  /* Set pointer to next entry to null */
  range_ptr->next_range_ptr = NULL;

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
  return(range_ptr);
}
/***********************************************************************

get_beta_turn  -  Get the beta turn residues

***********************************************************************/

void get_beta_turn(char *input_line,struct range **fst_range_ptr,
                   struct range **lst_range_ptr)
{
  char chain, res_num[4][6], res_name[4][4];
  char label[10], name[NAME_LEN + 1], range_name[NAME_LEN + 1];

  static char *colour[] = { "red", "orange", "green", "blue" };

  int ipos, ires;
  int sidechain;

  float wireframe;

  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Initialise variables */
  strcpy(range_name,"b_turn");
  strcpy(name,"beta turn");
  label[0] = '\0';
  sidechain = TRUE;
  wireframe = (float) 0.3;

  /* Get the chain-id */
  chain = input_line[0];
  ipos = 3;

  /* Loop over the 4 residues making up the turn to get names and numbers */
  for (ires = 0; ires < 4; ires++)
    {
      /* Get the info on this residue */
      strncpy(res_num[ires],input_line+ipos,5);
      res_num[ires][5] = '\0';
      strncpy(res_name[ires],input_line+ipos+7,3);
      res_name[ires][3] = '\0';

      /* Set position for next residue */
      ipos = ipos + 12;
    }

  /* Create a range record covering the whole turn */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[0],res_name[0],res_num[3],
                                  res_name[3],"white",label,sidechain,
                                  wireframe,TRUE,range_name,name);

  /* Reset identifiers */
  range_name[0] = '\0';
  name[0] = '\0';

  /* Create a range record for each residue in the turn */
  for (ires = 0; ires < 4; ires++)
    {
      /* Create this range record */
      range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                      chain,res_num[ires],res_name[ires],
                                      res_num[ires],res_name[ires],
                                      colour[ires],label,sidechain,
                                      wireframe,FALSE,range_name,name);
    }

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
}
/***********************************************************************

get_gamma_turn  -  Get the gamma turn residues

***********************************************************************/

void get_gamma_turn(char *input_line,struct range **fst_range_ptr,
                    struct range **lst_range_ptr)
{
  char chain, res_num[3][6], res_name[3][4];
  char label[10], name[NAME_LEN + 1], range_name[NAME_LEN + 1];

  static char *colour[] = { "red", "green", "blue" };

  int ipos, ires;
  int sidechain;

  float wireframe;

  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Initialise variables */
  strcpy(range_name,"g_turn");
  strcpy(name,"gamma turn");
  label[0] = '\0';
  sidechain = TRUE;
  wireframe = (float) 0.3;

  /* Get the chain-id */
  chain = input_line[0];
  ipos = 3;

  /* Loop over the 3 residues making up the turn to get names and numbers */
  for (ires = 0; ires < 3; ires++)
    {
      /* Get the info on this residue */
      strncpy(res_num[ires],input_line+ipos,5);
      res_num[ires][5] = '\0';
      strncpy(res_name[ires],input_line+ipos+7,3);
      res_name[ires][3] = '\0';

      /* Set position for next residue */
      ipos = ipos + 12;
    }

  /* Create a range record covering the whole turn */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[0],res_name[0],res_num[2],
                                  res_name[2],"white",label,sidechain,
                                  wireframe,TRUE,range_name,name);

  /* Reset identifiers */
  range_name[0] = '\0';
  name[0] = '\0';

  /* Create a range record for each residue in the turn */
  for (ires = 0; ires < 3; ires++)
    {
      /* Create this range record */
      range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                      chain,res_num[ires],res_name[ires],
                                      res_num[ires],res_name[ires],
                                      colour[ires],label,sidechain,
                                      wireframe,FALSE,range_name,name);
    }

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
}
/***********************************************************************

get_bulge  -  Get the beta bulge residues

***********************************************************************/

void get_bulge(char *input_line,struct range **fst_range_ptr,
               struct range **lst_range_ptr)
{
  char chain, res_num[5][6], res_name[5][4];
  char label[10], name[NAME_LEN + 1], range_name[NAME_LEN + 1];

  char number_string[20];

  static char bulge_pos[] = "X1234";
  static char *colour[] = { "red", "green", "blue", "purple", "purple" };

  int ipos, ires, strand;
  int sidechain;

  float wireframe;

  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Initialise variables */
  label[0] = '\0';
  sidechain = TRUE;
  wireframe = (float) 0.3;

  /* Get start and end residues used in beta-bulge plot */
  ipos = 199;
  for (strand = 0; strand < 2; strand++)
    {
      /* Get this strand's chain identifier */
      chain = input_line[ipos];
      ipos = ipos + 3;

      /* Loop over the start and end residues */
      for (ires = 0; ires < 2; ires++)
        {
          /* Get the residue number */
          strncpy(res_num[ires],input_line+ipos,5);
          res_num[ires][5] = '\0';

          /* Don't have the residue name */
          res_name[ires][0] = '\0';

          /* Set position for next residue */
          ipos = ipos + 7;
        }

      /* Form range name */
      strcpy(range_name,"strand");
      strcpy(name,"strand ");
      sprintf(number_string,"%d",strand + 1);
      strcat(range_name,number_string);
      strcat(name,number_string);

      /* Create a range record covering this strand */
      range_ptr
        = create_range_record(&first_range_ptr,&last_range_ptr,
                              chain,res_num[0],res_name[0],res_num[1],
                              res_name[1],"white",label,sidechain,
                              wireframe,TRUE,range_name,name);
    }

  /* Initialise start position for bulge residues */
  ipos = 3;

  /* Loop over the 5 residues making up the bulge to get names and numbers */
  for (ires = 0; ires < 5; ires++)
    {
      /* Get the chain-id */
      chain = input_line[ipos];
      ipos = ipos + 3;

      /* Get the info on this residue */
      strncpy(res_num[ires],input_line+ipos,5);
      res_num[ires][5] = '\0';
      strncpy(res_name[ires],input_line+ipos+7,3);
      res_name[ires][3] = '\0';

      /* Set position for next residue */
      ipos = ipos + 12;

      /* If have a valid bulge residue, write out */
      if (strcmp(res_name[ires],"   "))
        {
          /* Form range name */
          strcpy(range_name,"res");
          strcpy(name,"bulge residue ");
          sprintf(number_string,"%c",bulge_pos[ires]);
          strcat(range_name,number_string);
          strcat(name,number_string);

          /* Create this range record */
          range_ptr
            = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[ires],res_name[ires],
                                  res_num[ires],res_name[ires],
                                  colour[ires],label,sidechain,
                                  wireframe,FALSE,range_name,name);
        }
    }

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
}
/***********************************************************************

get_hairpin  -  Get the beta hairpin residues

***********************************************************************/

void get_hairpin(char *input_line,struct range **fst_range_ptr,
                 struct range **lst_range_ptr,int motif_type)
{
  char chain, res_num[2][2][6], res_name[2][2][4];
  char loop_num[2][6], loop_name[2][4];
  char label[10], name[NAME_LEN + 1], range_name[NAME_LEN + 1];

  int ipos, ires, strand;
  int sidechain;

  float wireframe;

  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Initialise variables */
  if (motif_type == BETA_HAIRPINS)
    {
      strcpy(range_name,"hairpin");
      strcpy(name,"beta hairpin");
    }
  else if (motif_type == BAB_MOTIFS)
    {
      strcpy(range_name,"bab_motif");
      strcpy(name,"beta-alpha-beta motif");
    }
  else if (motif_type == PSI_LOOPS)
    {
      strcpy(range_name,"psi_loop");
      strcpy(name,"psi loop");
    }
  label[0] = '\0';
  sidechain = TRUE;
  wireframe = (float) 0.3;

  /* Get the chain-id */
  chain = input_line[0];
  ipos = 3;

  /* Loop over the start- and end ranges of the two strands */
  for (strand = 0; strand < 2; strand++)
    {
      /* Loop over this strand's start and end residues */
      for (ires = 0; ires < 2; ires++)
        {
          /* Get the info on this residue */
          strncpy(res_num[strand][ires],input_line+ipos,5);
          res_num[strand][ires][5] = '\0';
          strncpy(res_name[strand][ires],input_line+ipos+7,3);
          res_name[strand][ires][3] = '\0';

          /* Set position for next residue */
          ipos = ipos + 12;
        }
    }
     
  /* Create a range record covering the whole hairpin */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[0][0],res_name[0][0],
                                  res_num[1][1],res_name[1][1],"white",
                                  label,sidechain,wireframe,TRUE,
                                  range_name,name);

  /* Reset identifiers */
  range_name[0] = '\0';
  name[0] = '\0';

  /* Create a range covering strand 1 */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[0][0],res_name[0][0],
                                  res_num[0][1],res_name[0][1],
                                  "blue",label,sidechain,
                                  wireframe,FALSE,range_name,name);

  /* Create a range covering strand 2 */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[1][0],res_name[1][0],
                                  res_num[1][1],res_name[1][1],
                                  "green",label,sidechain,
                                  wireframe,FALSE,range_name,name);

  /* Determine the loop's start residue number */
  strcpy(loop_num[0],res_num[0][1]);
  loop_num[0][4] = '\0';
  ires = atoi(loop_num[0]);
  ires++;
  sprintf(loop_num[0],"%4d ",ires);

  /* Repeat for its end residue number */
  strcpy(loop_num[1],res_num[1][0]);
  loop_num[1][4] = '\0';
  ires = atoi(loop_num[1]);
  ires--;
  sprintf(loop_num[1],"%4d ",ires);

  /* Set residue names as blank (as not currently needed) */
  loop_name[0][0] = '\0';
  loop_name[1][0] = '\0';

  /* Create a range covering the loop */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,loop_num[0],loop_name[0],
                                  loop_num[1],loop_name[1],
                                  "red",label,sidechain,
                                  wireframe,FALSE,range_name,name);

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
}
/***********************************************************************

format_resname  -  Join residue name and number

***********************************************************************/

void format_resname(char *res_name,char *res_num,char *out_name)
{
  char tmp_num[6];

  int len;

  /* Initialise variables */
  strcpy(out_name,res_name);

  /* Convert last 2 chars of residue name to lower case */
  to_lower(out_name+1,2);
  out_name[3] = '\0';

  /* Remove leading spaces from residue number */
  strcpy(tmp_num,res_num);
  remove_leading_blanks(tmp_num);

  /* Add on to residue name */
  strcat(out_name,tmp_num);

  /* Truncate the string */
  len = strlen(out_name);
  string_truncate(out_name,len+1);
}
/***********************************************************************

get_helix  -  Get the helix residues

***********************************************************************/

void get_helix(char *input_line,struct range **fst_range_ptr,
               struct range **lst_range_ptr,int helint)
{
  char chain, res_num[2][6], res_name[2][4], out_name[9];
  char label[10], name[NAME_LEN + 1], number_string[20];
  char range_name[NAME_LEN + 1];

  static char *colour[] = { "white", "red", "blue" };

  int helix, ipos, ires;
  int sidechain;

  float wireframe;

  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Initialise variables */
  strcpy(range_name,"helix");
  strcpy(name,"helix number ");
  label[0] = '\0';
  sidechain = TRUE;
  wireframe = (float) 0.3;

  /* Get the helix number */
  strncpy(number_string,input_line+3,3);
  number_string[3] = '\0';
  helix = atoi(number_string);
  sprintf(number_string,"%d",helix);
  strcat(range_name,number_string);
  strcat(name,number_string);

  /* Get the chain-id */
  chain = input_line[0];
  ipos = 8;

  /* Loop over the helix's start- and end residues */
  for (ires = 0; ires < 2; ires++)
    {
      /* Get the info on this residue */
      strncpy(res_num[ires],input_line+ipos,5);
      res_num[ires][5] = '\0';
      strncpy(res_name[ires],input_line+ipos+7,3);
      res_name[ires][3] = '\0';

      /* Set position for next residue */
      ipos = ipos + 12;
    }

  /* Add helix residue range to name */
  strcat(name," (");
  format_resname(res_name[0],res_num[0],out_name);
  strcat(name,out_name);
  strcat(name,"-");
  format_resname(res_name[1],res_num[1],out_name);
  strcat(name,out_name);
  strcat(name,")");
     
  /* Create a range record covering the whole helix */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain,res_num[0],res_name[0],
                                  res_num[1],res_name[1],colour[helint],
                                  label,sidechain,wireframe,TRUE,
                                  range_name,name);

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
}
/***********************************************************************

get_strand  -  Get the start and end strand residues

***********************************************************************/

void get_strand(char *input_line,struct range **fst_range_ptr,
                struct range **lst_range_ptr)
{
  char chain[2], res_num[2][6], res_name[2][4], out_name[9];
  char label[10], name[NAME_LEN + 1], number_string[20];
  char range_name[NAME_LEN + 1];

  int strand, ipos, ires;
  int sidechain;

  float wireframe;

  struct range *range_ptr;
  struct range *first_range_ptr, *last_range_ptr;

  /* Initialise range pointers */
  range_ptr = NULL;
  first_range_ptr = *fst_range_ptr;
  last_range_ptr = *lst_range_ptr;

  /* Initialise variables */
  strcpy(range_name,"strand");
  strcpy(name,"strand number ");
  label[0] = '\0';
  sidechain = TRUE;
  wireframe = (float) 0.3;

  /* Get the strand number */
  strncpy(number_string,input_line+3,3);
  number_string[3] = '\0';
  strand = atoi(number_string);
  sprintf(number_string,"%d",strand);
  strcat(range_name,number_string);
  strcat(name,number_string);

  /* Get the chain-id */
  chain[0] = input_line[0];
  chain[1] = '\0';
  ipos = 8;

  /* Loop over the strand's start- and end residues */
  for (ires = 0; ires < 2; ires++)
    {
      /* Get the info on this residue */
      strncpy(res_num[ires],input_line+ipos,5);
      res_num[ires][5] = '\0';
      strncpy(res_name[ires],input_line+ipos+7,3);
      res_name[ires][3] = '\0';

      /* Set position for next residue */
      ipos = ipos + 12;
    }

  /* Add strand residue range to name */
  strcat(name," (");
  format_resname(res_name[0],res_num[0],out_name);
  strcat(name,out_name);
  if (chain[0] != ' ')
    {
      strcat(name,chain);
    }
  strcat(name,"-");
  format_resname(res_name[1],res_num[1],out_name);
  strcat(name,out_name);
  if (chain[0] != ' ')
    {
      strcat(name,chain);
    }
  strcat(name,")");
     
  /* Create a range record covering the whole strand */
  range_ptr = create_range_record(&first_range_ptr,&last_range_ptr,
                                  chain[0],res_num[0],res_name[0],
                                  res_num[1],res_name[1],"red",
                                  label,sidechain,wireframe,FALSE,
                                  range_name,name);

  /* Return current pointers */
  *fst_range_ptr = first_range_ptr;
  *lst_range_ptr = last_range_ptr;
}
/***********************************************************************

get_chains  -  Get which protein chains are to be read in and output

***********************************************************************/

void get_chains(char chain,struct range *first_range_ptr,
                struct parameters *params)
{
  int ichain, nchains;
  int found;

  struct range *range_ptr;

  /* Initialise variables */
  nchains = 0;

  /* Store first chain as definitely wanted */
  params->chains[nchains] = chain;
  nchains++;

  /* Initialise pointer to the first record */
  range_ptr = first_range_ptr;

  /* Loop through all the range records to get the chains */
  while (range_ptr != NULL)
    {
      /* If this is not the default chains, see whether it is a new
         one that needs to be stored */
      if (range_ptr->chain != chain)
        {
          /* Initialise flag */
          found = FALSE;

          /* Loop through the stored chains to see if we already have
             this one */
          for (ichain = 1; ichain < nchains; ichain++)
            if (params->chains[ichain] == range_ptr->chain)
              found = TRUE;

          /* If chain not found, add to list */
          if (found == FALSE && nchains < MAX_CHAINS)
            {
              params->chains[nchains] = range_ptr->chain;
              nchains++;
            }
        }

      /* Get pointer to the next range record */
      range_ptr = range_ptr->next_range_ptr;
    }

  /* Store number of chains */
  params->nchains = nchains;
}
/***********************************************************************

get_motif  -  Pick up the required motif from the promotif.dat file

***********************************************************************/

int get_motif(struct sequence *sequence_ptr,char *file_name,
              struct range **fst_range_ptr,FILE *outfile_ptr,
              struct parameters *params)
{
  char chain, chains[2], input_line[LINELEN + 1], number_string[20];
  char search_motif[30];

  static char *motif_name[NMOTIFS] = {
    "BETA_TURNS", "GAMMA_TURNS", "BETA_BULGES", "BETA_HAIRPINS",
    "HELICES", "BAB_MOTIFS", "PSI_LOOPS", "STRANDS",
    "HELICES"
  };

  int helint, helix_id, helix_no[2], i, nmotif, strand_id;
  int motif_no;
  int done, found, started, wanted;

  struct range *first_range_ptr, *last_range_ptr;
   
  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  found = FALSE;
  helint = 0;
  motif_no = params->motif_no;
  nmotif = 0;
  started = FALSE;
  wanted = FALSE;
  first_range_ptr = NULL;
  last_range_ptr = NULL;
  strcpy(search_motif,motif_name[params->motif_type]);
  strand_id = 0;
  for (i = 0; i < 2; i++)
    {
      chains[i] = '\0';
      helix_no[i] = 0;
    }

  /* If looking for a specific helix interaction, need to identify
     which helices we're after */
  if (params->motif_type == HELIX_INTERACTIONS)
    get_helices(sequence_ptr->chain,motif_no,file_name,chains,
                helix_no,outfile_ptr);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      fprintf(outfile_ptr,"echo \"*** ERROR. Unable to open file ");
      fprintf(outfile_ptr,"[%s] \"\n",file_name);
      fprintf(outfile_ptr,"echo \" \"\n");
      return(FALSE);
    }

  /* Loop through the file until we pick up the motif we're after */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Check whether this is a START record */
      if (!strncmp(input_line,"START:",6))
        {
          /* If have already started, then this mark the end of our search */
          if (started == TRUE)
            done = TRUE;

          /* Otherwise, check whether this is the right motif */
          else if (!strcmp(input_line+6,search_motif))
            {
              /* Set flag that we have found the right motif section */
              started = TRUE;
            }
        }

      /* Otherwise, if this is a motif record and we have started,
         see if this is the one we're after */
      else if (started == TRUE && !strncmp(input_line+1,"::",2))
        {
          /* Initialise flag */
          wanted = FALSE;

          /* If reading the strands that make up the given sheet, get
             this strand's sheet id */
          if (params->motif_type == SHEETS &&
              !strcmp(search_motif,"STRANDS"))
            {
              /* Extract the strand id */
              strncpy(number_string,input_line+32,3);
              number_string[3] = '\0';
              strand_id = atoi(number_string);

              /* If this is the right sheet, then want this strand */
              if (strand_id == motif_no)
                wanted = TRUE;
            }

          /* If this is a helix interaction, check whether this is
             one of the required helices */
          else if (params->motif_type == HELIX_INTERACTIONS)
            {
              /* Extract the chain and helix id */
              chain = input_line[0];
              strncpy(number_string,input_line+3,3);
              number_string[3] = '\0';
              helix_id = atoi(number_string);

              /* Check whether it is either helix */
              for (i = 0; i < 2; i++)
                {
                  if (chain == chains[i] && helix_id == helix_no[i])
                    {
                      helint = i + 1;
                      wanted = TRUE;
                    }
                }
            }

          /* If this is the right chain-id, increment motif count */
          else if (input_line[0] == sequence_ptr->chain)
            {
              /* Increment motif count */
              nmotif++;

              /* Determine whether this record is wanted */
              if (nmotif == motif_no)
                {
                  /* Set flag */
                  wanted = TRUE;

                  /* If this is a sheet record, extract the sheet
                     identifer */
                }
            }

          /* If this is the correct motif, extract the relevant details
             and set flag */
          if (wanted == TRUE)
            {
              /* Extract motif details */
              if (params->motif_type == BETA_TURNS)
                get_beta_turn(input_line,&first_range_ptr,
                              &last_range_ptr);
              else if (params->motif_type == GAMMA_TURNS)
                get_gamma_turn(input_line,&first_range_ptr,
                               &last_range_ptr);
              else if (params->motif_type == BETA_BULGES)
                get_bulge(input_line,&first_range_ptr,
                          &last_range_ptr);
              else if (params->motif_type == BETA_HAIRPINS ||
                       params->motif_type == BAB_MOTIFS ||
                       params->motif_type == PSI_LOOPS)
                get_hairpin(input_line,&first_range_ptr,
                            &last_range_ptr,params->motif_type);
              else if (params->motif_type == HELICES ||
                       params->motif_type == HELIX_INTERACTIONS)
                get_helix(input_line,&first_range_ptr,&last_range_ptr,
                          helint);
              else if (params->motif_type == SHEETS)
                get_strand(input_line,&first_range_ptr,&last_range_ptr);

              /* Set flag that motif found */
              found = TRUE;

              /* Unless need to find more motif records, set flags that
                 we're done */
              if (params->motif_type != SHEETS &&
                  params->motif_type != HELIX_INTERACTIONS)
                done = TRUE;
            }
        }
    }

  /* Close input file */
  fclose(file_ptr);

  /* Get all the chains that need to be read in from the PDB file */
  get_chains(sequence_ptr->chain,first_range_ptr,params);

  /* Return this motif's first range */
  *fst_range_ptr = first_range_ptr;

  /* Return whether required motif was found */
  return(found);
}
/***********************************************************************

get_site  -  Read through the PDB file to pick up the residues belonging
             to the given SITE

***********************************************************************/

void get_site(struct sequence *sequence_ptr,char *pdb_name,
              FILE *outfile_ptr,struct parameters *params)
{
  char input_line[LINELEN + 1];
  char chain, res_num[6], res_name[4], site_name[4];

  int ichain, ipos, ires, len, nchains;
  int done, eof, found, gzipped, keep_reading;

  struct list *list_ptr;
    
#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  eof = FALSE;
  keep_reading = TRUE;
  nchains = 0;

  /* Initialise pointers */
  list_ptr = NULL;

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
    }
  else
    eof = TRUE;

  /* If file found, then read in the data */
  if (eof == TRUE)
    {
      fprintf(outfile_ptr,"echo \"*** ERROR. Unable to open PDB file ");
      fprintf(outfile_ptr,"[%s] \"\n",pdb_name);
      fprintf(outfile_ptr,"echo \" \"\n");
      return;
    }

  /* Initialise flag */
  keep_reading = TRUE;

  /* Loop through the PDB file, picking up all the atoms that are
     required */
  while (eof == FALSE && keep_reading == TRUE)
    {
      /* Read in the next record */
      if (gzipped == TRUE)
        {
#ifdef TEST
#else
          if (gzgets(gzfile_ptr,input_line,R_LINELEN) == NULL)
            eof = TRUE;
#endif
        }
      else
        {
          if (fgets(input_line,R_LINELEN,file_ptr) == NULL)
            eof = TRUE;
        }

      /* If record read in OK, then process it */
      if (eof == FALSE)
        {
          /* Check for the start of the ATOM records */
          if (!strncmp("ATOM  ",input_line,6) ||
              !strncmp("HETATM",input_line,6))
            keep_reading = FALSE;

          /* If this is a SITE record, determine whether this is the one */
          else if (!strncmp(input_line,"SITE  ",6))
            {
              /* Truncate the string to strip off the newline */
              string_truncate(input_line,LINELEN);

              /* Add blank to end of line to ensure that truncated
                 lines are correctly interpreted */
              strcat(input_line,"                                      ");

              /* Get the site name */
              strncpy(site_name,input_line+11,3);
              site_name[3] = '\0';

              /* If this is the correct site, get the residues from this
                 line */
              if (!strcmp(params->site_id,site_name))
                {
                  /* Initialise the record position for reading in
                     the residues in this record */
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
                      if (!strcmp(res_name,"   ") ||
                          !strcmp(res_num,"     "))
                        done = TRUE;

                      /* Otherwise, create a reslist record to store
                         residue */
                      else
                        {
                          /* Create the record */
                          create_list_record(&list_ptr,sequence_ptr,res_name,
                                             res_num,chain);

                          /* Set colour */
                          list_ptr->colour = 1;

                          /* If this is the first record, store the chain
                             identifier */
                          if (nchains == 0)
                            {
                              params->chains[nchains] = chain;
                              nchains++;
                            }

                          /* Otherwise, if this is a different chain,
                             add to list */
                          else if (chain != params->chains[0])
                            {
                              /* Initialise flag */
                              found = FALSE;

                              /* Loop through the stored chains to see if
                                 we already have this one */
                              for (ichain = 1; ichain < nchains; ichain++)
                                if (params->chains[ichain] == chain)
                                  found = TRUE;

                              /* If chain not found, add to list */
                              if (found == FALSE && nchains < MAX_CHAINS)
                                {
                                  params->chains[nchains] = chain;
                                  nchains++;
                                }
                            }
                        }

                      /* Increment position for next residue */
                      ipos = ipos + 11;
                    }
                }
            }
        }
    }

  /* Close the PDB file */
  if (gzipped == TRUE)
    {
#ifdef TEST
#else
      gzclose(gzfile_ptr);
#endif
    }
  else
    fclose(file_ptr);

  /* Store number of chains */
  params->nchains = nchains;
}
/***********************************************************************

create_hbond_entry  -  Create a hbond entry

***********************************************************************/

struct hbond *create_hbond_entry(struct hbond **fst_hbond_ptr,
                                 struct hbond **lst_hbond_ptr,
                                 struct atom *protein_atom_ptr,
                                 struct atom *ligand_atom_ptr,
                                 double length,int donor_first)
{
  struct hbond *hbond_ptr;
  struct hbond *first_hbond_ptr, *last_hbond_ptr;


  /* Initialise hbond pointers */
  hbond_ptr = NULL;
  first_hbond_ptr = *fst_hbond_ptr;
  last_hbond_ptr = *lst_hbond_ptr;

  /* Create hbond entry */
  hbond_ptr = (struct hbond*)malloc(sizeof(struct hbond));
  if (hbond_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct hbond\n");
      printf("***        Program romlas terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the data for this entry */
  if (donor_first == TRUE)
    {
      hbond_ptr->donor_atom_ptr = protein_atom_ptr;
      hbond_ptr->acceptor_atom_ptr = ligand_atom_ptr;
    }
  else
    {
      hbond_ptr->acceptor_atom_ptr = protein_atom_ptr;
      hbond_ptr->donor_atom_ptr = ligand_atom_ptr;
    }
  hbond_ptr->length = length;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_hbond_ptr == NULL)
    first_hbond_ptr = hbond_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_hbond_ptr->next_hbond_ptr = hbond_ptr;
                      
  /* Save the current residue's pointer */
  last_hbond_ptr = hbond_ptr;

  /* Set pointer to next entry to null */
  hbond_ptr->next_hbond_ptr = NULL;

  /* Return current pointers */
  *fst_hbond_ptr = first_hbond_ptr;
  *lst_hbond_ptr = last_hbond_ptr;
  return(hbond_ptr);
}
/***********************************************************************

get_hbonds  -  Read in the H-bond information from the appropriate
               grow.out file

***********************************************************************/

int get_hbonds(char file_name[FILENAME_LEN],struct sequence *sequence_ptr)
{
  char chain[2], atom_name[2][5], res_name[2][4], res_num[2][6];
  char input_line[LINELEN + 1], number_string[20];

  int line, nhbonds;
  int donor_first, have_file;

  double length;

  struct atom *atom_ptr, *ligand_atom_ptr, *protein_atom_ptr;
  struct hbond *first_hbond_ptr, *last_hbond_ptr, *hbond_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nhbonds = 0;

  /* Initialise hbond pointers */
  hbond_ptr = NULL;
  first_hbond_ptr = sequence_ptr->first_hbond_ptr;
  last_hbond_ptr = NULL;

  /* Initialise flag indicating whether GROW file present */
  have_file = TRUE;
  line = 0;

  /* Open the GROW file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If have the file, read in all its data */
  if (have_file == TRUE)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Only interested in protein-ligand or DNA-ligand H-bonds */
          if (!strncmp(input_line,"H P L",5) ||
              !strncmp(input_line,"H N L",5))
            {
              /* Extract the protein atom and residue details */
              strncpy(atom_name[0],input_line+24,4);
              atom_name[0][4] = '\0';
              strncpy(res_name[0],input_line+14,3);
              res_name[0][3] = '\0';
              strncpy(res_num[0],input_line+9,5);
              res_num[0][5] = '\0';
              chain[0] = input_line[6];

              /* Extract the protein atom and residue details */
              strncpy(atom_name[1],input_line+33,4);
              atom_name[1][4] = '\0';
              strncpy(res_name[1],input_line+44,3);
              res_name[1][3] = '\0';
              strncpy(res_num[1],input_line+49,5);
              res_num[1][5] = '\0';
              chain[1] = input_line[54];

              /* Initialise pointer to first atom */
              atom_ptr = sequence_ptr->first_atom_ptr;
              ligand_atom_ptr = NULL;
              protein_atom_ptr = NULL;

              /* Loop through all the atoms to check if it is one of the
                 atoms involved in the current H-bond */
              while (atom_ptr != NULL && (ligand_atom_ptr == NULL ||
                                          protein_atom_ptr == NULL))
                {
                  /* Get this atom's residue pointer */
                  residue_ptr = atom_ptr->residue_ptr;

                  /* Check whether it is the protein atom */
                  if (protein_atom_ptr == NULL &&
                      !strcmp(atom_ptr->atom_name,atom_name[0]) &&
                      !strcmp(residue_ptr->res_name,res_name[0]) &&
                      !strcmp(residue_ptr->res_num,res_num[0]) &&
                      (residue_ptr->chain == chain[0] ||
                       (residue_ptr->chain == ' ' && chain[0] == 'A')))
                    {
                      /* Have a match, so store pointer */
                      protein_atom_ptr = atom_ptr;
                    }

                  /* Check whether it is the ligand atom */
                  if (ligand_atom_ptr == NULL &&
                      !strcmp(atom_ptr->atom_name,atom_name[1]) &&
                      !strcmp(residue_ptr->res_name,res_name[1]) &&
                      !strcmp(residue_ptr->res_num,res_num[1]) &&
                      (residue_ptr->chain == chain[1] ||
                       (residue_ptr->chain == ' ' && chain[1] == 'A')))
                    {
                      /* Have a match, so store pointer */
                      ligand_atom_ptr = atom_ptr;
                    }

                  /* Get pointer to next atom in linked-list */
                  atom_ptr = atom_ptr->next_atom_ptr;
                }

              /* If have located both atoms, then can store this bond */
              if (ligand_atom_ptr != NULL && protein_atom_ptr != NULL)
                {
                  /* Get the bond length */
                  strncpy(number_string,input_line+55,6);
                  number_string[6] = '\0';
                  length = atof(number_string);

                  /* Determine which are the donor and acceptor atoms */
                  if (!strncmp(input_line+29,"-->",3))
                    donor_first = TRUE;
                  else
                    donor_first = FALSE;

                  /* Create an H-bond record for this bond */
                  hbond_ptr
                    = create_hbond_entry(&first_hbond_ptr,
                                         &last_hbond_ptr,protein_atom_ptr,
                                         ligand_atom_ptr,length,
                                         donor_first);

                  /* Increment count of H-bonds stored */
                  nhbonds++;
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Store current pointer */
  sequence_ptr->first_hbond_ptr = first_hbond_ptr;

  /* Return number of H-bonds stored */
  return(nhbonds);
}
/***********************************************************************

process_conect_record  -  Check whether the given CONECT record applies
                          to one of the ATOM or HETATM records already
                          written out and hence is required

***********************************************************************/

void process_conect_record(struct sequence *sequence_ptr,
                           char input_line[LINELEN + 1])
{
  char atmnum[10];

  int atom_number, iconect, ipos, nconect;
  int done, wanted;

  struct atom *atom_ptr, *save_atom_ptr, *store_atom_ptr[MAXCONECT];
  struct conect *conect_ptr;

  /* Initialise variables */
  done = FALSE;
  ipos = 6;
  nconect = 0;

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

          /* Check whether this atom number is one of those written
             out */
          wanted = FALSE;

          /* Initialise pointer to the first record */
          atom_ptr = sequence_ptr->first_atom_ptr;

          /* Loop through all the atom records */
          while (atom_ptr != NULL && wanted == FALSE)
            {
              /* Consider atom only if it is one of the ones to be
                 displayed */
              if (atom_ptr->display_atom == TRUE)
                {
                  /* If have a match, then set flag */
                  if (atom_number == atom_ptr->atom_number)
                    {
                      /* Save the atom pointer */
                      wanted = TRUE;
                      save_atom_ptr = atom_ptr;
                    }
                }

              /* Get pointer to the next atom record */
              atom_ptr = atom_ptr->next_atom_ptr;
            }

          /* If record is wanted then store */
          if (wanted == TRUE)
            {
              /* Get the saved atom pointer */
              atom_ptr = save_atom_ptr;

              /* Add atom number to the list of stored numbers */
              store_atom_ptr[nconect] = atom_ptr;
              nconect++;
            }

          /* If not wanted then check whether it is the first atom, in
             which case we don't want this CONECT record at all */
          else if (iconect == 0)
            done = TRUE;
        }

      /* Shift position along input line */
      ipos = ipos + 5;
    }

  /* If have stored at least one pair of connect-numbers, then create
     a CONECT record in which to store them */
  if (nconect > 1)
    {
      conect_ptr = (struct conect *) malloc(sizeof(struct conect));
      if (conect_ptr == NULL)
        {
          printf("*** Can't allocate memory for struct conect\n");
          exit (1);
        }

      /* Add link from previous CONECT record to the current one */
      if (last_conect_ptr != NULL)
        last_conect_ptr->next_conect_ptr = conect_ptr;
      else
        first_conect_ptr = conect_ptr;
      last_conect_ptr = conect_ptr;

      /* Initialise all the connectivity numbers */
      for (iconect = 0; iconect < MAXCONECT; iconect++)
        conect_ptr->atom_ptr[iconect] = NULL;

      /* Fill the connectivity numbers with the data stored */
      for (iconect = 0; iconect < nconect; iconect++)
        conect_ptr->atom_ptr[iconect] = store_atom_ptr[iconect];

      /* Initialise pointer to next record */
      conect_ptr->next_conect_ptr = NULL;

      /* If this is the sequence's first one, then store pointer */
      if (sequence_ptr->first_conect_ptr == NULL)
        sequence_ptr->first_conect_ptr = conect_ptr;

      /* Increment count of sequence's conect records */
      sequence_ptr->nconect++;
    }
}
/***********************************************************************

add_residue_range  -  Add current residue range to output line

***********************************************************************/

int add_residue_range(FILE *outfile_ptr,char object[15],
                      char output_line[LINELEN + 1],int from_atom,
                      int to_atom,int nranges)
{
  char atmnum[10];

  /* Check current length of string */
  output_line[0] = '\0';

  /* If the first range, then start the define */
  if (object[0] != '\0')
    {
      if (nranges == 0)
        sprintf(output_line,"define %s ",object);

      /* Otherwise, continue the define */
      else
        sprintf(output_line,"define %s %s | ",object,object);
    }

  /* First atom of range */
  strcat(output_line," (atomno >= ");

  /* Form atom number and insert in output line */
  sprintf(atmnum,"%d",from_atom);
  strcat(output_line,atmnum);

  /* Second atom of range */
  strcat(output_line," & atomno <= ");

  /* Form atom number and insert in output line */
  sprintf(atmnum,"%d",to_atom);
  strcat(output_line,atmnum);

  /* Close brackets */
  strcat(output_line,")");

  /* Increment count of ranges */
  nranges++;

  /* Write definition out to file */
  if (object[0] != '\0')
    fprintf(outfile_ptr,"%s\n",output_line);

  /* Return ranges count */
  return(nranges);
}
/***********************************************************************

write_script_residues  -  Write out the script commands to define the
                          ligand and non-ligand residues

***********************************************************************/

void write_script_residues(FILE *outfile_ptr,int ifile,
                           struct sequence *sequence_ptr,
                           int nligands,double rmsd,int nfit,
                           int in_ligand,char *oname,
                           struct parameters params)
{
  char chain_str[2], domain_str[5], object[15], object_type[10];
  char full_chain[5], label[10], residue_label[10], rmsd_out[10];
  char output_line[LINELEN + 1], number_string[3];
  char last_residue[4], last_res_number[6], last_chain;
  char ch, chain, pdb_code[5];
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];

  static char *col_name[] = {
    "purple", "yellow", "blue", "red", "green",
    "brown", "pink", "greenblue", "magenta",
    "cyan", "redorange", "violet", "orange",
    "white"
  };

  static char *col_nameb[] = {
    "purple", "blue", "red", "green",
    "brown", "pink", "greenblue", "magenta",
    "cyan", "redorange", "violet", "orange",
    "yellow", "white"
  };

  static char *col_numbers[] = {
    "[178,51,255]", "[255,255,0]", "[0,0,255]", "[255,0,0]", "[0,255,0]",
    "[128,0,0]", "[255,128,255]", "[51,153,102]", "[255,0,255]",
    "[128,255,255]", "[255,102,0]", "[238,130,238]", "[204,128,0]",
    "[255,255,255]"
  };

  static char *col_numbersb[] = {
    "[178,51,255]", "[0,0,255]", "[255,0,0]", "[0,255,0]",
    "[128,0,0]", "[255,128,255]", "[51,153,102]", "[255,0,255]",
    "[128,255,255]", "[255,102,0]", "[238,130,238]", "[204,128,0]",
    "[255,255,0]", "[255,255,255]"
  };

  int colour, domain_no, from_atom, ipos, jfile, len, nranges, to_atom;
  int label_atoms, label_residues, last_in_ligand, new_residue, type;
  int factor, i, lower_b, upper_b;
  int match;

  static int nlig = 0;
  static int nprot = 0;

  struct atom *atom_ptr;
  struct residue *residue_ptr;
  struct sequence *other_sequence_ptr;

  /* Initialise variables */
  chain = sequence_ptr->chain;
  chain_str[0] = sequence_ptr->chain;
  chain_str[1] = '\0';
  domain_no = sequence_ptr->domain_no;
  label[0] = '\0';
  label_atoms = FALSE;
  label_residues = FALSE;
  last_residue[0] = '\0';
  last_res_number[0] = '\0';
  last_chain = '\0';
  last_in_ligand = in_ligand;
  from_atom = 0;
  to_atom = 0;
  nranges = 0;
  residue_label[0] = '\0';
  strcpy(pdb_code,sequence_ptr->pdb_code);
  type = params.pdb_type;
  if (params.jmol == TRUE)
    factor = 1;
  else
    factor = 100;

  /* Initialise object name */
  if (in_ligand == TRUE)
    {
      /* Form object names */
      strcpy(object,"lig");
      strcpy(object_type,"ligand");

      /* For binding sites, just use a sequential number */
      if (params.binding_sites == TRUE)
        {
          /* Increment ligand counter, and add to object name */
          nlig++;
          sprintf(number_string,"%2d",nlig);
          if (number_string[0] == ' ')
            number_string[0] = '0';
          strcat(object,number_string);
        }

      /* Determine labels required, if any */
      if (params.label_atoms == TRUE)
        label_atoms = TRUE;
      if (params.label_ligand == TRUE)
        label_residues = TRUE;
    }
  else
    {
      /* Form object names */
      strcpy(object,"prot");
      strcpy(object_type,"protein");

      /* For binding sites, just use a sequential number */
      if (params.binding_sites == TRUE)
        {
          /* Increment protein-residues counter, and add to object name */
          nprot++;
          sprintf(number_string,"%2d",nprot);
          if (number_string[0] == ' ')
            number_string[0] = '0';
          strcat(object,number_string);
        }

      /* Determine labels required, if any */
      if (params.label_prot_residues == TRUE)
        label_residues = TRUE;
    }

  /* Add PDB code to object name, except for binding-sites database */
  if (params.binding_sites == FALSE)
    {
      if (!strcmp(sequence_ptr->pdb_code,"uplo"))
        strcat(object,"query");
      else
        {
          strcat(object,sequence_ptr->pdb_code);

          /* See if this is an Alpha Fold model */
          ch = sequence_ptr->pdb_code[0];
          if (ch >= 'A' && ch <= 'Z')
            type = ALPHA_FOLD_PDB;
        }

      if (sequence_ptr->chain != ' ')
        strcat(object,chain_str);

      /* Add domain identifier, if there is one */
      if (sequence_ptr->domain_no != 0)
        {
          /* Convert domain number into character-string and add on */
          sprintf(domain_str,"_d%d",sequence_ptr->domain_no);
          strcat(object,domain_str);
        }
    }

  /* Set appropriate atom-colour */

  /* SAS colouring scheme */
  if (params.sas_colour == TRUE)
    {
      /* If plotting ligands, then make colour cpk */
      if (in_ligand == TRUE)
        {
          strcpy(sequence_ptr->colour_name,"cpk");
          strcpy(sequence_ptr->colour_numbers,"cpk");
        }

      /* Check whether this sequence is all grey */
      else if (sequence_ptr->grey == TRUE)
        {
          strcpy(sequence_ptr->colour_name,"grey");
          strcpy(sequence_ptr->colour_numbers,"[151,151,151]");
        }

      /* Otherwise, set default colour as white */
      else
        {
          strcpy(sequence_ptr->colour_name,"white");
          strcpy(sequence_ptr->colour_numbers,"white");
        }
    }

  /* Separate colour for each structure */
  else if (params.colour_by_structure == TRUE)
    {
      /* Check whether we have a previous file for the same PDB
         code, chain and domain-number */
      match = FALSE;

      /* Loop through all previous entries */
      for (jfile = 0; jfile < ifile && match == FALSE; jfile++)
        {
          /* Get the sequence pointer for the other sequence */
          other_sequence_ptr = stored_sequence_ptr[jfile];

          /* Check for match */
          if (!strcmp(other_sequence_ptr->pdb_code,pdb_code) &&
              other_sequence_ptr->chain == chain &&
              other_sequence_ptr->domain_no == domain_no)
            {
              /* Have found a match */
              match = TRUE;

              /* Transfer the colour across */
              strcpy(sequence_ptr->colour_name,
                     other_sequence_ptr->colour_name);
              strcpy(sequence_ptr->colour_numbers,
                     other_sequence_ptr->colour_numbers);

              /* Blank out label so that separate entry not shown in the
                 key */
              sequence_ptr->label[0] = '\0';
            }
        }

      /* If no prior match, then get the colour for this file */
      if (match == FALSE)
        {
          colour = ifile % MAX_OBJ_COL;
          if (params.binding_sites == TRUE)
            {
              strcpy(sequence_ptr->colour_name,col_nameb[colour]);
              strcpy(sequence_ptr->colour_numbers,col_numbersb[colour]);
            }
          else
            {
              strcpy(sequence_ptr->colour_name,col_name[colour]);
              strcpy(sequence_ptr->colour_numbers,col_numbers[colour]);
            }
        }
    }

  /* PROSITE pattern or PROMOTIF so set default colour to white */
  else if (params.prosite_pattern == TRUE || params.promotif == TRUE ||
           params.residue_list[0] != '\0')
    {
      /* Set default colour to white */
      strcpy(sequence_ptr->colour_name,"white");
      strcpy(sequence_ptr->colour_numbers,"white");
    }
  else
    {
      strcpy(sequence_ptr->colour_name,"cpk");
      strcpy(sequence_ptr->colour_numbers,"cpk");
    }

  /* Write out heading */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Definition for %s\n",object);

  /* Start the output line */
  strcpy(output_line,"define ");
  strcat(output_line,object);

  /* Initialise pointer to the first record */
  atom_ptr = sequence_ptr->first_atom_ptr;

  /* Loop through all the atom records */
  while (atom_ptr != NULL)
    {
      /* Get this atom's residue pointer */
      residue_ptr = atom_ptr->residue_ptr;

      /* Check whether atom belongs to the current object type */
      if (residue_ptr->is_ligand == in_ligand)
        {
          /* If this is the first atom of this range, then store */
          if (from_atom == 0)
            from_atom = atom_ptr->out_number;

          /* Store end-atom of range */
          to_atom = atom_ptr->out_number;

          /* Check if residue has changed */
          if (strncmp(residue_ptr->res_name,last_residue,3) ||
              strncmp(residue_ptr->res_num,last_res_number,5) ||
              residue_ptr->chain != last_chain)
            {
              /* Set flag */
              new_residue = TRUE;

              /* If labelling atoms or residues, form the residue-label */
              if (label_atoms == TRUE || (label_residues == TRUE &&
                                          residue_ptr->in_list != NOT_IN_LIST))
                {
                  /* Form label for the current atom */
                  strcpy(residue_label,"%n%r");
                  if (residue_ptr->chain != ' ')
                    {
                      sprintf(full_chain,"(%%c)");
                      strcat(residue_label,full_chain);
                    }
                }
            }
          else
            new_residue = FALSE;

          /* If labelling all atoms, then label this one */
          if (label_atoms == TRUE)
            {
              /* Form label for the current atom, adding residue name
                 to first atom of the residue */
              if (new_residue == TRUE)
                {
                  strcpy(label,"-%a:");
                  strcat(label,residue_label);
                }
              else
                strcpy(label,"-%a");

              /* Select current atom and apply label to it */
              if (params.pymol == TRUE)
                {
                  fprintf(outfile_ptr,"create label%d, id %d\n",
                          atom_ptr->out_number,atom_ptr->out_number);
                  fprintf(outfile_ptr,"label label%d, \"%s\"\n",
                          atom_ptr->out_number,label);
                }
              else
                {
                  fprintf(outfile_ptr,"select atomno = %d\n",
                          atom_ptr->out_number);
                  fprintf(outfile_ptr,"label %s\n",label);
                }
            }

          /* Otherwise, if labelling residues, and this is a new one, then
             label the atom */
          else if (new_residue == TRUE && label_residues == TRUE &&
                   residue_ptr->in_list != NOT_IN_LIST)
            {
              if (params.pymol == TRUE)
                {
                  /* Select current atom and apply label to it */
                  fprintf(outfile_ptr,"create label%d, id %d\n",
                          atom_ptr->out_number,atom_ptr->out_number);
                  fprintf(outfile_ptr,"label label%d, \"%s\"\n",
                          atom_ptr->out_number,residue_label);
                }
              else
                {
                  /* Select current atom and apply label to it */
                  fprintf(outfile_ptr,"select atomno = %d\n",
                          atom_ptr->out_number);
                  fprintf(outfile_ptr,"label %s\n",residue_label);
                }
            }

          /* Store the current residue details */
          strcpy(last_residue,residue_ptr->res_name);
          strcpy(last_res_number,residue_ptr->res_num);
          last_chain = residue_ptr->chain;
        }

      /* Otherwise check if have changed */
      else if (residue_ptr->is_ligand != last_in_ligand && from_atom != 0)
        {
          /* Add this residue range to the end of the string */
          nranges = add_residue_range(outfile_ptr,object,output_line,
                                      from_atom,to_atom,nranges);

          /* Reset atom ranges */
          from_atom = 0;
          to_atom = 0;
        }

      /* Store whether atom in ligand or not */
      last_in_ligand = residue_ptr->is_ligand;

      /* Get pointer to the next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }

  /* Add final atom-range to list of ranges, if there is one */
  if (from_atom != 0)
    nranges = add_residue_range(outfile_ptr,object,output_line,
                                from_atom,to_atom,nranges);

  /* If have defined current object, then define and select it */
  if (nranges > 0)
    {
      /* Select and display the object according to type */
      fprintf(outfile_ptr,"select %s\n",object);

      /* If only a single ligand, then display as thick bonds or as
         ball-and-stick */
      if (in_ligand == TRUE && nligands < 2)
        {
          /* Thick bonds */
          if (params.display_ligand_only == TRUE ||
              params.display_ligand_and_residues == TRUE)
            {
              fprintf(outfile_ptr,"wireframe 0.4\n");
              if (params.metal == TRUE)
                fprintf(outfile_ptr,"spacefill 0.4\n");
            }

          /* Thick ball-and-stick */
          else
            {
              fprintf(outfile_ptr,"wireframe 0.3\n");
              fprintf(outfile_ptr,"spacefill 0.6\n");
            }
        }

      /* If multiple ligands, display each one as thin bonds or 
         ball-and-stick */
      else if (in_ligand == TRUE && nligands > 1)
        {
          /* Thin bonds */
          if (params.display_ligand_only == TRUE ||
              params.display_ligand_and_residues == TRUE)
            fprintf(outfile_ptr,"wireframe 0.3\n");

          /* Thin ball-and-stick */
          else
            {
              fprintf(outfile_ptr,"wireframe 0.22\n");
              fprintf(outfile_ptr,"spacefill 0.4\n");
            }
        }

      /* If only a single ligand, and these are the protein residues,
         display as thin sticks */
      else if (in_ligand == FALSE)
        {
          if (sequence_ptr->display_option == ALL_ATOMS ||
              sequence_ptr->display_option == MAIN_CHAIN)
            {
              if (params.pymol == TRUE)
                fprintf(outfile_ptr,"show lines, all\n");
              else
                fprintf(outfile_ptr,"wireframe on\n");
            }

          else if (sequence_ptr->display_option == C_ALPHA)
            {
              if (params.pymol == TRUE)
                fprintf(outfile_ptr,
                        "show lines, (name ca or name c or name n)\n");
              else
                fprintf(outfile_ptr,"backbone on\n");
            }

          else if (sequence_ptr->display_option == BINDING_SITE  ||
                   sequence_ptr->display_option == INT_RESIDUES ||
                   sequence_ptr->display_option == LIGAND_AND_RES)
            fprintf(outfile_ptr,"wireframe 0.1\n");

          else if (sequence_ptr->display_option == CARTOON)
            {
              if (params.pymol == TRUE)
                fprintf(outfile_ptr,"show cartoon, all\n");
              else
                fprintf(outfile_ptr,"cartoon on\n");
            }
        }

      /* If this an Alpha Fold structure, colour by B-factor */
      if (type == ALPHA_FOLD_PDB)
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
              fprintf(outfile_ptr,"select %s and (temperature > %d and "
                      "temperature <= %d)\n",object,lower_b,upper_b);
              fprintf(outfile_ptr,"colour [%d,%d,%d]\n",
                      AF_COLOUR_IRGB[i][0],
                      AF_COLOUR_IRGB[i][1],
                      AF_COLOUR_IRGB[i][2]);

              /* Define the colour name */
              strcpy(colour_name,"coloured by pred confidence");
            }
        }

      /* Otherwise, use standard colouring */
      else
        {
          /* Get default colour-name */
          strcpy(colour_name,sequence_ptr->colour_name);
          strcpy(colour_numbers,sequence_ptr->colour_numbers);
          if ((params.prosite_pattern == TRUE ||
               params.promotif == TRUE ||
               params.residue_list[0] != '\0') && in_ligand == TRUE)
            {
              strcpy(colour_name,"cpk");
              strcpy(colour_numbers,"cpk");
            }

          /* Write out the colour */
          fprintf(outfile_ptr,"colour atoms %s\n",colour_numbers);
        }

      /* For RasMol script, show the defiitions */
      if (params.pymol == FALSE)
        {
          /* Make length of object name a uniform length */
          len = strlen(object);
          for (ipos = 0; ipos < 12 - len; ipos++)
            strcat(object," ");

          /* Show the definition on screen */
          fprintf(outfile_ptr,"echo \"  %s ",object);
          if (in_ligand == TRUE)
            fprintf(outfile_ptr," ");
          if (!strcmp(sequence_ptr->pdb_code,"uplo"))
            fprintf(outfile_ptr,"= %s: query",object_type);
          else
            fprintf(outfile_ptr,"= %s: %s",object_type,sequence_ptr->pdb_code);
          if (sequence_ptr->chain != ' ')
            fprintf(outfile_ptr," chain %c",sequence_ptr->chain);
          if (sequence_ptr->domain_no != 0)
            fprintf(outfile_ptr,", domain %d",sequence_ptr->domain_no);
          if (params.sas_colour == FALSE)
            fprintf(outfile_ptr," (%s)",colour_name);
          
          /* If have a RMSD value for this structure, show it */
          if (params.cath_alignment == TRUE && ifile > 0 &&
              !strncmp(object_type,"protein",7))
            {
              if (rmsd < 1.0)
                sprintf(rmsd_out,"%5.3f",rmsd);
              else if (rmsd < 10.0)
                sprintf(rmsd_out,"%4.2f",rmsd);
              else if (rmsd < 100.0)
                sprintf(rmsd_out,"%5.2f",rmsd);
              else if (rmsd < 1000.0)
                sprintf(rmsd_out,"%6.2f",rmsd);
              else
                sprintf(rmsd_out,"%7.2f",rmsd);
              fprintf(outfile_ptr,"  RMSD = %s (%d atoms)",rmsd_out,nfit);
            }

          /* Close off the current echo statement */
          fprintf(outfile_ptr,"\"\n");
        }
    }

  /* Otherwise, reduce appropriate counter */
  else
    {
      if (in_ligand == TRUE)
        nlig--;
      else
        nprot--;
    }

  /* Store current object name */
  strcpy(oname,object);
}
/***********************************************************************

get_exon_colour  -  Retrieve the colour corresponding to this exon

***********************************************************************/

void get_exon_colour(char *colour,int iexon)
{
#define NEXCOL 14

  int icol;

  static char *col_name[] = {
    "purple", "yellow", "red", "blue", "brown", "green", "pink", "cyan",
    "lilac", "grey", "orange", "white", "sky_blue", "cream"
  };

  /* Initialise variables */
  icol = iexon % NEXCOL;

  /* Get the corresponding colour colour */
  strcpy(colour,col_name[icol]);
}
/***********************************************************************

get_rasmol_colour  -  Retrieve the corresponding RasMol colour triplet
                      from the supplied colour name

***********************************************************************/

void get_rasmol_colour(char *colour_name,char *colour_numbers)
{
#define NRASMOL_COLOURS 23

  int icol;
  int match;

  static char *col_name[] = {
    "black", "white", "cream", "paleblue",
    "purple", "yellow", "blue", "red", "green",
    "brown", "pink", "greenblue", "magenta",
    "cyan", "redorange", "violet", "orange",
    "peach", "sky_blue", "grey", "dark_grey", "lime_green", "lilac"
  };

  static char *col_numbers[] = {
    "[0,0,0]", "[255,255,255]", "[255,255,179]", "[173,216,230]",
    "[178,51,255]", "[255,255,0]", "[0,0,255]", "[255,0,0]", "[0,255,0]",
    "[128,0,0]", "[255,128,255]", "[51,153,102]", "[255,0,255]",
    "[128,255,255]", "[255,102,0]", "[238,130,238]", "[204,128,0]",
    "[255,204,153]", "[0,153,255]", "[128,128,128]", "[51,51,51]",
    "[128,255,0]", "[255,0,255]"
  };

  /* Initialise variables */
  match = FALSE;
  strcpy(colour_numbers,col_numbers[0]);

  /* Loop until find colour */
  for (icol = 0; icol < NRASMOL_COLOURS && match == FALSE; icol++)
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

write_script_reslist  -  Write out the script commands to define the
                         user-defined residue lists

***********************************************************************/

void write_script_reslist(FILE *outfile_ptr,struct parameters params,
                          char *object)
{
  char range[50];
  char c_name[COLNAM_LEN], colour_numbers[COLNO_LEN];

  int ipos, label_atom, len;
  int iatom, natoms, from_atom, to_atom;
  int done, first;

  struct atom *atom_ptr;
  struct residue *residue_ptr;
  struct reslist *reslist_ptr;

  /* Initiailise variables */
  ipos = 0;
  done = FALSE;
  first = TRUE;

  /* If have at least one residue start-end pair, write out heading */
  if (params.first_reslist_ptr != NULL)
    {
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"# Residue list\n");
    }

  /* Start at first residue range */
  reslist_ptr = params.first_reslist_ptr;

  /* Loop through the residue list */
  while (reslist_ptr != NULL)
    {
      /* If have start and end residues, then proceed */
      if (reslist_ptr->residue1_ptr != NULL && reslist_ptr->residue2_ptr)
        {
          /* Get the first residue */
          residue_ptr = reslist_ptr->residue1_ptr;

          /* Get first atom */
          from_atom = residue_ptr->first_atom_no;

          /* Get the last residue */
          residue_ptr = reslist_ptr->residue2_ptr;

          /* Get last atom */
          to_atom = residue_ptr->last_atom_no;

          /* Define residue range */
          sprintf(range,"(atomno >= %d & atomno <= %d)",from_atom,
                  to_atom);

          /* If first time, define residue range using given name */
          if (first == TRUE)
            fprintf(outfile_ptr,"define %s %s\n",object,range);
          else
            fprintf(outfile_ptr,"define %s %s | %s\n",object,object,range);
          first = FALSE;

          /* Select the residue range */
          fprintf(outfile_ptr,"select %s\n",range);

          /* Colour residues */
          if (!strcmp(params.rlist_colour,"exons"))
            get_exon_colour(c_name,reslist_ptr->number);
          else
            strcpy(c_name,params.rlist_colour);
          if (!strcmp(c_name,"cpk"))
            fprintf(outfile_ptr,"colour cpk\n");
          else
            {
              get_rasmol_colour(c_name,colour_numbers);
              fprintf(outfile_ptr,"colour %s\n",colour_numbers);
            }

          /* Show as thicker bonds */
          if (params.display_backbone_only == TRUE)
            fprintf(outfile_ptr,"backbone 0.3\n");
          else if (params.display_mainchain_only == TRUE)
            fprintf(outfile_ptr,"wireframe 0.3\n");

          /* If residues to be labelled, label first residue */
          if (params.csa == TRUE)
            {
              /* Initialise atom number for label atom */
              label_atom = 0;

              /* Get the first residue */
              residue_ptr = reslist_ptr->residue1_ptr;

              /* Get this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              natoms = residue_ptr->natoms;
              iatom = 0;
              done = FALSE;

              /* Loop through the residue's atoms */
              while (atom_ptr != NULL && iatom < natoms && done == FALSE)
                {
                  /* Save the atom number */
                  label_atom = atom_ptr->out_number;

                  /* If this is the mainchain nitrogen, end here */
                  if (!strcmp(atom_ptr->atom_name," N  "))
                    done = TRUE;

                  /* Get pointer to next atom in linked-list */
                  atom_ptr = atom_ptr->next_atom_ptr;
                  iatom++;
                }
              
              /* If have a label atom, write out the label */
              if (label_atom != 0)
                {
                  /* Select the atom */
                  fprintf(outfile_ptr,"select atomno = %d\n",label_atom);
                  if (reslist_ptr->chain != ' ')
                    fprintf(outfile_ptr,"label %%n%%r(%%c)\n");
                  else
                    fprintf(outfile_ptr,"label %%n%%r\n");
                }
            }
        }

      /* Go to next range */
      reslist_ptr = reslist_ptr->next_reslist_ptr;
    }
}
/***********************************************************************

write_script_fitlist  -  Write out the script commands to highlight
                         residues on which structures were fitted

***********************************************************************/

void write_script_fitlist(struct sequence *sequence_ptr,
                          FILE *outfile_ptr,char *pname,
                          int ifile,struct parameters params)
{
  char range[10];

  struct fitlist *fitlist_ptr;

  /* Write out heading */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Fitted residues list\n");

  /* Get the first of this structure's fitting ranges */
  fitlist_ptr = sequence_ptr->first_fitlist_ptr;

  /* Loop over all the fit-lists for this structure */
  while (fitlist_ptr != NULL)
    {
      /* Remove all blanks from this list */
      strcpy(range,fitlist_ptr->range);
      remove_blanks(range);

      /* Select the residue range */
      fprintf(outfile_ptr,"select %s & %s\n",pname,range);

      /* Show as thicker bonds */
      if (params.display_backbone_only == TRUE)
        fprintf(outfile_ptr,"backbone 0.3\n");
      else if (params.display_mainchain_only == TRUE)
        fprintf(outfile_ptr,"wireframe 0.3\n");

      /* Get the next fitlist record */
      fitlist_ptr = fitlist_ptr->next_fitlist_ptr;
    }
}
/***********************************************************************

form_residue_def  -  Form the residue definition in RasMol script format
                     (using atom-numbers as current versions of RasMol
                     cannot handle residue insertion codes)

***********************************************************************/

void form_residue_def(char chain,char res_name[4],char res_num[6],
                      char residue_def[16])
{
  char ch, temp_name[10];

  int ipos, fst_pos, len, lst_pos;
  int number;

  /* Initialise variables */
  residue_def[0] = '\0';

  /* Check whether the residue name contains any
     non-alphabetic characters (eg numbers as in SO4
     groups) */
  number = FALSE;
  fst_pos = -1;
  lst_pos = 0;
  for (ipos = 0; ipos < 3; ipos++)
    {
      /* Check if character falls outside alphabet */
      ch = res_name[ipos];
      if (ch < 'A' || ch > 'Z')
        number = TRUE;

      /* Check for blanks */
      if (ch != ' ')
        {
          lst_pos = ipos;
          if (fst_pos == -1)
            fst_pos = ipos;
        }
    }

  /* Form just the residue name, less any blanks */
  len = lst_pos - fst_pos + 1;
  strncpy(temp_name,res_name+fst_pos,len);
  temp_name[len] = '\0';

  /* Add residue name to definition, with square
     brackets, if necessary */
  if (number == TRUE)
    strcat(residue_def,"[");
  strcat(residue_def,temp_name);
  if (number == TRUE)
    strcat(residue_def,"]");

  /* Form the residue number */
  fst_pos = 0;
  if (res_num[0] == ' ')
    fst_pos = 1;
  if (res_num[1] == ' ')
    fst_pos = 2;
  if (res_num[2] == ' ')
    fst_pos = 3;
  strcpy(temp_name,res_num+fst_pos);

  /* Currently RasMol can't cope with residue 
     insertion codes (interpreting them as the chain-id),
     chain-id), so exclude the insertion code */
  /*   if (res_num[4] == ' ')
       temp_name[4 - fst_pos] = '\0';
       else
       temp_name[5 - fst_pos] = '\0'; */
  temp_name[4 - fst_pos] = '\0';

  /* Add the residue number to the residue definition */
  strcat(residue_def,temp_name);

  /* If there is a chain-id, then add that */
  if (chain != ' ')
    {
      temp_name[0] = chain;
      temp_name[1] = '\0';
      strcat(residue_def,":");
      strcat(residue_def,temp_name);
    }
}
/***********************************************************************

write_script_sascolour  -  Write out the SAS colouring for each residue
                           in the current structure

***********************************************************************/

void write_script_sascolour(FILE *outfile_ptr,struct sequence *sequence_ptr)
{
  char chain, res_name[4], res_num[6];
  char residue_def[16];

  static char *colour_name[] = { 
    "white", "red", "green", "blue", "yellow",
    "purple", "brown", "skyblue", "orange",
    "cyan", "grey", "black"
  };

  static char *colour_numbers[] = {
    "[255,255,255]", "[255,0,0]", "[0,255,0]", "[0,0,255]", "[255,255,0]",
    "[178,51,255]", "[128,0,0]", "[77,204,255]", "[204,128,0]",
    "[128,255,255]", "[151,151,151]", "[0,0,0]"
  };

  int colour;

  static int ncolours = 12;

  struct atom *atom_ptr, *next_atom_ptr;
  struct residue *next_residue_ptr, *residue_ptr;

  /* Write out heading */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# SAS colouring of residues\n");

  /* Initialise pointer to the first record */
  atom_ptr = sequence_ptr->first_atom_ptr;

  /* Loop through all the atom records */
  while (atom_ptr != NULL)
    {
      /* Get this atom's residue pointer */
      residue_ptr = atom_ptr->residue_ptr;

      /* Get the next atom and its residue pointer */
      next_atom_ptr = atom_ptr->next_atom_ptr;
      if (next_atom_ptr != NULL)
        next_residue_ptr = next_atom_ptr->residue_ptr;
      else
        next_residue_ptr = NULL;

      /* If this is the last atom of the current residue, then write out
         colour of residue */
      if (residue_ptr != next_residue_ptr && residue_ptr != NULL)
        {
          /* Check that the colour isn't the default */
          if (residue_ptr->colour != 0)
            {
              /* Get the residue's name, number and chain-id */
              strcpy(res_name,residue_ptr->res_name);
              strcpy(res_num,residue_ptr->res_num);
              chain = residue_ptr->chain;

              /* Form the residue definition in RasMol script format */
              form_residue_def(chain,res_name,res_num,residue_def);

              /* Select the given residue */
              fprintf(outfile_ptr,"select %s\n",residue_def);

              /* Check that colour number doesn't exceed number of
                 colours we have here */
              colour = residue_ptr->colour;
              if (colour < 0)
                colour = - colour;
              if (colour > ncolours - 1)
                colour = ncolours - 1;

              /* Write out the appropriate colour */
              fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
                      colour_numbers[colour],colour_name[colour]);
            }
        }

      /* Get pointer to the next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

write_script_sidechain  -  Write out the selected sidechains in the
                           given colour

***********************************************************************/

void write_script_sidechain(FILE *outfile_ptr,struct sequence *sequence_ptr,
                            struct parameters params)
{
  char atom_range[LINELEN + 1], colour_numbers[COLNO_LEN];
  char c_name[COLNAM_LEN];
  char object[15];

  static char *prosite_col_name[] = {
    "white", "red", "redorange", "orange",
    "yellow", "green", "greenblue", "purple",
    "skyblue", "blue"
  };

  static char *prosite_col_numbers[] = {
    "[255,255,255]", "[255,0,0]", "[255,102,0]", "[204,128,0]",
    "[255,255,0]", "[0,255,0]", "[51,153,102]", "[178,51,255]",
    "[77,204,255]", "[0,0,255]"
  };

  static char *sas_col_name[] = { 
    "white", "red", "green", "blue", "yellow",
    "purple", "brown", "skyblue", "orange",
    "cyan", "grey", "black", "magenta"
  };

  static char *sas_col_numbers[] = {
    "[255,255,255]", "[255,0,0]", "[0,255,0]", "[0,0,255]", "[255,255,0]",
    "[178,51,255]", "[128,0,0]", "[77,204,255]", "[204,128,0]",
    "[128,255,255]", "[151,151,151]", "[0,0,0]", "[240,0,108]"
  };

  int colour, from_atom, last_atom, previous_first, to_atom;

  struct residue *residue_ptr;

  /* Initialise variables */
  last_atom = 0;
  object[0] = '\0';
  previous_first = 0;

  /* Write out the het-groups in spacefill mode */
  if (params.prosite_pattern == TRUE)
    {
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"# Ligands\n");
      fprintf(outfile_ptr,"select hetero\n");
      fprintf(outfile_ptr,"colour atoms cpk\n");
      fprintf(outfile_ptr,"spacefill 1.2\n");
      fprintf(outfile_ptr,"\n");
    }

  /* Write out heading */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Residue sidechains\n");

  /* Get pointer to the first of this structure's residues */
  residue_ptr = sequence_ptr->first_residue_ptr;

  /* Loop through all the currently held residues for this structure to
     check whether we already have this residue */
  while (residue_ptr != NULL)
    {
      /* Determine whether this residue's sidechain is to be written out */
      if (residue_ptr->is_ligand == FALSE &&
          residue_ptr->in_list != NOT_IN_LIST)
        {
          /* Get the atom-range defining this residue */
          from_atom = residue_ptr->first_atom_no;
          to_atom = residue_ptr->last_atom_no;

          /* Select this residue using its atom-number range */
          atom_range[0] = '\0';
          add_residue_range(outfile_ptr,object,atom_range,from_atom,
                            to_atom,0);

          /* Write out the select command */
          fprintf(outfile_ptr,"select %s\n",atom_range);

          /* Write out the appropriate colour */
          if ((params.prosite_pattern == TRUE || params.side_sas == TRUE) &&
              residue_ptr->colour != 0)
            {
              /* Get the appropriate colour */
              colour = residue_ptr->colour;
              if (params.prosite_pattern == TRUE && colour > 9)
                colour = 9;
              if (params.side_sas == TRUE && colour > 13)
                colour = 13;

              /* Write out the colour for PROSITE pattern */
              if (params.prosite_pattern == TRUE)
                {
                  fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
                          prosite_col_numbers[colour],
                          prosite_col_name[colour]);

                  /* Set appropriate display mode */
                  fprintf(outfile_ptr,"wireframe 0.2\n");
                }
              else
                {
                  /* Bonds coloured according to SAS */
                  fprintf(outfile_ptr,"colour bonds %s  # (%s)\n",
                          sas_col_numbers[colour],sas_col_name[colour]);
                  fprintf(outfile_ptr,"wireframe 0.2\n");

                  /* Atoms in spacefill and coloured in CPK */
                  fprintf(outfile_ptr,"spacefill 0.5\n");
                  fprintf(outfile_ptr,"colour cpk\n");

                  /* Show carbons in dark grey */
                  fprintf(outfile_ptr,"select %s & carbon\n",atom_range);
                  fprintf(outfile_ptr,"colour %s  # (%s)\n",
                          sas_col_numbers[10],sas_col_name[10]);

                  /* For glycine, as only the CA atom is shown, pick it
                     out in the SAS colour */
                  if (!strcmp(residue_ptr->res_name,"GLY"))
                    fprintf(outfile_ptr,"colour %s  # (%s)\n",
                            sas_col_numbers[colour],sas_col_name[colour]);
                }

              /* If last residue was also shown, then join the two */
              if (last_atom == from_atom - 1)
                {
                  /* Select previous residue, plus this residue, using
                     their combined atom-number range */
                  atom_range[0] = '\0';
                  add_residue_range(outfile_ptr,object,atom_range,
                                    previous_first,to_atom,0);

                  /* Select the two atoms and join them */
                  fprintf(outfile_ptr,"select %s\n",atom_range);
                  fprintf(outfile_ptr,"wireframe 0.2\n");
                }

              /* Save last atom written out */
              previous_first = from_atom;
              last_atom = to_atom;
            }

          /* For residue list, show sidechains in given colour with
             thick bonds */
          else if (params.residue_list[0] != '\0' &&
                   residue_ptr->in_list != NOT_IN_LIST)
            {
              /* Get the RasMol numbers for the given colour */
              if (!strcmp(params.rlist_colour,"cpk"))
                strcpy(colour_numbers,"cpk");
              else
                {
                  //if (!strcmp(params.rlist_colour,"exons"))
                  //  get_exon_colour(c_name,iexon);
                  //else
                  strcpy(c_name,params.rlist_colour);
                  get_rasmol_colour(c_name,colour_numbers);
                }

              /* Write out as thick stick */
              fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
                      colour_numbers,c_name);
              fprintf(outfile_ptr,"wireframe 0.2\n");
            }

          /* Default */
          else
            {
              /* Write out default colour and display mode */
              fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
                      sequence_ptr->colour_numbers,
                      sequence_ptr->colour_name);
              fprintf(outfile_ptr,"wireframe on\n");
            }
        }

      /* Get pointer to the next residue record */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int check_for_metal(char atom_name[5])
{
  int imetal;
  int metal;

  static int nmetals = 25;

  static char *metal_name[] = {
    "CA  ", "FE  ", "ZN  ", " U  ", "CU  ",
    "MN  ", "CL  ", "MG  ", "NA  ", "HG  ",
    "CO  ", " K  ", "CD  ", "AU  ", "AG  ",
    "IN  ", "HO  ", " F  ", "PB  ", "AL  ",
    "YB  ", "BR  ", "....", "....", "...."
  };

  /* Initialise variables */
  metal = FALSE;

  /* Loop through all the metals */
  for (imetal = 0; imetal < nmetals && metal == FALSE; imetal++)
    {
      /* If name matches, then have a metal */
      if (!strncmp(atom_name,metal_name[imetal],2))
        metal = TRUE;
    }

  /* Return the answer */
  return(metal);
}
/***********************************************************************

write_script_metals  -  Write out the script commands to show any metal
                        ions interacting with the ligand

***********************************************************************/

void write_script_metals(FILE *outfile_ptr,struct sequence *sequence_ptr)
{
  int natoms;
  int is_metal;

  struct atom *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first record */
  atom_ptr = sequence_ptr->first_atom_ptr;

  /* Loop through all the atom records */
  while (atom_ptr != NULL)
    {
      /* Get this atom's residue pointer */
      residue_ptr = atom_ptr->residue_ptr;

      /* Get this residue's number of atoms */
      natoms = residue_ptr->natoms;

      /* If not a ligand residue, and has only a single atom, it could
         be a metal */
      if (residue_ptr->is_ligand == FALSE && natoms == 1 &&
          strncmp(atom_ptr->atom_name," CA ",4) &&
          strncmp(residue_ptr->res_name,"HOH",3))
        {
          /* Check if this is a metal ion */
          is_metal = check_for_metal(atom_ptr->atom_name);

          /* Select current atom and display as spacefill */
          if (is_metal == TRUE)
            {
              fprintf(outfile_ptr,"select atomno = %d\n",
                      atom_ptr->out_number);
              fprintf(outfile_ptr,"spacefill 0.6\n");
            }
        }

      /* Get pointer to the next atom record */
      atom_ptr = atom_ptr->next_atom_ptr;
    }
}
/***********************************************************************

write_script_hbonds  -  Write out the script commands to show the
                        hydrogen bonds

***********************************************************************/

int write_script_hbonds(FILE *outfile_ptr,struct sequence *sequence_ptr,
                        int nbonds)
{
  struct hbond *hbond_ptr;

  static char *hbond_colour[] = { "[0,255,255]" };

  /* Write out heading */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# H-bond definitions for %s",
          sequence_ptr->pdb_code);
  if (sequence_ptr->chain != ' ')
    fprintf(outfile_ptr,"(%c)",sequence_ptr->chain);
  fprintf(outfile_ptr,"\n");

  /* Get pointer to the first H-bond */
  hbond_ptr = sequence_ptr->first_hbond_ptr;

  /* Loop through all the H-bonds to label each one with its length */
  while (hbond_ptr != NULL)
    {
      /* Increment count of bonds read */
      nbonds++;

      /* Show all H-bonds as a sky-blue line */
      fprintf(outfile_ptr,"select HB%d\n",nbonds);
      fprintf(outfile_ptr,"colour %s\n",hbond_colour[0]);
      fprintf(outfile_ptr,"wireframe on\n");

      /* Select current bond and assign label to the mid-point atom */
      fprintf(outfile_ptr,"select HB%d.C\n",nbonds);
      fprintf(outfile_ptr,"label %4.2f\n",hbond_ptr->length);

      /* Get pointer to next H-bond */
      hbond_ptr = hbond_ptr->next_hbond_ptr;
    }

  /* Return number of bonds */
  return(nbonds);
}
/***********************************************************************

colour_chains  -  Write the script commands to colour the chains

***********************************************************************/

void colour_chains(FILE *outfile_ptr,struct parameters params)
{
#define NCOLS 12

  char object[15], range_name[NAME_LEN + 1];

  int ichain;

  static char *col_name[NCOLS] = {
    "white", "green", "yellow", "purple",
    "brown", "pink", "greenblue", "magenta",
    "cyan", "redorange", "violet", "orange"
  };

  static char *col_numbers[NCOLS] = {
    "[255,255,255]", "[0,255,0]", "[255,255,0]", "[178,51,255]",
    "[128,0,0]", "[255,128,255]", "[51,153,102]", "[255,0,255]",
    "[128,255,255]", "[255,102,0]", "[238,130,238]", "[204,128,0]"
  };

  /* Loop over the chains */
  for (ichain = 0; ichain < params.nchains; ichain++)
    {
      /* Give a name to this chain */
      if (params.chains[ichain] == ' ')
        strcpy(object,"chain_");
      else
        sprintf(object,"chain%c",params.chains[ichain]);
      sprintf(range_name," *:%c",params.chains[ichain]);

      /* Write out the definition */
      fprintf(outfile_ptr,"\n");
      fprintf(outfile_ptr,"# Definition for %s\n",object);

      /* Define object */
      fprintf(outfile_ptr,"define %s %s\n",object,range_name);

      /* Select it */
      fprintf(outfile_ptr,"select %s\n",object);

      /* Colour and display it */
      fprintf(outfile_ptr,"colour %s\n",col_numbers[ichain % NCOLS]);
      fprintf(outfile_ptr,"wireframe on\n");

      /* Show the definition on screen */
      strcat(object,"      ");
      fprintf(outfile_ptr,"echo \"  %s = chain %c (%s)\"\n",object,
              params.chains[ichain],col_name[ichain % NCOLS]);
    }
}
/***********************************************************************

show_ligands  -  Show any ligands that are in the structure

***********************************************************************/

void show_ligands(FILE *outfile_ptr,struct rasdef *first_rasdef_ptr)
{
  char lobject[15], mobject[15], range_name[NAME_LEN + 1];

  int nlig, nmet;

  struct rasdef *rasdef_ptr;

  /* Initialise variables */
  nlig = 0;
  nmet = 0;

  /* Form name of the site */
  strcpy(lobject,"lig");
  strcpy(mobject,"metals");

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      /* Check for a match to a ligand */
      if (rasdef_ptr->type == LIGAND && rasdef_ptr->count == 1)
        {
          /* Add to ligands definition */
          if (nlig == 0)
            fprintf(outfile_ptr,"define %s %s\n",lobject,
                    rasdef_ptr->rasmol_ligand_def);
          else
            fprintf(outfile_ptr,"define %s %s | %s\n",lobject,lobject,
                    rasdef_ptr->rasmol_ligand_def);

          /* Increment count of ligands stored */
          nlig++;
        }

      /* Check for a match to a metal */
      else if (rasdef_ptr->type == METAL && rasdef_ptr->count == 1)
        {
          /* Add to metals definition */
          if (nmet == 0)
            fprintf(outfile_ptr,"define %s %s\n",mobject,
                    rasdef_ptr->res_num);
          else
            fprintf(outfile_ptr,"define %s %s | %s\n",mobject,mobject,
                    rasdef_ptr->res_num);

          /* Increment count of metals stored */
          nmet++;
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }

  /* If have any ligands, display in spacefill */
  if (nlig > 0)
    {
      fprintf(outfile_ptr,"select %s\n",lobject);
      fprintf(outfile_ptr,"colour cpk\n");
      fprintf(outfile_ptr,"spacefill 1.2\n");
    }

  /* If have any metals, display in spacefill */
  if (nmet > 0)
    {
      fprintf(outfile_ptr,"select %s\n",mobject);
      fprintf(outfile_ptr,"colour cpk\n");
      fprintf(outfile_ptr,"spacefill on\n");
    }
}
/***********************************************************************

write_script_promotif  -  Write out the script commands to show the given
                          PROMOTIF secondary structure motif

***********************************************************************/

void write_script_promotif(FILE *outfile_ptr,
                           struct range *first_range_ptr,
                           struct parameters params)
{
#define NGROUP 3

  char res_num[2][6];
  char colour_numbers[COLNO_LEN];
  char object[15], helix_name[NAME_LEN + 1], range_name[NAME_LEN + 1];

  int i, igroup, ipos, len, nrange, nstrands;

  struct range *range_ptr;

  static char *colour[NGROUP] = {
    "[0,255,0]", "[255,0,0]", "[0,0,255]"
  };
  static char *group[NGROUP] = {
    "(ALA | VAL | PHE | PRO | MET | ILE | LEU | GLY)",
    "(ASP | GLU | LYS | ARG)",
    "(SER | THR | TYR | HIS | CYS | ASN | GLN | TRP)"
  };

  /* Initialise variables */
  helix_name[0] = '\0';
  nrange = 0;
  nstrands = 0;

  /* Display the chains in different colours */
  colour_chains(outfile_ptr,params);

  /* Initialise pointer to the first record */
  range_ptr = first_range_ptr;

  /* Loop through all the range records */
  while (range_ptr != NULL)
    {
      /* Format the residue numbers to remove leading spaces and the
         insertion code (which RasMol can't handle) */
      for (i = 0; i < 2; i++)
        {
          /* Save the residue number and remove insertion code */
          strcpy(res_num[i],range_ptr->res_num[i]);
          res_num[i][4] = '\0';

          /* Remove leading blanks */
          remove_leading_blanks(res_num[i]);
        }

      /* If have no identifier for this range, increment count of ranges
         and generate an identifier */
      if (range_ptr->range_name[0] == '\0')
        {
          /* Increment range */
          nrange++;

          /* Create range identifier */
          sprintf(range_name,"range%d",nrange);
        }

      /* Otherwise, use given range name and display */
      else
        {
          /* Store the name */
          strcpy(range_name,range_ptr->range_name);

          /* Make length of object name a uniform length */
          strcpy(object,range_ptr->range_name);
          len = strlen(object);
          for (ipos = 0; ipos < 12 - len; ipos++)
            strcat(object," ");

          /* Show the definition on screen */
          fprintf(outfile_ptr,"echo \"  %s = %s\"\n",object,
                  range_ptr->name);

          /* If showing single helix, save range name */
          if (params.motif_type == HELICES)
            strcpy(helix_name,range_name);
        }

      /* Define the given residue range */
      fprintf(outfile_ptr,"define %s %s-%s",range_name,res_num[0],
              res_num[1]);
      if (range_ptr->chain != ' ')
        fprintf(outfile_ptr," & *:%c",range_ptr->chain);
      if (range_ptr->sidechain == FALSE)
        fprintf(outfile_ptr," & mainchain");
      fprintf(outfile_ptr,"\n");

      /* Select the given residue range */
      fprintf(outfile_ptr,"select %s\n",range_name);

      /* Define colour */
      get_rasmol_colour(range_ptr->colour,colour_numbers);
      fprintf(outfile_ptr,"colour %s\n",colour_numbers);

      /* Define thickness */
      fprintf(outfile_ptr,"wireframe %3.1f\n",range_ptr->wireframe);

      /* If labelling the residue, select the N atoms of each residue */
      if (range_ptr->label[0] != '\0')
        {
          fprintf(outfile_ptr,"select range%d & *.N\n",nrange);
          fprintf(outfile_ptr,"label %s\n",range_ptr->label);
        }

      /* If H-bonds to be switched on, then do so */
      if (range_ptr->hbonds == TRUE)
        {
          fprintf(outfile_ptr,"select %s\n",range_name);
          fprintf(outfile_ptr,"hbonds on\n");
          fprintf(outfile_ptr,"colour hbonds [128,255,255]\n");
        }
      fprintf(outfile_ptr,"\n");

      /* For beta sheets, add current strand to sheet definition */
      if (params.motif_type == SHEETS)
        {
          /* If this is the first strand, start the sheet definition */
          if (nstrands == 0)
            fprintf(outfile_ptr,"define bsheet %s\n",range_name);

          /* Otherwise, join to previous definition */
          else
            fprintf(outfile_ptr,"define bsheet bsheet | %s\n",range_name);

          /* Increment strand count */
          nstrands++;
        }

      /* Get pointer to the next range record */
      range_ptr = range_ptr->next_range_ptr;
    }

  /* For beta bulge, show H-bonds between the two strands */
  if (params.motif_type == BETA_BULGES)
    {
      /* Switch on the H-bonds between the two strands */
      fprintf(outfile_ptr,"select strand1 | strand2\n");
      fprintf(outfile_ptr,"hbonds on\n");
      fprintf(outfile_ptr,"colour hbonds [128,255,255]\n");
    }

  /* For helix, colour the residues according to the helical wheel
     and helical net plots */
  if (params.motif_type == HELICES && helix_name[0] != '\0')
    {
      /* Loop over the residue groups */
      for (igroup = 0; igroup < NGROUP; igroup++)
        {
          /* Select this group */
          fprintf(outfile_ptr,"select %s & %s\n",helix_name,group[igroup]);

          /* Define colour */
          fprintf(outfile_ptr,"colour %s\n",colour[igroup]);
        }
    }

  /* For beta sheets, define sheet and show H-bonds */
  if (params.motif_type == SHEETS)
    {
      /* Select the beta sheet */
      fprintf(outfile_ptr,"select bsheet\n");

      /* Show the definition on screen */
      fprintf(outfile_ptr,
              "echo \"  bsheet       = %d-strand beta sheet (red)\"\n",
              nstrands);

      /* Switch on the H-bonds */
      fprintf(outfile_ptr,"hbonds on\n");
      fprintf(outfile_ptr,"colour hbonds [128,255,255]\n");
    }
}
/***********************************************************************

write_script_site  -  Write out the script commands to show the given
                      SITE residues

***********************************************************************/

void write_script_site(FILE *outfile_ptr,struct sequence *sequence_ptr,
                       struct rasdef *first_rasdef_ptr,
                       struct parameters params)
{
  char chain_str[2], object[15], list_name[NAME_LEN + 1], res_num[7];

  int ipos, len, nres;

  struct list *list_ptr;

  /* Initialise variables */
  nres = 0;

  /* Display the chains in different colours */
  colour_chains(outfile_ptr,params);

  /* Form name of the site */
  strcpy(object,"site");
  strcat(object,params.site_id);

  /* Write out the definition */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Definition for active site residues (%s)\n",object);

  /* Initialise pointer to the first record */
  list_ptr = sequence_ptr->first_list_ptr;

  /* Loop through all the list records */
  while (list_ptr != NULL)
    {
      /* Form residue identifier */
      strcpy(res_num,list_ptr->res_num);
      string_truncate(res_num,5);
      chain_str[0] = list_ptr->chain;
      chain_str[1] = '\0';
      strcat(res_num,chain_str);

      /* Add residue to current definition of the site */
      if (nres == 0)
        fprintf(outfile_ptr,"define %s %s\n",object,res_num);
      else
        fprintf(outfile_ptr,"define %s %s | %s\n",object,object,res_num);
      nres++;

      /* Get pointer to the next list record */
      list_ptr = list_ptr->next_list_ptr;
    }

  /* Select the active site residues */
  fprintf(outfile_ptr,"select %s\n",object);

  /* Make length of object name a uniform length */
  len = strlen(object);
  for (ipos = 0; ipos < 12 - len; ipos++)
    strcat(object," ");

  /* Show the definition on screen */
  fprintf(outfile_ptr,
          "echo \"  %s = active site %s residues (thick bonds)\"\n",
          object,params.site_id);

  /* Show as thick bonds */
  fprintf(outfile_ptr,"wireframe 0.3\n");
  fprintf(outfile_ptr,"colour cpk\n");

  /* Show any ligands that are in the structure */
  show_ligands(outfile_ptr,first_rasdef_ptr);

  /* Show any waters or metals as spacefill */
  fprintf(outfile_ptr,"select %s & !protein\n",object);
  fprintf(outfile_ptr,"spacefill on\n");
}
/***********************************************************************

get_chn_colour  -  Determine the colour for the given chain

***********************************************************************/

int get_chn_colour(char chain,char *colour_name,char *colour_numbers)
{
#define MAX_CHAIN_COL  8

  int colour, ichain, nchains;

  static char chain_list[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890";
  static char *col_name[] = {
    "purple", "red", "brown", "pink", "blue", "green", "cyan", "yellow"
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
          colour = colour % MAX_CHAIN_COL;
          strcpy(colour_name,col_name[colour]);
          strcpy(colour_numbers,col_numbers[colour]);
        }
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

render_objects  -  Generate script commands for each of the objects
                   described by the rasdef records

***********************************************************************/

void render_objects(FILE *outfile_ptr,struct sequence *sequence_ptr,
                    struct rasdef *first_rasdef_ptr,
                    struct parameters params)
{
  char object[20], object_type[10], range[50];
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char metal_name[4];

  int atom_no, colour, from_atom, ipos, len, ligno, to_atom;
  int first;

  struct atom *atom_ptr;
  struct rasdef *rasdef_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  ligno = 0;

  /* Get pointer to the first RasMol definition record */
  rasdef_ptr = first_rasdef_ptr;

  /* Loop through all the definition records */
  while (rasdef_ptr != NULL)
    {
      /* Check that we have a first and last atom defined */
      if (rasdef_ptr->first_atom_ptr != NULL &&
          rasdef_ptr->last_atom_ptr != NULL)
        {
          /* If this is a protein chain, set colour */
          if (rasdef_ptr->type == PROTEIN || rasdef_ptr->type == CA_ONLY ||
              rasdef_ptr->type == NUCLEIC_ACID)
            {
              /* Get first and last atoms */
              from_atom = rasdef_ptr->first_atom_ptr->out_number;
              to_atom = rasdef_ptr->last_atom_ptr->out_number;

              /* Define atom range */
              sprintf(range,"(atomno >= %d & atomno <= %d)",from_atom,
                      to_atom);

              /* Define object name */
              sprintf(object,"chn%c",rasdef_ptr->chain);

              /* Write out the definition */
              fprintf(outfile_ptr,"\n");
              if (rasdef_ptr->type == PROTEIN)
                {
                  fprintf(outfile_ptr,"# Protein chain %c\n",
                          rasdef_ptr->chain);
                  strcpy(object_type,"protein");
                }
              else if (rasdef_ptr->type == CA_ONLY)
                {
                  fprintf(outfile_ptr,"# CA-only chain %c\n",
                          rasdef_ptr->chain);
                  strcpy(object_type,"CA-only");
                }
              else
                {
                  fprintf(outfile_ptr,"# DNA/RNA chain %c\n",
                          rasdef_ptr->chain);
                  strcpy(object_type,"DNA/RNA");
                }

              /* Define the object using atom-number range */
              fprintf(outfile_ptr,"define %s %s\n",object,range);

              /* Select and colour according to chain colour */
              fprintf(outfile_ptr,"select %s\n",object);

              /* Display according to display option */
              if (params.display_cartoon == TRUE)
                fprintf(outfile_ptr,"cartoon on\n");
              else if (params.display_backbone_only == TRUE ||
                  rasdef_ptr->type == CA_ONLY)
                fprintf(outfile_ptr,"backbone on\n");
              else if (params.display_mainchain_only == TRUE &&
                       params.csa == FALSE)
                fprintf(outfile_ptr,"wireframe on\n");
              else
                fprintf(outfile_ptr,"cartoon on\n");

              /* Get this chain's colour */
              colour = get_chn_colour(rasdef_ptr->chain,colour_name,
                                      colour_numbers);

              /* Define the colour */
              fprintf(outfile_ptr,"colour %s\n",colour_numbers);

              /* Make length of object name a uniform length */
              len = strlen(object);
              for (ipos = 0; ipos < 12 - len; ipos++)
                strcat(object," ");

              /* Show object definition on screen */
              fprintf(outfile_ptr,"echo \"  %s ",object);
              fprintf(outfile_ptr,"= %s chain %c (%s)\"\n",object_type,
                      rasdef_ptr->chain,colour_name);
            }

          /* If this is a ligand, display */
          else if (rasdef_ptr->type == LIGAND)
            {
              /* If ligand hasn't been removed, can display */
              if (rasdef_ptr->first_atom_ptr->display_atom == TRUE)
                {
                  /* Increment ligand number */
                  ligno++;

                  /* Get first and last atoms */
                  from_atom = rasdef_ptr->first_atom_ptr->out_number;
                  to_atom = rasdef_ptr->last_atom_ptr->out_number;

                  /* Define atom range */
                  sprintf(range,"(atomno >= %d & atomno <= %d)",
                          from_atom,to_atom);

                  /* Define object name */
                  sprintf(object,"ligand%d",ligno);
                  fprintf(outfile_ptr,"\n");
                  fprintf(outfile_ptr,"# Ligand %s\n",object);
                  fprintf(outfile_ptr,"define %s %s\n",object,range);

                  /* Select and display */
                  fprintf(outfile_ptr,"select %s\n",object);
                  fprintf(outfile_ptr,"colour cpk\n");
                  fprintf(outfile_ptr,"spacefill 1.2\n");

                  /* Get ligand's first and last residues */
                  residue1_ptr = rasdef_ptr->first_atom_ptr->residue_ptr;
                  residue2_ptr = rasdef_ptr->last_atom_ptr->residue_ptr;

                  /* Get name of ligand's first residue */
                  if (residue1_ptr->chain != ' ')
                    sprintf(range,"%s.%s(%c)",residue1_ptr->res_name,
                            residue1_ptr->res_num,residue1_ptr->chain);
                  else
                    sprintf(range,"%s.%s",residue1_ptr->res_name,
                            residue1_ptr->res_num);
                  remove_blanks(range);

                  /* If last residue is different, add its name */
                  if (residue2_ptr != residue1_ptr)
                    {
                      /* Get name of ligand's last residue */
                      if (residue2_ptr->chain != ' ')
                        sprintf(range,"%s-%s.%s(%c)",range,
                                residue2_ptr->res_name,
                                residue2_ptr->res_num,
                                residue2_ptr->chain);
                      else
                        sprintf(range,"%s-%s.%s",range,
                                residue2_ptr->res_name,
                                residue2_ptr->res_num);
                      remove_blanks(range);
                    }

                  /* Make length of object name a uniform length */
                  len = strlen(object);
                  for (ipos = 0; ipos < 12 - len; ipos++)
                    strcat(object," ");

                  /* Write out this definition */
                  fprintf(outfile_ptr,"echo \"  %s ",object);
                  fprintf(outfile_ptr,"= ligand %s (spacefill)",range);
                  fprintf(outfile_ptr,"\"\n");
                }
            }

          /* If this is a metal, get all atom numbers and display */
          else if (rasdef_ptr->type == METAL)
            {
              /* Define object name */
              sprintf(object,"met_%s",rasdef_ptr->res_num);
              strcpy(metal_name,rasdef_ptr->res_num);
              metal_name[2] = '\0';
              to_lower(metal_name+1,1);
              first = TRUE;

              /* Write heading */
              fprintf(outfile_ptr,"\n");
              fprintf(outfile_ptr,"# Metal %s\n",object);

              /* Get pointer to this structure's first residue */
              residue_ptr = sequence_ptr->first_residue_ptr;

              /* Loop through all the residues to pick up all metals
                 belonging to current rasdef entry */
              while (residue_ptr != NULL)
                {
                  /* If is residue corresponds to the current metal,
                     get atom number */
                  if (residue_ptr->rasdef_ptr == rasdef_ptr)
                    {
                      /* Get the atom number */
                      if (residue_ptr->first_atom_ptr != NULL &&
                          residue_ptr->first_atom_ptr->display_atom == TRUE)
                        {
                          atom_no
                            = residue_ptr->first_atom_ptr->out_number;

                          /* Add atom number to metal definition */
                          if (first == TRUE)
                            fprintf(outfile_ptr,
                                    "define %s atomno = %d\n",
                                    object,atom_no);
                          else
                            fprintf(outfile_ptr,
                                    "define %s %s | atomno = %d\n",
                                    object,object,atom_no);

                          /* Set flag */
                          first = FALSE;
                        }
                    }

                  /* Get pointer to the next residue record */
                  residue_ptr = residue_ptr->next_residue_ptr;
                }

              /* If any definition records written out, display metal
                 in appropriate colour */
              if (first == FALSE)
                {
                  /* Get this metal's colour */
                  strcpy(colour_name,rasdef_ptr->colour);

                  /* Show metal */
                  fprintf(outfile_ptr,"select %s\n",object);
                  fprintf(outfile_ptr,"spacefill 1.5\n");
                  fprintf(outfile_ptr,"colour %s\n",colour_name);

                  /* Make length of object name a uniform length */
                  len = strlen(object);
                  for (ipos = 0; ipos < 12 - len; ipos++)
                    strcat(object," ");

                  /* Write out this definition */
                  fprintf(outfile_ptr,"echo \"  %s ",object);
                  fprintf(outfile_ptr,"= metal ion %s (%s sphere)",
                          metal_name,colour_name);
                  fprintf(outfile_ptr,"\"\n");
                }
            }
        }

      /* Get next definition record */
      rasdef_ptr = rasdef_ptr->next_rasdef_ptr;
    }
}
/***********************************************************************

write_script_csa  -  Write out the script commands to show the given
                     CSA residues

***********************************************************************/

void write_script_csa(FILE *outfile_ptr,struct sequence *sequence_ptr,
                      struct rasdef *first_rasdef_ptr,
                      struct parameters params)
{
  char chain_str[2], object[15], list_name[NAME_LEN + 1], res_num[7];

  int ipos, len, nres;

  struct range *range_ptr;

  /* Initialise variables */
  nres = 0;

  /* Render the different items in the rasdef records */
  render_objects(outfile_ptr,sequence_ptr,first_rasdef_ptr,params);

  /* Form name of the site */
  strcpy(object,"CSA");

  /* Write out the definition */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Definition for catalytic residues (%s)\n",object);

  /* Write out the catalytic residues */
  write_script_reslist(outfile_ptr,params,object);

  /* Make length of object name a uniform length */
  len = strlen(object);
  for (ipos = 0; ipos < 12 - len; ipos++)
    strcat(object," ");

  /* Show the definition on screen */
  fprintf(outfile_ptr,
          "echo \"  %s = catalytic residues (thick bonds, plus labels)\"\n",
          object);
}
/***********************************************************************

read_surface_vertices  -  Read through all the required .vtx and .pNN files
                          to get the total number of vertex files and
                          their minimum and maximum coordinates

***********************************************************************/

int read_surface_vertices(int out_number,int nfiles,char *protein_family,
                          struct limits *limit,
                          struct c_file_data file_name_ptr,
                          struct parameters params)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char root_name[FILENAME_LEN], pdbsum_dir[FILENAME_LEN];
  char chain, chain_id[2], domain, domain_id[3], pdb_code[5];

  char *fname;

  int atom_number;
  int dir_type, transform;

  double grid_spacing, x, x_new, y, y_new, z, z_new;

  double matrix[4][4];

  struct surface *surface_ptr;
  struct vertex *last_vertex_ptr, *vertex_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  atom_number = out_number;
  last_vertex_ptr = NULL;
  strcpy(pdbsum_dir,"./");

  /* Set mask pointer to the first mask */
  surface_ptr = first_surface_ptr;

  /* Loop over all the masks to be written out */
  while (surface_ptr != NULL)
    {
      /* Get the PDB code and chain-id for this mask */
      strncpy(pdb_code,surface_ptr->site_surface_code,4);
      pdb_code[4] = '\0';
      chain = surface_ptr->site_surface_code[4];
      domain = surface_ptr->site_surface_code[5];

      /* Form the name of this PDB entry's PDBsum directory */
      dir_type = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                                 params.pdb_type,pdbsum_dir);

      /* Form the root name of the .mat, .vtx and .pNN files */
      strcpy(root_name,protein_family);

      /* Add the chain-identifier, if there is one */
      if (chain != '0')
        {
          chain_id[0] = chain;
          chain_id[1] = '\0';
          strcat(root_name,chain_id);
        }

      /* Add domain number to file name, if not zero */
      if (domain != '0')
        {
          sprintf(domain_id,"_%c",domain);
          strcat(root_name,domain_id);
        }

      /* If may need to transform the vertex coords, form the filename
         of the .mat file */
      if ((nfiles > 1 || params.nsite_surface_files > 1)
          && protein_family[0] != '\0')
        {
          /* Form the filename */
          strcpy(file_name,pdbsum_dir);
          strcat(file_name,root_name);
          strcat(file_name,".mat");

          /* Pick up the transformation matrix, if there is one */
          transform = get_transformation(file_name,matrix);
        }

      /* Otherwise, no transformation required */
      else
        transform = FALSE;

      /* Form the filename for this mask file */
      strcpy(file_name,pdbsum_dir);
      strcat(file_name,root_name);

      /* If this is an X-SITE distribution, then add .pNN extension */
      if (surface_ptr->xsite_probe[0] != '\0')
        {
          /* If showing the spherical X-SITE representations add .s 
             part of extension */
          if (params.xsite_spheres == TRUE)
            strcat(file_name,".s");

          /* Otherwise, for ordinary X-SITE distribs, add .p part of
             extension */
          else
            strcat(file_name,".p");

          /* Add probe number */
          strcat(file_name,surface_ptr->xsite_probe);
        }

      /* Otherwise, add .vtx extension for the binding-site mask file */
      else
        strcat(file_name,".vtx");

      /* Open the .vtx file */
      if ((file_ptr = fopen(file_name,"r")) !=  NULL)
        {
          /* Re-initialise last-vertex pointer */
          last_vertex_ptr = NULL;

          /* Loop through the .vtx file, picking up all the vertex coords  */
          while (fgets(input_line,LINELEN,file_ptr)!= NULL)
            {
              /* If line starts with a #, check for grid-spacing */
              if (input_line[0] == '#')
                {
                  /* Check for #G line */
                  if (input_line[1] == 'G')
                    {
                      /* Retrieve the grid-spacing */
                      sscanf(input_line+2,"%lf",&grid_spacing);

                      /* Store it for this surface */
                      surface_ptr->grid_spacing = grid_spacing;
                    }
                }

              /* Otherwise, read in the x-, y- and z-coords */
              else
                {
                  /* Get the x-, y- and z-coordinates */
                  sscanf(input_line,"%lf %lf %lf",&x,&y,&z);

                  /* Increment atom number */
                  atom_number++;

                  /* Write out in PDB format */
                  if (atom_number < 100000)
                    {
                      /* If have a transformation, the apply to these
                         coordinates */
                      if (transform == TRUE)
                        {
                          /* Translate */
                          x_new = x - matrix[3][0];
                          y_new = y - matrix[3][1];
                          z_new = z - matrix[3][2];

                          /* Rotate */
                          x = x_new * matrix[0][0]
                            + y_new * matrix[0][1]
                            + z_new * matrix[0][2];
                          y = x_new * matrix[1][0]
                            + y_new * matrix[1][1]
                            + z_new * matrix[1][2];
                          z = x_new * matrix[2][0]
                            + y_new * matrix[2][1]
                            + z_new * matrix[2][2];
                        }

                      /* Allocate memory for structure to hold vertex info */
                      vertex_ptr = (struct vertex *)
                        malloc(sizeof(struct vertex));
                      if (vertex_ptr == NULL)
                        {
                          printf("*** Can't allocate memory for struct");
                          printf(" vertex\n");
                          exit (1);
                        }

                      /* Update the mask links to this vertex */
                      if (surface_ptr->first_vertex_ptr == NULL)
                        surface_ptr->first_vertex_ptr = vertex_ptr;

                      /* Add link from previous vertex to the current one */
                      if (last_vertex_ptr != NULL)
                        last_vertex_ptr->next_vertex_ptr = vertex_ptr;
                      last_vertex_ptr = vertex_ptr;

                      /* Fill in the elements of the structure */
                      vertex_ptr->atom_number = atom_number;
                      vertex_ptr->x = x;
                      vertex_ptr->y = y;
                      vertex_ptr->z = z;
                      vertex_ptr->next_vertex_ptr = NULL;

                      /* Update coordinate limits, if necessary */

                      /* If this is the first coordinate encountered, store it
                         as the present max and min values */
                      if (limit->first == TRUE)
                        {
                          limit->min_x = x;
                          limit->max_x = x;
                          limit->min_y = y;
                          limit->max_y = y;
                          limit->min_z = z;
                          limit->max_z = z;
                          limit->first = FALSE;
                        }

                      /* Otherwise, update the limits if coordinates are
                         outside them */
                      else
                        {
                          if (x < limit->min_x)
                            limit->min_x = x;
                          if (y < limit->min_y)
                            limit->min_y = y;
                          if (z < limit->min_z)
                            limit->min_z = z;
                          if (x > limit->max_x)
                            limit->max_x = x;
                          if (y > limit->max_y)
                            limit->max_y = y;
                          if (z > limit->max_z)
                            limit->max_z = z;
                        }
                    }
                }
            }
        }

      /* If file not found, then display warning message */
      else if (params.write_to_file == TRUE)
        {
          printf("\n*** Warning. Unable to open the .vtx vertices file ");
          printf("[%s]\n",file_name);
        }

      /* Get pointer to next mask record */
      surface_ptr = surface_ptr->next_surface_ptr;
    }

  /* Close mask vertices file */
  fclose(file_ptr);

  /* Return the current atom-number */
  return(atom_number);
}
/***********************************************************************

write_script_mask  -  Write script records to define and display
                      the binding-site mask

***********************************************************************/

void write_script_mask(FILE *outfile_ptr)
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char object[OBJ_NAME_LEN];

  static char *col_name[] = {
    "green", "cyan", "pink", "orange",
    "yellow",  "brown", "greenblue", "magenta", "blue",
    "red", "redorange", "violet", "purple",
    "white"
  };

  static char *col_numbers[] = {
    "[0,255,0]", "[128,255,255]", "[255,128,255]", "[204,128,0]",
    "[255,255,0]", "[128,0,0]", "[51,153,102]", "[255,0,255]", "[0,0,255]",
    "[255,0,0]", "[255,102,0]", "[238,130,238]", "[178,51,255]",
    "[255,255,255]"
  };

  int colour;

  struct surface *surface_ptr;

  /* Set mask pointer to the first mask */
  surface_ptr = first_surface_ptr;

  /* Loop over all the masks to be written out */
  while (surface_ptr != NULL)
    {
      /* Throw line */
      fprintf(outfile_ptr,"\n");

      /* Write comment line */
      fprintf(outfile_ptr,"# Mask %d\n",surface_ptr->surface_number);

      /* Define this surface */
      sprintf(object,"surf%d",surface_ptr->surface_number);
      fprintf(outfile_ptr,"define %s MAP%d:1\n",object,
              surface_ptr->surface_number);

      /* Set mask colour */
      colour = (surface_ptr->surface_number -1 ) % MAX_OBJ_COL;
      strcpy(colour_name,col_name[colour]);
      strcpy(colour_numbers,col_numbers[colour]);

      /* Define the chain colour */
      fprintf(outfile_ptr,"select %s\n",object);
      fprintf(outfile_ptr,"wireframe on\n");
      fprintf(outfile_ptr,"colour bonds %s  # (%s)\n",colour_numbers,
              colour_name);

      /* Save the mask colour */
      strcpy(surface_ptr->colour_name,colour_name);
      strcpy(surface_ptr->colour_numbers,colour_numbers);

      /* Show mask definition */
      fprintf(outfile_ptr,"echo \"  %s ",object);
      if (surface_ptr->xsite_probe[0] == '\0')
        fprintf(outfile_ptr," = binding site mask for ");
      else
        fprintf(outfile_ptr," = X-SITE distribution for ");
      fprintf(outfile_ptr,"%s\"\n",surface_ptr->site_surface_code);

      /* Get pointer to next mask record */
      surface_ptr = surface_ptr->next_surface_ptr;
    }
}
/***********************************************************************

check_for_hydrogen  -  Check whether this is a hydrogen or a pseudo atom
                       from any of the many possible representations of
                       such in the PDB

***********************************************************************/

int check_for_hydrogen(char atname[5])
{
  char atom_name[5];

  int is_hydrogen;

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

  return is_hydrogen;
}
/***********************************************************************

extract_domain_data - Extract the domain data from the free-format line
                      of the CATH domains file and store in a fixed-
                      format array

**********************************************************************/         

int extract_domain_data(char input_line[LINELEN + 1],
                        char domain_details[C_MAXDOMCHAIN][16],
                        int *ndomains,int *nfrags,
                        char error_message[LINELEN + 1],int verbose)
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
      token_len = get_token(input_line,LINELEN,&lpos,token_string,
                            TOKEN_LEN,&done);

      /* If in verbose mode, then print current token */
      if (verbose == TRUE)
        printf("Token [%s]\n",token_string);

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
          sprintf(domain_details[ndetails],"D%2d",idomain + 1);
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

read_in_domain_data  -  Read in the domain data on contacts with ligands
                        for the given structure, if present

***********************************************************************/

int read_in_domain_data(char pdb_code[5],char chain,
                        char cath_domains_file[FILENAME_LEN],
                        char domain_details[C_MAXDOMCHAIN][16],
                        char error_message[LINELEN + 1],
                        int verbose)
{
  char chain_id;
  char cath_code[C_CATH_LEN], wolf_code[7];
  char input_line[LINELEN + 1];

  int line, line_length;
  int ndetails, ndomains, nfrags;
  int got_rec;

  FILE *file_ptr;

  /* Initialise variables */
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
  while (fgets(input_line,LINELEN,file_ptr) != NULL && got_rec == FALSE)
    {
      /* Increment line-coint */
      line++;

      /* Get the line length */
      line_length = strlen(input_line);

      /* Get the chain-id */
      if (line_length > 14)
        {
          chain_id = input_line[4];

          /* Check whether this record applies to the PDB code and chain
             that we want */
          if (!strncmp(input_line,pdb_code,4) &&
              (chain_id == chain ||
               (chain_id == '0' && chain == ' ')))
            {
              /* Set flag to indicate we have found the record we're
                 after */
              got_rec = TRUE;

              /* Extract the domain data from the free-format line and
                 store in a fixed-format array for ease of use */
              /* Old version for old format CATH files
              ndetails = extract_domain_data(input_line,domain_details,
                                             &ndomains,&nfrags,
                                             error_message,verbose);
              */
              ndetails
                = extract_domain_info(input_line,wolf_code,cath_code,
                                      domain_details);

              /* Extract the Wolf code from the front of the line */
              //strncpy(wolf_code,input_line,6);
              //wolf_code[6] = '\0';

              /* If have an error message, display it */
              if (error_message[0] != '\0' && verbose == TRUE)
                {
                  /* Print the error message */
                  printf("*** ERROR in domains file at line %d for ",line);
                  printf("entry [%s]\n",wolf_code);
                  printf("***       %s\n",error_message);
                }
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Show number of entries read */
  if (verbose == TRUE)
    {
      printf("\n");
      printf("Number of entries read: %d\n",line);
      printf("\n");
    }

  /* Return number of detail items found */
  return(ndetails);
}
/***********************************************************************

get_domain  -  Determine which domain the given residue belongs to

***********************************************************************/

void get_domain(char res_num[6],char chain,int *in_range,int *domain_no,
                int *dom_detail,int ndetails,
                char domain_details[C_MAXDOMCHAIN][16])
{
  char number_string[20];

  int domain, domain_detail, idetail;
  int in_ran;

  /* Initialise variables */
  domain = *domain_no;
  domain_detail = *dom_detail;
  in_ran = *in_range;

  /* If are currently in a domain-range, then check for end of range */
  if (in_ran == TRUE)
    {
      /* Check whether this is the end of the residue-range */
      if (!strncmp(domain_details[domain_detail]+10,res_num,5) &&
          (domain_details[domain_detail][9] == chain ||
           (domain_details[domain_detail][9] == '0' && chain == ' ')))
        in_ran = FALSE;
    }

  /* Otherwise, look for domain that this residue might be the start of */
  else
    {
      /* Initialise */
      domain_detail = -1;
      domain = 0;

      /* Loop through the domain detail entries */
      for (idetail = 0; idetail < ndetails && domain_detail == -1; idetail++)
        {
          /* Check whether this is the start of the residue
             range for current domain segment */
          if (!strncmp(domain_details[idetail]+4,res_num,5) &&
              (domain_details[idetail][3] == chain ||
               (domain_details[idetail][3] == '0' && chain == ' ')))
            {
              /* Residue is start of a domain, so save details */
              in_ran = TRUE;
              domain_detail = idetail;

              /* Store the domain number */
              strncpy(number_string,domain_details[idetail]+1,2);
              number_string[2] = '\0';
              domain = atoi(number_string);

              /* Check whether this is also the end of the residue-range
                 (ie a single-residue domain!) */
              if (!strncmp(domain_details[domain_detail]+10,res_num,5) &&
                  (domain_details[domain_detail][9] == chain ||
                   (domain_details[domain_detail][9] == '0' && chain == ' ')))
                in_ran = FALSE;
            }
        }
    }

  /* Return the domain number and flag indicating whether in a residue
     range defined by a domain */
  *in_range = in_ran;
  *domain_no = domain;
  *dom_detail = domain_detail;
}
/***********************************************************************

check_if_in_pattern  -  Check whether the given residue is in one of
                        the stored PROSITE patterns

***********************************************************************/

int check_if_in_pattern(struct prosite_pattern 
                        *first_prosite_pattern_ptr,char chain,char aa_code,
                        char res_name[4],char res_num[6],int *position,
                        int *colour)
{
  int ipos, len;
  int wanted;

  struct prosite_pattern *prosite_pattern_ptr;

  /* Initialise variables */
  *colour = 0;
  ipos = 0;
  wanted = FALSE;

  /* Start at first pattern */
  prosite_pattern_ptr = first_prosite_pattern_ptr;

  /* Loop through all the patterns */
  while (prosite_pattern_ptr != NULL && *colour == 0)
    {
      /* Check whether have the right chain */
      if (prosite_pattern_ptr->chain == chain)
        {
          /* If current residue is the first residue of the pattern,
             then reset position to ensure we are in sync */
          if (!strcmp(prosite_pattern_ptr->first_res_num,res_num) &&
              prosite_pattern_ptr->residue_string[0] == aa_code)
            {
              /* Reset position */
              *position = prosite_pattern_ptr->position;

              /* Have a match at the first position */
              wanted = TRUE;
              ipos = 0;
            }

          /* Otherwise, check whether within the current pattern */
          else
            {
              /* Get length of pattern */
              len = strlen(prosite_pattern_ptr->residue_string);

              /* Get relative position */
              ipos = *position - prosite_pattern_ptr->position;

              /* If within the length, then have a hit if residue-id
                 matches */
              if (ipos >= 0 && ipos < len)
                {
                  if (prosite_pattern_ptr->residue_string[ipos] == aa_code)
                    wanted = TRUE;
                }
            }

          /* If have a match, then retrieve the residue colour */
          if (wanted == TRUE && *colour == 0)
            {
              /* Get length of colour string */
              len = strlen(prosite_pattern_ptr->colours_string);

              /* If position within it is valid, retrieve the colour */
              if (ipos < len)
                {
                  /* Convert character to a colour-number */
                  *colour = prosite_pattern_ptr->colours_string[ipos] - '0';
                }
            }
        }

      /* Go to next pattern */
      prosite_pattern_ptr = prosite_pattern_ptr->next_prosite_pattern_ptr;
    }

  /* Return whether pattern is wanted */
  return(wanted);
}
/***********************************************************************

check_if_in_list  -  Check whether the given residue is in one of the
                     ranges in the residue list ... NOT USED

***********************************************************************/

int check_if_in_list(struct reslist *first_reslist_ptr,char chain,
                     char *res_num,char *range_end)
{
  char res_chain[7];

  int found, in_range;

  struct reslist *reslist_ptr;

  /* Initialise variables */
  in_range = FALSE;
  res_chain[0] = chain;
  strcpy(res_chain+1,res_num);

  /* If have a range open, then check whether this is the end residue */
  if (range_end[0] != '\0')
    {
      /* Mark residue as being within the current range */
      in_range = TRUE;

      /* Check whether this is the end residue */
      if (!strcmp(res_chain,range_end))
        range_end[0] = '\0';
    }

  /* Otherwise, see if this is the first residue in a new range */
  else
    {
      /* Initialise flag */
      found = FALSE;

      /* Start at first residue range */
      reslist_ptr = first_reslist_ptr;

      /* Loop until all the residue ranges have been checked */
      while (reslist_ptr != NULL && found == FALSE)
        {
          /* If this is the start-residue of the range, then set
             flag and store range end residue */
          if (chain == reslist_ptr->chain &&
              !strcmp(res_num,reslist_ptr->res_num1))
            {
              /* Set flag that we're in-range */
              in_range = TRUE;

              /* Store range-end */
              range_end[0] = reslist_ptr->chain;
              strcpy(range_end+1,reslist_ptr->res_num2);

              /* Set flag that range found */
              found = TRUE;
            }

          /* Go to next range */
          reslist_ptr = reslist_ptr->next_reslist_ptr;
        }
    }

  /* Return whether residue is in one of the specified residue ranges */
  return(in_range);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

void create_residue_record(struct residue **new_residue_ptr,
                           struct sequence *sequence_ptr,char res_name[4],
                           char res_num[6],char chain,int in_list,
                           int from_alignment,int colour,int is_ligand,
                           int interacting_residue,int type,
                           struct rasdef *rasdef_ptr)
{
  struct residue *residue_ptr;

  /* Allocate memory for structure to hold residue info */
  residue_ptr = (struct residue *) malloc(sizeof(struct residue));
  if (residue_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct");
      printf(" residue\n");
      exit (1);
    }

  /* Update the sequence links to this residue */
  if (sequence_ptr->first_residue_ptr == NULL)
    sequence_ptr->first_residue_ptr = residue_ptr;

  /* Add link from previous residue to the current one */
  if (sequence_ptr->last_residue_ptr != NULL)
    sequence_ptr->last_residue_ptr->next_residue_ptr = residue_ptr;
  sequence_ptr->last_residue_ptr = residue_ptr;

  /* Initialise the elements of the structure */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->natoms = 0;

  /* Fill in currently-known details */
  residue_ptr->type = type;
  residue_ptr->in_list = in_list;
  residue_ptr->from_alignment = from_alignment;
  residue_ptr->colour = colour;
  residue_ptr->is_ligand = is_ligand;
  residue_ptr->interacting_residue = interacting_residue;
  residue_ptr->check_distance = FALSE;
  residue_ptr->rasdef_ptr = rasdef_ptr;

  residue_ptr->first_atom_no = 0;
  residue_ptr->last_atom_no = 0;
  residue_ptr->position = 0;
  residue_ptr->rmsd = 0;
  residue_ptr->fit = FALSE;
  residue_ptr->fit_used = FALSE;
  residue_ptr->first_atom_ptr = NULL;
  residue_ptr->N_atom_ptr = NULL;
  residue_ptr->CA_atom_ptr = NULL;
  residue_ptr->C_atom_ptr = NULL;
  residue_ptr->O_atom_ptr = NULL;
  residue_ptr->reference_residue_ptr = NULL;
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointer */
  *new_residue_ptr = residue_ptr;
}
/***********************************************************************

check_if_want_residue  -  Check whether this residue is required

***********************************************************************/

struct residue *check_if_want_residue(struct sequence *sequence_ptr,
                                      char *res_num,char *res_name,
                                      char chain,int ndetails,
                                      char domain_details[C_MAXDOMCHAIN][16],
                                      int *in_rng,int *dom_no,
                                      int *dom_det,int *position,
                                      char at_het,struct rasdef *rasdef_ptr,
                                      struct parameters params)
{
  char aa_code;

  static char range_end[6];

  int colour, domain_no, dom_detail, ichain, pos;
  int from_alignment, in_list, in_range, type, wanted;
  int fit, interacting_residue, is_ligand;

  static int first = TRUE;

  struct list *list_ptr, *save_list_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  in_list = NOT_IN_LIST;
  from_alignment = FALSE;
  is_ligand = FALSE;
  interacting_residue = FALSE;
  if (first == TRUE)
    {
      range_end[0] = '\0';
      first = FALSE;
    }
  if (rasdef_ptr == NULL)
    type = UNKNOWN;
  else
    type = rasdef_ptr-> type;

  /* Store returnable variables */
  fit = FALSE;
  in_range = *in_rng;
  domain_no = *dom_no;
  dom_detail = *dom_det;
  pos = *position;

  /* Assume residue is not wanted */
  colour = 0;
  residue_ptr = NULL;
  save_list_ptr = NULL;
  wanted = FALSE;

  /* Get this residue's single-letter code */
  aa_code = get_aa_code(res_name);

  /* Initialise search pointer to first residue in the linked list */
  list_ptr = sequence_ptr->first_list_ptr;
  
  /* Loop through all the stored residues to check for a match */
  while (list_ptr != NULL && wanted == FALSE)
    {
      /* Check whether this is the correct residue */
      if (!strncmp(res_num,list_ptr->res_num,5) &&
          !strncmp(res_name,list_ptr->res_name,3) &&
          (chain == list_ptr->chain ||
           (chain == 'A' && list_ptr->chain == ' ')))
        {
          /* Have a match, so residue is wanted */
          wanted = TRUE;

          /* Indicate that came from list */
          in_list = IN_LIST;
          
          /* Save colour */
          colour = list_ptr->colour;
          
          /* If from a CORA-format alignment file, then indicate as such */
          if (params.cath_alignment == TRUE)
            {
              /* Set flag */
              from_alignment = TRUE;
              
              /* Only mark as coming from a list is residue is flagged */
              if (colour == 0)
                in_list = NOT_IN_LIST;
            }
          
          /* Save whether ligand or interacting residue */
          is_ligand = list_ptr->is_ligand;
          interacting_residue = list_ptr->interacting_residue;

          /* Store whether a fitting residue */
          fit = list_ptr->fit;

          /* Store pointer from list item to residue created from it */
          save_list_ptr = list_ptr;
        }

      /* Get pointer to the next residue in the linked list */
      list_ptr = list_ptr->next_list_ptr;
    }

  /* If displaying PROSITE patterns, check if this residue belongs to
     one of those required */
  if (params.prosite_pattern == TRUE)
    {
      /* Check if this residue belongs to a PROSITE pattern */
      wanted = check_if_in_pattern(sequence_ptr->first_prosite_pattern_ptr,
                                   chain,aa_code,res_name,res_num,&pos,
                                   &colour);

      /* Indicate that came from list */
      in_list = IN_LIST;
    }

  /* If have a residue list, check whether this residue belongs */
  // if (params.residue_list[0] != '\0')
  //  {
  //    /* Determine whether residue is to be shown */
  //    in_list = check_if_in_list(params.first_reslist_ptr,chain,
  //                               res_num,range_end);

  //    /* If residue is in the list, mark as wanted */
  //    if (in_list == TRUE)
  //      wanted = TRUE;
  //  }

  /* If residue not selected, check whether all residues from this
     chain/domain required anyway */
  if (wanted == FALSE)
    {

      /* If have a residue list, accept all residues for now */
      if (params.residue_list[0] != '\0')
        wanted = TRUE;

      /* For PROMOTIF motifs, active sites or CSA catalytic residues
         retain all required chains */
      else if (params.promotif == TRUE || params.pdb_site == TRUE ||
          params.csa == TRUE)
        {
          /* Check whether chain is wanted */
          for (ichain = 0; ichain < params.nchains; ichain++)
            if (chain == params.chains[ichain])
              wanted = TRUE;

          /* Check whether this is a ligand or metal */
          if ((params.pdb_site == TRUE || params.csa == TRUE) &&
              (type == LIGAND || type == METAL))
            {
              /* Mark as ligand */
              wanted = TRUE;
              is_ligand = TRUE;
            }
        }

      /* Retain HET groups for case of PROSITE patterns */
      else if (params.prosite_pattern == TRUE &&
               (at_het == 'H' || type == LIGAND || type == METAL))
        {
          /* Check whether wanted from chain identifier */
          if (params.all_chains == TRUE ||
              sequence_ptr->chain == ' ' || sequence_ptr->chain == chain)
            {
              /* Mark as ligand */
              wanted = TRUE;
              is_ligand = TRUE;
            }
        }

      /* Keep HET groups if required for SAS alignments */
      else if ((at_het == 'H' && type != PROTEIN) || type == LIGAND ||
               type == METAL)
        {
          if (params.include_hets == TRUE)
            {
              wanted = TRUE;
              is_ligand = TRUE;
            }
          else
            wanted = FALSE;
        }

      /* If all chains required, regardless, then flag as wanted */
      else if (params.all_chains == TRUE)
        wanted = TRUE;

      /* If all of the given chain is required, then select */
      else if (params.single_chain == TRUE)
        {
          if (chain == ' ' || chain == sequence_ptr->chain)
            wanted = TRUE;
          else
            wanted = FALSE;
       }

      /* If a single domain, check whether this residue belongs to it */
      else if (params.single_domain == TRUE)
        {
          /* If this is a single-domain chain, then treat as single-chain
             case */
          if (sequence_ptr->domain_no == 0)
            {
              /* If chain matches, then want this residue */
              if (chain == ' ' || chain == sequence_ptr->chain)
                wanted = TRUE;
           }

          /* Otherwise, get this residue's domain */
          else
            {
              /* Get the domain number */
              get_domain(res_num,chain,&in_range,&domain_no,&dom_detail,
                         ndetails,domain_details);

              /* Check whether we have the right domain */
              if (domain_no == sequence_ptr->domain_no)
                wanted = TRUE;
            }
        }

      /* If this is the correct chain, then it is always wanted */
      else if (params.whole_protein == TRUE &&
               (chain == ' ' || chain == sequence_ptr->chain))
        wanted = TRUE;

      /* If this is a water, then don't want it */
      if (!strncmp(res_name,"HOH",3))
        wanted = FALSE;
    }

  /* If residue is wanted, then create residue record and store all
     the details */
  if (wanted == TRUE)
    {
      /* Have a match, so create residue */
      create_residue_record(&residue_ptr,sequence_ptr,res_name,res_num,
                            chain,in_list,from_alignment,colour,
                            is_ligand,interacting_residue,type,
                            rasdef_ptr);
      /* Store whether a fitting residue */
      residue_ptr->fit = fit;

      /* If residue came from a list, then update the pointers */
      if (save_list_ptr != NULL)
        {
          /* Store pointer from list item to residue created from it */
          save_list_ptr->residue_ptr = residue_ptr;

          /* If the list item referred to a reference residue, track down
             its pointer and store */
          if (save_list_ptr->reference_list_ptr != NULL)
            residue_ptr->reference_residue_ptr
              = save_list_ptr->reference_list_ptr->residue_ptr;
        }
    }

  /* Return returnable variables */
  *in_rng = in_range;
  *dom_no = domain_no;
  *dom_det = dom_detail;
  *position = pos;

  /* Return the residue pointer */
  return(residue_ptr);
}
/***********************************************************************

check_ranges  -  Check whether residue is within the given range

***********************************************************************/

int check_ranges(struct residue *residue_ptr,
                 struct range *first_range_ptr)
{
  char number[6];

  int ires, res_num, range_num[2]; 
  int in_motif;

  struct range *range_ptr;

  /* Initialise variables */
  in_motif = FALSE;
  strcpy(number,residue_ptr->res_num);
  number[4] = '\0';
  res_num = atoi(number);

  /* Initialise pointer to the first record */
  range_ptr = first_range_ptr;

  /* Loop through all the range records */
  while (range_ptr != NULL && in_motif == FALSE)
    {
      /* Check chain matches */
      if (residue_ptr->chain == range_ptr->chain)
        {
          /* Get the start and end residue numbers for this range */
          for (ires = 0; ires < 2; ires++)
            {
              /* Extract the number */
              strcpy(number,range_ptr->res_num[ires]);
              number[4] = '\0';
              range_num[ires] = atoi(number);
            }

          /* If our residue lies within this range it is in the motif */
          if (res_num >= range_num[0] && res_num <= range_num[1])
            in_motif = TRUE;
        }

      /* Get pointer to the next range record */
      range_ptr = range_ptr->next_range_ptr;
    }

  /* Return whether this residue is in a motif */
  return(in_motif);
}
/***********************************************************************

check_reslist  -  Check whether residue corresponds to a start or end
                  residue for any reslist record

***********************************************************************/

void check_reslist(struct residue *residue_ptr,
                   struct reslist *first_reslist_ptr)
{
  char chain, res_num[6];

  int ires, reslist_num[2]; 

  struct reslist *reslist_ptr;

  /* Initialise variables */
  chain = residue_ptr->chain;
  strcpy(res_num,residue_ptr->res_num);

  /* Initialise pointer to the first record */
  reslist_ptr = first_reslist_ptr;

  /* Loop through all the reslist records */
  while (reslist_ptr != NULL)
    {
      /* Check chain matches */
      if (chain == reslist_ptr->chain)
        {
          /* Check the start residue */
          if (!strcmp(reslist_ptr->res_num1,res_num))
            {
              /* Save residue pointer as first residue */
              reslist_ptr->residue1_ptr = residue_ptr;
              residue_ptr->in_list = IN_LIST;
            }

          /* Repeat for end residue */
          if (!strcmp(reslist_ptr->res_num2,res_num))
            {
              /* Save residue pointer as end residue */
              reslist_ptr->residue2_ptr = residue_ptr;
              residue_ptr->in_list = IN_LIST;
            }
        }

      /* Get pointer to the next reslist record */
      reslist_ptr = reslist_ptr->next_reslist_ptr;
    }
}
/***********************************************************************

create_atom_record  -  Create and initialise a new atom record

***********************************************************************/

void create_atom_record(struct atom **new_atom_ptr,int atom_number,
                        char atom_name[5],int out_number,char at_het,
                        int alt_occup,double x,double y,double z,
                        char occup_bvalue[13],int display_atom,int fit_atom,
                        int check_distance,struct residue *residue_ptr)
{
  struct atom *atom_ptr;
  struct rasdef *rasdef_ptr;

  /* Allocate memory for structure to hold atom info */
  atom_ptr = (struct atom *) malloc(sizeof(struct atom));
  if (atom_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct atom\n");
      exit (1);
    }

  /* Store given details */
  atom_ptr->atom_number = atom_number;
  strncpy(atom_ptr->atom_name,atom_name,4);
  atom_ptr->atom_name[4] = '\0';
  atom_ptr->residue_ptr = residue_ptr;
  atom_ptr->out_number = out_number;

  /* Store whether atom comes from ATOM or HETATM record */
  atom_ptr->at_het = at_het;

  /* Get the atom's coordinates and store */
  atom_ptr->x = x;
  atom_ptr->y = y;
  atom_ptr->z = z;

  /* Get the occupancy and B-value */
  strcpy(atom_ptr->occup_bvalue,occup_bvalue);

  /* Store whether this atom to be displayed */
  atom_ptr->display_atom = display_atom;

  /* Store whether this is a fitting atom */
  atom_ptr->fit_atom = fit_atom;

  /* Store whether distance-check required */
  if (check_distance == TRUE)
    residue_ptr->check_distance = TRUE;

  /* Update residue's record of which atoms it holds */
  if (residue_ptr->natoms == 0)
    {
      residue_ptr->first_atom_no = out_number;
      residue_ptr->first_atom_ptr = atom_ptr;
    }
  residue_ptr->last_atom_no = out_number;

  /* Increment count of this residue's atoms */
  residue_ptr->natoms++;

  /* Update corresponding rasdef record, if there is one */
  if (residue_ptr->rasdef_ptr != NULL)
    {
      /* Get the rasdef record */
      rasdef_ptr = residue_ptr->rasdef_ptr;

      /* If this is the first atom, store */
      if (rasdef_ptr->first_atom_ptr == NULL)
        rasdef_ptr->first_atom_ptr = atom_ptr;

      /* Store as last atom */
      rasdef_ptr->last_atom_ptr = atom_ptr;
    }

  /* Initialise pointer to next atom */
  atom_ptr->next_atom_ptr = NULL;

  /* Initialise details */
  atom_ptr->written = FALSE;

  /* Return current pointer */
  *new_atom_ptr = atom_ptr;
}
/***********************************************************************

get_atom_details  -  Check whether this atom is required

***********************************************************************/

void get_atom_details(char *input_line,char *atom_name,
                      int *atom_number,char *at_het,int *alt_occup,
                      char *occup_bvalue,double *x,double *y,double *z)
{
  char occup[4];
  char number_string[20];

  int len;

  /* Get the length of the input line */
  len = strlen(input_line);

  /* Get the atom type */
  strncpy(atom_name,input_line+12,4);
  atom_name[4] = '\0';

  /* Get the atom-number */
  strncpy(number_string,input_line+6,5);
  number_string[5] = '\0';
  *atom_number = atoi(number_string);

  /* Get whether atom comes from ATOM of HETATM record */
  *at_het = input_line[0];

  /* Determine whether atom is an alternate occupancy atom */
  *alt_occup = FALSE;
  strncpy(occup,input_line+56,3);
  occup[3] = '\0';
  if (input_line[16] != ' ' || strncmp(occup,"1.0",3))
    *alt_occup = TRUE;

  /* Retrieve the atomic coordinates, occupancy and B-value */
  strncpy(number_string,input_line+30,8);
  number_string[8] = '\0';
  *x = atof (number_string);
  strncpy(number_string,input_line+38,8);
  number_string[8] = '\0';
  *y = atof (number_string);
  strncpy(number_string,input_line+46,8);
  number_string[8] = '\0';
  *z = atof (number_string);

  /* Get the occupancy and B-value */
  occup_bvalue[0] = '\0';
  if (len > 65)
    {
      strncpy(occup_bvalue,input_line+54,12);
      occup_bvalue[12] = '\0';
    }
}
/***********************************************************************

assign_to_object  -  Check which of the elements in the rasmol.dfn file
                     the given residue should be assigned to

***********************************************************************/

struct rasdef *assign_to_object(struct rasdef *first_rasdef_ptr,
                                char *atom_name,char *res_name,
                                char *res_num,char chain)
{
  int type;
  int have_exact;

  struct rasdef *rasdef_ptr, *assigned_rasdef_ptr;

  /* Initialise variables */
  assigned_rasdef_ptr = NULL;
  have_exact = FALSE;
  type = UNKNOWN;

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

  /* Return the matched rasdef pointer */
  return(assigned_rasdef_ptr);
}
/***********************************************************************

check_if_want_atom  -  Check whether this atom is required

***********************************************************************/

struct atom *check_if_want_atom(char atom_name[5],int atom_number,
                                char at_het,int alt_occup,
                                double x,double y,double z,
                                char occup_bvalue[13],
                                struct residue *residue_ptr,
                                struct atom *first_atom_ptr,
                                int display_option,int *out_number,
                                struct atom **lst_atom_ptr,
                                struct sequence *sequence_ptr,
                                int *atom_counter,int in_motif,
                                struct parameters params)
{
  int fit_atom;
  int check_distance, display_atom, hydrogen, wanted;

  struct atom *atom_ptr, *check_atom_ptr, *last_atom_ptr;
  struct residue *res_ptr;

  /* Initialise variables */
  check_distance = FALSE;

  /* Assume atom is not wanted */
  atom_ptr = NULL;

  /* Get the pointer to the last atom */
  last_atom_ptr = *lst_atom_ptr;

  /* Check whether it is a hydrogen, in which case we don't want
     it */
  hydrogen = check_for_hydrogen(atom_name);
  if (hydrogen == TRUE)
    wanted = FALSE;
  else
    wanted = TRUE;

  /* If still wanted, make further checks */
  if (wanted == TRUE)
    {
      /* Set flag to store this atom */
      fit_atom = -1;
      display_atom = FALSE;
      wanted = FALSE;

      /* If this is a fitting residue only, then store only if it
         is a fitting atom */
      if (residue_ptr->from_alignment == TRUE)
        {
          /* Assume atom is wanted */
          wanted = TRUE;

          /* If fitting to backbone, only want the atom if it
             is a CA-alpha */
          if (params.fit_to_backbone == TRUE)
            {
              if (!strncmp(atom_name," CA ",4))
                fit_atom = 0;
              else
                wanted = FALSE;
            }

          /* If fitting to mainchain, then store if this is
             a mainchain atom */
          else if (params.fit_to_mainchain == TRUE)
            {
              if (!strncmp(atom_name," CA ",4))
                fit_atom = 0;
              else if (!strncmp(atom_name," N  ",4))
                fit_atom = 1;
              else if (!strncmp(atom_name," C  ",4))
                fit_atom = 2;
              else if (!strncmp(atom_name," O  ",4))
                fit_atom = 3;
              else
                wanted = FALSE;
            }
        }

      /* If all atoms to be displayed, then include */
      if (display_option == ALL_ATOMS)
        display_atom = TRUE;

      /* If a HET atom, then set for distance-check */
      if (params.include_hets == TRUE &&
          ((at_het == 'H' && residue_ptr->type != PROTEIN) ||
           residue_ptr->type == LIGAND ||
           residue_ptr->type == METAL) &&
          strncmp(residue_ptr->res_name,"HOH",3))
        {
          /* Mark atom as requiring a distance-check */
          display_atom = TRUE;
          check_distance = TRUE;
        }

      /* If this is a ligand residue, the check for ligand-display
         options */
      if (residue_ptr->is_ligand == TRUE)
        {
          /* Check for one of the options that displays the ligand */
          if (display_option == LIGAND_ONLY ||
              display_option == LIGAND_AND_RES)
            display_atom = TRUE;
        }

      /* Otherwise check the various other possibilities */
      else
        {
          /* Check for display options involving the whole
             residue */
          if (display_option == LIGAND_AND_RES ||
              display_option == BINDING_SITE ||
              display_option == INT_RESIDUES)
            display_atom = TRUE;

          /* Otherwise, check that atom type corresponds to the
             required display mode */
          else
            {
              /* For CSA catalytic residues, want every atom of the
                 residue */
              if (params.csa == TRUE && residue_ptr->in_list == IN_LIST)
                display_atom = TRUE;

              /* Main-chain or cartoon */
              else if (display_option == MAIN_CHAIN ||
                       display_option == CARTOON)
                {
                  /* Check whether this is a mainchain atom */
                  if (!strncmp(atom_name," CA ",4) ||
                      !strncmp(atom_name," N  ",4) ||
                      !strncmp(atom_name," C  ",4) ||
                      !strncmp(atom_name," O  ",4))
                    display_atom = TRUE;
                }

              /* C-alpha option */
              else if (display_option == C_ALPHA)
                {
                  /* Check whether this is a C-alpha atom */
                  if (!strncmp(atom_name," CA ",4))
                    display_atom = TRUE;
                }
            }
        }

      /* Display options for PROSITE patterns, or SAS-coloured structures
         with coloured sidechains shown */
      if (params.prosite_pattern == TRUE || params.side_sas == TRUE)
        {
          /* If this is one of the coloured PROSITE pattern residues,
             then display whole residue */
          if (params.prosite_pattern == TRUE && residue_ptr->colour > 0)
            display_atom = TRUE;

          /* If this is one of the coloured SAS residues, then if wanting
             to show the whole sidechain, grab any sidechain atoms */
          else if (params.side_sas == TRUE && residue_ptr->colour > 0)
            {
              /* If not a main-chain atom, grab it */
              if (strncmp(atom_name," CA ",4) &&
                  strncmp(atom_name," N  ",4) &&
                  strncmp(atom_name," C  ",4) &&
                  strncmp(atom_name," O  ",4))
                display_atom = TRUE;
            }

          /* Otherwise, if a HET group, then save for later checking */
          else if ((at_het == 'H' || residue_ptr->type == LIGAND ||
                    residue_ptr->type == METAL)
                   && strncmp(residue_ptr->res_name,"HOH",3))
            {
              /* Mark atom as requiring a distance-check */
              display_atom = TRUE;
              check_distance = TRUE;
            }

          /* Otherwise, for PROSITE patterns, save only mainchain atoms */
          else if (params.prosite_pattern == TRUE)
            {
              if (!strncmp(atom_name," CA ",4) ||
                  !strncmp(atom_name," N  ",4) ||
                  !strncmp(atom_name," C  ",4) ||
                  !strncmp(atom_name," O  ",4))
                display_atom = TRUE;
            }
        }

      /* If showing residue lists, get all atoms for the in-list
         residues */
      if (params.residue_list[0] != '\0')
        {
          /* If residue is in the list, then want this atom */
          if (residue_ptr->in_list != NOT_IN_LIST)
            {
              /* Main-chain or cartoon */
              if (display_option == MAIN_CHAIN ||
                  display_option == CARTOON)
                {
                  /* Check whether this is a mainchain atom */
                  if (!strncmp(atom_name," CA ",4) ||
                      !strncmp(atom_name," N  ",4) ||
                      !strncmp(atom_name," C  ",4) ||
                      !strncmp(atom_name," O  ",4))
                    display_atom = TRUE;
                }

              /* C-alpha option */
              else if (display_option == C_ALPHA)
                {
                  /* Check whether this is a C-alpha atom */
                  if (!strncmp(atom_name," CA ",4))
                    display_atom = TRUE;
                }

              /* Otherwise, assume all atoms required */
              else
                display_atom = TRUE;
            }
        }

      /* If this is a water, then display only if we are displaying
         binding-site waters and this is one of them */
      if (!strncmp(residue_ptr->res_name,"HOH",3))
        {
          if (params.pdb_site == TRUE ||
              (residue_ptr->is_ligand == TRUE &&
               params.display_binding_site_waters == TRUE))
            display_atom = TRUE;
          else
            {
              display_atom = FALSE;
              wanted = FALSE;
            }
        }

      /* If displaying a PROMOTIF motif, test whether current residue
         belongs to one of the motifs */
      if (params.promotif == TRUE && in_motif == TRUE)
        display_atom = TRUE;

      /* If residue belongs to a SITE residue, then want this atom */
      if (params.pdb_site == TRUE && residue_ptr->colour > 0)
        display_atom = TRUE;

      /* If atom to be displayed, then want to store it */
      if (display_atom == TRUE)
        wanted = TRUE;

      /* If this is an alternate occupancy atom check whether we already
         have it */
      if (alt_occup == TRUE)
        {
          /* Initialise pointer to first atom */
          check_atom_ptr = first_atom_ptr;
          
          /* Loop through all the atoms to check whether we
             already have this one */
          while (check_atom_ptr != NULL && wanted == TRUE)
            {
              /* Get this atom's residue pointer */
              res_ptr = check_atom_ptr->residue_ptr;

              /* If we already have this atom, then don't want
                 the second copy */
              if (!strncmp(res_ptr->res_name,residue_ptr->res_name,3) &&
                  !strncmp(res_ptr->res_num,residue_ptr->res_num,5) &&
                  res_ptr->chain == residue_ptr->chain &&
                  !strncmp(check_atom_ptr->atom_name,atom_name,4))
                wanted = FALSE;

              /* Get pointer to next atom in linked-list */
              check_atom_ptr = check_atom_ptr->next_atom_ptr;
            }
        }

      /* If atom wanted, then store its details */
      if (wanted == TRUE)
        {
          /* Increment count of atoms written out and store */
          (*out_number)++;

          /* Create a new atom record */
          create_atom_record(&atom_ptr,atom_number,atom_name,
                             *out_number,at_het,alt_occup,x,y,z,
                             occup_bvalue,display_atom,fit_atom,
                             check_distance,residue_ptr);

          /* Add link from previous atom to the current one */
          if (last_atom_ptr != NULL)
            last_atom_ptr->next_atom_ptr = atom_ptr;
          else
            sequence_ptr->first_atom_ptr = atom_ptr;
          last_atom_ptr = atom_ptr;

          /* If it is a fitting atom, then store a pointer to it from
             the current residue */
          if (fit_atom == 0)
            residue_ptr->CA_atom_ptr = atom_ptr;
          else if (fit_atom == 1)
            residue_ptr->N_atom_ptr = atom_ptr;
          else if (fit_atom == 2)
            residue_ptr->C_atom_ptr = atom_ptr;
          else if (fit_atom == 3)
            residue_ptr->O_atom_ptr = atom_ptr;

          /* Increment count of atoms read in */
          (*atom_counter)++;
        }
    }

  /* Return atom pointer */
  *lst_atom_ptr = last_atom_ptr;
  return(atom_ptr);
}
/***********************************************************************

get_atom_records  -  Read through the PDB file, writing out all the
                     atom and CONECT records for the residue-list read
                     in from ligplot.res

***********************************************************************/

int get_atom_records(struct sequence *sequence_ptr,char *pdb_name,
                     int *next_atom,struct parameters params,
                     char cath_domains_file[FILENAME_LEN],
                     struct range *first_range_ptr,
                     struct rasdef *first_rasdef_ptr,
                     FILE *outfile_ptr)
{
  char input_line[LINELEN + 1];
  char at_het, atom_name[5];
  char res_num[6], res_name[4], last_res_name[4], last_res_num[6];
  char occup_bvalue[13];
  char chain, last_chain;
  char domain_details[C_MAXDOMCHAIN][16], error_message[LINELEN + 1];

  char *fname;

  int keep_reading, wanted;
  int atom_counter, atom_number, len;
  int domain_no, dom_detail, ndetails;
  int out_number, position, type;
  int alt_occup, eof, gzipped, verbose, domain_OK, in_motif, in_range;

  double x, y, z;

  struct atom *atom_ptr, *last_atom_ptr;
  struct rasdef *rasdef_ptr;
  struct residue *residue_ptr;
    
#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  atom_counter = 0;
  verbose = FALSE;
  domain_OK = TRUE;
  dom_detail = 0;
  eof = FALSE;
  keep_reading = TRUE;
  in_motif = FALSE;
  in_range = FALSE;
  last_atom_ptr = NULL;
  last_res_name[0] = '\0';
  last_res_num[0] = '\0';
  last_chain = '\0';
  ndetails = 0;
  position = 0;
  residue_ptr = NULL;
  type = UNKNOWN;

  /* Initialise pointers */
  residue_ptr = NULL;

  /* Get next atom number */
  out_number = *next_atom;

  /* Read in the domain details for this PDB file, if any, and if required */
  if (sequence_ptr->domain_no != 0)
    {
      ndetails = read_in_domain_data(sequence_ptr->pdb_code,
                                     sequence_ptr->chain,
                                     cath_domains_file,domain_details,
                                     error_message,verbose);
    }

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
    }
  else
    eof = TRUE;

  /* If file found, then read in the data */
  if (eof == TRUE)
    {
      fprintf(outfile_ptr,"echo \"*** ERROR. Unable to open PDB file ");
      fprintf(outfile_ptr,"[%s] \"\n",pdb_name);
      fprintf(outfile_ptr,"echo \" \"\n");
      return(atom_counter);
    }

  /* Initialise flag */
  keep_reading = TRUE;

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
          /* Initialise flag */
          wanted = FALSE;

          /* Check for the end of an NMR structure model */
          if (!strncmp("ENDMDL",input_line,6))
            keep_reading = FALSE;

          /* If this is an ATOM or HETATM record, process it */
          else if (keep_reading == TRUE &&
                   (!strncmp(input_line,"ATOM",4) ||
                    !strncmp(input_line,"HETATM",6) ||
                    !strncmp(input_line,"ATZERO",6) ||
                    !strncmp(input_line,"HEZERO",6)))
            {
              /* Truncate the string to strip off the newline */
              string_truncate(input_line,LINELEN);

              /* Get the residue name, number, chain-ID and atom type */
              strncpy(res_name,input_line+17,3);
              res_name[3] = '\0';
              strncpy(res_num,input_line+22,5);
              res_num[5] = '\0';
              chain = input_line[21];

              /* Determine whether ATOM or HETATM record */
              at_het = input_line[0];

              /* Retrieve the atom details from the input line */
              get_atom_details(input_line,atom_name,&atom_number,&at_het,
                               &alt_occup,occup_bvalue,&x,&y,&z);

              /* Check whether atom belongs to a new residue */
              if (chain != last_chain || strncmp(res_num,last_res_num,5) ||
                  strncmp(res_name,last_res_name,3))
                {
                  /* If this is a new chain, re-initialise residue-position
                     counter */
                  if (chain != last_chain)
                    position = 0;

                  /* Increment residue position */
                  position++;

                  /* Determine residue type using the definitions read in
                     from the rasmol.dfn file */
                  rasdef_ptr
                    = assign_to_object(first_rasdef_ptr,atom_name,res_name,
                                       res_num,chain);

                  /* Determine whether this residue is required */
                  residue_ptr
                    = check_if_want_residue(sequence_ptr,res_num,res_name,
                                            chain,ndetails,domain_details,
                                            &in_range,&domain_no,&dom_detail,
                                            &position,at_het,rasdef_ptr,
                                            params);

                  /* If showing PROMOTIF motifs, check whether this
                     residue belongs to one of the motifs to be displayed */
                  if (residue_ptr != NULL && params.promotif == TRUE)
                    {
                      /* Check the motif residue ranges */
                      in_motif = check_ranges(residue_ptr,first_range_ptr);
                    }

                  /* If have a residue list, check whether this residue is
                     one of the start or end residues */
                  if (residue_ptr != NULL &&
                      params.first_reslist_ptr != NULL)
                    check_reslist(residue_ptr,params.first_reslist_ptr);

                  /* Save the current residue number and chain */
                  strncpy(last_res_name,res_name,3);
                  last_res_name[3] = '\0';
                  strncpy(last_res_num,res_num,5);
                  last_res_num[5] = '\0';
                  last_chain = chain;
                }

              /* If have a valid, wanted residue, then process this atom */
              if (residue_ptr != NULL)
                {
                  /* Check whether atom is wanted */
                  atom_ptr = check_if_want_atom(atom_name,atom_number,at_het,
                                                alt_occup,x,y,z,occup_bvalue,
                                                residue_ptr,
                                                sequence_ptr->first_atom_ptr,
                                                sequence_ptr->display_option,
                                                &out_number,&last_atom_ptr,
                                                sequence_ptr,&atom_counter,
                                                in_motif,params);
                }
            }

          /* If this is a CONECT record, then see if it is wanted */
          else if (!strncmp(input_line,"CONECT",6))
            process_conect_record(sequence_ptr,input_line);
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

  /* Save the last-used atom-number */
  *next_atom = out_number;

  /* Return the number of atoms read in */
  return(atom_counter);
}
/***********************************************************************

check_distance  -  Reject all HET groups if they're too far from anything

***********************************************************************/

void check_distance(int ifile,int *nligands)
{
  int iatom, natoms;
  int display;

  double dist_2, x, x1, y, y1, z, z1;

  struct atom *atom_ptr, *other_atom_ptr;
  struct residue *residue_ptr, *other_residue_ptr;
  struct sequence *sequence_ptr;

  /* Get the sequence pointer for the current PDB entry */
  sequence_ptr = stored_sequence_ptr[ifile];

  /* Initialise pointer to first residue */
  residue_ptr = sequence_ptr->first_residue_ptr;

  /* Loop through all the residues to calculate the minimum
     distance of this one from those already stored */
  while (residue_ptr != NULL)
    {
      /* If this residue requires a distance-check, then continue */
      if (residue_ptr->check_distance == TRUE)
        {
          /* Assume residue to be switched off */
          display = FALSE;

          /* Get this residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;
          iatom = 0;

          /* Loop through the residue's atoms */
          while (atom_ptr != NULL && iatom < natoms && display == FALSE)
            {
              /* Get atom's coordinates */
              x1 = atom_ptr->x;
              y1 = atom_ptr->y;
              z1 = atom_ptr->z;

              /* Initialise other pointer to first atom */
              other_atom_ptr = sequence_ptr->first_atom_ptr;

              /* Loop through all atoms to check distances */
              while (other_atom_ptr != NULL && display == FALSE)
                {
                  /* Get this atom's residue */
                  other_residue_ptr = other_atom_ptr->residue_ptr;

                  /* Consider only if not a distance-check atom */
                  if (other_residue_ptr->check_distance == FALSE)
                    {
                      /* Get this atom's coordinates */
                      x = other_atom_ptr->x;
                      y = other_atom_ptr->y;
                      z = other_atom_ptr->z;

                      /* Calculate distance */
                      dist_2 = (x - x1) * (x - x1)
                        + (y - y1) * (y - y1)
                        + (z - z1) * (z - z1);
                          
                      /* If within 5 Angstroms, then keep */
                      if (dist_2 < 25.0)
                        display = TRUE;
                    }

                  /* Get pointer to next atom in linked-list */
                  other_atom_ptr = other_atom_ptr->next_atom_ptr;
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }

          /* Start at first of resaidue's atoms again */
          atom_ptr = residue_ptr->first_atom_ptr;
          natoms = residue_ptr->natoms;
          iatom = 0;

          /* Loop through to set their display states */
          while (atom_ptr != NULL && iatom < natoms)
            {
              /* Set atom's display state */
              atom_ptr->display_atom = display;

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next_atom_ptr;
              iatom++;
            }

          /* If residue to be displayed as ligand, increment count */
          if (display == TRUE)
            (*nligands)++;
        }

      /* Get pointer to next residue in linked-list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

get_max_min_coords  -  Get the maximum and minimum x-, y-, and z-coords
                       of all atoms

***********************************************************************/

void get_max_min_coords(int nfiles,struct limits *limit)
{
  int ifile;

  double x, y, z;

  struct atom *atom_ptr;
  struct sequence *sequence_ptr;

  /* Loop over all the files */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* Initialise pointer to first atom of this structure */
      atom_ptr = sequence_ptr->first_atom_ptr;

      /* Loop through all the atoms to check whether we
         already have this one */
      while (atom_ptr != NULL)
        {
          /* Include atom only if it is to be displayed */
          if (atom_ptr->display_atom == TRUE)
            {
              /* Get the atom's coordinates */
              x = atom_ptr->x;
              y = atom_ptr->y;
              z = atom_ptr->z;

              /* If this is the first atom encountered, store its
                 coords as the present max and min values */
              if (limit->first == TRUE)
                {
                  limit->min_x = x;
                  limit->max_x = x;
                  limit->min_y = y;
                  limit->max_y = y;
                  limit->min_z = z;
                  limit->max_z = z;
                  limit->first = FALSE;
                }

              /* Otherwise, update the limits if coordinates are
                 outside them */
              else
                {
                  if (x < limit->min_x)
                    limit->min_x = x;
                  if (y < limit->min_y)
                    limit->min_y = y;
                  if (z < limit->min_z)
                    limit->min_z = z;
                  if (x > limit->max_x)
                    limit->max_x = x;
                  if (y > limit->max_y)
                    limit->max_y = y;
                  if (z > limit->max_z)
                    limit->max_z = z;
                }
            }

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next_atom_ptr;
        }
    }
}
/***********************************************************************

write_script_file_labels  -  Label the dummy atoms (to be written out
                             later) with each file name

***********************************************************************/

void write_script_file_labels(FILE *outfile_ptr,int nfiles,int out_number)
{
  int atom_number, ifile;

  struct surface *surface_ptr;
  struct sequence *sequence_ptr;

  /* Store next available atom-number */
  atom_number = out_number;
  atom_number++;
  if (atom_number > MAX_ATOMNO)
    return;

  /* Show comment line */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Dummy atom-labels\n");
  fprintf(outfile_ptr,"select atomno = %d\n",atom_number);
  fprintf(outfile_ptr,"label \"Key:\"\n");
  fprintf(outfile_ptr,"colour atoms red\n");
  fprintf(outfile_ptr,"\n");

  /* Loop through all the files */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* Process if this entry has a label */
      if (sequence_ptr->label[0] != '\0')
        {
          /* Increment atom number */
          atom_number++;

          /* Select this atom */
          fprintf(outfile_ptr,"# Label for file %s\n",sequence_ptr->pdb_code);
          fprintf(outfile_ptr,"select atomno = %d\n",atom_number);

          /* Colour according to colour for this file */
          fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
                  sequence_ptr->colour_numbers,sequence_ptr->colour_name);

          /* Show as small sphere */
          fprintf(outfile_ptr,"spacefill 0.2\n");

          /* Add the label */
          fprintf(outfile_ptr,"label %s\n",sequence_ptr->label);
          fprintf(outfile_ptr,"\n");
        }
    }

  /* Add labels for the binding-site masks */

  /* Set mask pointer to the first mask */
  surface_ptr = first_surface_ptr;

  /* Loop over all the masks to be written out */
  while (surface_ptr != NULL)
    {
      /* Increment atom number */
      atom_number++;

      /* Select this atom */
      fprintf(outfile_ptr,"# Label for mask %d\n",surface_ptr->surface_number);
      fprintf(outfile_ptr,"select atomno = %d\n",atom_number);

      /* Colour according to colour for this file */
      fprintf(outfile_ptr,"colour atoms %s  # (%s)\n",
              surface_ptr->colour_numbers,surface_ptr->colour_name);

      /* Show as small sphere */
      fprintf(outfile_ptr,"spacefill 0.2\n");

      /* Add the label */
      fprintf(outfile_ptr,"label %s(mask)\n",surface_ptr->site_surface_code);
      fprintf(outfile_ptr,"\n");

      /* Get pointer to next mask record */
      surface_ptr = surface_ptr->next_surface_ptr;
    }
}
/***********************************************************************

write_script_tail  -  Write all the final records to the output RasMol
                      script file

***********************************************************************/

void write_script_tail(FILE *outfile_ptr,int nfiles,
                       struct parameters params)
{
  /* Close off script */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"select all\n");
  fprintf(outfile_ptr,"set picking distance\n");

  /* Switch hydrogen bonds on */
  /*  fprintf(outfile_ptr,"colour hbonds cyan\n");
      fprintf(outfile_ptr,"hbonds on\n"); */

  /* If structures fitted using FASTA alignment from SAS, the show
     message */
  if (params.cath_alignment == TRUE && nfiles > 1)
    {
      /* Show basis of structural fitting */
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"echo \"* Structures fitted on ");
      if (params.fit_to_backbone == TRUE)
        fprintf(outfile_ptr,"C-alpha atoms ");
      else
        fprintf(outfile_ptr,"N, CA, C, O atoms ");
      if (params.from_fasta == TRUE)
        fprintf(outfile_ptr,"using FASTA alignment\"\n");
      else if (params.upload == TRUE)
        fprintf(outfile_ptr,"using sequence alignment\"\n");
      else
        fprintf(outfile_ptr,"using CATH alignment\"\n");
    }

  /* Show the colouring scheme used */
  if (params.sas_colour == TRUE)
    {
      /* Residues coloured by colouring in SAS alignment */
      fprintf(outfile_ptr,"echo \"* Residues coloured according to ");
      fprintf(outfile_ptr,"colouring in SAS alignment\"\n");
    }
  else if (params.colour_by_structure == TRUE)
    {
      /* Each structure a different colour */
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"echo \"* Each structure shown in a different ");
      fprintf(outfile_ptr,"colour\"\n");
    }

  /* Give some hints */
  fprintf(outfile_ptr,"echo \" \"\n");
  if (params.colour_by_structure == TRUE)
    {
      fprintf(outfile_ptr,"echo \"* To reset atom colours:\"\n");
      fprintf(outfile_ptr,"echo \"       RasMol> select all\"\n");
      fprintf(outfile_ptr,"echo \"       RasMol> colour atoms cpk\"\n");
      fprintf(outfile_ptr,"echo \" \"\n");
    }

  fprintf(outfile_ptr,"exit\n");
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

write_conect_records  -  Write out all the stored CONECT records, or
                         just those for the given structure

***********************************************************************/

void write_conect_records(FILE *outfile_ptr,struct sequence *sequence_ptr)
{
  int iconect, nconect, nwritten, ncon;
  int wanted;

  struct conect *conect_ptr;

  /* Initialise variables */
  nconect = INT_MAX;
  ncon = 0;

  /* Initialise pointer to first CONECT record */
  if (sequence_ptr != NULL)
    {
      conect_ptr = sequence_ptr->first_conect_ptr;
      nconect = sequence_ptr->nconect;
    }
  else
    conect_ptr = first_conect_ptr;

  /* Loop through all the CONECT records, writing each one out */
  while (conect_ptr != NULL && ncon < nconect)
    {
      /* Check whether this conect record is wanted */
      wanted = TRUE;
      if (conect_ptr->atom_ptr[0] == NULL ||
          conect_ptr->atom_ptr[0]->written == FALSE)
        wanted = FALSE;
      nwritten = 0;
      for (iconect = 1; iconect < MAXCONECT && wanted == TRUE; iconect++)
        {
          /* If this atom was written out, then increment count */
          if (conect_ptr->atom_ptr[iconect] != NULL &&
              conect_ptr->atom_ptr[iconect]->written == TRUE)
            nwritten++;
        }
      if (nwritten == 0)
        wanted = FALSE;

      /* If conect details are wanted, then write out */
      if (wanted == TRUE)
        {
          /* Start this record */
          fprintf(outfile_ptr,"CONECT");

          /* Loop through all the atom numbers in this record */
          for (iconect = 0; iconect < MAXCONECT; iconect++)
            {
              /* If this atom number if not zero, then write out */
              if (conect_ptr->atom_ptr[iconect] != NULL &&
                  conect_ptr->atom_ptr[iconect]->out_number < MAX_ATOMNO)
                fprintf(outfile_ptr,"%5d",
                        conect_ptr->atom_ptr[iconect]->out_number);
            }

          /* Finish this record */
          fprintf(outfile_ptr,"\n");
        }

      /* Increment count of this sequences conect records */
      ncon++;

      /* Get pointer to next conect in linked-list */
      conect_ptr = conect_ptr->next_conect_ptr;
    }
}
/***********************************************************************

write_atom_records  -  Write out all the stored ATOM records

***********************************************************************/

int write_atom_records(FILE *outfile_ptr,int nfiles,
                       struct parameters params)
{
  char atom_name[5], ch, chain, file_name[FILENAME_LEN], last_chain;
  char chain_code[2];

  int ifile, out_number;

  struct atom *atom_ptr;
  struct residue *residue_ptr;
  struct sequence *sequence_ptr;

  FILE *attribfile_ptr, *splitfile_ptr;

  /* Initialise variables */
  out_number = 0;

  /* If writing split PDB files, open Astex-style attribute file and
     write header records to it */
  if (params.split_pdb == TRUE)
    {
      /* Form name of output PDB file */
      if (params.work_dir[0] != '\0')
        {
          strcpy(file_name,params.work_dir);
          strcat(file_name,DIR_SEP);
        }
      else
        file_name[0] = '\0';

      /* Add file name */
      strcat(file_name,"astex.atr");

      /* Open the file */
      attribfile_ptr = open_output_file(file_name,TRUE);

      /* Write out the header records */
      fprintf(attribfile_ptr,"Atribute file Version 1.0\n");
      fprintf(attribfile_ptr,"\n");
      fprintf(attribfile_ptr,"Group \"group\" ");
    }

  /* Loop over all the files */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* If also writing to a separate PDB file, open the file */
      if (params.split_pdb == TRUE)
        {
          /* Form name of output PDB file */
          if (params.work_dir[0] != '\0')
            {
              strcpy(file_name,params.work_dir);
              strcat(file_name,DIR_SEP);
            }
          else
            file_name[0] = '\0';

          /* Add PDB code and chain id if present */
          strcat(file_name,sequence_ptr->pdb_code);
          if (sequence_ptr->chain != ' ')
            {
              sprintf(chain_code,"%c",sequence_ptr->chain);
              strcat(file_name,chain_code);
            }
          strcat(file_name,".pdb");

          /* Open the file */
          splitfile_ptr = open_output_file(file_name,TRUE);

          /* Add identifier to the attribute file */
          fprintf(attribfile_ptr," %s",sequence_ptr->pdb_code);
          if (sequence_ptr->chain != ' ')
            fprintf(attribfile_ptr,"%c",sequence_ptr->chain);
        }

      /* Write out the REMARK records giving the PDB name */
      fprintf(outfile_ptr,"PDBFIL %s\n",sequence_ptr->pdb_code);

      /* If resetting atom numbering for each structure, re-initialise */
      if (params.atnum_reset == TRUE)
        out_number = 0;

      /* Initialise pointer to first atom of this structure */
      atom_ptr = sequence_ptr->first_atom_ptr;
      last_chain = '\0';

      /* Loop through all the atoms, writing each one out */
      while (atom_ptr != NULL)
        {
          /* Get this atom's residue record */
          residue_ptr = atom_ptr->residue_ptr;

          /* Include only if this atom is to be displayed */
          if (atom_ptr->display_atom == TRUE &&
              (params.fit_out == FALSE ||
               residue_ptr->fit_used == TRUE ||
               residue_ptr->type == LIGAND ||
               residue_ptr->type == METAL ||
               residue_ptr->type == NUCLEIC_ACID))
            {
              /* Get the atom's out-number */
              if (params.atnum_reset == TRUE)
                {
                  out_number++;
                  atom_ptr->out_number = out_number;
                }
              out_number = atom_ptr->out_number;

              /* Get this atom type and chain-id */
              strcpy(atom_name,atom_ptr->atom_name);
              chain = residue_ptr->chain;

              /* If this is a new chain, then end previous chain with
                 a TER record */
              if (chain != last_chain && last_chain != '\0')
                fprintf(outfile_ptr,"TER\n");
              last_chain = chain;

              /* Check for unusual atom-naming of some residues which cause
                 Rasmol to go awry */
              if (atom_name[0] == 'A' || atom_name[0] == 'E' ||
                  atom_name[0] == 'N' || atom_name[0] == 'U')
                {
                  /* If this is C, N, O, or P, then blank out the initial
                     A or N */
                  ch = atom_name[1];
                  if (ch == 'C' || ch == 'N' || ch == 'O' || ch == 'P')
                    atom_name[0] = ' ';
                }

              /* Write whether this is an ATOM or HETATM record */
              if (atom_ptr->at_het == 'A')
                fprintf(outfile_ptr,"ATOM  ");
              else
                fprintf(outfile_ptr,"HETATM");

              /* Write the atom and residue details out */
              fprintf(outfile_ptr,"%5d %s %s %c%s   %8.3f%8.3f%8.3f%s\n",
                      out_number,atom_name,residue_ptr->res_name,
                      chain,residue_ptr->res_num,
                      atom_ptr->x,atom_ptr->y,atom_ptr->z,
                      atom_ptr->occup_bvalue);

              /* Mark atom as written */
              atom_ptr->written = TRUE;

              /* Write to split file, if required */
              if (params.split_pdb == TRUE)
                {
                  /* Write record identifier */
                  if (atom_ptr->at_het == 'A')
                    fprintf(splitfile_ptr,"ATOM  ");
                  else
                    fprintf(splitfile_ptr,"HETATM");

                  /* Write the atom and residue details out */
                  fprintf(splitfile_ptr,"%5d %s %s %c%s   "
                          "%8.3f%8.3f%8.3f%s\n",
                          out_number,atom_name,residue_ptr->res_name,
                          chain,residue_ptr->res_num,
                          atom_ptr->x,atom_ptr->y,atom_ptr->z,
                          atom_ptr->occup_bvalue);
                }
            }

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next_atom_ptr;
        }

      /* Throw a TER separator between files */
      fprintf(outfile_ptr,"TER\n");

      /* If writing each structure in its own block, write out this
         one's CONECT records */
      if (params.atnum_reset == TRUE)
        write_conect_records(outfile_ptr,sequence_ptr);

      /* If also writing to a separate PDB file, close the file */
      if (params.split_pdb == TRUE)
        fclose(splitfile_ptr);
    }

  /* If also writing to a separate PDB file, close the attribute file */
  if (params.split_pdb == TRUE)
    {
      /* Write out the closing line to the attribure file */
      fprintf(attribfile_ptr," {\n");
      fprintf(attribfile_ptr," }\n");

      /* Close the files */
      fclose(attribfile_ptr);
    }

  /* Return last atom-number used */
  return(out_number);
}
/***********************************************************************

write_surface_vertices  -  Write vertices from .vtx file out as HETATM
                           records in PDB format for display by RasMol

***********************************************************************/

void write_surface_vertices(FILE *outfile_ptr)
{
  int atom_number;
  int written;

  double diagonal, dist_sqrd, grid_spacing_sqrd;
  double x, y, z, other_x, other_y, other_z;

  struct surface *surface_ptr;
  struct vertex *vertex_ptr, *other_vertex_ptr;

  /* Set mask pointer to the first mask */
  surface_ptr = first_surface_ptr;

  /* Loop over all the masks to be written out */
  while (surface_ptr != NULL)
    {
      /* Get this surface's grid-spacing and calculate maximum distance
         that lines joining vertices can span */
      diagonal = 1.1 * surface_ptr->grid_spacing;
      grid_spacing_sqrd = diagonal * diagonal;

      /* Get the pointer to this mask's first vertex */
      vertex_ptr = surface_ptr->first_vertex_ptr;

      /* Loop through all this mask's vertices */
      while (vertex_ptr != NULL)
        {
          /* Get this vertex's transformed coordinates */
          x = vertex_ptr->x;
          y = vertex_ptr->y;
          z = vertex_ptr->z;
          atom_number = vertex_ptr->atom_number;

          /* Initialise flag determining whether to write this
             atom out */
          written = FALSE;

          /* Initialise search for all neighbours of current vertex */
          other_vertex_ptr = vertex_ptr->next_vertex_ptr;

          /* Loop through all subsequent vertices */
          while (other_vertex_ptr != NULL)
            {
              /* Get this vertex's transformed coordinates */
              other_x = other_vertex_ptr->x;
              other_y = other_vertex_ptr->y;
              other_z = other_vertex_ptr->z;

              /* Check whether this vertex is close enough to be a
                 neighbour */
              dist_sqrd = (x - other_x) * (x - other_x)
                + (y - other_y) * (y - other_y)
                + (z - other_z) * (z - other_z);

              /* If within the grid-spacing distance, then have a neighbour */
              if (dist_sqrd <= grid_spacing_sqrd &&
                  atom_number < MAX_ATOMNO)
                {
                  /* If haven't yet written out the current vertex, then
                     do so now */
                  if (written == FALSE)
                    {
                      /* Write out the coordinates of current vertex */
                      fprintf(outfile_ptr,"HETATM%5d  C   MAP 1  %2d    ",
                              atom_number,surface_ptr->surface_number);
                      fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",
                              x,y,z);

                      /* Set flag that vertex written */
                      written = TRUE;
                    }

                  /* Write out the coordinates of neighbour vertex */
                  fprintf(outfile_ptr,"HETATM%5d  C   MAP 1  %2d    ",
                          other_vertex_ptr->atom_number,
                          surface_ptr->surface_number);
                  fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",
                          other_x,other_y,other_z);
                }

              /* Get pointer to the next one */
              other_vertex_ptr = other_vertex_ptr->next_vertex_ptr;
            }

          /* If current vertex was written out, then add TER record */
          if (written == TRUE)
            fprintf(outfile_ptr,"TER\n");

          /* Get pointer to the next vertex */
          vertex_ptr = vertex_ptr->next_vertex_ptr;
        }

      /* Get pointer to next mask record */
      surface_ptr = surface_ptr->next_surface_ptr;
    }
}
/***********************************************************************

write_hbonds  -  Write out the hydrogen bonds, each as 3 atoms: donor,
                 acceptor, and mid-point for label showing H-bond length

***********************************************************************/

int write_hbonds(FILE *file_ptr,int nfiles,int out_number,
                 struct parameters params)
{
  int ifile, nbond;

  double mid_x, mid_y, mid_z, x1, x2, y1, y2, z1, z2;
  double mid_x1, mid_y1, mid_z1, mid_x2, mid_y2, mid_z2;

  struct atom *atom_ptr, *other_atom_ptr;
  struct hbond *hbond_ptr;
  struct sequence *sequence_ptr;

  FILE *outfile_ptr;

  /* Initialise variables */
  nbond = 0;

  /* Get the file pointer */
  outfile_ptr = file_ptr;

  /* @@@ */

  /* Loop through all the files, writing out a dummy atom for each */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* Get pointer to the first of this sequence's H-bonds */
      hbond_ptr = sequence_ptr->first_hbond_ptr;

      /* Loop through all the H-bonds */
      while (hbond_ptr != NULL)
        {
          /* Increment bond number */
          nbond++;

          /* Get the two atoms involved in this bond */
          atom_ptr = hbond_ptr->acceptor_atom_ptr;
          other_atom_ptr = hbond_ptr->donor_atom_ptr;

          /* Get their x-, y-, z-coords */
          x1 = atom_ptr->x;
          y1 = atom_ptr->y;
          z1 = atom_ptr->z;
          x2 = other_atom_ptr->x;
          y2 = other_atom_ptr->y;
          z2 = other_atom_ptr->z;

          /* Calculate coords of mid-point */
          mid_x = (x1 + x2) / 2.0;
          mid_y = (y1 + y2) / 2.0;
          mid_z = (z1 + z2) / 2.0;

          /* Calculate coords of quarter- and three-quarter points */
          mid_x1 = (x1 + mid_x) / 2.0;
          mid_y1 = (y1 + mid_y) / 2.0;
          mid_z1 = (z1 + mid_z) / 2.0;
          mid_x2 = (x2 + mid_x) / 2.0;
          mid_y2 = (y2 + mid_y) / 2.0;
          mid_z2 = (z2 + mid_z) / 2.0;

          /* Write out the three atoms */
          out_number++;
          fprintf(outfile_ptr,"HETATM%5d  A   HB   %4d    ",out_number,
                  nbond);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",x1,y1,z1);

          out_number++;
          fprintf(outfile_ptr,"HETATM%5d  B   HB   %4d    ",out_number,
                  nbond);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",
                  mid_x1,mid_y1,mid_z1);

          out_number++;
          fprintf(outfile_ptr,"HETATM%5d  C   HB   %4d    ",out_number,
                  nbond);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",
                  mid_x,mid_y,mid_z);

          out_number++;
          fprintf(outfile_ptr,"HETATM%5d  D   HB   %4d    ",out_number,
                  nbond);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",
                  mid_x2,mid_y2,mid_z2);

          out_number++;
          fprintf(outfile_ptr,"HETATM%5d  E   HB   %4d    ",out_number,
                  nbond);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",x2,y2,z2);

          /* Write out a TER record */
          fprintf(outfile_ptr,"TER   \n");

          /* Get pointer to next H-bond */
          hbond_ptr = hbond_ptr->next_hbond_ptr;
        }
    }

  /* Return last atom-number used */
  return(out_number);
}
/***********************************************************************

write_dummy_atoms  -  Write out the dummy atoms corresponding to each
                      file

***********************************************************************/

void write_dummy_atoms(FILE *outfile_ptr,int nfiles,int out_number,
                       struct limits limit)
{
  int atom_number, ifile;

  double x, y, z;

  struct surface *surface_ptr;
  struct sequence *sequence_ptr;

  /* Store next available atom-number */
  atom_number = out_number;

  /* Initialise starting coordinates */
  x = limit.max_x + 5.0;
  y = limit.max_y;
  z = (limit.min_z + limit.max_z) / 2.0;

  /* Write out the key */
  atom_number++;
  if (atom_number < MAX_ATOMNO)
    return;
  fprintf(outfile_ptr,"HETATM%5d  C   XXX 1   1    ",atom_number);
  fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",x,y,z);

  /* Initialise starting coordinates */
  y = limit.max_y - 2.0;

  /* Loop through all the files, writing out a dummy atom for each */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* Process if this entry has a label */
      if (sequence_ptr->label[0] != '\0' && atom_number < MAX_ATOMNO)
        {
          /* Increment atom number */
          atom_number++;

          /* Write out this dummy atom */
          fprintf(outfile_ptr,"HETATM%5d  C   XXX 1   1    ",atom_number);
          fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",x,y,z);

          /* Drop the y-coordinate for the next dummy atom */
          y = y - 2.0;
        }
    }

  /* Repeat for the binding-site masks */

  /* Set mask pointer to the first mask */
  surface_ptr = first_surface_ptr;

  /* Loop over all the masks to be written out */
  while (surface_ptr != NULL && atom_number < MAX_ATOMNO)
    {
      /* Increment atom number */
      atom_number++;

      /* Write out this dummy atom */
      fprintf(outfile_ptr,"HETATM%5d  C   XXX 1   1    ",atom_number);
      fprintf(outfile_ptr,"%8.3f%8.3f%8.3f  1.00 10.00\n",x,y,z);

      /* Drop the y-coordinate for the next dummy atom */
      y = y - 2.0;

      /* Get pointer to next mask record */
      surface_ptr = surface_ptr->next_surface_ptr;
    }
}
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char out_name[FILENAME_LEN];
  char lname[15], pname[15], rname[15];
  char pdb_name[FILENAME_LEN], reslist_name[FILENAME_LEN];
  char cath_name[FILENAME_LEN];
  char matfile_name[FILENAME_LEN], protein_family[FILENAME_LEN];
  char dfn_name[FILENAME_LEN];
  char grow_out_file[FILENAME_LEN], promotif_dat[FILENAME_LEN];

  int ifile, ifit, last_number, nalign_resid;
  int natoms, nbonds, nfiles, nresid, nseq, out_number;
  int nfit, nhbonds, nligs, nligands;
  int fit_flag, fit_flagged, have_dfn, transform;

  double fit_cog[3], matrix[4][4], rmsd;

  struct c_file_data file_name_ptr;
  struct limits limit;
  struct parameters params;

  struct range *first_range_ptr;
  struct rasdef *first_rasdef_ptr;
  struct sequence *sequence_ptr;

  FILE *outfile_ptr;

  /* Initialise variables */
  verbose = FALSE;
  fit_cog[0] = 0.0;
  fit_cog[1] = 0.0;
  fit_cog[2] = 0.0;
  fit_flagged = FALSE;
  first_conect_ptr = NULL;
  last_conect_ptr = NULL;
  limit.first = TRUE;
  nalign_resid = 0;
  nbonds = 0;
  nfit = 0;
  nligs = 0;
  nligands = 0;
  nresid = 0;
  nseq = 0;
  out_name[0] = '\0';
  out_number = 0;
  rmsd = 0.0;
  strcpy(rname,"reslist");
  dfn_name[0] = '\0';

  /* Initialise pointers */
  first_rasdef_ptr = NULL;

  /* Interpret the command-line parameters */
  nfiles = get_command_arguments(argv,argc - 1,protein_family,cath_name,
                                 &params);

  /* Sort out the display options for all the sequences */
  assign_display_options(nfiles,params,&nligands);

  /* If writing to file, then set file name */
  if (params.write_to_file == TRUE)
    strcpy(out_name,"romlas.out");

  /* Read in all the directory paths from the parameter file, CATHPARAM */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Open the output file */
  outfile_ptr = open_output_file(out_name,TRUE);

  /* If writing  PyMOL script, write header records */
  if (params.pymol == TRUE)
    write_pymol_headers(outfile_ptr,nfiles,file_name_ptr,params);

  /* If writing RasMol script, write header records for script */
  else if (params.rasmol_script == TRUE)
    write_script_headers(outfile_ptr,nfiles,file_name_ptr,params);

  /* If have a CATH alignment file, read in to get the residue 
     equivalences in the alignment, for use when fitting */
  if (params.cath_alignment == TRUE)
    fit_flagged
      = get_cora_alignment(cath_name,&nseq,&nalign_resid,&nfiles,params);

  /* Loop over all the PDB codes to get all the ATOM records and
     associated CONECT records */
  for (ifile = 0; ifile < nfiles; ifile++)
    {
      /* Initialise variables for this file */
      nhbonds = 0;
      if (ifile == 0)
        fit_flag = FALSE;
      else
        fit_flag = fit_flagged;

      /* Get the sequence pointer for the current PDB entry */
      sequence_ptr = stored_sequence_ptr[ifile];

      /* Get the input file names for this PDB code */
      form_file_names(sequence_ptr,pdb_name,dfn_name,reslist_name,
                      protein_family,matfile_name,grow_out_file,
                      promotif_dat,ifile,file_name_ptr,&params);

      /* Read in the RasMol definitions */
      have_dfn
        = get_dfn_data(sequence_ptr->pdb_code,dfn_name,&first_rasdef_ptr);

      /* If looking for the residues matching the given PROSITE
         pattern, pick these up */
      if (params.prosite_pattern == TRUE)
        get_prosite_patterns(sequence_ptr,file_name_ptr,params);

      /* Otherwise, if showing ligand, or ligand and interacting residues,
         read in the list of relevant residues from the ligplot.res or
         family.res files */
      else if (sequence_ptr->display_option == LIGAND_ONLY ||
               sequence_ptr->display_option == INT_RESIDUES ||
               sequence_ptr->display_option == LIGAND_AND_RES ||
               sequence_ptr->display_option == LIGAND_AND_BOX ||
               sequence_ptr->display_option == BINDING_SITE)
        get_residue_list(sequence_ptr,reslist_name,&nresid,params);

      /* If showing PROMOTIF motif, get its details from the promotif.dat
         file in the PDBsum directory */
      if (params.promotif == TRUE)
        get_motif(sequence_ptr,promotif_dat,&first_range_ptr,outfile_ptr,
                  &params);

      /* If picking up the SITE residues for a given site, read through
         the PDB file to identify the list of residues involved */
      if (params.pdb_site == TRUE)
        get_site(sequence_ptr,pdb_name,outfile_ptr,&params);

      /* Read through the PDB file to pick up all the required residues
         and associated CONECT records */
      natoms = get_atom_records(sequence_ptr,pdb_name,&out_number,
                                params,file_name_ptr.other_dir[CATH_DOMAINS],
                                first_range_ptr,first_rasdef_ptr,
                                outfile_ptr);

      /* If have at least one atom, proceed */
      if (natoms > 0)
        {
          /* Keep/discard het groups requiring distance checks */
          check_distance(ifile,&nligs);

          /* Read in the transformation for this structure from the
             appropriate family.mat file, if present and if required */
          if ((nfiles > 1 || params.nsite_surface_files > 1)
              && protein_family[0] != '\0')
            transform = get_transformation(matfile_name,matrix);

          /* If this is a CATH alignment, then calculate the transformation
             by fitting onto the first structure */
          else if (nfiles > 1 && params.cath_alignment == TRUE)
            transform = calculate_fit(sequence_ptr,matrix,fit_cog,&rmsd,
                                      &nfit,fit_flag,params);

          /* Otherwise, no transformation required */
          else
            transform = FALSE;

          /* Apply the transformation to all atoms */
          if (transform == TRUE)
            apply_transformation(sequence_ptr,matrix,fit_cog);

          /* For LIGPLOT display, pick up all the hydrogen bonds from
             the grow.out file */
          if (sequence_ptr->display_option == LIGAND_AND_RES)
            nhbonds = get_hbonds(grow_out_file,sequence_ptr);

          /* If writing RasMol script, define the ligand and protein
             residues */
          if (params.rasmol_script == TRUE)
            {
              /* Define protein chains for all cases but PROMOTIF */
              if (params.promotif == FALSE && params.pdb_site == FALSE &&
                  params.csa == FALSE)
                {
                  /* Define and display the ligand residues */
                  write_script_residues(outfile_ptr,ifile,sequence_ptr,
                                        nligands,rmsd,nfit,TRUE,lname,
                                        params);

                  /* Define and display the protein residues */
                  write_script_residues(outfile_ptr,ifile,sequence_ptr,
                                        nfiles,rmsd,nfit,FALSE,pname,
                                        params);
                }

              /* If have user-defined residue lists, highlight these */
              if (params.residue_list[0] != '\0' && params.csa == FALSE)
                write_script_reslist(outfile_ptr,params,rname);

              /* If have a range of residues on which fit performed
                 highlight this range */
              if (sequence_ptr->first_fitlist_ptr != NULL)
                write_script_fitlist(sequence_ptr,outfile_ptr,pname,
                                     ifile,params);

              /* If colouring by the SAS colouring scheme, colour the
                 individual residues */
              if (params.sas_colour == TRUE)
                write_script_sascolour(outfile_ptr,sequence_ptr);

              /* If have to add sidechains onto the current file, then do
                 so */
              if (params.add_sidechains == TRUE)
                write_script_sidechain(outfile_ptr,sequence_ptr,params);

              /* If ligand interacts with any metals, then display these as
                 spheres */
              else
                write_script_metals(outfile_ptr,sequence_ptr);

              /* If have any H-bonds to display, then write these out */
              if (nhbonds > 0)
                nbonds
                  = write_script_hbonds(outfile_ptr,sequence_ptr,nbonds);

              /* If showing a PROMOTIF secondary structure motif, write
                 the script commands to display it */
              if (params.promotif == TRUE)
                write_script_promotif(outfile_ptr,first_range_ptr,params);

              /* If showing the SITE residues, write appropriate script */
              if (params.pdb_site == TRUE)
                write_script_site(outfile_ptr,sequence_ptr,
                                  first_rasdef_ptr,params);

              /* If showing the CSA residues, write appropriate script */
              if (params.csa == TRUE)
                write_script_csa(outfile_ptr,sequence_ptr,first_rasdef_ptr,
                                 params);
            }
        }
    }

  /* Calculate the maximum and minimum coords */
  get_max_min_coords(nfiles,&limit);
  last_number = out_number;

  /* If have any binding-site masks, then write out the script definitions
     for them */
  if (params.nsite_surface_files > 0)
    {
      /* Read through the mask vertices to count their number and to
         adjust minimum and maximum coordinates */
      last_number = read_surface_vertices(out_number,nfiles,protein_family,
                                          &limit,file_name_ptr,params);

      /* If writing RasMol script, define the mask files */
      if (params.rasmol_script == TRUE)
        write_script_mask(outfile_ptr);
    }

  /* If writing RasMol script, then close off the script part */
  if (params.rasmol_script == TRUE)
    {
      /* If multiple files, each coloured differently, add
         dummy atoms with labels showing the PDB entries */
      if (params.colour_by_structure == TRUE && nfiles > 1)
        write_script_file_labels(outfile_ptr,nfiles,last_number);

      /* Close the script off */
      write_script_tail(outfile_ptr,nfiles,params);
    }

  /* Write out all the ATOM records */
  if (params.pdb_type != PDBSUM1 || params.pdb_file != &null)
    out_number = write_atom_records(outfile_ptr,nfiles,params);

  /* If have any H-bonds to write out, then do so */
  if (nbonds > 0)
    out_number = write_hbonds(outfile_ptr,nfiles,out_number,params);

  /* If have any mask-files to write out, then process these */
  if (params.nsite_surface_files > 0)
    write_surface_vertices(outfile_ptr);

  /* If have dummy atoms to write out, write these out */
  if (params.colour_by_structure == TRUE && nfiles > 1)
    write_dummy_atoms(outfile_ptr,nfiles,last_number,limit);

  /* Write out all the CONECT records, unless we have already done so */
  if (params.atnum_reset == FALSE)
    write_conect_records(outfile_ptr,NULL);

  return(0);
}
