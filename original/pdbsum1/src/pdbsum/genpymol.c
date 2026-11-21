/***********************************************************************

  genpymol.c - Program to generate PyMOL script for given PDB structure
               or part structure

************************************************************************

  Date:-         28 Sep 2016

  Written by:-   Roman Laskowski

Testing:

Why are protein cheins being marked as C-alpha?

   ~/src/pdbsum/genpymol 6f1t g -u 10278 -rasmol -labels -v O14576 -range 268 584 -n

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      Description
-----------      -----------
alignments.dat   Input file giving mapping from UniProt sequence to
                 closest PDB residue
CATHPARAM        Input parameter file
seq001.conserv   Inout residue conservation file
genpymol.pdb     Output coords in PDB format
genpymol.pml     Output PyMOL script
grow.out         Input list of protein-ligand/DNA/metal interactions
igroq.out.gz     Input list of protein-protein interactions
map1d_x?.dat     Input map of original image, giving mutation data
pdb.align        Input file of PDB mappings, created by disstruc
pdbXXXX.ent      Input PDB file
rasmol.dfn       Input PDBsum definition file
resdata.dat      Input file in PDBsum directory giving residue
                 annotations
uniplot.dat      Input file in PDBsum directory giving Pfam domain
                 definitions and colours
variants.dat     Output file listing the mapped variants for use on
                 the web page

Compilation
-----------

cc -c genpymol.c
cc -c common.c
cc -c reapar.c
cc -c readpdb.c
cc -c readsum.c
cc -o genpymol genpymol.o common.o reapar.o readpdb.o readsum.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c genpymol.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -O2 -funroll-loops -c readsum.c
cc -o genpymol genpymol.o common.o reapar.o readpdb.o readsum.o -lm -lz

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> extract_res_name (in common.c)
               -> to_upper (in common.c)
               -> get_aa_code (in common.c)
         -> create_addmut_record
         -> get_start_colour
   -> get_paths
         -> get_param (in reapar.c)
         -> store_string
   -> open_output_file (in common.c)
   -> find_pdbsum_dir (in common.c)
   -> get_pdbsum1_dir (in common.c)
   -> read_chain_equivs
         -> string_truncate (in common.c)
         -> create_chain_record
   -> get_data_dir
   -> get_closest_pdb
         -> string_chomp
   -> get_res_conservation
         -> string_chomp (in common.c)
         -> get_tokens (in common.c)
         -> create_residue_record
         -> free_tokens (in common.c)
   -> create_residue_index
   -> get_pdb_mappings
         -> get_chain_string
         -> string_chomp (in common.c)
         -> get_tokens (in common.c)
         -> to_upper (in common.c)
         -> free_tokens (in common.c)
   -> read_database_dat
         -> string_chomp (in common.c)
   -> r_get_pdb_data (in readpdb.c)
   -> s_get_resdata (in readsum.c)
         -> find_pdbsum_dir (in common.c)
         -> string_truncate (in common.c)
   -> s_transfer_resdata (in readsum.c)
         -> r_find_residue (in readpdb.c)
   -> map_to_pdb
   -> mark_local_residues
         -> check_cutoff
   -> delete_residues
   -> get_data_dir
   -> read_seq_fasta
         -> parse_uniprot_name (in common.c)
         -> store_string (in common.c)
   -> initialise_residues
         -> create_residue_record
               -> get_res_name (in common.c)
   -> create_residue_index
   -> get_mutations
         -> string_chomp (in common.c)
         -> find_disease_record
         -> create_disease_record
         -> find_desc_entry
         -> create_dlist_record
         -> create_annot_entry
         -> create_alist_record
         -> create_aamut_record
   -> sort_variants
         -> variant_sort
   -> added_variants
         -> create_desc_entry
         -> create_annot_entry
         -> create_alist_record
         -> create_aamut_record
   -> map_residues
         -> init_residues
         -> get_seq_start
               -> string_chomp (in common.c)
         -> string_chomp (in common.c)
         -> tokenize (in common.c)
         -> get_residue_details
               -> get_res_name (in common.c)
         -> to_upper (in common.c)
         -> r_find_residue (in readpdb.c)
   -> write_variants
         -> get_last_atom
   -> read_uniplot_dat
         -> string_chomp (in common.c)
         -> create_domain_record
   -> map_domains_to_pdb
   -> identify_neighbours
         -> find_pdbsum_dir (in common.c)
         -> s_get_grow_data (in readsum.c)
               -> create_s_grow_record (in readsum.c)
         -> add_covalent_bonds
               -> get_covalent_bonds
                     -> string_truncate (in common.c)
                     -> create_s_grow_record (in readsum.c)
         -> s_init_sort  (in readsum.c)
         -> s_sort_grow_records (in readsum.c)
         -> filter_neighbours
               -> check_residue
         -> s_free_grow_data (in readsum.c)
   -> write_pymol_pml
         -> form_filename (in common.c)
         -> open_output_file (in common.c)
         -> write_pymol_headers
               -> get_chain_col
                     -> get_colour_rgb (in common.c)
         -> write_pymol_chains
               -> remove_residue_blanks
               -> create_pobject_record
               -> get_chain_colname
               -> find_pdbsum_dir (in common.c)
               -> colour_alpha_fold
                     -> get_c_scores (in common.c)
                           -> create_c_score_entry (in common.c)
               -> store_string (in common.c)
         -> write_pymol_ligands
               -> remove_residue_blanks
               -> create_pobject_record
               -> store_string (in common.c)
         -> write_pymol_metals
               -> remove_residue_blanks
               -> create_pobject_record
               -> store_string (in common.c)
         -> write_pymol_variants
               -> get_max_variants
               -> to_lower (in common.c)
               -> remove_blanks (in common.c)
         -> write_pymol_surface
         -> write_pymol_domcol
   -> write_rasmol_script
         -> open_output_file (in common.c)
         -> write_rasmol_headers
               -> get_rasmol_colour
                     -> get_colour_rgb (in common.c)
         -> write_rasmol_protein
               -> get_chain_col
               -> find_pdbsum_dir (in common.c)
               -> colour_alpha_fold_rasmol
                     -> get_c_scores (in common.c)
                           -> create_c_score_entry (in common.c)
         -> write_rasmol_nucleic
               -> get_chain_col
         -> check_ligand_ranges
         -> write_rasmol_ligands
               -> display_ligand
         -> write_rasmol_metals
              -> get_rasmol_colour
         -> write_rasmol_rescons
         -> write_rasmol_domcol
              -> replace_char (in common.c)
         -> write_rasmol_variants
              -> get_max_variants
              -> remove_blanks (in common.c)
              -> write_def
         -> write_coords
               -> r_write_pdb_coords (in readpdb.c)
   -> write_coords
         -> r_write_pdb_coords (in readpdb.c)

*/

#include "common.h"
#include "readpdb.h"
#include "readsum.h"
#include "dcommon.h"
#include "genpymol.h"


/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int extract_res_name(char *residue,char *aac,char *aam,char **aa_ex,
                     char **mut_ex,int *typ,int *mut_len,char **add_var);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
char *form_filename(char *pdb_code,char *filename_template);
void free_tokens(char **token,int ntokens);
int get_c_scores(char *pdb_code,char *pdbsum_dir,
                 struct c_score **fst_c_score_ptr);
void get_colour_rgb(char *c_name,double rgb[3]);
void get_pdbsum1_dir(char *pdir,char *pdb_code,
                     struct c_file_data file_name_ptr);
void get_res_name(char aa_code,char *res_name);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
void parse_uniprot_name(char *string,char *protein,char *species,
                        char *uniprot_acc,char *uniprot_id);
void remove_blanks(char *string);
void replace_char(char *string,char replace,char by);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);
int tokenize(char *string,char **token,int max_tokens,char token_sep);
void to_lower(char *search_string,int length);
void to_upper(char *search_string,int length);

/* In reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/* In readpdb.c */
struct r_residue *r_find_residue(struct r_residue *first_r_residue_ptr,
                                 char *res_name,char *res_num,char chain);
struct r_pdb *r_get_pdb_data(char *pdb_code,
                             struct c_file_data file_name_ptr,
                             int pdb_type,int data_rqrd);
int r_write_pdb_coords(FILE *file_ptr,
                       struct r_residue *first_r_residue_ptr,int nresid,
                       int atom_number,int omit_hydrogens);

/* In readsum.c */
void s_append_grow_records(struct s_grow **fst_s_grow_ptr,
                           struct s_grow *second_s_grow_ptr);
struct s_grow *create_s_grow_record(struct s_grow **fst_s_grow_ptr,
                                    struct s_grow **lst_s_grow_ptr);
void s_free_grow_data(struct s_grow *first_s_grow_ptr);
int s_get_grow_data(char *pdb_code,char *pdbsum_dir,
                    struct s_grow **fst_s_grow_ptr,int type);
int s_get_resdata(char *pdb_code,char chain,int pdb_type,
                  struct r_residue **fst_r_residue_ptr,
                  struct c_file_data file_name_ptr);
void s_init_sort(struct s_grow *first_s_grow_ptr);
void s_sort_grow_records(struct s_grow **fst_s_grow_ptr,int ngrow);
void s_transfer_resdata(struct r_residue *first_r_residue_ptr,
                        struct r_residue *first_pdb_r_residue_ptr);


/***********************************************************************

create_addmut_record  -  Create a new addmut record

***********************************************************************/

struct addmut *create_addmut_record(struct addmut **fst_addmut_ptr,
                                    struct addmut **lst_addmut_ptr,
                                    int seq_no,char aa_code,char aa_mut,
                                    char *added_variant,char *aa_extra,
                                    char *mut_extra,int mut_length,
                                    int type)
{
  int i;

  struct addmut *addmut_ptr, *first_addmut_ptr, *last_addmut_ptr;

  /* Initialise addmut pointers */
  addmut_ptr = NULL;
  first_addmut_ptr = *fst_addmut_ptr;
  last_addmut_ptr = *lst_addmut_ptr;

  /* Create addmut entry */
  addmut_ptr = (struct addmut*)malloc(sizeof(struct addmut));
  if (addmut_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct addmut\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  addmut_ptr->seq_no = seq_no;
  addmut_ptr->aa_code = aa_code;
  addmut_ptr->aa_mut = aa_mut;
  store_string(&(addmut_ptr->added_variant),added_variant);
  addmut_ptr->aa_extra = aa_extra;
  addmut_ptr->mut_extra = mut_extra;
  addmut_ptr->type = type;
  addmut_ptr->mut_length = mut_length;

  /* Initialise pointer to next addmut record */
  addmut_ptr->next_addmut_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of addmuted addmut */
  if (first_addmut_ptr == NULL)
    first_addmut_ptr = addmut_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_addmut_ptr->next_addmut_ptr = addmut_ptr;
                      
  /* Save the current addmut's pointer */
  last_addmut_ptr = addmut_ptr;

  /* Return current pointers */
  *fst_addmut_ptr = first_addmut_ptr;
  *lst_addmut_ptr = last_addmut_ptr;
  return(addmut_ptr);
}
/***********************************************************************

get_start_colour  -  Get the domain start colour

***********************************************************************/

int get_start_colour(char *uniprot_acc)
{
  char u_acc[UNIPROT_LEN];

  char *string_ptr;

  int i, prot_no, start_col;

  /* Replace hyphen with underscore, if required */
  strcpy(u_acc,uniprot_acc);
  string_ptr = strchr(u_acc,'-');
  if (string_ptr != NULL)
    string_ptr[0] = '_';

  /* Get the starting colour for this protein */
  prot_no = -1;
  for (i = 0; i < NCOVID_PROTEINS && prot_no == -1; i++)
    {
      /* If match found, save the protein number */
      if (!strcmp(u_acc,COVID_PROTEIN[i]))
        prot_no = i;
    }

  /* If protein not found, use default colour */
  if (prot_no == -1)
    prot_no = 0;

  /* Get the corresponding starting colour */
  start_col = 3 * prot_no;
  start_col = start_col % NPFAM_COLOURS;

  /* Return the starting colour */
  return(start_col);
}
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int ntoken,
                           struct parameters *params)
{
  char aa_code, aa_mut, ch, string[3][20];

  char *aa_extra, *added_variant, *mut_extra;
  
  int istring, itoken, len, mut_length, seq_no, type;

  struct addmut *addmut_ptr, *last_addmut_ptr;

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->uniprot_acc[0] = '\0';
  params->chain = '*';
  params->covid19 = FALSE;
  params->cutoff_dist = CUTOFF_DIST;
  params->disastr = FALSE;
  params->pml_file = &null;
  params->neighbours = FALSE;
  params->labels = FALSE;
  params->start_colour = 0;
  params->start = 0;
  params->end = 0;
  params->pdb_type = STANDARD_PDB;
  params->output_type = PYMOL_SCRIPT;
  params->user_id = 0;
  params->variants = FALSE;
  params->res_cons = FALSE;
  params->colour_by_domain = FALSE;
  params->run_64 = TRUE;
  params->naddmut = 0;
  params->first_addmut_ptr = last_addmut_ptr = NULL;

  /* If -help entered on the command line, show options available */
  if (ntoken < 1 || (ntoken == 1 && (!strcmp(strptrs[1],"-help") ||
                                     !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n");
      printf("\n");
      printf("    genpymol  pdb_code  [chain]  [-n]  [-c]  "
             "[-rasmol|-jmol|-JSmol]  [-disastr]\n"
             "              [-labels]  [-v uniprot_acc]  "
             "[-range start end]  [-var XnY]\n"
             "              [-u user_id]  [-covid19 uniprot_acc] "
             "[-cons]  [-domcol]  [-pml filename]\n"
             "              [UPLOAD|MCSG]  [-pdbsum1]\n");
      printf("where\n");
      printf("     * pdb_code = PDB code for file to be processed\n");
      printf("     * chain    = Chain identifier\n");
      printf("     * -disastr = Flag indicating that running for VarSite\n");
      printf("     * -n       = Retain neighbours of given chain only\n");
      printf("     * -jmol    = Output a Jmol script and coords\n");
      printf("     * -JSmol   = Output a RasMol script and coords to "
             "separate files\n");
      printf("     * -rasmol  = Output a RasMol script\n");
      printf("     * -labels  = Label the variants\n");
      printf("     * -c       = Output the coords only\n");
      printf("     * -u       = User-id indicating tmp dir with coords\n");
      printf("     * -v uniprot_acc = Show variants from corresponding "
             "VarSite UniProt entry\n");
      printf("     * -range start end = Range of UniProt residues mapped "
             "to the PDB entry\n");
      printf("     * -var XnY = user-added variant(s); n=seq no, "
             "X,Y=optional variant aa codes\n"
             "                eg  80, L34R, KV101QT, G99Ter\n");
      printf("     * UPLOAD   = Uploaded PDB file for PDBsum generation\n");
      printf("     * MCSG     = MCSG structure\n");
      printf("     * -varsite uniprot_acc = VarSite UniProt accession\n");
      printf("     * -covid19 uniprot_acc = Covid19 protein\n");
      printf("     * -cut dist = Cut-off dist for local residues\n");
      printf("     * -cons    = Colour by residue conservation\n");
      printf("     * -domcol  = Colour by Pfam domains\n");
      printf("     * -pdbsum1 = PDBsum1 structure\n");
      printf("\n");
      exit(1);
    }

  /* Get the PDB code from the first token */
  strcpy(params->pdb_code,strptrs[1]);
  params->pdb_code[4] = '\0';

  /* Check whether this is an Alpha Fold structure (first character is
     an upper-case letter) */
  ch = params->pdb_code[0];
  if (ch >= 'A' && ch <= 'Z')
    params->pdb_type = ALPHA_FOLD_PDB;

  /* Loop through any additional command arguments */
  for (itoken = 2; itoken < ntoken + 1; itoken++)
    {
      /* Check for the -n option */
      if (!strcmp(strptrs[itoken],"-n"))
        params->neighbours = TRUE;

      /* Check for the -c option */
      else if (!strcmp(strptrs[itoken],"-c"))
        params->output_type = ATOM_COORDS;

      /* Check for the -rasmol option */
      else if (!strcmp(strptrs[itoken],"-rasmol"))
        params->output_type = RASMOL_SCRIPT;

      /* Check for the -jmol option */
      else if (!strcmp(strptrs[itoken],"-jmol"))
        params->output_type = JMOL_SCRIPT;

      /* Check for the -JSmol option */
      else if (!strcmp(strptrs[itoken],"-JSmol"))
        params->output_type = JSMOL;

      /* Check for the -labels option */
      else if (!strcmp(strptrs[itoken],"-labels"))
        params->labels = TRUE;

      /* Check for the -cons option */
      else if (!strcmp(strptrs[itoken],"-cons"))
        params->res_cons = TRUE;

      /* Check for the -domcol option */
      else if (!strcmp(strptrs[itoken],"-domcol"))
        params->colour_by_domain = TRUE;

      /* Check for the -covid19 option */
      else if (!strcmp(strptrs[itoken],"-covid19"))
        {
          /* Get the UniProt accession from the next parameter */
          if (itoken < ntoken)
            {
              strcpy(params->uniprot_acc,strptrs[itoken + 1]);
              itoken++;
              params->covid19 = TRUE;

              /* Get the corresponding start colour if plotting by
                 domains */
              params->start_colour
                = get_start_colour(params->uniprot_acc);
            }
        }

      /* Check for the cut-off distance */
      else if (!strcmp(strptrs[itoken],"-cut"))
        {
          /* Get the distance from the next parameter */
          if (itoken < ntoken)
            {
              params->cutoff_dist = atof(strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the name of the output .pml file */
      else if (!strcmp(strptrs[itoken],"-pml"))
        {
          /* Get the file name from the next parameter */
          if (itoken < ntoken)
            {
              store_string(&(params->pml_file),strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for VarSite UniProt accession */
      else if (!strcmp(strptrs[itoken],"-varsite"))
        {
          /* Get the UniProt accession from the next parameter */
          if (itoken < ntoken)
            {
              strcpy(params->uniprot_acc,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for addition of variants from VarSite */
      else if (!strcmp(strptrs[itoken],"-v"))
        {
          /* Get the UniProt accession from the next parameter */
          if (itoken < ntoken)
            {
              strcpy(params->uniprot_acc,strptrs[itoken + 1]);
              itoken++;
              params->variants = TRUE;
            }          
        }

      /* Check for user-id giving reference to tmp directory */
      else if (!strcmp(strptrs[itoken],"-u"))
        {
          /* Get the user id from the next parameter */
          if (itoken < ntoken)
            {
              params->user_id = atoi(strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for residue range */
      else if (!strcmp(strptrs[itoken],"-range"))
        {
          /* Get the start and end residues from the next two parameters */
          if (itoken < ntoken - 1)
            {
              params->start = atoi(strptrs[itoken + 1]);
              itoken++;
              params->end = atoi(strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for user-supplied variant */
      else if (!strcmp(strptrs[itoken],"-var"))
        {
          /* Get the variant (XnY or XxxnYyy) from the next parameter */
          if (itoken < ntoken)
            {
              /* Extract the sequence number and two residue names */
              seq_no = extract_res_name(strptrs[itoken + 1],&aa_code,
                                        &aa_mut,&aa_extra,&mut_extra,
                                        &type,&mut_length,&added_variant);

              /* If have a valid variant, save */
              if (seq_no > 0)
                {
                  /* Create a record to store this variant */
                  addmut_ptr
                    = create_addmut_record(&(params->first_addmut_ptr),
                                           &last_addmut_ptr,seq_no,
                                           aa_code,aa_mut,added_variant,
                                           aa_extra,mut_extra,mut_length,
                                           type);

                  /* Increment count of added variants */
                  params->naddmut++;
                }
            }
        }

      /* Check for the -disastr option */
      else if (!strcmp(strptrs[itoken],"-disastr"))
        params->disastr = TRUE;

      /* Check for the UPLOAD option */
      else if (!strcmp(strptrs[itoken],"UPLOAD"))
        params->pdb_type = UPLOAD_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the PDBsum1 option */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        params->pdb_type = PDBSUM1;

      /* Is single-character parameter, take it to be the chain id */
      else if (strlen(strptrs[itoken]) == 1)
        params->chain = strptrs[itoken][0];
    }
}
/***********************************************************************

get_data_dir  -  Get the directory where the data are

***********************************************************************/

void get_data_dir(char **ddir,char *root_dir,char *uniprot_acc)
{
  char file_name[C_FILENAME_LEN], parent_dir[C_FILENAME_LEN];

  char *data_dir;

  int have_file, len, loop;
  int error;

  /* Get the length of the UniProt account */
  len = strlen(uniprot_acc);

  /* Form the name of the parent directory */
  strcpy(file_name,root_dir);
  strcat(file_name,"/");
  strcat(file_name,uniprot_acc + len - 2);
  strcat(file_name,"/");
  strcpy(parent_dir,file_name);
  strcat(file_name,uniprot_acc);
  strcat(file_name,"/");
  store_string(&data_dir,file_name);

  /* Return the data directory */
  *ddir = data_dir;
}
/***********************************************************************

get_closest_pdb  -  Get the closest PDB entry and corresponding residue

***********************************************************************/

int get_closest_pdb(char *pdbcode,char *chn,char *uniprot_acc,
                    int res_no,char *code_dir)
{
  char file_name[C_FILENAME_LEN], input_line[LINELEN + 1];
  char aa_code, chain, pdb_code[5], res_num[6], res_name[4];
  char closest_chain, closest_pdb_code[5];

  char *string_ptr;

  int closest_res_no, dist2, len, min_dist2, seq_pos;
  int found;

  FILE *file_ptr;

  /* Initialise variables */
  aa_code = ' ';
  *chn = chain = closest_chain = ' ';
  closest_res_no = 0;
  dist2 = min_dist2 = 0;
  found = FALSE;
  pdbcode[0] = closest_pdb_code[0] = '\0';
  res_num[0] = res_name[0] = '\0';

  /* Form the name of the alignments.dat file */
  strcpy(file_name,code_dir);
  strcat(file_name,"/align/alignments.dat");

  /* Open the alignment file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* Show warning */
      printf("*** ERROR. File [%s] not found\n",file_name);
      exit(1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LONG_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the line at the last non-blank character */
      len = string_chomp(input_line);

      /* Check for name of PDB entry */
      if (!strncmp(input_line,"ALIGN_PDB_CODE[",15))
        {
          /* Extract the PDB code */
          string_ptr = strchr(input_line,' ');
          if (string_ptr != NULL)
            strcpy(pdb_code,string_ptr + 1);
        }

      /* If this is an alignment record, read it in */
      else if (strncmp(input_line,"NALIGN",6) &&
               strchr(input_line,'_') == NULL &&
               strchr(input_line,'[') == NULL)
        {
          /* Get the sequence position */
          string_ptr = strchr(input_line,' ');

          /* If have a number, read in the rest of the line */
          if (string_ptr != NULL)
            {
              /* Get the length of the number */
              len = string_ptr - input_line;

              /* Extract the amino acid code and PDB residue */
              aa_code = input_line[len + 1];
              strncpy(res_name,input_line + len + 5,3);
              res_name[3] = '\0';
              strncpy(res_num,input_line + len + 9,5);
              res_num[5] = '\0';
              chain = input_line[len + 15];

              /* Get the sequence number */
              string_ptr[len] = '\0';
              seq_pos = atoi(input_line);

              /* If this is our residue, then save PDB id and chain */
              if (seq_pos == res_no)
                {
                  /* Save details */
                  strcpy(pdbcode,pdb_code);
                  *chn = chain;

                  /* Set flag that residue found */
                  found = TRUE;
                  printf("FOUND residue %d with PDB: %s %c\n",
                         res_no,pdbcode,chain);
                  fflush(stdout);
                }

              /* If not our residue, then see if it is the closest so
                 far */
              else
                {
                  /* Calculate the distance between this residue and
                     the one we're looking for */
                  dist2 = (seq_pos - res_no) * (seq_pos - res_no);

                  /* If we don't have a closest, yet, or this is the
                     closest so far, save as such */
                  if (closest_pdb_code[0] == '\0' ||
                      dist2 < min_dist2)
                    {
                      /* Save closest residue and PDB details */
                      closest_res_no = seq_pos;
                      strcpy(closest_pdb_code,pdb_code);
                      closest_chain = chain;

                      /* Save the minimum distance so far */
                      min_dist2 = dist2;
                    }
                }
            }
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* If our residue not found, then use the closest residue */
  if (found == FALSE)
    {
      /* If no closest residue (ie no structure) show error message and
         give up */
      if (closest_pdb_code[0] == '\0')
        {
          printf("*** ERROR. No corresponding PDB entry found\n");
          exit(1);
        }

      /* Show that closest residue being used */
      printf("Residue not found. Using closest %d PDB: %s %c\n",
             closest_res_no,closest_pdb_code,closest_chain);
      fflush(stdout);

      /* Save the closest details */
      res_no = closest_res_no;
      strcpy(pdbcode,closest_pdb_code);
      *chn = closest_chain;
    }

  /* Return the closest residue number */
  return(res_no);
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      int seq_no,char aa_code)
{
  char residue_name[4];

  int i, type;

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

  /* Store the amino acid code */
  residue_ptr->aa_code = aa_code;
      
  /* Get the corresponding residue */
  get_res_name(aa_code,residue_ptr->res_name);
  residue_ptr->seq_no = seq_no;

  /* Initialise remaining fields */
  residue_ptr->pdb_code[0] = '\0';
  residue_ptr->res_name[0] = '\0';
  residue_ptr->res_num[0] = '\0';
  residue_ptr->chain = ' ';
  residue_ptr->pdb_res_name[0] = '\0';
  residue_ptr->pdb_res_num[0] = '\0';
  residue_ptr->pdb_chain = ' ';
  residue_ptr->sec_struc = ' ';
  residue_ptr->conservation = 0;
  residue_ptr->last_atom_name[0] = '\0';
  residue_ptr->r_residue_ptr = NULL;
  residue_ptr->first_alist_ptr = NULL;
  residue_ptr->last_alist_ptr = NULL;

  /* If this is the very first, save its pointer */
  if (first_residue_ptr == NULL)
    first_residue_ptr = residue_ptr;

  /* Add link from previous residue to the current one */
  if (last_residue_ptr != NULL)
    last_residue_ptr->next_residue_ptr = residue_ptr;
  last_residue_ptr = residue_ptr;

  /* Initialise pointer to the next residue */
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return new residue pointer */
  return(residue_ptr);
}
/***********************************************************************

get_res_conservation  -  Read in the residue conservation scores from
                         the sas001.conserv file

***********************************************************************/

int get_res_conservation(char *code_dir,struct residue **fst_residue_ptr,
                         float *mx_con,FILE *logfile_ptr)
{
  char input_line[LONG_LINELEN + 1], file_name[FILENAME_LEN + 1];
  char number_string[20];
  char aa_code, chain, res_name[4], res_num[6];

  char *string_ptr, *token[MAX_TOKEN];

  int ipos, iresid, len, line, nresid, nscores, ntokens, seq_no;
  int warned;

  float conservation, max_cons;

  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  *mx_con = max_cons = 0;
  nresid = 0;
  nscores = 0;
  warned = FALSE;

  /* Initialise residue pointers */
  residue_ptr = first_residue_ptr = last_residue_ptr = NULL;

  /* Form the name of the data file */
  strcpy(file_name,code_dir);
  strcat(file_name,"/seq001.conserv");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",
             file_name);
      return(nresid);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LONG_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);
      line++;
      
      /* Extract the information from the line */
      aa_code = input_line[0];
      strncpy(res_name,input_line + 3,3);
      res_name[3] = '\0';
      strncpy(res_num,input_line + 7,5);
      res_num[5] = '\0';
      seq_no = atoi(res_num);
      chain = input_line[13];
      if (chain == ' ')
        chain = 'A';

      /* Get the conservation score */
      strncpy(number_string,input_line + 16,8);
      number_string[8] = '\0';
      conservation = atof(number_string);

      /* Create this residue */
      residue_ptr
        = create_residue_record(&first_residue_ptr,
                                &last_residue_ptr,seq_no,aa_code);

      /* Store the residue details */
      strcpy(residue_ptr->res_name,res_name);
      strcpy(residue_ptr->res_num,res_num);
      residue_ptr->chain = chain;

      /* Store the residue conservation */
      residue_ptr->conservation = conservation;

      /* If highest value so far, then store */
      if (residue_ptr->conservation > max_cons)
        max_cons = residue_ptr->conservation;

      /* Increment count of residues */
      nresid++;
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Return pointer to first residue */
  *fst_residue_ptr = first_residue_ptr;

  /* Show the number of residues created */
  fprintf(logfile_ptr,"Number of residues in sequence:     %8d\n",
          nresid);

  /* Return the maximum conservation score */
  *mx_con = max_cons;

  /* Return number of residues read in */
  return(nresid);
}
/***********************************************************************

get_chain_string  -  Create the string of equivalent chains

***********************************************************************/

int get_chain_string(struct chain *first_chain_ptr,char *equiv_chains,
                     char chain)
{
  char ref_chain;

  int nchains;

  struct chain *chain_ptr, *our_chain_ptr;

  /* Initialise variables */
  nchains = 0;
  ref_chain = chain;

  /* Get the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over the chain equivalences */
  while (chain_ptr != NULL)
    {
      /* If this is our chain, save the pointer */
      if (chain_ptr->chain_id == chain)
        {
          /* If this is a copy, then save the reference chain */
          if (chain_ptr->copy_of != '\0')
            ref_chain = chain_ptr->copy_of;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Put the reference chain in the string */
  equiv_chains[nchains] = ref_chain;
  nchains++;

  /* Get the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over the chain equivalences */
  while (chain_ptr != NULL)
    {
      /* If this is our chain, save the pointer */
      if (chain_ptr->copy_of == ref_chain)
        {
          /* Add to the string */
          equiv_chains[nchains] = chain_ptr->chain_id;
          nchains++;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Close off the chain string */
  equiv_chains[nchains] = '\0';

  /* Return number of equivalenced chains */
  return(nchains);
}
/***********************************************************************

get_pdb_mappings  -  Read in all the PDB mappings from the pdb.align
                     file

***********************************************************************/

struct residue *get_pdb_mappings(struct residue **residue_index_ptr,
                                 int nresid,char *data_dir,
                                 char *pdb_code,char chain,int res_no,
                                 struct chain *first_chain_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char aa_code, align_chain, align_pdb_code[5], res_name[4], res_num[6];
  char equiv_chains[MAX_CHAINS];

  char *string_ptr, *token[MAX_TOKEN];

  int align_start, align_end, len, line, nmatch, ntokens, seq_no;
  int domain, entry_no, iseq, nchains;
  int ok;

  struct residue *residue_ptr, *found_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  nmatch = 0;
  ok = TRUE;

  /* Create the string of equivalent chains */
  nchains = get_chain_string(first_chain_ptr,equiv_chains,chain);

  /* Initialise pointers */
  found_residue_ptr = NULL;

  /* Form the name of the data file */
  strcpy(file_name,data_dir);
  strcat(file_name,"/pdb.align");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",
             file_name);

      /* Return flag indicating we don't have PDB data */
      return(FALSE);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Increment ipos */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);
      
      /* Split the line up into its tokens */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,'\t');

      /* If have right number of tokens, store the data */
      if (ntokens == NALIGN_FIELDS)
        {
          /* Get the PDB code and chain id */
          strcpy(align_pdb_code,token[ALIGN_PDB_CODE]);
          align_chain = token[ALIGN_CHAIN][0];

          /* If this is the right PDB code and chain id, extract the
             alignment data */
          if (!strcmp(pdb_code,align_pdb_code) &&
              strchr(equiv_chains,align_chain) != NULL)
            {
              /* Extract the remaining fields */
              aa_code = token[ALIGN_AA_CODE][0];
              seq_no = atoi(token[ALIGN_SEQ_NO]);
              entry_no = atoi(token[ALIGN_ENTRY_NO]);
              domain = atoi(token[ALIGN_DOMAIN]);
              iseq = atoi(token[ALIGN_ISEQ]);
              strcpy(res_name,token[ALIGN_RES_NAME]);
              to_upper(res_name,3);
              strcpy(res_num,token[ALIGN_RES_NUM]);
              align_start = atoi(token[ALIGN_START]);
              align_end = atoi(token[ALIGN_END]);

              /* If this is a valid sequence, get the corresponding residue
                 and store the details */
              if (seq_no > 0 && seq_no < nresid + 1)
                {
                  /* Get the corresponding residue */
                  residue_ptr = residue_index_ptr[seq_no - 1];

                  /* Store the details from the PDB entry */
                  strcpy(residue_ptr->pdb_res_name,res_name);
                  strcpy(residue_ptr->pdb_res_num,res_num);
                  strcpy(residue_ptr->pdb_code,pdb_code);
                  residue_ptr->pdb_chain = chain;

                  /* If this is our residue, save its pointer */
                  if (seq_no == res_no)
                    found_residue_ptr = residue_ptr;

                  /* Increment count of residues matched */
                  nmatch++;
                }
            }
        }
      
      /* Free up the memory taken by the tokens */
      free_tokens(token,ntokens);
    }

  /* Close the input file */
  fclose(file_ptr);

  /* If no matches, show warning message */
  if (nmatch == 0)
    printf("*** Warning. No residues matched for %s %c in  %s\n",
           pdb_code,chain,file_name);

  /* Return the key residue */
  return(found_residue_ptr);
}
/***********************************************************************

read_database_dat  -  Read in the start residue for the given AlphaFold
                      model from the database.dat file

***********************************************************************/

int read_database_dat(char *pdb_code,char *pdbsum_dir)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char number_string[20];

  char *string_ptr;

  int len, add_res;
  int done;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  add_res = 0;

  /* Form the name of the PDBsum database.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"database.dat");

  /* Open the database.dat file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. No database.dat file for %s",pdb_code);
      return(add_res);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* Check for first residue */
      if (!strncmp(input_line,"AFFRAG_ADD_RES ",15))
        {
          /* Get the residue number */
          strcpy(number_string,input_line + 15);
          add_res = atoi(number_string);

          /* Set flag that done */
          done = TRUE;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the starting residue */
  return(add_res);
}
/***********************************************************************

map_to_pdb  -  Map to the PDB structure just read in

***********************************************************************/

void map_to_pdb(struct residue **residue_index_ptr,int nresid,
                struct r_residue *first_r_residue_ptr,int add_res,
                float max_cons,int res_cons)
{
  int iatom, ires, last_res, natoms, nmapped, resno;
  int found, wanted;

  float conservation, factor;

  struct residue *residue_ptr;
  struct r_atom *r_atom_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  if (max_cons > 0)
    factor = 80.0 / max_cons;
  last_res = -1;
  nmapped = 0;

  /* Get pointer to the first residue */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop over all the residues to initialise */
  while (r_residue_ptr != NULL)
    {
      /* Get number of atoms in this residue */
      natoms = r_residue_ptr->natoms;

      /* Initialise pointer to first atom of this r_residue */
      r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
      iatom = 0;

      /* If starting residue not 1, adjust accordingly */
      if (add_res > 0)
        {
          /* Get the current residue number */
          resno = atoi(r_residue_ptr->res_num) + add_res;
          sprintf(r_residue_ptr->res_num,"%4d ",resno);
        }

      /* Loop through all the residue's atoms to initialise B-value */
      while (iatom < natoms && r_atom_ptr != NULL)
        {
          /* Set B-value field to zero */
          if (res_cons == TRUE)
            r_atom_ptr->bvalue = 0;

          /* Get pointer to next atom in linked-list */
          r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
          iatom++;
        }

      /* Loop until this residue is matched */
      found = FALSE;
      for (ires = last_res + 1; ires < nresid && found == FALSE; ires++)
        {
          /* Get this residue */
          residue_ptr = residue_index_ptr[ires];

          /* If we have a match, store the conservation score in the
             B-value field */
          wanted = FALSE;
          if (residue_ptr->pdb_chain == r_residue_ptr->chain &&
              !strcmp(residue_ptr->pdb_res_num,r_residue_ptr->res_num) &&
              !strcmp(residue_ptr->pdb_res_name,r_residue_ptr->res_name))
            wanted = TRUE;
          else if (residue_ptr->chain == r_residue_ptr->chain &&
                   !strcmp(residue_ptr->res_num,r_residue_ptr->res_num) &&
                   !strcmp(residue_ptr->res_name,r_residue_ptr->res_name))
            wanted = TRUE;

          /* Set flag that residue found */
          if (wanted == TRUE)
            found = TRUE;

          /* If wanted, store the conservation score */
          if (wanted == TRUE && res_cons == TRUE)
            {
              /* Convert the conservation score to a B-value between
                 zero and 80 */
              conservation = 80 * residue_ptr->conservation;

              /* Get number of atoms in this residue */
              natoms = r_residue_ptr->natoms;

              /* Initialise pointer to first atom of this r_residue */
              r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
              iatom = 0;
  
              /* Loop through all the residue's atoms to store B-value */
              while (iatom < natoms && r_atom_ptr != NULL)
                {
                  /* Save the convation score in the B-value field */
                  r_atom_ptr->bvalue = conservation;

                  /* Get pointer to next atom in linked-list */
                  r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
                  iatom++;
                }
            }

          /* If wanted, store the pointer */
          if (wanted == TRUE)
            {
              /* Save the pointer */
              residue_ptr->r_residue_ptr = r_residue_ptr;

              /* Save the last residue matched */
              last_res = ires;
            }
        }

      /* If mapped, increment count */
      if (found == TRUE)
        nmapped++;

      /* Get pointer to the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }

  /* Number of residues mapped */
  printf("Number of residues mapped = %d\n",nmapped);
  fflush(stdout);
}
/***********************************************************************

check_cutoff  -  See if residue is within the cut-off distance of the
                 current residue

***********************************************************************/

int check_cutoff(struct r_residue *r_residue_ptr,
                 struct r_residue *other_r_residue_ptr,double cutoff_dist)
{
  int i, iatom, niatoms, jatom, njatoms;
  int wanted, toofar;

  double cutoff_dist_sqrd, dist_sqrd;

  struct r_atom *r_atom_ptr, *other_r_atom_ptr;

  /* Initialise variables */
  cutoff_dist_sqrd = cutoff_dist * cutoff_dist;
  wanted = FALSE;

  /* Check whether these residue are within reach */
  toofar = FALSE;
  for (i = 0; i < 3 && toofar == FALSE; i++)
    {
      if (r_residue_ptr->coords_min[i]
          - other_r_residue_ptr->coords_max[i] > cutoff_dist)
        toofar = TRUE;
      if (other_r_residue_ptr->coords_min[i]
          - r_residue_ptr->coords_max[i] > cutoff_dist)
        toofar = TRUE;
    }
      
  /* If residues within range calculate all atom-atom distances */
  if (toofar == FALSE)
    {
      /* Get the first of the first residue's atoms */
      r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
      iatom = 0;
      niatoms = r_residue_ptr->natoms;
  
      /* Loop through all the residue's atoms */
      while (r_atom_ptr != NULL && iatom < niatoms && wanted == FALSE)
        {
          /* Check whether this atom within range of the other
             residue */
          toofar = FALSE;
          for (i = 0; i < 3 && toofar == FALSE; i++)
            {
              if (r_atom_ptr->coords[i]
                  - other_r_residue_ptr->coords_max[i] > cutoff_dist ||
                  other_r_residue_ptr->coords_min[i]
                  - r_atom_ptr->coords[i] > cutoff_dist)
                toofar = TRUE;
            }

          /* If atom not too far from the other residue, check its
             distance from all that residue's atoms */
          if (toofar == FALSE)
            {
              /* Get the first of the other residue's atoms */
              other_r_atom_ptr = other_r_residue_ptr->first_r_atom_ptr;
              jatom = 0;
              njatoms = other_r_residue_ptr->natoms;
  
              /* Loop through all the residue's atoms */
              while (jatom < njatoms && other_r_atom_ptr != NULL &&
                     wanted == FALSE)
                {
                  /* Calculate the squared distance between these
                     two atoms */
                  dist_sqrd = 0.0;
                  for (i = 0; i < 3; i++)
                    dist_sqrd
                      = dist_sqrd
                      + (r_atom_ptr->coords[i]
                         - other_r_atom_ptr->coords[i])
                      * (r_atom_ptr->coords[i]
                         - other_r_atom_ptr->coords[i]);

                  /* If within covalent bonding distance, have a
                     bond between the two residues */
                  if (dist_sqrd < cutoff_dist_sqrd)
                    wanted = TRUE;

                  /* Get pointer to next atom in linked-list */
                  other_r_atom_ptr = other_r_atom_ptr->next_r_atom_ptr;
                  jatom++;
                }
            }

          /* Get pointer to next atom in linked-list */
          r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
          iatom++;
        }
    }

  /* Return whether residue is wanted */
  return(wanted);
}
/***********************************************************************

mark_local_residues  -  Remove all but the local residues around our
                        reference residue

***********************************************************************/

int mark_local_residues(struct r_residue *first_r_residue_ptr,
                        struct r_residue *ref_r_residue_ptr,
                        double cutoff_dist)
{
  int nwanted;

  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  nwanted = 0;

  /* Mark the reference residue as wanted */
  ref_r_residue_ptr->wanted = TRUE;

  /* Get pointer to the first residue */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop over all the residues to determine which are close enough to
     the reference to be included */
  while (r_residue_ptr != NULL)
    {
      /* See if this residue is within the cut-off distance */
      if (r_residue_ptr != ref_r_residue_ptr)
        r_residue_ptr->wanted
          = check_cutoff(r_residue_ptr,ref_r_residue_ptr,cutoff_dist);

      /* If residue wanted, then increment count */
      if (r_residue_ptr->wanted == TRUE)
        nwanted++;

      /* Get pointer to the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }

  /* Return the number of wanted residues */
  return(nwanted);
}
/***********************************************************************

delete_residues  -  Delete the unwanted residue records

***********************************************************************/

int delete_residues(struct r_pdb *r_pdb_ptr)
{
  int iresid, nres, nresid, nsaved;

  struct r_molec *r_molec_ptr;
  struct r_residue *r_residue_ptr, *next_r_residue_ptr, *prev_r_residue_ptr;

  /* Initialise variables */
  nres = 0;

  /* Get pointer to the first molecule in the PDB entry */
  r_molec_ptr = r_pdb_ptr->first_r_molec_ptr;

  /* Loop through all the molecules to calc coord limits */
  while (r_molec_ptr != NULL)
    {
      /* Initialise pointers */
      prev_r_residue_ptr = NULL;

      /* Get the first r_residue */
      r_residue_ptr = r_molec_ptr->first_r_residue_ptr;

      r_molec_ptr->first_r_residue_ptr = NULL;
      nresid = r_molec_ptr->nresid;
      iresid = 0;
      nsaved = 0;

      /* Loop through the residue entries to deleted any unwanted ones */
      while (r_residue_ptr != NULL && iresid < nresid)
        {
          /* Get the next residue */
          next_r_residue_ptr = r_residue_ptr->next_r_residue_ptr;

          /* If r_residue not required, then delete */
          if (r_residue_ptr->wanted == FALSE)
            {
              /* If have a previous residue, then get it to point
                 around this one */
              if (prev_r_residue_ptr != NULL)
                prev_r_residue_ptr->next_r_residue_ptr = next_r_residue_ptr;

              /* Otherwise, this must be the first residue, so make
                 the next one the new first */
              else
                r_molec_ptr->first_r_residue_ptr = next_r_residue_ptr;

              /* Free up the memory taken by the deleted r_residue */
              free(r_residue_ptr);
            }

          /* Otherwise, save as the new previous record */
          else
            {
              /* Save current pointer */
              prev_r_residue_ptr = r_residue_ptr;

              /* Increment count of undeleted entries */
              nsaved++;
              nres++;
            }

          /* Get the next residue */
          r_residue_ptr = next_r_residue_ptr;
          iresid++;
        }

      /* Save the number of residues now in the molecule */
      r_molec_ptr->nresid = nsaved;

      /* Get the next molecule */
      r_molec_ptr = r_molec_ptr->next_r_molec_ptr;
    }

  /* Save the new number of residues */
  r_pdb_ptr->nresid = nres;

  /* Return number of remaining entries */
  return(nres);
}
/***********************************************************************

read_seq_fasta  -  Read in the sequence from the original FASTA sequence
                   file

***********************************************************************/

int read_seq_fasta(char *dir,char *fname,char **sname,char **seq)
{
  char file_name[C_FILENAME_LEN], input_line[LONG_LINELEN + 1];
  char protein[STRING_LEN + 1], species[STRING_LEN + 1];
  char uniprot_acc[STRING_LEN + 1], uniprot_id[STRING_LEN + 1];
  char *seq_name, *sequence, *tmp;

  int seq_len, len;

  FILE *file_ptr;

  /* Initialise variables */
  protein[0] = '\0';
  seq_len = 0;
  seq_name = *sname = &null;
  sequence = *seq = &null;
  species[0] = '\0';
  uniprot_acc[0] = '\0';
  uniprot_id[0] = '\0';

  /* Form the name of the sprot.fasta file */
  strcpy(file_name,dir);
  strcat(file_name,"/");
  strcat(file_name,fname);

  /* Open the sequence file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* Show warning */
      printf("*** Warning. File [%s] not found\n",file_name);
      return(seq_len);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LONG_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the line at the last non-blank character */
      len = string_truncate(input_line,LONG_LINELEN);

      /* Check for sequence name */
      if (input_line[0] == '>')
        {
          /* If input line too long, then truncate it */
          if (strlen(input_line) > STRING_LEN)
            input_line[STRING_LEN - 1] = '\0';

          /* Extract the sequence name */
          parse_uniprot_name(input_line,protein,species,uniprot_acc,
                             uniprot_id);
          store_string(&seq_name,protein);
        }

      /* Otherwise, take this to be the sequence, so add to end */
      else if (len > 0)
        {
          /* Get the new sequence length */
          seq_len = strlen(sequence) + len;

          /* Allocate memory for the string */
          tmp = (char *) malloc(sizeof(char)*(seq_len + 1));

          /* Copy the old sequence across */
          strcpy(tmp,sequence);

          /* Add the new sequence */
          strcat(tmp,input_line);

          /* Free up memoryt from previous array */
          if (sequence != &null)
            free(sequence);
          sequence = &null;

          /* Store the current sequence */
          store_string(&sequence,tmp);

          /* Free up temporary memory */
          free(tmp);

          /* Get length of current sequence */
          seq_len = strlen(sequence);
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the sequence and sequence name */
  *seq = sequence;
  *sname = seq_name;

  /* Return the sequence length */
  return(seq_len);
}
/***********************************************************************

initialise_residues  -  Create a residue record for every position in
                        the query sequence

***********************************************************************/

int initialise_residues(struct residue **fst_residue_ptr,
                        char *query_sequence,FILE *logfile_ptr)
{
  char aa_code;

  int ires, nresid, seqlen;

  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  /* Initialise variables */
  nresid = 0;
  seqlen = strlen(query_sequence);

  /* Initialise residue pointers */
  residue_ptr = first_residue_ptr = last_residue_ptr = NULL;

  /* Loop through the sequence */
  for (ires = 0; ires < seqlen; ires++)
    {
      /* Get this residue */
      aa_code = query_sequence[ires];

      /* Increment residue count */
      nresid++;

      /* Create this residue */
      residue_ptr
        = create_residue_record(&first_residue_ptr,&last_residue_ptr,
                                nresid,aa_code);
    }

  /* Return pointer to first residue */
  *fst_residue_ptr = first_residue_ptr;

  /* Show the number of residues created */
  fprintf(logfile_ptr,"Number of residues in sequence:     %8d\n",
          nresid);

  /* Return the number of residues */
  return(nresid);
}
/***********************************************************************

create_residue_index  -  Create an index of the residues in the sequence

***********************************************************************/

int create_residue_index(struct residue *first_residue_ptr,int nresid,
                         struct residue ***res_index_ptr)
{
  int iresid;

  struct residue *residue_ptr;

  struct residue **residue_index_ptr;

  /* Initialise variables */
  iresid = 0;

  /* Create an array of residue pointers */
  residue_index_ptr
    = (struct residue **) malloc ((nresid + 1) * sizeof(struct residue *));
  if (residue_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for residue index\n");
      exit(1);
    }

  /* Get the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all the residues to store their pointers in the array */
  while (residue_ptr != NULL)
    {
      /* Store this residue's pointer */
      if (iresid < nresid)
        {
          /* Store and increment count */
          residue_index_ptr[iresid] = residue_ptr;
          iresid++;
        }

      /* Get the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the index array */
  *res_index_ptr = residue_index_ptr;

  /* Return number of pointers stored */
  return(iresid);
}
/***********************************************************************

create_disease_record  -  Create a new disease record

***********************************************************************/

struct disease *create_disease_record(struct disease **fst_disease_ptr,
                                      struct disease **lst_disease_ptr,
                                      char *id)
{
  int i;

  struct disease *disease_ptr;
  struct disease *first_disease_ptr, *last_disease_ptr;

  /* Initialise disease pointers */
  disease_ptr = NULL;
  first_disease_ptr = *fst_disease_ptr;
  last_disease_ptr = *lst_disease_ptr;

  /* Create disease entry */
  disease_ptr = (struct disease*)malloc(sizeof(struct disease));
  if (disease_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct disease\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  store_string(&(disease_ptr->id),id);

  /* Initialise remaining fields */
  disease_ptr->name = &null;
  disease_ptr->description = &null;
  disease_ptr->omim_id = &null;
  disease_ptr->type = 0;
  disease_ptr->ncodes = 0;
  disease_ptr->nmut = 0;
  disease_ptr->is_note = FALSE;

  /* Initialise pointer to next disease record */
  disease_ptr->next_disease_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of diseases */
  if (first_disease_ptr == NULL)
    first_disease_ptr = disease_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_disease_ptr->next_disease_ptr = disease_ptr;
                      
  /* Save the current disease's pointer */
  last_disease_ptr = disease_ptr;

  /* Return current pointers */
  *fst_disease_ptr = first_disease_ptr;
  *lst_disease_ptr = last_disease_ptr;
  return(disease_ptr);
}
/***********************************************************************

create_desc_entry  -  Create a domain decription record

***********************************************************************/

struct desc *create_desc_entry(struct desc **fst_desc_ptr,
                               struct desc **lst_desc_ptr,char *annot_id)
{
  char *string_ptr;

  struct desc *desc_ptr;
  struct desc *first_desc_ptr, *last_desc_ptr;

  /* Initialise pointers */
  desc_ptr = NULL;
  first_desc_ptr = *fst_desc_ptr;
  last_desc_ptr = *lst_desc_ptr;

  /* Create entry */
  desc_ptr = (struct desc *) malloc(sizeof(struct desc));
  if (desc_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct desc\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  store_string(&(desc_ptr->annot_id),annot_id);

  /* Initialise remaining fields */
  desc_ptr->name = &null;
  desc_ptr->description = &null;
  desc_ptr->type = UNASSIGNED;
  desc_ptr->colour_num = 0;
  desc_ptr->nannot = 0;

  /* Initialise pointer to next entry */
  desc_ptr->next_desc_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_desc_ptr == NULL)
    first_desc_ptr = desc_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_desc_ptr->next_desc_ptr = desc_ptr;
                      
  /* Save the current residue's pointer */
  last_desc_ptr = desc_ptr;

  /* Return current pointers */
  *fst_desc_ptr = first_desc_ptr;
  *lst_desc_ptr = last_desc_ptr;
  return(desc_ptr);
}
/***********************************************************************

find_desc_entry  -  Find the given pfam

***********************************************************************/

struct desc *find_desc_entry(struct desc *first_desc_ptr,char *annot_id)
{
  struct desc *desc_ptr, *found_desc_ptr;

  /* Initialise variables */
  found_desc_ptr = NULL;

  /* Get pointer to the first description record */
  desc_ptr = first_desc_ptr;

  /* Loop over all the pfams to find the one just read in */
  while (desc_ptr != NULL && found_desc_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (!strcmp(annot_id,desc_ptr->annot_id))
        found_desc_ptr = desc_ptr;

      /* Get pointer to the next pfam */
      desc_ptr = desc_ptr->next_desc_ptr;
    }

  /* Return the matching pfam pointer */
  return(found_desc_ptr);
}
/***********************************************************************

create_annot_entry  -  Create an annotation record

***********************************************************************/

struct annot *create_annot_entry(struct annot **fst_annot_ptr,
                                 struct annot **lst_annot_ptr,
                                 struct desc *desc_ptr)
{
  int i;

  struct annot *annot_ptr;
  struct annot *first_annot_ptr, *last_annot_ptr;

  /* Initialise domain pointers */
  annot_ptr = NULL;
  first_annot_ptr = *fst_annot_ptr;
  last_annot_ptr = *lst_annot_ptr;

  /* Create domain entry */
  annot_ptr = (struct annot *) malloc(sizeof(struct annot));
  if (annot_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct annot\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  annot_ptr->desc_ptr = desc_ptr;

  /* Initialise remaining fields */
  annot_ptr->domain = 0;
  annot_ptr->start_res = 0;
  annot_ptr->pre_start = 0;
  annot_ptr->end_res = 0;
  annot_ptr->post_end = 0;
  annot_ptr->annot_no = 0;
  annot_ptr->length = 0;
  annot_ptr->sequence = &null;
  annot_ptr->deleted = FALSE;
  annot_ptr->nrefs = 0;
  annot_ptr->serial_no = 0;
  annot_ptr->aamut_ptr = NULL;

  /* Initialise pointer to next entry */
  annot_ptr->next_annot_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_annot_ptr == NULL)
    first_annot_ptr = annot_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_annot_ptr->next_annot_ptr = annot_ptr;
                      
  /* Save the current residue's pointer */
  last_annot_ptr = annot_ptr;

  /* Return current pointers */
  *fst_annot_ptr = first_annot_ptr;
  *lst_annot_ptr = last_annot_ptr;
  return(annot_ptr);
}
/***********************************************************************

create_alist_record  -  Create a new alist record

***********************************************************************/

struct alist *create_alist_record(struct alist **fst_alist_ptr,
                                  struct alist **lst_alist_ptr,
                                  struct annot *annot_ptr)
{
  struct alist *alist_ptr;
  struct alist *first_alist_ptr, *last_alist_ptr;

  /* Initialise alist pointers */
  alist_ptr = NULL;
  first_alist_ptr = *fst_alist_ptr;
  last_alist_ptr = *lst_alist_ptr;

  /* Create alist entry */
  alist_ptr = (struct alist*)malloc(sizeof(struct alist));
  if (alist_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct alist\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  alist_ptr->annot_ptr = annot_ptr;

  /* Initialise pointer to next alist record */
  alist_ptr->next_alist_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of alisted alist */
  if (first_alist_ptr == NULL)
    first_alist_ptr = alist_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_alist_ptr->next_alist_ptr = alist_ptr;
                      
  /* Save the current alist's pointer */
  last_alist_ptr = alist_ptr;

  /* Return current pointers */
  *fst_alist_ptr = first_alist_ptr;
  *lst_alist_ptr = last_alist_ptr;
  return(alist_ptr);
}
/***********************************************************************

create_aamut_record  -  Create a new aamut record

***********************************************************************/

struct aamut *create_aamut_record(struct aamut **fst_aamut_ptr,
                                  struct aamut **lst_aamut_ptr,
                                  int seq_no,struct annot *annot_ptr)
{
  int i;

  struct aamut *aamut_ptr;
  struct aamut *first_aamut_ptr, *last_aamut_ptr;

  /* Initialise aamut pointers */
  aamut_ptr = NULL;
  first_aamut_ptr = *fst_aamut_ptr;
  last_aamut_ptr = *lst_aamut_ptr;

  /* Create aamut entry */
  aamut_ptr = (struct aamut*)malloc(sizeof(struct aamut));
  if (aamut_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct aamut\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  aamut_ptr->seq_no = seq_no;
  aamut_ptr->annot_ptr = annot_ptr;

  /* Initialise remaining records */
  aamut_ptr->aa_code = '\0';
  aamut_ptr->aa_mut = '\0';
  aamut_ptr->nwritten = 0;
  aamut_ptr->ndiseases = 0;
  aamut_ptr->var_type = -1;
  aamut_ptr->aa_extra = &null;
  aamut_ptr->mut_extra = &null;
  aamut_ptr->mut_length = 0;
  aamut_ptr->disease_ptr = NULL;

  /* Initialise pointer to next aamut record */
  aamut_ptr->next_aamut_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of aamuted aamut */
  if (first_aamut_ptr == NULL)
    first_aamut_ptr = aamut_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_aamut_ptr->next_aamut_ptr = aamut_ptr;
                      
  /* Save the current aamut's pointer */
  last_aamut_ptr = aamut_ptr;

  /* Return current pointers */
  *fst_aamut_ptr = first_aamut_ptr;
  *lst_aamut_ptr = last_aamut_ptr;
  return(aamut_ptr);
}
/***********************************************************************

find_disease_record  -  Find the given script entry

***********************************************************************/

struct disease *find_disease_record(struct disease *first_disease_ptr,
                                    char *id)
{
  struct disease *disease_ptr, *found_disease_ptr;

  /* Initialise variables */
  found_disease_ptr = NULL;
  if (id[0] == '\0')
    return(found_disease_ptr);

  /* Get pointer to the first script */
  disease_ptr = first_disease_ptr;

  /* Loop over all the scripts to find the one just read in */
  while (disease_ptr != NULL && found_disease_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (!strcmp(disease_ptr->id,id))
        found_disease_ptr = disease_ptr;

      /* Get pointer to the next script */
      disease_ptr = disease_ptr->next_disease_ptr;
    }

  /* Return the matching script pointer */
  return(found_disease_ptr);
}
/***********************************************************************

get_mutations  -  Read in the mutation details from the uprotein.dat
                  file

***********************************************************************/

int get_mutations(char *uniprot_acc,struct residue **residue_index_ptr,
                  int nresid,char *code_dir,FILE *logfile_ptr)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN + 1];
  char disease_id[LINELEN + 1], number_string[20];
  char aa_code, aa_mut, annot_id[20], res_name[4], mut_name[4];

  char *string_ptr;

  int annot_no, error, len, line, mut_type, res_no, seq_no;
  int nannot, ndesc, ndisease, nmut, nvar_res, var_type;
  int done;

  long offset;

  struct aamut *aamut_ptr, *first_aamut_ptr, *last_aamut_ptr;
  struct alist *alist_ptr;
  struct annot *annot_ptr, *first_annot_ptr, *last_annot_ptr;
  struct desc *desc_ptr, *first_desc_ptr, *last_desc_ptr;
  struct disease *disease_ptr, *first_disease_ptr, *last_disease_ptr;
  struct dlist *dlist_ptr;
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  aa_code = aa_mut = res_name[0] = mut_name[0] = '\0';
  annot_id[0] = '\0';
  annot_no = 0;
  done = FALSE;
  mut_type = UNASSIGNED;
  nannot = 0;
  ndesc = 0;
  ndisease = 0;
  nvar_res = 0;
  nmut = 0;
  offset = -1;
  seq_no = -1;
  res_no = -1;
  var_type = -1;

  /* Initialise pointers */
  residue_ptr = NULL;
  aamut_ptr = first_aamut_ptr = last_aamut_ptr = NULL;
  annot_ptr = first_annot_ptr = last_annot_ptr = NULL;
  desc_ptr = first_desc_ptr = last_desc_ptr = NULL;
  disease_ptr = first_disease_ptr = last_disease_ptr = NULL;

  /* Form the name of the data file */
  strcpy(file_name,code_dir);
  strcat(file_name,"/uprotein.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",
             file_name);
      done = TRUE;
    }

  /* Loop while reading in records from the file */
  while (done == FALSE && fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* Get the blank character after the field name */
      string_ptr = strchr(input_line,' ');

      /* Get the variant type */
      if (!strncmp(input_line,"VARIANT_NAME[",13))
        {
          /* Initialise variant type */
          var_type = -1;

          /* Get the variant type from the first subscript */
          string_ptr = strchr(input_line,']');
          if (string_ptr != NULL)
            {
              /* Extract the number */
              string_ptr[0] = '\0';
              strcpy(number_string,input_line + 13);
              var_type = atoi(number_string) - 1;

              /* Check we have a valid number */
              if (var_type < 0 || var_type > NVAR_TYPES - 1)
                var_type = -1;
            }

          /* Initialise disease record */
          disease_ptr = NULL;
        }

      /* If this is the disease id, then store */
      else if (!strncmp(input_line,"DISEASE_ID[",11))
        {
          /* Extract the disease id */
          if (string_ptr != NULL)
            {
              /* Get the id */
              strcpy(disease_id,string_ptr + 1);

              /* See if we already have this disease */
              disease_ptr
                = find_disease_record(first_disease_ptr,disease_id);

              /* If not found, create a new record */
              if (disease_ptr == NULL)
                {
                  /* Create a disease record */
                  disease_ptr
                    = create_disease_record(&first_disease_ptr,
                                            &last_disease_ptr,disease_id);

                  /* Determine whether this is a note */
                  if (var_type == DISEASE)
                    disease_ptr->is_note = FALSE;
                  else
                    disease_ptr->is_note = TRUE;

                  /* Increment count of diseases */
                  ndisease++;
                }
            }
        }

      /* If this is the disease name, then store */
      else if (!strncmp(input_line,"DISEASE_NAME[",13))
        {
          /* Store the name */
          if (string_ptr != NULL && disease_ptr != NULL)
            store_string(&(disease_ptr->name),string_ptr);
        }

      /* Ditto for disease description */
      else if (!strncmp(input_line,"DISEASE_DESCRIPTION[",20))
        {
          /* Store the name */
          if (string_ptr != NULL && disease_ptr != NULL)
            store_string(&(disease_ptr->description),string_ptr);
        }

      /* OMIM id */
      else if (!strncmp(input_line,"OMIM_ID[",8))
        {
          /* Store the name */
          if (string_ptr != NULL && disease_ptr != NULL)
            store_string(&(disease_ptr->omim_id),string_ptr);
        }

      /* Amino acid code */
      else if (!strncmp(input_line,"AA_CODE[",8) && string_ptr != NULL)
        aa_code = string_ptr[1];

      /* Residue name */
      else if (!strncmp(input_line,"RES_NAME[",9) && string_ptr != NULL)
        strcpy(res_name,string_ptr + 1);

      /* If this is the id, create annotation record */
      else if (!strncmp(input_line,"ANNOT_ID[",9))
        {
          /* Get the annot id */
          annot_id[0] = '\0';
          string_ptr = strchr(input_line,' ');
          if (string_ptr != NULL)
            strcpy(annot_id,string_ptr + 1);
 
          /* If not found, create a desc entry */
          if (desc_ptr == NULL)
            {
              /* Create a annot decription record */
              desc_ptr
                = create_desc_entry(&first_desc_ptr,&last_desc_ptr,
                                    annot_id);

              /* Increment count of annot descriptions */
              ndesc++;
            }
        }

      /* Mutation code */
      else if (!strncmp(input_line,"AA_MUT[",7) && string_ptr != NULL &&
               aamut_ptr != NULL)
        {
          /* Get the variants */
          aa_mut = string_ptr[1];

          /* Store it */
          aamut_ptr->aa_mut = aa_mut;
        }

      /* Mutation name */
      else if (!strncmp(input_line,"MUT_NAME[",9) && string_ptr != NULL)
        strcpy(mut_name,string_ptr + 1);

      /* If sequence number, then save mutation details */
      else if (!strncmp(input_line,"SEQ_NO[",7) && string_ptr != NULL)
        {
          /* Initialise pointers */
          residue_ptr = NULL;
          aamut_ptr = NULL;
          annot_ptr = NULL;

          /* Get the residue number */
          if (string_ptr != NULL)
            {
              /* Get the residue number */
              seq_no = atoi(string_ptr + 1);
              res_no = seq_no - 1;

              /* Get the corresponding residue record */
              if (res_no > -1 && res_no < nresid)
                {
                  /* Get the corresponding residue record */
                  residue_ptr = residue_index_ptr[res_no];

                  /* If invalid, then unset residue number */
                  if (residue_ptr == NULL)
                    res_no = -1;
                }

              /* Show warning */
              else
                printf("*** Warning. Residue number off end of sequence "
                       "%d (%d)\n",res_no,nresid);
            }
          else
            res_no = -1;

          /* If have a valid residue, create the mutation record
             for it */
          if (res_no > -1)
            {
              /* See if we have this annot already */
              desc_ptr = find_desc_entry(first_desc_ptr,annot_id);

              /* Create an annotation record */
              annot_ptr
                = create_annot_entry(&first_annot_ptr,&last_annot_ptr,
                                     desc_ptr);

              /* Inrement count of annotations */
              nannot++;

              /* If this is this residue's first variant, increment
                 count of residues with variants */
              if (residue_ptr->first_alist_ptr == NULL)
                nvar_res++;

              /* Create a alist record to attach this annotation to this
                 residue */
              alist_ptr
                = create_alist_record(&(residue_ptr->first_alist_ptr),
                                      &(residue_ptr->last_alist_ptr),
                                      annot_ptr);

              /* Create a mutation record */
              aamut_ptr
                = create_aamut_record(&first_aamut_ptr,&last_aamut_ptr,
                                      seq_no,annot_ptr);

              /* Attach mutation to the annotation record */
              annot_ptr->aamut_ptr = aamut_ptr;

              /* Save the known mutation details */
              aamut_ptr->aa_code = aa_code;
              aamut_ptr->var_type = var_type;
              aamut_ptr->disease_ptr = disease_ptr;

              /* Increment count of mutations */
              nmut++;
            }
        }

      /* If have hit the end, we are done */
      else if (!strncmp(input_line,"NVAR_TYPES ",11))
        done = TRUE;
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Show number of mutations read in */
  fprintf(logfile_ptr,"Number of UniProt mutations read in:%8d\n",nmut);

  /* Return count of residues with variants read in */
  return(nvar_res);
}
/***********************************************************************

added_variants  -  Add any added variants then to the relevant residues

***********************************************************************/

int added_variants(struct residue **residue_index_ptr,int nresid,
                   int nvar_res,struct parameters params)
{
  int res_no, seq_no;

  struct aamut *aamut_ptr, *first_aamut_ptr, *last_aamut_ptr;
  struct addmut *addmut_ptr, *last_addmut_ptr;
  struct alist *alist_ptr;
  struct annot *annot_ptr, *first_annot_ptr, *last_annot_ptr;
  struct desc *desc_ptr, *first_desc_ptr, *last_desc_ptr;
  struct residue *residue_ptr;

  /* Initialise pointers */
  residue_ptr = NULL;
  aamut_ptr = first_aamut_ptr = last_aamut_ptr = NULL;
  annot_ptr = first_annot_ptr = last_annot_ptr = NULL;
  desc_ptr = first_desc_ptr = last_desc_ptr = NULL;

  /* If we don't have any added variants, then return */
  if (params.naddmut == 0)
    return(nvar_res);

  /* Create a desc entry */
  desc_ptr
    = create_desc_entry(&first_desc_ptr,&last_desc_ptr,"Added variant");

  /* Get the first added variant */
  addmut_ptr = params.first_addmut_ptr;

  /* Loop over all the variants */
  while (addmut_ptr != NULL)
    {
      /* Initialise */
      residue_ptr = NULL;

      /* Get the residue number */
      seq_no = addmut_ptr->seq_no;
      res_no = seq_no - 1;

      /* Get the corresponding residue record */
      if (res_no > -1 && res_no < nresid)
        {
          /* Get the corresponding residue record */
          residue_ptr = residue_index_ptr[res_no];

          /* If invalid, then unset residue number */
          if (residue_ptr == NULL)
            res_no = -1;
        }

      /* If outside valid range, show warning */
      else
        {
          /* Show warning */
          printf("*** Warning. Residue number off end of sequence "
                 "%d (%d)\n",res_no,nresid);

          /* Unset residue number */
          res_no = -1;
        }

      /* If residue located, save this variant */
      if (residue_ptr != NULL)
        {
          /* Create an annotation record */
          annot_ptr
            = create_annot_entry(&first_annot_ptr,&last_annot_ptr,
                                 desc_ptr);

          /* If this is this residue's first variant, increment
             count of residues with variants */
          if (residue_ptr->first_alist_ptr == NULL)
            nvar_res++;

          /* Create a alist record to attach this annotation to this
             residue */
          alist_ptr
            = create_alist_record(&(residue_ptr->first_alist_ptr),
                                  &(residue_ptr->last_alist_ptr),
                                  annot_ptr);

          /* Create a mutation record */
          aamut_ptr
            = create_aamut_record(&first_aamut_ptr,&last_aamut_ptr,
                                  seq_no,annot_ptr);

          /* Attach mutation to the annotation record */
          annot_ptr->aamut_ptr = aamut_ptr;

          /* Save the known mutation details */
          aamut_ptr->aa_code = addmut_ptr->aa_code;
          aamut_ptr->aa_mut = addmut_ptr->aa_mut;
          aamut_ptr->var_type = ADDED_VAR;
          aamut_ptr->aa_extra = addmut_ptr->aa_extra;
          aamut_ptr->mut_extra = addmut_ptr->mut_extra;
          aamut_ptr->mut_length = addmut_ptr->mut_length;
        }

      /* Get the next variant */
      addmut_ptr = addmut_ptr->next_addmut_ptr;
    }

  /* Return the number of residues with variants */
  return(nvar_res);
}
/***********************************************************************

variant_sort  -  Sort the variants by sequence number using a simple
                  bubble sort

***********************************************************************/

void variant_sort(struct alist **fst_alist_ptr)
{
  int i, next_type, swap, swapped, type;

  struct alist *alist_ptr, *first_alist_ptr;
  struct alist *next_alist_ptr, *prev_alist_ptr, *swap_alist_ptr;
  struct annot *annot_ptr, *next_annot_ptr;

  /* Get pointer to the first variant record */
  first_alist_ptr = *fst_alist_ptr;
  swapped = TRUE;

  /* Loop until Pfam domains in required order */
  while (swapped == TRUE)
    {
      /* Initialise swapped flag */
      swapped = FALSE;
          
      /* Get pointer to first variant record */
      alist_ptr = first_alist_ptr;
      prev_alist_ptr = NULL;

      /* Loop over all the variants, adjusting their order where
         necessary */
      while (alist_ptr != NULL)
        {
          /* Get pointer to the next record */
          next_alist_ptr = alist_ptr->next_alist_ptr;

          /* Determine whether current two records are in the right
             order */
          swap = FALSE;
          if (next_alist_ptr != NULL)
            {
              /* Get the corresponding annotation records */
              annot_ptr = alist_ptr->annot_ptr;
              next_annot_ptr = next_alist_ptr->annot_ptr;

              /* Get the variant types */
              type = annot_ptr->aamut_ptr->var_type + 1;
              next_type = next_annot_ptr->aamut_ptr->var_type + 1;

              /* Check seq_no order */
              if (VAR_PRIORITY[type] < VAR_PRIORITY[next_type])
                swap = TRUE;
            }

          /* If need to swap, perform the swap */
          if (swap == TRUE)
            {
              /* If this is the first variant record, adjust first pointer */
              if (alist_ptr == first_alist_ptr)
                first_alist_ptr = next_alist_ptr;

              /* Otherwise, get the previous variant pointer to point
                 to next one */
              else
                prev_alist_ptr->next_alist_ptr = next_alist_ptr;

              /* Swap pointers within the records */
              alist_ptr->next_alist_ptr
                = next_alist_ptr->next_alist_ptr;
              next_alist_ptr->next_alist_ptr = alist_ptr;

              /* Swap the pointers themselves */
              swap_alist_ptr = alist_ptr;
              alist_ptr = next_alist_ptr;
              next_alist_ptr = swap_alist_ptr;

              /* Set flag that swap made */
              swapped = TRUE;
            }

          /* Save current pointer */
          prev_alist_ptr = alist_ptr;

          /* Get the next variant record */
          alist_ptr = next_alist_ptr;
        }

      /* Save the first variant pointer */
      *fst_alist_ptr = first_alist_ptr;
    }
}
/***********************************************************************

sort_variants  -  Sort each residue's variants by importance

***********************************************************************/

void sort_variants(struct residue **residue_index_ptr,int nresid)
{
  int iresid;
  
  struct residue *residue_ptr;

  /* Loop over all the residues in the sequence */
  for (iresid = 0; iresid < nresid; iresid++)
    {
      /* Get this residue */
      residue_ptr = residue_index_ptr[iresid];

      /* If it has any variants, sort them */
      if (residue_ptr->first_alist_ptr != NULL)
        variant_sort(&(residue_ptr->first_alist_ptr));
    }
}
/***********************************************************************

get_residue_details  -  Interpret the residue details from the CORA-file
                        segment

***********************************************************************/

int get_residue_details(char *residue_details,int version,char *aa_code,
                        char *res_name,char *res_num,
                        char *similarity_marker,int *similarity,
                        char *sec_struc)
{
  char number_string[20];

  int a, aa_num, number, z;
  int blank;

  /* Initialise variables */
  a = 'a';
  z = 'z';
  blank = FALSE;
  *aa_code = '0';
  *similarity = 0;
  *sec_struc = ' ';

  /* Get the residue name, number and chain-id from the
     residue details */
  if (version == 0)
    *aa_code = residue_details[6];
  else
    *aa_code = residue_details[7];

  /* Check for lower-case letters (which represent
     cysteines involved in disulphide bonds) */
  aa_num = *aa_code;
  if (aa_num >= a && aa_num <= z)
    *aa_code = 'C';

  /* If this is a gap, identify as a blank record */
  if (*aa_code == '0')
    blank = TRUE;

  /* Update this residue's details */
  else
    {
      /* Get the corresponding three-letter code */
      get_res_name(*aa_code,res_name);

      /* Get the residue number and chain id */
      if (version == 0)
        strncpy(res_num,residue_details,5);
      else
        strncpy(res_num,residue_details+1,5);
      res_num[5] = '\0';

      /* Get the similarity and FASTA similarity marker */
      if (version == 0)
        strncpy(number_string,residue_details+7,2);
      else
        strncpy(number_string,residue_details+8,2);
      number_string[2] = '\0';
      number = atoi(number_string);

      /* Extract the similarity and FASTA similarity marker
         from the number */
      if (number >= 60)
        {
          *similarity_marker = 'X';
          *similarity = number - 60;
        }
      else if (number >= 40)
        {
          *similarity_marker = ':';
          *similarity = number - 40;
        }
      else if (number >= 20)
        {
          *similarity_marker = '.';
          *similarity = number - 20;
        }
      else
        {
          *similarity_marker = ' ';
          *similarity = number;
        }

      /* Get the secondary structure assignment */
      if (version == 0)
        *sec_struc = residue_details[9];
      else
        *sec_struc = residue_details[10];
      if (*sec_struc == 'G')
        *sec_struc = 'H';
      if (*sec_struc == '.')
        *sec_struc = ' ';
    }

  /* Return whether residue corresponds to a blank position */
  return(blank);
}
/***********************************************************************

init_residues  -  Initialise all the PDB residue variant flags

***********************************************************************/

void init_residues(struct r_residue *first_r_residue_ptr)
{
  struct r_residue *r_residue_ptr;

  /* Get pointer to the first residue */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop over all the residues to initialise */
  while (r_residue_ptr != NULL)
    {
      /* Initialise */
      r_residue_ptr->wanted = NO_VARIANT;

      /* Get pointer to the next residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
    }
}
/***********************************************************************

get_seq_start  -  Get the starting sequence position for this domain

***********************************************************************/

int get_seq_start(char *dir,int start,int end)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char tmp[LINELEN + 1];

  char *string_ptr;

  int ipos, len, seq_end, seq_start, seq_len;
  int done, overlap;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  seq_end = seq_start = seq_len = 0;

  /* Form the name of the align.head file */
  strcpy(file_name,dir);
  strcat(file_name,"/align.head");

  /* Open the file to see if it exists */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(seq_start);

  /* Loop through the file, line by line */
  while (done == FALSE &&
         fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_chomp(input_line);

      /* If this is the name line, check if have whole sequence */
      if (!strncmp(input_line,"NAME: ",6) &&
          (string_ptr = strstr(input_line,"Whole sequence")) != NULL)
        {
          /* Start residue if the first residue */
          seq_start = 1;
          done = TRUE;

          /* Get the length of the sequence */
          strcpy(tmp,string_ptr + 17);
          string_ptr = strchr(tmp,' ');
          if (string_ptr != NULL)
            string_ptr[0] = '\0';
          seq_len = atoi(tmp);

          /* Get the last residue */
          seq_end = seq_start + seq_len - 1;
          string_ptr = NULL;
        }

      /* Otherwise, get the residue range */
      else if (!strncmp(input_line,"NAME: ",6))
        {
          /* If an unassigned domain, gets start and end residues */
          if (strstr(input_line,"Unassigned") != NULL)
            {
              /* Get the start residue */
              strcpy(tmp,input_line + 18);
              string_ptr = strchr(tmp,' ');
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              seq_start = atoi(tmp);

              /* Find start of last residue */
              string_ptr = strstr(input_line + 18," to ");
            }

          /* Otherwise, take to be a named domain */
          else
            {
              /* Find start of seq range */
              string_ptr = strstr(input_line," Seq ");
              if (string_ptr != NULL)
                {
                  /* Save this position */
                  ipos = string_ptr - input_line;
                  
                  /* Get the first residue */
                  strcpy(tmp,string_ptr + 5);
                  string_ptr = strchr(tmp,' ');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  seq_start = atoi(tmp);
              
                  /* Find start of last residue */
                  string_ptr = strstr(input_line + ipos + 5," to ");
                }
            }
             
          /* Get the length of the sequence */
          if (string_ptr != NULL)
            {
              /* Get the last residue */
              strcpy(tmp,string_ptr + 4);
              string_ptr = strchr(tmp,' ');
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              seq_end = atoi(tmp);
              
              /* Get the length */
              seq_len = seq_end - seq_start + 1;
            }

          /* Set flag that we're done */
          done = TRUE;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* See if the range overlaps our range of interest */
  overlap = TRUE;
  if (seq_start > end || seq_end < start)
    overlap = FALSE;

  /* If no overlap, then don't want to read the corresponding .aln file */
  if (overlap == FALSE)
    seq_start = 0;

  /* Return the first residue */
  return(seq_start);
}
/***********************************************************************

map_residues  -  Scan the align.aln files to map the residues onto our
                 PDB structure

***********************************************************************/

int map_residues(char *dir_name,struct residue **residue_index_ptr,
                 int nresid,struct r_pdb *r_pdb_ptr,char chain,
                 int nvar_res,int start, int end)
{
  char domain_dir[FILENAME_LEN + 1];
  char file_name[FILENAME_LEN + 1], input_line[LONG_LINELEN + 1];
  char aa_code[2], domain_id[10], res_name[2][4], res_num[2][6];
  char residue_details[12], sec_struc, similarity_marker, stype;
  char *string_ptr, *token[MAX_TOKENS];

  int align_length, idom, itoken, len, line, nmapped, nseqs, ntokens;
  int ipos, ires, iresid, iseq, loop, nres, seq_pos, similarity;
  int column_width, next_resid, start_pos, version;
  int seq_start, start_res, end_res;
  int blank[2], check, done, have_file, save, wanted;

  int *mapped;

  struct residue *residue_ptr;
  struct r_residue *r_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  check = TRUE;
  done = FALSE;
  have_file = FALSE;
  idom = 0;
  nmapped = 0;
  nseqs = 0;
  save = TRUE;
  seq_start = 0;
  start_res = end_res = 0;

  /* Create array with flags recording which residues have been mapped */
  nres = end - start + 1;
  mapped = (int *) malloc(nres * sizeof(int));
  for (ires = 0; ires < nres; ires++)
    mapped[ires] = FALSE;

  /* Initialise all the PDB residue variant flags */
  init_residues(r_pdb_ptr->first_r_residue_ptr);

  /* Loop through all the possible domain directories */
  for (idom = 0; idom < MAX_DOM && done == FALSE; idom++)
    {
      /* Re-intialise */
      line = 0;
      seq_pos = -1;

      /* Initialise CORA-format parameters */
      column_width = 10;
      start_pos = 13;
      version = 0;

      /* Form the domain identifier */
      sprintf(domain_id,"%3.3d",idom);

      /* Form the name of the domain directory */
      strcpy(domain_dir,dir_name);
      strcat(domain_dir,"/domains/domain");
      strcat(domain_dir,domain_id);

      /* Get the starting sequence position for this domain */
      seq_start = get_seq_start(domain_dir,start,end);

      /* If seq-start is zero, then don't need to read the corresponding
         .aln file */
      if (seq_start > 0)
        {
          /* Form the name of the corresponding align.aln file */
          strcpy(file_name,domain_dir);
          strcat(file_name,"/align.aln");

          /* Initialise flag */
          have_file = TRUE;
          wanted = TRUE;

          /* Open the file to see if it exists */
          if ((file_ptr = fopen(file_name,"r")) ==  NULL)
            {
              /* Set flags */
              have_file = FALSE;
              done = TRUE;
              wanted = FALSE;
            }
        }
      else
        have_file = wanted = FALSE;

      /* Initialise alignment's start and end residues */
      start_res = end_res = 0;

      /* Loop through the file, line by line */
      while (done == FALSE && wanted == TRUE && 
             fgets(input_line,LONG_LINELEN,file_ptr) != NULL)
        {
          /* Truncate the input line */
          len = string_chomp(input_line);

          /* Check for new CORA version */
          if (!strncmp(input_line,"#FM CORA_FORMAT",15))
            {
              /* Set CORA-format parameters */
              column_width = 11;
              start_pos = 14;
              version = 1;
            }

          /* Process only if this is not a comment line */
          else if (input_line[0] != '#')
            {
              /* If this is the second record, read in the PDB codes and
                 check they match our sequence records  */
              if (line == 1)
                {
                  /* Pick up all the tokens from the input line */
                  ntokens = tokenize(input_line,token,MAX_TOKENS,' ');
                  nseqs = ntokens;

                  /* Check that last token isn't empty */
                  if (strlen(token[ntokens - 1]) < 6)
                    nseqs--;

                  /* Initialise flag */
                  wanted = FALSE;

                  /* Loop over the tokens to see if our PDB chain is
                     amongst them */
                  for (itoken = 0; itoken < ntokens; itoken++)
                    {
                      /* Check for a match */
                      if (!strncmp(token[itoken],r_pdb_ptr->pdb_code,4) &&
                          token[itoken][4] == chain)
                        {
                          /* Set flag that this file is wanted */
                          wanted = TRUE;

                          /* Save the column of the match */
                          seq_pos = itoken;
                       }
                    }
                }

              /* If this is the third record, read in the total number of 
                 residues in the alignment */
              else if (line == 2)
                align_length = atoi(input_line);

              /* Otherwise, get the alignment at this residue position */
              else if (line > 2)
                {
                  /* Initialise the sequence position */
                  iseq = 0;

                  /* Loop twice to process the query residue and our PDB
                     residue */
                  for (loop = 0; loop < 2; loop++)
                    {
                      /* Get the start position for this residue */
                      ipos = start_pos + iseq * column_width;

                      /* Extract the residue details for this residue */
                      strncpy(residue_details,input_line+ipos,column_width);
                      residue_details[column_width] = '\0';

                      /* Interpret the residue details */
                      blank[loop]
                        = get_residue_details(residue_details,version,
                                              &(aa_code[loop]),
                                              res_name[loop],
                                              res_num[loop],
                                              &similarity_marker,&similarity,
                                              &sec_struc);

                      /* Set the sequence to the PDB sequence */
                      iseq = seq_pos;
                    }
                  
                  /* If neither position is blank, save the mapping */
                  if (blank[0] == FALSE && blank[1] == FALSE)
                    {
                      /* Get the residue number */
                      iresid = atoi(res_num[0]) + seq_start - 1;

                      /* Update start and end residues */
                      if (start_res == 0)
                        start_res = iresid;
                      end_res = iresid;

                      /* Get the corresponding residue in the UniProt
                         sequence */
                      residue_ptr = residue_index_ptr[iresid - 1];

                      /* Convert residue name to upper case */
                      to_upper(res_name[1],3);

                      /* Find the corresponding PDB residue */
                      r_residue_ptr
                        = r_find_residue(r_pdb_ptr->first_r_residue_ptr,
                                         res_name[1],res_num[1],
                                         chain);

                      /* If we have both residues, then save the mapping */
                      if (residue_ptr != NULL && r_residue_ptr != NULL)
                        {
                          /* Save the mapping */
                          residue_ptr->r_residue_ptr = r_residue_ptr;

                          /* Work out the relative position of this
                             residue in the range to be mapped */
                          ires = iresid - start + 1;

                          /* If not already mapped, then set flag */
                          if (ires > -1 && ires < nres &&
                              mapped[ires] == FALSE)
                            {
                              /* Set flag */
                              mapped[ires] = TRUE;

                              /* Increment count of mapped residues */
                              nmapped++;
                            }
                        }
                    }
                }

              /* Increment line count */
              line++;
            }
        }

      /* Close the input file */
      if (have_file == TRUE)
        fclose(file_ptr);

      /* If we have mapped all the residues we needed to, can end here */
      if (nmapped == nres)
        done = TRUE;
    }

  /* Return the number of variants mapped */
  return(nmapped);
}
/***********************************************************************

get_last_atom  -  Get this residue's last-defined atom

***********************************************************************/

int get_last_atom(struct r_residue *r_residue_ptr,char *last_atom_name)
{
  int iatom, natoms;

  struct r_atom *r_atom_ptr;

  /* Initialise variables */
  last_atom_name[0] = '\0';

  /* Get number of atoms */
  natoms = r_residue_ptr->natoms;

  /* Initialise pointer to first atom of this r_residue */
  r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
  iatom = 0;

  /* Loop through all the residue's atoms and create corresponding
     atom records for our residue */
  while (iatom < natoms && r_atom_ptr != NULL)
    {
      /* Save the name of this atom, if not a hydrogen */
      if (strcmp(r_atom_ptr->element," H"))
        strcpy(last_atom_name,r_atom_ptr->atom_name);

      /* Get pointer to next atom in linked-list */
      r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
      iatom++;
    }

  /* Return last atom */
  return(iatom);
}
/***********************************************************************

write_variants  -  Write the variants to the variants.dat and
                   variants.map files

***********************************************************************/

void write_variants(struct residue **residue_index_ptr,int nresid,
                    struct parameters params)
{
  char file_name[FILENAME_LEN + 1], last_atom_name[5], top_mut;
  
  static char *variant_type[NVAR_TYPES] = { 
    "disease variant", "disease note", "mutagenesis expt",
    "DDD variant", "Added variant"
  };

  int i, iresid, max_type, nmut_type[NVAR_TYPES], nres, nvar, type;
  int added_var, have_added;
  
  struct addmut *addmut_ptr;
  struct aamut *aamut_ptr;
  struct annot *annot_ptr;
  struct alist *alist_ptr;
  struct disease *disease_ptr;
  struct residue *residue_ptr;
  struct r_residue *r_residue_ptr;

  FILE *datfile_ptr, *mapfile_ptr;

  /* Initialise variants */
  added_var = FALSE;
  nres = nvar = 0;
  for (i = 0; i < NVAR_TYPES; i++)
    nmut_type[i] = 0;

  /* Form name of the output .map file */
  strcpy(file_name,"variants.map");

  /* Open the output file */
  mapfile_ptr = open_output_file(file_name,TRUE);
  if (mapfile_ptr == NULL)
    {
      printf("*** ERROR. Unable to open output file %s\n",file_name);
      exit(1);
    }

  /* Form name of the output .dat file */
  strcpy(file_name,"variants.dat");

  /* Open the output file */
  datfile_ptr = open_output_file(file_name,TRUE);
  if (datfile_ptr == NULL)
    {
      printf("*** ERROR. Unable to open output file %s\n",file_name);
      exit(1);
    }

  /* Loop over all the residues in the sequence */
  for (iresid = 0; iresid < nresid; iresid++)
    {
      /* Get this residue */
      residue_ptr = residue_index_ptr[iresid];

      /* If it has any variants, and has a mapping to the PDB file, then
         write out its details to both output files */
      if (residue_ptr->first_alist_ptr != NULL &&
          residue_ptr->r_residue_ptr != NULL)
        {
          /* Increment count of residues written out */
          nres++;

          /* Initialise variant count */
          have_added = FALSE;
          max_type = -1;
          nvar = 0;
          top_mut = 'X';
          
          /* Write out the UniProt residue details */
          fprintf(datfile_ptr,"SEQ_NO[%d] %d\n",nres,iresid + 1);
          fprintf(datfile_ptr,"AA_CODE[%d] %c\n",nres,residue_ptr->aa_code);

          /* Get the corresponding PDB residue */
          r_residue_ptr = residue_ptr->r_residue_ptr;
          
          /* Write out the PDB residue details */
          fprintf(datfile_ptr,"RES_NAME[%d] %s\n",nres,
                  r_residue_ptr->res_name);
          fprintf(datfile_ptr,"RES_NUM[%d] %s\n",nres,
                  r_residue_ptr->res_num);
          fprintf(datfile_ptr,"RES_CHAIN[%d] %c\n",nres,
                  r_residue_ptr->chain);

          /* Get the first variant */
          alist_ptr = residue_ptr->first_alist_ptr;

          /* Loop over all the variants */
          while (alist_ptr != NULL)
            {
              /* Get the corresponding annotation record */
              annot_ptr = alist_ptr->annot_ptr;
              
              /* Increment variant count */
              nvar++;

              /* Get the mutation type and corresponding mutation record */
              aamut_ptr = annot_ptr->aamut_ptr;

             /* Write out the variant */
              if (aamut_ptr != NULL)
                {
                  /* Get the mutation type */
                  type = aamut_ptr->var_type;

                  /* If this is an added variant, set flag */
                  if (type == ADDED_VAR)
                    {
                      /* Set flag */
                      have_added = TRUE;
                      added_var = TRUE;

                      /* Write out data */
                      fprintf(datfile_ptr,"NAME[%d][%d] %s\n",nres,nvar,
                          VARIANT_TYPE[type]);
                      fprintf(datfile_ptr,"IS_NOTE[%d][%d] TRUE\n",
                              nres,nvar);
                    }

                  /* If the first, then save */
                  if (max_type < 0)
                    {
                      /* Save this mutation's details */
                      max_type = type;
                      top_mut = aamut_ptr->aa_mut;
                    }

                  /* Write out the mutation */
                  fprintf(datfile_ptr,"MUT_TYPE[%d][%d] %d\n",nres,nvar,
                          type + 1);
                  fprintf(datfile_ptr,"MUT_TYPE_NAME[%d][%d] %s\n",nres,nvar,
                          VARIANT_TYPE[type]);
                  fprintf(datfile_ptr,"AA_MUT[%d][%d] %c\n",nres,nvar,
                          aamut_ptr->aa_mut);

                  /* Increment count of mutation types */
                  nmut_type[type]++;

                  /* Get the corresponding disease record */
                  disease_ptr = aamut_ptr->disease_ptr;

                  /* Write out the disease details */
                  if (disease_ptr != NULL)
                    {
                      fprintf(datfile_ptr,"DISEASE_ID[%d][%d] %s\n",
                              nres,nvar,disease_ptr->id);
                      fprintf(datfile_ptr,"NAME[%d][%d] %s\n",
                              nres,nvar,disease_ptr->name);
                      fprintf(datfile_ptr,"DESCRIPTION[%d][%d] %s\n",
                              nres,nvar,disease_ptr->description);
                      fprintf(datfile_ptr,"OMIM_ID[%d][%d] %s\n",
                              nres,nvar,disease_ptr->omim_id);
                      if (disease_ptr->is_note == TRUE)
                        fprintf(datfile_ptr,"IS_NOTE[%d][%d] TRUE\n",
                                nres,nvar);
                      else
                        fprintf(datfile_ptr,"IS_NOTE[%d][%d] FALSE\n",
                                nres,nvar);
                    }
                }
              
              /* Get the next variant */
              alist_ptr = alist_ptr->next_alist_ptr;
            }

          /* If residue has any associated variants, write it out to the
             .map file */
          if (nvar > 0)
            {
              /* Get this residue's last-defined atom */
              get_last_atom(r_residue_ptr,last_atom_name);

              /* Save the last atom */
              strcpy(residue_ptr->last_atom_name,last_atom_name);

              /* Write out line to the .map file */
              fprintf(mapfile_ptr,"%d %c -> %s %s %c :%c%c %s\n",
                      iresid + 1,residue_ptr->aa_code,
                      r_residue_ptr->res_name,
                      r_residue_ptr->res_num,
                      r_residue_ptr->chain,
                      top_mut,VARIANT_CHAR[max_type],last_atom_name);
              
              /* Write out this residue's colour */
              fprintf(datfile_ptr,"RES_COLOUR[%d] %d\n",nres,
                      max_type + 1);

              /* If we have an added variant at this residue position,
                 set the highlight flag */
              if (have_added == TRUE)
                fprintf(datfile_ptr,"HIGHLIGHT_ADDED[%d] TRUE\n",nres);
              else
                fprintf(datfile_ptr,"HIGHLIGHT_ADDED[%d] FALSE\n",nres);
            }

          /* Write out the number of this residue's variants */
          fprintf(datfile_ptr,"NVARIANTS[%d] %d\n",nres,nvar);
        }
    }
  
  /* Write out the number of residues to the .dat file */
  fprintf(datfile_ptr,"NRES %d\n",nres);

  /* If an added variant included, set flag */
  if (added_var == TRUE)
    {
      /* Set flag */
      fprintf(datfile_ptr,"HAVE_ADDED TRUE\n");

      /* Write out the variant itself */
      if (params.first_addmut_ptr != NULL)
        fprintf(datfile_ptr,"ADDED_VARIANT %s\n",
                params.first_addmut_ptr->added_variant);
    }
  else
    fprintf(datfile_ptr,"HAVE_ADDED FALSE\n");
  
  /* Initialise count of variant type written out */
  nvar = 0;

  /* Loop over the variant types, writing out number of each */
  for (i = 0; i < NVAR_TYPES; i++)
    {
      /* If not zero, then write out */
      if (nmut_type[i] > 0)
        {
          /* Increment count of variants written out */
          nvar++;

          /* Write out variant type and count */
          fprintf(datfile_ptr,"VARIANT_TYPE[%d] %d\n",nvar,i + 1);
          fprintf(datfile_ptr,"VARIANT_DESC[%d] %s\n",nvar,
                  variant_type[i]);
          fprintf(datfile_ptr,"VARIANT_COUNT[%d] %d\n",nvar,
                  nmut_type[i]);
        }
    }

  /* Write out the number of variant types written out */
  fprintf(datfile_ptr,"NVAR_TYPES %d\n",nvar);

  /* Close both ouput files */
  fclose(mapfile_ptr);
  fclose(datfile_ptr);
}
/***********************************************************************

create_domain_record  -  Create a new domain record

***********************************************************************/

struct domain *create_domain_record(struct domain **fst_domain_ptr,
                                    struct domain **lst_domain_ptr,
                                    char *pfam_id)
{
  int i;

  struct domain *domain_ptr;
  struct domain *first_domain_ptr, *last_domain_ptr;

  /* Initialise domain pointers */
  domain_ptr = NULL;
  first_domain_ptr = *fst_domain_ptr;
  last_domain_ptr = *lst_domain_ptr;

  /* Create domain entry */
  domain_ptr = (struct domain*)malloc(sizeof(struct domain));
  if (domain_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct domain\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  store_string(&(domain_ptr->pfam_id),pfam_id);

  /* Initialise other fields */
  domain_ptr->name = &null;
  domain_ptr->description = &null;
  domain_ptr->type = UNKNOWN;
  domain_ptr->colour = 0;
  domain_ptr->start = 0;
  domain_ptr->end = 0;
  domain_ptr->start_r_residue_ptr = NULL;;
  domain_ptr->end_r_residue_ptr = NULL;;

  /* Initialise pointer to next domain record */
  domain_ptr->next_domain_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of domained domain */
  if (first_domain_ptr == NULL)
    first_domain_ptr = domain_ptr;

  /* Otherwise, make previous entry point to the current one */
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

read_uniplot_dat  -  Read through the given PDB entry's uniplot.dat file
                     to pick up the Pfam domains

***********************************************************************/

int read_uniplot_dat(char *pdb_code,char chain,char *pdbsum_dir,
                     struct domain **fst_domain_ptr,char **seq,
                     char **pdb_seq)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];

  char *pdb_sequence, *sequence, *string_ptr;

  int len, ndomains;
  int done, wanted;

  struct domain *domain_ptr, *first_domain_ptr, *last_domain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  ndomains = 0;
  sequence = *seq = &null;
  pdb_sequence = *pdb_seq = &null;
  wanted = FALSE;

  /* Initialise pointers */
  *fst_domain_ptr = first_domain_ptr = last_domain_ptr = NULL;

  /* Form name of uniplot.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"uniplot.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file, line by line */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             done == FALSE)
        {
          /* Truncate the input line */
          len = string_chomp(input_line);

          /* Locate space */
          string_ptr = strchr(input_line,' ');

          /* If none found, the add two */
          if (string_ptr == NULL)
            {
              strcat(input_line,"  ");
              string_ptr = strchr(input_line,' ');
            }

          /* If this is the chain id, see if it is the one we want */
          if (!strncmp(input_line,"UNIPLOT_CHAIN[",14))
            {
              /* If it's the one we want, set flag */
              if (string_ptr[1] == chain || chain == '*')
                wanted = TRUE;

              /* Otherwise, see if we are done */
              else
                {
                  /* Unset flag */
                  wanted = FALSE;

                  /* If we already have the domains, can end here */
                  if (ndomains > 0)
                    done = TRUE;
                }
            }

          /* If this is a copy chain, see if it's the one we want */
          else if (!strncmp(input_line,"UNIPLOT_COPY[",13) &&
                   wanted == FALSE)
            {
              /* If it's the one we want, set flag */
              if (string_ptr[1] == chain || chain == '*')
                wanted = TRUE;
            }

          /* If this is our UniProt sequence, then save */
          else if (!strncmp(input_line,"UNIPLOT_UNIPROT_SEQ[",20))
            {
              /* Save the sequence */
              store_string(&sequence,string_ptr + 1);
            }

          /* If this is the aligned PDB sequence, then save */
          else if (!strncmp(input_line,"UNIPLOT_PDB_SEQ[",16))
            {
              /* Save the sequence */
              store_string(&pdb_sequence,string_ptr + 1);
            }

          /* If we are reading in the appropriate sequence's Pfam
             domains save domain record */
          else if (!strncmp(input_line,"UPLOT_PFAM_ID",13) &&
                   wanted == TRUE)
            {
              /* Create a Pfam domain record */
              domain_ptr
                = create_domain_record(&first_domain_ptr,
                                       &last_domain_ptr,
                                       string_ptr + 1);

              /* Increment count of domain records */
              ndomains++;
            }

          /* If have domain name, then save */
          else if (!strncmp(input_line,"UPLOT_DOMAIN_NAME",17) &&
                   wanted == TRUE && domain_ptr != NULL)
            store_string(&(domain_ptr->name),string_ptr+1);

          /* If have domain description, then save */
          else if (!strncmp(input_line,"UPLOT_DOMAIN_DESC",17) &&
                   wanted == TRUE && domain_ptr != NULL)
            store_string(&(domain_ptr->description),string_ptr+1);

          /* If Pfam A domain, then mark as such */
          else if (!strncmp(input_line,"UPLOT_PFAMA",11) &&
                   wanted == TRUE && domain_ptr != NULL)
            domain_ptr->type = PFAMA_DOMAIN;

          /* If Pfam B domain, then mark as such */
          else if (!strncmp(input_line,"UPLOT_PFAMB",11) &&
                   wanted == TRUE && domain_ptr != NULL)
            domain_ptr->type = PFAMB_DOMAIN;

          /* Save the domain start */
          else if (!strncmp(input_line,"UPLOT_PFAM_START",16) &&
                   wanted == TRUE && domain_ptr != NULL)
            domain_ptr->start = atoi(string_ptr + 1);

          /* Save the domain end */
          else if (!strncmp(input_line,"UPLOT_PFAM_END",14) &&
                   wanted == TRUE && domain_ptr != NULL)
            domain_ptr->end = atoi(string_ptr + 1);

          /* Get this domain's colour */
          else if (!strncmp(input_line,"UPLOT_PFAM_DOMAIN_COLOUR",23) &&
                   wanted == TRUE && domain_ptr != NULL)
            domain_ptr->colour = atoi(string_ptr + 1);
        }

      /* Close file */
      fclose(file_ptr);
    }

  /* Return current domain pointer */
  *fst_domain_ptr = first_domain_ptr;

  /* Return the UniProt and PDB sequences */
  *seq = sequence;
  *pdb_seq = pdb_sequence;

  /* Return number of domains */
  return(ndomains);
}
/***********************************************************************

map_domains_to_pdb  -  For covid-19 structure, colour according to
                         residue conservation (ie B-factor)

***********************************************************************/

void map_domains_to_pdb(struct domain *first_domain_ptr,char *sequence,
                        char *pdb_sequence,struct r_pdb *r_pdb_ptr,
                        char chain)
{
  char aa_start, aa_end, ch, pdb_aa_start, pdb_aa_end;

  int i, ires, len, pdb_len, pdb_start, pdb_end, seq_no, start, end;
  int done, started; 

  struct domain *domain_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  aa_start = aa_end = pdb_aa_start = pdb_aa_end = '?';
  done = FALSE;
  len = strlen(sequence);
  pdb_len = strlen(pdb_sequence);


  /* Get pointer to first domain */
  domain_ptr = first_domain_ptr;

  /* Loop over all the domains to colour */
  while (domain_ptr != NULL)
    {
      /* Initialise UniProt start and end positions */
      done = FALSE;
      seq_no = 0;
      start = end = -1;

      /* Loop over the UniProt sequence to find the start and end of
         this domain */
      for (i = 0; i < len && done == FALSE; i++)
        {
          /* If not a hyphen, then increment sequence position */
          ch = sequence[i];
          if (ch != '-')
            {
              /* Increment the sequence position */
              seq_no++;

              /* If this is the start of the domain, then save start
                 position */
              if (seq_no == domain_ptr->start)
                {
                  /* Save the start position and amino acid */
                  aa_start = ch;
                  start = i;
                }

              /* If we have started, save this as the domain end */
              if (start > -1)
                {
                  /* Save the end position and amino acid */
                  aa_end = ch;
                  end = i;
                }

              /* If this is the end of the domain, save it */
              if (seq_no == domain_ptr->end)
                {
                  /* Set flag that we're done */
                  done = TRUE;
                }
            }
        }

      /* Initialise PDB start and end positions */
      done = FALSE;
      seq_no = 0;
      pdb_start = pdb_end = -1;
      started = FALSE;

      /* If don't have valid start and end positions, then give up */
      if (start == -1 || end == -1)
        done = TRUE;

      /* Loop over the PDB sequence to find the start and end positions
         for this domain */
      for (i = 0; i < pdb_len && done == FALSE; i++)
        {
          /* See if this is the start-position of the aligned domain */
          if (i == start)
            started = TRUE;

          /* If not a hyphen, then increment sequence position */
          ch = pdb_sequence[i];
          if (ch != '-')
            {
              /* Increment the sequence position */
              seq_no++;

              /* If domain has started, then take this to be the first
                 PDB residue of the domain */
              if (started == TRUE)
                {
                  /* If the first, save this residue number */
                  if (pdb_start == -1)
                    {
                      pdb_start = seq_no;
                      pdb_aa_start = ch;
                    }

                  /* Save as the last PDB residue */
                  pdb_end = seq_no;
                  pdb_aa_end = ch;
                }
            }
             
          /* If this is the end-position unset flag */
          if (i == end)
            done = TRUE;
        }

      /* Initialise residue count */
      done = FALSE;
      ires = 0;
      started = FALSE;
      if (pdb_start == -1 || pdb_end == -1)
        done = TRUE;

      /* Get pointer to the first residue */
      r_residue_ptr = r_pdb_ptr->first_r_residue_ptr;

      /* Loop over all the residues to identify the start and end ones of
         the domain */
      while (r_residue_ptr != NULL && done == FALSE)
        {
          /* Increment count if this is the right chain */
          if (r_residue_ptr->chain == chain)
            ires++;

          /* If this is the start residue, then save */
          if (ires == pdb_start)
            {
              /* Save the starting residue */
              domain_ptr->start_r_residue_ptr = r_residue_ptr;

              /* Set flag */
              started = TRUE;
            }
              
          /* If we've started, save this as the end residue */
          if (started == TRUE)
            domain_ptr->end_r_residue_ptr = r_residue_ptr;

          /* If it is the end residue, we are done */
          if (ires == pdb_end)
            done = TRUE;

          /* Get pointer to the next residue */
          r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
        }

      /* Get the next domain record */
      domain_ptr = domain_ptr->next_domain_ptr;
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

int read_chain_equivs(struct chain **fst_chain_ptr,char *pdbsum_dir,
                      int pdb_type)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain_id, copy_of;

  int len, nchains, nlines;

  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;
  nlines = 0;

  /* Initialise chain pointers */
  chain_ptr = NULL;
  first_chain_ptr = NULL;
  last_chain_ptr = NULL;

  /* Form filename of the chains.dat file */
  strcpy(file_name,pdbsum_dir);
  if (pdb_type == PDBSUM1)
    strcpy(file_name,"../newsec/");
  strcat(file_name,"chains.dat");

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

          /* If have just a blank line, then take this to be chain
             blank */
          if (len == 0)
            input_line[0] = ' ';

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

          /* Increment chain count */
          nchains++;
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Return the number of chains stored */
  return(nchains);
}
/***********************************************************************

get_covalent_bonds  -  Read in any covalent bonds from the .res file and
                       add to the grow data

***********************************************************************/

int get_covalent_bonds(struct r_pdb *r_pdb_ptr,
                       struct s_grow **fst_s_grow_ptr,
                       struct s_grow **lst_s_grow_ptr,
                       struct r_rasdef *r_rasdef_ptr,char *pdbsum_dir)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char atom_name[2][5], chain[2], res_name[2][4], res_num[2][6];
  char lig_type;

  int i, j, match, ngrow;

  struct s_grow *s_grow_ptr, *first_s_grow_ptr, *last_s_grow_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  if (r_rasdef_ptr->type == R_METAL)
    lig_type = 'M';
  else
    lig_type = 'L';
  ngrow = 0;

  /* Initialise s_grow pointers */
  s_grow_ptr = NULL;
  first_s_grow_ptr = *fst_s_grow_ptr;
  last_s_grow_ptr = *lst_s_grow_ptr;

  /* Form the name of the .res file */
  strcpy(file_name,pdbsum_dir);
  if (r_rasdef_ptr->type == R_LIGAND)
    strcat(file_name,"ligplot");
  else
    strcat(file_name,"metplot");
  strcat(file_name,r_rasdef_ptr->ligplot_id);
  strcat(file_name,".res");

  /* Open the .res file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(ngrow);

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* If this is the target's accession code, check if it
         matches our protein */
      if (!strncmp(input_line,"#Covalent links: ",17))
        {
          string_truncate(input_line,LINELEN);

          /* Extract the atoms involved in the bond */
            strncpy(atom_name[0],input_line + 17,4);
            atom_name[0][4] = '\0';
            strncpy(res_name[0],input_line + 22,3);
            res_name[0][3] = '\0';
            strncpy(res_num[0],input_line + 26,5);
            res_num[0][5] = '\0';
            chain[0] = input_line[32];
            strncpy(atom_name[1],input_line + 37,4);
            atom_name[1][4] = '\0';
            strncpy(res_name[1],input_line + 42,3);
            res_name[1][3] = '\0';
            strncpy(res_num[1],input_line + 46,5);
            res_num[1][5] = '\0';
            chain[1] = input_line[52];

            /* See if either atom belongs to our rasdef residue */
            match = -1;
            if (!strcmp(res_num[0],r_rasdef_ptr->res_num) &&
                !strcmp(res_name[0],r_rasdef_ptr->res_name) &&
                chain[0] == r_rasdef_ptr->chain)
              match = 0;
            else if (!strcmp(res_num[1],r_rasdef_ptr->res_num) &&
                     !strcmp(res_name[1],r_rasdef_ptr->res_name) &&
                     chain[1] == r_rasdef_ptr->chain)
              match = 1;

            /* If found, then create a special grow record for this
               bond */
            if (match > -1)
              {
                /* Create a grow record */
                s_grow_ptr
                  = create_s_grow_record(&first_s_grow_ptr,
                                         &last_s_grow_ptr);

                /* Increment count of s_grow records stored */
                ngrow++;

                /* Save the imteracting atoms */
                for (i = 0; i < 2; i++)
                  {
                    /* Flip, if necessary to have the protein atom first */
                    if (match == 0)
                      j = 1 - i;
                    else
                      j = i;

                    /* Save the details of this atom */
                    strcpy(s_grow_ptr->atom_name[i],atom_name[j]);
                    strcpy(s_grow_ptr->res_name[i],res_name[j]);
                    strcpy(s_grow_ptr->res_num[i],res_num[j]);
                    s_grow_ptr->chain[i] = chain[j];
                  }

                /* Set the interaction type */
                s_grow_ptr->type = S_COVALENT;

                /* Save the residue types */
                s_grow_ptr->res_type[0] = 'P';
                s_grow_ptr->res_type[1] = lig_type;
              }
        }
    }
  
  /* Close the file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_s_grow_ptr = first_s_grow_ptr;
  *lst_s_grow_ptr = last_s_grow_ptr;

  /* Return number of grow records created */
  return(ngrow);
}
/***********************************************************************

add_covalent_bonds  -  Add any protein-ligand or protein-metal covalent
                       bonds, as found in the corresponding .res file

***********************************************************************/

void add_covalent_bonds(struct r_pdb *r_pdb_ptr,char *pdbsum_dir,
                        struct s_grow **fst_s_grow_ptr)
{
  int line;

  struct r_rasdef *r_rasdef_ptr;
  struct s_grow *s_grow_ptr, *first_s_grow_ptr, *last_s_grow_ptr;

  /* Initialise variables */
  line = 0;

  /* Inituialise pointers */
  last_s_grow_ptr = NULL;

  /* Get the first of the grow records */
  first_s_grow_ptr = *fst_s_grow_ptr;
  s_grow_ptr = first_s_grow_ptr;

  /* Loop through the grow records to attach those involving ligands
     to the corresponding ligand records */
  while (s_grow_ptr != NULL)
    {
      /* Save as the last record */
      last_s_grow_ptr = s_grow_ptr;

      /* Get pointer to next grow record */
      s_grow_ptr = s_grow_ptr->next_s_grow_ptr;
    }

  /* Get the first of the rasdef records */
  r_rasdef_ptr = r_pdb_ptr->first_r_rasdef_ptr;

  /* Loop through all the definition records */
  while (r_rasdef_ptr != NULL)
    {
      /* If this is a ligand or metal record, then check the .res file */
      if (r_rasdef_ptr->type == R_LIGAND || r_rasdef_ptr->type == R_METAL)
        {
          /* If we have a LIGPLOT id, can read the relevant file */
          if (strlen(r_rasdef_ptr->ligplot_id) > 4)
            {
              /* Read in any covalent bonds from the .res file and add
                 to the grow data */
              get_covalent_bonds(r_pdb_ptr,&first_s_grow_ptr,
                                 &last_s_grow_ptr,r_rasdef_ptr,
                                 pdbsum_dir);
            }
        }

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

  /* Return current pointer */
  *fst_s_grow_ptr = first_s_grow_ptr;
}
/***********************************************************************

check_residue  -  Check whether the residue from the grow record
                  corresponds to the metal residue or one of the
                  ligand's residues

***********************************************************************/

int check_residue(struct r_molec *r_molec_ptr,struct s_grow *s_grow_ptr,
                  int ires)
{
  int iresid, nresid;
  int found;

  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  found = FALSE;

  /* Get pointer to this molecule's first residue */
  r_residue_ptr = r_molec_ptr->first_r_residue_ptr;
  nresid = r_molec_ptr->nresid;
  iresid = 0;

  /* Loop through all the residues to find the one we're looking  for */
  while (r_residue_ptr != NULL && iresid < nresid && found == FALSE)
    {
      /* If this is the same residue as in the grow interaction, then
         set flag */
      if (r_residue_ptr->chain == s_grow_ptr->chain[ires] &&
          !strcmp(r_residue_ptr->res_num,s_grow_ptr->res_num[ires]) &&
          !strcmp(r_residue_ptr->res_name,s_grow_ptr->res_name[ires]))
        found = TRUE;

      /* Get the next r_residue */
      r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
      iresid++;
    }

  /* Return flag indicating if resdieu was found */
  return(found);
}
/***********************************************************************

filter_neighbours  -  Filter out all but the interacting molecules

***********************************************************************/

void filter_neighbours(struct r_pdb *r_pdb_ptr,char chain,
                       struct s_grow *first_s_grow_ptr)
{
  char inch, mol_type, res_name[4], res_num[6];

  int ires, nmolec;
  int check, wanted;

  struct r_molec *r_molec_ptr;
  struct r_residue *r_residue_ptr;
  struct s_grow *s_grow_ptr;

  /* Initialise variables */
  nmolec = 0;

  /* Get pointer to the first molecule in the PDB entry */
  r_molec_ptr = r_pdb_ptr->first_r_molec_ptr;

  /* Loop through all the molecules to calc coord limits */
  while (r_molec_ptr != NULL)
    {
      /* Increment molecule count */
      nmolec++;

      /* Initialise wanted flag */
      check = wanted = FALSE;
      mol_type = '\0';

      /* Get this molecule's first residue */
      r_residue_ptr = r_molec_ptr->first_r_residue_ptr;

      /* If this is a protein chain, then need to check interaction */
      if (r_molec_ptr->type == R_PROTEIN ||
          r_molec_ptr->type == R_CALPHA_ONLY)
        {
          /* If have a first residue, check chain id */
          if (r_residue_ptr != NULL)
            {
              /* If this is our chain, then don't need to perform check */
              if (r_residue_ptr->chain == chain)
                wanted = TRUE;
              else
                {
                  /* Set flags and molecule type */
                  check = TRUE;
                  mol_type = 'P';
                }
            }
        }

      /* If a DNA chain, ligand or metal, then need to check interaction */
      else if (r_molec_ptr->type == R_NUCLEIC_ACID ||
               r_molec_ptr->type == R_LIGAND ||
               r_molec_ptr->type == R_METAL)
        {
          /* Set flags and molecule type */
          check = TRUE;
          if (r_molec_ptr->type == R_NUCLEIC_ACID)
            mol_type = 'N';
          else if (r_molec_ptr->type == R_LIGAND)
            mol_type = 'L';
          else if (r_molec_ptr->type == R_METAL)
            mol_type = 'M';
        }

      /* If check required, then read through all the grow files to
         perform it */
      if (check == TRUE)
        {
          /* Get the first grow record */
          s_grow_ptr = first_s_grow_ptr;

          /* Loop through the grows to see if we have an interaction
             between the current molecule and the protein chain in
             question */
          while (s_grow_ptr != NULL && wanted == FALSE)
            {
              /* Initialise residue to be checked */
              ires = -1;

              /* If either interaction corresponds to our chain,
                 check the other interaction */
              if (s_grow_ptr->res_type[0] == 'P' &&
                  s_grow_ptr->chain[0] == chain &&
                  s_grow_ptr->res_type[1] == mol_type)
                ires = 1;
              else if (s_grow_ptr->res_type[1] == 'P' &&
                  s_grow_ptr->chain[1] == chain &&
                  s_grow_ptr->res_type[0] == mol_type)
                ires = 0;

              /* If we need to check the interaction, then do so */
              if (ires != -1 && s_grow_ptr->chain[ires]
                  == r_residue_ptr->chain)
                {
                  /* If we are checking protein or DNA, then we
                     have a match */
                  if (mol_type == 'P' || mol_type == 'N')
                    wanted = TRUE;

                  /* Otherwise, for ligand or metal, need to check
                     the residue name and number */
                  else
                    wanted = check_residue(r_molec_ptr,s_grow_ptr,ires);
                }

              /* Get the next grow record */
              s_grow_ptr = s_grow_ptr->next_s_grow_ptr;
            }
        }

      /* If the molecule is not wanted, change its type */
      if (wanted == FALSE)
        r_molec_ptr->type = R_UNKNOWN;
  
      /* Get the next molecule */
      r_molec_ptr = r_molec_ptr->next_r_molec_ptr;
    }
}
/***********************************************************************

identify_neighbours  -  Identify all neighbouring molecules that interact
                        with the given protein chain

***********************************************************************/

void identify_neighbours(struct r_pdb *r_pdb_ptr,char chain,
                         struct c_file_data file_name_ptr,
                         struct parameters params)
{
  char pdbsum_dir[FILENAME_LEN];

  int dir_type, ngrow, ndna, nlig, nprot;
  
  struct s_grow *first_s_grow_ptr;
  struct s_grow *dna_s_grow_ptr, *prot_s_grow_ptr;
  struct r_residue *first_r_residue_ptr;

  /* Find the PDBsum directory for this entry */
  dir_type
    = find_pdbsum_dir(r_pdb_ptr->pdb_code,file_name_ptr.pdbsum_template,
                      params.pdb_type,pdbsum_dir);

  /* Change directory for PDBsum1 */
  if (params.pdb_type == PDBSUM1)
    sprintf(pdbsum_dir,"../newsec/");

  /* Read in the data from the appropriate grow.out file */
  nlig = s_get_grow_data(r_pdb_ptr->pdb_code,pdbsum_dir,
                         &first_s_grow_ptr,S_GROW);

  /* Change directory for PDBsum1 */
  if (params.pdb_type == PDBSUM1)
    sprintf(pdbsum_dir,"../ligplot/");

  /* Add any protein-ligand or protein-metal covalent bonds,
     as found in the corresponding .res file */
  add_covalent_bonds(r_pdb_ptr,pdbsum_dir,&first_s_grow_ptr);

  /* Change directory for PDBsum1 */
  if (params.pdb_type == PDBSUM1)
    sprintf(pdbsum_dir,"../newsec/");

  /* Get the protein-protein interaction data */
  nprot = s_get_grow_data(r_pdb_ptr->pdb_code,pdbsum_dir,
                          &prot_s_grow_ptr,S_IGROW);

  /* Append all the interactions together */
  s_append_grow_records(&first_s_grow_ptr,prot_s_grow_ptr);

  /* Initialise the grow records prior to sorting */
  s_init_sort(first_s_grow_ptr);

  /* Sort the grow records by chain, residue number and interaction
     type */
  //s_sort_grow_records(&first_s_grow_ptr,ngrow);

  /* Filter out any molecules that do not interact with the given
     protein chain */
  filter_neighbours(r_pdb_ptr,chain,first_s_grow_ptr);

  /* Free up the memory taken by the grow data */
  s_free_grow_data(first_s_grow_ptr);
}
/***********************************************************************

get_chain_col  -  Get the given chain's colour

***********************************************************************/

int get_chain_col(char chain,double *rgb,double *rgb_hel,
                  double *rgb_str,int split,char *colour_name)
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

  /* Get the colour for this chain */
  colour = ichain % NCHAIN_COLOURS;
  strcpy(colour_name,CHAIN_COLOUR[colour]);

  /* Get rgb value for this colour */
  get_colour_rgb(colour_name,rgb);

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

write_pymol_headers  -  Write out the standard PyMOL script header lines

***********************************************************************/

void write_pymol_headers(FILE *file_ptr,char *pdb_file_name,
                         struct parameters params)
{
  char colour_name[COLNAM_LEN], rgb_string[20];

  int i, irgb[3], colour;
  int split;
  
  double rgb[3], rgb_hel[3], rgb_str[3];

  /* Initialise variables */
  split = TRUE;

  /* Write out the PDB code and initial parameters */
  fprintf(file_ptr,"#\n");
  fprintf(file_ptr,"# PDBsum PyMol script for PDB entry %s\n",
          params.pdb_code);
  fprintf(file_ptr,"#\n");
  fprintf(file_ptr,"load %s, %s, format=pdb\n",pdb_file_name,
          params.pdb_code);
  fprintf(file_ptr,"bg_color white\n");
  fprintf(file_ptr,"set depth_cue=0\n");
  fprintf(file_ptr,"set ray_trace_fog=0\n");
  fprintf(file_ptr,"set cartoon_ring_mode,1\n");
  fprintf(file_ptr,"hide all\n");
  fprintf(file_ptr,"set auto_color, 0\n");

  /* If not colouring by domain, define all the colours */
  if (params.colour_by_domain == FALSE)
    {
      /* Define all the RGB colours */
      for (colour = 0; colour < NCOLOURS; colour++)
        fprintf(file_ptr,"set_color %s, [ %6.3f, %6.3f, %6.3f ]\n",
                COLOUR_LIST[colour],RGB_TRIPLET[colour][0],
                RGB_TRIPLET[colour][1],RGB_TRIPLET[colour][2]);

      /* Define all the major and minor secondary structure colours */
      for (colour = 0; colour < NCHAIN_COLOURS; colour++)
        {
          /* Get this chain's colour */
          get_chain_col(CHAIN_LIST[colour],rgb,rgb_hel,rgb_str,split,
                        colour_name);

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

      /* If have an Alpha Fold structure, colour by model confidence */
      if (params.pdb_type == ALPHA_FOLD_PDB)
        {
          /* Define all the RGB colours */
          for (colour = 0; colour < NAF_RANGES; colour++)
            fprintf(file_ptr,"set_color %s, [ %s ]\n",
                    AF_COLOUR_CONF[colour],AF_COLOUR_RGB[colour]);
        }
    }

  /* If colouring by domain, define the domain colours */
  else
    {
      /* Define all the RGB colours */
      for (colour = 0; colour < NPFAM_COLOURS; colour++)
        {
          /* Get the RGB components of this colour */
          strcpy(rgb_string,COLOUR_DEFN[colour][2]);

          /* Extract the components */
          sscanf(rgb_string,"%d %d %d",&irgb[0],&irgb[1],&irgb[2]);

          /* Convert to floating point numbers */
          for (i = 0; i < 3; i++)
            rgb[i] = (double) irgb[i] / 255.0;

          /* Write out the colour definition */
          fprintf(file_ptr,"set_color %s, [ %6.3f, %6.3f, %6.3f ]\n",
                  COLOUR_DEFN[colour][0],rgb[0],rgb[1],rgb[2]);
        }
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

colour_alpha_fold  -  Use the AF_ranges file to assign the confidence
                      colours

***********************************************************************/

void colour_alpha_fold(FILE *pymol_script_ptr,char *pdb_code,
                       char *pdbsum_dir)
{
  int i, ires, first_res, last_range, last_res, len, nscores, range;

  struct c_score *c_score_ptr, *first_c_score_ptr;

  /* Initialise variables */
  last_range = -1;
  first_res = last_res = -1;

  /* Get the confidence scores */
  nscores = get_c_scores(pdb_code,pdbsum_dir,&first_c_score_ptr);

  /* Get first score record */
  c_score_ptr = first_c_score_ptr;

  /* Loop over  the scores */
  while (c_score_ptr != NULL)
    {
      /* Get this residue */
      ires = atoi(c_score_ptr->res_num);

      /* See if this is a new range */
      if (c_score_ptr->range != last_range)
        {
          /* If have a previous range, then write out */
          if (last_range != -1)
            fprintf(pymol_script_ptr,"color %s, (resi %d-%d)\n",
                    AF_COLOUR_CONF[last_range],first_res,last_res);

          /* Save this residue as the first of the new range */
          first_res = ires;
        }

      /* Save this residue as the last of the current range */
      last_res = ires;

      /* Save the range */
      last_range = c_score_ptr->range;

      /* Get the next record */
      c_score_ptr = c_score_ptr->next_c_score_ptr;
    }

  /* Write out the final range, if there is one */
  if (last_range != -1)
    fprintf(pymol_script_ptr,"color %s, (resi %d-%d)\n",
            AF_COLOUR_CONF[last_range],first_res,last_res);
}
/***********************************************************************

write_pymol_chains  -  Create PyMol objects for each protein and DNA
                       chain in the structure

***********************************************************************/

void write_pymol_chains(FILE *pymol_script_ptr,struct r_pdb *r_pdb_ptr,
                        struct pobject **fst_pobject_ptr,
                        struct pobject **lst_pobject_ptr,
                        int colour_as_orig,int transparent_copies,
                        int biomol,struct c_file_data file_name_ptr,
                        struct parameters params)
{
  char chain, object[TOKEN_LEN + 1], chain_name[TOKEN_LEN + 1];
  char colour_name[COLNAM_LEN], display_type[20];
  char colour_major[COLNAM_LEN + 6], colour_minor[COLNAM_LEN + 6];
  char colour_coil[COLNAM_LEN + 5], number_string[20];
  char name[3], res_num1[6], res_num2[6], residue_range[13];
  char pdb_code[5], pdbsum_dir[FILENAME_LEN];

  int colour, iresid, metal_type, nresid, nsite, number, type, valid;
  int dir_type, transparency;
  int show, split, wanted;

  double rgb[3], rgb_hel[3], rgb_str[3];

  struct r_chain *r_chain_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  strcpy(colour_name,"white");
  strcpy(colour_coil,"white");
  strcpy(colour_major,"white");
  strcpy(colour_minor,"white");
  show = TRUE;
  split = TRUE;
  transparency = PARTIAL;
  strcpy(pdb_code,r_pdb_ptr->pdb_code);

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Get pointer to the first chain */
  r_chain_ptr = r_pdb_ptr->first_r_chain_ptr;

  /* Loop through all the chains to define the various objects */
  while (r_chain_ptr != NULL)
    {
      /* Get chain type */
      type = r_chain_ptr->type;

      /* Determine if this chain is wanted */
      wanted = FALSE;
      if (type == R_PROTEIN || type == R_NUCLEIC_ACID ||
          type == R_CALPHA_ONLY)
        wanted = TRUE;

      /* If colouring by domain, show all chains except ours */
      if (params.colour_by_domain == TRUE &&
          r_chain_ptr->chain_id == params.chain)
        wanted = FALSE;

      /* For PDBsum1, only show if chain selected */
      if (params.pdb_type == PDBSUM1 && wanted == TRUE)
        {
          /* Only include if chain matches */
          if (params.chain == '*' ||
              r_chain_ptr->chain_id == params.chain)
            wanted = TRUE;
          else
            wanted = FALSE;
        }

      /* If protein or DNA, create a name for this chain */
      if (wanted == TRUE)
        {
          /* Show comment */
          fprintf(pymol_script_ptr,"# Definition for chain [%c] - ",
                  r_chain_ptr->chain_id);
          if (type == R_PROTEIN)
            fprintf(pymol_script_ptr,"PROTEIN\n");
          else if (type == R_NUCLEIC_ACID)
            fprintf(pymol_script_ptr,"DNA\n");
          else
            fprintf(pymol_script_ptr,"CALPHA_ONLY\n");

          /* Determine whether this chain is to be shown */
          show = TRUE;
          if (params.res_cons == TRUE &&
              r_chain_ptr->chain_id == params.chain)
            show = FALSE;

          /* Initialise the first residue and last residues */
          res_num1[0] = '\0';
          res_num2[0] = '\0';

          /* Get this chain's first residue */
          r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
          nresid = r_chain_ptr->nresid;
          iresid = 0;
          valid = FALSE;

          /* Loop through this chain's residues, to find start and
             end residues */
          while (r_residue_ptr != NULL && iresid < nresid)
            {
              /* If this is the first residue, save number */
              if (res_num1[0] == '\0')
                strcpy(res_num1,r_residue_ptr->res_num);

              /* Save number as last residue */
              strcpy(res_num2,r_residue_ptr->res_num);

              /* Get the next residue */
              r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
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
          if (r_chain_ptr->chain_id == ' ')
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
              sprintf(object,"chain_%c",r_chain_ptr->chain_id);

              /* If this is a duplicate biomol chain, include
                 residue definition in its object name */
              if (biomol == TRUE && r_chain_ptr->symm_chain == TRUE &&
                  valid == TRUE)
                fprintf(pymol_script_ptr,
                        "create %s, (resi %s and chain %c and not het) \n",
                        object,residue_range,r_chain_ptr->chain_id);

              /* Otherwise write out the definition using chain name
                 only */
              else
                fprintf(pymol_script_ptr,
                        "create %s, (chain %c and not het) \n",
                        object,r_chain_ptr->chain_id);
            }

          /* Create a pymol object record */
          pobject_ptr
            = create_pobject_record(&first_pobject_ptr,&last_pobject_ptr,
                                    object,type,r_chain_ptr->chain_id);

          /* Get the chain */
          chain = r_chain_ptr->chain_id;
          if (biomol == TRUE && colour_as_orig == TRUE)
            chain = r_chain_ptr->link_r_chain_ptr->chain_id;

          /* Define the chain colour */
          if (params.colour_by_domain == FALSE)
            get_chain_colname(chain,colour_name,colour_coil,colour_major,
                              colour_minor);

          /* If this is an Alpha Fold model, use the AF_ranges file to
             assign the confidence colours */
          if (params.pdb_type == ALPHA_FOLD_PDB)
            {
              /* Get the PDBsum directory of this protein */
              dir_type
                = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                                  params.pdb_type,pdbsum_dir);

              /* Colour the structure */
              colour_alpha_fold(pymol_script_ptr,pdb_code,pdbsum_dir);

              /* Unset the show flag */
              show = FALSE;
            }

          /* If colouring, do so according to type */
          if (show == TRUE)
            {
              if (type == R_NUCLEIC_ACID)
                fprintf(pymol_script_ptr,"color %s, %s\n",colour_name,
                        object);
              else
                fprintf(pymol_script_ptr,"color %s, %s\n",colour_coil,
                        object);

              /* Define secondary structure colours according to which
                 secondary structure type is the "dominant" one */
              if (type != R_NUCLEIC_ACID &&
                  r_chain_ptr->dominant_sst == R_HELIX)
                {
                  fprintf(pymol_script_ptr,"color %s, (%s and ss h)\n",
                          colour_major,object);
                  fprintf(pymol_script_ptr,"color %s, (%s and ss s)\n",
                          colour_minor,object);
                }
              else if (type != R_NUCLEIC_ACID)
                {
                  fprintf(pymol_script_ptr,"color %s, (%s and ss h)\n",
                          colour_minor,object);
                  fprintf(pymol_script_ptr,"color %s, (%s and ss s)\n",
                          colour_major,object);
                }
            }

          /* Define display type for protein */
          if (type == R_PROTEIN)
            strcpy(display_type,"cartoon");

          /* Cartoon for CA */
          else if (type == R_CALPHA_ONLY)
            {
              fprintf(pymol_script_ptr,"set cartoon_trace, 1, %s\n",
                      object);
              strcpy(display_type,"cartoon");
            }

          /* Sticks for DNA/RNA */
          else if (type == R_NUCLEIC_ACID)
            strcpy(display_type,"cartoon");

          /* Determine whether chain to be shown transparent */
          if (biomol == TRUE && r_chain_ptr->symm_chain == TRUE &&
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

      /* If this is a protein chain, check for any residues belonging
         to active sites */
      if (type == R_PROTEIN)
        {
          /* Save the name of the chain */
          strcpy(chain_name,object);

          /* Define name of site */
          if (r_chain_ptr->chain_id == ' ')
            strcpy(object,"sites_");
          else
            sprintf(object,"sites_%c",r_chain_ptr->chain_id);

          /* Initialiase residue counter */
          nsite = 0;

          /* Get pointer to the first of this chain's residues */
          r_residue_ptr = r_chain_ptr->first_r_residue_ptr;
          nresid = r_chain_ptr->nresid;
          iresid = 0;

          /* Loop through this chain's residues to find any active site
             residues */
          while (r_residue_ptr != NULL && iresid < nresid)
            {
              /* If this is an active site residue, add to definition */
              if (r_residue_ptr->site[0] != '\0')
                {
                  /* Get the residue number and check that valid */
                  strcpy(number_string,r_residue_ptr->res_num);
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
                                  r_chain_ptr->chain_id);
                          fprintf(pymol_script_ptr,
                                  "create %s, "
                                  "(%s and resi %s and not het)\n",
                                  object,chain_name,r_residue_ptr->res_num);
                        }

                      /* Otherwise, add residue to current definition */
                      else
                        {
                          fprintf(pymol_script_ptr,
                                  "create %s, (%s or "
                                  "(%s and resi %s and not het))\n",
                                  object,object,chain_name,
                                  r_residue_ptr->res_num);
                        }

                      /* Increment count of site residues in this chain */
                      nsite++;
                    }
                }

              /* Get the next residue */
              r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
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
                                        r_chain_ptr->chain_id);
              store_string(&(pobject_ptr->display_type),"sticks");


              /* Flag chain as having active sites */
              pobject_ptr->got_sites = TRUE;
            }
        }

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }

  /* Return current pobject pointers */
  *fst_pobject_ptr = first_pobject_ptr;
  *lst_pobject_ptr = last_pobject_ptr;
}
/***********************************************************************

write_pymol_ligands  -  Create PyMol objects for each ligand in the
                        structure

***********************************************************************/

void write_pymol_ligands(FILE *pymol_script_ptr,struct r_pdb *r_pdb_ptr,
                         struct pobject **fst_pobject_ptr,
                         struct pobject **lst_pobject_ptr,
                         int transparent_copies,int biomol)
{
  char chain, range_name[20], residue_range[20];

  int ilig, nligands;

  struct r_ligand *r_ligand_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct r_residue *r_residue_ptr;
  struct s_grow *s_grow_ptr;

  /* Initialise variables */
  nligands = 0;

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Get pointer to the first ligand in the PDB entry */
  r_ligand_ptr = r_pdb_ptr->first_r_ligand_ptr;

  /* Loop through all the ligands to calc coord limits */
  while (r_ligand_ptr != NULL)
    {
      /* If ligand, then show as spheres */
      if (r_ligand_ptr->metal == FALSE)
        {
          /* Increment ligand count */
          nligands++;

          /* If first ligand, write out comment */
          if (nligands == 1)
            fprintf(pymol_script_ptr,"# Ligands\n");

          /* Form the ligand name */
          sprintf(range_name,"lrange_%d",nligands);

          /* Get the chain and start- and end-residues for ligand */
          chain = r_ligand_ptr->first_residue[8];
          strncpy(residue_range,r_ligand_ptr->first_residue + 3,5);
          residue_range[5] = '\0';
          strcat(residue_range,"-");
          strncat(residue_range,r_ligand_ptr->last_residue + 3,5);
          residue_range[11] = '\0';

          /* Remove leading spaces from ligand number */
          remove_residue_blanks(residue_range);

          /* Write out its definition */
          fprintf(pymol_script_ptr,
                  "create %s, (resi %s and chain %c and het)\n",
                  range_name,residue_range,chain);
        }
      
      /* Get the next ligand */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }

  /* If have any ligands, write out the commands to render them */
  if (nligands > 0)
    {
      /* Form the ligands definition */
      fprintf(pymol_script_ptr,"create ligands, (");

      /* Loop over the ligand ranges to write out */
      for (ilig = 0; ilig < nligands; ilig++)
        {
          /* If not the first, write out " or " */
          if (ilig > 0)
            fprintf(pymol_script_ptr," or ");

          /* Write out this ligand range name */
          fprintf(pymol_script_ptr,"lrange_%d",ilig + 1);
        }

      /* Close off the list of ranges */
      fprintf(pymol_script_ptr,")\n");

      /* Write out the rendering commands */
      fprintf(pymol_script_ptr,"show spheres, ligands\n");
      fprintf(pymol_script_ptr,"util.cbag ligands\n");
      fprintf(pymol_script_ptr,"color black, (ligands and elem C)\n");
    }
}
/***********************************************************************

write_pymol_metals  -  Create PyMol objects for each metal ion in the
                        structure

***********************************************************************/

void write_pymol_metals(FILE *pymol_script_ptr,struct r_pdb *r_pdb_ptr,
                        struct pobject **fst_pobject_ptr,
                        struct pobject **lst_pobject_ptr,
                        int transparent_copies,int biomol)
{
  char chain, residue_range[20];
  char range_name[20];

  int nmetals;

  struct r_ligand *r_ligand_ptr;
  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
  struct r_residue *r_residue_ptr;
  struct s_grow *s_grow_ptr;

  /* Initialise variables */
  nmetals = 0;

  /* Initialise pointers */
  first_pobject_ptr = *fst_pobject_ptr;
  last_pobject_ptr = *lst_pobject_ptr;

  /* Get pointer to the first ligand in the PDB entry */
  r_ligand_ptr = r_pdb_ptr->first_r_ligand_ptr;

  /* Loop through all the ligands to calc coord limits */
  while (r_ligand_ptr != NULL)
    {
      /* If ligand, then show as spheres */
      if (r_ligand_ptr->metal == TRUE)
        {
          /* Increment ligand count */
          nmetals++;

          /* If first metal, write out comment */
          if (nmetals == 1)
            fprintf(pymol_script_ptr,"# Metals\n");

          /* Create the metal definition */
          fprintf(pymol_script_ptr,"create metal%s, (elem %s)\n",
                  r_ligand_ptr->metal_name,r_ligand_ptr->metal_name);

          /* Define its colour */
          fprintf(pymol_script_ptr,"colour %s, metal%s\n",
                  R_METAL_COLOUR[r_ligand_ptr->metal_colour],
                  r_ligand_ptr->metal_name);

          /* Define its size */
          fprintf(pymol_script_ptr,"alter metal%s, vdw=1.5\n",
                  r_ligand_ptr->metal_name);

          /* Render the metal */
          fprintf(pymol_script_ptr,"rebuild\n");
          fprintf(pymol_script_ptr,"show spheres, metal%s\n",
                  r_ligand_ptr->metal_name);
        }
      
      /* Get the next ligand */
      r_ligand_ptr = r_ligand_ptr->next_r_ligand_ptr;
    }
}
/***********************************************************************

get_max_variants  -  Get each residue's 'highest' type of variant

***********************************************************************/

int get_max_variants(struct residue **residue_index_ptr,int nresid,
                     int **vtype)
{
  int iresid, max_type, nres, nvar, type;

  int *var_type;

  struct aamut *aamut_ptr;
  struct annot *annot_ptr;
  struct alist *alist_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  nres = nvar = 0;

  /* Create an array of variant types */
  var_type = (int *) malloc (nresid * sizeof(int));
  if (var_type == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for var_type\n");
      exit(1);
    }

  /* Initialise the array */
  for (iresid = 0; iresid < nresid; iresid++)
    var_type[iresid] = -1;

  /* Loop over all the residues in the sequence */
  for (iresid = 0; iresid < nresid; iresid++)
    {
      /* Get this residue */
      residue_ptr = residue_index_ptr[iresid];

      /* If it has any variants, then save 'highest' type */
      if (residue_ptr->first_alist_ptr != NULL &&
          residue_ptr->r_residue_ptr != NULL)
        {
          /* Initialise variant count */
          max_type = -1;
          nvar = 0;

          /* Get the first variant */
          alist_ptr = residue_ptr->first_alist_ptr;

          /* Loop over all the variants */
          while (alist_ptr != NULL)
            {
              /* Get the corresponding annotation record */
              annot_ptr = alist_ptr->annot_ptr;
              
              /* Increment variant count */
              nvar++;

             /* Get the mutation type and corresponding mutation record */
              aamut_ptr = annot_ptr->aamut_ptr;

             /* Get the variant type */
              if (aamut_ptr != NULL)
                {
                  /* Get the mutation type */
                  type = aamut_ptr->var_type;
                  
                  /* If the first, then save */
                  if (max_type < 0)
                    {
                      /* Save this mutation's details */
                      max_type = type;
                    }
                }
              
              /* Get the next variant */
              alist_ptr = alist_ptr->next_alist_ptr;
            }

          /* If residue has any associated variants, store the highest */
          if (nvar > 0)
            {
              var_type[iresid] = max_type;

              /* Increment count of residues with variants */
              nres++;
            }
        }
    }

  /* Return the array of maximum variant types */
  *vtype = var_type;

  /* Return count of residues with variants */
  return(nres);
}
/***********************************************************************

write_pymol_variants  -  Write the variants to the PyMOL script

***********************************************************************/

void write_pymol_variants(FILE *pymol_script_ptr,
                          struct residue **residue_index_ptr,int nresid,
                          int labels)
{
  char defn_string[LONG_LINELEN + 1], residue_def[30];
  char label[30], res_name[4];

  int iresid, nlabel, nres, nvar_res, type;

  int *var_type;

  struct residue *residue_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  nlabel = 0;

  /* Determine each residue's 'highest' variant type */
  nvar_res = get_max_variants(residue_index_ptr,nresid,&var_type);

  /* Loop over the variant types, rendering the list of residues
     belonging to each type */
  for (type = 0; type < NVAR_TYPES && nvar_res > 0; type++)
    {
      /* Initialise count of residues belonging to this type */
      nres = 0;

      /* Initialise string holding list of this type's residues */
      defn_string[0] = '\0';

      /* Loop over all the residues to pick out the ones of this type */
      for (iresid = 0; iresid < nresid; iresid++)
        {
          /* If residue of the correct type, then add to definition */
          if (var_type[iresid] == type)
            {
              /* Get this residue */
              residue_ptr = residue_index_ptr[iresid];

              /* Get the corresponding PDB residue */
              r_residue_ptr = residue_ptr->r_residue_ptr;
          
              /* Convert residue number into PyMOL format */
              if (r_residue_ptr->chain != ' ')
                sprintf(residue_def,"(resi %s & chain %c)",
                        r_residue_ptr->res_num,r_residue_ptr->chain);
              else
                sprintf(residue_def,"resi %s",r_residue_ptr->res_num);

              /* If first one, start definition */
              if (nres == 0)
                sprintf(defn_string,"create var_type%d, (",type + 1);
              else
                strcat(defn_string," | ");

              /* Select the residue */
              strcat(defn_string,residue_def);

              /* Increment count of residues of this type */
              nres++;

              /* If labelling variant, write out the label */
              if (labels == TRUE)
                {
                  /* Increment label count */
                  nlabel++;
                  
                  /* Convert residue name to lower case */
                  strcpy(res_name,r_residue_ptr->res_name);
                  to_lower(res_name + 1,2);

                  /* Form name of residue */
                  sprintf(label,"%s%s",res_name,r_residue_ptr->res_num);

                  /* Remove any blanks from the label */
                  remove_blanks(label);
      
                  /* Select the atom */
                  fprintf(pymol_script_ptr,
                          "create label%d, name %s and resi %s and chain %c\n",
                          nlabel,residue_ptr->last_atom_name,
                          r_residue_ptr->res_num,
                          r_residue_ptr->chain);

                  /* Show label */
                  fprintf(pymol_script_ptr,"label label%d, \"%s\"\n",
                          nlabel,label);

                  /* Define label colour */
                  fprintf(pymol_script_ptr,
                          "set label_color, %s, label%d\n",
                          PYMOL_COLOUR[type],nlabel);
                }
            }
        }

      /* If we have any residues of this type, render them */
      if (nres > 0)
        {
          /* Close off the definition string */
          strcat(defn_string,")");

          /* Write out the definition */
          fprintf(pymol_script_ptr,"%s\n",defn_string);
      
          /* Define this variant's colour */
          fprintf(pymol_script_ptr,"color %s, var_type%d\n",
                  PYMOL_COLOUR[type],type + 1);

          /* Show as sticks */
          fprintf(pymol_script_ptr,"set stick_radius, 0.2, var_type%d\n",
                  type + 1);
          fprintf(pymol_script_ptr,"show sticks, var_type%d\n",
                  type + 1);
          fprintf(pymol_script_ptr,"rebuild\n");
        }
    }

  /* Free the var_type array */
  free(var_type);
}
/***********************************************************************

write_pymol_surface  -  Show selected chain as a surface and colour by
                        residue conservation

***********************************************************************/

void write_pymol_surface(FILE *pymol_script_ptr,struct parameters params)
{
  char chain, chain_name[TOKEN_LEN + 1];

  /* Initialise variables */
  chain = params.chain;

  /* Show heading */
  fprintf(pymol_script_ptr,"\n");
  fprintf(pymol_script_ptr,"# Residue conservation\n");

  /* Show chain as surface */
  fprintf(pymol_script_ptr,"show surface, chain_%c\n",chain);

  /* Colour by B-factor */
  fprintf(pymol_script_ptr,"spectrum b, rainbow, chain_%c, "
          "minimum=0, maximum=80\n",chain);

  /* Define transparency */
  fprintf(pymol_script_ptr,"set transparency, 0.2\n");
}
/***********************************************************************

write_pymol_domcol  -  Colour according to Pfam domain

***********************************************************************/

void write_pymol_domcol(FILE *pymol_script_ptr,
                        struct domain *first_domain_ptr,
                        struct parameters params)
{
  char chain,  colour_name[20], colour_string[20], res_range[20];
  char object[OBJ_NAME_LEN];

  int icol, ndomain, start_col;

  struct domain *domain_ptr;

  /* Initialise variables */
  chain = params.chain;
  ndomain = 0;

  /* Write out heading */
  fprintf(pymol_script_ptr,"\n");
  fprintf(pymol_script_ptr,"# Colour chain %c by domain\n",chain);

  /* Show the given chain as a white cartoon */
  fprintf(pymol_script_ptr,"show cartoon, (chain %c and not het)\n",
          chain);
  fprintf(pymol_script_ptr,"color Black, (chain %c and not het)\n",
          chain);
  fprintf(pymol_script_ptr,"\n");

  /* Get the starting domain colour */
  start_col = PFAMA_STARTCOL;

  /* If this is a covid-19 protein, define corresponding start colour */
  if (params.start_colour > 0)
    start_col = params.start_colour;

  /* Get pointer to first domain */
  domain_ptr = first_domain_ptr;

  /* Loop over all the domains to colour */
  while (domain_ptr != NULL)
    {
      /* If have sytart and end residues, define the domain and its
         colour */
      if (domain_ptr->start_r_residue_ptr != NULL &&
          domain_ptr->end_r_residue_ptr != NULL)
        {
          /* Increment domain number */
          ndomain++;

          /* Create the domain name */
          sprintf(object,"domain%d",ndomain);

          /* Determine this domain's colour */
          icol = (start_col + domain_ptr->colour) % (NPFAM_COLOURS - 2);
          strcpy(colour_name,COLOUR_DEFN[icol][0]);

          /* Write out the domain name and colour */
          fprintf(pymol_script_ptr,"# Domain %d - %s = %s (%s)\n",
                  ndomain,domain_ptr->pfam_id,domain_ptr->name,
                  colour_name);

          /* Form the names of the start and end residues */
          sprintf(res_range,"%s-%s",
                  domain_ptr->start_r_residue_ptr->res_num,
                  domain_ptr->end_r_residue_ptr->res_num);

          /* Remove any blanks from the definition */
          remove_blanks(res_range);

          /* Create the domain */
          fprintf(pymol_script_ptr,
                  "create %s, (resi %s & chain %c and not het)\n",
                  object,res_range,chain);

          /* Set the colour */
          fprintf(pymol_script_ptr,
                  "color %s, (resi %s & chain %c and not het)\n",
                  colour_name,res_range,chain);
          fprintf(pymol_script_ptr,"\n");
        }

      /* Get the next domain record */
      domain_ptr = domain_ptr->next_domain_ptr;
    }
}
/***********************************************************************

write_pymol_pml  -  Write command to PyMol script to define all objects
                    and their representation

***********************************************************************/

void write_pymol_pml(struct r_pdb *r_pdb_ptr,FILE **pymol_script_ptr,
                     struct pobject **fst_pobject_ptr,
                     struct residue **residue_index_ptr,int nresid,
                     struct domain *first_domain_ptr,char *pdir,
                     struct c_file_data file_name_ptr,
                     struct parameters params)
{
  char file_name[FILENAME_LEN + 1], number_string[20];
  char pdb_file_name[FILENAME_LEN + 1];
  char colour_name[COLNAM_LEN], display_type[20];
  char colour_coil[COLNAM_LEN + 5];

  char *pdbsum_dir;

  int iresid, metal_type, nsite, number, type, valid;
  int colour_as_orig, transparent_copies, biomol, wanted;

  struct pobject *pobject_ptr, *first_pobject_ptr, *last_pobject_ptr;
 
  FILE *file_ptr;

  /* Initialise variables */
  colour_as_orig = TRUE;
  transparent_copies = biomol = FALSE;

  /* Initialise pointers */
  pobject_ptr = first_pobject_ptr = last_pobject_ptr = NULL;

  /* Form name of pymol script file, pymol.pml */
  if (params.pml_file != &null)
    strcpy(file_name,params.pml_file);
  else
    strcpy(file_name,"genpymol.pml");
  strcpy(pdb_file_name,"genpymol.pdb");

  /* For PDBsum1, put full file name */
  if (params.pdb_type == PDBSUM1)
    {
      /* Get the name of the PDBsum directory */
      pdbsum_dir
        = form_filename(params.pdb_code,
                        file_name_ptr.other_dir[PDBSUM_TEMPLATE]);

      /* Form the name of the PDB file */
      sprintf(pdb_file_name,"%s/%s.new",pdbsum_dir,params.pdb_code);
    }

  /* If running on the web, pick up the coordinates via the appropriate
     URL */
  if (params.user_id > 0)
    {
      /* Form the name of the URL for retrieving the PDB coordinates */
      sprintf(pdb_file_name,"%s%d",
              file_name_ptr.other_dir[GET_PYMOL_FILE],
              params.user_id);
    }

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);
  if (file_ptr == NULL)
    {
      printf("*** ERROR. Unable to open output file %s\n",file_name);
      exit(1);
    }

  /* Write out the header records */
  write_pymol_headers(file_ptr,pdb_file_name,params);

  /* Write out the definitions for all the protein and DNA chains */
  write_pymol_chains(file_ptr,r_pdb_ptr,&first_pobject_ptr,
                     &last_pobject_ptr,colour_as_orig,transparent_copies,
                     biomol,file_name_ptr,params);

  /* If showing surface coloured by residue conservation, write out
     commands */
  if (params.res_cons == TRUE)
    write_pymol_surface(file_ptr,params);

  /* If colouring by domain, use the domain definitions */
  if (params.colour_by_domain == TRUE)
    write_pymol_domcol(file_ptr,first_domain_ptr,params);

  /* For PDBsum1, if single chain selected, don't show ligands or
     metals */
  wanted = TRUE;
  if (params.pdb_type == PDBSUM1 && params.chain != '*')
    wanted = FALSE;

  /* If wanted, write all the ligands and metals */
  if (wanted == TRUE)
    {
      /* Write out the definitions for all the ligands */
      write_pymol_ligands(file_ptr,r_pdb_ptr,&first_pobject_ptr,
                          &last_pobject_ptr,transparent_copies,biomol);

      /* Write out the definitions for all the metals */
      write_pymol_metals(file_ptr,r_pdb_ptr,&first_pobject_ptr,
                         &last_pobject_ptr,transparent_copies,biomol);
    }

  /* Write out any variants */
  if (params.variants == TRUE)
    write_pymol_variants(file_ptr,residue_index_ptr,nresid,
                         params.labels);

  /* Rotate by 180 about the x-axis */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"# Rotate 180 about the x-axis\n");
  fprintf(file_ptr,"turn x, 180\n");
  

  /* Zoom in on final picture */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"# Zoom in on the structure\n");
  if (params.chain == '*')
    fprintf(file_ptr,"zoom all, complete=1, buffer=2\n");
  else
    fprintf(file_ptr,"zoom chain_%c, complete=1, buffer=2\n",
            params.chain);

  /* Return first pobject pointer */
  *fst_pobject_ptr = first_pobject_ptr;

  /* Return pointer to script file */
  *pymol_script_ptr = file_ptr;

  /* Close the PyMol script file */
  fclose(file_ptr);
}
/***********************************************************************

write_coords  -  Write out the required coordinates

***********************************************************************/

void write_coords(struct r_pdb *r_pdb_ptr,FILE *file_ptr,int use_atomno)
{
  char file_name[FILENAME_LEN + 1], number_string[20];

  int atom_number, nmolec, type;
  int new_file;

  struct r_atom *r_atom_ptr;
  struct r_molec *r_molec_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  atom_number = nmolec = 0;
  new_file = FALSE;

  /* If we don't already have an open file, create */
  if (file_ptr == NULL)
    {
      /* Form name of the output file */
      strcpy(file_name,"genpymol.pdb");

      /* Open the output file */
      file_ptr = open_output_file(file_name,TRUE);
      if (file_ptr == NULL)
        {
          printf("*** ERROR. Unable to open output file %s\n",file_name);
          exit(1);
        }

      /* Set flag that new */
      new_file = TRUE;
    }

  /* Get pointer to the first molecule in the PDB entry */
  r_molec_ptr = r_pdb_ptr->first_r_molec_ptr;

  /* Loop through all the molecules to calc coord limits */
  while (r_molec_ptr != NULL)
    {
      /* Get the molecule type */
      type = r_molec_ptr->type;

      /* If this molecule is one of the ones we want, then write out */
      if (type == R_PROTEIN || type == R_CALPHA_ONLY ||
          type == R_NUCLEIC_ACID || type == R_LIGAND ||
          type == R_METAL)
        {
          /* Increment molecule count */
          nmolec++;

          /* If using the stored atom numbers, get the first */
          if (use_atomno == TRUE)
            {
              /* Get the first residue */
              r_residue_ptr = r_molec_ptr->first_r_residue_ptr;

              /* Get its first atom */
              if (r_residue_ptr != NULL)
                {
                  r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
                  atom_number = r_atom_ptr->atom_number - 1;
                }
            }

          /* Write out this molecule's coords */
          if (r_molec_ptr->nresid > 0)
            atom_number
              = r_write_pdb_coords(file_ptr,
                                   r_molec_ptr->first_r_residue_ptr,
                                   r_molec_ptr->nresid,atom_number,TRUE);
        }

      /* Get the next molecule */
      r_molec_ptr = r_molec_ptr->next_r_molec_ptr;
    }

  /* Close the output file */
  if (new_file == TRUE)
    fclose(file_ptr);
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

write_rasmol_headers  -  Write all the header records to the output
                         RasMol script file

***********************************************************************/

void write_rasmol_headers(FILE *outfile_ptr,char *pdb_code,char chain,
                          char *uniprot_acc,struct parameters params)
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  
  /* Write out the starting lines */
  fprintf(outfile_ptr,"#!rasmol -script\n");
  fprintf(outfile_ptr,
          "# Generated by genpymol - R A Laskowski, Jan 2017\n");
  fprintf(outfile_ptr,"\n");
  
  /* Show PDB code and PDB file name */
  fprintf(outfile_ptr,"# PDB code:   %s(%c)\n",pdb_code,chain);
  fprintf(outfile_ptr,"# UniProt id: %s\n",uniprot_acc);
  fprintf(outfile_ptr,"\n");
  
  /* Initilise and indicate that ATOM records to follow */
  if (params.output_type == RASMOL_SCRIPT)
    {
      fprintf(outfile_ptr,"zap\n");
      fprintf(outfile_ptr,"load inline\n");
     }
 
  /* Set background colour */
  strcpy(colour_name,"cream");
  
  /* Get the corresponding RasMol colour triplet */
  get_rasmol_colour(colour_name,colour_numbers);
  
  /* Write out the script line defining the background colour */
  fprintf(outfile_ptr,"background %s\n",colour_numbers);

  /* Other rendering parameters */
  if (params.output_type == RASMOL_SCRIPT)
    {
      fprintf(outfile_ptr,"set ambient 60\n");
      fprintf(outfile_ptr,"set specular off\n");
      fprintf(outfile_ptr,"slab off\n");
  
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
     }
 
  /* Show molecule details */
  fprintf(outfile_ptr,"# PDB details\n");
  fprintf(outfile_ptr,"echo \" \"\n");
  fprintf(outfile_ptr,"echo \"* PDB entry: %s(%c)\"\n",pdb_code,chain);
  if (uniprot_acc[0] != '\0')
    {
      fprintf(outfile_ptr,"echo \" \"\n");
      fprintf(outfile_ptr,"echo \"  showing variants from UniProt id %s\"\n",
              uniprot_acc);
    }
  fprintf(outfile_ptr,"echo \" \"\n");
  
  /* Switch off display */
  fprintf(outfile_ptr,"select all\n");
  fprintf(outfile_ptr,"wireframe off\n");
  fprintf(outfile_ptr,"spacefill off\n");
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"\n");
}
/***********************************************************************

colour_alpha_fold_rasmol  -  Use the AF_ranges file to assign the
                             confidence colours

***********************************************************************/

void colour_alpha_fold_rasmol(FILE *outfile_ptr,char *object,
                              struct parameters params)
{
  int factor, i, lower_b, upper_b;

  /* Initialise variables */
  if (params.output_type == JSMOL || params.output_type == JMOL_SCRIPT)
    factor = 1;
  else
    factor = 100;

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
    }
}
/***********************************************************************

write_rasmol_protein  -  Write script records defining protein chain

***********************************************************************/

int write_rasmol_protein(FILE *outfile_ptr,struct r_pdb *r_pdb_ptr,
                         struct c_file_data file_name_ptr,
                         struct parameters params)
{
  char add_on[11], chain, object[OBJ_NAME_LEN];
  char pdb_code[5], pdbsum_dir[FILENAME_LEN];
  char colour_coil[COLNAM_LEN + 5], colour_name[COLNAM_LEN];
  char major[10], minor[10];

  int atom_number, colour, dir_type, i, ichain, type;
  int first_atom_number, iatom, last_atom_number, natoms;
  int iresid, nresid;
  int split, use_atomno, wanted;

  double rgb[3], rgb_major[3], rgb_minor[3];

  struct r_atom *r_atom_ptr;
  struct r_chain *r_chain_ptr;
  struct r_molec *r_molec_ptr;
  struct r_residue *r_residue_ptr;


  /* Initialise variables */
  atom_number = 0;
  strcpy(pdb_code,r_pdb_ptr->pdb_code);
  split = TRUE;
  use_atomno = TRUE;

 /* Initialise strings */
  object[0] = '\0';
  add_on[0] = '\0';

  /* Get pointer to the first molecule in the PDB entry */
  r_molec_ptr = r_pdb_ptr->first_r_molec_ptr;

  /* Loop through all the molecules to calc coord limits */
  while (r_molec_ptr != NULL)
    {
      /* Get this molecule type */
      type = r_molec_ptr->type;

      /* If this is a protein chain, then need to check interaction */
      if (type == R_PROTEIN || type == R_CALPHA_ONLY)
        {
          /* Get this molecule's first residue */
          r_residue_ptr = r_molec_ptr->first_r_residue_ptr;
          nresid = r_molec_ptr->nresid;
          iresid = 0;

          /* Get this chain id */
          chain = r_residue_ptr->chain;
          r_chain_ptr = r_residue_ptr->r_chain_ptr;

          /* Initialise atom numbers */
          first_atom_number = last_atom_number = 0;

          /* Loop over the residues to renumber the atoms */
          while (r_residue_ptr != NULL && iresid < nresid)
            {
              /* Get this residue's first atom */
              r_atom_ptr = r_residue_ptr->first_r_atom_ptr;
              natoms = r_residue_ptr->natoms;
              iatom = 0;

              /* Loop through all the residue's atoms to renumber them */
              while (iatom < natoms && r_atom_ptr != NULL)
                {
                  /* Increment atom counter */
                  atom_number++;

                  /* If this is the chain's first atom, then save */
                  if (iresid == 0 && iatom == 0)
                    first_atom_number = atom_number;

                  /* Save chain's last astom number */
                  last_atom_number = atom_number;

                  /* Save the atom number */
                  r_atom_ptr->atom_number = atom_number;

                  /* Get pointer to next atom in linked-list */
                  r_atom_ptr = r_atom_ptr->next_r_atom_ptr;
                  iatom++;
                }

              /* Get the next residue */
              r_residue_ptr = r_residue_ptr->next_r_residue_ptr;
              iresid++;
            }

          /* Define chain name */
          sprintf(object,"chn%c",chain);

          /* If CA-only, add to name */
          add_on[0] = '\0';
          if (type == R_CALPHA_ONLY)
            strcpy(add_on," (CA-only)");

          /* Write comment line */
          fprintf(outfile_ptr,"# Chain %c - protein%s\n",chain,
                  add_on);

          /* Define this chain and then select it */
          fprintf(outfile_ptr,
                  "define %s atomno >= %d && atomno <= %d\n",object,
                  first_atom_number,last_atom_number);
          fprintf(outfile_ptr,"select %s & protein\n",object);

          /* Switch off wireframe */
          fprintf(outfile_ptr,"wireframe off\n");

          /* Show as cartoon */
          fprintf(outfile_ptr,"cartoon on\n");

          /* If not AlphaFold model, colour by chain */
          if (params.pdb_type != ALPHA_FOLD_PDB)
            {
              /* Get this chain's colour */
              get_chain_col(chain,rgb,rgb_major,rgb_minor,split,
                            colour_name);

              /* Define default chain colour */
              fprintf(outfile_ptr,"colour [%d,%d,%d] # (%s)\n",
                      (int) (255 * rgb[0]),(int) (255 * rgb[1]),
                      (int) (255 * rgb[2]),colour_name);

              /* Define secondary structure colours according to which
                 secondary structure type is the "dominant" one */
              if (r_chain_ptr->dominant_sst == R_STRAND)
                {
                  strcpy(major,"sheet");
                  strcpy(minor,"helix");
                }
              else
                {
                  strcpy(major,"helix");
                  strcpy(minor,"sheet");
                }

              /* Show dominant colour */
              fprintf(outfile_ptr,"select %s & protein & %s\n",object,
                      major);
              fprintf(outfile_ptr,"colour [%d,%d,%d] # (%s_major)\n",
                      (int) (255 * rgb_major[0]),
                      (int) (255 * rgb_major[1]),
                      (int) (255 * rgb_major[2]),colour_name);

              /* Repeat for non-dominant colour */
              fprintf(outfile_ptr,"select %s & protein & %s\n",object,
                      minor);
              fprintf(outfile_ptr,"colour [%d,%d,%d] # (%s_minor)\n",
                      (int) (255 * rgb_minor[0]),
                      (int) (255 * rgb_minor[1]),
                      (int) (255 * rgb_minor[2]),colour_name);

              /* Blank line */
              fprintf(outfile_ptr,"\n");
              fflush(outfile_ptr);
            }

          /* If this is an Alpha Fold model, use the AF_ranges file to
             assign the confidence colours */
          else
            colour_alpha_fold_rasmol(outfile_ptr,object,params);
        }

      /* Get the next molecule */
      r_molec_ptr = r_molec_ptr->next_r_molec_ptr;
    }

  /* Return flag */
  return(use_atomno);
}
/***********************************************************************

write_rasmol_nucleic  -  Write script records defining nucleic acid
                         chains

***********************************************************************/

void write_rasmol_nucleic(FILE *outfile_ptr,struct r_pdb *r_pdb_ptr)
{
  char chain, object[OBJ_NAME_LEN];
  char colour_coil[COLNAM_LEN + 5], colour_name[COLNAM_LEN];

  int colour, i, ichain, type;
  int split, wanted;

  double rgb[3], rgb_major[3], rgb_minor[3];

  struct r_chain *r_chain_ptr;
  struct r_residue *r_residue_ptr;

  /* Initialise variables */
  split = TRUE;

  /* Initialise strings */
  object[0] = '\0';

  /* Get pointer to the first chain */
  r_chain_ptr = r_pdb_ptr->first_r_chain_ptr;

  /* Loop through all the chains to render the protein chains */
  while (r_chain_ptr != NULL)
    {
      /* Get chain id and type */
      chain = r_chain_ptr->chain_id;
      type = r_chain_ptr->type;

      /* If protein or DNA, create a name for this chain */
      if (type == R_NUCLEIC_ACID)
        {
          /* Write comment line */
          fprintf(outfile_ptr,"# Chain %c - nucleic acid\n",chain);
          sprintf(object,"chn%c",chain);

          /* Define this chain and then select it */
          fprintf(outfile_ptr,"define %s *:%c\n",object,chain);
          fprintf(outfile_ptr,"select %s & !hydrogen\n",object);

          /* Show as wireframe */
          fprintf(outfile_ptr,"wireframe 0.2\n");
          fprintf(outfile_ptr,"colour atoms cpk\n");

          /* Select phosphate backbone */
          fprintf(outfile_ptr,"select %s & *.P\n",object);

          /* Get this chain's colour */
          get_chain_col(chain,rgb,rgb_major,rgb_minor,split,colour_name);

          /* Define default chain colour */
          fprintf(outfile_ptr,"colour atoms [%d,%d,%d] # (%s)\n",
                  (int) (255 * rgb[0]),(int) (255 * rgb[1]),
                  (int) (255 * rgb[2]),colour_name);

          /* Switch backbone mode on */
          fprintf(outfile_ptr,"backbone on\n");
        }

      /* Get the next chain */
      r_chain_ptr = r_chain_ptr->next_r_chain_ptr;
    }
}
/***********************************************************************

check_ligand_ranges  -  Check for excessively long ligand ranges and
                        split if necessary

***********************************************************************/

void check_ligand_ranges(struct r_pdb *r_pdb_ptr)
{
  int add_number, last_number, len, nligands, nres, number, outlen;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  add_number = 0;
  last_number = -1;
  len = outlen = 0;
  nligands = 0;
  nres = 0;

  /* Get the first rasdef record */
  r_rasdef_ptr = r_pdb_ptr->first_r_rasdef_ptr;

  /* Loop over all the rasdef records to pick out all the ligand
     residue sequences */
  while (r_rasdef_ptr != NULL)
    {
      /* If ligand residue, then write out  */
      if (r_rasdef_ptr->type == R_LIGAND)
        {
          /* Increment residue count */
          nres++;

          /* Get the ligand number */
          number = r_rasdef_ptr->number;

          /* Get the length of the ligand range */
          len = strlen(r_rasdef_ptr->ligand_range) + 2;

          /* If this will make the output string too long for RasMol,
             need to break up */
          if (len + outlen > MAX_LIGRANGE)
            {
              /* Increment addition */
              add_number++;

              /* Reset the output length */
              outlen = 0;
            }

          /* Update the range number */
          r_rasdef_ptr->number = r_rasdef_ptr->number + add_number;

          /* Update the output length */
          outlen = outlen + len;
        }

      /* Get the next record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }
}
/***********************************************************************

display_ligand  -  Write out the RasMol commands to display given ligand

***********************************************************************/

void display_ligand(FILE *outfile_ptr,char *object)
{
  /* Select and display the ligand */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"select %s & !hydrogen\n",object);
  fprintf(outfile_ptr,"wireframe 0.2\n");
  fprintf(outfile_ptr,"spacefill on\n");
  fprintf(outfile_ptr,"colour atoms cpk\n");
  fprintf(outfile_ptr,"select %s & carbon\n",object);
  fprintf(outfile_ptr,"colour [56,56,70]\n");
  fprintf(outfile_ptr,"\n");

  /* Re-initialise object name */
  object[0] = '\0';
}
/***********************************************************************

write_rasmol_ligands  -  Write script records to display ligands

***********************************************************************/

void write_rasmol_ligands(FILE *outfile_ptr,struct r_pdb *r_pdb_ptr)
{
  char chain, object[OBJ_NAME_LEN], range_name[20], residue_range[20];

  int last_number, nligands, nres, number;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  last_number = -1;
  nligands = 0;
  nres = 0;
  object[0] = '\0';

  /* Get the first rasdef record */
  r_rasdef_ptr = r_pdb_ptr->first_r_rasdef_ptr;

  /* Loop over all the rasdef records to pick out all the ligand
     residue sequences */
  while (r_rasdef_ptr != NULL)
    {
      /* If ligand residue, then write out  */
      if (r_rasdef_ptr->type == R_LIGAND)
        {
          /* Get the ligand number */
          number = r_rasdef_ptr->number;

          /* If this is a new number, close off the previous ligand
             and start a new one */
          if (number != last_number)
            {
              /* If have a previous ligand, then display it */
              if (object[0] != '\0')
                display_ligand(outfile_ptr,object);

              /* Show ligand name */
              fprintf(outfile_ptr,"# Ligand %d - %s\n",number,
                      r_rasdef_ptr->ligand_name);

              /* Create the object record for the new ligand */
              sprintf(object,"lig%4.4d",number);

              /* Commence its definition */
              fprintf(outfile_ptr,"define %s ",object);

              /* Save the ligand number */
              last_number = number;

              /* Initialise count of residues in the ;igand */
              nres = 0;
            }

          /* If not first residue of the ligand, append a comma */
          if (nres > 0)
            fprintf(outfile_ptr,", ");

          /* Add current residue to the ligand definition */
          fprintf(outfile_ptr,"%s",r_rasdef_ptr->ligand_range);

          /* Increment residue count */
          nres++;
        }

      /* Get the next record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

  /* If have a last ligand, then display it */
  if (object[0] != '\0')
    display_ligand(outfile_ptr,object);
}
/***********************************************************************

write_rasmol_metals  -  Write script records to display metals

***********************************************************************/

void write_rasmol_metals(FILE *outfile_ptr,struct r_pdb *r_pdb_ptr)
{
  char chain, colour_numbers[COLNO_LEN], object[OBJ_NAME_LEN];
  char full_name[8], metal_name[5], tmp_name[5];

  int ipos, nmetals;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  nmetals = 0;
  object[0] = '\0';

  /* Get the first rasdef record */
  r_rasdef_ptr = r_pdb_ptr->first_r_rasdef_ptr;

  /* Loop over all the rasdef records to pick out all the metal
     residue sequences */
  while (r_rasdef_ptr != NULL)
    {
      /* If metal residue, then write out  */
      if (r_rasdef_ptr->type == R_METAL)
        {
          /* Increment count of metals */
          nmetals++;

          /* Get the residue name */
          strcpy(metal_name,r_rasdef_ptr->ligand_name);

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

          /* Form the metal name */
          strcpy(full_name,"*.");
          strcat(full_name,metal_name);

          /* Deal with special cases: Ca and Cd */
          if (!strncmp(metal_name,"CA",2))
            strcpy(full_name,"CA.CA");
          if (!strncmp(metal_name,"Ca",2))
            strcpy(full_name,"CA.Ca");
          if (!strncmp(metal_name,"CD",2))
            strcpy(full_name,"CD");

          /* Special case of single-letter metals such as K, U, etc */
          if (metal_name[0] == '_')
            {
              full_name[0] = metal_name[1];
              full_name[1] = '\0';
            }

          /* Show metal name */
          fprintf(outfile_ptr,"# Metal %d - %s\n",nmetals,
                  r_rasdef_ptr->ligand_name);

          /* Create the object record for this metal */
          sprintf(object,"metal%s",metal_name);

          /* Define metal */
          fprintf(outfile_ptr,"define %s %s\n",object,full_name);

          /* Get the corresponding RasMol colour triplet */
          get_rasmol_colour(r_rasdef_ptr->colour,colour_numbers);
  
          /* Select and display the metal */
          fprintf(outfile_ptr,"select %s & !hydrogen\n",object);
          fprintf(outfile_ptr,"colour atoms %s\n",colour_numbers);
          fprintf(outfile_ptr,"spacefill 1.5\n");

          /* Blank line */
          fprintf(outfile_ptr,"\n");
        }

      /* Get the next record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }
}
/***********************************************************************

write_def  -  Write the variant definition to the RasMol script

***********************************************************************/

void write_def(FILE *outfile_ptr,char *defn_string,int type,int num)
{

  /* Close off the definition string */
  strcat(defn_string,")");

  /* Write out the definition */
  fprintf(outfile_ptr,"%s\n",defn_string);
  fprintf(outfile_ptr,"select var_type%d_%d\n",type + 1,num);
      
  /* Define this variant's colour */
  fprintf(outfile_ptr,"colour %s\n",PYMOL_COLOUR[type]);

  /* Show as sticks */
  fprintf(outfile_ptr,"wireframe 0.2\n");
}
/***********************************************************************

write_rasmol_variants  -  Write the variants to the RasMol script

***********************************************************************/

void write_rasmol_variants(FILE *outfile_ptr,
                           struct residue **residue_index_ptr,int nresid,
                           int labels)
{
  char defn_string[RASMOL_LINELEN + 1], residue_def[30], string[30];
  char atom_def[30];

  int iresid, len, nres, num, nvar_res, type;

  int *var_type;

  struct residue *residue_ptr;
  struct r_residue *r_residue_ptr;

  /* Write out comment */
  fprintf(outfile_ptr,"# Variants\n");

  /* Determine each residue's 'highest' variant type */
  nvar_res = get_max_variants(residue_index_ptr,nresid,&var_type);

  /* Loop over the variant types, rendering the list of residues
     belonging to each type */
  for (type = 0; type < NVAR_TYPES && nvar_res > 0; type++)
    {
      /* Initialise count of residues belonging to this type */
      nres = num = 0;

      /* Initialise string holding list of this type's residues */
      defn_string[0] = '\0';

      /* Loop over all the residues to pick out the ones of this type */
      for (iresid = 0; iresid < nresid; iresid++)
        {
          /* If residue of the correct type, then add to definition */
          if (var_type[iresid] == type)
            {
              /* Check that the definition string is not too long */
              len = strlen(defn_string);
              if (len > RASMOL_LINELEN - 10)
                {
                  /* Write out the definition we have so far */
                  write_def(outfile_ptr,defn_string,type,num);

                  /* Re-initialise */
                  nres = 0;
                  defn_string[0] = '\0';
                }
                  
              /* Get this residue */
              residue_ptr = residue_index_ptr[iresid];

              /* Get the corresponding PDB residue */
              r_residue_ptr = residue_ptr->r_residue_ptr;
          
              /* Convert residue number into RasMol format */
              sprintf(residue_def,"%s",r_residue_ptr->res_num);

              /* If first one, start definition */
              if (nres == 0)
                {
                  /* Increment count for this definition */
                  num++;

                  /* Start definition */
                  sprintf(defn_string,"define var_type%d_%d ",type + 1,num);

                  /* If have chain is, include that in the definition */
                  if (r_residue_ptr->chain != ' ')
                    {
                      /* Include chain id */
                      sprintf(string,"*:%c & ",r_residue_ptr->chain);
                      strcat(defn_string,string);
                    }
                  
                  /* Open brackets for list of residue numbers */
                  strcat(defn_string,"(");
                }
              else
                strcat(defn_string," | ");

              /* Add residue to list */
              strcat(defn_string,residue_def);

              /* Increment count of residues of this type */
              nres++;

              /* If labelling variant, write out the label */
              if (labels == TRUE)
                {
                  /* Form name of atom */
                  sprintf(atom_def,"%s%s:%c.%s",
                          r_residue_ptr->res_name,
                          r_residue_ptr->res_num,
                          r_residue_ptr->chain,
                          residue_ptr->last_atom_name);

                  /* Remove any blanks from the definition */
                  remove_blanks(atom_def);
      
                  /* Select the atom */
                  fprintf(outfile_ptr,"select %s\n",atom_def);
                  fprintf(outfile_ptr,"label %%n%%r\n");
                }
            }
        }

      /* If we have any residues of this type, render them */
      if (nres > 0)
        {
          /* Write out the definition */
          write_def(outfile_ptr,defn_string,type,num);
        }
    }

  /* Free the var_type array */
  free(var_type);
}
/***********************************************************************

write_rasmol_rescons  -  Colour according to residue conservation (ie
                         B-factor)

***********************************************************************/

void write_rasmol_rescons(FILE *outfile_ptr,struct parameters params)
{
#define NRANGES       7

  char chain;

  static char *RANGE_COLOUR[NRANGES] = { 
    "purple", "blue", "cyan", "green", "yellow", "orange", "red" };

  int diff, irange, value;

  /* Initialise variables */
  value = 0;
  diff = (int) (8000.0 / NRANGES);

  /* Write out the commands to colour by B-value */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Residue conservation\n");

  /* Select the appropriate chain and set to wireframe */
  chain = params.chain;
  fprintf(outfile_ptr,"select *:%c\n",chain);
  fprintf(outfile_ptr,"cartoon off\n");
  fprintf(outfile_ptr,"wireframe on\n");
  fprintf(outfile_ptr,"backbone 0.2\n");

  /* Loop over the residue conservation ranges to colour accordingly */
  for (irange = 0; irange < NRANGES; irange++)
    {
      /* Write out the colour for this range */
      fprintf(outfile_ptr,"select *:%c & temperature > %d\n",chain,
              value);
      fprintf(outfile_ptr,"colour %s\n",RANGE_COLOUR[irange]);

      /* On highest range, make the residues thicker */
      if (irange == NRANGES - 1)
        fprintf(outfile_ptr,"wireframe 0.1\n");

      /* Increment the value */
      value = value + diff;
    }

  /* Show dot surface */
  fprintf(outfile_ptr,"select *:%c & protein\n",chain);
  fprintf(outfile_ptr,"dots on\n");
}
/***********************************************************************

write_rasmol_domcol  -  Colour according to Pfam domain

***********************************************************************/

void write_rasmol_domcol(FILE *outfile_ptr,struct domain *first_domain_ptr,
                         struct parameters params)
{
  char chain, colour_string[20], object[OBJ_NAME_LEN], res_range[20];
  char chain_select[20], colour_name[20];

  int icol, ndomain, start_col;

  struct domain *domain_ptr;

  /* Initialise variables */
  ndomain = 0;
  chain_select[0] = '\0';
  if (params.chain != '*')
    sprintf(chain_select," & *:%c",params.chain);

  /* Get the starting domain colour */
  start_col = PFAMA_STARTCOL;

  /* If this is a covid-19 protein, define corresponding start colour */
  if (params.start_colour > 0)
    start_col = params.start_colour;

  /* Write out heading */
  fprintf(outfile_ptr,"\n");
  fprintf(outfile_ptr,"# Colour by domain\n");

  /* Select the appropriate chain and colour white */
  chain = params.chain;
  fprintf(outfile_ptr,"select protein\n");
  fprintf(outfile_ptr,"colour white\n");
  fprintf(outfile_ptr,"select *:%c & protein\n",params.chain);
  fprintf(outfile_ptr,"colour grey\n");
  fprintf(outfile_ptr,"\n");

  /* Get pointer to first domain */
  domain_ptr = first_domain_ptr;

  /* Loop over all the domains to colour */
  while (domain_ptr != NULL)
    {
      /* If have sytart and end residues, define the domain and its
         colour */
      if (domain_ptr->start_r_residue_ptr != NULL &&
          domain_ptr->end_r_residue_ptr != NULL)
        {
          /* Increment domain number */
          ndomain++;

          /* Form domain name */
          sprintf(object,"domain%4.4d",ndomain);

          /* Determine this domain's colour */
          icol = (start_col + domain_ptr->colour) % (NPFAM_COLOURS - 2);
          sprintf(colour_string,"[%s]",COLOUR_DEFN[icol][2]);
          strcpy(colour_name,COLOUR_DEFN[icol][0]);

          /* Replace spaces by commas */
          replace_char(colour_string,' ',',');

          /* Write out the domain name and colour */
          fprintf(outfile_ptr,"# Domain %d - %s = %s (%s)\n",
                  ndomain,domain_ptr->pfam_id,domain_ptr->name,
                  colour_name);

          /* Form the names of the start and end residues */
          sprintf(res_range,"%s-%s",
                  domain_ptr->start_r_residue_ptr->res_num,
                  domain_ptr->end_r_residue_ptr->res_num);

          /* Remove any blanks from the definition */
          remove_blanks(res_range);

          /* Show domain definition */
          fprintf(outfile_ptr,"echo \"Domain %d: %s - %s (%s)\"\n",
                  ndomain,domain_ptr->pfam_id,domain_ptr->name,
                  colour_name);

          /* Define this domain */
          fprintf(outfile_ptr,"define %s %s%s\n",object,res_range,
                  chain_select);

          /* Select it */
          fprintf(outfile_ptr,"select %s\n",object);

          /* Set the colour */
          fprintf(outfile_ptr,"colour %s # (%s)\n",colour_string,
                  colour_name);
          fprintf(outfile_ptr,"\n");
        }

      /* Get the next domain record */
      domain_ptr = domain_ptr->next_domain_ptr;
    }
}
/***********************************************************************

write_rasmol_script  -  Write command to RasMol script to define all objects
                     and their representation

***********************************************************************/

void write_rasmol_script(struct r_pdb *r_pdb_ptr,
                         struct residue **residue_index_ptr,int nresid,
                         struct residue *ref_residue_ptr,
                         struct domain *first_domain_ptr,
                         struct c_file_data file_name_ptr,
                         struct parameters params)
{
  char file_name[FILENAME_LEN + 1], number_string[20];
  char colour_name[COLNAM_LEN], display_type[20];
  char colour_coil[COLNAM_LEN + 5];

  int iresid, metal_type, nsite, number, type, valid;
  int colour_as_orig, use_atomno;

  FILE *file_ptr;

  /* Initialise variables */
  colour_as_orig = TRUE;

  /* Form name of rasmol script file, rasmol.pml */
  strcpy(file_name,"genpymol.script");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);
  if (file_ptr == NULL)
    {
      printf("*** ERROR. Unable to open output file %s\n",file_name);
      exit(1);
    }

  /* Write out the header records */
  write_rasmol_headers(file_ptr,params.pdb_code,params.chain,
                       params.uniprot_acc,params);

  /* Write out the definitions for all the protein chains */
  use_atomno
    = write_rasmol_protein(file_ptr,r_pdb_ptr,file_name_ptr,params);

  /* Write out the definitions for all the DNA/RNA chains */
  write_rasmol_nucleic(file_ptr,r_pdb_ptr);

  /* Check for excessively long ligand ranges and split if necessary */
  check_ligand_ranges(r_pdb_ptr);

  /* Write out the definitions for all the ligands */
  write_rasmol_ligands(file_ptr,r_pdb_ptr);

  /* Write out the definitions for all the metals */
  write_rasmol_metals(file_ptr,r_pdb_ptr);

  /* Write out any variants */
  if (params.variants == TRUE)
    write_rasmol_variants(file_ptr,residue_index_ptr,nresid,
                          params.labels);

  /* Colour according to residue conservation, ir required */
  if (params.res_cons == TRUE)
    write_rasmol_rescons(file_ptr,params);

  /* If colouring by domain, use the domain definitions */
  if (params.colour_by_domain == TRUE)
    write_rasmol_domcol(file_ptr,first_domain_ptr,params);

  /* Rotate by 180 about x-axis to be in same orientation as in PDB
     file */
  if (params.output_type != RASMOL_SCRIPT)
    fprintf(file_ptr,"rotate x 180\n");

  /* Close off the script commands */
  fprintf(file_ptr,"\n");
  fprintf(file_ptr,"exit\n");
  fprintf(file_ptr,"\n");

  /* Append the atom coordinates to the file */
  if (params.output_type == RASMOL_SCRIPT)
    write_coords(r_pdb_ptr,file_ptr,use_atomno);

  /* Close the Rasmol script file */
  fclose(file_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char chain, log_name[FILENAME_LEN], pdbsum_dir[FILENAME_LEN];
  char pdir[FILENAME_LEN];

  char *code_dir, *pdb_sequence, *seq_name, *sequence;

  int nchains, npdb_resid, nresid, nwanted, nvar, nvar_res, seq_len;
  int dir_type, ndomains, add_res;
  int have_conservation, use_atomno;

  float max_cons;

  struct c_file_data file_name_ptr;
  struct parameters params;

  struct chain *first_chain_ptr;
  struct domain *first_domain_ptr;
  struct pobject *first_pobject_ptr;
  struct r_pdb *r_pdb_ptr;
  struct residue *first_residue_ptr, *ref_residue_ptr;
  struct r_residue *first_r_residue_ptr;

  struct residue **residue_index_ptr;

  FILE *logfile_ptr, *pymol_script_ptr;

  /* Initialise variables */
  chain = '\0';
  strcpy(log_name,"genpymol.out");
  have_conservation = FALSE;
  max_cons = 0;
  ndomains = 0;
  pdir[0] = '\0';
  pymol_script_ptr = NULL;
  add_res = 0;
  use_atomno = FALSE;

  /* Initialise pointers */
  first_domain_ptr = NULL;
  first_residue_ptr = ref_residue_ptr = NULL;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Read in all the directory paths from the parameter file, CATHPARAM */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Open the log file */
  logfile_ptr = open_output_file(log_name,TRUE);

  /* For PDBsum1, get the PDBsum1 directory for this entry */
  if (params.pdb_type == PDBSUM1)
    get_pdbsum1_dir(pdir,params.pdb_code,file_name_ptr);

  /* Find the PDBsum directory for this entry */
  else
    dir_type
      = find_pdbsum_dir(params.pdb_code,file_name_ptr.pdbsum_template,
                        params.pdb_type,pdbsum_dir);

  /* Read in the equivalent chains */
  nchains = read_chain_equivs(&first_chain_ptr,pdbsum_dir,
                              params.pdb_type);

  /* If have a UniProt accession, pick up the associated data */
  if (params.uniprot_acc[0] != '\0')
    {
      /* Get the VarSite directory */
      get_data_dir(&code_dir,
                   file_name_ptr.other_dir[VARSITE_REAL_UNIPROT_DIR],
                   params.uniprot_acc);

      /* Get the chain idebtifier */
      chain = params.chain;

      /* Pick up the closest PDB file for the given residue */
      if (params.pdb_code[0] == '\0')
        params.start
          = get_closest_pdb(params.pdb_code,&chain,params.uniprot_acc,
                            params.start,code_dir);

      /* Read in the residue conservation scores */
      nresid
        = get_res_conservation(code_dir,&first_residue_ptr,&max_cons,
                               logfile_ptr);

      /* Store the residue pointers in an index array */
      nresid = create_residue_index(first_residue_ptr,nresid,
                                    &residue_index_ptr);

      /* Get all the residue mappings to this PDB file */
      if (params.pdb_type != ALPHA_FOLD_PDB)
        ref_residue_ptr
          = get_pdb_mappings(residue_index_ptr,nresid,code_dir,
                             params.pdb_code,chain,params.start,
                             first_chain_ptr);

      /* For AlphaFold model, get the start residue from its database.dat
         file */
      //else
      //add_res = read_database_dat(params.pdb_code,pdbsum_dir);
      add_res = 0;

      /* Show reference residue found */
      if (ref_residue_ptr != NULL)
        {
          fprintf(logfile_ptr,"Reference residue: %d %c  PDB: %s %c - %s %s\n",
                  ref_residue_ptr->seq_no,
                  ref_residue_ptr->aa_code,
                  ref_residue_ptr->pdb_code,
                  ref_residue_ptr->pdb_chain,
                  ref_residue_ptr->pdb_res_name,
                  ref_residue_ptr->pdb_res_num);
          fflush(logfile_ptr);
        }
    }

  /* Read in the required PDB entry */
  r_pdb_ptr = r_get_pdb_data(params.pdb_code,file_name_ptr,params.pdb_type,
                             COORDS_ONLY);

  /* If PDB file not found, then abort */
  if (r_pdb_ptr == NULL)
    {
      fprintf(logfile_ptr,"*** ERROR. PDB file for %s not found\n",
              params.pdb_code);
      exit(1);
    }


  /* Read in the secondary structure assignments from the resdata.dat file */
  npdb_resid = s_get_resdata(params.pdb_code,chain,params.pdb_type,
                             &first_r_residue_ptr,file_name_ptr);

  /* Transfer assignments to the PDB residues */
  s_transfer_resdata(r_pdb_ptr->first_r_residue_ptr,first_r_residue_ptr);

  /* If this is a VarSite or covid-19 protein, map the residues to the
     PDB structure */
  if (params.uniprot_acc[0] != '\0')
    {
      /* Map to the PDB structure just read in */
      map_to_pdb(residue_index_ptr,nresid,r_pdb_ptr->first_r_residue_ptr,
                 add_res,max_cons,params.res_cons);

      /* If we have a reference residue, delete all but the neighbouring
         residues - not wanted at present */
      ref_residue_ptr = NULL;
      if (ref_residue_ptr != NULL)
        {
          /* Remove all but the local residues around our reference */
          nwanted = mark_local_residues(r_pdb_ptr->first_r_residue_ptr,
                                        ref_residue_ptr->r_residue_ptr,
                                        params.cutoff_dist);

          /* Delete unwanted residues */
          npdb_resid = delete_residues(r_pdb_ptr);
        }

      /* If showing VarSite variants, then read them in from the appropriate
         map file */
      if (params.variants == TRUE)
        {
          /* If seq start and end not defined, use first and last residues */
          if (params.start == 0 || params.end == 0)
            {
              params.start = 1;
              params.end = seq_len;
            }

          /* No longer needed */

          /* Get the UniProt sequence from the sprot.fasta file */
          //seq_len
          //  = read_seq_fasta(code_dir,"sprot.fasta",&seq_name,&sequence);

          /* Read in the residue conservation scores */
          //          nresid
          //  = get_res_conservation(code_dir,&first_residue_ptr,&max_cons,
          //                         logfile_ptr);

          /* Create a residue record for every position in the query
             sequence */
          //if (nresid == 0)
          //  nresid
          //    = initialise_residues(&first_residue_ptr,sequence,logfile_ptr);

          /* Store the residue pointers in an index array */
          //nresid
          //  = create_residue_index(first_residue_ptr,nresid,
          //                         &residue_index_ptr);

          /* Pick up the mutation details from the uprotein.dat file */
          nvar_res
            = get_mutations(params.uniprot_acc,residue_index_ptr,nresid,
                            code_dir,logfile_ptr);

          /* If we have any added variants, add then to the relevant
             residues */
          nvar_res
            = added_variants(residue_index_ptr,nresid,nvar_res,params);

          /* Sort each residue's variants by importance */
          sort_variants(residue_index_ptr,nresid);

          /* No longer needed */

          /* Scan the align.aln files to map the residues to our PDB
             structure */
          //          map_residues(code_dir,residue_index_ptr,nresid,r_pdb_ptr,
          //             params.chain,nvar_res,params.start,params.end);

          /* Write the variants out to the variants.map and variants.dat
             files */
          write_variants(residue_index_ptr,nresid,params);
        }
    }

  /* If colouring by domain, pick up all the Pfam domains for this
     sequence from the PDB entry's corresponding uniplot.dat file */
  if (params.colour_by_domain == TRUE)
    {
      /* Read in the Pfam domains */
      ndomains
        = read_uniplot_dat(params.pdb_code,params.chain,pdbsum_dir,
                           &first_domain_ptr,&sequence,&pdb_sequence);

      /* Get the start and end PDB residues for each domain */
      map_domains_to_pdb(first_domain_ptr,sequence,pdb_sequence,
                         r_pdb_ptr,params.chain);
    }

  /* If only want neighbouring molecules, delete any not in contact with
     the given chain */
  if (params.neighbours == TRUE)
    identify_neighbours(r_pdb_ptr,params.chain,file_name_ptr,params);

  /* Write the PyMOL script */
  if (params.output_type == PYMOL_SCRIPT)
    write_pymol_pml(r_pdb_ptr,&pymol_script_ptr,&first_pobject_ptr,
                    residue_index_ptr,nresid,first_domain_ptr,pdir,
                    file_name_ptr,params);

  /* Write the RasMol script */
  else if (params.output_type == RASMOL_SCRIPT ||
           params.output_type == JMOL_SCRIPT ||
           params.output_type == JSMOL)
    write_rasmol_script(r_pdb_ptr,residue_index_ptr,nresid,
                        ref_residue_ptr,first_domain_ptr,file_name_ptr,
                        params);

  /* Write out the coordinates */
  if (params.output_type != RASMOL_SCRIPT &&
      params.pdb_type != PDBSUM1)
    write_coords(r_pdb_ptr,NULL,use_atomno);

  /* Close the log file */
  fclose(logfile_ptr);
  
  exit(1);
}
