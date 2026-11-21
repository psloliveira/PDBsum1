/***********************************************************************

  uniplot.c - Program to plot the parts of the whole UniProt sequence
              that correspond to the solved protein chain in PDB, adding
              PFam domains where known

************************************************************************


Testing in $SCRATCH_DIR_NEW/run1 on 1aog

Date:-         21 Sep 2006
Dir update:-   19 Nov 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name       Description
-----------       -----------
appdrugs_pdb.dat  Input file containing the Het Group to DrugPort mappings
chains.dat        Input chains file listing unique and duplicate
                  protein chains
dataset.lst       Input list of PDB codes to be processed
letters.dat       Input file containing all the character definitions
lig2pfam.dat      Output file listing ligands and the Pfam domains to
                  which they are bound in the PDB
lig2pfam.idx      Output index file for the above
lig2pfam.txt      Output text file for the above
pfam2lig.dat      Output file listing Pfam domains and the ligands that
                  bind to them in the PDB
pfam2lig.idx      Output index file for the above
pfam2lig.txt      Output text file for the above
pfam_domains.parsed.sorted Input file giving short names and descriptions
                  of all PfamA entries
rasmol.dfn        Input definition file for each PDB entry
resdata.dat       Input file listing all the protein residues and
                  their annotations
seqXXX.fasta      Input FASTA-format file for each unique sequence of
                  every PDB entry
seqinfo.dat       Input file in PDBsum directory giving mapping of each
                  protein chain to the corresponding UniProt sequence
swisspfam.gz      Input file mapping Pfam domains onto UniProt sequences
uniprot_sprot.fasta.gz  Gzipped UniProt sequence database
uniprot_trembl.fasta.gz    "       "       "        "
uniprot.fasta     Output UniProt sequence for comparison against
                  corresponding PDB sequence
grow.out          Input file of interactions made by protein
pdbsum.dr         Output file of DR lines for each UniProt code to PDB
                  chain equivalence
uniplot.dat       Output file giving alignment to UniProt seq and any
                  Pfam domains

Compilation
-----------

cc -c uniplot.c
cc -c common.c
cc -c align.c
cc -c reapar.c
cc -c readpdb.c
cc -o uniplot uniplot.o common.o align.o reapar.o readpdb.o -lm -lz

Optimization flags:-

cc -O2 -funroll-loops -c uniplot.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c align.c
cc -O2 -funroll-loops -c reapar.c
cc -O2 -funroll-loops -c readpdb.c
cc -o uniplot uniplot.o common.o align.o reapar.o readpdb.o -lm -lz

For testing on laptop, without read of gzipped files, use:
    gcc -DTEST -c uniplot.c

For testing without system calls, use:
    gcc -DTEST -DNOSYST -c uniplot.c

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)
   -> get_letters
         -> string_truncate (in common.c)
   -> get_hetdrug_mappings
         -> get_tokens (in common.c)
         -> create_hetdrug_entry
   -> get_pdb_codes
         -> string_truncate (in common.c)
         -> create_pdb_entry
   -> process_codes
         -> find_pdbsum_dir (in common.c)
         -> get_chains_dat
         -> get_chains_dat
               -> string_truncate (in common.c)
               -> create_chain_record
         -> get_seqinfo_dat
               -> string_truncate (in common.c)
               -> find_code_record
               -> create_code_record
               -> create_clist_record
               -> create_smatch_record
   -> read_old_codes
         -> string_truncate (in common.c)
   -> create_code_index
         -> code_sort
   -> read_swisspfam
         -> string_chomp (in common.c)
         -> binary_search
         -> get_tokens (in common.c)
         -> find_pfam_entry
         -> create_pfam_entry
         -> create_domain_entry
         -> free_tokens (in common.c)
   -> get_missing_domains
         -> perform_splice (in common.c)
         -> empty_check (in common.c)
         -> run_wget (in common.c)
               -> empty_check (in common.c)
         -> read_pfam_xml -- No longer used since Pfam shut down
               -> store_string (in common.c)
               -> find_pfam_entry
               -> create_pfam_entry
               -> create_domain_entry
               -> get_tokens (in common.c)
               -> free_tokens (in common.c)
         -> get_pfam_tsv
               -> read_pfam_tsv
                     -> string_chomp (in common.c)
                     -> get_tokens (in common.c)
                     -> store_string (in common.c)
                     -> create_domain_entry
                     -> free_tokens (in common.c)
               -> find_pfam_entry
               -> create_pfam_entry
   -> create_pfam_index
         -> pfam_sort
   -> get_pfam_names
         -> string_chomp (in common.c)
         -> tokenize (in common.c)
   -> download_uniprot_seqs
         -> empty_check (in common.c)
         -> run_wget (in common.c)
               -> empty_check (in common.c)
         -> align_seqs
               -> get_sequence
                     -> string_chomp (in common.c)
                     -> get_uniprot_data
                     -> store_string (in common.c)
                     -> append_string (in common.c)
               -> align_sequences (in align.c)
   -> process_alignments
         -> compare_seqs (as below)
         -> finish_entries (as below)
   -> delete_fasta_files
   -> process_uniprot_seqs
         -> string_truncate (in common.c)
         -> compare_seqs
               -> find_pdbsum_dir (in common.c)
               -> colour_pfam_domains
                     -> create_domain_entry
                     -> order_domains
                           -> domain_sort
               -> run_fasta
               -> read_fasta_out
                     -> store_seq_pair
               -> read_sequence
               -> get_resdata
                     -> create_residue_record
                     -> get_tokens (in common.c)
                     -> create_ligand_record
                           -> store_string (in common.c)
                     -> free_tokens (in common.c)
               -> get_c_scores (in common.c)
                     -> create_c_score_entry (in common.c)
               -> match_to_residues
               -> r_get_dfn_data (in readpdb.c)
                     -> form_filename (in common.c)
                           -> store_string (in common.c)
                     -> string_chop (in common.c)
                     -> create_r_rasdef_entry (in readpdb.c)
                     -> r_get_ligand_residues (in readpdb.c)
                           -> r_get_token (in readpdb.c)
                           -> r_format_res_number (in readpdb.c)
                           -> r_remove_leading_blanks (in readpdb.c)
                           -> create_r_rasdef_entry (in readpdb.c)
                           -> store_string (in common.c) (in readpdb.c)
                     -> r_to_lower (in readpdb.c)
                     -> store_string (in common.c)
                     -> r_number_ligands (in readpdb.c)
                     -> map_to_ligplots (in readpdb.c)
                           -> string_chomp (in common.c)
                           -> get_tokens (in common.c)
                           -> split_metal_range (in readpdb.c)
                           -> create_r_rasdef_entry (in readpdb.c)
                           -> r_remove_leading_blanks (in common.c)
                           -> store_string (in common.c)
               -> identify_ligands
                     -> add_ligand
                           -> find_liglist_entry
                           -> create_liglist_entry
                     -> find_pfamlig_entry
                     -> create_pfamlig_entry
               -> find_residue
               -> flag_ligands
               -> unique_ligands
               -> write_ligands
               -> multiple_domains
                     -> create_pfamlig_entry
                     -> create_plist_entry
                     -> sort_domains
               -> free_rasdefs
               -> identify_cath_domains
                     -> find_cath_record
                     -> create_cath_record
                     -> create_domain_entry
               -> match_cath_domains
                     -> create_domain_entry
               -> plot_miniwir
                     -> get_domain_colour
                     -> generate_image
                           -> get_chain_colour (in common.c)
                                 -> get_colour_rgb (in common.c)
                           -> get_pfam_colour
                           -> open_output_file (in common.c)
               -> free_residue
               -> create_smatch_record
                     -> store_string (in common.c)
         -> finish_entries
               -> write_uniplot_dat
                     -> write_int_ligands
                           -> check_for_aa
                           -> store_string (in common.c)
                           -> binary_hetdrug_search
                           -> to_lower (in common.c)
                           -> sort_ligands
                           -> find_ligand
                     -> open_output_file (in common.c)
                     -> check_other_chains
                           -> add_name
                     -> perform_splice (in common.c)
                     -> write_pfam_details
                           -> order_domains (as above)
                     -> write_domain_details
               -> free_chain
                     -> free_domain
               -> free_clist
         -> binary_search
         -> open_output_file
   -> write_ligpfam_analyses
         -> sort_pfamlig_entries
               -> pfamlig_sort1
               -> pfamlig_sort2
         -> add_old_pfamlig_data
               -> create_pdb_index
                     -> pdb_sort
               -> create_pfam_index
                     -> pfam_sort
               -> binary_pdb_search
               -> binary_pfam_search
               -> create_plist_entry
               -> create_pfamlig_record
         -> write_lig2pfam_dat
         -> write_index
         -> write_pfam2lig_dat


*/

#include "common.h"
#include "readpdb.h"
#include "uniplot.h"

/* Prototypes
   ---------- */

/* In common.c */
void append_string(char **string_ptr,char *string);
void call_system(char *command,int dos);
int check_file(char *file_name);
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int empty_check(char *file_name);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
void free_tokens(char **token,int ntokens);
int get_c_scores(char *pdb_code,char *pdbsum_dir,
                 struct c_score **fst_c_score_ptr);
int get_chain_colour(char chain,double *rgb,double *rgb_light);
void get_colour_rgb(char *c_name,double rgb[3]);
int get_tokens(char *line,char **token,int max_token,char token_sep);
FILE *open_output_file(char *out_name,int terminate);
int perform_splice(char *line,char *string,char *new_string);
int run_wget(char *wget_exe,char *url,char *out_file,int dos);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);
int tokenize(char *string,char **token,int max_tokens,char token_sep);
void to_lower(char *search_string,int length);
void to_upper(char *search_string,int length);

/* In align.c */
float align_sequences(char *sequence1,char *sequence2,
                      char *alignment[2],int *aln_len,int *nident);

/* In readpdb.c */
int r_get_dfn_data(char *pdb_code,int pdb_type,char *pdbsum_template,
                   struct r_rasdef **fst_r_rasdef_ptr,
                   struct r_ligand **fst_r_ligand_ptr,int *nlig,
                   int *nmet);
int r_list_rasdef_entries(struct r_rasdef *first_r_rasdef_ptr);

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

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->pdb_type = STANDARD_PDB;
  params->pdbsum_type = PDBSUM_CURRENT;
  params->run_type = STANDARD;
  params->covid19 = FALSE;
  params->wget = FALSE;
  params->run_64 = FALSE;

  /* If -halp option entered, show options available */
  if (num == 1 && (!strcmp(strptrs[1],"-help") ||
                   !strcmp(strptrs[1],"-h")))
    {
      printf("\n");
      printf("Usage is ...\n"
             "\n"
             "  uniplot  [UPLOAD | MCSG | AF]  [-weekly]  [-covid19]  "
             "[-pdbsum1 pdb_code]  [-wget]\n"
             "\n"
             "where\n"
             "   * UPLOAD            = Uploaded PDB file for PDBsum "
             "generation\n"
             "   * MCSG              = MCSG structure\n"
             "     * AF              = Alpha Fold structure\n"
             "   * -weekly           = Indicates weekly update run\n"
             "   * -covid19          = processing COVID-19 structures\n"
             "   * -pdbsum1 pdb_code = PDBsum1 structure\n"
             "   * -wget             = Indicates Pfam data to be retrieved "
             "using wget\n"
             "   * -64               = run on 64-bit processor\n"
             "\n"
             "\n");
      exit(1);
    }

  /* Loop through the command arguments */
  for (itoken = 1; itoken < num + 1; itoken++)
    {
      /* Check for the UPLOAD option */
      if (!strcmp(strptrs[itoken],"UPLOAD"))
        {
          params->pdb_type = UPLOAD_PDB;
          params->pdbsum_type = PDBSUM_UPLOAD;
        }

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        {
          params->pdb_type = MCSG_PDB;
          params->pdbsum_type = PDBSUM_UPLOAD;
        }

      /* Check for the Alpha Fold option */
      else if (!strcmp(strptrs[itoken],"AF"))
        params->pdb_type = ALPHA_FOLD_PDB;

      /* Check for the weekly update flag */
      else if (!strcmp(strptrs[itoken],"-weekly"))
        params->run_type = WEEKLY_UPDATE;

      /* Check for the covid19 option */
      else if (!strcmp(strptrs[itoken],"-covid19"))
        params->covid19 = TRUE;

      /* Check for the wget option */
      else if (!strcmp(strptrs[itoken],"-wget"))
        params->wget = TRUE;

      /* Check for the PDBsum1 option */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        {
          /* Get the PDB code from the next token */
          if (itoken < num)
            {
              /* Get the PDB code */
              strcpy(params->pdb_code,strptrs[itoken + 1]);

              /* Increment token count */
              itoken++;
            }

          /* Set flag */
          params->pdb_type = PDBSUM1;
        }

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
/***********************************************************************

get_letters  -  Pick up the character definitions from the
                       raslet.dat file

***********************************************************************/

int get_letters(char *character,struct letter *letter_ptr,
                struct c_file_data file_name_ptr)
{
  char ch, file_name[FILENAME_LEN], input_line[LINELEN + 1];

  int icol, iletter, ipos, irow, len, line, ncols, nletters;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  ncols = 0;
  nletters = 0;
  for (iletter = 0; iletter < NLETTERS + 1; iletter++)
    character[iletter] = '\0';
  iletter = -1;

  /* Form full name of letters.dat file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"letters.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file processing each record */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment line-count */
          line++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Check whether this is the start of a new character */
          if (len > 1 && input_line[0] == '#')
            {
              /* Get the character */
              ch = input_line[1];

              /* If cannot store this letter end here */
              if (nletters > NLETTERS - 1)
                {
                  printf("*** ERROR. Arrays not large enough to store "
                         "letters: %d\n",NLETTERS);
                  exit(1);
                }

              /* Increment count of letters */
              iletter = nletters;
              nletters++;

              /* Store this character in the array */
              character[iletter] = ch;

              /* Initialise row counter */
              irow = 0;
            }

          /* Otherwise, store the current row of the character definition */
          else if (irow < LETTER_ROWS)
            {
              /* If this is the first row, determine its length and
                 grab memory for array */
              if (irow == 0)
                {
                  /* Get number of characters in this row */
                  ncols = len;

                  /* Allocate memory to store array */
                  letter_ptr[iletter].array
                    = (char *) malloc(ncols * LETTER_ROWS * sizeof(char));
                  if (letter_ptr[iletter].array == NULL)
                    {
                      printf("*** ERROR. Unable to allocate memory for "
                             "letter array\n");
                      exit(1);
                    }

                  /* Store letter's dimensions */
                  letter_ptr[iletter].width = len;
                  letter_ptr[iletter].height = LETTER_ROWS;
                }

              /* Store current row */
              for (icol = 0; icol < ncols; icol++)
                {
                  /* Get the location in the array */
                  ipos = irow * ncols + icol;

                  /* Get current character, replacing dot by space */
                  ch = input_line[icol];
                  if (ch == '.')
                    ch = ' ';

                  /* Copy current character to array */
                  letter_ptr[iletter].array[ipos] = ch;
                }

              /* Increment row count */
              irow++;
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show warning that couldn't find the label file */
  else
    {
      printf("*** Warning. Unable to find letters.dat file [%s]\n",
             file_name);
    }

  /* Return the number of characters read in */
  return(nletters);
}
/***********************************************************************

create_hetdrug_entry  -  Create a Het Group to drug-id entry

***********************************************************************/

struct hetdrug *create_hetdrug_entry(struct hetdrug **fst_hetdrug_ptr,
                                     struct hetdrug **lst_hetdrug_ptr,
                                     char *het_name,char *drug_id,
                                     char *generic_name,char *drug_desc)
{
  struct hetdrug *hetdrug_ptr;
  struct hetdrug *first_hetdrug_ptr, *last_hetdrug_ptr;

  /* Initialise hetdrug pointers */
  hetdrug_ptr = NULL;
  first_hetdrug_ptr = *fst_hetdrug_ptr;
  last_hetdrug_ptr = *lst_hetdrug_ptr;

  /* Create hetdrug entry */
  hetdrug_ptr = (struct hetdrug*)malloc(sizeof(struct hetdrug));
  if (hetdrug_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct hetdrug\n");
      printf("***        Program terminated with");
      printf(" error.\n");
      exit (1);
    }

  /* Store the UniProt accession code and nunmber of sequences */
  strcpy(hetdrug_ptr->het_name,het_name);
  store_string(&(hetdrug_ptr->drug_id),drug_id);
  store_string(&(hetdrug_ptr->generic_name),generic_name);
  store_string(&(hetdrug_ptr->drug_desc),drug_desc);

  /* Initialise pointer to next entry */
  hetdrug_ptr->next_hetdrug_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_hetdrug_ptr == NULL)
    first_hetdrug_ptr = hetdrug_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_hetdrug_ptr->next_hetdrug_ptr = hetdrug_ptr;
                      
  /* Save the current residue's pointer */
  last_hetdrug_ptr = hetdrug_ptr;

  /* Return current pointers */
  *fst_hetdrug_ptr = first_hetdrug_ptr;
  *lst_hetdrug_ptr = last_hetdrug_ptr;
  return(hetdrug_ptr);
}
/***********************************************************************

get_hetdrug_mappings  -  Read in the Het Groups that correspond to drug
                         molecules

***********************************************************************/

int get_hetdrug_mappings(struct hetdrug **fst_hetdrug_ptr,
                         struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char dir_name[3], journal[LINELEN + 1], search_name[6];
  char drug_id[20], drug_desc[30], generic_name[LINELEN], het_name[4];

  char *string_ptr;

  int i, len, nhetdrug;
  int wanted;

  struct hetdrug *hetdrug_ptr, *first_hetdrug_ptr, *last_hetdrug_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  drug_id[0] = drug_desc[0] = generic_name[0] = '\0';
  nhetdrug = 0;

  /* Initialise hetdrug pointers */
  hetdrug_ptr = first_hetdrug_ptr = last_hetdrug_ptr = NULL;

  /* Form full name of the appdrugs_pdb.dat file */
  strcpy(file_name,file_name_ptr.other_dir[DRUGPORT_REAL_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"appdrugs_pdb.dat");
  printf("Reading %s ...\n",file_name);
  fflush(stdout);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file, line by line */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the input line */
          len = string_truncate(input_line,LINELEN);

          /* Separate the field name from the field value */
          string_ptr = strstr(input_line," ");
          if (string_ptr == NULL)
            {
              strcat(input_line," ");
              string_ptr = strstr(input_line," ");
            }

          /* If this is the generic name, save it */
          if (!strncmp(input_line,"GENERIC_NAME[",13))
            strcpy(generic_name,string_ptr+1);

          /* If this is the DrugPort identifier, save it */
          else if (!strncmp(input_line,"DRUGNAME_ID[",12))
            strcpy(drug_id,string_ptr+1);

          /* If this is the drug description, save it */
          else if (!strncmp(input_line,"DRUG_DESC[",10))
            strcpy(drug_desc,string_ptr+1);

          /* If this is the Het Group name, create a new hetdrug record */
          else if (!strncmp(input_line,"DRUG_HET_NAME[",14) &&
                   drug_id[0] != '\0')
            {
              /* Get the Het Group name */
              strcpy(het_name,string_ptr+1);

              /* Assume it is wanted */
              wanted = TRUE;

              /* Check it is not one of the small molecules that aren't
                 really drugs */
              sprintf(search_name,"|%s|",het_name);
              if (strstr(NOT_DRUGS_LIST,search_name) != NULL)
                wanted = FALSE;

              /* If wanted, create a hetdrug record */
              if (wanted == TRUE)
                {
                  /* Create new hetdrug record to store the mapping */
                  hetdrug_ptr
                    = create_hetdrug_entry(&first_hetdrug_ptr,
                                           &last_hetdrug_ptr,het_name,
                                           drug_id,generic_name,drug_desc);
                  
                  /* Increment count of records created */
                  nhetdrug++;
                }

              /* Re-initialise drug id and generic name */
              drug_id[0] = generic_name[0] = '\0';
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return current pointer */
  *fst_hetdrug_ptr = first_hetdrug_ptr;

  /* Return number of hetdrug records created */
  printf("Number of Het Group -> drug mappings = %8d\n",nhetdrug);
  return(nhetdrug);
}
/***********************************************************************

create_hetdrug_index  -  Create the array holding the sort-indices for the
                      hetdrug data

***********************************************************************/

void create_hetdrug_index(struct hetdrug *first_hetdrug_ptr,
                          struct hetdrug ***hetdrug_ind_ptr,
                          int nhetdrug)
{
  int nentry;

  struct hetdrug *hetdrug_ptr;

  struct hetdrug **hetdrug_index_ptr;

  /* Initialise variables */
  nentry = 0;

  /* Create the index array as one large array */
  hetdrug_index_ptr
    = (struct hetdrug **) malloc(nhetdrug * sizeof(struct hetdrug *));

  /* Check that sufficient memory was allocated */
  if (hetdrug_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for array\n");
      exit(1);
    }

  /* Get pointer to the first hetdrug entry */
  hetdrug_ptr = first_hetdrug_ptr;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (hetdrug_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      hetdrug_index_ptr[nentry] = hetdrug_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next record */
      hetdrug_ptr = hetdrug_ptr->next_hetdrug_ptr;
    }

  /* Return the array pointer */
  *hetdrug_ind_ptr = hetdrug_index_ptr;
}
/***********************************************************************

sort_hetdrug_entries  -  Sort the Het Group hetdrug entries

***********************************************************************/

void sort_hetdrug_entries(struct hetdrug **hetdrug_index_ptr,int nhetdrug)
{
  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;

  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;

  double dble, dble1;

  struct hetdrug *hetdrug_ptr, *other_hetdrug_ptr, *swap_hetdrug_ptr;

  /*--Initialise values */
  n = nhetdrug;
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

              /* Get the two NIST records given by the current indices */
              hetdrug_ptr = hetdrug_index_ptr[indi];
              other_hetdrug_ptr = hetdrug_index_ptr[indl];

              /* If names in wrong order, then swap sort-pointers */
              if (strcmp(hetdrug_ptr->het_name,
                         other_hetdrug_ptr->het_name) > 0)
                {
                  swap_hetdrug_ptr = hetdrug_index_ptr[indi];
                  hetdrug_index_ptr[indi] = hetdrug_index_ptr[indl];
                  hetdrug_index_ptr[indl] = swap_hetdrug_ptr;

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

create_pdb_entry  -  Create a pdb entry

***********************************************************************/

struct pdb *create_pdb_entry(struct pdb **fst_pdb_ptr,
                             struct pdb **lst_pdb_ptr,char *pdb_code)
{
  struct pdb *pdb_ptr;
  struct pdb *first_pdb_ptr, *last_pdb_ptr;

  /* Initialise pdb pointers */
  pdb_ptr = NULL;
  first_pdb_ptr = *fst_pdb_ptr;
  last_pdb_ptr = *lst_pdb_ptr;

  /* Create pdb entry */
  pdb_ptr = (struct pdb *) malloc(sizeof(struct pdb));
  if (pdb_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct pdb\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the PDB code */
  strcpy(pdb_ptr->pdb_code,pdb_code);

  /* Initialise remaining fields */
  pdb_ptr->nchains = 0;
  pdb_ptr->nunique = 0;
  pdb_ptr->have_cath = FALSE;
  pdb_ptr->have_pfam = FALSE;
  pdb_ptr->coloured = FALSE;
  pdb_ptr->done = FALSE;
  pdb_ptr->first_chain_ptr = NULL;

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

get_pdb_codes  -  Read through the dataset.lst file to read in the PDB
                  codes to be processed

***********************************************************************/

int get_pdb_codes(char *dataset_name,struct pdb **fst_pdb_ptr,
                  char *single_pdb_code)
{
  char pdb_code[5];
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];

  int len, ncodes;

  struct pdb *pdb_ptr;
  struct pdb *first_pdb_ptr, *last_pdb_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  pdb_ptr = NULL;
  first_pdb_ptr = last_pdb_ptr = NULL;
  ncodes = 0;

  /* If have a single PDB code, then save that and return */
  if (single_pdb_code[0] != '\0')
    {
      /* Create a new pdb entry */
      pdb_ptr = create_pdb_entry(&first_pdb_ptr,&last_pdb_ptr,
                                 single_pdb_code);

      /* Increment count */
      ncodes++;

      /* Return pointer to first pdb entry */
      *fst_pdb_ptr = first_pdb_ptr;

      /* Return the number of sequences read in */
      return(ncodes);
    }

  /* Form name of dataset.lst file */
  strcpy(file_name,dataset_name);

  /* Open the dataset.lst file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** ERROR. Unable to open file [%s]\n",file_name);
      exit(1);
    }

  /* Loop through the file, line by line */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_truncate(input_line,LINELEN);

      /* Check line-length */
      if (len >= 4)
        {
          /* Get the PDB code */
          strncpy(pdb_code,input_line,4);
          pdb_code[4] = '\0';

          /* Create a new pdb entry */
          pdb_ptr = create_pdb_entry(&first_pdb_ptr,&last_pdb_ptr,
                                     pdb_code);

          /* Increment count */
          ncodes++;
        }
    }

  /* Close file */
  fclose(file_ptr);

  /* Return pointer to first pdb entry */
  *fst_pdb_ptr = first_pdb_ptr;

  /* Return the number of sequences read in */
  printf("Number of PDB codes read in:       %8d\n",ncodes);
  return(ncodes);
}
/***********************************************************************

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  struct pdb *pdb_ptr,char chain_id,
                                  char copy_of,int chain_no)
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
  chain_ptr->chain_no = chain_no;
  chain_ptr->chain_id = chain_id;
  chain_ptr->copy_of = copy_of;
  chain_ptr->pdb_ptr = pdb_ptr;
  chain_ptr->next_chain_ptr = NULL;

  /* Initialise remaining fields */
  chain_ptr->ncodes = 0;
  chain_ptr->nresid = 0;
  chain_ptr->image_height = 0;
  chain_ptr->image_width = 0;
  chain_ptr->last_width = 0;
  chain_ptr->ncath_dom = 0;
  chain_ptr->nligs = 0;
  chain_ptr->nligands = 0;
  for (i = 0; i < MAX_CODES; i++)
    {
      chain_ptr->nimages[i] = 0;
      chain_ptr->code_ptr[i] = NULL;
      chain_ptr->plot_no[i] = 0;
      chain_ptr->smatch_ptr[i] = NULL;
    }
  chain_ptr->align_len = 0;
  chain_ptr->nidentical = 0;
  chain_ptr->identity = 0;
  chain_ptr->alignment[0]= chain_ptr->alignment[1] = &null;
  chain_ptr->first_cath_domain_ptr = NULL;
  chain_ptr->first_ligand_ptr = NULL;
  chain_ptr->last_ligand_ptr = NULL;
  chain_ptr->first_liglist_ptr = NULL;
  chain_ptr->last_liglist_ptr = NULL;

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return new chain pointer */
  return(chain_ptr);
}
/***********************************************************************

get_chains_dat  -  Read in the chain equivalences from the chains.dat
                   file

***********************************************************************/

void get_chains_dat(struct pdb *pdb_ptr,char *pdbsum_dir)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char chain_id, copy_of;

  int len, nlines;

  struct chain *chain_ptr, *last_chain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nlines = 0;

  /* Initialise chain pointers */
  chain_ptr = NULL;
  last_chain_ptr = NULL;

  /* Form the name of the chains.dat file in the PDBsum directory */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"chains.dat");

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment count of lines read in */
          nlines++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

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
              pdb_ptr->nunique++;
            }

          /* Create a chain record */
          chain_ptr
            = create_chain_record(&(pdb_ptr->first_chain_ptr),
                                  &last_chain_ptr,pdb_ptr,chain_id,
                                  copy_of,pdb_ptr->nunique);

          /* Increment count of this PDB entry's protein chains */
          pdb_ptr->nchains++;
        }

      /* Close the input file */
      fclose(file_ptr);
    }
}
/***********************************************************************

find_code_record  -  Find the given code entry

***********************************************************************/

struct code *find_code_record(struct code *first_code_ptr,
                              char *uniprot_acc)
{
  struct code *code_ptr, *found_code_ptr;

  /* Initialise variables */
  found_code_ptr = NULL;
  if (uniprot_acc[0] == '\0')
    return(found_code_ptr);

  /* Get pointer to the first code */
  code_ptr = first_code_ptr;

  /* Loop over all the codes to find the one just read in */
  while (code_ptr != NULL && found_code_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (!strcmp(code_ptr->uniprot_acc,uniprot_acc))
        found_code_ptr = code_ptr;

      /* Get pointer to the next code */
      code_ptr = code_ptr->next_code_ptr;
    }

  /* Return the matching code pointer */
  return(found_code_ptr);
}
/***********************************************************************

create_code_record  -  Create a new code record

***********************************************************************/

struct code *create_code_record(struct code **fst_code_ptr,
                                struct code **lst_code_ptr,
                                char *uniprot_acc)
{
  struct code *code_ptr;
  struct code *first_code_ptr, *last_code_ptr;

  /* Initialise code pointers */
  code_ptr = NULL;
  first_code_ptr = *fst_code_ptr;
  last_code_ptr = *lst_code_ptr;

  /* Create code entry */
  code_ptr = (struct code*)malloc(sizeof(struct code));
  if (code_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct code\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  strcpy(code_ptr->uniprot_acc,uniprot_acc);

  /* Initialise other fields */
  code_ptr->uniprot_id[0] = '\0';
  code_ptr->name = &null;
  code_ptr->seq_len = 0;
  code_ptr->pfam_seq_len = 0;
  code_ptr->have_seq = FALSE;
  code_ptr->done = FALSE;
  code_ptr->covid_protein = -1;
  code_ptr->first_clist_ptr = NULL;
  code_ptr->last_clist_ptr = NULL;
  code_ptr->first_domain_ptr = NULL;

  /* Initialise pointer to next code record */
  code_ptr->next_code_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of codeed code */
  if (first_code_ptr == NULL)
    first_code_ptr = code_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_code_ptr->next_code_ptr = code_ptr;
                      
  /* Save the current code's pointer */
  last_code_ptr = code_ptr;

  /* Return current pointers */
  *fst_code_ptr = first_code_ptr;
  *lst_code_ptr = last_code_ptr;
  return(code_ptr);
}
/***********************************************************************

create_clist_entry  -  Create a clist entry

***********************************************************************/

struct clist *create_clist_entry(struct clist **fst_clist_ptr,
                                 struct clist **lst_clist_ptr,
                                 struct chain *chain_ptr,int icode)
{
  struct clist *clist_ptr;
  struct clist *first_clist_ptr, *last_clist_ptr;

  /* Initialise clist pointers */
  clist_ptr = NULL;
  first_clist_ptr = *fst_clist_ptr;
  last_clist_ptr = *lst_clist_ptr;

  /* Create clist entry */
  clist_ptr = (struct clist *) malloc(sizeof(struct clist));
  if (clist_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct clist\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  clist_ptr->chain_ptr = chain_ptr;
  clist_ptr->icode = icode;

  /* Initialise pointer to next entry */
  clist_ptr->next_clist_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_clist_ptr == NULL)
    first_clist_ptr = clist_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_clist_ptr->next_clist_ptr = clist_ptr;
                      
  /* Save the current residue's pointer */
  last_clist_ptr = clist_ptr;

  /* Return current pointers */
  *fst_clist_ptr = first_clist_ptr;
  *lst_clist_ptr = last_clist_ptr;
  return(clist_ptr);
}
/***********************************************************************

create_smatch_record  -  Create and initialise a new smatch record

***********************************************************************/

struct smatch *create_smatch_record(struct code *code_ptr)
{
  int i;

  struct smatch *smatch_ptr;

  /* Initialise smatch pointers */
  smatch_ptr = NULL;

  /* Allocate memory for structure to hold smatch info */
  smatch_ptr = (struct smatch *) malloc(sizeof(struct smatch));
  if (smatch_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct smatch\n");
      exit (1);
    }

  /* Store known data */
  smatch_ptr->code_ptr = code_ptr;

  /* Initialise remaining fields */
  smatch_ptr->align_len = 0;
  smatch_ptr->nmismatch = 0;
  smatch_ptr->uniprot_seq = &null;
  smatch_ptr->pdb_seq = &null;
  smatch_ptr->pdb_hyphens = &null;
  smatch_ptr->pdb_resno = &null;
  smatch_ptr->diffs = &null;
  smatch_ptr->first_domain_ptr = NULL;

  /* Return new smatch pointer */
  return(smatch_ptr);
}
/***********************************************************************

get_seqinfo_dat  -  Read in all the references, from the pnames.dat file

***********************************************************************/

int get_seqinfo_dat(struct pdb *pdb_ptr,struct code **fst_code_ptr,
                    struct code **lst_code_ptr,int *ncod,
                    char *pdbsum_dir,struct code *dummy_code_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char number_string[20], uniprot_acc[UNIPROT_LEN];

  char *string_ptr;

  int i, icode, ichain, iunique, len, ncodes, plot_no;
  int found, have_info;

  struct chain *chain_ptr;
  struct code *code_ptr;
  struct code *first_code_ptr, *last_code_ptr;
  struct smatch *smatch_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  have_info = FALSE;
  ncodes = *ncod;
  plot_no = 0;

  /* Initialise pointers */
  code_ptr = NULL;
  first_code_ptr = *fst_code_ptr;
  last_code_ptr = *lst_code_ptr;

  /* Form the name of the input seqinfo.dat file in the PDBsum
     directory */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"seqinfo.dat");

  /* Open the file to see if it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop through the file, line by line */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the input line */
          len = string_truncate(input_line,LINELEN);

          /* Separate the field name from the field value */
          string_ptr = strstr(input_line," ");
          if (string_ptr == NULL)
            {
              strcat(input_line," ");
              string_ptr = strstr(input_line," ");
            }

          /* If this is the start of a new reference, then create a
             new record to store details */
          if (!strncmp(input_line,"SWISSPROT_ID",12) &&
              strchr(input_line,'-') == NULL)
            {
              /* Get the UniProt identifier */
              strcpy(uniprot_acc,string_ptr+1);

              /* See if we already have this UniProt code */
              code_ptr = find_code_record(first_code_ptr,uniprot_acc);

              /* If don't yet have this record, create a new one */
              if (code_ptr == NULL)
                {
                  /* Create a UniProt code record */
                  code_ptr
                    = create_code_record(&first_code_ptr,&last_code_ptr,
                                         uniprot_acc);

                  /* Increment count of code records */
                  ncodes++;
                }

              /* Get which unique chain this applies to */
              plot_no++;
              ichain = 0;
              string_ptr = strstr(input_line,"[");
              if (string_ptr != NULL)
                {
                  strcpy(number_string,string_ptr+1);
                  string_ptr = strchr(number_string,']');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Get the chain number */
                  ichain = atoi(number_string);
                }

              /* If have a valid chain assignment, store the associated
                 UniProt code */
              if (ichain > 0)
                {
                  /* Initialise flag */
                  found = FALSE;
                  iunique = 0;

                  /* Get the first of this PDB entry's chains */
                  chain_ptr = pdb_ptr->first_chain_ptr;

                  /* Loop through all the chain records until find
                     the correct one */
                  while (chain_ptr != NULL && found == FALSE)
                    {
                      /* If this is a unique chain, increment unique
                         chains count */
                      if (chain_ptr->copy_of == '\0')
                        iunique++;
                      
                      /* If this is the chain for which we have just
                         picked up the UniProt code, then store */
                      if (iunique == ichain)
                        {
                          /* Initialise */
                          i = 0;

                          /* If not too many codes stored already,
                             store this one */
                          if (chain_ptr->ncodes < MAX_CODES)
                            {
                              /* Store the current code */
                              i = chain_ptr->ncodes;
                              chain_ptr->code_ptr[i] = code_ptr;

                              /* Increment code count */
                              chain_ptr->ncodes++;

                              /* Store count of chains with UniProt ids */
                              chain_ptr->plot_no[i] = plot_no;

                              /* Create new clist record to link
                                 current UniProt code to the current
                                 chain */
                              create_clist_entry(&(code_ptr->first_clist_ptr),
                                                 &(code_ptr->last_clist_ptr),
                                                 chain_ptr,i);

                              /* Create a new smatch record to store the
                                 details of this UniProt code */
                              smatch_ptr
                                = create_smatch_record(code_ptr);

                              /* Store for this chain */
                              chain_ptr->smatch_ptr[i] = smatch_ptr;
                            }

                          /* Set flag that we're done */
                          found = TRUE;

                        }

                      /* Get pointer to the next chain record */
                      chain_ptr = chain_ptr->next_chain_ptr;
                    }
                }

              /* Set flag that we have an associated UniProt id */
              have_info = TRUE;

            }

          /* If this is the UniProt id, save it */
          else if (!strncmp(input_line,"SWISSPROT_CODE",14) &&
                   code_ptr != NULL)
            strcpy(code_ptr->uniprot_id,string_ptr+1);
        }

      /* Close it */
      fclose(file_ptr);
    }

  /* Re-initialise plot number */
  plot_no = 0;

  /* Loop over this PDB entry's unique chains to link any without
     UniProt codes to the dummy code */
  chain_ptr = pdb_ptr->first_chain_ptr;

  /* Loop through all the chain records to check which don't have a
     UniProt code */
  while (chain_ptr != NULL)
    {
      /* If this is a unique chain with no UniProt code, then assign
         to the blank */
      if (chain_ptr->copy_of == '\0' && chain_ptr->ncodes == 0)
        {
          /* Add the dummy code to the current chain */
          chain_ptr->code_ptr[0] = dummy_code_ptr;

          /* Increment code count */
          chain_ptr->ncodes++;

          /* Set sequence counter from seqinfo.dat to zero */
          chain_ptr->plot_no[0] = 0;

          /* Create new clist record to link dummy UniProt code to
             the current chain */
          create_clist_entry(&(dummy_code_ptr->first_clist_ptr),
                             &(dummy_code_ptr->last_clist_ptr),
                             chain_ptr,0);

          /* Create a new smatch record to store the
             details of this UniProt code */
          smatch_ptr = create_smatch_record(dummy_code_ptr);

          /* Store for this chain */
          chain_ptr->smatch_ptr[0] = smatch_ptr;
        }

      /* If this is a unique chain with a UniProt code, then reset the
         plot number */
      if (chain_ptr->copy_of == '\0')
        {
          /* Loop over this chain's UniProt codes */
          for (icode = 0; icode < chain_ptr->ncodes; icode++)
            {
              /* Increment plot number */
              plot_no++;

              /* Save this plot number */
              chain_ptr->plot_no[icode] = plot_no;
            }
        }

      /* Get pointer to the next chain record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return current pointers */
  *fst_code_ptr = first_code_ptr;
  *lst_code_ptr = last_code_ptr;

  /* Return current number of unique UniProt codes */
  *ncod = ncodes;

  /* Return whether any UniProt assignments found */
  return(have_info);
}
/***********************************************************************

process_codes  -  Process all the codes in turn

***********************************************************************/

int process_codes(struct pdb *first_pdb_ptr,struct code **fst_code_ptr,
                  struct code **lst_code_ptr,int npdb,
                  struct c_file_data file_name_ptr,
                  struct parameters params)
{
  char pdbsum_dir[FILENAME_LEN + 1];

  int dir_type, ncodes, nread;
  int have_info;

  struct code *code_ptr, *first_code_ptr, *last_code_ptr;
  struct pdb *pdb_ptr;

  FILE *auditfile_ptr;

  /* Initialise variables */
  ncodes = 0;
  nread = 0;
  first_code_ptr = NULL;
  last_code_ptr = NULL;

  /* Create a dummy UniProt record for all protein chains that have no
     UniProt assignment */
  code_ptr = create_code_record(&first_code_ptr,&last_code_ptr,"NONE");
  ncodes++;

  /* Get the first of the PDB entries */
  pdb_ptr = first_pdb_ptr;

  /* Loop through the PDB entries */
  while (pdb_ptr != NULL)
    {
      /* Form the name of the PDBsum directory */
      if (params.pdb_type == PDBSUM1)
        {
          strcpy(pdbsum_dir,"../newsec/");
          dir_type = PDBSUM1;
        }
      else
        dir_type = find_pdbsum_dir(pdb_ptr->pdb_code,
                                   file_name_ptr.pdbsum_template,
                                   params.pdbsum_type,pdbsum_dir);

      /* If have PDBsum directory, then continue */
      if (dir_type != NOT_FOUND)
        {
          /* Pick up the chains from the chains.dat file, if present */
          get_chains_dat(pdb_ptr,pdbsum_dir);

          /* If PDB file contains any protein chains, then process */
          if (pdb_ptr->nchains > 0)
            {
              /* Form the name of the PDBsum directory */
              if (params.pdb_type == PDBSUM1)
                strcpy(pdbsum_dir,"../seqdata/");

              /* Pick up the UniProt code corresponding to each chain
                 from the seqinfo.dat file, if present */
              have_info
                = get_seqinfo_dat(pdb_ptr,&first_code_ptr,&last_code_ptr,
                                  &ncodes,pdbsum_dir,code_ptr);
            }
        }

      /* Display error as PDBsum directory not found */
      else
        printf("*** Warning. Unable to locate PDBsum directory for %s\n",
               pdb_ptr->pdb_code);

      /* Increment count */
      nread++;

      /* Write progress to audit file */
      if ((nread % NAUDIT_PDB) == 0)
        {
          auditfile_ptr = open_output_file("audit.out",TRUE);
          fprintf(auditfile_ptr,"Processing PDB files. %d of %d done ...\n",
                  nread,npdb);
          fclose(auditfile_ptr);
        }

      /* Get the next of the PDB entries */
      pdb_ptr = pdb_ptr->next_pdb_ptr;
    }

  /* Return current pointers */
  *fst_code_ptr = first_code_ptr;
  *lst_code_ptr = last_code_ptr;

  /* Return number of unique UniProt ids stored */
  return(ncodes);
}
/***********************************************************************

read_old_codes  -  Pick up the existing data in the lig2pfam.dat
                   file and create a record for any data we do not
                   have

***********************************************************************/

int read_old_codes(struct code **fst_code_ptr,struct code *last_code_ptr,
                   int ncodes,struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char uniprot_acc[UNIPROT_LEN];

  char *string_ptr;

  int len, line;

  struct code *code_ptr, *first_code_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  uniprot_acc[0] = '\0';

  /* Initialise pointers */
  first_code_ptr = *fst_code_ptr;

  /* Form name of input file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"lig2pfam.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. File not found: %s\n",file_name);
      return(ncodes);
    }

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* If this is the UniProt accession code, retrieve it */
      if (!strncmp(input_line,"L2P_UNIPROT_ACC",15))
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);
          line++;

          /* Separate the field name from the field value */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            {
              /* Get the UniProt accession code */
              strcpy(uniprot_acc,string_ptr+1);

              /* See if we already have this UniProt code */
              code_ptr = find_code_record(first_code_ptr,uniprot_acc);

              /* If don't yet have this record, create a new one */
              if (code_ptr == NULL)
                {
                  /* Create a UniProt code record */
                  code_ptr
                    = create_code_record(&first_code_ptr,&last_code_ptr,
                                         uniprot_acc);

                  /* Increment count of code records */
                  ncodes++;
                }
            }
        }

      /* If this is the UniProt id, save it */
      else if (!strncmp(input_line,"L2P_UNIPROT_ID",14) &&
               code_ptr != NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Separate the field name from the field value */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            {
              /* Store the UniProt id */
              strcpy(code_ptr->uniprot_id,string_ptr+1);
            }
        }
    }

  /* Close file */
  fclose(file_ptr);

  /* Return first code pointer */
  *fst_code_ptr = first_code_ptr;

  /* Return number of unique UniProt codes */
  return(ncodes);
}
/***********************************************************************

code_sort  -  Shell-sort which sorts the UniProt codes into ascending
              order

***********************************************************************/

void code_sort(struct code **code_index_ptr,int nindex)
{
  char code_name1[UNIPROT_LEN], code_name2[UNIPROT_LEN];

  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct code *code1_ptr, *code2_ptr, *swap_code_ptr;

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

              /* Get the two code records */
              code1_ptr = code_index_ptr[indi];
              code2_ptr = code_index_ptr[indl];

              /* Form full code names on which sort is based */
              strcpy(code_name1,code1_ptr->uniprot_acc);
              strcpy(code_name2,code2_ptr->uniprot_acc);

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(code_name1,code_name2) > 0)
                {
                  swap_code_ptr = code_index_ptr[indi];
                  code_index_ptr[indi] = code_index_ptr[indl];
                  code_index_ptr[indl] = swap_code_ptr;

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

create_code_index  -  Create index of UniProt records for use in binary
                      search

***********************************************************************/

int create_code_index(struct code *first_code_ptr,
                      struct code ***code_index_ptr,int ncodes)
{
  int ncovid, noncovid, nentry;
  int entry_types;

  struct code *code_ptr;

  struct code **code_ind_ptr;

  /* Initialise variables */
  entry_types = NO_COVID;
  ncovid = noncovid = 0;
  nentry = 0;

  /* Create the index array as one large array */
  code_ind_ptr = (struct code **) malloc(ncodes * sizeof(struct code *));

  /* Check that sufficient memory was allocated */
  if (code_ind_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Get pointer to the first code entry */
  code_ptr = first_code_ptr;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (code_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      if (nentry < ncodes)
        code_ind_ptr[nentry] = code_ptr;

      /* If this is a COVID-19 UniProt code, set flag */
      if (strchr(code_ptr->uniprot_acc,'_') != NULL ||
          strstr(code_ptr->uniprot_id,"_SARS2") != NULL)
        ncovid++;
      else if (strcmp(code_ptr->uniprot_acc,"NONE"))
        noncovid++;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      code_ptr = code_ptr->next_code_ptr;
    }

  /* Sort the code pointers into order of UniProt code */
  if (ncodes > 1)
    code_sort(code_ind_ptr,ncodes);

  /* Return the array pointer */
  *code_index_ptr = code_ind_ptr;

  /* Determine entry types */
  if (ncovid > 0 && noncovid > 0)
    entry_types = MIXED;
  else if (ncovid > 0)
    entry_types = COVID_ONLY;
  else
    entry_types = NO_COVID;

  /* Return flag indicating entry types */
  return(entry_types);
}
/***********************************************************************

binary_search  -  Search for UniProt record

***********************************************************************/

struct code *binary_search(struct code **code_index_ptr,int ncodes,
                           char *uniprot_acc)
{
  int bottom_point, ipoint, new_point, top_point;
  int found;

  struct code *code_ptr, *found_code_ptr, *prev_code_ptr;

  /* Initialise variables */
  found = FALSE;
  found_code_ptr = NULL;

  /* Make starting point half-way into the index */
  ipoint = ncodes / 2;
  bottom_point = 0;
  top_point = ncodes - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated GO pointer */
      code_ptr = code_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(uniprot_acc,code_ptr->uniprot_acc))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_code_ptr = code_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(uniprot_acc,code_ptr->uniprot_acc) < 0)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_code_ptr = code_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_code_ptr->uniprot_acc,uniprot_acc) < 0)
                found = TRUE;

              /* Otherwise, need to search earlier in the list */
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

      /* Otherwise, have to look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= ncodes - 1)
            {
              /* Set flag as done */
              found = TRUE;
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

  /* Return the pointer */
  return(found_code_ptr);
}
/***********************************************************************

find_pfam_entry  -  Find the given pfam

***********************************************************************/

struct pfam *find_pfam_entry(struct pfam *first_pfam_ptr,char *pfam_id)
{
  struct pfam *pfam_ptr, *found_pfam_ptr;

  /* Initialise variables */
  found_pfam_ptr = NULL;

  /* Get pointer to the first pfam */
  pfam_ptr = first_pfam_ptr;

  /* Loop over all the pfams to find the one just read in */
  while (pfam_ptr != NULL && found_pfam_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (!strcmp(pfam_id,pfam_ptr->pfam_id))
        found_pfam_ptr = pfam_ptr;

      /* Get pointer to the next pfam */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Return the matching pfam pointer */
  return(found_pfam_ptr);
}
/***********************************************************************

create_pfam_entry  -  Create a Pfam record holding name and decription
                      of Pfam domain

***********************************************************************/

struct pfam *create_pfam_entry(struct pfam **fst_pfam_ptr,
                               struct pfam **lst_pfam_ptr,char *pfam_id)
{
  struct pfam *pfam_ptr;
  struct pfam *first_pfam_ptr, *last_pfam_ptr;

  /* Initialise pfam pointers */
  pfam_ptr = NULL;
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_ptr = *lst_pfam_ptr;

  /* Create pfam entry */
  pfam_ptr = (struct pfam *) malloc(sizeof(struct pfam));
  if (pfam_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct pfam\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  store_string(&(pfam_ptr->pfam_id),pfam_id);
  if (!strncmp(pfam_id,"PF",2))
    pfam_ptr->type = PFAMA;
  else if (!strncmp(pfam_id,"PB",2))
    pfam_ptr->type = PFAMB;
  else if (!strncmp(pfam_id,"COV",2))
    pfam_ptr->type = COVID19;
  else
    pfam_ptr->type = UNKNOWN;

  /* Initialise remaining fields */
  pfam_ptr->name = &null;
  pfam_ptr->description = &null;

  /* Initialise pointer to next entry */
  pfam_ptr->next_pfam_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_pfam_ptr == NULL)
    first_pfam_ptr = pfam_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_pfam_ptr->next_pfam_ptr = pfam_ptr;
                      
  /* Save the current residue's pointer */
  last_pfam_ptr = pfam_ptr;

  /* Return current pointers */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;
  return(pfam_ptr);
}
/***********************************************************************

create_domain_entry  -  Create a domain record being a domain attached to
                        a specific UniProt code

***********************************************************************/

struct domain *create_domain_entry(struct domain **fst_domain_ptr,
                                   struct domain **lst_domain_ptr,
                                   int type,struct pfam *pfam_ptr,
                                   char *id,int start_res,int end_res)
{
  int i;

  struct domain *domain_ptr;
  struct domain *first_domain_ptr, *last_domain_ptr;

  /* Initialise domain pointers */
  domain_ptr = NULL;
  first_domain_ptr = *fst_domain_ptr;
  last_domain_ptr = *lst_domain_ptr;

  /* Create domain entry */
  domain_ptr = (struct domain *) malloc(sizeof(struct domain));
  if (domain_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct domain\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the id */
  domain_ptr->id = &null;
  if (id[0] != '\0')
    store_string(&(domain_ptr->id),id);

  /* Store the known data */
  domain_ptr->pfam_ptr = pfam_ptr;
  domain_ptr->start_res = start_res;
  domain_ptr->end_res = end_res;
  domain_ptr->type = type;

  /* Initialise remaining fields */
  domain_ptr->cath_ptr = NULL;
  domain_ptr->domain_no = 0;
  domain_ptr->colour = 0;
  for (i = 0; i < 4; i++)
    domain_ptr->map[i] = -1;

  /* Initialise pointer to next entry */
  domain_ptr->next_domain_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_domain_ptr == NULL)
    first_domain_ptr = domain_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_domain_ptr->next_domain_ptr = domain_ptr;
                      
  /* Save the current residue's pointer */
  last_domain_ptr = domain_ptr;

  /* Return current pointers */
  *fst_domain_ptr = first_domain_ptr;
  *lst_domain_ptr = last_domain_ptr;
  return(domain_ptr);
}
/***********************************************************************

get_pfam_names  -  Read in the pfam names from the
                   pfam_domains.parsed.sorted file

***********************************************************************/

void get_pfam_names(struct pfam **pfam_index_ptr,int npfam,
                    struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char last_pfam_id[PFAM_LEN], pfam_id[PFAM_LEN];
  char pfam_name[NAME_LEN], short_name[NAME_LEN];

  char *string_ptr, *token[MAX_TOKENS];

  int ipfam, itoken, len, nlines, ntokens, nmatched;
  int done, next_code, next_record;

  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  ipfam = 0;
  last_pfam_id[0] = '\0';
  next_code = TRUE;
  next_record = TRUE;
  nlines = nmatched = 0;

  /* Get pointer to the first pfam record */
  if (npfam == 0)
    return;
  pfam_ptr = pfam_index_ptr[ipfam];
  ipfam++;
  next_code = FALSE;

  /* Form name of input file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,"/");
  strcat(file_name,"pfam_names.txt");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      printf("Warning. Unable to open data file [%s]\n",file_name);
      return;
    }

  /* Loop until all codes have been assigned a protein name */
  while (done == FALSE)
    {
      /* If need to go to the next code, then do so */
      if (next_code == TRUE)
        {
          /* Get pointer to the next pfam record */
          if (ipfam < npfam)
            {
              pfam_ptr = pfam_index_ptr[ipfam];
              ipfam++;
            }
          else
            done = TRUE;

          /* Check for end */
          if (pfam_ptr == NULL)
            done = TRUE;

          /* Set flag */
          next_code = FALSE;
        }

      /* If need to read to next record, then loop until we have the
         next useful one */
      while (done == FALSE && next_record == TRUE)
        {
          /* Read in the next record from the file */
          if (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Increment line count */
              nlines++;

              /* Get the tab-delimited tokens on this line */
              ntokens = tokenize(input_line,token,MAX_TOKENS,'\t');

              /* If have 3 tokens, process them */
              if (ntokens == 3)
                {
                  /* Check that pfam id not too long */
                  len = strlen(token[0]);
                  if (len > PFAM_LEN - 1)
                    {
                      printf("*** Warning. Invalid Pfam id: [%s]\n",
                             input_line);
                      return;
                    }

                  /* Get the Pfam id */
                  strcpy(pfam_id,token[0]);

                  /* Check that pfam ids in the file are in sorted
                     order */
                  if (strcmp(pfam_id,last_pfam_id) < 0)
                    {
                      printf("*** Warning. Codes not sorted. "
                              "Current [%s] lt previous [%s]\n",
                              pfam_id,last_pfam_id);
                      return;
                    }

                  /* Save the long and short Pfam names */
                  strcpy(pfam_name,token[1]);
                  strcpy(short_name,token[2]);
                  string_chomp(short_name);

                  /* If this is less than our current pfam record,
                     need to keep reading */
                  if (strcmp(pfam_id,pfam_ptr->pfam_id) < 0)
                    next_record = TRUE;
                  else
                    next_record = FALSE;
                }

              /* Save the current pfam id */
              strcpy(last_pfam_id,pfam_id);
            }

          /* Have reached the end of the file, so we're done */
          else
            done = TRUE;
        }

      /* If not yet done, check the match */
      if (done == FALSE)
        {
          /* If this matches our current pfam record, store the name */
          if (!strcmp(pfam_id,pfam_ptr->pfam_id))
            {
              /* Store the name */
              store_string(&(pfam_ptr->name),pfam_name);
              store_string(&(pfam_ptr->description),short_name);

              /* Increment count of matches */
              nmatched++;

              /* Set both flags to go to the next record */
              next_code = TRUE;
              next_record = TRUE;
            }

          /* Otherwise, if it is higher, need to get the next code */
          else if (strcmp(pfam_id,pfam_ptr->pfam_id) > 0)
            {
              /* Set flags */
              next_code = TRUE;
              next_record = FALSE;
            }

          /* Otherwise, must be lower, so need to read the next
             record */
          else
              next_record = TRUE;
        }
    }

  /* Close the input data file */
  fclose(file_ptr);
}
/***********************************************************************

get_uniprot_data  -  Retrieve the UniProt accession code, id and
                     partial sequence name from the input line

***********************************************************************/

int get_uniprot_data(char *input_line,char *uniprot_acc,
                     char *uniprot_id,char *protein_name)
{
  char number_string[20];
  char *string_ptr;

  int end_pos, i, ipos, jpos, len, seq_len, start_pos, tlen;
  int new_format, newest_format;

  /* Initialise variables */
  seq_len = 0;
  uniprot_acc[0] = uniprot_id[0] = protein_name[0] = '\0';
  newest_format = new_format = FALSE;

  /* Check for the latest format */
  if (!strncmp(input_line,">sp|",4) || !strncmp(input_line,">tr|",4))
    newest_format = TRUE;

  /* Determine whether code is in new or old format */
  else if ((string_ptr = strchr(input_line,'|')) != NULL &&
           (string_ptr - input_line) < 16)
    new_format = TRUE;
  else
    new_format = FALSE;

  /* Latest format is >sp|id|code name or >tr|id|code name */
  if (newest_format == TRUE)
    {
      /* Store the accession number */
      strncpy(uniprot_acc,input_line+4,UNIPROT_LEN);
      uniprot_acc[UNIPROT_LEN]= '\0';
      string_ptr = strchr(uniprot_acc,'|');
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Extract the UniProt id */
      string_ptr = strchr(input_line+4,'|');
      if (string_ptr != NULL)
        {
          strncpy(uniprot_id,string_ptr+1,UNIPROT_LEN);
          uniprot_id[UNIPROT_LEN]= '\0';
          string_ptr = strstr(uniprot_id," ");
          if (string_ptr != NULL)
            string_ptr[0] = '\0';
        }
      else
        uniprot_id[0]= '\0';

      /* Get the sequence name */
      string_ptr = strstr(input_line," ");
      if (string_ptr != NULL)
        strcpy(protein_name,string_ptr + 1);

      /* Get the length of the sequence */
      if ((string_ptr = strstr(protein_name," aa)")) != NULL)
        {
          /* Get position of " aa" */
          end_pos = string_ptr - protein_name;
          start_pos = 0;

          /* Loop to find space prior to start of number */
          for (i = 0; i < end_pos; i++)
            {
              /* If this is a space, stroe its position */
              if (protein_name[i] == ' ')
                start_pos = i + 1;
            }

          /* Get the length of the number */
          start_pos++;
          len = end_pos - start_pos;

          /* Retrieve the number */
          if (len < 14)
            {
              strncpy(number_string,protein_name+start_pos,len);
              number_string[len] = '\0'; 
              seq_len = atoi(number_string);
            }
        }
    }

  /* If new format, code is in form >id|code name */
  else if (new_format == TRUE)
    {
      /* Store the accession number */
      strncpy(uniprot_acc,input_line+1,UNIPROT_LEN);
      uniprot_acc[UNIPROT_LEN]= '\0';
      string_ptr = strchr(uniprot_acc,'|');
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Extract the UniProt id */
      string_ptr = strchr(input_line,'|');
      if (string_ptr != NULL)
        {
          strncpy(uniprot_id,string_ptr+1,UNIPROT_LEN);
          uniprot_id[UNIPROT_LEN]= '\0';
          string_ptr = strstr(uniprot_id," ");
          if (string_ptr != NULL)
            string_ptr[0] = '\0';
        }
      else
        uniprot_id[0]= '\0';

      /* Get the protein name */
      string_ptr = strchr(input_line,' ');
      if (string_ptr != NULL)
        strcpy(protein_name,string_ptr+1);
    }

  /* Otherwise, got old format, code is >code (id) */
  else
    {
      /* Extract the UniProt id */
      strncpy(uniprot_id,input_line+1,UNIPROT_LEN);
      uniprot_id[UNIPROT_LEN]= '\0';
      string_ptr = strstr(uniprot_id," ");
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Get UniProt accession (between brackets) */
      string_ptr = strstr(input_line," (");
      if (string_ptr != NULL)
        {
          /* Get start position of UniProt id */
          ipos = string_ptr - input_line + 2;

          /* Get end position */
          string_ptr = strstr(input_line + ipos,")");
          if (string_ptr != NULL)
            {
              /* Get length of id */
              jpos = string_ptr - input_line;
              tlen = jpos - ipos;

              /* Extract the UniProt id */
              strncpy(uniprot_acc,input_line + ipos,tlen);
              uniprot_acc[len] = '\0';

              /* Extract the protein name */
              if (jpos + 2 < len)
                strcpy(protein_name,string_ptr+2);
            }
        }
    }

  /* Return the sequence length */
  return(seq_len);
}
/***********************************************************************

get_sequence  -  Get the sequence from the given .fasta file

***********************************************************************/

int get_sequence(char *file_name,char **seq,char *name)
{
  char input_line[LINELEN + 1];
  char uniprot_acc[UNIPROT_LEN], uniprot_id[UNIPROT_LEN];

  char *sequence;

  int len, line, seqlen;

  FILE *file_ptr;

  /* Initialise variables */
  line = seqlen = 0;
  *seq = sequence = &null;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    return(seqlen);

  /* Loop while reading in from file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Increment line count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* If first line, retrieve the protein name from it */
      if (line == 1)
        get_uniprot_data(input_line,uniprot_acc,uniprot_id,name);

      /* If second line, then save sequence */
      else if (line == 2)
        store_string(&sequence,input_line);

      /* If subsquent, then append */
      else if (line > 2)
        append_string(&sequence,input_line);
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the sequence read in */
  *seq = sequence;

  /* Get its length */
  seqlen = strlen(sequence);

/* Return the sequence length */
return(seqlen);
}
/***********************************************************************

align_seqs  -  Read in the PDB and UniProt sequence files and align the
               two sequences

***********************************************************************/

void align_seqs(struct code *code_ptr,char *file_name,int pdb_type,
                int pdbsum_type,struct c_file_data file_name_ptr)
{
  char seq_file[FILENAME_LEN], name[LINELEN + 1];
  char pdbsum_dir[FILENAME_LEN + 1];

  char *alignment[2], *sequence[2];

  int dir_type, pdb_len, uniprot_len;

  struct chain *chain_ptr;
  struct clist *clist_ptr;
  struct pdb *pdb_ptr;

  /* Get the corresponding UniProt sequence */
  printf("In align_seqs for %s ...\n",code_ptr->uniprot_acc);
  fflush(stdout);
  uniprot_len = get_sequence(file_name,&(sequence[UNIPROT_SEQ]),name);

  /* Save the protein name and sequence length */
  store_string(&(code_ptr->name),name);
  code_ptr->seq_len = uniprot_len;
  code_ptr->have_seq = TRUE;

  /* Get pointer to the first clist entry for this UniProt code */
  clist_ptr = code_ptr->first_clist_ptr;

  /* Loop through all entries to compare PDB chain sequence against
     UniProt code sequence */
  while (clist_ptr != NULL)
    {
      /* Get the corresponding chain */
      chain_ptr = clist_ptr->chain_ptr;

      /* Get the corresponding PDB record */
      pdb_ptr = chain_ptr->pdb_ptr;

      /* Form name of seqXXX.fasta file for this chain */
      if (pdb_type == PDBSUM1)
        sprintf(seq_file,"../newsec/seq%3.3d.fasta",chain_ptr->chain_no);
      else
        {
          /* Get the PDBsum directory */
          dir_type = find_pdbsum_dir(pdb_ptr->pdb_code,
                                     file_name_ptr.pdbsum_template,
                                     pdbsum_type,pdbsum_dir);

          /* Form the neame of the FASTA sequence file */
          sprintf(seq_file,"%sseq%3.3d.fasta",pdbsum_dir,
                  chain_ptr->chain_no);
        }

      /* Get the PDB sequence */
      pdb_len = get_sequence(seq_file,&(sequence[PDB_SEQ]),name);

      /* If neither sequence is empty, then perform alignment */
      if (pdb_len > 0 && uniprot_len > 0)
        {
          /* Align the two sequences */
          chain_ptr->identity
            = align_sequences(sequence[PDB_SEQ],sequence[UNIPROT_SEQ],
                              chain_ptr->alignment,
                              &(chain_ptr->align_len),
                              &(chain_ptr->nidentical));
        }

      /* Get pointer to the next entry in the linked list */
      clist_ptr = clist_ptr->next_clist_ptr;
    }
  printf("End of align_seqs\n");
}
/***********************************************************************

download_uniprot_seqs  -  Use wget to retrieve the required UniProt
                          FASTA files and align the sequences to the PDB
                          chain sequences

***********************************************************************/

void download_uniprot_seqs(struct code **code_index_ptr,int ncodes,
                           int pdb_type,int pdbsum_type,
                           struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], url[FILENAME_LEN];
  char command[FILENAME_LEN + 1];

  char *pdbsum_dir, *string_ptr;

  int dir_type, icode, loop;
  int dos, have_file, itry, ntries;
  int ok;

  struct chain *chain_ptr;
  struct code *code_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  ntries = 5;
  ok = FALSE;

  /* Loop over all the codes */
  for (icode = 0; icode < ncodes; icode++)
    {
      /* Initialise flag */
      ok = FALSE;

      /* Get the corresponding code */
      code_ptr = code_index_ptr[icode];

      /* If not the dummy record, download the sequence */
      if (strcmp(code_ptr->uniprot_acc,"NONE"))
        {
          /* Form the name of the output .fasta file */
          sprintf(file_name,"%s.fasta",code_ptr->uniprot_acc);

          /* See if we already have this file */
          have_file = empty_check(file_name);
          if (have_file == FALSE)
            {
              /* Form the URL of the required .fasta file */
              strcpy(url,file_name_ptr.other_dir[GET_UNIPROT_DATA]);
              strcat(url,code_ptr->uniprot_acc);
              strcat(url,".fasta");

              /* Initialise flag */
              ok = FALSE;

              /* Loop while trying to download the file */
              for (itry = 0; itry < ntries && ok == FALSE; itry++)
                {
                  /* Run the wget command to retrieve the Pfam data */
                  printf("Running wget %d ...\n",itry + 1);
                  fflush(stdout);
                  ok = run_wget(file_name_ptr.exe_name[WGET],url,
                                file_name,dos);
                  printf("... done with ok = %s\n",T_F[ok]);
                  fflush(stdout);

                  /* If not OK, pause before trying again */
                  if (ok == FALSE)
                    {
                      /* Pause for 5 seconds */
                      sprintf(command,"%s",SLEEP_COMMAND);
                      for (loop = 0; loop < 5; loop++)
                        call_system(command,dos);
                    }
                }
            }
          else
            ok = TRUE;

          /* If OK, align the PDB sequence to the UniProt sequence */
          if (ok == TRUE)
            {
              /* Perform the alignments */
              align_seqs(code_ptr,file_name,pdb_type,pdbsum_type,
                         file_name_ptr);

              /* Delete the downloaded .fasta file */
              sprintf(command,"%s %s",RM_COMMAND,file_name);
            }
        }
    }
}
/***********************************************************************

read_swisspfam_data  -  Pick up the Pfam domains associated with each
                        UniProt sequences from the swisspfam_data.txt.gz
                        file

***********************************************************************/

int read_swisspfam_data(struct code **code_index_ptr,int ncodes,
                        struct pfam **fst_pfam_ptr,
                        struct pfam **lst_pfam_ptr,int entry_types,
                        struct parameters params,
                        struct c_file_data file_name_ptr)
{
  char ch, file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char tmp_line[LINELEN + 1];
  char number_string[20];
  char pfam_id[PFAM_LEN];
  char uniprot_acc[UNIPROT_LEN], uniprot_id[UNIPROT_LEN];

  char *string_ptr, *token[MAX_TOKEN];

  int dotpos, end_res, err_no, i, icopy, idpos, ipos, len;
  int ncopy, ndom, ndots, npfam, ntokens, seq_len, start_res;
  int floop, last_loop, loop;
  int numeric;

  long line, ndone, nseqs;

  struct code *code_ptr;
  struct domain *domain_ptr, *last_domain_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  FILE *auditfile_ptr;
  gzFile gzfile_ptr;

  /* Initialise variables */
  line = 0;
  ndom = 0;
  ndone = 0;
  npfam = 0;
  nseqs = 0;
  seq_len = 0;
  uniprot_id[0] = '\0';

  /* Initialise pointers */
  code_ptr = NULL;
  last_domain_ptr = NULL;
  pfam_ptr = first_pfam_ptr = last_pfam_ptr = NULL;

  /* Form name of input file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"swisspfam_data.txt.gz");
  printf("Reading Pfam domains file %s ...\n",file_name);
  fflush(stdout);

  /* Open the file, if it exists */
  gzfile_ptr = gzopen(file_name,"rb");
  if (gzfile_ptr ==  NULL)
    {
      /* Check for insufficiency of memory */
      if(err_no)
        {
          printf("*** ERROR. Insufficient memory opening file [%s]\n",
                 file_name);
          return(npfam);
        }
      
      /* Otherwise, error must be due to missing file */
      printf("\n*** Warning. Unable to open file [%s]\n",file_name);
      return(npfam);
    }

  /* Loop while reading in records from the file */
  while (gzgets(gzfile_ptr,input_line,LINELEN)!= NULL && ndone < ncodes)
    {
      /* Increment line count */
      line++;

      /* Truncate the input line */
      len = string_chomp(input_line);
              
      /* Increment count of sequence read in */
      nseqs++;

      /* Write to audit file, if required */
      if ((nseqs % NAUDIT) == 0)
        {
          /* Write sequence count to audit file to allow
             check on progress */
          auditfile_ptr = open_output_file("audit.out",TRUE);
          fprintf(auditfile_ptr,"Pfam Seqs %ld\n",nseqs);
          fclose(auditfile_ptr);
        }
          
      /* Split line up into space-separated tokens */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

      /* If have enough tokens, check if this UniProt record is wanted */
      if (ntokens > 4)
        {
          /* Get the UniProt accession */
          strcpy(uniprot_acc,token[0]);

          /* See if it is one of the ones we're interested in */
          code_ptr = binary_search(code_index_ptr,ncodes,uniprot_acc);

          /* If record found, save the sequence length as given in the
             Pfam record */
          if (code_ptr != NULL)
            {
              /* Store the details */
              strcpy(uniprot_id,token[1]);
              strcpy(pfam_id,token[2]);

              /* Find whether we already have this Pfam record */
              pfam_ptr = find_pfam_entry(first_pfam_ptr,pfam_id);

              /* If not found, create a new entry */
              if (pfam_ptr == NULL)
                {
                  /* Create record */
                  pfam_ptr
                    = create_pfam_entry(&first_pfam_ptr,&last_pfam_ptr,
                                        pfam_id);

                  /* Increment count of Pfam domains */
                  npfam++;
                }

              /* Get the start and end residues of this domain */
              start_res = atoi(token[3]);
              end_res = atoi(token[4]);

              /* Create a domain record */
              domain_ptr
                = create_domain_entry(&(code_ptr->first_domain_ptr),
                                      &last_domain_ptr,PFAM,
                                      pfam_ptr,pfam_ptr->pfam_id,
                                      start_res,end_res);

              /* Increment count of domains */
              ndom++;
            }
        }
              
      /* Free the tokens */
      free_tokens(token,ntokens);
    }

  /* Close file */
  gzclose(gzfile_ptr);

  /* Show stats */
  printf("Number of Pfam records read in:    %8ld\n",nseqs);
  printf("Number of Pfam records stored:     %8d\n",npfam);
  printf("Number of domain records created:  %8d\n",ndom);
  fflush(stdout);

  /* Return pointers to Pfam records */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;

  /* Return number of Pfam domains created */
  return(npfam);
}
/***********************************************************************

read_swisspfam  -  Pick up the Pfam domains associated with each
                   UniProt sequences from the swisspfam.gz file

***********************************************************************/

int read_swisspfam(struct code **code_index_ptr,int ncodes,
                   struct pfam **fst_pfam_ptr,
                   struct pfam **lst_pfam_ptr,int entry_types,
                   struct parameters params,
                   struct c_file_data file_name_ptr)
{
  char ch, file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char tmp_line[LINELEN + 1];
  char number_string[20];
  char pfam_id[PFAM_LEN];
  char uniprot_acc[UNIPROT_LEN], uniprot_id[UNIPROT_LEN];

  char *string_ptr, *token[MAX_TOKEN];

  int dotpos, end_res, err_no, i, icopy, idpos, ipos, len;
  int ncopy, ndom, ndots, npfam, ntokens, seq_len, start_res;
  int floop, last_loop, loop;
  int numeric;

  long line, ndone, nseqs;

  struct code *code_ptr;
  struct domain *domain_ptr, *last_domain_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  FILE *auditfile_ptr;
#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  line = 0;
  ndom = 0;
  ndone = 0;
  npfam = 0;
  nseqs = 0;
  seq_len = 0;
  uniprot_id[0] = '\0';

  /* Initialise pointers */
  code_ptr = NULL;
  last_domain_ptr = NULL;
  pfam_ptr = first_pfam_ptr = last_pfam_ptr = NULL;

  /* Determine whether COVID-19 swisspfam file needs to be read */
  floop = 1;
  last_loop = 2;
  if (params.covid19 == TRUE || entry_types == COVID_ONLY)
    {
      floop = 0;
      last_loop = 1;
    }
  else if (entry_types == MIXED)
    {
      floop = 0;
      last_loop = 2;
    }

  /* Loop over the swisspfam files to be read */
  for (loop = floop; loop < last_loop; loop++)
    {

#ifdef TEST
      /* Form name of input file */
      if (loop == 0)
        strcpy(file_name,"swisspfam_covid19");
      else
        strcpy(file_name,"swisspfam");
      printf("Reading Pfam domains file %s ...\n",file_name);

      /* Open the file, if it exists */
      if ((file_ptr = fopen(file_name,"r")) ==  NULL)
        {
          printf("\n*** Warning. Unable to open file [%s]\n",file_name);
          return(npfam);
        }

      /* Loop while reading in records from input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL && ndone < ncodes)
#else
      /* Form name of input file */
      strcpy(file_name,file_name_ptr.other_dir[PFAM_DATA_DIR]);
      strcat(file_name,DIR_SEP);
      if (loop == 0)
        strcat(file_name,"swisspfam_covid19");
      else
        strcat(file_name,"swisspfam.gz");
      printf("Reading Pfam domains file %s ...\n",file_name);
      fflush(stdout);

      /* Open the file, if it exists */
      gzfile_ptr = gzopen(file_name,"rb");
      if (gzfile_ptr ==  NULL)
        {
          /* Check for insufficiency of memory */
          if(err_no)
            {
              printf("*** ERROR. Insufficient memory opening file [%s]\n",
                     file_name);
              return(npfam);
            }

          /* Otherwise, error must be due to missing file */
          printf("\n*** Warning. Unable to open file [%s]\n",file_name);
          return(npfam);
        }

      /* Loop while reading in records from the file */
      while (gzgets(gzfile_ptr,input_line,LINELEN)!= NULL && ndone < ncodes)
#endif
        {
          /* Increment line count */
          line++;

          /* If this is the start of a new UniProt entry, get code and check
             whether it is one of ours */
          if (input_line[0] == '>')
            {
              /* Truncate the input line */
              len = string_chomp(input_line);
              
              /* Increment count of sequence read in */
              nseqs++;

              /* Write to audit file, if required */
              if ((nseqs % NAUDIT) == 0)
                {
                  /* Write sequence count to audit file to allow
                     check on progress */
                  auditfile_ptr = open_output_file("audit.out",TRUE);
                  fprintf(auditfile_ptr,"Pfam Seqs %ld\n",nseqs);
                  fclose(auditfile_ptr);
                }

              /* If we've just done a code record, increment count */
              if (code_ptr != NULL)
                ndone++;

              /* Initialise UniProt id */
              uniprot_acc[0] = '\0';
              uniprot_id[0] = '\0';
              code_ptr = NULL;
              seq_len = 0;

              /* Extract the UniProt id */
              strcpy(tmp_line,input_line+1);
              string_ptr = strstr(tmp_line," ");
              if (string_ptr != NULL)
                string_ptr[0]= '\0';
              len = strlen(tmp_line);
              if (len < UNIPROT_LEN)
                strcpy(uniprot_id,tmp_line);

              /* Initialise variables */
              last_domain_ptr = NULL;

              /* Find location of UniProt accession code */
              string_ptr = strstr(input_line,"=| ");
              if (string_ptr != NULL)
                {
                  /* Extract the sequence length */
                  strcpy(tmp_line,string_ptr+3);
                  string_ptr = strstr(tmp_line," a.a.");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  string_ptr = strchr(tmp_line,' ');
                  if (string_ptr != NULL && strlen(string_ptr) < 15)
                    {
                      /* Extract length of sequence */
                      strcpy(number_string,string_ptr);
                      seq_len = atoi(number_string);

                      /* Terminate the string just after Pfam number */
                      string_ptr[0] = '\0';
                    }

                  /* Extract the UniProt accession code */
                  string_ptr = strchr(tmp_line,'.');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  strcpy(uniprot_acc,tmp_line);
                }

              /* If have managed to extract the UniProt accession, see if
                 it is one of the ones we're interested in */
              if (uniprot_acc[0] != '\0' && ncodes > 0)
                code_ptr = binary_search(code_index_ptr,ncodes,uniprot_acc);
              else
                code_ptr = NULL;

              /* If record found, save the sequence length as given in the
                 Pfam record */
              if (code_ptr != NULL)
                {
                  /* If have a UniProt id, store it */
                  if (uniprot_id[0] != '\0')
                    strcpy(code_ptr->uniprot_id,uniprot_id);

                  /* Check that we haven't already processed this record */
                  if (code_ptr->first_domain_ptr == NULL)
                    code_ptr->pfam_seq_len = seq_len;

                  /* Otherwise, leave it */
                  else
                    code_ptr = NULL;
                }
            }

          /* If this is the start of a new UniProt entry, get code and check
             whether it is one of ours */
          else if (code_ptr != NULL && input_line[0] != '\0')
            {
              /* Truncate the input line */
              len = string_chomp(input_line);

              /* Split line up into space-separated tokens */
              ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');
              if (ntokens < 0)
                printf("*** Warning. Pfam entry for %s has too many repeats!\n",
                       code_ptr->uniprot_acc);

              /* Get number of copies of this domain */
              if (ntokens > 1)
                ncopy = atoi(token[1]);

              /* Loop through the tokens to identify the one containing the
                 Pfam id */
              idpos = -1;
              for (i = 0; i < ntokens && idpos == -1; i++)
                {
                  /* Check for PF or PB at start and rest numeric with a
                     single dot */
                  if (!strncmp(token[i],"PF",2) || !strncmp(token[i],"PB",2))
                    {
                      /* Initialise flags */
                      ndots = 0;
                      dotpos = 0;
                      numeric = TRUE;

                      /* Check remaining characters are numeric and there
                         is a single dot */
                      for (ipos = 2; ipos < strlen(token[i]); ipos++)
                        {
                          /* Get this character */
                          ch = token[i][ipos];

                          /* If dot, then increment count of dots and
                             store position */
                          if (ch == '.')
                            {
                              ndots++;
                              dotpos = ipos;
                            }

                          /* Check if numeric */
                          else if (ch < '0' || ch > '9')
                            numeric = FALSE;
                        }

                      /* If valid, store position of Pfam id */
                      if (ndots < 2 && numeric == TRUE && dotpos < PFAM_LEN)
                        idpos = i;
                    }

                  /* Check for COVID-19 domain */
                  else if (!strncmp(token[i],"COV",3))
                    {
                      /* Initialise flags */
                      ndots = 0;
                      dotpos = 0;

                      /* Store position of Pfam id */
                      idpos = i;
                    }
                }
              
              /* If Pfam id identified, extract domain information */
              if (idpos > -1)
                {
                  /* Store the Pfam id */
                  if (dotpos > 0)
                    token[idpos][dotpos] = '\0';
                  strcpy(pfam_id,token[idpos]);

                  /* Find whether we already have this Pfam record */
                  pfam_ptr = find_pfam_entry(first_pfam_ptr,pfam_id);

                  /* If not found, create a new entry */
                  if (pfam_ptr == NULL)
                    {
                      /* Create record */
                      pfam_ptr
                        = create_pfam_entry(&first_pfam_ptr,&last_pfam_ptr,
                                            pfam_id);

                      /* Increment count of Pfam domains */
                      npfam++;

                      /* Store domain's short name */
                      store_string(&(pfam_ptr->name),token[0]);

                      /* Form the long name and store as the description */
                      tmp_line[0] = '\0';
                      for (i = idpos + 1; i < ntokens - ncopy; i++)
                        {
                          /* If not first token, add space */
                          if (tmp_line[0] != '\0')
                            strcat(tmp_line," ");

                          /* Add current token to description */
                          strcat(tmp_line,token[i]);
                        }

                      /* Store description */
                      if (tmp_line[0] != '\0')
                        store_string(&(pfam_ptr->description),tmp_line);
                    }

                  /* Loop over the copies of this domain to create domain
                     entries for current UniProt code */
                  for (icopy = 0; icopy < ncopy; icopy++)
                    {
                      /* Get position of residue range for this copy of the
                         domain */
                      i = ntokens - ncopy + icopy;

                      /* Extract the start and end residue numbers */
                      if (strlen(token[i]) < 14)
                        {
                          /* Initialise start and end residue numbers */
                          start_res = 0;
                          end_res = 0;

                          /* Extract start and end residues */
                          string_ptr = strchr(token[i],'-');
                          if (string_ptr != NULL)
                            {
                              /* Extract end residue */
                              strcpy(number_string,string_ptr+1);
                              end_res = atoi(number_string);

                              /* Extract start residue */
                              string_ptr[0] = '\0';
                              start_res = atoi(token[i]);
                            }

                          /* Create a domain record */
                          domain_ptr
                            = create_domain_entry(&(code_ptr->first_domain_ptr),
                                                  &last_domain_ptr,PFAM,
                                                  pfam_ptr,pfam_id,
                                                  start_res,end_res);

                          /* Increment count of domains */
                          ndom++;
                        }
                    }
                }
              
              /* Free the tokens */
              free_tokens(token,ntokens);
            }
        }

      /* Close file */
#ifdef TEST
      fclose(file_ptr);
#else
      gzclose(gzfile_ptr);
#endif
    }

  /* Show stats */
  printf("Number of Pfam records read in:    %8ld\n",nseqs);
  printf("Number of Pfam records stored:     %8d\n",npfam);
  printf("Number of domain records created:  %8d\n",ndom);
  fflush(stdout);

  /* Return pointers to Pfam records */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;

  /* Return number of Pfam domains created */
  return(npfam);
}
/***********************************************************************

read_pfam_xml  -  Read the output XML file, pfam.xml

***********************************************************************/

int read_pfam_xml(struct pfam **fst_pfam_ptr,struct pfam **lst_pfam_ptr,
                  int npfam,struct code *code_ptr,char *file_name)
{
  char input_line[LINELEN + 1];
  char pfam_id[PFAM_LEN + 1], number_string[20], tmp[LINELEN + 1];

  char *sequence, *string_ptr, *token[MAX_TOKEN];

  int itoken, len, line, ndomains, ntokens, type;
  int end_res, start_res, start_seq;

  struct domain *domain_ptr, *last_domain_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  sequence = &null;

  /* Initialise pointers */
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_ptr = *lst_pfam_ptr;
  last_domain_ptr = NULL;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. Unable to find file [%s]\n",file_name);
      return(npfam);
    }

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment line-count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Check for the protein sequence */
      if ((string_ptr = strstr(input_line,"<sequence ")) != NULL)
        {
          /* Find where the sequence starts */
          string_ptr = strchr(string_ptr,'>');

          /* If found, then save the start-point */
          if (string_ptr != NULL)
            start_seq = string_ptr - input_line + 1;
          else
            start_seq = 0;

          /* Remove the end tag */
          string_ptr = strchr(input_line + start_seq,'<');
          if (string_ptr != NULL)
            string_ptr[0] = '\0';

          /* Store the sequence */
          store_string(&sequence,input_line + start_seq);

          /* Save the sequence length */
          code_ptr->pfam_seq_len = strlen(sequence);

          /* Free up memory taken by the sequence */
          free(sequence);
          sequence = &null;
        }

      /* Check for a domain-match line */
      else if ((string_ptr = strstr(input_line,"<match ")) != NULL)
        {
          /* Initialise data */
          pfam_id[0] = tmp[0] = '\0';
          type = PFAMB;
          pfam_ptr = NULL;

          /* Break line up into tokens */
          ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

          /* Loop over the tokens to extract the domain details */
          for (itoken = 0; itoken < ntokens; itoken++)
            {
              /* If this is the Pfam id, then store */
              if (!strncmp(token[itoken],"accession=",10))
                {
                  /* Store the Pfam id */
                  strcpy(pfam_id,token[itoken] + 11);

                  /* Remove the final double-quote */
                  len = strlen(pfam_id);
                  pfam_id[len - 1] = '\0';
                }

              /* If it's the Pfam name, then store */
              else if (!strncmp(token[itoken],"id=",3))
                {
                  /* Store the description */
                  strcpy(tmp,token[itoken] + 4);

                  /* Remove the final double-quote */
                  len = strlen(tmp);
                  tmp[len - 1] = '\0';
                }

              /* If it's the Pfam type, then get the type */
              else if (!strncmp(token[itoken],"type=",5))
                {
                  /* Get the type */
                  if (strstr(token[itoken],"Pfam-A") != NULL)
                    type = PFAMA;
                  else
                    type = PFAMB;
                }
            }

          /* See if we already have this Pfam domain */
          pfam_ptr = find_pfam_entry(first_pfam_ptr,pfam_id);

          /* If not found, create Pfam record */
          if (pfam_ptr == NULL && pfam_id[0] != '\0')
            {
              /* Create a Pfam record */
              pfam_ptr
                = create_pfam_entry(&first_pfam_ptr,&last_pfam_ptr,
                                    pfam_id);

              /* Store name */
              store_string(&(pfam_ptr->name),tmp);

              /* Increment coiunt of Pfam domains stored */
              npfam++;
            }

          /* Free the tokens */
          free_tokens(token,ntokens);
        }

      /* Check for a residue-range line */
      else if ((string_ptr = strstr(input_line,"<location ")) != NULL)
        {
          /* Initialise domain range */
          start_res = end_res = -1;

          /* Break line up into tokens */
          ntokens = get_tokens(input_line,token,MAX_TOKEN,' ');

          /* Loop over the tokens to extract the domain details */
          for (itoken = 0; itoken < ntokens; itoken++)
            {
              /* If this is the start of the range, then save */
              if (!strncmp(token[itoken],"start=",6))
                {
                  /* Extract the start */
                  strcpy(number_string,token[itoken] + 7);

                  /* Remove the final double-quote */
                  len = strlen(number_string);
                  number_string[len - 1] = '\0';

                  /* Get the start */
                  start_res = atoi(number_string);
                }

              /* If this is the end of the range, then save */
              else if (!strncmp(token[itoken],"end=",4))
                {
                  /* Extract the end */
                  strcpy(number_string,token[itoken] + 5);

                  /* Remove the final double-quote */
                  len = strlen(number_string);
                  number_string[len - 1] = '\0';

                  /* Get the end */
                  end_res = atoi(number_string);
                }
            }

          /* Have a valid residue range, then create a domain record */
          if (end_res > start_res && start_res != -1)
            {
              /* Create a domain record */
              domain_ptr
                = create_domain_entry(&(code_ptr->first_domain_ptr),
                                      &last_domain_ptr,PFAM,pfam_ptr,
                                      pfam_id,start_res,end_res);

              /* Increment coiunt of Domain domains stored */
              ndomains++;
            }

          /* Free the tokens */
          free_tokens(token,ntokens);
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;

  /* Return the number of domains */
  return(npfam);
}
/***********************************************************************

read_pfam_tsv  -  Read the output .tsv file returned by InterPro api

***********************************************************************/

int read_pfam_tsv(char *file_name,struct domain **fst_domain_ptr,
                  char *input_line)
{
  char pfam_id[PFAM_LEN + 1], number_string[20], pfam_name[LINELEN];
  char accession[LINELEN];

  char *string_ptr, *token[MAX_TOKEN];

  int itoken, jtoken, len, line, ndomains, ntokens;
  int end_res, start_res, start_seq;
  int done, locations;

  struct domain *domain_ptr, *first_domain_ptr, *last_domain_ptr;
  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  locations = FALSE;
  ndomains = line = 0;

  /* Initialise pointers */
  *fst_domain_ptr = first_domain_ptr = last_domain_ptr = NULL;
  pfam_ptr = NULL;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. Unable to find file [%s]\n",file_name);
      return(ndomains);
    }

  /* Loop through the file processing each record */
  while (fgets(input_line,LONG_LINELEN,file_ptr) != NULL)
    {
      /* Increment line-count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* Break line up into tokens */
      ntokens = get_tokens(input_line,token,MAX_TOKEN,'"');

      /* Initialise flag */
      locations = FALSE;
      start_res = end_res = -1;
      pfam_id[0] = '\0';
      pfam_name[0] = '\0';

      /* Loop over the tokens to extract the domain details */
      for (itoken = 0; itoken < ntokens; itoken++)
        {
          /* Check for Pfam id and name */
          if (!strcmp(token[itoken],"accession") && itoken + 6 < ntokens)
            {
              /* Check if this refers to a Pfam doma1n */
              if ((!strncmp(token[itoken + 2],"PF",2) ||
                  !strncmp(token[itoken + 2],"PB",2)) &&
                  !strncmp(token[itoken + 6],"pfam",4))
                {
                  /* Get the Pfam id */
                  strcpy(pfam_id,token[itoken + 2]);
                }
              else
                pfam_id[0] = '\0';

              /* Unset locations flag */
              locations = FALSE;
              start_res = end_res = -1;
            }

          /* Check for start of locations */
          else if (!strcmp(token[itoken],"locations") &&
                   pfam_id[0] != '\0')
            locations = TRUE;
          
          /* Check for domain start */
          else if (!strcmp(token[itoken],"start") && itoken + 1 < ntokens)
            {
              /* Extract the start-position */
              strcpy(number_string,token[itoken + 1] + 1);

              /* Remove the final comma */
              len = strlen(number_string);
              number_string[len - 1] = '\0';

              /* Get the start */
              start_res = atoi(number_string);
            }
          
          /* If this is the end of the range, then save */
          else if (!strcmp(token[itoken],"end") && itoken + 1 < ntokens &&
                   locations == TRUE)
            {
              /* Extract the end */
              strcpy(number_string,token[itoken + 1] + 1);

              /* Remove the final comma */
              len = strlen(number_string);
              number_string[len - 1] = '\0';
              
              /* Get the end */
              end_res = atoi(number_string);

              /* Create a domain record */
              domain_ptr
                = create_domain_entry(&first_domain_ptr,
                                      &last_domain_ptr,PFAM,pfam_ptr,
                                      pfam_id,start_res,end_res);

              /* Increment count of domains stored */
              ndomains++;
            }
        }

      /* Free the tokens */
      free_tokens(token,ntokens);
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_domain_ptr = first_domain_ptr;

  /* Return the number of domains */
  return(ndomains);
}
/***********************************************************************

delete_domains  -  Delete all the temporary domain records

***********************************************************************/

void delete_domains(struct domain *first_domain_ptr)
{
  struct domain *domain_ptr, *next_domain_ptr;

  /* Get first domain record */
  domain_ptr = first_domain_ptr;

  /* Loop through the domains until done */
  while (domain_ptr != NULL)
    {
      /* Get pointer to the next domain record */
      next_domain_ptr = domain_ptr->next_domain_ptr;

      /* Free up memory taken by current domain */
      free(domain_ptr);

      /* Go to the next domain */
      domain_ptr = next_domain_ptr;
    }
}
/***********************************************************************

get_pfam_tsv  -  Read the output .tsv file returned by InterPro api

***********************************************************************/

int get_pfam_tsv(struct pfam **fst_pfam_ptr,struct pfam **lst_pfam_ptr,
                 int npfam,struct code *code_ptr,char *file_name,
                 int *nlines,char *long_line)
{
  char pfam_id[PFAM_LEN + 1], number_string[20], pfam_name[LINELEN];
  char accession[LINELEN];

  char *sequence, *string_ptr, *token[MAX_TOKEN];

  int itoken, jtoken, len, line, ndomains, nread, ntokens, type;
  int end_res, start_res, start_seq;
  int done;

  struct domain *domain_ptr, *first_domain_ptr, *last_domain_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  *nlines = line = nread = 0;
  sequence = &null;

  /* Initialise pointers */
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_ptr = *lst_pfam_ptr;
  domain_ptr = first_domain_ptr = last_domain_ptr = NULL;

  /* Read in the Pfam domains */
  ndomains = read_pfam_tsv(file_name,&first_domain_ptr,long_line);

  /* Get pointer to the first domain */
  domain_ptr = first_domain_ptr;

  /* Loop over all the domains to link to Pfam records */
  while (domain_ptr != NULL)
    {
      /* Get the Pfam id */
      strcpy(pfam_id,domain_ptr->id);
      
      /* See if we already have this Pfam domain */
      pfam_ptr = find_pfam_entry(first_pfam_ptr,pfam_id);

      /* If not found, create Pfam record */
      if (pfam_ptr == NULL && pfam_id[0] != '\0')
        {
          /* Create a Pfam record */
          pfam_ptr
            = create_pfam_entry(&first_pfam_ptr,&last_pfam_ptr,
                                pfam_id);

          /* Increment coiunt of Pfam domains stored */
          npfam++;
          nread++;
        }
      
      /* Create a domain record */
      create_domain_entry(&(code_ptr->first_domain_ptr),
                          &last_domain_ptr,PFAM,pfam_ptr,pfam_id,
                          domain_ptr->start_res,domain_ptr->end_res);

      /* Get pointer to the next domain */
      domain_ptr = domain_ptr->next_domain_ptr;
    }

  /* Free uo memory taken by the temporary domains */
  delete_domains(first_domain_ptr);
  
  /* Return current pointers */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;

  /* Return the line count */
  *nlines = line;

  /* Show number of Pfam domains foubd */
  printf("Number of Pfam domains: %d\n",nread);
  fflush(stdout);

  /* Return the number of domains */
  return(npfam);
}
/***********************************************************************

get_missing_domains  -  For any UniProt sequences that have no, or
                        mismatched, Pfam domain info, use the Pfam web
                        server to retrive the domain data

***********************************************************************/

int get_missing_domains(struct code **code_index_ptr,int ncodes,
                        struct pfam **fst_pfam_ptr,
                        struct pfam *last_pfam_ptr,int npfam,
                        struct c_file_data file_name_ptr)
{
  char command[FILENAME_LEN], file_name[FILENAME_LEN], url[FILENAME_LEN];

  char *long_line;

  int icode, loop, nlines;
  int dos, have_file, itry, ntries;
  int ok;

  struct code *code_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  ntries = 5;
  ok = TRUE;

  /* Allocate memory for very long input lines */
  long_line = (char *) malloc(sizeof(char)*(LONG_LINELEN + 1));

  /* Initialise pointers */
  first_pfam_ptr = *fst_pfam_ptr;
  
  /* Loop over all the UniProt entries */
  for (icode = 0; icode < ncodes; icode++)
    {
      /* Get the UniProt record */
      code_ptr = code_index_ptr[icode];

      /* If UniProt and Pfam sequence lengths don't match, pick up Pfam
         data */
      if (code_ptr->pfam_seq_len == 0 &&
          strcmp(code_ptr->uniprot_acc,"NONE"))
        {
          /* Reset Pfam details */
          code_ptr->pfam_seq_len = 0;
          code_ptr->first_domain_ptr = NULL;

          /* Form the output file name */
          sprintf(file_name,"%s.tsv",code_ptr->uniprot_acc);
          ok = TRUE;

          /* See if we already have this file */
          have_file = empty_check(file_name);
          if (have_file == FALSE)
            {
              printf("Retrieving Pfam data for %s\n",code_ptr->uniprot_acc);
              fflush(stdout);

              /* Form the URL of the required UniProt accession number */
              strcpy(url,"\"");
              strcat(url,file_name_ptr.other_dir[URL_PFAM]);

              /* Append the UniProt accession, and the extra features flag */
              strcat(url,code_ptr->uniprot_acc);
              strcat(url,"?extra_features=true");
              strcat(url,"\"");

              /* Initialise flag */
              ok = FALSE;

              /* Loop while trying to download the file */
              for (itry = 0; itry < ntries && ok == FALSE; itry++)
                {
                  /* Run the wget command to retrieve the Pfam data */
                  printf("Running wget %d ...\n",itry + 1);
                  fflush(stdout);
                  ok
                    = run_wget(file_name_ptr.exe_name[WGET],url,
                               file_name,dos);
                  printf("... done with ok = %s\n",T_F[ok]);
                  fflush(stdout);

                  /* Pause for 5 seconds */
                  sprintf(command,"%s",SLEEP_COMMAND);
                  for (loop = 0; loop < 2; loop++)
                    call_system(command,dos);
                }
            }
          
          /* If OK, interpret the returned XML file */
          if (ok == TRUE)
            {
              /* Read in the Pfam domain assignments */
              npfam = get_pfam_tsv(&first_pfam_ptr,&last_pfam_ptr,
                                   npfam,code_ptr,file_name,
                                   &nlines,long_line);

              /* Delete the downloaded .tsv file */
              sprintf(command,"%s %s",RM_COMMAND,file_name);
              call_system(command,dos);
            }
        }
    }

  /* Return pointers to Pfam records */
  *fst_pfam_ptr = first_pfam_ptr;

  /* Return number of Pfam domains created */
  return(npfam);
}
/***********************************************************************

domain_sort  -  Simple bubble-sort to order the domains by starting
                residue

***********************************************************************/

void domain_sort(struct domain **domain_index_ptr,int *index,
                 int ndomains)
{
  int i, idomain, j;
  int swapped;

  /* Initialise swap flag */
  swapped = TRUE;

  /* Keep looping while any elements are still being swapped */
  while (swapped == TRUE)
    {
      /* Unset the swap flag */
      swapped = FALSE;

      /* Loop through all the domainures */
      for (idomain = 0; idomain < ndomains - 1; idomain++)
        {
          /* Get the two domainure indexes */
          i = index[idomain];
          j = index[idomain + 1];

          /* If the scores are the wrong way round, then need to
             swap them */
          if (domain_index_ptr[i]->start_res >
              domain_index_ptr[j]->start_res)
            {
              /* Swap the corresponding indices */
              index[idomain] = j;
              index[idomain + 1] = i;

              /* Set flag to indicate that a swap took place in
                 this loop */
              swapped = TRUE;
            }
        }
    }
}
/***********************************************************************

order_domains  -  Sort the domains by start residue

***********************************************************************/

int order_domains(struct domain **fst_domain_ptr)
{
  int iorder, ipos, ndomains;
  int *index;

  struct domain *domain_ptr, *first_domain_ptr, *last_domain_ptr;
  struct domain **domain_index_ptr;

  /* Initialise variables */
  ndomains = 0;

  /* Initialise pointers */
  domain_ptr = first_domain_ptr = *fst_domain_ptr;
  last_domain_ptr = NULL;

  /* Loop through domains to count */
  while (domain_ptr != NULL)
    {
      /* Increment domain counter */
      ndomains++;

      /* Get next domain record */
      domain_ptr = domain_ptr->next_domain_ptr;
    }

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(ndomains * sizeof(int));
  domain_index_ptr
    = (struct domain **) malloc(ndomains * sizeof(struct domain *));

  /* Get pointer to first domain record */
  domain_ptr = first_domain_ptr;
  ipos = 0;

  /* Loop through domain records to initialise the sort index */
  while (domain_ptr != NULL && ipos < ndomains)
    {
      /* Add domain to sort index */
      index[ipos] = ipos;
      domain_index_ptr[ipos] = domain_ptr;

      /* Increment domain counter */
      ipos++;

      /* Get next domain record */
      domain_ptr = domain_ptr->next_domain_ptr;
    }

  /* Save number of domains in index */
  ndomains = ipos;

  /* Sort the domains */
  if (ndomains > 1)
    domain_sort(domain_index_ptr,index,ndomains);

  /* Having sorted the domains, rearrange order of the domain-record
     pointers in ascending order of E-value */
  for (iorder = 0; iorder < ndomains; iorder++)
    {
      /* Get the array position of this domain record */
      ipos = index[iorder];

      /* Get the current domain record */
      domain_ptr = domain_index_ptr[ipos];
 
      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_domain_ptr = domain_ptr;

      /* Otherwise, get the previous domain record to point to this one */
      else
        last_domain_ptr->next_domain_ptr = domain_ptr;

      /* Blank out current record's next domain pointer */
      domain_ptr->next_domain_ptr = NULL;

      /* Store current record as last encountered */
      last_domain_ptr = domain_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(domain_index_ptr);

  /* Return current pointer */
  *fst_domain_ptr = first_domain_ptr;

  /* Return number of sorted domains */
  return(ndomains);
}
/***********************************************************************

colour_pfam_domains  -  Check all this PDB entry's Pfam domains to
                        assign a consistent colouring scheme to them

***********************************************************************/

void colour_pfam_domains (struct pdb *pdb_ptr)
{
  char pfam_id[PFAM_LEN];
  char stored_id[NDOMAIN_TYPES][MAX_DOMAINS][PFAM_LEN];

  int icode, idomain, itype, ipos, ndomains[NDOMAIN_TYPES], type;
  int found;

  struct chain *chain_ptr;
  struct code *code_ptr;
  struct domain *chain_domain_ptr, *domain_ptr, *last_domain_ptr;
  struct smatch *smatch_ptr;

  /* Initialise variables */
  for (itype = 0; itype < NDOMAIN_TYPES; itype++)
    {
      ndomains[itype] = 0;
      for (ipos = 0; ipos < MAX_DOMAINS; ipos++)
        stored_id[itype][ipos][0] = '\0';
    }

  /* Get the first of this PDB entry's chains */
  chain_ptr = pdb_ptr->first_chain_ptr;

  /* Loop through all the chain records to identify and colour the
     domains */
  while (chain_ptr != NULL)
    {
      /* Consider only if this is a unique chain  */
      if (chain_ptr->copy_of == '\0' && chain_ptr->ncodes > 0)
        {
          /* Loop over this chain's UniProt codes */
          for (icode = 0; icode < chain_ptr->ncodes; icode++)
            {
              /* Get the corresponding smatch record */
              smatch_ptr = chain_ptr->smatch_ptr[icode];

              /* Initialise first domain pointer */
              smatch_ptr->first_domain_ptr = NULL;
              last_domain_ptr = NULL;

              /* Get the UniProt record */
              code_ptr = chain_ptr->code_ptr[icode];

              /* Sort the domains in order of start residue */
              order_domains(&(code_ptr->first_domain_ptr));

              /* Get the first of this sequence's Pfam domains */
              domain_ptr = code_ptr->first_domain_ptr;

              /* Loop over all the domains to determine their colours */
              while (domain_ptr != NULL)
                {
                  /* Create a copy of this domain record, attached to
                     current chain record */
                  chain_domain_ptr
                    = create_domain_entry(&(smatch_ptr->first_domain_ptr),
                                          &last_domain_ptr,PFAM,
                                          domain_ptr->pfam_ptr,
                                          domain_ptr->pfam_ptr->pfam_id,
                                          domain_ptr->start_res,
                                          domain_ptr->end_res);

                  /* Get this domain's Pfam id */
                  strcpy(pfam_id,domain_ptr->pfam_ptr->pfam_id);

                  /* Update PDB records flag to indicate that we have
                     Pfam ids */
                  pdb_ptr->have_pfam = TRUE;

                  /* Get domain type */
                  type = domain_ptr->pfam_ptr->type;

                  /* If have a valid type see if we already have this
                     domain */
                  if (type != UNKNOWN)
                    {
                      /* Check whether we have encountered this domain */
                      found = FALSE;
                      for (idomain = 0;
                           idomain < ndomains[type] && found == FALSE;
                           idomain++)
                        {
                          /* If domain is in the stored list, then save its
                             position */
                          if (!strcmp(pfam_id,stored_id[type][idomain]))
                            {
                              /* Store as the domain's colour */
                              chain_domain_ptr->colour = idomain;
                              chain_domain_ptr->domain_no = idomain;

                              /* Set flag that domain found */
                              found = TRUE;
                            }
                        }

                      /* If domain not in list, add it */
                      if (found == FALSE && ndomains[type] < MAX_DOMAINS)
                        {
                          /* Get next domain number */
                          idomain = ndomains[type];

                          /* Store domain's Pfam id */
                          strcpy(stored_id[type][idomain],pfam_id);

                          /* Store as colour */
                          chain_domain_ptr->colour = idomain;

                          /* Increment count */
                          ndomains[type]++;

                          /* Store domain number */
                          chain_domain_ptr->domain_no = idomain;
                        }
                    }

                  /* Get the next domain */
                  domain_ptr = domain_ptr->next_domain_ptr;
                }
            }
        }

      /* Get pointer to the next chain record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Set flag to indicate that this PDB entry has been coloured */
  pdb_ptr->coloured = TRUE;
}
/***********************************************************************

run_fasta  -  Run FASTA to compare PDB sequence against corresponding
              UniProt sequence

***********************************************************************/

void run_fasta(char *fasta_exe,char *seq_file,char *uniprot_file,
               int run_no,int dos)
{
  char command[FILENAME_LEN + 1], fasta_out[FILENAME_LEN + 1];

  /* Form name out output file */
  if (run_no == 0)
    strcpy(fasta_out,"fasta.out");
  else
    sprintf(fasta_out,"fasta%3.3d.out",run_no);

  /* Form the parameters to run the FASTA program */
  sprintf(command,"%s -Q -a %s %s > %s",fasta_exe,
          seq_file,uniprot_file,fasta_out);
/* debug */
  printf("Run %d. COMMAND: %s\n",run_no,command);
  fflush(stdout);
/* debug */

  /* Run the FASTA program */
#ifdef NOSYST
#else
  call_system(command,dos);
#endif
}
/***********************************************************************

store_sequence_pair  -  Store the alignment fragments just read in for
                        the two aligned sequences

***********************************************************************/

void store_sequence_pair(char seq_frag[2][FRAG_LEN + 1],
                         char sequence[2][MAX_SEQLEN],char *uniprot_acc)
{
  int i, j, len, max_len, nblanks, seq_len;

  /* Initialise variables */
  max_len = 0;
  seq_len = 0;

  /* Loop over the two sequence fragments (from target and matched
     sequences) to determine which is the longest */
  for (i = 0; i < 2; i++)
    {
      /* Get the length of sequence */
      len = strlen(sequence[i]);

      /* If this is the longest so far, store */
      if (len > seq_len)
        seq_len = len;

      /* Get the length of this fragment */
      len = strlen(seq_frag[i]);

      /* If this is the longest so far, store */
      if (len > max_len)
        max_len = len;
    }

  /* If not too long, add fragments onto current sequences */
  if (seq_len + max_len < MAX_SEQLEN - 1)
    {
      /* Loop over the two fragments again to add to the existing sequence
         fragments */
      for (i = 0; i < 2; i++)
        {
          /* Store the string */
          strcat(sequence[i],seq_frag[i]);

          /* Add required number of blank spaces onto the end */
          nblanks = max_len - strlen(seq_frag[i]);
          for (j = 0; j < nblanks; j++)
            strcat(sequence[i]," ");

          /* Reinitialise sequence fragment */
          seq_frag[i][0] = '\0';
        }
    }

  /* Otherwise, show warning massage */
  else
    {
      printf("*** Warning. Sequence %s exceeds array length: %d\n",
             uniprot_acc,MAX_SEQLEN);
    }
}
/***********************************************************************

read_fasta_out  -  Read in the alignment from the FASTA run

***********************************************************************/

void read_fasta_out(struct chain *chain_ptr,char *fasta_id,
                    char sequence[2][MAX_SEQLEN],char *uniprot_acc,
                    int run_no)
{
  char input_line[LINELEN + 1], fasta_out[FILENAME_LEN + 1];
  char seq_name[LINELEN + 1];
  char alignment_id[2][7];
  char e_value_string[20];
  char chain, pdb_code[5];
  char seq_frag[2][FRAG_LEN + 1];

  int i, len, nalign;
  int nresidues[2], overlap, sw_score;
  int line_type;
  int have_uniprot_seq, have_pdb_seq, markers_next, reading_alignment;
  int first;

  float identity;
  float e_value, z_score;

  FILE *file_ptr;

  /* Initialise variables */
  chain = ' ';
  first = TRUE;
  have_uniprot_seq = FALSE;
  have_pdb_seq = FALSE;
  nalign = 0;
  e_value = 0.0;
  strcpy(e_value_string,"        ");
  identity = 0.0;
  markers_next = FALSE;
  nresidues[0] = 0;
  nresidues[1] = 0;
  overlap = 0;
  pdb_code[0] = '\0';
  reading_alignment = FALSE;
  for (i = 0; i < 2; i++)
    seq_frag[i][0] = '\0';
  sw_score = 0;
  seq_name[0] = '\0';
  z_score = 0.0;

  /* Define sequence identifiers in alignment */
  sprintf(alignment_id[PDB_SEQ],"%s:%c",chain_ptr->pdb_ptr->pdb_code,
          chain_ptr->chain_id);
  strcpy(alignment_id[UNIPROT_SEQ],fasta_id);

  /* Form file name */
  if (run_no == 0)
    strcpy(fasta_out,"fasta.out");
  else
    sprintf(fasta_out,"fasta%3.3d.out",run_no);

  /* Open the fasta.out file */
  if ((file_ptr = fopen(fasta_out,"r")) ==  NULL)
    {
      /* If can't find fasta.out file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",fasta_out);
      return;
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If have start of new alignment or have hit end, finish off
         previous alignment if there is one */
      if (!strncmp(input_line,">>",2) ||
          !strncmp(input_line,"Function",8))
        {
          /* If new alignment, store the UniProt alignment id */
          if (!strncmp(input_line,">>",2))
            {
              strncpy(alignment_id[UNIPROT_SEQ],input_line + 2,6);
              alignment_id[UNIPROT_SEQ][6] = '\0';
            }
          
          /* If have a sequence to store, then do so */
          if (have_uniprot_seq == TRUE || have_pdb_seq == TRUE)
            store_sequence_pair(seq_frag,sequence,uniprot_acc);

          /* Initialise all variables */
          e_value = 0.0;
          strcpy(e_value_string,"        ");
          have_uniprot_seq = FALSE;
          have_pdb_seq = FALSE;
          identity = 0.0;
          nresidues[1] = 0;
          overlap = 0;
          reading_alignment = FALSE;
          for (i = 0; i < 2; i++)
            seq_frag[i][0] = '\0';
          sw_score = 0;
          z_score = 0.0;
        }

      /* Check for line containing the z-score and E value */
      else if (!strncmp(input_line," initn:",7))
        {
          /* Extract all the relevant data from this line */
          //process_initn_line(input_line,&z_score,&e_value,
          //                   e_value_string);
        }

      /* Check for line that signals start of this sequence alignment */
      else if (!strncmp(input_line,"Smith-Waterman score:",21))
        {
          /* Extract Smith-Waterman score, chain identity and overlap */
          // process_sw_line(input_line,&sw_score,&identity,&overlap);

          /* Increment counts of sequences and alignments stored */
          nalign++;

          /* Set flag that can read in this alignment */
          reading_alignment = TRUE;
        }

      /* Otherwise, if reading in the alignment, check which line
         this is */
      else if (reading_alignment == TRUE)
        {
          /* Check for either sequence's identifier */
          line_type = NO_SEQUENCE;
          if (!strncmp(input_line,alignment_id[PDB_SEQ],6))
            {
              /* If have a previous part of the alignment to store,
                 then do so */
              if (have_uniprot_seq == TRUE || have_pdb_seq == TRUE)
                store_sequence_pair(seq_frag,sequence,uniprot_acc);

              /* Set flags */
              line_type = PDB_SEQ;
              markers_next = TRUE;
              have_pdb_seq = TRUE;
              have_uniprot_seq = FALSE;
            }
          else if (!strncmp(input_line,alignment_id[UNIPROT_SEQ],6))
            {
              /* If have a previous part of the alignment to store,
                 then do so */
              if (have_uniprot_seq == TRUE)
                store_sequence_pair(seq_frag,sequence,uniprot_acc);

              /* Set flags */
              line_type = UNIPROT_SEQ;
              have_pdb_seq = FALSE;
              have_uniprot_seq = TRUE;
            }
          else if (markers_next == TRUE)
            {
              line_type = SIMILARITY_MARKERS;
              markers_next = FALSE;
            }

          /* If have a sequence line, extract the sequence itself */
          if (line_type != NO_SEQUENCE &&
              line_type != SIMILARITY_MARKERS && len > 7)
            {
              /* Store this sequence */
              strcpy(seq_frag[line_type],input_line+7);
            }
        }
    }

  /* Close the fasta.out file */
  fclose(file_ptr);
}
/***********************************************************************

read_sequence  -  Read in the protein's sequence from the seqXXX.fasta
                  file

***********************************************************************/

void read_sequence(char *file_name,char sequence[2][MAX_SEQLEN])
{
  char input_line[LINELEN + 1];

  int len, line;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in from file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Increment line count */
          line++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If second line, then save sequence */
          if (line == 2)
            {
              /* Check that sequence not too long */
              if (len < MAX_SEQLEN)
                strcpy(sequence[PDB_SEQ],input_line);
              else
                {
                  strncpy(sequence[PDB_SEQ],input_line,MAX_SEQLEN - 1);
                  sequence[PDB_SEQ][MAX_SEQLEN - 1] = '\0';
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }
}
/***********************************************************************

create_residue_record  -  Create and initialise a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **fst_residue_ptr,
                                      struct residue **lst_residue_ptr,
                                      char chain,char *res_name,
                                      char *res_num,char aa_code,
                                      char sec_struc,int domain)
{
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

  /* If this is the very first, save its pointer */
  if (first_residue_ptr == NULL)
    first_residue_ptr = residue_ptr;

  /* Add link from previous residue to the current one */
  if (last_residue_ptr != NULL)
    last_residue_ptr->next_residue_ptr = residue_ptr;
  last_residue_ptr = residue_ptr;

  /* Get the single-letter amino acid code */
  residue_ptr->aa_code = aa_code;

  /* Initialise the elements of the structure */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;
  residue_ptr->sec_struc = sec_struc;
  residue_ptr->cat_res = FALSE;
  residue_ptr->beta_turn = FALSE;
  residue_ptr->gamma_turn = FALSE;
  residue_ptr->disulphide_id = 0;
  for (type = 0; type < CONTACT_TYPES; type++)
    {
      residue_ptr->nhbonds[type] = 0;
      residue_ptr->ncontacts[type] = 0;
    }
  residue_ptr->domain = domain;
  residue_ptr->nsites = 0;
  residue_ptr->npatterns = 0;
  for (i = 0; i < MAX_SITES; i++)
    residue_ptr->site_name[i] = &null;
  for (i = 0; i < MAX_PATTERNS; i++)
    {
      residue_ptr->pattern_id[i] = 0;
      residue_ptr->pattern_colour[i] = -1;
    }
  residue_ptr->similarity = 0;
  residue_ptr->pred_confidence = 0.0;
  residue_ptr->pfam_no = -1;
  residue_ptr->nligs = 0;
  residue_ptr->range = -1;
  residue_ptr->chain_ptr = NULL;
  residue_ptr->pfam_ptr = NULL;
  residue_ptr->first_liglist_ptr = NULL;
  residue_ptr->last_liglist_ptr = NULL;

  /* Initialise pointer to the next residue */
  residue_ptr->next_residue_ptr = NULL;

  /* Return current pointers */
  *fst_residue_ptr = first_residue_ptr;
  *lst_residue_ptr = last_residue_ptr;

  /* Return new residue pointer */
  return(residue_ptr);
}
/***********************************************************************

create_ligand_entry  -  Create a ligand entry

***********************************************************************/

struct ligand *create_ligand_entry(struct ligand **fst_ligand_ptr,
                                   struct ligand **lst_ligand_ptr,
                                   char *ligand_name,int type,int number)
{
  struct ligand *ligand_ptr;
  struct ligand *first_ligand_ptr, *last_ligand_ptr;

  /* Initialise ligand pointers */
  ligand_ptr = NULL;
  first_ligand_ptr = *fst_ligand_ptr;
  last_ligand_ptr = *lst_ligand_ptr;

  /* Create ligand entry */
  ligand_ptr = (struct ligand *) malloc(sizeof(struct ligand));
  if (ligand_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct ligand\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  store_string(&(ligand_ptr->ligand_name),ligand_name);
  ligand_ptr->type = type;
  ligand_ptr->number = number;

  /* Initialise remaining fields */
  ligand_ptr->description = &null;
  ligand_ptr->nhets = 0;
  ligand_ptr->first_het_ptr = NULL;

  /* Initialise pointer to next entry */
  ligand_ptr->next_ligand_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_ligand_ptr == NULL)
    first_ligand_ptr = ligand_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_ligand_ptr->next_ligand_ptr = ligand_ptr;
                      
  /* Save the current residue's pointer */
  last_ligand_ptr = ligand_ptr;

  /* Return current pointers */
  *fst_ligand_ptr = first_ligand_ptr;
  *lst_ligand_ptr = last_ligand_ptr;
  return(ligand_ptr);
}
/***********************************************************************

create_het_record  -  Create a new het record

***********************************************************************/

struct het *create_het_record(struct het **fst_het_ptr,
                              struct het **lst_het_ptr,
                              char *name,char *desc)
{
  struct het *het_ptr;
  struct het *first_het_ptr, *last_het_ptr;

  /* Initialise het pointers */
  last_het_ptr = *lst_het_ptr;
  first_het_ptr = *fst_het_ptr;

  /* Create a new het entry */
  het_ptr = (struct het*)malloc(sizeof(struct het));
  if (het_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct het\n");
      printf("***        Program resdata terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(het_ptr->name,name);
  store_string(&(het_ptr->desc),desc);

  /* Initialise pointer to next het record */
  het_ptr->next_het_ptr = NULL;

  /* If this is the first entry stored, then update pointer to
     head of linked list */
  if (first_het_ptr == NULL)
    first_het_ptr = het_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_het_ptr->next_het_ptr = het_ptr;
                      
  /* Save the current het's pointer */
  last_het_ptr = het_ptr;

  /* Return current pointers */
  *fst_het_ptr = first_het_ptr;
  *lst_het_ptr = last_het_ptr;

  return(het_ptr);
}
/***********************************************************************

get_resdata  -  Read in the given structure from the PDB file

***********************************************************************/

int get_resdata(char *pdb_code,char chain,struct chain *chain_ptr,
                int icode,char *pdbsum_dir,
                struct residue **fst_residue_ptr,int *nx)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char number_string[20];
  char aa_code, inchain, res_name[4], res_num[6], sec_struc;
  char description[LINELEN + 1], ligand_name[LINELEN + 1];

  static char subscript[] = { 'D', 'L', 'M' };

  char *string_ptr, *token[MAX_TOKENS];

  int ctype, domain, i, ipattn, last_domain, len, nox, nresid;
  int itoken, ntokens, number, type;
  int done, got_chain, reading_sequence, wanted;

  struct het *het_ptr, *last_het_ptr;
  struct ligand *ligand_ptr;
  struct residue *residue_ptr, *first_residue_ptr, *last_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  last_domain = -99;
  nresid = 0;
  done = FALSE;
  got_chain = FALSE;
  *nx = nox = 0;
  wanted = FALSE;
  reading_sequence = FALSE;

  /* Initialise pointers */
  residue_ptr = NULL;
  first_residue_ptr = *fst_residue_ptr = NULL;
  last_residue_ptr = NULL;
  ligand_ptr = NULL;

  /* Form name of the resdata.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"resdata.dat");

  /* Open the resdata.dat file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find resdata.dat file, then return empty-handed */
      printf("*** Warning. Unable to find resdata.dat file for %s\n",
             pdb_code);
      printf("***          %s\n",file_name);
      return(0);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,LINELEN);

      /* Check for line signalling start of residue sequence */
      if (!strncmp(input_line,"----------",10))
        reading_sequence = TRUE;

      /* If reading sequence, interpret line */
      else if (reading_sequence == TRUE)
        {
          /* If first character not a space, then have a residue */
          if (input_line[0] != ' ')
            {
              /* Get the chain identifier */
              inchain = input_line[12];

              /* If chain is wanted, process the record */
              if (inchain == chain)
                {
                  /* Set flag that we have the chain we're after */
                  wanted = TRUE;
                  got_chain = TRUE;

                  /* Extract the residue information on the line */
                  aa_code = input_line[0];
                  strncpy(res_name,input_line+2,3);
                  res_name[3] = '\0';
                  strncpy(res_num,input_line+6,5);
                  res_num[5] = '\0';
                  inchain = input_line[12];

                  /* Get the secondary structure assignment */
                  sec_struc = input_line[14];
                  
                  /* Get the domain number */
                  strcpy(number_string,input_line+16);
                  domain = atoi(number_string);

                  /* Create a residue record and store data just read in */
                  residue_ptr
                    = create_residue_record(&first_residue_ptr,
                                            &last_residue_ptr,inchain,
                                            res_name,res_num,aa_code,
                                            sec_struc,domain);

                  /* Add chain record */
                  residue_ptr->chain_ptr = chain_ptr;

                  /* Increment count of residues stored */
                  nresid++;

                  /* If aa code for residue is 'X', increment count */
                  if (aa_code == 'X')
                    nox++;
                }

              /* Otherwise, check whether we're done */
              else
                {
                  /* Set flag that chain not wanted */
                  wanted = FALSE;

                  /* If we already have the chain we were after, we're
                     done */
                  if (got_chain == TRUE)
                    done = TRUE;
                }
            }

          /* Otherwise, if item belongs to a chain of interest, store
             for last-read residue */
          else if (wanted == TRUE)
            {
              /* Disulphide bond */
              if (!strncmp(input_line," DISUL",6))
                {
                  /* Extract the disulphide identifier */
                  strcpy(number_string,input_line+6);
                  residue_ptr->disulphide_id = atoi(number_string);
                }

              /* Beta turn */
              else if (!strncmp(input_line," BTURN",6))
                residue_ptr->beta_turn = TRUE;

              /* Gamma turn */
              else if (!strncmp(input_line," GTURN",6))
                residue_ptr->gamma_turn = TRUE;

              /* Active site */
              else if (!strncmp(input_line," SITE",5))
                {
                  /* Store name of active site, if have room */
                  if (residue_ptr->nsites < MAX_SITES)
                    {
                      /* Store the name of the active site */
                      i = residue_ptr->nsites;
                      store_string(&(residue_ptr->site_name[i]),
                                   input_line+6);

                      /* Increment count of sites that this
                         residue belongs to  */
                      residue_ptr->nsites++;
                    }
                }

              /* PROSITE pattern */
              else if (!strncmp(input_line," PROSITE",8))
                {
                  /* If have room to store, do so */
                  if (residue_ptr->npatterns < MAX_PATTERNS)
                    {
                      /* Get the prosite pattern number and store */
                      strncpy(number_string,input_line+8,3);
                      number_string[3] = '\0';
                      ipattn = residue_ptr->npatterns;
                      residue_ptr->pattern_id[ipattn]
                        = atoi(number_string);

                      /* Get the residue colour */
                      strcpy(number_string,input_line+11);
                      residue_ptr->pattern_colour[ipattn]
                        = atoi(number_string);

                      /* Increment count of patterns that this
                         residue belongs to  */
                      residue_ptr->npatterns++;
                    }
                      
                }

              /* H-bonds */
              else if (!strncmp(input_line," HB_",4))
                {
                  /* Get the contact type */
                  ctype = -1;
                  for (i = 0; i < CONTACT_TYPES; i++)
                    if (input_line[4] == subscript[i])
                      ctype = i;

                  /* If valid, store the count of H-bonds */
                  if (ctype != -1)
                    {
                      strcpy(number_string,input_line+5);
                      residue_ptr->nhbonds[ctype]
                        = atoi(number_string);
                    }
                }

              /* Non-bonded contacts */
              else if (!strncmp(input_line," NB_",4))
                {
                  /* Get the contact type */
                  ctype = -1;
                  for (i = 0; i < CONTACT_TYPES; i++)
                    if (input_line[4] == subscript[i])
                      ctype = i;

                  /* If valid, store the count of H-bonds */
                  if (ctype != -1)
                    {
                      strcpy(number_string,input_line+5);
                      residue_ptr->ncontacts[ctype]
                        = atoi(number_string);
                    }
                }
            }
        }

      /* If ligand or metal, pick up the relevant data */
      else if (!strncmp(input_line,"LIGAND ",7) ||
               !strncmp(input_line,"METAL ",6))
        {
          /* Split into tokens */
          ntokens = get_tokens(input_line,token,MAX_TOKENS,' ');
 
          /* If have enough tokens, then process */
          if (ntokens > 2)
            {
              /* Get the ligand number */
              number = atoi(token[1]);

              /* Get the ligand/metal sequence */
              strcpy(ligand_name,token[2]);

              /* Initialise ligand description */
              description[0] = '\0';

              /* Define the ligand type */
              if (!strcmp(token[0],"LIGAND"))
                type = LIGAND;

              /* Otherwise, get the metal description */
              else
                {
                  /* Define type */
                  type = METAL;

                  /* Get the metal name */
                  description[0] = '\0';
                  for (itoken = 3; itoken < ntokens; itoken++)
                    {
                      /* If not first token, add space */
                      if (description[0] != '\0')
                        {
                          strcat(description," ");
                          strcat(description,token[itoken]);
                        }

                      /* Otherwise, copy current token */
                      else
                        strcpy(description,token[itoken]);
                    }

                  /* Change to lower case */
                  len = strlen(description);
                  to_lower(description,len);
                }

              /* Create a ligand record */
              ligand_ptr
                = create_ligand_entry(&(chain_ptr->first_ligand_ptr),
                                      &(chain_ptr->last_ligand_ptr),
                                      ligand_name,type,number);

              /* If this is a metal, then save the description */
              if (type == METAL)
                store_string(&(ligand_ptr->description),description);

              /* Increment ligand count */
              chain_ptr->nligands++;

              /* Initialise last Het Group record */
              last_het_ptr = NULL;
            }

          /* Free the tokens */
          free_tokens(token,ntokens);
        }

      /* If this is a Het Group definition, then save */
      else if (!strncmp(input_line," HET_GRP: ",10) &&
               ligand_ptr != NULL)
        {
          /* Split into tokens */
          ntokens = get_tokens(input_line + 1,token,MAX_TOKENS,' ');
 
          /* If have enough tokens, then process */
          if (ntokens > 2)
            {
              /* Get the Het Group name */
              description[0] = '\0';
              for (itoken = 2; itoken < ntokens; itoken++)
                {
                  /* If not first token, add space */
                  if (description[0] != '\0')
                    {
                      strcat(description," ");
                      strcat(description,token[itoken]);
                    }

                  /* Otherwise, copy current token */
                  else
                    strcpy(description,token[itoken]);
                }

              /* Create and store a Het Group record */
              het_ptr
                = create_het_record(&(ligand_ptr->first_het_ptr),
                                    &last_het_ptr,token[1],description);

              /* Increment count of Het Groups */
              ligand_ptr->nhets++;
            }

          /* Free the tokens */
          free_tokens(token,ntokens);
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return pointer to first residue */
  *fst_residue_ptr = first_residue_ptr;

  /* Return number of residues labelled as 'X' */
  *nx = nox;

  /* Return number of residues read in */
  return(nresid);
}
/***********************************************************************

match_to_residues  -  Match the sequence to the residue records

***********************************************************************/

int match_to_residues(struct chain *chain_ptr,int icode,
                      struct residue *first_residue_ptr,int nresid,
                      struct c_score *first_c_score_ptr,int nscores,
                      int nox,char sequence[2][MAX_SEQLEN],
                      char *sec_struc,char *res_range,int *cath_domain,
                      char pdb_resno[2][MAX_SEQLEN],char *uniprot_acc)
{
  char aa_code, number_string[20], res_num[6];

  int align_len, i, ipos, ires, last_pos, jpos, seq_len, strand_length;
  int end_res, num_len, num_start, pfam_no, res_no, start_res;
  int done, ok, warned;

  struct c_score *c_score_ptr;
  struct domain *domain_ptr;
  struct pfam *pfam_ptr;
  struct residue *residue_ptr;
  struct smatch *smatch_ptr;

  /* Initialise variables */
  ires = 0;
  last_pos = -1;
  ok = TRUE;
  pfam_ptr = NULL;
  seq_len = 0;
  start_res = end_res = -1;
  warned = FALSE;

  /* Get the alignment and sequence lengths */
  align_len = strlen(sequence[PDB_SEQ]);
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* If residue not blank and not a hyphen, increment seq len */
      aa_code = sequence[PDB_SEQ][ipos];
      if (aa_code != ' ' && aa_code != '-')
        seq_len++;

      /* Initialise secondary structure and CATH domain */
      sec_struc[ipos] = ' ';
      res_range[ipos] = ' ';
      cath_domain[ipos] = -99;
      if (aa_code != ' ' && aa_code != '-')
        pdb_resno[HYPHENS][ipos] = '-';
      else
        pdb_resno[HYPHENS][ipos] = ' ';
      pdb_resno[NUMBERS][ipos] = ' ';
    }
  pdb_resno[HYPHENS][align_len] = '\0';
  pdb_resno[NUMBERS][align_len] = '\0';

  /* Check whether number of stored residues corresponds to the
     known sequence length */
  if (nresid < seq_len || nresid > seq_len + nox)
    {
      /* Show warning message */
      printf("*** Warning. Mismatch in sequence length for PDB code");
      printf(" %s, chain '%c'\n",chain_ptr->pdb_ptr->pdb_code,
             chain_ptr->chain_id);
      printf("             Seq length = %d, no. resid in resdata.dat "
             "= %d\n",seq_len,nresid);

      /* Flag this sequence so that it isn't used */
      ok = FALSE;
      return(ok);
    }

  /* Store number of residues in this chain */
  chain_ptr->nresid = nresid;

  /* Get first residue */
  residue_ptr = first_residue_ptr;
  c_score_ptr = first_c_score_ptr;

  /* Loop through the sequence to match up the residues */
  for (ipos = 0; ipos < align_len && residue_ptr != NULL; ipos++)
    {
      /* Determine the residue number in the UniProt sequence */
      if (strcmp(uniprot_acc,"NONE"))
        {
          aa_code = sequence[UNIPROT_SEQ][ipos];
          if (aa_code != ' ' && aa_code != '-')
            ires++;
        }

      /* Get this amino acid code */
      aa_code = sequence[PDB_SEQ][ipos];

      /* If not blank and not a hyphen, then process */
      if (aa_code != ' ' && aa_code != '-')
        {
          /* If no domain, then see if this residue belongs to a Pfam
             domain */
          if (ires < start_res || ires > end_res || pfam_ptr == NULL)
            {
              /* Initialise Pfam pointer */
              pfam_ptr = NULL;
              pfam_no = 0;

              /* Get the first UniProt match to this chain */
              smatch_ptr = chain_ptr->smatch_ptr[icode];
              if (smatch_ptr != NULL)
                domain_ptr = smatch_ptr->first_domain_ptr;
              else
                domain_ptr = NULL;

              /* Initialise flag */
              done = FALSE;

              /* Loop through the domains to see if current residue
                 belongs to one of the Pfam domains */
              while (domain_ptr != NULL && done == FALSE)
                {
                  /* Get the corresponding Pfam pointer */
                  pfam_ptr = domain_ptr->pfam_ptr;

                  /* If this is a Pfam domain determine whether
                     the current residue belongs to it */
                  if (pfam_ptr != NULL)
                    {
                      /* Increment Pfam domain count */
                      pfam_no++;

                      /* If residue within the domain range,
                         then can assign it to this domain */
                      if (ires >= domain_ptr->start_res &&
                          ires <= domain_ptr->end_res)
                        {
                          /* Attach Pfam pointer to this residue */
                          residue_ptr->pfam_ptr = pfam_ptr;
                          residue_ptr->pfam_no = pfam_no;

                          /* Set flag that found */
                          done = TRUE;
                        }
                    }

                  /* Get the next domain record */
                  domain_ptr = domain_ptr->next_domain_ptr;
                }
            }

          /* If have an AF score, copy across the range */
          if (c_score_ptr != NULL &&
              !strcmp(c_score_ptr->res_name,residue_ptr->res_name) &&
              !strcmp(c_score_ptr->res_num,residue_ptr->res_num) &&
              c_score_ptr->chain == residue_ptr->chain)
            {
              residue_ptr->range = c_score_ptr->range;
            }

          /* Check for a match with the residue record */
          if (aa_code == residue_ptr->aa_code ||
              residue_ptr->aa_code == 'X')
            {
              /* Store the secondary structure and CATH domain
                 assignments */
              if (ipos < MAX_SEQLEN - 1)
                {
                  sec_struc[ipos] = residue_ptr->sec_struc;
                  cath_domain[ipos] = residue_ptr->domain;

                  /* Repeat for AF score, if we have one */
                  if (residue_ptr->range > -1)
                    {
                      sprintf(number_string,"%d",residue_ptr->range);
                      res_range[ipos] = number_string[0];
                    }

                  /* Check whether to store this residue number in
                     the residue numbers line */
                  strcpy(res_num,residue_ptr->res_num);
                  res_num[4] = '\0';
                  res_no = atoi(res_num);

                  /* If number divisble by 10, consider writing it out */
                  if (res_no % 10 == 0)
                    {
                      /* Get residue number start point */
                      if (res_num[2] == ' ')
                        num_start = 3;
                      else if (res_num[1] == ' ')
                        num_start = 2;
                      else if (res_num[0] == ' ')
                        num_start = 1;
                      else
                        num_start = 0;

                      /* Get the length of the residue number */
                      num_len = 4 - num_start;
                      if (residue_ptr->res_num[4] != ' ')
                        num_len++;

                      /* Get the residue number */
                      strncpy(number_string,
                              residue_ptr->res_num+num_start,num_len);
                      number_string[num_len] = '\0';

                      /* Check if there is enough space to write residue
                         number out */
                      if (ipos - last_pos > num_len)
                        {
                          /* Transfer residue number across */
                          for (i = 0; i < num_len; i++)
                            {
                              /* Get position in residue number string */
                              jpos = i + ipos - num_len + 1;

                              /* Copy number across */
                              pdb_resno[NUMBERS][jpos]
                                = residue_ptr->res_num[num_start + i];
                            }
                          
                          /* Mark the hyphens line */
                          pdb_resno[HYPHENS][ipos] = '+';

                          /* Store the last number position */
                          last_pos = ipos;
                        }
                    }
                }

              /* If off end of sequence, show error message */
              else
                {
                  /* Show warning if we haven't already done so */
                  if (warned == FALSE)
                    printf("*** Warning. Sequence too long: %d ...\n",
                           seq_len);
                  warned = TRUE;
                }

              /* Get the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              if (c_score_ptr != NULL)
                c_score_ptr = c_score_ptr->next_c_score_ptr;
            }

          /* If don't have a match, have a problem */
          else
            {
              /* Show warning message */
              printf("*** Warning. Sequence mismatch for PDB code");
              printf(" %s, chain '%c'\n",chain_ptr->pdb_ptr->pdb_code,
                     chain_ptr->chain_id);
              printf("           %c  vs  %c  (%s %s %c)\n",aa_code,
                     residue_ptr->aa_code,residue_ptr->res_name,
                     residue_ptr->res_num,residue_ptr->chain);

              /* Mark chain as not OK */
              ok = FALSE;

              /* End here */
              residue_ptr = NULL;
            }
        }
    }

  /* If all OK store the secondary structure assignments */
  if (ok == TRUE)
    {
      /* Close off secondary structures and CATH domains strings */
      sec_struc[align_len] = '\0';
      cath_domain[align_len] = '\0';

      /* Check through the secondary structure assignments to extend
         strands where start or end in beta-turns */
      for (ipos = 0; ipos < align_len; ipos++)
        {
          /* If turn, check secondary structure either side */
          if (sec_struc[ipos] == 'T')
            {
              /* If have strand in previous residue, make this one strand */
              if (ipos > 0 && sec_struc[ipos - 1] == 'E')
                sec_struc[ipos] = 'E';

              /* Ditto for following residue */
              if (ipos < align_len - 1 && sec_struc[ipos + 1] == 'E')
                sec_struc[ipos] = 'E';
            }
        }

      /* Initialise strand length variable */
      strand_length = 0;

      /* Loop through again to delete any very short strands of 1 or 2
         residues in length */
      for (ipos = 0; ipos < align_len; ipos++)
        {
          /* If strand, then increment count */
          if (sec_struc[ipos] == 'E')
            strand_length++;

          /* Otherwise, if have just finished strand, check its length */
          else if (strand_length > 0)
            {
              /* If strand too short, then eliminate it */
              if (strand_length < 3)
                {
                  /* Loop over previous strand position and convert to
                     coil */
                  for (jpos = ipos - strand_length; jpos < ipos; jpos++)
                    sec_struc[jpos] = '.';
                }

              /* Re-initialise strand length variable */
              strand_length = 0;
            }
        }
    }

  /* Return whether match was successful */
  return(ok);
}
/***********************************************************************

find_residue  -  Find the given residue entry

***********************************************************************/

struct residue *find_residue(struct residue *first_residue_ptr,
                               char *res_name,char *res_num,char chain)
{
  struct residue *residue_ptr, *found_residue_ptr;

  /* Initialise variables */
  found_residue_ptr = NULL;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop over all the residues to find the one of interest */
  while (residue_ptr != NULL && found_residue_ptr == NULL)
    {
      /* If this is the right residue, store pointer */
      if (chain == residue_ptr->chain &&
          !strcmp(res_num,residue_ptr->res_num) &&
          !strcmp(res_name,residue_ptr->res_name))
        found_residue_ptr = residue_ptr;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the residueing residue pointer */
  return(found_residue_ptr);
}
/***********************************************************************

find_liglist_entry  -  Find the given liglist entry

***********************************************************************/

struct liglist *find_liglist_entry(struct liglist *first_liglist_ptr,
                                    struct r_rasdef *r_rasdef_ptr)
{
  struct liglist *liglist_ptr, *found_liglist_ptr;

  /* Initialise variables */
  found_liglist_ptr = NULL;
  if (r_rasdef_ptr == NULL)
    return(found_liglist_ptr);

  /* Get pointer to the first liglist */
  liglist_ptr = first_liglist_ptr;

  /* Loop over all the liglists to find the one just read in */
  while (liglist_ptr != NULL && found_liglist_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (liglist_ptr->r_rasdef_ptr == r_rasdef_ptr)
        found_liglist_ptr = liglist_ptr;

      /* Get pointer to the next liglist */
      liglist_ptr = liglist_ptr->next_liglist_ptr;
    }

  /* Return the matching liglist pointer */
  return(found_liglist_ptr);
}
/***********************************************************************

create_liglist_entry  -  Create a liglist entry

***********************************************************************/

struct liglist *create_liglist_entry(struct liglist **fst_liglist_ptr,
                                     struct liglist **lst_liglist_ptr,
                                     struct r_rasdef *r_rasdef_ptr)
{
  struct liglist *liglist_ptr;
  struct liglist *first_liglist_ptr, *last_liglist_ptr;

  /* Initialise liglist pointers */
  liglist_ptr = NULL;
  first_liglist_ptr = *fst_liglist_ptr;
  last_liglist_ptr = *lst_liglist_ptr;

  /* Create liglist entry */
  liglist_ptr = (struct liglist *) malloc(sizeof(struct liglist));
  if (liglist_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct liglist\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  liglist_ptr->r_rasdef_ptr = r_rasdef_ptr;

  /* Initialise pointer to next entry */
  liglist_ptr->next_liglist_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_liglist_ptr == NULL)
    first_liglist_ptr = liglist_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_liglist_ptr->next_liglist_ptr = liglist_ptr;
                      
  /* Save the current residue's pointer */
  last_liglist_ptr = liglist_ptr;

  /* Return current pointers */
  *fst_liglist_ptr = first_liglist_ptr;
  *lst_liglist_ptr = last_liglist_ptr;
  return(liglist_ptr);
}
/***********************************************************************

add_ligand  -  Find the given pfamlig entry

***********************************************************************/

void add_ligand(struct residue *residue_ptr,struct r_rasdef *r_rasdef_ptr)
{
  struct chain *chain_ptr;
  struct liglist *liglist_ptr;

  /* See if this rasdef entry already linked to this residue */
  liglist_ptr
    = find_liglist_entry(residue_ptr->first_liglist_ptr,r_rasdef_ptr);

  /* If not found, then create */
  if (liglist_ptr == NULL)
    {
      /* Add rasdef entry to current residue record */
      liglist_ptr
        = create_liglist_entry(&(residue_ptr->first_liglist_ptr),
                               &(residue_ptr->last_liglist_ptr),
                               r_rasdef_ptr);

      /* Increment count of this residue's ligands */
      residue_ptr->nligs++;

      /* Get the corresponding chain record */
      chain_ptr = residue_ptr->chain_ptr;

      /* If valid, add resdef entry to chain */
      if (chain_ptr != NULL)
        {
          /* See if this rasdef entry already linked to this chain */
          liglist_ptr
            = find_liglist_entry(chain_ptr->first_liglist_ptr,
                                 r_rasdef_ptr);

          /* If not found, then create */
          if (liglist_ptr == NULL)
            {
              /* Add rasdef entry to current chain record */
              liglist_ptr
                = create_liglist_entry(&(chain_ptr->first_liglist_ptr),
                                       &(chain_ptr->last_liglist_ptr),
                                       r_rasdef_ptr);

              /* Increment count of this chain's ligands */
              chain_ptr->nligs++;
            }
        }
    }
}
/***********************************************************************

find_pfamlig_record  -  Find the given pfamlig entry

***********************************************************************/

struct pfamlig *find_pfamlig_record(struct pfamlig *first_pfamlig_ptr,
                                    char *ligand_name,char *pdb_code,
                                    char chain,char *ligand_id,
                                    struct pfam *pfam_ptr,
                                    struct code *code_ptr,int domain_no)
{
  struct pfamlig *pfamlig_ptr, *found_pfamlig_ptr;

  /* Initialise variables */
  found_pfamlig_ptr = NULL;
  if (ligand_name[0] == '\0')
    return(found_pfamlig_ptr);

  /* Get pointer to the first pfamlig */
  pfamlig_ptr = first_pfamlig_ptr;

  /* Loop over all the pfamligs to find the one just read in */
  while (pfamlig_ptr != NULL && found_pfamlig_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (pfamlig_ptr->domain_no == domain_no &&
          pfamlig_ptr->pfam_ptr == pfam_ptr &&
          pfamlig_ptr->code_ptr == code_ptr &&
          !strcmp(pfamlig_ptr->pdb_code,pdb_code) &&
          pfamlig_ptr->chain == chain &&
          !strcmp(pfamlig_ptr->ligand_id,ligand_id) &&
          !strcmp(pfamlig_ptr->ligand_name,ligand_name))
        found_pfamlig_ptr = pfamlig_ptr;

      /* Get pointer to the next pfamlig */
      pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;
    }

  /* Return the matching pfamlig pointer */
  return(found_pfamlig_ptr);
}
/***********************************************************************

create_pfamlig_record  -  Create a new pfamlig record

***********************************************************************/

struct pfamlig *create_pfamlig_record(struct pfamlig **fst_pfamlig_ptr,
                                      struct pfamlig **lst_pfamlig_ptr,
                                      char *ligand_name,char *ligand_id,
                                      char *pdb_code,char chain,
                                      int seq_no,struct pfam *pfam_ptr,
                                      struct code *code_ptr,
                                      int domain_no,int type)
{
  struct pfamlig *pfamlig_ptr;
  struct pfamlig *first_pfamlig_ptr, *last_pfamlig_ptr;

  /* Initialise pfamlig pointers */
  pfamlig_ptr = NULL;
  first_pfamlig_ptr = *fst_pfamlig_ptr;
  last_pfamlig_ptr = *lst_pfamlig_ptr;

  /* Create pfamlig entry */
  pfamlig_ptr = (struct pfamlig*)malloc(sizeof(struct pfamlig));
  if (pfamlig_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct pfamlig\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  store_string(&(pfamlig_ptr->ligand_name),ligand_name);
  strcpy(pfamlig_ptr->ligand_id,ligand_id);
  strcpy(pfamlig_ptr->pdb_code,pdb_code);
  pfamlig_ptr->chain = chain;
  pfamlig_ptr->seq_no = seq_no;
  pfamlig_ptr->domain_no = domain_no;
  pfamlig_ptr->type = type;
  pfamlig_ptr->pfam_ptr = pfam_ptr;
  pfamlig_ptr->code_ptr = code_ptr;

  /* Initialise other fields */
  pfamlig_ptr->domain_list = &null;
  pfamlig_ptr->ncontacts = 0;
  pfamlig_ptr->max_contacts = 0;
  pfamlig_ptr->done = FALSE;
  pfamlig_ptr->first_plist_ptr = NULL;

  /* Initialise pointer to next pfamlig record */
  pfamlig_ptr->next_pfamlig_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of pfamliged pfamlig */
  if (first_pfamlig_ptr == NULL)
    first_pfamlig_ptr = pfamlig_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_pfamlig_ptr->next_pfamlig_ptr = pfamlig_ptr;
                      
  /* Save the current pfamlig's pointer */
  last_pfamlig_ptr = pfamlig_ptr;

  /* Return current pointers */
  *fst_pfamlig_ptr = first_pfamlig_ptr;
  *lst_pfamlig_ptr = last_pfamlig_ptr;
  return(pfamlig_ptr);
}
/***********************************************************************

identify_ligands  -  Identify which interact with the current chain

***********************************************************************/

int identify_ligands(char *pdb_code,struct chain *chain_ptr,
                     struct code *code_ptr,
                     struct r_rasdef *first_r_rasdef_ptr,
                     struct residue *first_residue_ptr,
                     struct pfamlig **fst_pfamlig_ptr,
                     struct pfamlig **lst_pfamlig_ptr,
                     struct pfamlig **strt_pfamlig_ptr,int *npflig,
                     char *pdbsum_dir)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, inchain, res_name[4], res_num[6];
  char other_inchain, other_res_name[4], other_res_num[6];
  char ligand_id[7], other_molecule_type, residue_name[15];
  char res_names[RESNAMES_LEN + 1];

  char *string_ptr;

  int chain_no, ipos, len;
  int ndefs, nligand, npfamlig, nres, nstored;
  int rlen, type;
  int wanted;

  struct pfam *pfam_ptr;
  struct pfamlig *pfamlig_ptr, *first_pfamlig_ptr, *last_pfamlig_ptr;
  struct pfamlig *start_pfamlig_ptr;
  struct r_rasdef *r_rasdef_ptr;
  struct r_rasdef *stored_r_rasdef_ptr[MAX_RESNAMES];
  struct residue *residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  chain = chain_ptr->chain_id;
  chain_no = chain_ptr->chain_no;
  len = 0;
  ndefs = 0;
  nligand = 0;
  npfamlig = *npflig;
  nres = 0;
  nstored = 0;
  res_names[0] = '\0';

  /* Initialise pointers */
  first_pfamlig_ptr = *fst_pfamlig_ptr;
  last_pfamlig_ptr = *lst_pfamlig_ptr;
  start_pfamlig_ptr = NULL;

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through all the definition records to store the ligand, and
     metal residues, plus DNA chains for quick scanning */
  while (r_rasdef_ptr != NULL && nstored < MAX_RESNAMES)
    {
      /* Initialise flag */
      wanted = FALSE;

      /* If this is a nucleic acid, store its chain */
      if (r_rasdef_ptr->type == NUCLEIC_ACID && len < LINELEN - 3)
        {
          /* Form chain name */
          sprintf(residue_name,"|%c|........",r_rasdef_ptr->chain);

          /* Set flag */
          wanted = TRUE;
        }

      /* Check for a match to a ligand */
      else if (r_rasdef_ptr->type == LIGAND && len < LINELEN - 11)
        {
          /* Form full residue name */
          sprintf(residue_name,"|%s%s%c|",r_rasdef_ptr->res_name,
                  r_rasdef_ptr->res_num,r_rasdef_ptr->chain);

          /* Set flag */
          wanted = TRUE;
        }

      /* Check for a match to a metal */
      else if (r_rasdef_ptr->type == METAL && len < LINELEN - 5)
        {
          /* Form residue name */
          sprintf(residue_name,"|%s|......",r_rasdef_ptr->res_name);

          /* Set flag */
          if (r_rasdef_ptr->ligplot_id[0] != '\0')
            wanted = TRUE;
        }

      /* If this definition wanted, store residue name and pointer */
      if (wanted == TRUE)
        {
          /* Add to end of residue names list */
          strcat(res_names,residue_name);

          /* Increase length of names list */
          len = len + 11;

          /* Store the pointer */
          stored_r_rasdef_ptr[nstored] = r_rasdef_ptr;

          /* Increment count of stored pointers */
          nstored++;
        }

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

  /* Form the name of the grow.out */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"grow.out");

  /* Open the GROW file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Initialise flag */
          rlen = 0;
          type = OTHER;
          r_rasdef_ptr = NULL;
          string_ptr = NULL;

          /* Truncate the string to strip off the newline */
          string_truncate(input_line,LINELEN);

          /* If this is an interaction to the current chain, check
             which molecule the interaction is with */
          if (input_line[0] != '#' && input_line[6] == chain)
            {
              /* Get the interacting molecule */
              other_molecule_type = input_line[4];

              /* Get the residue name of the interacting molecule */
              strncpy(other_res_name,input_line+44,3);
              other_res_name[3] = '\0';
              strncpy(other_res_num,input_line+49,5);
              other_res_num[5] = '\0';
              other_inchain = input_line[54];

              /* Form unique ligand id */
              ligand_id[0] = other_inchain;
              ligand_id[1] = '\0';
              strcat(ligand_id,other_res_num);

              /* If this is a ligand, form full residue name */
              if (other_molecule_type == 'L')
                {
                  /* Get type */
                  type = LIGAND;
                  rlen = 10;

                  /* Form search string */
                  sprintf(residue_name,"|%s%s%c|",other_res_name,
                          other_res_num,other_inchain);
                }

              /* Ditto for metal */
              else if (other_molecule_type == 'M')
                {
                  /* Get type */
                  type = METAL;
                  rlen = 5;

                  /* Form search string */
                  sprintf(residue_name,"|%s|",other_res_name);
                }

              /* Ditto for DNA */
              else if (other_molecule_type == 'N')
                {
                  /* Get type */
                  type = NUCLEIC_ACID;
                  rlen = 3;

                  /* Form search string */
                  sprintf(residue_name,"|%c|",other_inchain);
                }

              /* Check whether this residue is in our list */
              if (type != OTHER)
                string_ptr = strstr(res_names,residue_name);

              /* If in list, then need to flag the corresponding
                 ligand residue */
              if (string_ptr != NULL)
                {
                  /* Get the corresponding rasdef pointer */
                  ipos = (string_ptr - res_names) / 11;
                  if (ipos < nstored)
                    r_rasdef_ptr = stored_r_rasdef_ptr[ipos];
                }
            }

          /* If need to flag interacting molecule, do so */
          if (r_rasdef_ptr != NULL)
            {
              /* Set the wanted flag for this residue */
              if (r_rasdef_ptr->wanted == FALSE)
                nligand++;
              r_rasdef_ptr->wanted = TRUE;

              /* Get the residue name, number and chain */
              strncpy(res_name,input_line+14,3);
              res_name[3] = '\0';
              strncpy(res_num,input_line+9,5);
              res_num[5] = '\0';
              inchain = input_line[6];

              /* Identify the interacting residue */
              residue_ptr
                = find_residue(first_residue_ptr,res_name,res_num,
                               inchain);

              /* If residue found, and it belongs to a Pfam domain,
                 create a link between ligand and Pfam domain */
              if (residue_ptr != NULL)
                {
                  /* Add the rasdef record to the residue and chain */
                  add_ligand(residue_ptr,r_rasdef_ptr);

                  /* Get the corresponding Pfam domain */
                  pfam_ptr = residue_ptr->pfam_ptr;

                  /* If have a Pfam domain, create ligand-domain link
                     if don't already have one */
                  if (pfam_ptr != NULL)
                    {
                      /* Check previous links to see if we already have
                         this one */
                      pfamlig_ptr
                        = find_pfamlig_record(start_pfamlig_ptr,
                                              r_rasdef_ptr->ligand_name,
                                              pdb_code,inchain,
                                              ligand_id,pfam_ptr,code_ptr,
                                              residue_ptr->pfam_no);

                      /* If entry not found, create a new one */
                      if (pfamlig_ptr == NULL)
                        {
                          /* Create a new pfamlig entry */
                          pfamlig_ptr
                            = create_pfamlig_record(&first_pfamlig_ptr,
                                                    &last_pfamlig_ptr,
                                                    r_rasdef_ptr->ligand_name,
                                                    ligand_id,pdb_code,
                                                    inchain,chain_no,
                                                    pfam_ptr,code_ptr,
                                                    residue_ptr->pfam_no,
                                                    SINGLE);

                          /* Increment count of links */
                          npfamlig++;

                          /* Store the domain name in the domain list */
                          store_string(&(pfamlig_ptr->domain_list),
                                       pfam_ptr->pfam_id);

                          /* If this is the first one for the current
                             chain, store it */
                          if (start_pfamlig_ptr == NULL)
                            start_pfamlig_ptr = pfamlig_ptr;
                        }

                      /* Increment count of contacts */
                      pfamlig_ptr->ncontacts++;
                    }
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Return current pointers */
  *strt_pfamlig_ptr = start_pfamlig_ptr;
  *fst_pfamlig_ptr = first_pfamlig_ptr;
  *lst_pfamlig_ptr = last_pfamlig_ptr;

  /* Return count of ligand to Pfam domain links */
  *npflig = npfamlig;

  /* Retuirn number of ligands marked as wanted */
  return(nligand);
}
/***********************************************************************

create_plist_entry  -  Create a plist entry

***********************************************************************/

struct plist *create_plist_entry(struct plist **fst_plist_ptr,
                                 struct plist **lst_plist_ptr,
                                 struct pfam *pfam_ptr,int domain_no,
                                 int ncontacts)
{
  struct plist *plist_ptr;
  struct plist *first_plist_ptr, *last_plist_ptr;

  /* Initialise plist pointers */
  plist_ptr = NULL;
  first_plist_ptr = *fst_plist_ptr;
  last_plist_ptr = *lst_plist_ptr;

  /* Create plist entry */
  plist_ptr = (struct plist *) malloc(sizeof(struct plist));
  if (plist_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct plist\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  plist_ptr->pfam_ptr = pfam_ptr;
  plist_ptr->domain_no = domain_no;
  plist_ptr->ncontacts = ncontacts;

  /* Initialise pointer to next entry */
  plist_ptr->next_plist_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_plist_ptr == NULL)
    first_plist_ptr = plist_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_plist_ptr->next_plist_ptr = plist_ptr;
                      
  /* Save the current residue's pointer */
  last_plist_ptr = plist_ptr;

  /* Return current pointers */
  *fst_plist_ptr = first_plist_ptr;
  *lst_plist_ptr = last_plist_ptr;
  return(plist_ptr);
}
/***********************************************************************

sort_domains  -  Sort the domains using a simple bubble sort and
                 identify the one with the largest number of contacts to
                 the ligand

***********************************************************************/

void sort_domains(struct pfamlig *pfamlig_ptr)
{
  int max_contacts, swap, swapped;

  struct plist *plist_ptr, *first_plist_ptr;
  struct plist *next_plist_ptr, *prev_plist_ptr, *swap_plist_ptr;

  /* Initialise variables */
  max_contacts = 0;

  /* Get pointer to the first pfam list record */
  first_plist_ptr = pfamlig_ptr->first_plist_ptr;
  swapped = TRUE;

  /* Get pointer to first plist record */
  plist_ptr = first_plist_ptr;

  /* Loop over all the domains, to identify the one making most contacts
     with the ligand  */
  while (plist_ptr != NULL)
    {
      /* If maximum number of contacts so far, store the pointer */
      if (plist_ptr->ncontacts > max_contacts)
        {
          /* Store pointer and number of contacts */
          max_contacts = plist_ptr->ncontacts;
          pfamlig_ptr->pfam_ptr = plist_ptr->pfam_ptr;
          pfamlig_ptr->max_contacts = max_contacts;
        }

      /* Get the next plist record */
      plist_ptr = plist_ptr->next_plist_ptr;
    }

  /* Loop until Pfam domains in required order */
  while (swapped == TRUE)
    {
      /* Initialise swapped flag */
      swapped = FALSE;
          
      /* Get pointer to first plist record */
      plist_ptr = first_plist_ptr;
      prev_plist_ptr = NULL;

      /* Loop over all the domains, adjusting their order where
         necessary */
      while (plist_ptr != NULL)
        {
          /* Get pointer to the next record */
          next_plist_ptr = plist_ptr->next_plist_ptr;

          /* Determine whether current two records are in the right
             order according to domain number */
          swap = FALSE;
          if (next_plist_ptr != NULL &&
              plist_ptr->domain_no > next_plist_ptr->domain_no)
            swap = TRUE;

          /* If need to swap, perform the swap */
          if (swap == TRUE)
            {
              /* If this is the first plist record, adjust first pointer */
              if (plist_ptr == first_plist_ptr)
                first_plist_ptr = next_plist_ptr;

              /* Otherwise, get the previous plist pointer to point
                 to next one */
              else
                prev_plist_ptr->next_plist_ptr = next_plist_ptr;

              /* Swap pointers within the records */
              plist_ptr->next_plist_ptr
                = next_plist_ptr->next_plist_ptr;
              next_plist_ptr->next_plist_ptr = plist_ptr;

              /* Swap the pointers themselves */
              swap_plist_ptr = plist_ptr;
              plist_ptr = next_plist_ptr;
              next_plist_ptr = swap_plist_ptr;

              /* Set flag that swap made */
              swapped = TRUE;
            }

          /* Save current pointer */
          prev_plist_ptr = plist_ptr;

          /* Get the next plist record */
          plist_ptr = next_plist_ptr;
        }

      /* Save the first pfam list pointer */
      pfamlig_ptr->first_plist_ptr = first_plist_ptr;
    }
}
/***********************************************************************

multiple_domains  -  Identify cases where ligand binds to more than one
                     Pfam domain

***********************************************************************/

int multiple_domains(struct pfamlig *start_pfamlig_ptr,
                     struct pfamlig **fst_pfamlig_ptr,
                     struct pfamlig **lst_pfamlig_ptr,int npfamlig)
{
  char domain_list[LINELEN], pfam_id[LINELEN];

  char *string_ptr;

  int last_no, n;
  int multi_domain;

  struct pfam *pfam_ptr, *other_pfam_ptr;
  struct pfamlig *pfamlig_ptr, *new_pfamlig_ptr, *other_pfamlig_ptr;
  struct pfamlig *first_pfamlig_ptr, *last_pfamlig_ptr;
  struct plist *plist_ptr, *last_plist_ptr;

  /* Initialise pointers */
  first_pfamlig_ptr = *fst_pfamlig_ptr;
  last_pfamlig_ptr = *lst_pfamlig_ptr;

  /* Get pointer to first pfamlig record */
  pfamlig_ptr = start_pfamlig_ptr;

  /* Loop through pfamlig records to compare against all others */
  while (pfamlig_ptr != NULL)
    {
      /* Initialise flag */
      multi_domain = FALSE;
      new_pfamlig_ptr = NULL;
      last_plist_ptr = NULL;

      /* If have more than a minimum number of contacts, then take as
         binding the given ligand */
      if (pfamlig_ptr->ncontacts > MIN_CONTACTS &&
          pfamlig_ptr->done == FALSE)
        {
          /* Get the corresponding Pfam record */
          pfam_ptr = pfamlig_ptr->pfam_ptr;

          /* Initialise multiple Pfam id */
          strcpy(pfam_id,pfam_ptr->pfam_id);

          /* Get pointer to next pfamlig record */
          other_pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;

          /* Loop through pfamlig records and compare */
          while (other_pfamlig_ptr != NULL)
            {
              /* Initialise flag */
              multi_domain = FALSE;

              /* If have more than a minimum number of contacts, then
                 take as binding the given ligand */
              if (other_pfamlig_ptr->ncontacts > MIN_CONTACTS &&
                  other_pfamlig_ptr->done == FALSE)
                {
                  /* Check whether domains are different and bind the
                     same ligand */
                  if (pfamlig_ptr->domain_no !=
                      other_pfamlig_ptr->domain_no &&
                      !strcmp(pfamlig_ptr->ligand_id,
                              other_pfamlig_ptr->ligand_id))
                    {
                      /* Mark both domains as part of a multiple
                         domain */
                      pfamlig_ptr->type = PART_OF_MULTIPLE;
                      other_pfamlig_ptr->type = PART_OF_MULTIPLE;

                      /* Get the corresponding Pfam record */
                      other_pfam_ptr = other_pfamlig_ptr->pfam_ptr;

                      /* See if we already have this domain */
                      string_ptr
                        = strstr(pfam_id,other_pfam_ptr->pfam_id);

                      /* If not found, then add */
                      if (string_ptr == NULL)
                        {
                          /* Append current Pfam id */
                          strcat(pfam_id,".");
                          strcat(pfam_id,other_pfam_ptr->pfam_id);

                          /* Set flag that have a multi-domain entry */
                          multi_domain = TRUE;
                        }
                    }
                }

              /* If domain is part of a multi-domain contact, store
                 the appropriate details */
              if (multi_domain == TRUE)
                {
                  /* Mark the current record as done */
                  other_pfamlig_ptr->done = TRUE;

                  /* If this is a new multiple, create the
                     pfamlig record to hold the details */
                  if (new_pfamlig_ptr == NULL)
                    {
                      /* Create a new pfamlig record */
                      new_pfamlig_ptr
                        = create_pfamlig_record(&first_pfamlig_ptr,
                                                &last_pfamlig_ptr,
                                                pfamlig_ptr->ligand_name,
                                                pfamlig_ptr->ligand_id,
                                                pfamlig_ptr->pdb_code,
                                                pfamlig_ptr->chain,
                                                pfamlig_ptr->seq_no,
                                                pfam_ptr,
                                                pfamlig_ptr->code_ptr,
                                                pfamlig_ptr->domain_no,
                                                MULTIPLE);

                      /* Increment count of records */
                      npfamlig++;

                      /* Add first Pfam domain to it */
                      plist_ptr
                        = create_plist_entry(&(new_pfamlig_ptr->first_plist_ptr),
                                             &last_plist_ptr,pfam_ptr,
                                             pfamlig_ptr->domain_no,
                                             pfamlig_ptr->ncontacts);
                    }

                  /* Add the current Pfam domain */
                  plist_ptr
                    = create_plist_entry(&(new_pfamlig_ptr->first_plist_ptr),
                                         &last_plist_ptr,other_pfam_ptr,
                                         other_pfamlig_ptr->domain_no,
                                         other_pfamlig_ptr->ncontacts); 
               }

              /* Get next pfamlig record */
              other_pfamlig_ptr = other_pfamlig_ptr->next_pfamlig_ptr;
            }
        }

      /* If we created a new pfamlig record, bubble sort the domains */
      if (new_pfamlig_ptr != NULL)
        {
          /* Sort the domains and identify the one with the largest number
             of contacts to the ligand */
          sort_domains(new_pfamlig_ptr);

          /* Initialise domain-list and count */
          domain_list[0] = '\0';
          last_no = -1;

          /* Get first of the Pfam domains */
          plist_ptr = new_pfamlig_ptr->first_plist_ptr;

          /* Loop through the Pfam domains to list */
          while (plist_ptr != NULL)
            {
              /* Determine number of domains since last one */
              if (last_no > -1)
                {
                  /* Determine number of domains since last one */
                  // n = plist_ptr->domain_no - last_no - 1;

                  /* Add separating dot */
                  strcat(domain_list,".");
                }

              /* Add domain to domain string */
              strcat(domain_list,plist_ptr->pfam_ptr->pfam_id);

              /* Store this domain number */
              last_no = plist_ptr->domain_no;

              /* Get next plist record */
              plist_ptr = plist_ptr->next_plist_ptr;
            }

          /* Store the domain list */
          store_string(&(new_pfamlig_ptr->domain_list),domain_list);
        }

      /* Get next pfamlig record */
      pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;
    }

  /* Return current pointers */
  *fst_pfamlig_ptr = first_pfamlig_ptr;
  *lst_pfamlig_ptr = last_pfamlig_ptr;

  /* Return count of npfamlig records */
  return(npfamlig);
}
/***********************************************************************

flag_ligands  -  Mark a single residue record for each wanted ligand

***********************************************************************/

void flag_ligands(struct r_rasdef *first_r_rasdef_ptr)
{
  int ires, ligand_no, type;

  struct r_rasdef *r_rasdef_ptr, *other_rasdef_ptr;

  /* Initialise variables */

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop to flag first residue of each wanted ligand */
  while (r_rasdef_ptr != NULL)
    {
      /* If residue wanted, mark first residue belonging to the parent
         ligand */
      if (r_rasdef_ptr->wanted == TRUE && 
          (r_rasdef_ptr->type == LIGAND ||
           r_rasdef_ptr->type == METAL ||
           r_rasdef_ptr->type == NUCLEIC_ACID))
        {
          /* Get the ligand number and type */
          ligand_no = r_rasdef_ptr->number;
          type = r_rasdef_ptr->type;
          ires = 0;

          /* Get pointer to the first RasMol definition again */
          other_rasdef_ptr = first_r_rasdef_ptr;
              
          /* Loop to write out the current ligand */
          while (other_rasdef_ptr != NULL)
            {
              /* If this is the right ligand number, write out */
              if (other_rasdef_ptr->type == type &&
                  other_rasdef_ptr->number == ligand_no)
                {
                  /* If not the first member of this ligand, then
                     mark as not wanted */
                  if (ires == 0)
                    other_rasdef_ptr->wanted = TRUE;
                  else
                    other_rasdef_ptr->wanted = FALSE;
                  ires++;
                }

              /* Get next definition record */
              other_rasdef_ptr = other_rasdef_ptr->next_r_rasdef_ptr;
            }
        }

      /* Otherwise mark as unwanted */
      else
        r_rasdef_ptr->wanted = FALSE;

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }
}
/***********************************************************************

unique_ligands  -  Combine ligands of same name

***********************************************************************/

void unique_ligands(struct r_rasdef *first_r_rasdef_ptr)
{
  int type;

  struct r_rasdef *r_rasdef_ptr, *other_rasdef_ptr;

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop to combine duplicate ligands */
  while (r_rasdef_ptr != NULL)
    {
      if (r_rasdef_ptr->wanted == TRUE)
        {
          /* Get ligand type */
          type = r_rasdef_ptr->type;

          /* Set count for this ligand to 1 */
          r_rasdef_ptr->ncopy = 1;

          /* Get pointer to the next RasMol definition */
          other_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
              
          /* Loop to write out the current ligand */
          while (other_rasdef_ptr != NULL)
            {
              /* If ligand is wanted and has same  */
              if (other_rasdef_ptr->type == type &&
                  other_rasdef_ptr->wanted == TRUE &&
                  !strcmp(other_rasdef_ptr->ligand_name,
                          r_rasdef_ptr->ligand_name))
                {
                  /* Mark second copy as not wanted and increment
                     count of first */
                  other_rasdef_ptr->wanted = FALSE;
                  r_rasdef_ptr->ncopy++;
                }

              /* Get next definition record */
              other_rasdef_ptr = other_rasdef_ptr->next_r_rasdef_ptr;
            }
        }

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }
}
/***********************************************************************

write_ligands  -  Write out all ligands/metals/DNA that interact
                  with the given chain

***********************************************************************/

int write_ligands(struct r_rasdef *first_r_rasdef_ptr,FILE *outfile_ptr)
{
  int nout;

  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  nout = 0;

  /* Get pointer to the first RasMol definition record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop to write out all the interacting ligands, metals and DNA */
  while (r_rasdef_ptr != NULL)
    {
      if (r_rasdef_ptr->wanted == TRUE)
        {
          /* If not first, precede by comma */
          if (nout == 0)
            fprintf(outfile_ptr,"; ");
          else
            fprintf(outfile_ptr,", ");
          nout++;

          /* Write out */
          fprintf(outfile_ptr,"%s",r_rasdef_ptr->ligand_name);

          /* If more than a single copy, write out number */
          if (r_rasdef_ptr->ncopy > 1)
            fprintf(outfile_ptr,"x%d",r_rasdef_ptr->ncopy);
        }

      /* Get next definition record */
      r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;
    }

  /* Return number of ligands written out */
  return(nout);
}
/***********************************************************************

free_rasdefs  -  Free up the memory taken up by the current entry

***********************************************************************/

void free_rasdefs(struct r_rasdef *first_r_rasdef_ptr)
{
  struct r_rasdef *r_rasdef_ptr, *next_r_rasdef_ptr;

  /* Get first rasdef record */
  r_rasdef_ptr = first_r_rasdef_ptr;

  /* Loop through the rasdefs until done */
  while (r_rasdef_ptr != NULL)
    {
      /* Get pointer to the next rasdef record */
      next_r_rasdef_ptr = r_rasdef_ptr->next_r_rasdef_ptr;

      /* Free up memory taken by current rasdef */
      free(r_rasdef_ptr);

      /* Go to the next rasdef */
      r_rasdef_ptr = next_r_rasdef_ptr;
    }
}
/***********************************************************************

free_residue  -  Free up the memory taken up by the current entry's data
                details

***********************************************************************/

void free_residue(struct residue **fst_residue_ptr)
{
  struct residue *residue_ptr, *next_residue_ptr;

  /* Get first residue record */
  residue_ptr = *fst_residue_ptr;

  /* Loop through the residues until done */
  while (residue_ptr != NULL)
    {
      /* Get pointer to the next residue record */
      next_residue_ptr = residue_ptr->next_residue_ptr;

      /* Free up memory taken by current residue */
      free(residue_ptr);

      /* Go to the next residue */
      residue_ptr = next_residue_ptr;
    }

  /* Re-initialise first record */
  *fst_residue_ptr = NULL;
}
/***********************************************************************

find_cath_record  -  Find the given cath entry

***********************************************************************/

struct cath *find_cath_record(struct cath *first_cath_ptr,
                              char *cath_code)
{
  struct cath *cath_ptr, *found_cath_ptr;

  /* Initialise variables */
  found_cath_ptr = NULL;
  if (cath_code[0] == '\0')
    return(found_cath_ptr);

  /* Get pointer to the first cath */
  cath_ptr = first_cath_ptr;

  /* Loop over all the caths to find the one just read in */
  while (cath_ptr != NULL && found_cath_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (!strcmp(cath_ptr->cath_code,cath_code))
        found_cath_ptr = cath_ptr;

      /* Get pointer to the next cath */
      cath_ptr = cath_ptr->next_cath_ptr;
    }

  /* Return the matching cath pointer */
  return(found_cath_ptr);
}
/***********************************************************************

create_cath_record  -  Create a new cath record

***********************************************************************/

struct cath *create_cath_record(struct cath **fst_cath_ptr,
                                struct cath **lst_cath_ptr,
                                char *cath_code)
{
  struct cath *cath_ptr;
  struct cath *first_cath_ptr, *last_cath_ptr;

  /* Initialise cath pointers */
  cath_ptr = NULL;
  first_cath_ptr = *fst_cath_ptr;
  last_cath_ptr = *lst_cath_ptr;

  /* Create cath entry */
  cath_ptr = (struct cath*)malloc(sizeof(struct cath));
  if (cath_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct cath\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store info that we have so far */
  store_string(&(cath_ptr->cath_code),cath_code);

  /* Initialise other fields */
  cath_ptr->class = &null;
  cath_ptr->architecture = &null;

  /* Initialise pointer to next cath record */
  cath_ptr->next_cath_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of cathed cath */
  if (first_cath_ptr == NULL)
    first_cath_ptr = cath_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_cath_ptr->next_cath_ptr = cath_ptr;
                      
  /* Save the current cath's pointer */
  last_cath_ptr = cath_ptr;

  /* Return current pointers */
  *fst_cath_ptr = first_cath_ptr;
  *lst_cath_ptr = last_cath_ptr;
  return(cath_ptr);
}
/***********************************************************************

identify_cath_domains  -  Get the CATH domain assignments from the
                          database.dat file

***********************************************************************/

void identify_cath_domains(struct pdb *pdb_ptr,struct chain *chain_ptr,
                           char *pdbsum_dir)
{
  char chain, input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char cath_code[30], number_string[20], pdb_code[5];

  char *string_ptr;

  int domain, len;
  int done, have_chain, have_file;

  struct cath *cath_ptr, *first_cath_ptr, *last_cath_ptr;
  struct domain *cath_domain_ptr, *last_cath_domain_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  strcpy(pdb_code,pdb_ptr->pdb_code);
  chain = chain_ptr->chain_id;
  have_chain = FALSE;

  /* Initialise pointers */
  chain_ptr->first_cath_domain_ptr = NULL;
  chain_ptr->ncath_dom = 0;
  last_cath_domain_ptr = NULL;
  first_cath_ptr = NULL;
  last_cath_ptr = NULL;

  /* Form the name of the PDBsum database.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"database.dat");
  have_file = FALSE;

  /* Open the database.dat file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL &&
             done == FALSE)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* Separate the field name from the field value */
          string_ptr = strstr(input_line," ");

          /* If this is the CHAIN record, check whether it is the
             right chain */
          if (!strncmp(input_line,"CHAIN[",6))
            {
              /* If we already have the right chain, must have processed
                 it, so can end here */
              if (have_chain == TRUE)
                done = TRUE;

              /* If no chain-id, check iof our chain is blank */
              else if (string_ptr == NULL && chain == ' ')
                have_chain = TRUE;

              /* Otherwise, check whether chain matches */
              else if (string_ptr != NULL && string_ptr[1] == chain)
                have_chain = TRUE;

              /* Otherwise, this is not our chain */
              else
                have_chain = FALSE;

              /* initialise variables */
              cath_ptr = NULL;
              domain = -1;
            }

          /* If reading through our chain, pick up the relevant
             CATH details */
          else if (have_chain == TRUE && string_ptr != NULL)
            {
              /* If this is the wolf code, get the domain number */
              if (!strncmp(input_line,"WOLF_CODE[",10))
                {
                  /* Extract the domain number */
                  strcpy(number_string,string_ptr+6);
                  domain = atoi(number_string);

                  /* Id domain number is zero, renumber it as 1 */
                  if (domain == 0)
                    domain = 1;
                }

              /* If CATH code, find/create CATH record */
              else if (!strncmp(input_line,"CATH_CODE[",10))
                {
                  /* Get the CATH code */
                  strcpy(cath_code,string_ptr+1);

                  /* See if we already have this record */
                  cath_ptr = find_cath_record(first_cath_ptr,cath_code);
                  
                  /* If don't yet have this record, create a new one */
                  if (cath_ptr == NULL)
                    {
                      /* Create a new CATH record */
                      cath_ptr
                        = create_cath_record(&first_cath_ptr,
                                             &last_cath_ptr,cath_code);
                    }

                  /* Create a CATH domain record for this domain */
                  cath_domain_ptr
                    = create_domain_entry(&(chain_ptr->first_cath_domain_ptr),
                                          &last_cath_domain_ptr,CATH,
                                          NULL,cath_code,0,0);

                  /* Store pointer to CATH record */
                  cath_domain_ptr->cath_ptr = cath_ptr;
                  cath_domain_ptr->domain_no = domain;

                  /* Increment count of domains for this chain */
                  chain_ptr->ncath_dom++;

                  /* Update PDB records flag to indicate that we have
                     CATH domains */
                  pdb_ptr->have_cath = TRUE;
                }

              /* If CATH class, store in CATH record */
              else if (!strncmp(input_line,"CATH_CLASS[",11) &&
                       cath_ptr != NULL && cath_ptr->class == &null)
                {
                  /* Store the CATH class */
                  store_string(&(cath_ptr->class),string_ptr+1);
                }

              /* If CATH architecture, store in CATH record */
              else if (!strncmp(input_line,"CATH_ARCH[",10) &&
                       cath_ptr != NULL &&
                       cath_ptr->architecture == &null)
                {
                  /* Store the CATH architecture */
                  store_string(&(cath_ptr->architecture),string_ptr+1);
                }

            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* Otherwise, if have no database.dat file, show warning */
  else
    printf("*** Warning. File missing: %s\n",file_name);
}
/***********************************************************************

match_cath_domains  -  Renumber domains for each residue so that they
                       match the CATH domain numbers

***********************************************************************/

void match_cath_domains(struct chain *chain_ptr,int *cath_domain,
                        char *sequence,int align_len)
{
  char cath_code[30];
  
  int domain, ipos, ires, last_domain, mindom, maxdom;
  int ok;

  struct cath *cath_ptr;
  struct domain *cath_domain_ptr, *domain_ptr, *last_cath_domain_ptr;

  /* Initialise variables */
  cath_code[0] = '\0';
  mindom = maxdom = cath_domain[0];
  ok = TRUE;

  /* Loop over the domain assignments to count how many there are */
  for (ipos = 0; ipos < align_len; ipos++)
    {
      /* Update minimum and maximim domain numbers */
      if (cath_domain[ipos] <  mindom)
        mindom = cath_domain[ipos];
      if (cath_domain[ipos] >  maxdom)
        maxdom = cath_domain[ipos];
    }

  /* If number of domains does not match the top number, then need to
     renumber the domains */
  if (maxdom != chain_ptr->ncath_dom)
    {
      /* If all domains numbered as -1 or zero, then renumber as 1 */
      if ((maxdom == -1 || maxdom == 0) && chain_ptr->ncath_dom == 1)
        {
          /* Loop over all the alignment positions */
          for (ipos = 0; ipos < align_len; ipos++)
            {
              /* If this corresponds to the top domain, renumber it */
              if (cath_domain[ipos] == maxdom)
                cath_domain[ipos] = 1;
            }
        }

      /* Otherwise, can't reconcile domain assignments */
      else
        {
          /* Show warning message and unset flag */
          printf("*** Warning. Can't reconcile CATH domains for");
          printf(" %s, chain '%c'\n",chain_ptr->pdb_ptr->pdb_code,
                 chain_ptr->chain_id);
          ok = FALSE;
        }
    }

  /* If all OK, create any additional domain records that may be
     required */
  if (ok == TRUE)
    {
      /* Initialise pointer to last domain record */
      last_cath_domain_ptr = NULL;

      /* Get the first of this sequence's CATH domains */
      cath_domain_ptr = chain_ptr->first_cath_domain_ptr;

      /* Loop through domains to record last one */
      while (cath_domain_ptr != NULL)
        {
          /* Store pointer as last encountered */
          last_cath_domain_ptr = cath_domain_ptr;

          /* Get the next domain */
          cath_domain_ptr = cath_domain_ptr->next_domain_ptr;
        }

      /* Initialise variables */
      cath_domain_ptr = NULL;
      last_domain = -99;
      ires = 0;

      /* Loop over the residues, storing start and end residue numbers
         and creating any additional domain records for split domains */
      for (ipos = 0; ipos < align_len; ipos++)
        {
          /* If residue not blank and not a hyphen, increment residue no */
          if (sequence[ipos] != ' ' && sequence[ipos] != '-')
            ires++;

          /* Get the domain at this alignment position */
          domain = cath_domain[ipos];

          /* If this is the same domain, update end-residue position */
          if (domain == last_domain && cath_domain_ptr != NULL)
            cath_domain_ptr->end_res = ires;

          /* If this is a new domain, match to the CATH domain record */
          else if (domain != last_domain && domain > 0)
            {
              /* Get the first of the CATH domain records */
              cath_domain_ptr = NULL;
              cath_ptr = NULL;
              domain_ptr = chain_ptr->first_cath_domain_ptr;

              /* Loop through records to find corresponding domain
                 record */
              while (domain_ptr != NULL && cath_domain_ptr == NULL)
                {
                  /* If domain matches number and hasn't already been
                     updated, then this is the one we want */
                  if (domain_ptr->domain_no == domain &&
                      domain_ptr->start_res == 0)
                    {
                      /* Save this domain pointer */
                      cath_domain_ptr = domain_ptr;
                    }

                  /* If record is the correct domain, but from a
                     previous segment, store the CATH pointer */
                  else if (domain_ptr->domain_no == domain)
                    {
                      cath_ptr = domain_ptr->cath_ptr;
                    }

                  /* Get the next domain record */
                  domain_ptr = domain_ptr->next_domain_ptr;
                }

              /* If clean domain record not found, then create one */
              if (cath_domain_ptr == NULL)
                {
                  /* Create new CATH domain record */
                  cath_domain_ptr
                    = create_domain_entry(&(chain_ptr->first_cath_domain_ptr),
                                          &last_cath_domain_ptr,CATH,
                                          NULL,cath_code,0,0);

                  /* Store pointer to CATH record */
                  cath_domain_ptr->cath_ptr = cath_ptr;
                  cath_domain_ptr->domain_no = domain;

                  /* Increment count of domains for this chain */
                  chain_ptr->ncath_dom++;
                }

              /* Store start residue position */
              cath_domain_ptr->start_res = ires;

              /* Store the last domain */
              last_domain = domain;
            }
        }
    }
}
/***********************************************************************

free_clist  -  Free up the memory taken up by the current entry's data
               details

***********************************************************************/

void free_clist(struct clist **fst_clist_ptr)
{
  struct clist *clist_ptr, *next_clist_ptr;

  /* Get first clist record */
  clist_ptr = *fst_clist_ptr;

  /* Loop through the clists until done */
  while (clist_ptr != NULL)
    {
      /* Get pointer to the next clist record */
      next_clist_ptr = clist_ptr->next_clist_ptr;

      /* Free up memory taken by current clist */
      free(clist_ptr);

      /* Go to the next clist */
      clist_ptr = next_clist_ptr;
    }

  /* Re-initialise first record */
  *fst_clist_ptr = NULL;
}
/***********************************************************************

get_pfam_colour  -  Get the colour defined by the Pfam domain code

***********************************************************************/

int get_pfam_colour(char code,int *rgb,int *trgb,int dpos)
{
  char *string_ptr;

  int colour, i, red, green, blue;

  float dist;

  /* Initialise variables */
  colour = WHITE_COL;

  /* Determine which colour this code represents */
  string_ptr = strchr(COLOUR_CODE,code);
  if (string_ptr != NULL)
    colour = string_ptr - COLOUR_CODE;

  /* Get the corresponding RGB values for this colour */
  sscanf(COLOUR_DEFN[colour][2],"%d %d %d",&red,&green,&blue);

  /* Store the RGB values */
  rgb[0] = red;
  rgb[1] = green;
  rgb[2] = blue;

  /* If have a valid domain colour (ie not white), calculate colour of
     text label giving domain name to go on the bar */
  if (colour != WHITE_COL)
    {
      /* Determine corresponding text colour */
      dist = red * red + green * green + blue * blue;
      if (dist > 0)
        dist = sqrt(dist);
      if (dist > HALF_TEXT)
        trgb[0] = trgb[1] = trgb[2] = 0;
      else
        trgb[0] = trgb[1] = trgb[2] = 255;
    }

  /* Adjust shade according to row-position */
  if (dpos >= 0 && dpos < PFAMA_HEIGHT)
    {
      /* Adjust each component of the colour */
      for (i = 0; i < 3; i++)
        {
          /* Check for the highlight row */
          if (dpos == 3)
            rgb[i] = 255;

          /* Adjust percentage drop to zero */
          else if (rgb[i] > 0)
            rgb[i] = (int) ((float) rgb[i]
                            * (float) DOMAIN_SHADE[dpos] / 100.0);
        }
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

get_domain_colour  -  Get the colour of the given residue, according to
                      its CATH domain assignment

***********************************************************************/

int get_domain_colour(int domain)
{
#define NDOMAIN_COLOURS   8

  int dcol, icol, colour;

  static char *colour_name[NDOMAIN_COLOURS] = {
    "Black", "Red", "Blue", "Green", "SaddleBrown", "Orange",
    "Cyan", "HotPink" };

  /* Initialise colour */
  colour = 0;

  /* Determine which colour this corresponds to */
  dcol = domain % NDOMAIN_COLOURS;

  /* Loop over the colour definitions to get the corresponding
     colour number */
  for (icol = 0; icol < NPFAM_COLOURS; icol++)
    {
      /* If colour name matches, then store corresponding colour
         number */
      if (!strcmp(colour_name[dcol],COLOUR_DEFN[icol][0]))
        colour = icol;
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

generate_image  -  Generate a .ppm image file holding the secondary
                   structure assignments

***********************************************************************/

void generate_image(char *array[PIC_HEIGHT],int width,
                    struct pdb *pdb_ptr,struct chain *chain_ptr,
                    int icode,int plot_no,int image_height,
                    char *pdbsum_dir,int image_split,char *res_range,
                    struct c_file_data file_name_ptr,
                    struct parameters params)
{
  char ch, cnum[2], last_ch;
  char command[FILENAME_LEN], file_name[FILENAME_LEN];
  char ppmfile_name[FILENAME_LEN], number_string[20];

  int dpos, i, image, image_width, inum, irgb[3], iwidth, last_width;
  int argb[3], image_no, irange, lrgb[3], prgb[3], trgb[3];
  int first_image, height, nimages, ntrip, xpos, ypos, ystart, yend;
  int dos, first;

  double rgb[3], rgb_light[3];

  FILE *file_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  last_ch = '\0';
  for (i = 0; i < 3; i++)
    prgb[i] = 0;
  trgb[0] = trgb[1] = trgb[2] = 0;

  /* Get the chain colour */
  get_chain_colour(chain_ptr->chain_id,rgb,rgb_light);

  /* Convert RGB value to integers */
  for (i = 0; i < 3; i++)
    {
      irgb[i] = (int) (255 * rgb[i]);
      lrgb[i] = (int) (255 * rgb_light[i]);
    }

  /* Determine how many images we need to have */
  nimages = width / IMAGE_WIDTH + 1;
  last_width = image_width = IMAGE_WIDTH;
  first_image = 0;

  /* If final image particularly small, tag onto end of previous image */
  if (nimages > 1 && (width % IMAGE_WIDTH) < IMAGE_ALLOW)
    {
      /* Reduce number of images by one and set width of last image */
      nimages--;
      last_width = image_width + width % IMAGE_WIDTH;
    }

  /* Otherwise, calculate width of final image */
  else
    last_width = width - (nimages - 1) * IMAGE_WIDTH;

  /* For COVID-19 protein page, have just a single image */
  if (params.covid19 == TRUE)
    {
      nimages = 1;
      image_width = width;
      last_width = width;
      first_image = -1;
    }

  /* Store number of images and their dimensions */
  chain_ptr->nimages[icode] = nimages;

  chain_ptr->image_width = image_width;
  chain_ptr->last_width = last_width;

  /* Loop over the images to be generated */
  for (image = first_image; image < nimages; image++)
    {
      /* If this is the last image, adjust width */
      if (image == nimages - 1)
        image_width = last_width;

      /* Get the image number */
      image_no = image;
      if (image < 0)
        image_no = 0;

      /* Form name of .ppm file */
      strcpy(ppmfile_name,"miniwir.ppm");

      /* Open the output file */
      file_ptr = open_output_file(ppmfile_name,TRUE);

      /* Define the start end rows of this image */
      height = image_height;
      ystart = image_height - 1;
      yend = -1;

      /* For COVID-19 protein page, define start and end depending on
         whether generating domain image or the secondary structure
         image */
      if (params.covid19 == TRUE && image_split > 0)
        {
          /* Get the start and end row, depening on which image this is */
          if (image == first_image)
            {
              /* First image is just the domain diagram */
              ystart = image_height - 1;
              yend = image_height - image_split;
            }
          else
            {
              /* Second image is just the secondary structure diagram */
              ystart = image_height - image_split;
              yend = -1;
            }

          /* Get the corresponding image height */
          height = ystart - yend;
        }

      /* Write out header records */
      fprintf(file_ptr,"P3\n");
      fprintf(file_ptr,"# miniwir.ppm\n");
      fprintf(file_ptr,"%d %d\n",image_width,height);
      fprintf(file_ptr,"255\n");

      /* Initialise triplet count and domain-shade position */
      ntrip = 0;
      dpos = -1;

      /* Loop through the rows making up the image */
      for (ypos = ystart; ypos > yend; ypos--)
        {
          /* Initialise flag indicating when first Pfam domain
             encountered in row */
          first = TRUE;
          last_ch = '\0';

          /* Loop across the width of the image */
          for (iwidth = 0; iwidth < image_width; iwidth++)
            {
              /* Get the x-position in array */
              xpos = iwidth + image_no * IMAGE_WIDTH;

              /* If this pixel blank, then write out white */
              if (array[ypos][xpos] == ' ')
                fprintf(file_ptr,"255 255 255 ");

              /* If a sequence-break dot, then write grey */
              else if (array[ypos][xpos] == '.')
                fprintf(file_ptr,"%d %d %d ",RGB_BREAK[0],RGB_BREAK[1],
                        RGB_BREAK[2]);

              /* If a black mark show in black */
              else if (array[ypos][xpos] == '#')
                fprintf(file_ptr,"0 0 0 ");

              /* If cross, show grey */
              else if (array[ypos][xpos] == '+')
                fprintf(file_ptr,"128 128 128 ");

              /* If part of a label character, show in text colour */
              else if (array[ypos][xpos] == '$')
                fprintf(file_ptr,"%d %d %d ",trgb[0],trgb[1],trgb[2]);

              /* If asterisk, write out in this chain's colour */
              else if (array[ypos][xpos] == '*')
                fprintf(file_ptr,"%d %d %d ",irgb[0],irgb[1],irgb[2]);

              /* If a number, write out a modified version of the chain
                 colour */
              else if (array[ypos][xpos] >= '0' &&
                       array[ypos][xpos] <= '9')
                {
                  /* Get the number */
                  cnum[0] = array[ypos][xpos];
                  cnum[1] = '\0';
                  inum = atoi(cnum);

                  /* If have an AF score, change the chain colour according to
                     the confidence range */
                  if (res_range[xpos] != ' ')
                    {
                      /* Extract the range number */
                      number_string[0] = res_range[xpos];
                      number_string[1] = '\0';
                      irange = atoi(number_string);

                      /* Copy across the range colour components */
                      irgb[0] = AF_COLOUR_IRGB[irange][0];
                      irgb[1] = AF_COLOUR_IRGB[irange][1];
                      irgb[2] = AF_COLOUR_IRGB[irange][2];
                    }

                  /* Compute the adjusted values */
                  for (i = 0; i < 3; i++)
                    {
                      /* If lighter colour, make lighter */
                      if (ADJUST_COL[inum] <= 1.0)
                        {
                          argb[i]
                            = (int) (irgb[i]
                            + (255.0 - irgb[i])
                            * (1.0 - ADJUST_COL[inum]));
                          if (argb[i] > 255)
                            argb[i] = 255;
                        }

                      /* Otherwise, make darker */
                      else
                        {
                          argb[i]
                            = (int) ((float) irgb[i]
                                     * (2.0 - ADJUST_COL[inum]));
                        }
                    }

                  /* Write out */
                  fprintf(file_ptr,"%d %d %d ",argb[0],argb[1],argb[2]);
                }

              /* If '@' symbol, write out in lighter version of chain's
                 colour */
              else if (array[ypos][xpos] == '@')
                fprintf(file_ptr,"%d %d %d ",lrgb[0],lrgb[1],lrgb[2]);

              /* Otherwise get Pfam domain colour */
              else
                {
                  /* If this if the first on this row, increment
                     domain-column position and unset flag */
                  if (first == TRUE)
                    {
                      dpos++;
                      first = FALSE;
                    }

                  /* Get the colour character */
                  ch = array[ypos][xpos];

                  /* If this is a new colour, need to get its RGB
                     values */
                  if (ch != last_ch)
                    {
                      /* Get the Pfam colour corresponding to this code */
                      get_pfam_colour(ch,prgb,trgb,dpos);

                      /* Save the colour character */
                      last_ch = ch;
                    }

                  /* Write out colour */
                  fprintf(file_ptr,"%d %d %d ",prgb[0],prgb[1],prgb[2]);
                }

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

      /* Form name of output GIF file */
      strcpy(file_name,pdbsum_dir);
      if (params.covid19 == TRUE)
        strcat(file_name,"covplot");
      else
        strcat(file_name,"uniplot");
      sprintf(number_string,"%3.3d",plot_no);
      strcat(file_name,number_string);
      strcat(file_name,"_");
      sprintf(number_string,"%2.2d",image + 1);
      strcat(file_name,number_string);
      strcat(file_name,".gif");

      /* Convert to a .gif file */
      strcpy(command,file_name_ptr.exe_name[CONVERT_EXE]);
      strcat(command," ");
      strcat(command,ppmfile_name);
      strcat(command," ");
      strcat(command,file_name);
      strcat(command," > convert.log");

      /* Run convert */
#ifdef NOSYST
#else
      call_system(command,dos);
#endif
    }
}
/***********************************************************************

plot_miniwir  -  Plot the mini wiring diagram

***********************************************************************/

int plot_miniwir(struct chain *chain_ptr,int icode,
                 char sequence[2][MAX_SEQLEN],char *sec_struc,
                 char *res_range,struct pdb *pdb_ptr,char *character,
                 struct letter *letter_ptr,int nletters,
                 char *pdbsum_dir,struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char ch, range, sst;
  char label_array[LETTER_ROWS][MAX_LABEL_WIDTH];

  char *array[PIC_HEIGHT], *blank;

  int align_len, cpos, i, icol, ipos, ires, irow, j, jpos, len, start_col;
  int helix_len, ncols, next_j, nmismatch, nrows, plot, last_plot;
  int run_len, start_pos, strand_len;
  int colour, end_res, start_res, type, xpos, ypos;
  int domain, domain_height, domain_width, iletter, height, width;
  int image_height, image_split, len_defn, plot_no, prot_no;
  int endy, label_width, label_height, starty, xoffset, yoffset;
  int plot_domains, used[LETTER_ROWS];

  struct code *code_ptr;
  struct domain *domain_ptr, *cath_domain_ptr, *last_domain_ptr;
  struct smatch *smatch_ptr;

  /* Initialise variables */
  helix_len = 0;
  image_split = 0;
  last_plot = COIL;
  nmismatch = 0;
  run_len = 0;
  strand_len = 0;
  last_domain_ptr = NULL;

  /* Get the corresponding smatch and code records */
  code_ptr = chain_ptr->code_ptr[icode];
  smatch_ptr = chain_ptr->smatch_ptr[icode];

  /* Get length of UniProt vs PDB sequence alignment */
  align_len = strlen(sequence[PDB_SEQ]);
  smatch_ptr->align_len = align_len;

  /* Get the corresponding UniProt sequence record */
  plot_no = chain_ptr->plot_no[icode];
  // plot_no = chain_ptr->chain_no;

  /* If length of sequence as given in Pfam matches that in UniProt,
     then can annotate sequence with UniProt domains */
  plot_domains = TRUE;
  //if (abs(code_ptr->seq_len - code_ptr->pfam_seq_len) < 4)
  //  plot_domains = TRUE;
  //else
  //  {
      /* Show warning and set flag */
  //    printf("*** Warning. Seqlen mismatch UniProt & Pfam (%s)"
  //           " for %s(%c): %d & %d\n",code_ptr->uniprot_acc,pdb_ptr->pdb_code,
  //           chain_ptr->chain_id,code_ptr->seq_len,
  //           code_ptr->pfam_seq_len);
  //    plot_domains = FALSE;
  //  }
  if (code_ptr->first_domain_ptr == NULL)
    plot_domains = FALSE;

  /* Determine size of array required for storing mini wiring diagram
     image */
  width = align_len;

  /* Create the array to hold image */
  for (ypos = 0; ypos < PIC_HEIGHT; ypos++)
    {
      array[ypos] = (char *) malloc((width + 1) * sizeof(char));

      /* Check that sufficient memory was allocated */
      if (array[ypos] == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for image "
                 "array\n");
          exit(1);
        }

      /* Initialise this line of the array */
      for (xpos = 0; xpos < width; xpos++)
        array[ypos][xpos] = ' ';

      /* Terminate line */
      array[ypos][width] = '\0';
    }

  /* Determine the actual picture height */
  image_height = WIR_HEIGHT + 1;
  if (smatch_ptr->first_domain_ptr != NULL && pdb_ptr->have_cath == TRUE)
    image_height = PIC_HEIGHT;
  else if (smatch_ptr->first_domain_ptr != NULL)
    image_height = image_height + WIR_GAP + SEQ_HEIGHT;
  else if (pdb_ptr->have_cath == TRUE)
    image_height = image_height + CATH_HEIGHT + WIR_GAP;
  else if (sequence[UNIPROT_SEQ][0] != '\0')
    image_height = image_height + SEQ_HEIGHT;

  /* Create array to hold info on which pixels are in use and initialise */
  blank = (char *) malloc((width + 1) * sizeof(char));
  for (xpos = 0; xpos < width; xpos++)
    blank[xpos] = ' ';

  /* Calculate y-position of coil midpoint */
  cpos = WIR_HEIGHT / 2;
  if (pdb_ptr->have_cath == TRUE)
    cpos = cpos + CATH_HEIGHT;

  /* Loop over the positions in the array */
  for (xpos = 0; xpos < width; xpos++)
    {
      /* Initialise the secondary structure assignments */
      sst = ' ';
      ires = xpos;

      /* Initialise plot type */
      plot = GAP_REGION;

      /* Get the secondary structure assignment for current residue */
      sst = sec_struc[ires];
      range = res_range[ires];

      /* Secondary structure is helix */
      if (sst == 'H')
        plot = HELIX;

      /* Or strand */
      else if (sst == 'E')
        plot = STRAND;

      /* Otherwise, if not a blank, plot as coil */
      else if (sst != ' ')
        plot = COIL;

      /* If have strand-end, plot full length of strand */
      if ((plot != STRAND && strand_len > 0) ||
          (plot == STRAND && xpos == (width - 1)))
        {
          /* Start strand at previous residue position */
          start_pos = xpos - 1;

          /* If this is the last residue position, increment strand
             length and start strand at current position */
          if (plot == STRAND && xpos == (width - 1))
            {
              strand_len++;
              start_pos = xpos;
            }

          /* Plot the full length of the strand */
          for (i = 0; i < STRAND_HEIGHT; i++)
            {
              /* Get y-position of current point */
              ypos = cpos - STRAND_HEIGHT / 2 + i;
              j = LEN_STRAND_DEFN;

              /* Loop over the strand length */
              for (jpos = 0; jpos < strand_len; jpos++)
                {
                  /* Determine which part of the strand to plot */
                  j--;
                  if (j < 0)
                    j = 0;

                  /* Mark this point on the array */
                  array[ypos][start_pos - jpos] = STRAND_DEFN[i][j];
                }
            }

          /* Re-initialise strand length */
          strand_len = 0;
        }

      /* If plotting coil, mark as single dot */
      if (plot == COIL)
        {
          /* Loop over coil height */
          for (i = 0; i < COIL_HEIGHT; i++)
            {
              /* Get y-position of current point */
              ypos = cpos - COIL_HEIGHT / 2 + i;

              /* Mark this point on the array */
              array[ypos][xpos] = COIL_DEFN[i][0];
            }
        }

      /* For helix, plot appropriate section */
      else if (plot == HELIX)
        {
          /* Plot line in current direction */
          for (i = 0; i < HELIX_HEIGHT; i++)
            {
              /* Get y-position of current point */
              ypos = cpos - HELIX_HEIGHT / 2 + i;

              /* Mark this point on the array */
              array[ypos][xpos] = HELIX_DEFN[i][helix_len];
            }

          /* Increment helix length */
          helix_len++;

          /* If longer than helix definition, reset length */
          if (helix_len > LEN_HELIX_DEFN - 1)
            helix_len = 0;
        }

      /* For strand plot appropriate section */
      else if (plot == STRAND)
        {
          /* Increment strand length */
          strand_len++;
        }

      /* Otherwise, mark this pixel position as blank */
      else
        blank[xpos] = 'Y';

      /* If sec struc plotted was not a helix, re-initialise strand
         length */
      if (plot != HELIX)
        helix_len = 0;

      /* Save current plot */
      last_plot = plot;
    }

  /* Loop over the pixel positions to draw any breaks in structure */
  for (xpos = 0; xpos < width; xpos++)
    {
      /* If have break in structure, draw break */
      if (blank[xpos] == 'Y')
        {
          /* Increment length of blank run */
          run_len++;

          /* If at break-point, draw break line */
          if ((xpos > 0 && blank[xpos - 1] == ' ') ||
              (xpos < width - 1 && blank[xpos + 1] == ' '))
            {
              /* Loop over height of break line */
              for (i = 0; i < BREAK_HEIGHT; i++)
                {
                  /* Get y-position of current point */
                  ypos = cpos - BREAK_HEIGHT / 2 + i;

                  /* Loop over the break length */
                  for (j = 0; j < LEN_BREAK_DEFN; j++)
                    {
                      /* Get position on array */
                      jpos = xpos - LEN_BREAK_DEFN / 2 + j;

                      /* Mark this point on the array if not off end */
                      if (jpos > -1 && jpos < width)
                        array[ypos][jpos] = BREAK_DEFN[i][j];
                    }
                }
            }
        }

      /* If blank run just ended, draw dotted line */
      if (run_len > 0 && (blank[xpos] == ' ' || xpos == width - 1))
        {
          /* Draw dot every 3rd position */
          for (i = 2; i < run_len - 1; i = i + 3)
            {
              /* Get the corresponding position for this dot */
              jpos = xpos - run_len + i;
              if (blank[xpos] == 'Y' && xpos == width - 1)
                jpos++;
              ypos = cpos;

              /* Mark this point on the array if not off end */
              if (jpos > -1 && jpos < width)
                array[ypos][jpos] = '.';
            }
        }

      /* If not blank, then re-initialise blank run count */
      if (blank[xpos] == ' ')
        run_len = 0;
    }

  /* Show the bar representing the UniProt sequence at the top */
  if (plot_domains == TRUE)
    ypos = cpos + WIR_HEIGHT / 2 + WIR_GAP + SEQ_HEIGHT / 2;
  else if (sequence[UNIPROT_SEQ][0] != '\0')
    ypos = cpos + WIR_HEIGHT / 2 + WIR_GAP + 1;
  else
    ypos = cpos;
  cpos = ypos;

  /* Plot line across length of sequence */
  if (sequence[UNIPROT_SEQ][0] != '\0')
    {
      for (xpos = 0; xpos < width; xpos++)
        {
          /* If not a gap in the UniProt sequence, plot bar */
          if (sequence[UNIPROT_SEQ][xpos] != ' ')
            array[ypos][xpos] = '*';
        }
    }

  /* Get the first of this sequence's Pfam domains */
  if (smatch_ptr != NULL)
    domain_ptr = smatch_ptr->first_domain_ptr;
  else
    domain_ptr = NULL;

  /* Loop over all the domains to plot */
  while (plot_domains == TRUE && domain_ptr != NULL)
    {
      /* Get the start- and end-residues for this domain */
      start_res = domain_ptr->start_res;
      end_res = domain_ptr->end_res;

      /* Get type and colour */
      type = domain_ptr->pfam_ptr->type;
      colour = domain_ptr->colour;

      /* Initialise residue counter */
      ires = 0;
      next_j = 0;

      /* Get height of domain bar, according to type */
      if (type == PFAMA || type == COVID19)
        {
          /* Determine the domain height */
          height = PFAMA_HEIGHT;
          len_defn = LEN_PFAMA_DEFN;

          /* Define the starting colour */
          if (params.covid19 == TRUE || type == COVID19)
            {
              /* Get the starting colour for this protein */
              prot_no = -1;
              for (i = 0; i < NCOVID_PROTEINS && prot_no == -1; i++)
                {
                  /* If match found, save the protein number */
                  if (!strcmp(code_ptr->uniprot_acc,COVID_PROTEIN[i]))
                    prot_no = i;
                }

              /* Save the protein number */
              code_ptr->covid_protein = prot_no;

              /* If protein not found, use default colour */
              if (prot_no == -1)
                prot_no = 0;

              /* Get the corresponding starting colour */
              start_col = 3 * prot_no;
              start_col = start_col % NPFAM_COLOURS;
            }

          /* Define the starting colour */
          else if (type == PFAMA)
            start_col = PFAMA_STARTCOL;
        }
      else
        {
          /* Determine the PFAMB domain height */
          height = PFAMB_HEIGHT;
          len_defn = LEN_PFAMB_DEFN;

          /* Define the starting colour */
          start_col = PFAMB_STARTCOL;
        }

      /* Scan line, plotting domain in appropriate position */
      for (xpos = 0; xpos < width && type != UNKNOWN; xpos++)
        {
          /* If not a gap in the UniProt sequence, increment
             residue number */
          if (sequence[UNIPROT_SEQ][xpos] != ' ')
            ires++;

          /* If residue within range of current domain, then
             plot appropriate part of domain diagram */
          if (sequence[UNIPROT_SEQ][xpos] != ' ' &&
              ires >= start_res && ires <= end_res)
            {
              /* Determine which column to plot from */
              if (next_j > len_defn - 1)
                next_j = 0;
              j = next_j;
              next_j++;

              /* Loop to plot the vertical part of the bar */
              for (i = 0; i < height; i++)
                {
                  /* Determine y position in array */
                  ypos = cpos - height / 2 + i;

                  /* DetermonFill in array cell with appropriate pixel
                     definition */
                  if (type == PFAMA || type == COVID19)
                    ch = PFAMA_DEFN[height - i - 1][j];
                  else
                    ch = PFAMB_DEFN[height - i - 1][j];

                  /* Determine colour to use */
                  if (ch == '*')
                    {
                      icol = (start_col + domain_ptr->colour)
                        % (NPFAM_COLOURS - 2);
                      ch = COLOUR_CODE[icol];
                    }

                  /* Store appropriate colour code */
                  array[ypos][xpos] = ch;

                  /* Update map coordinates for this domain */
                  if (domain_ptr->map[XMIN] == -1)
                    {
                      domain_ptr->map[XMIN] = xpos;
                      domain_ptr->map[XMAX] = xpos;
                      domain_ptr->map[YMIN] = ypos;
                      domain_ptr->map[YMAX] = ypos;
                    }
                  else
                    {
                      if (xpos < domain_ptr->map[XMIN])
                        domain_ptr->map[XMIN] = xpos;
                      if (xpos > domain_ptr->map[XMAX])
                        domain_ptr->map[XMAX] = xpos;
                      if (ypos < domain_ptr->map[YMIN])
                        domain_ptr->map[YMIN] = ypos;
                      if (ypos > domain_ptr->map[YMAX])
                        domain_ptr->map[YMAX] = ypos;
                    }
                }
            }
        }

      /* Initialise character positions */
      label_height = 0;
      label_width = 0;
      len = strlen(domain_ptr->pfam_ptr->name);
      domain_width = domain_ptr->map[XMAX] - domain_ptr->map[XMIN] + 1;
      domain_height = domain_ptr->map[YMAX] - domain_ptr->map[YMIN] + 1;

      /* Initialise label array */
      for (irow = 0; irow < LETTER_ROWS; irow++)
        {
          /* Initialise flags indicating which rows are in use */
          used[irow] = FALSE;

          /* Loop over the columns to blank out */
          for (icol = 0; icol < MAX_LABEL_WIDTH; icol++)
            label_array[irow][icol] = ' ';
        }

      /* Loop over the letters of the domain name to write into 
         letter array */
      for (i = 0; i < len; i++)
        {
          /* Get the current letter */
          ch = domain_ptr->pfam_ptr->name[i];

          /* If this is a space, add appropriate number of blanks */
          if (ch == ' ')
            {
              /* Increment label width */
              if (label_width + SPACE_WIDTH < MAX_LABEL_WIDTH)
                label_width = label_width + SPACE_WIDTH;
            }

          /* Otherwise, find the corresponding letter definition */
          else
            {
              /* Find letter */
              iletter = -1;
              for (j = 0; j < nletters && iletter == -1; j++)
                {
                  if (ch == character[j])
                    iletter = j;
                }

              /* If letter found and it can be fitted into letter array,
                 copy across */
              if (iletter > -1 &&
                  (label_width
                   + letter_ptr[iletter].width) < MAX_LABEL_WIDTH)
                {
                  /* Get this letter's width */
                  nrows = letter_ptr[iletter].height;
                  ncols = letter_ptr[iletter].width;

                  /* Loop over the rows of the letter */
                  for (irow = 0; irow < nrows; irow++)
                    for (icol = 0; icol < ncols; icol++)
                      {
                        /* Get position in array */
                        ipos = irow * ncols + icol;

                        /* Update label pixel */
                        label_array[irow][label_width + icol]
                          = letter_ptr[iletter].array[ipos];

                        /* If this is an 'X', update this label row
                           as used */
                        if (letter_ptr[iletter].array[ipos] == 'X')
                          used[irow] = TRUE;
                      }

                  /* Increment current label width */
                  label_width = label_width + ncols;
                }
            }
        }

      /* See if domain's name can be written into domain box */
      if (label_width > 0 && label_width < domain_width)
        {
          /* Determine the x-offset */
          xoffset = (domain_width - label_width) / 2;

          /* Determine range of rows in the label array that contain
             character parts */
          starty = 0;
          endy = LETTER_ROWS;
          for (i = 0; i < LETTER_ROWS && used[i] == FALSE; i++)
            starty = i + 1;
          for (i = LETTER_ROWS - 1; i > -1 && used[i] == FALSE; i--)
            endy = i;
          label_height = endy - starty;

          /* Determine y-offset */
          yoffset = (domain_height - label_height) / 2;

          /* Loop over the rows of the label */
          for (irow = starty; irow < endy; irow++)
            {
              /* Get row in image array */
              ypos = domain_ptr->map[YMAX] - yoffset - irow; 

              /* Loop over the label columns */
              for (icol = 0; icol < label_width; icol++)
                {
                  /* If label pixel belongs to a character, add to
                     image array */
                  if (label_array[irow][icol] == 'X')
                    {
                      /* Get column in image array */
                      xpos = domain_ptr->map[XMIN] + xoffset + icol;

                      /* Store pixel */
                      if (domain_ptr->pfam_ptr->type == PFAMA ||
                          domain_ptr->pfam_ptr->type == COVID19)
                        array[ypos][xpos] = '$';
                      else
                        array[ypos][xpos] = '#';
                    }
                }
            }
        }

      /* Get the next domain */
      domain_ptr = domain_ptr->next_domain_ptr;
    }

  /* Get position for black crosses */
  cpos = WIR_HEIGHT + 1;
  if (pdb_ptr->have_cath == TRUE)
    cpos = cpos + CATH_HEIGHT;

  /* For COVID-19 page, save row at which to split the domain image
     from the secondary structure image */
  if (params.covid19 == TRUE && sequence[UNIPROT_SEQ][0] != '\0')
    image_split = cpos + 2;

  /* Compare the two sequences and mark any differences */
  if (sequence[UNIPROT_SEQ][0] != '\0')
    {
      for (xpos = 0; xpos < width; xpos++)
        {
          /* If the sequences differ, show cross */
          if (blank[xpos] == ' ' && sequence[UNIPROT_SEQ][xpos] != ' ' &&
              sequence[UNIPROT_SEQ][xpos] != sequence[PDB_SEQ][xpos])
            {
              /* Plot the centre of the cross */
              array[cpos][xpos] = '#';
          
              /* Plot cross lines */
              if (xpos > 0)
                {
                  array[cpos + 1][xpos - 1] = '#';
                  array[cpos - 1][xpos - 1] = '#';
                }
              if (xpos < width - 2)
                {
                  array[cpos + 1][xpos + 1] = '#';
                  array[cpos - 1][xpos + 1] = '#';
                }

              /* Increment count of mismatched residues */
              nmismatch++;
            }
        }
    }

  /* Determine mid-point for CATH domains, if any */
  cpos = CATH_HEIGHT / 2;

  /* Get the first of this sequence's CATH domains */
  cath_domain_ptr = chain_ptr->first_cath_domain_ptr;

  /* Loop over all the domains to plot */
  while (cath_domain_ptr != NULL)
    {
      /* Get the start- and end-residues for this domain */
      start_res = cath_domain_ptr->start_res;
      end_res = cath_domain_ptr->end_res;

      /* Get the domain colour */
      domain = cath_domain_ptr->domain_no;
      colour = get_domain_colour(domain);

      /* Initialise residue counter */
      ires = 0;

      /* Get height of domain bar */
      height = CATH_HEIGHT;

      /* Scan line, plotting domain in appropriate position */
      for (xpos = 0; xpos < width; xpos++)
        {
          /* If not a gap in the UniProt sequence, increment
             residue number */
          if (sequence[PDB_SEQ][xpos] != ' ' && blank[xpos] != 'Y')
            ires++;

          /* If residue within range of current domain, then plot
             appropriate part of domain diagram */
          if (sequence[PDB_SEQ][xpos] != ' ' && blank[xpos] != 'Y' &&
              ires >= start_res && ires <= end_res)
            {
              /* Determine which column to plot from */
              if (ires - start_res > -1 && ires - start_res < 3)
                j = ires - start_res;
              else if ((end_res - ires) < 3)
                j = 6 - (end_res - ires);
              else
                j = 3;

              /* Loop to plot the vertical part of the bar */
              for (i = 0; i < height; i++)
                {
                  /* Determine y position in array */
                  ypos = cpos - height / 2 + i;

                  /* Fill in array cell with appropriate pixel
                     definition */
                  ch = CATH_DEFN[i][j];

                  /* Determine colour to use */
                  if (ch == '*')
                    ch = COLOUR_CODE[colour];

                  /* Store appropriate colour code */
                  array[ypos][xpos] = ch;

                  /* Update map coordinates for this domain */
                  if (cath_domain_ptr->map[XMIN] == -1)
                    {
                      cath_domain_ptr->map[XMIN] = xpos;
                      cath_domain_ptr->map[XMAX] = xpos;
                      cath_domain_ptr->map[YMIN] = ypos;
                      cath_domain_ptr->map[YMAX] = ypos;
                    }
                  else
                    {
                      if (xpos < cath_domain_ptr->map[XMIN])
                        cath_domain_ptr->map[XMIN] = xpos;
                      if (xpos > cath_domain_ptr->map[XMAX])
                        cath_domain_ptr->map[XMAX] = xpos;
                      if (ypos < cath_domain_ptr->map[YMIN])
                        cath_domain_ptr->map[YMIN] = ypos;
                      if (ypos > cath_domain_ptr->map[YMAX])
                        cath_domain_ptr->map[YMAX] = ypos;
                    }
                }
            }
        }

      /* Get the next domain */
      cath_domain_ptr = cath_domain_ptr->next_domain_ptr;
    }

  /* Generate the image */
  generate_image(array,width,pdb_ptr,chain_ptr,icode,plot_no,image_height,
                 pdbsum_dir,image_split,res_range,file_name_ptr,params);

  /* Free up memory taken up by image array */
  for (ypos = 0; ypos < PIC_HEIGHT; ypos++)
    free(array[ypos]);
  free(blank);

  /* Store this chain's image height */
  chain_ptr->image_height = image_height;

  /* Return count of mismatched residues */
  return(nmismatch);
}
/***********************************************************************

compare_seqs  -  Compare all PDB chains for current UniProt code
                 against the UniProt sequence

***********************************************************************/

int compare_seqs(struct code *code_ptr,char *uniprot_fasta,
                 char *fasta_id,int seq_len,int run_no,char *character,
                 struct letter *letter_ptr,int nletters,
                 FILE *drfile_ptr,struct pfamlig **fst_pfamlig_ptr,
                 struct pfamlig **lst_pfamlig_ptr,int *npflig,
                 struct c_file_data file_name_ptr,
                 struct parameters params)
{
  char file_name[FILENAME_LEN], seq_name[FILENAME_LEN];
  char res_range[MAX_SEQLEN], sec_struc[MAX_SEQLEN], sequence[2][MAX_SEQLEN];
  char diffs[MAX_SEQLEN], pdb_resno[2][MAX_SEQLEN];
  char pdbsum_dir[FILENAME_LEN + 1], template[FILENAME_LEN + 1];

  int align_len, cath_domain[MAX_SEQLEN], i, icode, nligands, nscores;
  int dir_type, ndiffs, nligand, nmetal, nmismatch, nox, npfamlig, nresid;
  int dos, pdbsum_type;
  int ok;

  struct c_score *first_c_score_ptr;
  struct chain *chain_ptr;
  struct clist *clist_ptr;
  struct r_ligand *first_r_ligand_ptr;
  struct pdb *pdb_ptr;
  struct pfamlig *first_pfamlig_ptr, *last_pfamlig_ptr;
  struct pfamlig *start_pfamlig_ptr;
  struct r_rasdef *first_r_rasdef_ptr;
  struct residue *first_residue_ptr;
  struct smatch *smatch_ptr;

  /* Initialise variables */
  dos = file_name_ptr.dos;
  nligands = nox = 0;
  npfamlig = *npflig;
  nscores = 0;
  pdbsum_type = params.pdbsum_type;
  res_range[0] = '\0';

  /* Store this UniProt sequence length */
  code_ptr->seq_len = seq_len;

  /* Initialise pointers */
  first_c_score_ptr = NULL;
  first_pfamlig_ptr = *fst_pfamlig_ptr;
  last_pfamlig_ptr = *lst_pfamlig_ptr;

  /* Get pointer to the first clist entry for this UniProt code */
  clist_ptr = code_ptr->first_clist_ptr;

  /* Loop through all entries to compare PDB chain sequence against
     UniProt code sequence */
  while (clist_ptr != NULL)
    {
      /* Initialise PDB and UniProt sequences */
      sequence[UNIPROT_SEQ][0] = '\0';
      sequence[PDB_SEQ][0] = '\0';

      /* Initialise first residue */
      first_residue_ptr = NULL;
      nligands = 0;

      /* Get the corresponding chain and PDB entry */
      chain_ptr = clist_ptr->chain_ptr;
      icode = clist_ptr->icode;
      pdb_ptr = chain_ptr->pdb_ptr;
      smatch_ptr = chain_ptr->smatch_ptr[icode];

      /* Form the name of the PDBsum directory */
      if (params.pdb_type == PDBSUM1)
        {
          strcpy(pdbsum_dir,"../newsec/");
          dir_type = PDBSUM1;

          /* Form the name of the UniProt .fasta file */
          sprintf(uniprot_fasta,"%s.fasta",code_ptr->uniprot_acc);
        }
      else
        {
          dir_type = find_pdbsum_dir(pdb_ptr->pdb_code,
                                     file_name_ptr.pdbsum_template,
                                     pdbsum_type,pdbsum_dir);
          sprintf(uniprot_fasta,"%s.fasta",code_ptr->uniprot_acc);
        }
      
      /* If have PDBsum directory, then continue */
      if (dir_type != NOT_FOUND)
        {
          /* Show entry being processed */
          printf("    PDB entry %s\n",pdb_ptr->pdb_code);
          fflush(stdout);

          /* Check all this PDB code's Pfam domains to assign colours to
             each */
          if (pdb_ptr->coloured == FALSE)
            colour_pfam_domains(pdb_ptr);

          /* Process this chain */
          if (chain_ptr->ncodes > 0)
            {
              /* Initialise variables */
              nmismatch = 0;

              /* Form name of seqXXX.fasta file for this chain */
              sprintf(seq_name,"seq%3.3d.fasta",chain_ptr->chain_no);
              strcpy(file_name,pdbsum_dir);
              strcat(file_name,seq_name);

              /* Increment run number (for testing only) */
#ifdef TEST
              run_no++;
#endif

              /* If not dummy UniProt code, run FASTA to compare the two
                 sequences */
              if (strcmp(code_ptr->uniprot_acc,"NONE"))
                {
                  /* For PDBsum1, we already have the alignment */
                  if (params.pdb_type == PDBSUM1)
                    {
                      /* Retrieve the alignment */
                      strcpy(sequence[PDB_SEQ],
                             chain_ptr->alignment[PDB_SEQ]);
                      strcpy(sequence[UNIPROT_SEQ],
                             chain_ptr->alignment[UNIPROT_SEQ]);
                    }

                  /* Otherwise, run FASTA to get it */
                  else
                    {
                      /* Run FASTA to compare seqXXX.fasta file vs
                         the UniProt sequence */
                      run_fasta(file_name_ptr.exe_name[FASTA],
                                file_name,uniprot_fasta,run_no,dos);

                      /* Extract the sequence alignment from the
                         fasta.out file */
                      read_fasta_out(chain_ptr,fasta_id,
                                     sequence,code_ptr->uniprot_acc,run_no);
                    }
                }

              /* Otherwise, for dummy sequence, read in the sequuence from
                 the seqXXX.fasta file */
              else
                read_sequence(file_name,sequence);

              /* Read in the residue data for this chain from the PDB entry's
                 resdata.dat file in its PDBsum directory */
              nresid = get_resdata(pdb_ptr->pdb_code,chain_ptr->chain_id,
                                   chain_ptr,icode,pdbsum_dir,
                                   &first_residue_ptr,&nox);


              /* If this is an Alpha Fold structure, pick up the confidence
                 scores from the AF_scores.dat file */
              if (params.pdb_type == ALPHA_FOLD_PDB)
                nscores
                  = get_c_scores(pdb_ptr->pdb_code,pdbsum_dir,
                                 &first_c_score_ptr);

              /* Match the sequence to the residue records read in from
                 resdata.dat */
              ok = match_to_residues(chain_ptr,icode,first_residue_ptr,nresid,
                                     first_c_score_ptr,nscores,nox,
                                     sequence,sec_struc,res_range,cath_domain,
                                     pdb_resno,code_ptr->uniprot_acc);

              /* If all OK, get ligand interactions and plot the mini wiring
                 diagram */
              if (ok == TRUE)
                {
                  /* Write out start of DR line to pdbsum.dr */
                  fprintf(drfile_ptr,"%s",code_ptr->uniprot_acc);
                  for (i = 0; i < 10 - strlen(code_ptr->uniprot_acc); i++)
                    fprintf(drfile_ptr," ");
                  fprintf(drfile_ptr,"DR   PDBsum; %s;",pdb_ptr->pdb_code);

                  /* Show chain */
                  if (chain_ptr->chain_id != ' ')
                    fprintf(drfile_ptr," %c",chain_ptr->chain_id);
                  else
                    fprintf(drfile_ptr," @");

                  /* If have copy chains, add these */
                  //if (chain_list[0] != '\0')
                  //  fprintf(drfile_ptr,"%s",chain_list);

                  /* Initialise count of interacting ligands */
                  nligand = 0;

                  /* Get the PDBsum directory template */
                  strcpy(template,
                         file_name_ptr.pdbsum_template[pdbsum_type]);
                  if (params.pdb_type == PDBSUM1)
                    strcpy(template,"..");

                  /* Read in the definitions from the rasmol.dfn file */
                  r_get_dfn_data(pdb_ptr->pdb_code,params.pdb_type,
                                 template,&first_r_rasdef_ptr,
                                 &first_r_ligand_ptr,&nligand,&nmetal);
      
                  /* If have ligands, find which interact with the current
                     chain and each of its Pfam domains */
                  if (nligand > 0 || nmetal > 0)
                    nligand
                      = identify_ligands(pdb_ptr->pdb_code,chain_ptr,
                                         code_ptr,first_r_rasdef_ptr,
                                         first_residue_ptr,
                                         &first_pfamlig_ptr,
                                         &last_pfamlig_ptr,
                                         &start_pfamlig_ptr,&npfamlig,
                                         pdbsum_dir);

                  /* If any interacting ligands found, write to pdbsum.dr
                     file */
                  if (nligand > 0)
                    {
                      /* Flag single ligand record for each ligand of
                         interest */
                      flag_ligands(first_r_rasdef_ptr);

                      /* Combine ligands of same name */
                      unique_ligands(first_r_rasdef_ptr);

                      /* Write out the ligands that interact with this chain
                         to the pdbsum.dr file, using interaction data
                         from the grow.out file */
                      write_ligands(first_r_rasdef_ptr,drfile_ptr);

                      /* Identify ligands that bind to multiple Pfam
                         domains */
                      npfamlig
                        = multiple_domains(start_pfamlig_ptr,
                                           &first_pfamlig_ptr,
                                           &last_pfamlig_ptr,npfamlig);
                    }

                  /* Close off the DR line */
                  fprintf(drfile_ptr,".\n");
                  fflush(drfile_ptr);

                  /* Get the CATH domain assignments from the database.dat
                     file */
                  if (params.pdb_type != PDBSUM1)
                    identify_cath_domains(pdb_ptr,chain_ptr,pdbsum_dir);

                  /* If have any CATH domains for this chain, renumber
                     domains for each residue */
                  if (chain_ptr->ncath_dom > 0)
                    match_cath_domains(chain_ptr,cath_domain,
                                       sequence[PDB_SEQ],
                                       strlen(sequence[PDB_SEQ]));

                  /* Store the details of this sequence match in the
                     corresponding smatch record */
                  if (sequence[UNIPROT_SEQ][0] != '\0')
                    store_string(&(smatch_ptr->uniprot_seq),
                                 sequence[UNIPROT_SEQ]);
                  if (sequence[PDB_SEQ][0] != '\0')
                    store_string(&(smatch_ptr->pdb_seq),sequence[PDB_SEQ]);
                  store_string(&(smatch_ptr->pdb_hyphens),pdb_resno[HYPHENS]);
                  store_string(&(smatch_ptr->pdb_resno),pdb_resno[NUMBERS]);

                  /* For PDBsum1, redefine current directory as PDBsum
                     directory */
                  if (params.pdb_type == PDBSUM1)
                    strcpy(pdbsum_dir,"./");

                  /* Plot the mini wiring diagram */
                  nmismatch
                    = plot_miniwir(chain_ptr,icode,sequence,sec_struc,
                                   res_range,pdb_ptr,character,letter_ptr,
                                   nletters,pdbsum_dir,file_name_ptr,params);
                }

              /* Free up memory taken up by residue and ligand data */
              free_residue(&first_residue_ptr);

              /* Get the alignment length */
              align_len = strlen(sequence[PDB_SEQ]);
              ndiffs = 0;

              /* Form string giving differences */
              for (i = 0; i < align_len &&
                     sequence[UNIPROT_SEQ][0] != '\0'; i++)
                {
                  /* If sequences differ, then mark differences line */
                  if (sequence[UNIPROT_SEQ][i] != ' ' &&
                      sequence[PDB_SEQ][i] != ' ' &&
                      sequence[UNIPROT_SEQ][i] != sequence[PDB_SEQ][i])
                    {
                      /* Put marker at this alignment position */
                      diffs[i] = 'X';
                      ndiffs++;
                    }
                  else
                    diffs[i] = '.';
                }
              diffs[align_len] = '\0';

              /* If have any differences, then store */
              if (ndiffs > 0)
                store_string(&(smatch_ptr->diffs),diffs);

              /* Store number of mismatched residues */
              smatch_ptr->nmismatch = nmismatch;
            }
        }
      else
        printf("*** Warning. Unable to locate PDBsum directory for %s\n",
               pdb_ptr->pdb_code);

      /* Get pointer to the next entry in the linked list */
      clist_ptr = clist_ptr->next_clist_ptr;
    }

  /* Return current pointers */
  *fst_pfamlig_ptr = first_pfamlig_ptr;
  *lst_pfamlig_ptr = last_pfamlig_ptr;

  /* Return count of ligand to Pfam domain links */
  *npflig = npfamlig;

  /* Return current FASTA run number */
  return(run_no);
}
/***********************************************************************

write_pfam_details  -  Write Pfam domain details to the uniplot.dat
                       file

***********************************************************************/

void write_pfam_details(FILE *file_ptr,struct chain *chain_ptr,
                        struct smatch *smatch_ptr,char *sub1)
{
  int ncopy, ndomain, type;
  int done;

  struct domain *domain_ptr, *other_domain_ptr;
  struct pfam *pfam_ptr;

  /* Initialise variables */
  ndomain = 0;

  /* Sort the domains in order of start residue */
  order_domains(&(smatch_ptr->first_domain_ptr));

  /* Get the first of this sequence's Pfam domains */
  domain_ptr = smatch_ptr->first_domain_ptr;

  /* Identify domain number from list of current
     chain's domain */
  while (domain_ptr != NULL)
    {
      /* Get the corresponding Pfam pointer */
      pfam_ptr = domain_ptr->pfam_ptr;
      done = FALSE;

      /* If this is a Pfam domain determine whether to write out its
         details */
      if (pfam_ptr != NULL)
        {
          /* Get pointer to first domain again */
          other_domain_ptr = smatch_ptr->first_domain_ptr;
          done = FALSE;

          /* Loop to see whether we have already processed this Pfam
             domain */
          while (other_domain_ptr != domain_ptr &&
                 other_domain_ptr != NULL && done == FALSE)
            {
              /* If this is the same Pfam code, then we have already
                 done it */
              if (other_domain_ptr->pfam_ptr == domain_ptr->pfam_ptr)
                done = TRUE;

              /* Get the next domain */
              other_domain_ptr = other_domain_ptr->next_domain_ptr;
            }
        }

      /* If this is a valid Pfam domain which hasn't yet been processed,
         then process now */
      if (pfam_ptr != NULL && done == FALSE)
        {
          /* Increment domain counter */
          ndomain++;

          /* Get domain type */
          type = pfam_ptr->type;

          /* Write out the domain id */
          fprintf(file_ptr,"UPLOT_PFAM_ID%s[%d] %s\n",sub1,
                  ndomain,pfam_ptr->pfam_id);

          /* Write out PfamA or PfamB flag */
          if (type == PFAMA)
            fprintf(file_ptr,"UPLOT_PFAMA%s[%d] TRUE\n",sub1,
                    ndomain);
          else if (type == COVID19)
            fprintf(file_ptr,"UPLOT_COVID19%s[%d] TRUE\n",sub1,
                    ndomain);
          else
            fprintf(file_ptr,"UPLOT_PFAMB%s[%d] TRUE\n",sub1,
                    ndomain);
                      
          /* Write out the domain identifier and description */
          fprintf(file_ptr,"UPLOT_DOMAIN_NAME%s[%d] %s\n",sub1,
                  ndomain,pfam_ptr->name);
          fprintf(file_ptr,"UPLOT_DOMAIN_DESC%s[%d] %s\n",sub1,
                  ndomain,pfam_ptr->description);

          /* Initialise count of copies of this domain */
          ncopy = 0;

          /* Get pointer to first domain again */
          other_domain_ptr = smatch_ptr->first_domain_ptr;

          /* Loop to write out start and end residues, domain number
             and colour */
          while (other_domain_ptr != NULL)
            {
              /* If this is the same Pfam code, then write out
                 required details */
              if (other_domain_ptr->pfam_ptr == domain_ptr->pfam_ptr)
                {
                  /* Increment copy counter */
                  ncopy++;

                  /* Write out the start and end residues */
                  fprintf(file_ptr,"UPLOT_PFAM_START%s[%d][%d] %d\n",
                          sub1,ndomain,ncopy,
                          other_domain_ptr->start_res);
                  fprintf(file_ptr,"UPLOT_PFAM_END%s[%d][%d] %d\n",
                          sub1,ndomain,ncopy,
                          other_domain_ptr->end_res);

                  /* Write out domain number */
                  fprintf(file_ptr,"UPLOT_PFAM_DOMAIN_NO%s[%d][%d] "
                          "%d\n",sub1,ndomain,ncopy,
                          other_domain_ptr->domain_no);

                  /* Write out domain colour */
                  fprintf(file_ptr,"UPLOT_PFAM_DOMAIN_COLOUR%s[%d][%d] "
                          "%d\n",sub1,ndomain,ncopy,
                          other_domain_ptr->colour);
                }

              /* Get the next domain */
              other_domain_ptr = other_domain_ptr->next_domain_ptr;
            }

          /* Write out the number of copies of this domain within the
             current chain */
          fprintf(file_ptr,"UPLOT_PFAM_NCOPY%s[%d] %d\n",sub1,
                  ndomain,ncopy);
        }

      /* Get the next domain */
      domain_ptr = domain_ptr->next_domain_ptr;
    }

  /* Write out number of Pfam domains, if any */
  if (ndomain > 0)
    fprintf(file_ptr,"UPLOT_NPFAM%s %d\n",sub1,ndomain);
}
/***********************************************************************

write_domain_details  -  Write Pfam domain details to the uniplot.dat
                         file for each image

***********************************************************************/

void write_domain_details(FILE *file_ptr,struct chain *chain_ptr,
                          int icode,struct smatch *smatch_ptr,char *sub1,
                          int *have_pfam)
{
  char flag[8], id[20], name[128], description[512];

  int end_pos, i, image, image_height, image_width, loop, map[4];
  int ndomain, start_pos, type;
  int in_image;

  struct cath *cath_ptr;
  struct domain *domain_ptr;
  struct pfam *pfam_ptr;

  /* Initialise variables */
  image_height = chain_ptr->image_height;
  start_pos = 0;

  /* Loop over the images generated for this chain */
  for (image = 0; image < chain_ptr->nimages[icode]; image++)
    {
      /* Initialise count of domains in this image */
      ndomain = 0;

      /* If this is the last image, adjust width */
      if (image == chain_ptr->nimages[icode] - 1)
        image_width = chain_ptr->last_width;
      else
        image_width = chain_ptr->image_width;

      /* Get the last pixel position (relative to original unsplit
         image) */
      end_pos = start_pos + image_width - 1;

      /* Write coords for area of image that will click for sequence
         alignment */
      map[XMIN] = 0;
      map[XMAX] = image_width - 1;
      map[YMIN] = 0;
      map[YMAX] = image_height - 1;
      for (i = 0; i < 4; i++)
        fprintf(file_ptr,"UNIPLOT_SST_%s%s[%d] %d\n",
                MAP_TYPE[i],sub1,image + 1,map[i]);

      /* Loop over the two types of domain (Pfam and CATH) */
      for (loop = 0; loop < 2; loop++)
        {
          /* Get the first of this chain's Pfam or CATH domains */
          if (loop == PFAM)
            domain_ptr = smatch_ptr->first_domain_ptr;
          else
            domain_ptr = chain_ptr->first_cath_domain_ptr;

          /* Loop over all the domains to determine their colours */
          while (domain_ptr != NULL)
            {
              /* For Pfam domains check domain type */
              if (loop == PFAM)
                {
                  /* Get whether PfamA or PfamB */
                  type = domain_ptr->pfam_ptr->type;

                  /* Update domain type flag */
                  if (type != UNKNOWN)
                    have_pfam[type] = TRUE;
                }

              /* Initialise flag */
              in_image = FALSE;

              /* See if start of this domain appears in the current
                 image */
              if (domain_ptr->map[XMIN] >= start_pos &&
                  domain_ptr->map[XMIN] <= end_pos)
                {
                 /* Store map start- and end-points */
                  map[XMIN] = domain_ptr->map[XMIN] - start_pos;
                  if (domain_ptr->map[XMAX] > end_pos)
                    map[XMAX] = end_pos - start_pos;
                  else
                    map[XMAX] = domain_ptr->map[XMAX] - start_pos;

                  /* Set flag */
                  in_image = TRUE;
                }

              /* See if end of this domain appears in the current
                 image */
              else if ((domain_ptr->map[XMAX] >= start_pos &&
                        domain_ptr->map[XMAX] <= end_pos))
                {
                  /* Store map start- and end-points */
                  map[XMAX] = domain_ptr->map[XMAX] - start_pos;
                  if (domain_ptr->map[XMIN] < start_pos)
                    map[XMIN] = 0;
                  else
                    map[XMIN] = domain_ptr->map[XMIN] - start_pos;

                  /* Set flag */
                  in_image = TRUE;
                }

              /* Otherwise, check whether domain spans full length of
                 image */
              else if ((domain_ptr->map[XMIN] <= start_pos &&
                        domain_ptr->map[XMAX] >= end_pos))
                {
                  /* Store map start- and end-points */
                  map[XMIN] = 0;
                  map[XMAX] = image_width - 1;

                  /* Set flag */
                  in_image = TRUE;
                }

              /* If domain is in image, write out the details */
              if (in_image == TRUE)
                {
                  /* Get the y-range for domain in image */
                  map[YMAX] = image_height - domain_ptr->map[YMIN] - 1;
                  map[YMIN] = image_height - domain_ptr->map[YMAX] - 1;

                  /* Increment count of domains in this image */
                  ndomain++;

                  /* Get the corresponding domain details */
                  if (loop == PFAM)
                    {
                      /* Get Pfam record */
                      pfam_ptr = domain_ptr->pfam_ptr;

                      /* Store id, name and description */
                      strcpy(id,pfam_ptr->pfam_id);
                      strcpy(name,pfam_ptr->name);
                      strcpy(description,pfam_ptr->description);

                      /* Determine whether this is a Pfam id */
                      if (!strncmp(id,"PF",2) || !strncmp(id,"PB",2))
                        strcpy(flag,"TRUE");
                      else
                        strcpy(flag,"FALSE");

                      /* Write out the domain type */
                      fprintf(file_ptr,"UNIPLOT_PFAM%s[%d][%d] %s\n",
                              sub1,image + 1,ndomain,flag);

                      /* If this is a PfamB domainh, write out flag */
                      if (pfam_ptr->type == PFAMB)
                        fprintf(file_ptr,
                                "UNIPLOT_PFAMB%s[%d][%d] TRUE\n",
                                sub1,image + 1,ndomain);
                    }
                  else if (domain_ptr->cath_ptr != NULL)
                    {
                      /* Get CATH record */
                      cath_ptr = domain_ptr->cath_ptr;

                      /* Store id, name and description */
                      strcpy(id,cath_ptr->cath_code);
                      strcpy(name,cath_ptr->class);
                      strcpy(description,cath_ptr->architecture);

                      /* Write out the domain type */
                      fprintf(file_ptr,"UNIPLOT_CATH%s[%d][%d] TRUE\n",
                          sub1,image + 1,ndomain);
                    }

                  /* Write out the domain identifier and description */
                  fprintf(file_ptr,"UNIPLOT_DOMAIN_ID%s[%d][%d] %s\n",
                          sub1,image + 1,ndomain,id);
                  fprintf(file_ptr,"UNIPLOT_DOMAIN_NAME%s[%d][%d] %s\n",
                          sub1,image + 1,ndomain,name);
                  fprintf(file_ptr,"UNIPLOT_DOMAIN_DESC%s[%d][%d] %s\n",
                          sub1,image + 1,ndomain,description);

                  /* Write out the map range for this domain */
                  for (i = 0; i < 4; i++)
                    fprintf(file_ptr,"UNIPLOT_MAP_%s%s[%d][%d] %d\n",
                            MAP_TYPE[i],sub1,image + 1,ndomain,
                            map[i]);
                }

              /* Get the next domain */
              domain_ptr = domain_ptr->next_domain_ptr;
            }
        }

      /* Write out the number of domains in current image */
      fprintf(file_ptr,"UNIPLOT_NDOMAINS%s[%d] %d\n",sub1,
              image + 1,ndomain);

      /* Update next start-position */
      start_pos = end_pos + 1;
    }
}
/***********************************************************************

add_name  -  Add current protein name to list

***********************************************************************/

int add_name(char *list,char *name,int nprot)
{
  int len, lentot;


  /* Get the total length if this name is added to the total name */
  len = strlen(name);
  lentot = len + strlen(list) + 2;

  /* If OK to add, then do so */
  if (lentot < MAX_LABEL_WIDTH)
    {
      /* Check that name not already in the list */
      if (strstr(list,name) == NULL)
        {
          /* If list not empty, then append comma */
          if (list[0] != '\0')
            strcat(list,", ");
          
          /* Append name */
          strcat(list,name);

          /* Increment protein count */
          nprot++;
        }
    }

  /* Return number of unique proteins in the list */
  return(nprot);
}
/***********************************************************************

sort_ligands  -  Sort the ligands by their name using a simple bubble
                 sort

***********************************************************************/

void sort_ligands(struct chain *chain_ptr)
{
  int swap, swapped;

  struct liglist *liglist_ptr, *first_liglist_ptr;
  struct liglist *next_liglist_ptr, *prev_liglist_ptr, *swap_liglist_ptr;
  struct r_rasdef *r_rasdef_ptr, *next_r_rasdef_ptr;

  /* Get pointer to the first record */
  first_liglist_ptr = chain_ptr->first_liglist_ptr;
  swapped = TRUE;

  /* Get pointer to first liglist record */
  liglist_ptr = first_liglist_ptr;

  /* Loop until Pfam domains in required order */
  while (swapped == TRUE)
    {
      /* Initialise swapped flag */
      swapped = FALSE;
          
      /* Get pointer to first liglist record */
      liglist_ptr = first_liglist_ptr;
      prev_liglist_ptr = NULL;

      /* Loop over all the records, adjusting their order where
         necessary */
      while (liglist_ptr != NULL)
        {
          /* Get pointer to the next record */
          next_liglist_ptr = liglist_ptr->next_liglist_ptr;

          /* Get the corresponding rasdef records */
          r_rasdef_ptr = liglist_ptr->r_rasdef_ptr;
          if (next_liglist_ptr == NULL)
            next_r_rasdef_ptr = NULL;
          else
            next_r_rasdef_ptr = next_liglist_ptr->r_rasdef_ptr;

          /* Determine whether current two records are in the right
             order according to ligand name */
          swap = FALSE;
          if (next_r_rasdef_ptr != NULL)
            {
              if (strcmp(r_rasdef_ptr->ligand_name,
                         next_r_rasdef_ptr->ligand_name) > 0)
                swap = TRUE;
              else if (!strcmp(r_rasdef_ptr->ligand_name,
                               next_r_rasdef_ptr->ligand_name) &&
                       r_rasdef_ptr->chain > next_r_rasdef_ptr->chain)
                swap = TRUE;
              else if (!strcmp(r_rasdef_ptr->ligand_name,
                               next_r_rasdef_ptr->ligand_name) &&
                       r_rasdef_ptr->chain == next_r_rasdef_ptr->chain &&
                       strcmp(r_rasdef_ptr->res_num,
                              next_r_rasdef_ptr->res_num) > 0)
                swap = TRUE;
            }

          /* If need to swap, perform the swap */
          if (swap == TRUE)
            {
              /* If this is the first liglist record, adjust first pointer */
              if (liglist_ptr == first_liglist_ptr)
                first_liglist_ptr = next_liglist_ptr;

              /* Otherwise, get the previous liglist pointer to point
                 to next one */
              else
                prev_liglist_ptr->next_liglist_ptr = next_liglist_ptr;

              /* Swap pointers within the records */
              liglist_ptr->next_liglist_ptr
                = next_liglist_ptr->next_liglist_ptr;
              next_liglist_ptr->next_liglist_ptr = liglist_ptr;

              /* Swap the pointers themselves */
              swap_liglist_ptr = liglist_ptr;
              liglist_ptr = next_liglist_ptr;
              next_liglist_ptr = swap_liglist_ptr;

              /* Set flag that swap made */
              swapped = TRUE;
            }

          /* Save current pointer */
          prev_liglist_ptr = liglist_ptr;

          /* Get the next liglist record */
          liglist_ptr = next_liglist_ptr;
        }

      /* Save the first list pointer */
      chain_ptr->first_liglist_ptr = first_liglist_ptr;
    }
}
/***********************************************************************

find_ligand  -  Identify the current ligand

***********************************************************************/

struct ligand *find_ligand(struct chain *chain_ptr,char *ligand_name)
{
  struct ligand *ligand_ptr, *found_ligand_ptr;


  /* Initialise pointer */
  found_ligand_ptr = NULL;

  /* Get the first of this chain's ligand definitions */
  ligand_ptr = chain_ptr->first_ligand_ptr;

  /* Loop over the ligands to find the one we're after */
  while (ligand_ptr != NULL && found_ligand_ptr == NULL)
    {
      /* If this is the right ligand, then save pointer */
      if (!strcmp(ligand_ptr->ligand_name,ligand_name))
        found_ligand_ptr = ligand_ptr;

      /* Get next ligand record */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Return the found ligand */
  return(found_ligand_ptr);
}
/***********************************************************************

check_for_aa  -  Check whether any of the Het Groups are one of the
                 standard amino acid codes

***********************************************************************/

void check_for_aa(char *res_name,char *res_desc)
{
  int iamino, namino;

  /* Initialise variables */
  res_desc[0] = '\0';

  /* Get the number of amino acids */
  namino = strlen(AMINO);

  /* Loop through all the amino acids until we find the right one */
  for (iamino = 0; iamino < namino && res_desc[0] == '\0'; iamino++)
    {
      /* If this is right code, save the description */
      if (!strcmp(res_name,AMINO_CODE[iamino]))
        {
          strcpy(res_desc,"- ");
          strcat(res_desc,AMINO_NAME[iamino]);
        }
    }
}
/***********************************************************************

binary_hetdrug_search  -  Search for the given Het Group

***********************************************************************/

struct hetdrug *binary_hetdrug_search(struct hetdrug **hetdrug_index_ptr,
                                      int nhetdrug,char *het_name)
{
  int bottom_point, ipoint, new_point, top_point;
  int hetdrug, found;

  struct hetdrug *hetdrug_ptr;

  /* Make starting point half-way into the index */
  ipoint = nhetdrug / 2;
  bottom_point = 0;
  top_point = nhetdrug - 1;

  /* Initialise variables */
  hetdrug = FALSE;
  found = FALSE;
  hetdrug_ptr = NULL;
  if (nhetdrug == 0)
    hetdrug = TRUE;

  /* Loop until required array position found */
  while (hetdrug == FALSE)
    {
      /* Get the associated pointer */
      hetdrug_ptr = hetdrug_index_ptr[ipoint];

      /* Check if array value is equal to the value we want */
      if (!strcmp(het_name,hetdrug_ptr->het_name))
        {
          /* Set flag to indicate that match found */
          hetdrug = TRUE;
          found = TRUE;
        }

      /* If array value is higher than the one we want, then check
         whether value of entry before is lower */
      else if (strcmp(het_name,hetdrug_ptr->het_name) < 0)
        {
          /* If this is the first dictionary in the sorted list, then
             can't go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are hetdrug */
              hetdrug = TRUE;
            }

          /* Otherwise, check whether the previous entry's het_name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              hetdrug_ptr = hetdrug_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(het_name,hetdrug_ptr->het_name) > 0)
                hetdrug = TRUE;

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
          if (ipoint >= nhetdrug - 1)
            {
              /* Set flag as hetdrug */
              hetdrug = TRUE;
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
    hetdrug_ptr = NULL;
  return(hetdrug_ptr);
}
/***********************************************************************

write_int_ligands  -  Write out list of ligands interacting with the 
                      given chain

***********************************************************************/

int write_int_ligands(FILE *file_ptr,struct chain *chain_ptr,char *sub1,
                      struct hetdrug **hetdrug_index_ptr,int nhetdrug)
{
  char description[LINELEN], desc[LINELEN], last_ligname[LINELEN];
  char ligand_name[LINELEN], string[10];

  int len, level1, level2, ncopy, nhets, nligs;

  struct het *het_ptr;
  struct hetdrug *hetdrug_ptr;
  struct ligand *ligand_ptr;
  struct liglist *liglist_ptr;
  struct r_rasdef *r_rasdef_ptr;

  /* Initialise variables */
  last_ligname[0] = '\0';
  ncopy = nligs = 0;

  /* Sort the ligands attached to this chain */
  sort_ligands(chain_ptr);

  /* Get the first of this chain's ligand definitions */
  ligand_ptr = chain_ptr->first_ligand_ptr;

  /* Loop over the ligands */
  while (ligand_ptr != NULL)
    {
      /* Only process if not a metal */
      if (ligand_ptr->type != METAL)
        {
          /* Initialise description */
          description[0] = '\0';

          /* Get the first of this ligand's Het Groups */
          het_ptr = ligand_ptr->first_het_ptr;
          
          /* Loop over the Het Group to write out */
          while (het_ptr != NULL)
            {
              /* Format the Het Group name */
              strcpy(desc,het_ptr->desc);
              len = strlen(desc) - 1;
              to_lower(desc + 1,len);

              /* If we have only a single Het Group, save its name */
              if (ligand_ptr->nhets == 1)
                {
                  /* Save the description */
                  strcpy(description,"- ");
                  strcat(description,desc);
                }

              /* Otherwise, have a list of Het Groups in brackets */
              else
                {
                  /* If the first of this ligand's Het Group, then start
                     bracketed list */
                  if (description[0] == '\0')
                    {
                      strcpy(description,"(");
                      strcat(description,het_ptr->name);
                    }
                  else
                    {
                      strcat(description,", ");
                      strcat(description,het_ptr->name);
                    }

                  /* Add description of Het Group */
                  strcat(description," = ");
                  strcat(description,desc);
                }

              /* Get the next Het Group records */
              het_ptr = het_ptr->next_het_ptr;
            }

          /* If no description, see if this is a single amino acid */
          len = strlen(ligand_ptr->ligand_name);
          if (description[0] == '\0' && len == 3)
            check_for_aa(ligand_ptr->ligand_name,description);

          /* If we have a description, save it */
          if (description[0] != '\0')
            {
              /* If we had a bracketed list, then close bracket */
              if (ligand_ptr->nhets > 1)
                strcat(description,")");

              /* Save the description */
              store_string(&(ligand_ptr->description),description);
            }
        }

      /* Get next ligand record */
      ligand_ptr = ligand_ptr->next_ligand_ptr;
    }

  /* Get pointer to the first liglist record */
  liglist_ptr = chain_ptr->first_liglist_ptr;

  /* Loop to write out all the interacting ligands */
  while (liglist_ptr != NULL)
    {
      /* Get the corresponding rasdef record */
      r_rasdef_ptr = liglist_ptr->r_rasdef_ptr;

      /* If this is a new ligand, write out its details */
      if (strcmp(r_rasdef_ptr->ligand_name,last_ligname))
        {
          /* If have a previous ligand, write out its copy count */
          if (nligs > 0)
            fprintf(file_ptr,"LIG_NCOPY%s[%d] %d\n",sub1,nligs,ncopy);

          /* Increment ligand count */
          nligs++;

          /* Write out ligand name */
          fprintf(file_ptr,"LIG_NAME%s[%d] %s\n",sub1,nligs,
                  r_rasdef_ptr->ligand_name);

          /* Write out the type */
          fprintf(file_ptr,"LIG_TYPE%s[%d] %d\n",sub1,nligs,
                  r_rasdef_ptr->type);

          /* Determine whether this is a drug molecule */
          hetdrug_ptr
            = binary_hetdrug_search(hetdrug_index_ptr,nhetdrug,
                                    r_rasdef_ptr->ligand_name);

          /* If found, write out DrugPort id and generic
             name */
          if (hetdrug_ptr != NULL)
            {
              /* Write out drug flag */
              fprintf(file_ptr,"LIG_IS_DRUG%s[%d] TRUE\n",sub1,nligs);
              
              /* Write out drug id */
              fprintf(file_ptr,"LIG_DRUG_ID%s[%d] %s\n",sub1,nligs,
                      hetdrug_ptr->drug_id);
              
              /* Write out drug description */
              fprintf(file_ptr,"LIG_DRUG_DESC%s[%d] %s\n",sub1,nligs,
                      hetdrug_ptr->drug_desc);

              /* Write out generic name */
              fprintf(file_ptr,"LIG_GENERIC_NAME%s[%d] %s\n",
                      sub1,nligs,hetdrug_ptr->generic_name);
            }
          else
            fprintf(file_ptr,"LIG_IS_DRUG%s[%d] FALSE\n",sub1,nligs);

          /* Get the ligand name and convert to upper case */
          strcpy(ligand_name,r_rasdef_ptr->ligand_name);
          len = strlen(ligand_name);
          to_upper(ligand_name,len);

          /* Identify the corresponding ligand */
          ligand_ptr = find_ligand(chain_ptr,ligand_name);
          nhets = 0;

          /* If found, then write out name */
          if (ligand_ptr != NULL)
            {      
              /* Write out ligand description */
              if (ligand_ptr->description != &null)
                fprintf(file_ptr,"LIG_DESC%s[%d] %s\n",sub1,nligs,
                        ligand_ptr->description);
              else if (ligand_ptr->nhets == 0)
                fprintf(file_ptr,"LIG_DESC%s[%d] %s\n",sub1,nligs,
                        ligand_ptr->ligand_name);

              /* Get the first of this ligand's Het Groups */
              het_ptr = ligand_ptr->first_het_ptr;

              /* Loop over the Het Group to write out */
              while (het_ptr != NULL)
                {
                  /* Increment Het group count */
                  nhets++;

                  /* Write out this Het Group */
                  fprintf(file_ptr,"LIG_HET_NAME%s[%d][%d] %s\n",
                          sub1,nligs,nhets,het_ptr->name);
                  fprintf(file_ptr,"LIG_HET_DESC%s[%d][%d] %s\n",
                          sub1,nligs,nhets,het_ptr->desc);

                  /* Get the next Het Group records */
                  het_ptr = het_ptr->next_het_ptr;
                }

              /* Write out the number of Het Groups */
              fprintf(file_ptr,"LIG_NHETS%s[%d] %d\n",sub1,nligs,
                      nhets);
            }

          /* Initialise ligand count */
          ncopy = 0;

          /* Save the ligand name */
          strcpy(last_ligname,r_rasdef_ptr->ligand_name);
        }

      /* Increment copy count */
      ncopy++;

      /* Write out this ligand's detail */
      fprintf(file_ptr,"LIG_RES_NAME%s[%d][%d] %s\n",sub1,nligs,ncopy,
              r_rasdef_ptr->res_name);
      fprintf(file_ptr,"LIG_RES_NUM%s[%d][%d] %s\n",sub1,nligs,ncopy,
              r_rasdef_ptr->res_num);
      fprintf(file_ptr,"LIG_CHAIN%s[%d][%d] %c\n",sub1,nligs,ncopy,
              r_rasdef_ptr->chain);

      /* Write out the Ligplot id */
      fprintf(file_ptr,"LIG_LIGPLOT%s[%d][%d] %s\n",sub1,nligs,ncopy,
              r_rasdef_ptr->ligplot_id);

      /* If this is the first, then write out Ligplot id for ligand */
      if (ncopy == 1)
        {
          /* Show the LigPlot id */
          fprintf(file_ptr,"LIG_LIGPLOT_ID%s[%d] %s\n",sub1,nligs,
                  r_rasdef_ptr->ligplot_id);

          /* Extract the corresponding level numbers */
          level1 = level2 = 1;
          strncpy(string,r_rasdef_ptr->ligplot_id,2);
          string[2] = '\0';
          level1 = atoi(string);
          strncpy(string,r_rasdef_ptr->ligplot_id + 3,2);
          string[2] = '\0';
          level2 = atoi(string);

          /* Write out the level numbers */
          fprintf(file_ptr,"LIG_L_ID%s[%d] %d.%d\n",sub1,nligs,
                  level1,level2);
        }

      /* Get next definition record */
      liglist_ptr = liglist_ptr->next_liglist_ptr;
    }

  /* If have a last ligand, write out its count */
  if (nligs > 0)
    fprintf(file_ptr,"LIG_NCOPY%s[%d] %d\n",sub1,nligs,ncopy);

  /* Write out the number of different ligand types */
  fprintf(file_ptr,"LIG_NLIGANDS%s %d\n",sub1,nligs);
  if (nligs > 0)
    fprintf(file_ptr,"LIG_HAVE_LIGANDS%s TRUE\n",sub1);
  else
    fprintf(file_ptr,"LIG_HAVE_LIGANDS%s FALSE\n",sub1);

  /* Return number of ligands written out */
  return(nligs);
}
/***********************************************************************

check_other_chains  -  Check this protein's other chains and save what
                       they are

***********************************************************************/

void check_other_chains(FILE *file_ptr,char *sub1,struct chain *chain_ptr,
                        struct chain *first_chain_ptr)
{
  char human_list[MAX_LABEL_WIDTH + 1], other_list[MAX_LABEL_WIDTH + 1];
  char other_covid_list[MAX_LABEL_WIDTH + 1];
  char tmp[LINELEN];

  char *string_ptr;

  int subscript;
  int nhuman, nother, nother_covid;
  int human, other, other_covid;

  struct chain *other_chain_ptr;
  struct code *code_ptr, *other_code_ptr;

  /* Initialise variables */
  nother_covid = nhuman = nother = 0;
  other_covid_list[0] = human_list[0] = other_list[0] = '\0';

  /* Get the corresponding code record */
  code_ptr = chain_ptr->code_ptr[0];
  if (code_ptr == NULL)
    return;

  /* Get the first of the chain records */
  other_chain_ptr = first_chain_ptr;

  /* Loop over the other chains to see what they are */
  while (other_chain_ptr != NULL)
    {
      /* If this isn't our chain, see what it is */
      if (other_chain_ptr != chain_ptr &&
          other_chain_ptr->copy_of == '\0')
        {
          /* Get the corresponding code record */
          other_code_ptr = other_chain_ptr->code_ptr[0];

          /* Initialise flags */
          other_covid = human = other = FALSE;

          /* If don't have a UniProt code, classify as "other" */
          if (other_code_ptr == NULL)
            other = TRUE;

          /* If it's a human protein, set flag */
          else if (strstr(other_code_ptr->uniprot_id,"_HUMAN") != NULL)
            human = TRUE;

          /* Otherwise, see if it might be a different COVID-19 protein */
          else if (strcmp(other_code_ptr->uniprot_acc,
                          code_ptr->uniprot_acc))
            {
              /* Get the UniProt accession */
              strcpy(tmp,other_code_ptr->uniprot_acc);

              /* Remove subscript, if there is one */
              subscript = 0;
              string_ptr = strchr(tmp,'_');
              if (string_ptr != NULL)
                {
                  /* Extract the subscript */
                  subscript = atoi(string_ptr + 1);
                  string_ptr[0] = '\0';
                }

              /* If this is a COVID19 structure, set
                 flag */
              if (subscript > 0 ||
                  strstr(other_code_ptr->uniprot_id,"_SARS2") != NULL)
                other_covid = TRUE;
              else
                other = TRUE;
            }

          /* If we have any other proteins, save the name in the
             appropriate list */
          if (human == TRUE || other_covid == TRUE || other == TRUE)
            {
              /* Get the protein name */
              if (other_code_ptr != NULL &&
                  other_code_ptr->name != &null)
                strcpy(tmp,other_code_ptr->name);
              else
                sprintf(tmp,"Chain %c",other_chain_ptr->chain_id);

              /* Remove unnecessary text from the protein name */
              string_ptr = strstr(tmp," OS=");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              string_ptr = strstr(tmp," (");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
              string_ptr = strstr(tmp," - ");
              if (string_ptr != NULL)
                string_ptr[0] = '\0';

              /* Add to appropriate list */
              if (human == TRUE)
                nhuman = add_name(human_list,tmp,nhuman);
              else if (other_covid == TRUE)
                nother_covid
                  = add_name(other_covid_list,tmp,nother_covid);
              else
                nother = add_name(other_list,tmp,nother);
            }
        }

      /* Get pointer to the next chain record */
      other_chain_ptr
        = other_chain_ptr->next_chain_ptr;
    }

  /* Write out flags and names, where relevant */

  /* Other COVID-19 proteins */
  if (nother_covid > 0)
    {
      /* Write out flag  and count */
      fprintf(file_ptr,"UNIPLOT_OTHER_COVID%s TRUE\n",sub1);
      fprintf(file_ptr,"UNIPLOT_NOTHER_COVID%s %d\n",sub1,nother_covid);

      /* Write out list of proteins */
      fprintf(file_ptr,"UNIPLOT_NOTHER_COVID_PROT%s %s\n",sub1,
              other_covid_list);
    }

  /* Other human proteins */
  if (nother > 0)
    {
      /* Write out flag  and count */
      fprintf(file_ptr,"UNIPLOT_OTHER%s TRUE\n",sub1);
      fprintf(file_ptr,"UNIPLOT_NOTHER%s %d\n",sub1,nother);

      /* Write out list of proteins */
      fprintf(file_ptr,"UNIPLOT_OTHER_PROT%s %s\n",sub1,other_list);
    }

  /* Other proteins */
  if (nhuman > 0)
    {
      /* Write out flag  and count */
      fprintf(file_ptr,"UNIPLOT_HUMAN%s TRUE\n",sub1);
      fprintf(file_ptr,"UNIPLOT_NHUMAN%s %d\n",sub1,nhuman);

      /* Write out list of proteins */
      fprintf(file_ptr,"UNIPLOT_HUMAN_PROT%s %s\n",sub1,human_list);
    }
}
/***********************************************************************

write_uniplot_dat  -  Write all data to the uniplot.dat file

***********************************************************************/

void write_uniplot_dat(struct pdb *pdb_ptr,FILE *donefile_ptr,
                       char *pdbsum_dir,
                       struct hetdrug **hetdrug_index_ptr,int nhetdrug,
                       struct parameters params)
{
  char chain, file_name[FILENAME_LEN], resno[MAX_SEQLEN];
  char add_chain[3], chain_list[LINELEN], number_string[20];
  char tmp[LINELEN], sub1[5];

  char *string_ptr;

  int i, icode, ipos, ires, j, len, ncopies, plot_no, num_len, subscript;
  int fcode, nchains, ncodes, nunder;
  int have_cath, have_pfam[NDOMAIN_TYPES];

  struct chain *chain_ptr, *other_chain_ptr;
  struct code *code_ptr, *other_code_ptr;
  struct smatch *smatch_ptr;

  FILE *auditfile_ptr, *file_ptr;

  /* Initialise variables */
  have_cath = FALSE;
  for (i = 0; i < NDOMAIN_TYPES; i++)
    have_pfam[i] = FALSE;
  nchains = 0;
  sub1[0] = '\0';

  /* Write PDB code to audit file */
  auditfile_ptr = open_output_file("audit.out",TRUE);
  fprintf(auditfile_ptr,"Finishing PDB file %s\n",pdb_ptr->pdb_code);
  fclose(auditfile_ptr);

  /* Form the name of the output uniplot.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.covid19 == TRUE)
    strcat(file_name,"covplot.dat");
  else
    strcat(file_name,"uniplot.dat");
#ifdef TEST
  if (params.covid19 == TRUE)
    strcpy(file_name,"covplot.dat");
  else
    strcpy(file_name,"uniplot.dat");
#endif

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Write out PDB code */
  fprintf(file_ptr,"UNIPLOT_PDB_CODE %s\n",pdb_ptr->pdb_code);

  /* Write out flag indicating new version of SST diagrams */
  fprintf(file_ptr,"UNIPLOT_NEW TRUE\n");

  /* Get the first of this PDB entry's chains */
  chain_ptr = pdb_ptr->first_chain_ptr;
  plot_no = 0;

  /* Loop through all the chain records to write out info on
     Pfam and uniplot images */
  while (chain_ptr != NULL)
    {
      /* If this is a unique chain write out details */
      if (chain_ptr->copy_of == '\0')
        {
          /* Initialise first and last codes */
          nchains++;
          fcode = 0;
          ncodes = chain_ptr->ncodes;

          /* For covid-19 run, only record relevant UniProt code */
          if (params.covid19 == TRUE)
            {
              /* Only show one UniProt sequence */
              ncodes = 1;

              /* Write out the key record for this chain */
              fprintf(file_ptr,"KEY: %d\n",nchains);
            }

          /* Loop over the UniProt codes for this sequence */
          for (icode = fcode; icode < ncodes; icode++)
            {
              /* Get the corresponding UniProt code */
              code_ptr = chain_ptr->code_ptr[icode];
              smatch_ptr = chain_ptr->smatch_ptr[icode];

              /* If have CATH domain info, set flag */
              if (chain_ptr->first_cath_domain_ptr != NULL)
                have_cath = TRUE;

              /* Increment unique count */
              plot_no = chain_ptr->plot_no[icode];
              sprintf(sub1,"[%d]",plot_no);

              /* For covid-19 run, don't need subscript */
              if (params.covid19 == TRUE)
                sub1[0] = '\0';

              /* Get the chain id and write out */
              chain = chain_ptr->chain_id;
              fprintf(file_ptr,"UNIPLOT_CHAIN%s %c\n",sub1,
                      chain);

              /* Write out number of residues, if non-zero */
              if (chain_ptr->nresid > 0)
                fprintf(file_ptr,"UNIPLOT_NRESID%s %d\n",sub1,
                        chain_ptr->nresid);
              if (code_ptr->seq_len > 0)
                fprintf(file_ptr,"UNIPLOT_SEQLEN%s %d\n",sub1,
                        code_ptr->seq_len);

              /* Write UniProt code */
              if (code_ptr->uniprot_acc[0] != '\0')
                {
                  /* Get the UniProt accession */
                  strcpy(tmp,code_ptr->uniprot_acc);

                  /* Remove subscript, if there is one */
                  subscript = 0;
                  string_ptr = strchr(tmp,'_');
                  if (string_ptr != NULL)
                    {
                      /* Extract the subscript */
                      subscript = atoi(string_ptr + 1);
                      string_ptr[0] = '\0';
                    }

                  /* Write out the UniProt accession */
                  fprintf(file_ptr,"UNIPLOT_UNIPROT_ACC%s %s\n",
                          sub1,tmp);

                  /* If this is a COVID19 structure, write flag and
                     subscript */
                  if (subscript > 0 ||
                      strstr(code_ptr->uniprot_id,"_SARS2") != NULL)
                    {
                      /* Write out the COVID-19 flag and subscript */
                      fprintf(file_ptr,"UNIPLOT_COVID_19%s TRUE\n",
                              sub1);
                      fprintf(file_ptr,"UNIPLOT_SUBSCRIPT%s %d\n",
                              sub1,subscript);

                      /* Write out which covid-19 protein this is */
                      fprintf(file_ptr,"UNIPLOT_COVID_PROTEIN%s %d\n",
                              sub1,code_ptr->covid_protein);

                      /* Check the other chains to see if any are other
                         COVID-19 proteins or human proteins */
                      check_other_chains(file_ptr,sub1,chain_ptr,
                                         pdb_ptr->first_chain_ptr);
                    }
                }

              /* Write UniProt id */
              if (code_ptr->uniprot_id[0] != '\0')
                {
                  /* Get the UniProt id */
                  strcpy(tmp,code_ptr->uniprot_id);

                  /* Remove subscript, if there is one */
                  ipos = -1;
                  nunder = 0;
                  len = strlen(tmp);
                  for (i = 0; i < len; i++)
                    if (tmp[i] == '_')
                      {
                        ipos = i;
                        nunder++;
                      }
                  if (nunder > 1 && ipos > -1)
                    tmp[ipos] = '\0';

                  /* Write out the UniProt id */
                  fprintf(file_ptr,"UNIPLOT_UNIPROT_ID%s %s\n",
                          sub1,tmp);
                }

              /* Write out the protein name, if we have one */
              if (code_ptr->name != &null)
                {
                  /* Remove unnecessary text */
                  strcpy(tmp,code_ptr->name);
                  string_ptr = strstr(tmp," OX=");
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Replace "OS=" by "from " */
                  if (strstr(tmp,"OS=") != NULL)
                    perform_splice(tmp,"OS=","from ");

                  /* Write out the protein name */
                  fprintf(file_ptr,"UNIPLOT_UNIPROT_NAME%s %s\n",
                          sub1,tmp);
                }
                  
              /* Write out the list of interacting ligands */
              write_int_ligands(file_ptr,chain_ptr,sub1,
                                hetdrug_index_ptr,nhetdrug);

              /* If have match data, then write out */
              if (smatch_ptr != NULL)
                {
                  /* Write out the alignment length */
                  if (smatch_ptr->align_len > 0)
                    fprintf(file_ptr,"UNIPLOT_ALIGN_LEN%s %d\n",sub1,
                            smatch_ptr->align_len);

                  /* Write out the UniProt sequence */
                  if (smatch_ptr->uniprot_seq != &null)
                    {
                      /* Replace spaces by hyphens */
                      for (i = 0; i < strlen(smatch_ptr->uniprot_seq); i++)
                        if (smatch_ptr->uniprot_seq[i] == ' ')
                          smatch_ptr->uniprot_seq[i] = '-';

                      /* Write out */
                      fprintf(file_ptr,"UNIPLOT_UNIPROT_SEQ%s %s\n",
                              sub1,smatch_ptr->uniprot_seq);

                      /* Write out the hyphens for this sequence */
                      fprintf(file_ptr,"UNIPLOT_UNIPROT_HYPHENS%s ",
                              sub1);
                      ires = 0;
                      for (i = 0; i < strlen(smatch_ptr->uniprot_seq); i++)
                        if (smatch_ptr->uniprot_seq[i] != '-')
                          {
                            /* Increment residue counter and show + when 
                               divisible by 10 */
                            ires++;
                            if (ires % 10 == 0)
                              fprintf(file_ptr,"+");
                            else
                              fprintf(file_ptr,"-");
                          }
                        else
                          fprintf(file_ptr," ");
                      fprintf(file_ptr,"\n");
                      
                      /* Initialise variables */
                      ires = 0;
                      for (i = 0; i < smatch_ptr->align_len; i++)
                        resno[i] = '.';
                      resno[smatch_ptr->align_len] = '\0';

                      /* Loop through the sequence to find multiples
                         of 10 */
                      for (i = 0; i < strlen(smatch_ptr->uniprot_seq); i++)
                        if (smatch_ptr->uniprot_seq[i] != '-')
                          {
                            /* Increment residue counter and write 
                               out number when divisible by 10 */
                            ires++;
                            if (ires % 10 == 0)
                              {
                                /* Form the number */
                                sprintf(number_string,"%d",ires);
                                num_len = strlen(number_string);

                                /* Copy across */
                                if (i - num_len > -1)
                                  {
                                    for (j = 0; j < num_len; j++)
                                      {
                                        /* If have valid position,
                                           copy across */
                                        ipos = i + j - num_len + 1;
                                        if (ipos > -1)
                                          resno[ipos] = number_string[j];
                                      }
                                  }
                              }
                          }
                      
                      /* Write out the UniProt sequence numbers */
                      fprintf(file_ptr,"UNIPLOT_UNIPROT_RESNO%s %s\n",
                              sub1,resno);
                    }

                  /* Repeat for PDB sequence */
                  if (smatch_ptr->pdb_seq != &null)
                    {
                      /* Replace spaces by hyphens */
                      for (i = 0; i < strlen(smatch_ptr->pdb_seq); i++)
                        if (smatch_ptr->pdb_seq[i] == ' ')
                          smatch_ptr->pdb_seq[i] = '-';

                      /* Write out */
                      fprintf(file_ptr,"UNIPLOT_PDB_SEQ%s %s\n",sub1,
                              smatch_ptr->pdb_seq);
                    }
              
                  /* If have PDB residue numbers, write these out */
                  if (smatch_ptr->pdb_resno != &null)
                    {
                      /* Replace spaces by dots */
                      for (i = 0; i < strlen(smatch_ptr->pdb_resno); i++)
                        if (smatch_ptr->pdb_resno[i] == ' ')
                          smatch_ptr->pdb_resno[i] = '.';

                      /* Write out */
                      fprintf(file_ptr,"UNIPLOT_PDB_RESNO%s %s\n",sub1,
                              smatch_ptr->pdb_resno);
                    }

                  /* Repeat for hyphens */
                  if (smatch_ptr->pdb_hyphens != &null)
                    {
                      /* Replace spaces by dots */
                      for (i = 0; i < strlen(smatch_ptr->pdb_hyphens); i++)
                        if (smatch_ptr->pdb_hyphens[i] == ' ')
                          smatch_ptr->pdb_hyphens[i] = '.';

                      /* Write out */
                      fprintf(file_ptr,"UNIPLOT_PDB_HYPHENS%s %s\n",
                              sub1,smatch_ptr->pdb_hyphens);
                    }

                  /* If there are any differences, write out the
                     differences string */
                  if (smatch_ptr->diffs != &null &&
                      smatch_ptr->uniprot_seq != &null)
                    fprintf(file_ptr,"UNIPLOT_DIFFS%s %s\n",sub1,
                            smatch_ptr->diffs);

                  /* Write out number of mismatched residues, if
                     non-zero */
                  if (smatch_ptr->nmismatch > 0 &&
                      smatch_ptr->uniprot_seq != &null)
                    fprintf(file_ptr,"UNIPLOT_NMISMATCH%s %d\n",sub1,
                            smatch_ptr->nmismatch);
                }
              else
                fprintf(file_ptr,"UNIPLOT_NMISMATCH%s 0\n",sub1);

              /* Get next chain */
              other_chain_ptr = chain_ptr->next_chain_ptr;
              ncopies = 0;
              chain_list[0] = '\0';

              /* Loop to write out all equivalent chains */
              while (other_chain_ptr != NULL)
                {
                  /* If this is the right chain, write out */
                  if (other_chain_ptr->copy_of == chain)
                    {
                      /* Increment count of copies */
                      ncopies++;

                      fprintf(file_ptr,"UNIPLOT_COPY%s[%d] %c\n",
                              sub1,ncopies,
                              other_chain_ptr->chain_id);

                      /* Add to chain list */
                      sprintf(add_chain,"/%c",other_chain_ptr->chain_id);
                      strcat(chain_list,add_chain);
                    }

                  /* Get pointer to the next chain record */
                  other_chain_ptr
                    = other_chain_ptr->next_chain_ptr;
                }

              /* If have any copy chains, write out count */
              if (ncopies > 0)
                fprintf(file_ptr,"UNIPLOT_NCOPIES%s %d\n",
                        sub1,ncopies);

              /* If have code, write out the details */
              if (code_ptr->have_seq == TRUE)
                {
                  /* Write out flag */
                  fprintf(file_ptr,"UNIPLOT_HAVE_UNIPROT%s TRUE\n",
                          sub1);

                  /* Write out all of this chain's Pfam domains */
                  write_pfam_details(file_ptr,chain_ptr,smatch_ptr,
                                     sub1);

                  /* If have Pfam or CATH domain information, write out
                     the domain data */
                  if (smatch_ptr->first_domain_ptr != NULL ||
                      chain_ptr->first_cath_domain_ptr != NULL)
                    write_domain_details(file_ptr,chain_ptr,icode,
                                         smatch_ptr,sub1,have_pfam);
                }

              /* Write out number of image files */
              if (chain_ptr->nimages[icode] > 0)
                fprintf(file_ptr,"UNIPLOT_NIMAGES%s %d\n",
                        sub1,chain_ptr->nimages[icode]);
            }
        }

      /* Get pointer to the next chain record */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Write out the number of chain entries */
  fprintf(file_ptr,"NSPROT_ID %d\n",plot_no);

  /* Write out flags indicating which domains we have */
  if (have_pfam[PFAMA] == TRUE)
    fprintf(file_ptr,"UNIPLOT_HAVE_PFAMA TRUE\n");
  if (have_pfam[PFAMB] == TRUE)
    fprintf(file_ptr,"UNIPLOT_HAVE_PFAMB TRUE\n");
  if (have_pfam[COVID19] == TRUE)
    fprintf(file_ptr,"UNIPLOT_ENTRY_TYPES19 TRUE\n");

  /* If any CATH domains, write out flag indicating we have those */
  if (have_cath == TRUE)
    fprintf(file_ptr,"UNIPLOT_HAVE_CATH TRUE\n");

  /* Close the output file */
  fclose(file_ptr);

  /* Write PDB code to done file */
  fprintf(donefile_ptr,"%s\n",pdb_ptr->pdb_code);
}
/***********************************************************************

free_domain  -  Free up the memory taken up by the current entry's data
                details

***********************************************************************/

void free_domain(struct domain **fst_domain_ptr)
{
  struct domain *domain_ptr, *next_domain_ptr;

  /* Get first domain record */
  domain_ptr = *fst_domain_ptr;

  /* Loop through the domains until done */
  while (domain_ptr != NULL)
    {
      /* Get pointer to the next domain record */
      next_domain_ptr = domain_ptr->next_domain_ptr;

      /* Free up memory taken by current domain */
      free(domain_ptr);

      /* Go to the next domain */
      domain_ptr = next_domain_ptr;
    }

  /* Re-initialise first record */
  *fst_domain_ptr = NULL;
}
/***********************************************************************

free_chain  -  Free up the memory taken up by the current entry's data
               details

***********************************************************************/

void free_chain(struct chain **fst_chain_ptr)
{
  int icode;

  struct chain *chain_ptr, *next_chain_ptr;
  struct smatch *smatch_ptr;

  /* Get first chain record */
  chain_ptr = *fst_chain_ptr;

  /* Loop through the chains until done */
  while (chain_ptr != NULL)
    {
      /* Get pointer to the next chain record */
      next_chain_ptr = chain_ptr->next_chain_ptr;

      /* Free up sequence match data */
      for (icode = 0; icode < chain_ptr->ncodes; icode++)
        {
          /* Get the corresponding smatch record */
          smatch_ptr = chain_ptr->smatch_ptr[icode];

          /* If not null, free up memory */
          if (smatch_ptr != NULL)
            {
              /* Free up the sequences */
              if (smatch_ptr->uniprot_seq != &null)
                free(smatch_ptr->uniprot_seq);
              if (smatch_ptr->pdb_seq != &null)
                free(smatch_ptr->pdb_seq);
              if (smatch_ptr->pdb_hyphens != &null)
                free(smatch_ptr->pdb_hyphens);
              if (smatch_ptr->pdb_resno != &null)
                free(smatch_ptr->pdb_resno);
              if (smatch_ptr->diffs != &null)
                free(smatch_ptr->diffs);

              /* Free up the domain records */
              free_domain(&(smatch_ptr->first_domain_ptr));

              /* Clear the smatch record */
              free(smatch_ptr);
            }
        }

      /* Free up the memory taken by the Pfam and CATH domains */
      free_domain(&(chain_ptr->first_cath_domain_ptr));

      /* Free up memory taken by current chain */
      //      free(chain_ptr);

      /* Go to the next chain */
      chain_ptr = next_chain_ptr;
    }

  /* Re-initialise first record */
  *fst_chain_ptr = NULL;
}
/***********************************************************************

finish_entries  -  Write uniplot.dat files for all completed PDB
                   entries and free up memory taken by their data

***********************************************************************/

void finish_entries(struct code *code_ptr,FILE *donefile_ptr,
                    struct hetdrug **hetdrug_index_ptr,int nhetdrug,
                    struct c_file_data file_name_ptr,
                    struct parameters params)
{
  char pdbsum_dir[FILENAME_LEN + 1];

  int dir_type, i, icode;
  int done, multimap;

  struct chain *chain_ptr;
  struct clist *clist_ptr;
  struct pdb *pdb_ptr;

  /* Get pointer to the first clist entry for this UniProt code */
  clist_ptr = code_ptr->first_clist_ptr;

  /* Loop through all entries to process any completed PDB entries */
  while (clist_ptr != NULL)
    {

      /* Get the corresponding chain and PDB entry */
      chain_ptr = clist_ptr->chain_ptr;
      icode = clist_ptr->icode;
      pdb_ptr = chain_ptr->pdb_ptr;

      /* If this PDB entry has not yet been done, then process it */
      if (pdb_ptr->done == FALSE)
        {
          /* For PDBsum1, redefine current directory as PDBsum
             directory */
          if (params.pdb_type == PDBSUM1)
            {
              strcpy(pdbsum_dir,"./");
              dir_type = PDBSUM1;
            }

          /* Form the name of the PDBsum directory */
          else
            dir_type = find_pdbsum_dir(pdb_ptr->pdb_code,
                                       file_name_ptr.pdbsum_template,
                                       params.pdbsum_type,pdbsum_dir);

          /* Assume entry is complete and each chain maps onto a single
             UniProt code */
          done = TRUE;
          multimap = FALSE;

          /* Get this PDB entry's first chain */
          chain_ptr = pdb_ptr->first_chain_ptr;

          /* Check whether this PDB entry is complete */
          while (chain_ptr != NULL && done == TRUE)
            {
              /* If this is a unique chain check whether corresponding
                 UniProt sequence has been retrieved */
              if (chain_ptr->copy_of == '\0')
                {
                  /* Loop over the chain's UniProt codes */
                  for (i = 0; i < chain_ptr->ncodes; i++)
                    {
                      /* If chain has not been processed, then can't delete
                         the current PDB entry */
                      if (chain_ptr->code_ptr[i] != NULL &&
                          (chain_ptr->smatch_ptr[i] == NULL ||
                           chain_ptr->smatch_ptr[i]->align_len == 0))
                        done = FALSE;
                    }
                }

              /* Get pointer to the next chain record */
              chain_ptr = chain_ptr->next_chain_ptr;
            }

          /* If this PDB entry is done, then write out and free up its
             memory */
          if (done == TRUE && pdb_ptr->first_chain_ptr != NULL)
            {
              /* Write out the uniplot.dat file */
              write_uniplot_dat(pdb_ptr,donefile_ptr,pdbsum_dir,
                                hetdrug_index_ptr,nhetdrug,params);

              /* Free up memory taken by this entry's chains */
              free_chain(&(pdb_ptr->first_chain_ptr));

              /* Set flag */
              pdb_ptr->done = TRUE;
            }
        }

      /* Get pointer to the next entry in the linked list */
      clist_ptr = clist_ptr->next_clist_ptr;
    }

  /* Remove the clist pointers from the current UniProt code record */
  free_clist(&(code_ptr->first_clist_ptr));
}
/***********************************************************************

process_uniprot_seqs  -  Pick up the UniProt sequences from the UniProt
                         .fasta files

***********************************************************************/

int process_uniprot_seqs(struct code **code_index_ptr,int ncodes,
                         char *character,struct letter *letter_ptr,
                         int nletters,struct pfamlig **fst_pfamlig_ptr,
                         int *npflig,FILE *donefile_ptr,
                         struct code *first_code_ptr,int entry_types,
                         struct hetdrug **hetdrug_index_ptr,
                         int nhetdrug,struct c_file_data file_name_ptr,
                         struct parameters params)
{
#define NUNIPROT  3

  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char tmp_line[LINELEN + 1], name[LINELEN + 1];
  char fasta_id[7], uniprot_fasta[FILENAME_LEN];
  char uniprot_acc[UNIPROT_LEN], uniprot_id[UNIPROT_LEN];

  char *string_ptr;

  static char *FASTA_FILE[NUNIPROT] = {
    "COVID-19.fasta", "uniprot_sprot.fasta", "uniprot_trembl.fasta" };

  int err_no, floop, last_loop, ipos, jpos, len, line, loop;
  int npfamlig, nseqs, seq_len;
  int ndone, run_no, tlen;
  int dummy;

  struct code *code_ptr;
  struct pfamlig *first_pfamlig_ptr, *last_pfamlig_ptr;

  FILE *auditfile_ptr, *drfile_ptr, *outfile_ptr;
#ifdef TEST
  FILE *file_ptr;
#else
  FILE *file_ptr;
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  dummy = FALSE;
  fasta_id[0] = '\0';
  floop = 0;
  last_loop = NUNIPROT - 1;
  line = 0;
  ndone = 0;
  npfamlig = 0;
  nseqs = 0;
  outfile_ptr = NULL;
  run_no = 0;
  seq_len = 0;
  strcpy(uniprot_fasta,"uniprot.fasta");
  if (params.covid19 == TRUE || entry_types == COVID_ONLY)
    {
      floop = 0;
      last_loop = 1;
    }
  if (entry_types == NO_COVID)
    floop = 1;

  /* Initialise pointers */
  first_pfamlig_ptr = NULL;
  last_pfamlig_ptr = NULL;

  /* Open the pdbsum.dr file */
  drfile_ptr = open_output_file("pdbsum.dr",TRUE);

  /* Loop over the SWISS-PROT/TrEMBL data files */
  for (loop = floop; loop < last_loop && ndone < ncodes; loop++)
    {
      /* Initialise variables */
      code_ptr = NULL;

      /* On first loop, process the unassigned codes */
      if (loop == floop && first_code_ptr != NULL &&
          !strcmp(first_code_ptr->uniprot_acc,"NONE"))
        {
          code_ptr = first_code_ptr;
          dummy = TRUE;
        }

      /* Form name of input file */
      if (loop > 0)
        strcpy(file_name,file_name_ptr.other_dir[SWISSPROT_DIR]);
      else
        strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
      strcat(file_name,DIR_SEP);
      strcat(file_name,FASTA_FILE[loop]);
      if (loop > 0 && loop != 3)
        strcat(file_name,".gz");
      printf("Reading uniprot file %s ...\n",file_name);

      /* Open the file, if it exists */
      gzfile_ptr = gzopen(file_name,"rb");
      if (gzfile_ptr ==  NULL)
        {
          /* Check for insufficiency of memory */
          if(err_no)
            {
              printf("*** ERROR. Insufficient memory opening file [%s]\n",
                     file_name);
              return(nseqs);
            }

          /* Otherwise, error must be due to missing file */
          printf("\n*** Warning. Unable to open file [%s]\n",file_name);
          return(nseqs);
        }

      /* Loop while reading in records from the file */
      while (gzgets(gzfile_ptr,input_line,LINELEN)!= NULL && ndone < ncodes)
        {
          /* Increment line count */
          line++;

          /* If this is a new UniProt record, get code and check
             whether it is one of ours */
          if (input_line[0] == '>')
            {
              /* If have a previous sequence to process, do so now */
              if (code_ptr != NULL && code_ptr->done == FALSE)
                {
                  /* Close the output .fasta file */
                  if (dummy == FALSE)
                    fclose(outfile_ptr);

                  /* Find all PDB chains with this as their UniProt
                     sequence and run their sequences against it
                     using FASTA to determine overlap */
                  run_no = compare_seqs(code_ptr,uniprot_fasta,fasta_id,
                                        seq_len,run_no,character,
                                        letter_ptr,nletters,drfile_ptr,
                                        &first_pfamlig_ptr,
                                        &last_pfamlig_ptr,&npfamlig,
                                        file_name_ptr,params);

                  /* Write uniplot.dat files for all completed PDB
                     entries and free up memory taken by their data */
                  finish_entries(code_ptr,donefile_ptr,
                                 hetdrug_index_ptr,nhetdrug,
                                 file_name_ptr,params);

                  /* Flag code as done */
                  code_ptr->done = TRUE;

                  /* Increment count of codes done */
                  ndone++;
                }

              /* Truncate the input line */
              len = string_truncate(input_line,LINELEN);

              /* Increment count of sequence read in */
              nseqs++;

              /* Write to audit file, if required */
              if ((nseqs % NAUDIT) == 0)
                {
                  /* Write sequence count to audit file to allow
                     check on progress */
                  auditfile_ptr = open_output_file("audit.out",TRUE);
                  fprintf(auditfile_ptr,"Loop %d. Seqs %d\n",loop,nseqs);
                  fclose(auditfile_ptr);
                }

              /* Initialise UniProt id */
              uniprot_acc[0] = '\0';
              uniprot_id[0] = '\0';
              code_ptr = NULL;

              /* Initialise protein name */
              name[0] = '\0';

              /* Pick up the UniProt codes and protein name from the
                 input line */
              get_uniprot_data(input_line,uniprot_acc,uniprot_id,name);

              /* If have managed to extract the UniProt id, see if
                 it is one of the ones we're interested in */
              if (uniprot_acc[0] != '\0')
                {
                  /* See if this is one of our UniProt codes */
                  code_ptr
                    = binary_search(code_index_ptr,ncodes,uniprot_acc);

                  /* If valid code, create a .fasta file to hold
                     its UniProt sequence */
                  if (code_ptr != NULL)
                    {
                      /* Open a new file to store the sequence */
                      outfile_ptr = open_output_file(uniprot_fasta,TRUE);

                      /* Write current header line to it */
                      fprintf(outfile_ptr,"%s\n",input_line);

                      /* Extract sequence identifier as it will be
                         seen in the eventual fasta.out file */
                      strncpy(fasta_id,input_line+1,6);
                      fasta_id[6] = '\0';

                      /* Flag code as having a FASTA sequence */
                      code_ptr->have_seq = TRUE;

                      /* Store the protein name, if have it */
                      if (name[0] != '\0')
                        store_string(&(code_ptr->name),name);

                      /* Initialise sequence length */
                      seq_len = 0;

                      /* If have a UniProt id, store it */
                      if (uniprot_id[0] != '\0' &&
                          code_ptr->uniprot_id[0] == '\0')
                        strcpy(code_ptr->uniprot_id,uniprot_id);

                      /* Set flag */
                      dummy = FALSE;
                    }
                }
            }

          /* If sequence line belongs to a UniProt entry of interest,
             write out to .fasta file */
          else if (code_ptr != NULL && dummy == FALSE)
            {
              /* Truncate the input line */
              len = string_truncate(input_line,LINELEN);

              /* Increment UniProt sequence length */
              seq_len = seq_len + len;

              /* Write sequence line to .fasta file */
              fprintf(outfile_ptr,"%s\n",input_line);
            }
        }

      /* Close file */
#ifdef TEST
      fclose(file_ptr);
#else
      gzclose(gzfile_ptr);
#endif

      /* If have a previous sequence to process, do so now */
      if (code_ptr != NULL && code_ptr->done == FALSE)
        {
          /* Close the output .fasta file */
          if (outfile_ptr != NULL)
            {
              fclose(outfile_ptr);
              outfile_ptr = NULL;
            }

          /* Find all PDB chains with this as their UniProt
             sequence and run their sequences against it
             using FASTA to determine overlap */
          run_no
            = compare_seqs(code_ptr,uniprot_fasta,fasta_id,seq_len,
                           run_no,character,letter_ptr,nletters,
                           drfile_ptr,&first_pfamlig_ptr,
                           &last_pfamlig_ptr,&npfamlig,file_name_ptr,
                           params);

          /* Write uniplot.dat files for all completed PDB
             entries and free up memory taken by their data */
          finish_entries(code_ptr,donefile_ptr,hetdrug_index_ptr,nhetdrug,
                         file_name_ptr,params);

          /* Flag code as done */
          code_ptr->done = TRUE;

          /* Increment count of codes done */
          ndone++;
        }

      /* Return current pointers */
      *fst_pfamlig_ptr = first_pfamlig_ptr;

      /* Return count of ligand to Pfam domain links */
      *npflig = npfamlig;
    }

  /* Close the pdbsum.dr file */
  fclose(drfile_ptr);

  /* Show stats */
  printf("Number of UniProt seqs read:       %8d\n",nseqs);
  fflush(stdout);

  /* Return number of sequences */
  return(nseqs);
}
/***********************************************************************

process_alignments  -  Process the alginments of the PDB sequences vs
                       the equivalent UniProt sequences

***********************************************************************/

int process_alignments(struct code **code_index_ptr,int ncodes,
                       char *character,struct letter *letter_ptr,
                       int nletters,struct pfamlig **fst_pfamlig_ptr,
                       int *npflig,FILE *donefile_ptr,
                       struct code *first_code_ptr,int entry_types,
                       struct hetdrug **hetdrug_index_ptr,
                       int nhetdrug,struct c_file_data file_name_ptr,
                       struct parameters params)
{
  char fasta_id[7], uniprot_fasta[FILENAME_LEN];

  char *string_ptr;

  int icode, ndone, npfamlig, run_no, seq_len;
  int dummy;

  struct code *code_ptr;
  struct pfamlig *first_pfamlig_ptr, *last_pfamlig_ptr;

  FILE *drfile_ptr;

  /* Initialise variables */
  dummy = FALSE;
  ndone = 0;
  npfamlig = 0;
  run_no = 0;
  seq_len = 0;
  fasta_id[0] = '\0';
  uniprot_fasta[0] = '\0';
  printf("Processing alignments ...\n");
  fflush(stdout);

  /* Initialise pointers */
  first_pfamlig_ptr = last_pfamlig_ptr = NULL;

  /* Open the pdbsum.dr file */
  drfile_ptr = open_output_file("pdbsum.dr",TRUE);

  /* Loop over all the codes */
  for (icode = 0; icode < ncodes; icode++)
    {
      /* Get the corresponding code */
      code_ptr = code_index_ptr[icode];

      /* Perform the compare routines */
      run_no
        = compare_seqs(code_ptr,uniprot_fasta,fasta_id,
                       code_ptr->seq_len,run_no,character,letter_ptr,
                       nletters,drfile_ptr,&first_pfamlig_ptr,
                       &last_pfamlig_ptr,&npfamlig,file_name_ptr,
                       params);

      /* Write uniplot.dat files for all completed PDB
         entries and free up memory taken by their data */
      finish_entries(code_ptr,donefile_ptr,hetdrug_index_ptr,nhetdrug,
                     file_name_ptr,params);

      /* Increment count of codes done */
      ndone++;
    }

  /* Return current pointers */
  *fst_pfamlig_ptr = first_pfamlig_ptr;

  /* Return count of ligand to Pfam domain links */
  *npflig = npfamlig;
  
  /* Close the pdbsum.dr file */
  fclose(drfile_ptr);

  /* Return number of sequences processed */
  return(ndone);
}
/***********************************************************************

delete_fasta_files  -  Delete all the downloaded .fasta files

***********************************************************************/

void delete_fasta_files(struct code **code_index_ptr,int ncodes,int dos)
{
  char command[FILENAME_LEN], uniprot_fasta[FILENAME_LEN];

  int icode;
  int delete;

  struct code *code_ptr;

  /* Loop over all the codes */
  for (icode = 0; icode < ncodes; icode++)
    {
      /* Get the corresponding code */
      code_ptr = code_index_ptr[icode];

      /* If accession not NONE, then delete */
      if (strcmp(code_ptr->uniprot_acc,"NONE"))
        {
          /* Check for covid-19 files so these don't get deleted */
          delete = TRUE;
          if (!strncmp(code_ptr->uniprot_acc,"P0DTD1",6) ||
              !strncmp(code_ptr->uniprot_acc,"P0DTC",5))
            delete = FALSE;

          /* If to be deleted, then do it */
          if (delete == TRUE)
            {
              /* Form the name of the file */
              sprintf(uniprot_fasta,"%s.fasta",code_ptr->uniprot_acc);
          
              /* Delete it */
              sprintf(command,"%s %s",RM_COMMAND,uniprot_fasta);
#ifdef NOSYST
#else
              call_system(command,dos);
#endif
            }
        }
    }
}
/***********************************************************************

close_pdb_entries  -  Pick up the UniProt sequences from the UniProt
                     .fasta files

***********************************************************************/

void close_pdb_entries(struct pdb *first_pdb_ptr,FILE *donefile_ptr,
                       struct hetdrug **hetdrug_index_ptr,int nhetdrug,
                       struct c_file_data file_name_ptr,
                       struct parameters params)
{
  char pdbsum_dir[FILENAME_LEN + 1];

  int dir_type;

  struct pdb *pdb_ptr;

  /* Get pointer to the first PDB record */
  pdb_ptr = first_pdb_ptr;

  /* Loop through all entries to write out any incomplete entries */
  while (pdb_ptr != NULL)
    {
      /* If have any chains, then process */
      if (pdb_ptr->first_chain_ptr != NULL)
        {
          /* For PDBsum1, redefine current directory as PDBsum
             directory */
          if (params.pdb_type == PDBSUM1)
            {
              strcpy(pdbsum_dir,"./");
              dir_type = PDBSUM1;
            }

          /* Form the name of the PDBsum directory */
          else
            dir_type = find_pdbsum_dir(pdb_ptr->pdb_code,
                                       file_name_ptr.pdbsum_template,
                                       params.pdbsum_type,pdbsum_dir);

          /* Write out the uniplot.dat file */
          if (dir_type != NOT_FOUND)
            write_uniplot_dat(pdb_ptr,donefile_ptr,pdbsum_dir,
                              hetdrug_index_ptr,nhetdrug,params);
        }

      /* Get pointer to the next entry in the linked list */
      pdb_ptr = pdb_ptr->next_pdb_ptr;
    }
}
/***********************************************************************

pfamlig_sort1  -  Shell-sort the pfamlig entries by ligand and Pfam
                  domain

***********************************************************************/

void pfamlig_sort1(int *index,struct pfamlig **pfamlig_index_ptr,
                   int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  int swap;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct code *code1_ptr, *code2_ptr;
  struct pfamlig *pfamlig1_ptr, *pfamlig2_ptr;

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

              /* Get the two pfamlig records */
              pfamlig1_ptr = pfamlig_index_ptr[index[indi]];
              pfamlig2_ptr = pfamlig_index_ptr[index[indl]];

              /* Get the corresponding UniProt records */
              code1_ptr = pfamlig1_ptr->code_ptr;
              code2_ptr = pfamlig2_ptr->code_ptr;

              /* Determine whether pointers need to be swapped */
              swap = FALSE;
              if (strcmp(pfamlig1_ptr->ligand_name,
                         pfamlig2_ptr->ligand_name) > 0)
                swap = TRUE;
              else if (!strcmp(pfamlig1_ptr->ligand_name,
                               pfamlig2_ptr->ligand_name))
                {
                  /* Same name, so compare Pfam id */
                  if (strcmp(pfamlig1_ptr->domain_list,
                             pfamlig2_ptr->domain_list) > 0)
                    swap = TRUE;
                  else if (!strcmp(pfamlig1_ptr->domain_list,
                                   pfamlig2_ptr->domain_list))
                    {
                      /* Same Pfam id, so sort on UniProt code */
                      if (strcmp(code1_ptr->uniprot_acc,
                                 code2_ptr->uniprot_acc) > 0)
                        swap = TRUE;
                      else if (!strcmp(code1_ptr->uniprot_acc,
                                       code2_ptr->uniprot_acc))
                        {
                          /* Same UniProt code, so sort on PDB code */
                          if (strcmp(pfamlig1_ptr->pdb_code,
                                     pfamlig2_ptr->pdb_code) > 0)
                            swap = TRUE;
                          else if (!strcmp(pfamlig1_ptr->pdb_code,
                                           pfamlig2_ptr->pdb_code))
                            {
                              /* Same PDB code, so sort on chain-id */
                              if (pfamlig1_ptr->chain >
                                  pfamlig2_ptr->chain)
                                swap = TRUE;
                              else if (pfamlig1_ptr->chain ==
                                       pfamlig2_ptr->chain)
                                {
                                  /* Same chain, so sort on domain
                                     number */
                                  if (pfamlig1_ptr->domain_no >
                                      pfamlig2_ptr->domain_no)
                                    swap = TRUE;
                                  else if (pfamlig1_ptr->domain_no ==
                                           pfamlig2_ptr->domain_no)
                                    {
                                      /* Same domain number, so sort on
                                         ligand id */
                                      if (strcmp(pfamlig1_ptr->ligand_id,
                                                 pfamlig2_ptr->ligand_id) > 0)
                                        swap = TRUE;
                                    }
                                }
                            }
                        }
                    }
                }

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

pfamlig_sort2  -  Shell-sort the pfamlig entries by Pfam domain and
                  ligand

***********************************************************************/

void pfamlig_sort2(int *index,struct pfamlig **pfamlig_index_ptr,
                   int nindex)
{
  int endloop, i, indi, indl, iswap, j, k, l, lognb2, m, n, nn;
  int swap;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct code *code1_ptr, *code2_ptr;
  struct pfamlig *pfamlig1_ptr, *pfamlig2_ptr;

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

              /* Get the two pfamlig records */
              pfamlig1_ptr = pfamlig_index_ptr[index[indi]];
              pfamlig2_ptr = pfamlig_index_ptr[index[indl]];

              /* Get the corresponding UniProt records */
              code1_ptr = pfamlig1_ptr->code_ptr;
              code2_ptr = pfamlig2_ptr->code_ptr;

              /* Determine whether pointers need to be swapped */
              swap = FALSE;
              if (strcmp(pfamlig1_ptr->domain_list,
                         pfamlig2_ptr->domain_list) > 0)
                swap = TRUE;
              else if (!strcmp(pfamlig1_ptr->domain_list,
                               pfamlig2_ptr->domain_list))
                {
                  /* Same Pfam id, so compare ligand names */
                  if (strcmp(pfamlig1_ptr->ligand_name,
                             pfamlig2_ptr->ligand_name) > 0)
                    swap = TRUE;
                  else if (!strcmp(pfamlig1_ptr->ligand_name,
                                   pfamlig2_ptr->ligand_name))
                    {
                      /* Same ligand name, so sort on UniProt code */
                      if (strcmp(code1_ptr->uniprot_acc,
                                 code2_ptr->uniprot_acc) > 0)
                        swap = TRUE;
                      else if (!strcmp(code1_ptr->uniprot_acc,
                                       code2_ptr->uniprot_acc))
                        {
                          /* Same UniProt code, so sort on PDB code */
                          if (strcmp(pfamlig1_ptr->pdb_code,
                                     pfamlig2_ptr->pdb_code) > 0)
                            swap = TRUE;
                          else if (!strcmp(pfamlig1_ptr->pdb_code,
                                           pfamlig2_ptr->pdb_code))
                            {
                              /* Same PDB code, so sort on chain-id */
                              if (pfamlig1_ptr->chain >
                                  pfamlig2_ptr->chain)
                                swap = TRUE;
                              else if (pfamlig1_ptr->chain ==
                                       pfamlig2_ptr->chain)
                                {
                                  /* Same chain, so sort on domain
                                     number */
                                  if (pfamlig1_ptr->domain_no >
                                      pfamlig2_ptr->domain_no)
                                    swap = TRUE;
                                  else if (pfamlig1_ptr->domain_no ==
                                           pfamlig2_ptr->domain_no)
                                    {
                                      /* Same domain number, so sort on
                                         ligand id */
                                      if (strcmp(pfamlig1_ptr->ligand_id,
                                                 pfamlig2_ptr->ligand_id) > 0)
                                        swap = TRUE;
                                    }
                                }
                            }
                        }
                    }
                }

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

sort_pfamlig_entries  -  Sort the pfamligs by pfamlig name and containing Het
                       Group

***********************************************************************/

int sort_pfamlig_entries(struct pfamlig **fst_pfamlig_ptr,int npfamlig,
                         int order)
{
  int iorder, ipos, nentry;

  int *index;

  struct pfamlig *first_pfamlig_ptr, *last_pfamlig_ptr;
  struct pfamlig *pfamlig_ptr;

  struct pfamlig **pfamlig_index_ptr;

  /* Initialise pointers */
  first_pfamlig_ptr = *fst_pfamlig_ptr;
  last_pfamlig_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (int *) malloc(npfamlig * sizeof(int *));
  pfamlig_index_ptr
    = (struct pfamlig **) malloc(npfamlig * sizeof(struct pfamlig *));

  /* Get pointer to first pfamlig record */
  pfamlig_ptr = first_pfamlig_ptr;
  nentry = 0;

  /* Loop through pfamlig records to initialise the sort index */
  while (pfamlig_ptr != NULL)
    {
      /* If not off end off array, add to array */
      if (nentry < npfamlig)
        {
          /* Initialise sort index */
          index[nentry] = nentry;
          pfamlig_index_ptr[nentry] = pfamlig_ptr;

          /* Increment counter */
          nentry++;
        }

      /* Get next pfamlig record */
      pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;
    }

  /* Sort the pfamlig records in required order */
  if (order == LIG_PFAM)
    pfamlig_sort1(index,pfamlig_index_ptr,npfamlig);
  else
    pfamlig_sort2(index,pfamlig_index_ptr,npfamlig);

  /* Having sorted the records, rearrange in their sorted order */
  for (iorder = 0; iorder < npfamlig; iorder++)
    {
      /* Get the array position of this pfamlig record */
      ipos = index[iorder];

      /* Get the current pfamlig record */
      pfamlig_ptr = pfamlig_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_pfamlig_ptr = pfamlig_ptr;

      /* Otherwise, get the previous pfamlig record to point to this one */
      else
        last_pfamlig_ptr->next_pfamlig_ptr = pfamlig_ptr;

      /* Blank out current record's next pfamlig pointer */
      pfamlig_ptr->next_pfamlig_ptr = NULL;

      /* Store current record as last encountered */
      last_pfamlig_ptr = pfamlig_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);
  free(pfamlig_index_ptr);

  /* Return current pointer */
  *fst_pfamlig_ptr = first_pfamlig_ptr;

  /* Return number of sorted pfamligs */
  return(npfamlig);
}
/***********************************************************************

pdb_sort  -  Shell-sort which sorts the UniProt pdbs into ascending
              order

***********************************************************************/

void pdb_sort(struct pdb **pdb_index_ptr,int nindex)
{
  char pdb_name1[5], pdb_name2[5];

  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct pdb *pdb1_ptr, *pdb2_ptr, *swap_pdb_ptr;

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

              /* Get the two pdb records */
              pdb1_ptr = pdb_index_ptr[indi];
              pdb2_ptr = pdb_index_ptr[indl];

              /* Form full pdb names on which sort is based */
              strcpy(pdb_name1,pdb1_ptr->pdb_code);
              strcpy(pdb_name2,pdb2_ptr->pdb_code);

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(pdb_name1,pdb_name2) > 0)
                {
                  swap_pdb_ptr = pdb_index_ptr[indi];
                  pdb_index_ptr[indi] = pdb_index_ptr[indl];
                  pdb_index_ptr[indl] = swap_pdb_ptr;

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

create_pdb_index  -  Create index of PDB records for use in binary
                     search

***********************************************************************/

void create_pdb_index(struct pdb *first_pdb_ptr,
                      struct pdb ***pdb_ind_ptr,int npdb)
{
  int nentry;

  struct pdb *pdb_ptr;

  struct pdb **pdb_index_ptr;

  /* Initialise variables */

  /* Create the index array as one large array */
  pdb_index_ptr
    = (struct pdb **) malloc(npdb * sizeof(struct pdb *));

  /* Check that sufficient memory was allocated */
  if (pdb_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Get pointer to the first record */
  pdb_ptr = first_pdb_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (pdb_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      if (nentry < npdb)
        pdb_index_ptr[nentry] = pdb_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      pdb_ptr = pdb_ptr->next_pdb_ptr;
    }

  /* Sort the pdbs */
  if (npdb > 1)
    pdb_sort(pdb_index_ptr,npdb);

  /* Return the array pointer */
  *pdb_ind_ptr = pdb_index_ptr;
}
/***********************************************************************

binary_pdb_search  -  Search for the given PDB code

***********************************************************************/

struct pdb *binary_pdb_search(struct pdb **pdb_index_ptr,int npdb,
                              char *pdb_code)
{
  int bottom_point, ipoint, new_point, top_point;
  int found;

  struct pdb *pdb_ptr, *found_pdb_ptr, *prev_pdb_ptr;

  /* Initialise variables */
  found = FALSE;
  found_pdb_ptr = NULL;

  /* Make starting point half-way into the index */
  ipoint = npdb / 2;
  bottom_point = 0;
  top_point = npdb - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated pointer */
      pdb_ptr = pdb_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(pdb_code,pdb_ptr->pdb_code))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_pdb_ptr = pdb_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(pdb_code,pdb_ptr->pdb_code) < 0)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_pdb_ptr = pdb_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_pdb_ptr->pdb_code,pdb_code) < 0)
                found = TRUE;

              /* Otherwise, need to search earlier in the list */
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

      /* Otherwise, have to look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= npdb - 1)
            {
              /* Set flag as done */
              found = TRUE;
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

  /* Return the pointer */
  return(found_pdb_ptr);
}
/***********************************************************************

pfam_sort  -  Shell-sort which sorts the UniProt pfams into ascending
              order

***********************************************************************/

void pfam_sort(struct pfam **pfam_index_ptr,int nindex)
{
  char pfam_id1[PFAM_LEN], pfam_id2[PFAM_LEN];

  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct pfam *pfam1_ptr, *pfam2_ptr, *swap_pfam_ptr;

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

              /* Get the two pfam records */
              pfam1_ptr = pfam_index_ptr[indi];
              pfam2_ptr = pfam_index_ptr[indl];

              /* Form full pfam names on which sort is based */
              strcpy(pfam_id1,pfam1_ptr->pfam_id);
              strcpy(pfam_id2,pfam2_ptr->pfam_id);

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(pfam_id1,pfam_id2) > 0)
                {
                  swap_pfam_ptr = pfam_index_ptr[indi];
                  pfam_index_ptr[indi] = pfam_index_ptr[indl];
                  pfam_index_ptr[indl] = swap_pfam_ptr;

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

create_pfam_index  -  Create index of PFAM records for use in binary
                      search

***********************************************************************/

void create_pfam_index(struct pfam *first_pfam_ptr,
                       struct pfam ***pfam_ind_ptr,int npfam)
{
  int nentry;

  struct pfam *pfam_ptr;

  struct pfam **pfam_index_ptr;

  /* Initialise variables */

  /* Create the index array as one large array */
  pfam_index_ptr
    = (struct pfam **) malloc(npfam * sizeof(struct pfam *));

  /* Check that sufficient memory was allocated */
  if (pfam_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Get pointer to the first record */
  pfam_ptr = first_pfam_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (pfam_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      if (nentry < npfam)
        pfam_index_ptr[nentry] = pfam_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Sort the pfams */
  if (npfam > 1)
    pfam_sort(pfam_index_ptr,npfam);

  /* Return the array pointer */
  *pfam_ind_ptr = pfam_index_ptr;
}
/***********************************************************************

binary_pfam_search  -  Search for the given PFAM code

***********************************************************************/

struct pfam *binary_pfam_search(struct pfam **pfam_index_ptr,int npfam,
                              char *pfam_id)
{
  int bottom_point, ipoint, new_point, top_point;
  int found;

  struct pfam *pfam_ptr, *found_pfam_ptr, *prev_pfam_ptr;

  /* Initialise variables */
  found = FALSE;
  found_pfam_ptr = NULL;

  /* Make starting point half-way into the index */
  ipoint = npfam / 2;
  bottom_point = 0;
  top_point = npfam - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated pointer */
      pfam_ptr = pfam_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(pfam_id,pfam_ptr->pfam_id))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_pfam_ptr = pfam_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(pfam_id,pfam_ptr->pfam_id) < 0)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's name
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_pfam_ptr = pfam_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_pfam_ptr->pfam_id,pfam_id) < 0)
                found = TRUE;

              /* Otherwise, need to search earlier in the list */
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

      /* Otherwise, have to look further up the list */
      else
        {
          /* Save the current position as the lowest position to search
             on */
          bottom_point = ipoint;

          /* If have reached the top of the list, then our value is
             higher than any value in the list */
          if (ipoint >= npfam - 1)
            {
              /* Set flag as done */
              found = TRUE;
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

  /* Return the pointer */
  return(found_pfam_ptr);
}
/***********************************************************************

add_old_pfamlig_data  -  Pick up the existing data in the lig2pfam.dat
                         file and create a record for any data we do not
                         have

***********************************************************************/

int add_old_pfamlig_data(struct pfamlig **fst_pfamlig_ptr,
                         struct pdb *first_pdb_ptr,int npdb,
                         struct pfam *first_pfam_ptr,int npfam,
                         struct code **code_index_ptr,int ncodes,
                         struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, ligand_id[7], ligand_name[LINELEN + 1];
  char number_string[20];
  char domain_list[LINELEN + 1], pfam_id[LINELEN + 1];
  char uniprot_acc[UNIPROT_LEN], uniprot_id[UNIPROT_LEN];
  char last_pdb[PFAM_LEN], pdb_code[5];

  char *string_ptr;

  int domain_no, len, line, ncontacts, npfamlig, seq_no, type;
  int max_contacts, multi_domain_no;

  struct code *code_ptr;
  struct pdb *pdb_ptr;
  struct pfam *pfam_ptr;
  struct pfamlig *pfamlig_ptr, *first_pfamlig_ptr, *last_pfamlig_ptr;
  struct plist *plist_ptr, *first_plist_ptr, *last_plist_ptr;

  struct pdb **pdb_index_ptr;
  struct pfam **pfam_index_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  chain = '\0';
  domain_list[0] = '\0';
  domain_no = 0;
  last_pdb[0] = '\0';
  ligand_name[0] = '\0';
  ligand_id[0] = '\0';
  max_contacts = 0;
  multi_domain_no = 0;
  ncontacts = 0;
  npfamlig = 0;
  pdb_code[0] = '\0';
  pfam_id[0] = '\0';
  seq_no = 0;
  type = SINGLE;
  uniprot_acc[0] = '\0';
  uniprot_id[0] = '\0';

  /* Initialise pointers */
  last_pfamlig_ptr = NULL;
  first_plist_ptr = last_plist_ptr = NULL;

  /* Get pointer to first pfamlig record */
  pfamlig_ptr = first_pfamlig_ptr = *fst_pfamlig_ptr;

  /* Loop through pfamlig records to count them and get pointer to last
     one */
  while (pfamlig_ptr != NULL)
    {
      /* Store this record as the last one */
      last_pfamlig_ptr = pfamlig_ptr;

      /* Increment count */
      npfamlig++;

      /* Get next pfamlig record */
      pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;
    }

  /* Create an index of PDB codes to identify codes which have just
     been processed */
  create_pdb_index(first_pdb_ptr,&pdb_index_ptr,npdb);

  /* Create an index of Pfam domains */
  create_pfam_index(first_pfam_ptr,&pfam_index_ptr,npfam);

  /* Form name of input file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,DIR_SEP);
  strcat(file_name,"lig2pfam.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** Warning. File not found: %s\n",file_name);
      return(npfamlig);
    }

  /* Loop through the file processing each record */
  while (fgets(input_line,LINELEN,file_ptr)!= NULL)
    {
      /* Increment line-count */
      line++;

      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Separate the field name from the field value */
      string_ptr = strstr(input_line," ");
      if (string_ptr == NULL)
        {
          strcat(input_line," ");
          string_ptr = strstr(input_line," ");
        }

      /* If this is the KEY record, re-initialise */
      if (!strncmp(input_line,"KEY: ",5))
        type = SINGLE;

      /* If this is the ligand name, then store */
      else if (!strncmp(input_line,"L2P_LIGAND_NAME",15))
        strcpy(ligand_name,string_ptr+1);

      /* If this is the domain type, then store */
      else if (!strncmp(input_line,"L2P_DOMAIN_TYPE",15))
        {
          /* Get the domain type */
          strcpy(number_string,string_ptr+1);
          type = atoi(number_string);

          /* Re-initialise pointers */
          first_plist_ptr = last_plist_ptr = NULL;
        }

      /* Store Pfam id */
      else if (!strncmp(input_line,"L2P_PFAM_ID",11))
        {
          strcpy(pfam_id,string_ptr+1);
          strcpy(domain_list,string_ptr+1);
        }

      /* Pfam domain belonging to multi-domain */
      else if (!strncmp(input_line,"L2P_MULTI_DOMAIN[",17))
        strcpy(pfam_id,string_ptr+1);

      /* Domain number in multi-domain */
      else if (!strncmp(input_line,"L2P_MULTI_DOMAIN_NO",19))
        {
          /* Get the domain number */
          strcpy(number_string,string_ptr+1);
          multi_domain_no = atoi(number_string);
        }

      /* Number of this domain's contacts */
      else if (!strncmp(input_line,"L2P_NMULTI_CONTACTS",19))
        {
          /* Get the number of contacts */
          strcpy(number_string,string_ptr+1);
          ncontacts = atoi(number_string);

          /* Identify the corresponding Pfam record */
          pfam_ptr = binary_pfam_search(pfam_index_ptr,npfam,pfam_id);

          /* If have a valid Pfam record, then create plist
             record */
          if (pfam_ptr != NULL)
            {
              /* Create a Pfam domain list entry */
              plist_ptr
                = create_plist_entry(&first_plist_ptr,&last_plist_ptr,
                                     pfam_ptr,multi_domain_no,
                                     ncontacts); 
            }

          /* Otherwise, show warning message */
          else
            {
              if (pfam_ptr == NULL)
                printf("*** Warning. Unable to locate Pfam "
                       "domain [%s] ... (1)\n",pfam_id);
            }
        }

      /* Number of contacts from max-domain */
      else if (!strncmp(input_line,"L2P_MAX_NCONTACTS",17))
        {
          strcpy(number_string,string_ptr+1);
          max_contacts = atoi(number_string);
        }

      /* Domain making most contacts to ligand */
      else if (!strncmp(input_line,"L2P_MAX_CONTACT_DOMAIN",22))
        strcpy(pfam_id,string_ptr+1);

      /* Store corresponding UniProt accession code */
      else if (!strncmp(input_line,"L2P_UNIPROT_ACC",15))
        strcpy(uniprot_acc,string_ptr+1);

      /* Store corresponding UniProt id */
      else if (!strncmp(input_line,"L2P_UNIPROT_ID",14))
        strcpy(uniprot_id,string_ptr+1);

      /* Store PDB code */
      else if (!strncmp(input_line,"L2P_PDB_ID",10))
        strcpy(pdb_code,string_ptr+1);

      /* Store chain identifier */
      else if (!strncmp(input_line,"L2P_CHAIN",9))
        {
          if (string_ptr[1] != '\0')
            chain = string_ptr[1];
          else
            chain = ' ';
        }

      /* Store PDBsum sequence identier */
      else if (!strncmp(input_line,"L2P_SEQ_NO",10))
        {
          strcpy(number_string,string_ptr+1);
          seq_no = atoi(number_string);
        }

      /* Store sequential Pfam domain number within PDB entry */
      else if (!strncmp(input_line,"L2P_DOMAIN_NO",13))
        {
          strcpy(number_string,string_ptr+1);
          domain_no = atoi(number_string);
        }

      /* Store ligand residue number */
      else if (!strncmp(input_line,"L2P_RES_NUM",11))
        {
          strcpy(ligand_id," ");
          strcat(ligand_id,string_ptr+1);

          /* Check for missing insertion code */
          if (ligand_id[5] == '\0')
            strcat(ligand_id," ");
        }

      /* Store ligand chain id */
      else if (!strncmp(input_line,"L2P_LIG_CHAIN",13))
        {
          if (string_ptr[1] != '\0')
            ligand_id[0] = string_ptr[1];
        }

      /* Store number of ligand-domain contacts */
      else if (!strncmp(input_line,"L2P_NCONTACTS",13))
        {
          /* Store the number of contacts */
          strcpy(number_string,string_ptr+1);
          ncontacts = atoi(number_string);

          /* Check whether we have data relating to the current PDB
             code */
          pdb_ptr = binary_pdb_search(pdb_index_ptr,npdb,pdb_code);

          /* If we don't already have this PDB code, then create
             a new pfamlig record */
          if (pdb_ptr == NULL)
            {
              /* Identify the corresponding Pfam record */
              pfam_ptr
                = binary_pfam_search(pfam_index_ptr,npfam,pfam_id);

              /* Identify the corresponding UniProt record */
              code_ptr
                = binary_search(code_index_ptr,ncodes,uniprot_acc);

              /* If have a valid Pfam record, then create pfamlig
                 record */
              if (pfam_ptr != NULL && code_ptr != NULL)
                {
                  /* Create a new pfamlig entry */
                  pfamlig_ptr
                    = create_pfamlig_record(&first_pfamlig_ptr,
                                            &last_pfamlig_ptr,
                                            ligand_name,ligand_id,
                                            pdb_code,chain,seq_no,
                                            pfam_ptr,code_ptr,domain_no,
                                            type);

                  /* Store the domain name in the domain list */
                  store_string(&(pfamlig_ptr->domain_list),domain_list);

                  /* Store number of contacts */
                  pfamlig_ptr->ncontacts = ncontacts;

                  /* If this is a multi-domain entry, store pointer
                     to the first of the plist records */
                  if (type == MULTIPLE)
                    {
                      pfamlig_ptr->first_plist_ptr = first_plist_ptr;
                      pfamlig_ptr->max_contacts = max_contacts;
                    }

                  /* Increment count of links */
                  npfamlig++;
                }

              /* Otherwise, show warning message */
              else
                {
                  if (pfam_ptr == NULL && code_ptr == NULL)
                    printf("*** Warning. Unable to locate UniProt "
                           "code [%s] and Pfam domain [%s]\n",
                           uniprot_acc,pfam_id);
                  else if (code_ptr == NULL)
                    printf("*** Warning. Unable to locate UniProt "
                           "code [%s]\n",uniprot_acc);
                  else
                    printf("*** Warning. Unable to locate Pfam "
                           "domain [%s] ... (2)\n",
                           pfam_id);
                }
            }
        }
    }

  /* Close file */
  fclose(file_ptr);

  /* Return pointer to first pfamlig record */
  *fst_pfamlig_ptr = first_pfamlig_ptr;

  /* Return number of pfamlig records */
  return(npfamlig);
}
/***********************************************************************

write_lig2pfam_dat  -  Write out the lig2pfam.dat file listing ligands
                       in the PDB and the Pfam domains to which they are
                       bound

***********************************************************************/

void write_lig2pfam_dat(struct pfamlig *first_pfamlig_ptr)
{
  char file_name[FILENAME_LEN];
  char last_chain;
  char last_list[LINELEN], last_name[LINELEN], last_pdb[5];
  char last_uniprot_acc[UNIPROT_LEN];

  int last_domain_no, ncode, ndomain, ndomain_copy, nlist, npair, npdb;
  int i, last_no, last_type, n;

  struct code *code_ptr;
  struct pfam *pfam_ptr;
  struct pfamlig *pfamlig_ptr;
  struct plist *plist_ptr;

  FILE *file_ptr, *txtfile_ptr;

  /* Initialise variables */
  last_chain = '\0';
  last_list[0] = '\0';
  last_domain_no = -1;
  last_name[0] = '\0';
  last_pdb[0] = '\0';
  last_type = -1;
  last_uniprot_acc[0] = '\0';
  ncode = 0;
  ndomain = 0;
  ndomain_copy = 0;
  npair = 0;
  npdb = 0;

  /* Form the name of output file */
  strcpy(file_name,"lig2pfam.dat");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Form the name of output file */
  strcpy(file_name,"lig2pfam.txt");

  /* Open the output file */
  txtfile_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to first pfamlig record */
  pfamlig_ptr = first_pfamlig_ptr;

  /* Loop through pfamlig records to write out */
  while (pfamlig_ptr != NULL)
    {
      /* Get the corresponding Pfam and UniProt pointers */
      code_ptr = pfamlig_ptr->code_ptr;
      pfam_ptr = pfamlig_ptr->pfam_ptr;

      /* If new ligand or domain or PDB code or ligand instance,
         write out count of ligand-Pfam domain pairings in last PDB
         entry */
      if (strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(code_ptr->uniprot_acc,last_uniprot_acc) ||
          strcmp(pfamlig_ptr->pdb_code,last_pdb) ||
          pfamlig_ptr->chain != last_chain ||
          pfamlig_ptr->domain_no != last_domain_no)
        {
          /* Write out number of pairs */
          if (npair > 0)
            fprintf(file_ptr,"L2P_NPAIRS[%d][%d][%d][%d] %d\n",ndomain,
                    ncode,npdb,ndomain_copy,npair);

          /* Re-initialise count */
          npair = 0;
        }

      /* If new ligand or domain or PDB code, write out number of
         copies of the domain involved in ligand interactions in last
         PDB entry */
      if (strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(code_ptr->uniprot_acc,last_uniprot_acc) ||
          strcmp(pfamlig_ptr->pdb_code,last_pdb) ||
          pfamlig_ptr->chain != last_chain)
        {
          /* Write out number of domain copies */
          if (ndomain_copy > 0)
            fprintf(file_ptr,"L2P_NDOMAIN_COPY[%d][%d][%d] %d\n",ndomain,
                    ncode,npdb,ndomain_copy);

          /* Re-initialise count */
          ndomain_copy = 0;
        }

      /* If new ligand or domain, write out number of PDB entries
         that contained last domain */
      if (strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(code_ptr->uniprot_acc,last_uniprot_acc))
        {
          /* Write out number of PDB entries */
          if (npdb > 0)
            fprintf(file_ptr,"L2P_NPDB[%d][%d] %d\n",ndomain,ncode,npdb);

          /* Re-initialise count */
          npdb = 0;
        }

      /* If new ligand or domain, write out number of UniProt code
         that contained last domain */
      if (strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(pfamlig_ptr->domain_list,last_list))
        {
          /* Write out number of UniProt entries */
          if (ncode > 0)
            fprintf(file_ptr,"L2P_NUNIPROT[%d] %d\n",ndomain,ncode);

          /* Re-initialise count */
          ncode = 0;
        }

      /* If this is a new ligand, re-initialise */
      if (strcmp(pfamlig_ptr->ligand_name,last_name))
        {
          /* Close off previous ligand */
          if (ndomain > 0)
            fprintf(file_ptr,"L2P_NDOMAINS %d\n",ndomain);

          /* Write out key record */
          fprintf(file_ptr,"KEY: %s\n",pfamlig_ptr->ligand_name);

          /* Write out ligand name */
          fprintf(file_ptr,"L2P_LIGAND_NAME %s\n",
                  pfamlig_ptr->ligand_name);

          /* Re-initialise variables */
          ndomain = 0;
          last_list[0] = '\0';

          /* Save the ligand name */
          strcpy(last_name,pfamlig_ptr->ligand_name);
        }

      /* If this is a new Pfam domain, write out code */
      if (strcmp(pfamlig_ptr->domain_list,last_list))
        {
          /* Increment domain count */
          ndomain++;

          /* Write out whether domain is a single domain, part of a
             multiple-domain interaction group, or represents a
             multiple-domain interaction group */
          fprintf(file_ptr,"L2P_DOMAIN_TYPE[%d] %d\n",ndomain,
                  pfamlig_ptr->type);

          /* Write out this Pfam domain's id */
          fprintf(file_ptr,"L2P_PFAM_ID[%d] %s\n",ndomain,
                  pfamlig_ptr->domain_list);

          /* If not a multiple-domain, write out name and description */
          if (pfamlig_ptr->type != MULTIPLE)
            {
              /* Write out its name */
              fprintf(file_ptr,"L2P_PFAM_NAME[%d] %s\n",ndomain,
                      pfam_ptr->name);

              /* Write out its description */
              if (pfam_ptr->description != &null)
                fprintf(file_ptr,"L2P_PFAM_DESC[%d] %s\n",ndomain,
                        pfam_ptr->description);
              else
                fprintf(file_ptr,"L2P_PFAM_DESC[%d] NONE\n",ndomain);
            }

          /* If a multiple-domain, write out the interacting domains */
          if (pfamlig_ptr->type == MULTIPLE)
            {
              /* Initialise domain count */
              nlist = 0;

              /* Get first of the Pfam domains */
              plist_ptr = pfamlig_ptr->first_plist_ptr;

              /* Loop through the Pfam domains to list */
              while (plist_ptr != NULL)
                {
                  /* Increment count */
                  nlist++;

                  /* Get the corresponding Pfam domains */
                  pfam_ptr = plist_ptr->pfam_ptr;

                  /* Write out this domain */
                  fprintf(file_ptr,"L2P_MULTI_DOMAIN[%d][%d] %s\n",
                          ndomain,nlist,pfam_ptr->pfam_id);

                  /* Write out domain number */
                  fprintf(file_ptr,"L2P_MULTI_DOMAIN_NO[%d][%d] %d\n",
                          ndomain,nlist,plist_ptr->domain_no);

                  /* Write out number of contacts */
                  fprintf(file_ptr,"L2P_NMULTI_CONTACTS[%d][%d] %d\n",
                          ndomain,nlist,plist_ptr->ncontacts);

                  /* Get next plist record */
                  plist_ptr = plist_ptr->next_plist_ptr;
                }

              /* Write out count of domains */
              fprintf(file_ptr,"L2P_NMULTI_DOMAINS[%d] %d\n",ndomain,
                      nlist);

              /* Write out the number of contacts made by main domain */
              fprintf(file_ptr,"L2P_MAX_NCONTACTS[%d] %d\n",ndomain,
                      pfamlig_ptr->max_contacts);

              /* Write out the domain that makes the most contacts
                 to the ligand */
              fprintf(file_ptr,"L2P_MAX_CONTACT_DOMAIN[%d] %s\n",ndomain,
                      pfamlig_ptr->pfam_ptr->pfam_id);
            }

          /* Initialise count of UniProt codes in which this ligand-domain
             combination is found */
          ncode = 0;
          last_uniprot_acc[0] = '\0';

          /* Save the domain id */
          strcpy(last_list,pfamlig_ptr->domain_list);
        }

      /* If this is a new UniProt code, write out id and accession number */
      if (strcmp(code_ptr->uniprot_acc,last_uniprot_acc))
        {
          /* Increment code count */
          ncode++;

          /* Write out this UniProt accession number */
          fprintf(file_ptr,"L2P_UNIPROT_ACC[%d][%d] %s\n",ndomain,ncode,
                  code_ptr->uniprot_acc);

          /* Write out this UniProt id */
          fprintf(file_ptr,"L2P_UNIPROT_ID[%d][%d] %s\n",ndomain,ncode,
                  code_ptr->uniprot_id);

          /* Initialise count of PDB entries */
          npdb = 0;
          last_pdb[0] = '\0';

          /* Save the UniProt code */
          strcpy(last_uniprot_acc,code_ptr->uniprot_acc);
        }

      /* If this is a new PDB entry, write out code */
      if (strcmp(pfamlig_ptr->pdb_code,last_pdb))
        {
          /* Increment count of PDB entries */
          npdb++;

          /* Write out the PDB entry */
          fprintf(file_ptr,"L2P_PDB_ID[%d][%d][%d] %s\n",ndomain,ncode,
                  npdb,pfamlig_ptr->pdb_code);

          /* Initialise count of domain copies in this PDB entry */
          ndomain_copy = 0;
          last_chain = '\0';
          last_domain_no = -1;

          /* Save the PDB code */
          strcpy(last_pdb,pfamlig_ptr->pdb_code);
        }

      /* If this is a new domain or chain, write out chain id and
         domain number */
      if (pfamlig_ptr->domain_no != last_domain_no ||
          pfamlig_ptr->chain != last_chain)
        {
          /* Increment count of domain copies*/
          ndomain_copy++;

          /* Write out the domain number */
          fprintf(file_ptr,"L2P_CHAIN[%d][%d][%d][%d] %c\n",ndomain,
                  ncode,npdb,ndomain_copy,pfamlig_ptr->chain);
          fprintf(file_ptr,"L2P_SEQ_NO[%d][%d][%d][%d] %d\n",ndomain,
                  ncode,npdb,ndomain_copy,pfamlig_ptr->seq_no);
          fprintf(file_ptr,"L2P_DOMAIN_NO[%d][%d][%d][%d] %d\n",ndomain,
                  ncode,npdb,ndomain_copy,pfamlig_ptr->domain_no);

          /* Initialise pair count */
          npair = 0;

          /* Save the last domain number and chain id */
          last_domain_no = pfamlig_ptr->domain_no;
          last_chain = pfamlig_ptr->chain;
        }

      /* Increment ligand-domain pairing count */
      npair++;

      /* Write out the chain and residue number for the current ligand */
      fprintf(file_ptr,"L2P_RES_NUM[%d][%d][%d][%d][%d] %s\n",ndomain,
              ncode,npdb,ndomain_copy,npair,pfamlig_ptr->ligand_id+1);
      fprintf(file_ptr,"L2P_LIG_CHAIN[%d][%d][%d][%d][%d] %c\n",ndomain,
              ncode,npdb,ndomain_copy,npair,pfamlig_ptr->ligand_id[0]);

      /* Write out the number of ligand-domain contacts */
      fprintf(file_ptr,"L2P_NCONTACTS[%d][%d][%d][%d][%d] %d\n",ndomain,
              ncode,npdb,ndomain_copy,npair,pfamlig_ptr->ncontacts);

      /* Get the Pfam domain pointer */
      pfam_ptr = pfamlig_ptr->pfam_ptr;

      /* Write out this record to the .txt file */
      fprintf(txtfile_ptr,"%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%c\t"
              "%s\t%c\t%d\t%c\t%s\n",
              pfamlig_ptr->ligand_name,ndomain,pfam_ptr->pfam_id,
              pfam_ptr->name,ncode,code_ptr->uniprot_acc,
              code_ptr->uniprot_id,npdb,pfamlig_ptr->pdb_code,
              pfamlig_ptr->chain,pfamlig_ptr->ligand_id+1,
              pfamlig_ptr->ligand_id[0],pfamlig_ptr->ncontacts,
              DOMAIN_TYPE[pfamlig_ptr->type],pfamlig_ptr->domain_list);

      /* Get next pfamlig record */
      pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;
    }

  /* Close off last domain and ligand's details */
  if (npair > 0)
    fprintf(file_ptr,"L2P_NPAIRS[%d][%d][%d][%d] %d\n",ndomain,ncode,
            npdb,ndomain_copy,npair);
  if (ndomain_copy > 0)
    fprintf(file_ptr,"L2P_NDOMAIN_COPY[%d][%d][%d] %d\n",ndomain,ncode,
            npdb,ndomain_copy);
  if (npdb > 0)
    fprintf(file_ptr,"L2P_NPDB[%d][%d] %d\n",ndomain,ncode,npdb);
  if (ncode > 0)
    fprintf(file_ptr,"L2P_NUNIPROT[%d] %d\n",ndomain,ncode);
  if (ndomain > 0)
    fprintf(file_ptr,"L2P_NDOMAINS %d\n",ndomain);

  /* Close the files */
  fclose(file_ptr);
  fclose(txtfile_ptr);
}
/***********************************************************************

write_index  -  Read through the given file to pick up the offset 
                positions of each key term and write out to index file

***********************************************************************/

int write_index(char *fname)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];

  int len, nkeys, nlines;

  long offset;

  FILE *file_ptr, *outfile_ptr;

  /* Initialise variables */
  nkeys = 0;
  nlines = 0;
  offset = 0;

  /* Form name of the input data file */
  strcpy(file_name,fname);
  strcat(file_name,".dat");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      printf("\n*** Unable to open data file [%s]\n",file_name);
      return(nkeys);
    }
    
  /* Form name of the output index file */
  strcpy(file_name,fname);
  strcat(file_name,".idx");

  /* Open the output file */
  outfile_ptr = open_output_file(file_name,TRUE);

  /* Get current file position */
  offset = ftell(file_ptr);

  /* Read through the input file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Determine whether this is a key record */
      if (!strncmp(input_line,"KEY: ",5))
        {
          /* Truncate the string */
          len = string_truncate(input_line,LINELEN);

          /* Write out offset to index file */
          fprintf(outfile_ptr,"%s %ld\n",input_line+5,offset);

          /* Increment count of index entries */
          nkeys++;
        }

      /* Get current file position */
      offset = ftell(file_ptr);
    }

  /* Close both files */
  fclose(file_ptr);
  fclose(outfile_ptr);

  /* Return number of index entries encountered */
  return(nkeys);
}
/***********************************************************************

write_pfam2lig_dat  -  Write out the pfam2lig.dat file listing Pfam
                       domains in the PDB and the ligands that bind to
                       them

***********************************************************************/

void write_pfam2lig_dat(struct pfamlig *first_pfamlig_ptr)
{
  char file_name[FILENAME_LEN];
  char last_chain;
  char last_list[LINELEN], last_name[LINELEN], last_pdb[5];
  char last_uniprot_acc[UNIPROT_LEN];

  int last_domain_no, ncode, nligand, ndomain_copy, npair, npdb;

  struct code *code_ptr;
  struct pfam *pfam_ptr;
  struct pfamlig *pfamlig_ptr;

  FILE *file_ptr, *txtfile_ptr;

  /* Initialise variables */
  last_chain = '\0';
  last_list[0] = '\0';
  last_domain_no = -1;
  last_name[0] = '\0';
  last_pdb[0] = '\0';
  last_uniprot_acc[0] = '\0';
  ncode = 0;
  nligand = 0;
  ndomain_copy = 0;
  npair = 0;
  npdb = 0;

  /* Form the name of output file */
  strcpy(file_name,"pfam2lig.dat");

  /* Open the output file */
  file_ptr = open_output_file(file_name,TRUE);

  /* Form the name of output file */
  strcpy(file_name,"pfam2lig.txt");

  /* Open the output file */
  txtfile_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to first pfamlig record */
  pfamlig_ptr = first_pfamlig_ptr;

  /* Loop through pfamlig records to write out */
  while (pfamlig_ptr != NULL)
    {
      /* Get the corresponding Pfam and UniProt pointers */
      code_ptr = pfamlig_ptr->code_ptr;
      pfam_ptr = pfamlig_ptr->pfam_ptr;

      /* If new ligand or domain or PDB code or ligand instance,
         write out count of ligand-Pfam domain pairings in last PDB
         entry */
      if (strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(code_ptr->uniprot_acc,last_uniprot_acc) ||
          strcmp(pfamlig_ptr->pdb_code,last_pdb) ||
          pfamlig_ptr->chain != last_chain ||
          pfamlig_ptr->domain_no != last_domain_no)
        {
          /* Write out number of pairs */
          if (npair > 0)
            fprintf(file_ptr,"P2L_NPAIRS[%d][%d][%d][%d] %d\n",nligand,
                    ncode,npdb,ndomain_copy,npair);

          /* Re-initialise count */
          npair = 0;
        }

      /* If new ligand or domain or PDB code, write out number of
         copies of the domain involved in ligand interactions in last
         PDB entry */
      if (strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(code_ptr->uniprot_acc,last_uniprot_acc) ||
          strcmp(pfamlig_ptr->pdb_code,last_pdb) ||
          pfamlig_ptr->chain != last_chain)
        {
          /* Write out number of domain copies */
          if (ndomain_copy > 0)
            fprintf(file_ptr,"P2L_NDOMAIN_COPY[%d][%d][%d] %d\n",nligand,
                    ncode,npdb,ndomain_copy);

          /* Re-initialise count */
          ndomain_copy = 0;
        }

      /* If new ligand or domain, write out number of PDB entries
         that contained last domain */
      if (strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(pfamlig_ptr->ligand_name,last_name) ||
          strcmp(code_ptr->uniprot_acc,last_uniprot_acc))
        {
          /* Write out number of PDB entries */
          if (npdb > 0)
            fprintf(file_ptr,"P2L_NPDB[%d][%d] %d\n",nligand,ncode,npdb);

          /* Re-initialise count */
          npdb = 0;
        }

      /* If new ligand or domain, write out number of UniProt codes
         that contained last domain */
      if (strcmp(pfamlig_ptr->domain_list,last_list) ||
          strcmp(pfamlig_ptr->ligand_name,last_name))
        {
          /* Write out number of UniProt entries */
          if (ncode > 0)
            fprintf(file_ptr,"P2L_NUNIPROT[%d] %d\n",nligand,ncode);

          /* Re-initialise count */
          ncode = 0;
        }

      /* If this is a new Pfam domain, re-initialise */
      if (strcmp(pfamlig_ptr->domain_list,last_list))
        {
          /* Close off previous ligand */
          if (nligand > 0)
            fprintf(file_ptr,"P2L_NLIGANDS %d\n",nligand);

          /* Write out key record */
          fprintf(file_ptr,"KEY: %s\n",pfamlig_ptr->domain_list);

          /* Write out whether domain is a single domain, part of a
             multiple-domain interaction group, or represents a
             multiple-domain interaction group */
          fprintf(file_ptr,"P2L_DOMAIN_TYPE %d\n",pfamlig_ptr->type);

          /* Write out Pfam id */
          fprintf(file_ptr,"P2L_PFAM_ID %s\n",pfamlig_ptr->domain_list);

          /* If not a multiple-domain, write out domain name */
          if (pfamlig_ptr->type != MULTIPLE)
            {
              /* Write out its name */
              fprintf(file_ptr,"P2L_PFAM_NAME %s\n",pfam_ptr->name);

              /* Write out its description */
              if (pfam_ptr->description != &null)
                fprintf(file_ptr,"P2L_PFAM_DESC %s\n",
                        pfam_ptr->description);
              else
                fprintf(file_ptr,"P2L_PFAM_DESC NONE\n");
            }

          /* Otherwise, list the individual domains */
          else
            {

/* @@@ */
            }

          /* Re-initialise variables */
          nligand = 0;
          last_list[0] = '\0';

          /* Save the ligand name */
          strcpy(last_list,pfamlig_ptr->domain_list);
        }

      /* If this is a new ligand, write out */
      if (strcmp(pfamlig_ptr->ligand_name,last_name))
        {
          /* Increment domain count */
          nligand++;

          /* Write out this Pfam domain's id */
          fprintf(file_ptr,"P2L_LIGAND_NAME[%d] %s\n",nligand,
                  pfamlig_ptr->ligand_name);

          /* Initialise count of UniProt codes in which this ligand-domain
             combination is found */
          ncode = 0;
          last_uniprot_acc[0] = '\0';

          /* Save the domain id */
          strcpy(last_name,pfamlig_ptr->ligand_name);
        }

      /* If this is a new UniProt code, write out id and accession number */
      if (strcmp(code_ptr->uniprot_acc,last_uniprot_acc))
        {
          /* Increment code count */
          ncode++;

          /* Write out this UniProt accession number */
          fprintf(file_ptr,"P2L_UNIPROT_ACC[%d][%d] %s\n",nligand,ncode,
                  code_ptr->uniprot_acc);

          /* Write out this UniProt id */
          fprintf(file_ptr,"P2L_UNIPROT_ID[%d][%d] %s\n",nligand,ncode,
                  code_ptr->uniprot_id);

          /* Initialise count of PDB entries */
          npdb = 0;
          last_pdb[0] = '\0';

          /* Save the UniProt code */
          strcpy(last_uniprot_acc,code_ptr->uniprot_acc);
        }

      /* If this is a new PDB entry, write out code */
      if (strcmp(pfamlig_ptr->pdb_code,last_pdb))
        {
          /* Increment count of PDB entries */
          npdb++;

          /* Write out the PDB entry */
          fprintf(file_ptr,"P2L_PDB_ID[%d][%d][%d] %s\n",nligand,
                  ncode,npdb,pfamlig_ptr->pdb_code);

          /* Initialise count of domain copies in this PDB entry */
          ndomain_copy = 0;
          last_chain = '\0';
          last_domain_no = -1;

          /* Save the PDB code */
          strcpy(last_pdb,pfamlig_ptr->pdb_code);
        }

      /* If this is a new domain or chain, write out chain id and
         domain number */
      if (pfamlig_ptr->domain_no != last_domain_no ||
          pfamlig_ptr->chain != last_chain)
        {
          /* Increment count of domain copies*/
          ndomain_copy++;

          /* Write out the domain number */
          fprintf(file_ptr,"P2L_CHAIN[%d][%d][%d][%d] %c\n",nligand,
                  ncode,npdb,ndomain_copy,pfamlig_ptr->chain);
          fprintf(file_ptr,"P2L_SEQ_NO[%d][%d][%d][%d] %d\n",nligand,
                  ncode,npdb,ndomain_copy,pfamlig_ptr->seq_no);
          fprintf(file_ptr,"P2L_DOMAIN_NO[%d][%d][%d][%d] %d\n",nligand,
                  ncode,npdb,ndomain_copy,pfamlig_ptr->domain_no);

          /* Initialise pair count */
          npair = 0;

          /* Save the last domain number and chain id */
          last_domain_no = pfamlig_ptr->domain_no;
          last_chain = pfamlig_ptr->chain;
        }

      /* Increment ligand-domain pairing count */
      npair++;

      /* Write out the chain and residue number for the current ligand */
      fprintf(file_ptr,"P2L_RES_NUM[%d][%d][%d][%d][%d] %s\n",nligand,
              ncode,npdb,ndomain_copy,npair,pfamlig_ptr->ligand_id+1);
      fprintf(file_ptr,"P2L_LIG_CHAIN[%d][%d][%d][%d][%d] %c\n",nligand,
              ncode,npdb,ndomain_copy,npair,pfamlig_ptr->ligand_id[0]);

      /* Write out the number of ligand-domain contacts */
      fprintf(file_ptr,"P2L_NCONTACTS[%d][%d][%d][%d][%d] %d\n",nligand,
              ncode,npdb,ndomain_copy,npair,pfamlig_ptr->ncontacts);

      /* Get the Pfam pointer */
      pfam_ptr = pfamlig_ptr->pfam_ptr;

      /* Write out this record to the .txt file */
      fprintf(txtfile_ptr,"%s\t%c\t%d\t%s\t%d\t%s\t%s\t%d\t%s\t%c\t"
              "%s\t%c\t%s\t%s\t%d\n",
              pfamlig_ptr->domain_list,DOMAIN_TYPE[pfamlig_ptr->type],
              nligand,pfamlig_ptr->ligand_name,ncode,code_ptr->uniprot_acc,
              code_ptr->uniprot_id,npdb,pfamlig_ptr->pdb_code,
              pfamlig_ptr->chain,pfamlig_ptr->ligand_id+1,
              pfamlig_ptr->ligand_id[0],pfam_ptr->pfam_id,pfam_ptr->name,
              pfamlig_ptr->ncontacts);

      /* Get next pfamlig record */
      pfamlig_ptr = pfamlig_ptr->next_pfamlig_ptr;
    }

  /* Close off last domain and ligand's details */
  if (npair > 0)
    fprintf(file_ptr,"P2L_NPAIRS[%d][%d][%d][%d] %d\n",nligand,
            ncode,npdb,ndomain_copy,npair);
  if (ndomain_copy > 0)
    fprintf(file_ptr,"P2L_NDOMAIN_COPY[%d][%d][%d] %d\n",nligand,ncode,
            npdb,ndomain_copy);
  if (npdb > 0)
    fprintf(file_ptr,"P2L_NPDB[%d][%d] %d\n",nligand,ncode,npdb);
  if (ncode > 0)
    fprintf(file_ptr,"P2L_NUNIPROT[%d] %d\n",nligand,ncode);
  if (nligand > 0)
    fprintf(file_ptr,"P2L_NLIGANDS %d\n",nligand);

  /* Close the files */
  fclose(file_ptr);
  fclose(txtfile_ptr);
}
/***********************************************************************

write_ligpfam_analyses  -  Write out the ligand-Pfam domain analyses

***********************************************************************/

void write_ligpfam_analyses(struct pfamlig *first_pfamlig_ptr,
                            int npfamlig,struct pdb *first_pdb_ptr,
                            int npdb,struct pfam *first_pfam_ptr,
                            int npfam,struct code **code_index_ptr,
                            int ncodes,struct c_file_data file_name_ptr)
{
  /* Sort all the pfamlig entries by ligand name and Pfam domain */
  if (npfamlig > 0)
    sort_pfamlig_entries(&first_pfamlig_ptr,npfamlig,LIG_PFAM);

  /* Read in the data in the current lig2pfam.dat file and append to
     existing data */
  npfamlig = add_old_pfamlig_data(&first_pfamlig_ptr,first_pdb_ptr,
                                  npdb,first_pfam_ptr,npfam,
                                  code_index_ptr,ncodes,file_name_ptr);

  /* Sort once again to include all the entries just read in from
     file */
  if (npfamlig > 0)
    sort_pfamlig_entries(&first_pfamlig_ptr,npfamlig,LIG_PFAM);

  /* Write out the lig2pfam.dat file */
  write_lig2pfam_dat(first_pfamlig_ptr);

  /* Write out the index file, lig2pfam.idx */
  write_index("lig2pfam");

  /* Sort all the pfamlig entries by Pfam domain and ligand name */
  if (npfamlig > 0)
    sort_pfamlig_entries(&first_pfamlig_ptr,npfamlig,PFAM_LIG);

  /* Write out the pfam2lig.dat file */
  write_pfam2lig_dat(first_pfamlig_ptr);

  /* Write out the index file, lig2pfam.idx */
  write_index("pfam2lig");
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char dataset_name[FILENAME_LEN], character[NLETTERS + 1];

  int ncodes, nhetdrug, nletters, npdb, npfam, npfamlig;
  int dos, entry_types, os_type;

  struct c_file_data file_name_ptr;
  struct letter letter_ptr[NLETTERS];

  struct code *first_code_ptr, *last_code_ptr;
  struct hetdrug *first_hetdrug_ptr;
  struct pdb *first_pdb_ptr;
  struct pfam *first_pfam_ptr, *last_pfam_ptr;
  struct pfamlig *first_pfamlig_ptr;

  struct code **code_index_ptr;
  struct hetdrug **hetdrug_index_ptr;
  struct pfam **pfam_index_ptr;

  struct parameters params;

  FILE *donefile_ptr;

  /* Initialise variables */
  entry_types = FALSE;
  nhetdrug = npfam = 0;

  /* Initialise pointers */
  hetdrug_index_ptr = NULL;
  first_hetdrug_ptr = NULL;
  first_pfam_ptr = last_pfam_ptr = NULL;

  /* Read in the command-line parameters and store */
  get_command_arguments(argv,argc - 1,&params);

  /* Read in all the directory paths from the parameter file, CATHPARAM */
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

  /* Get all the character definitions from the letters.ppm file */
  nletters = get_letters(character,letter_ptr,file_name_ptr);

  /* Read in the Het Group to DrugPort mappings from the appdrugs_pdb.dat
     file */
  if (params.pdb_type != PDBSUM1)
    nhetdrug = get_hetdrug_mappings(&first_hetdrug_ptr,file_name_ptr);

  /* Create index of Het Group -> drug mappings */
  if (nhetdrug > 0)
    create_hetdrug_index(first_hetdrug_ptr,&hetdrug_index_ptr,nhetdrug);

  /* Sort the hetdrug records */
  if (nhetdrug > 1)
    sort_hetdrug_entries(hetdrug_index_ptr,nhetdrug);

  /* Define input filename */
  strcpy(dataset_name,"dataset.lst");

  /* Read in the PDB codes in the dataset.lst file */
  npdb = get_pdb_codes(dataset_name,&first_pdb_ptr,params.pdb_code);

  /* Process each of the codes in turn */
  if (npdb > 0)
    {
      /* Process the PDB codes in turn to pick up the unique chains for
         each entry and any UniProt codes associated with these */
      ncodes
        = process_codes(first_pdb_ptr,&first_code_ptr,&last_code_ptr,
                        npdb,file_name_ptr,params);

      /* Read through the current lig2pfam.dat file to pick up the
         UniProt codes it contains */
      if (params.run_type == WEEKLY_UPDATE)
        ncodes = read_old_codes(&first_code_ptr,last_code_ptr,ncodes,
                                file_name_ptr);

      /* Create an index to allow binary search through UniProt sequences
         file */
      if (ncodes > 0)
        entry_types
          = create_code_index(first_code_ptr,&code_index_ptr,ncodes);

      /* If have any UniProt codes (other than the dummy one) search Pfam */
      if (ncodes > 1)
        {
          /* Pick up the Pfam domains from the swisspfam_data.txt.gz file */
          if (params.covid19 == FALSE && params.pdb_type != PDBSUM1)
            npfam
              = read_swisspfam_data(code_index_ptr,ncodes,&first_pfam_ptr,
                                    &last_pfam_ptr,entry_types,params,
                                    file_name_ptr);

          /* Read the swisspfam.gz file to pick up the Pfam domains */
          if (params.pdb_type != PDBSUM1)
            npfam = read_swisspfam(code_index_ptr,ncodes,&first_pfam_ptr,
                                   &last_pfam_ptr,entry_types,params,
                                   file_name_ptr);

          /* For any UniProt codes not in swisspfam.gz, get the 
             information from their .xml files */
          if (params.covid19 == FALSE && params.wget == TRUE)
            npfam = get_missing_domains(code_index_ptr,ncodes,
                                        &first_pfam_ptr,last_pfam_ptr,
                                        npfam,file_name_ptr);

          /* Sort the Pfam records */
          if (npfam > 0 && params.covid19 == FALSE)
            create_pfam_index(first_pfam_ptr,&pfam_index_ptr,npfam);

          /* Pick up the Pfam-A short names and descriptions from the
             pfam_domains.parsed.sorted file */
          if (params.covid19 == FALSE)
            get_pfam_names(pfam_index_ptr,npfam,file_name_ptr);
        }

      /* Open the done.lst file */
      donefile_ptr = open_output_file("done.lst",TRUE);

      /* Locate each UniProt entry in the UniProt sequence file */
      if (ncodes > 0)
        {
          /* Download the required sequences and align against the
             PDB sequences */
          download_uniprot_seqs(code_index_ptr,ncodes,params.pdb_type,
                                params.pdbsum_type,file_name_ptr);

          /* Process the alignments */
          process_alignments(code_index_ptr,ncodes,character,
                             letter_ptr,nletters,&first_pfamlig_ptr,
                             &npfamlig,donefile_ptr,first_code_ptr,
                             entry_types,hetdrug_index_ptr,nhetdrug,
                             file_name_ptr,params);

          /* Delete all the downloaded .fasta files */
          delete_fasta_files(code_index_ptr,ncodes,dos);

          /* Process the sequences */
          //else
          //  process_uniprot_seqs(code_index_ptr,ncodes,character,
          //                       letter_ptr,nletters,&first_pfamlig_ptr,
          //                       &npfamlig,donefile_ptr,first_code_ptr,
          //                       entry_types,hetdrug_index_ptr,nhetdrug,
          //                       file_name_ptr,params);
        }

      /* Process any remaining incomplete PDB entries */
      close_pdb_entries(first_pdb_ptr,donefile_ptr,hetdrug_index_ptr,
                        nhetdrug,file_name_ptr,params);

      /* Close the done.lst file */
      fclose(donefile_ptr);

      /* If weekly update, write out weekly files */
      if (params.run_type == WEEKLY_UPDATE)
        {
          /* Write out the ligand-Pfam domain analyses */
          write_ligpfam_analyses(first_pfamlig_ptr,npfamlig,first_pdb_ptr,
                                 npdb,first_pfam_ptr,npfam,code_index_ptr,
                                 ncodes,file_name_ptr);
        }
    }
/* debug */
  printf("END\n");
  fflush(stdout);
/* debug */

  return(0);
}
