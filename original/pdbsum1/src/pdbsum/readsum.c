/***********************************************************************

  readsum.c - Program to read in data from PDBsum data files

************************************************************************

Date:-         23 Nov 2009

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name  Description
-----------  -----------
database.dat Input file in PDBsum directory giving header information
             on PDBsum entry
grow.out     Input file in PDBsum directory giving protein-ligand,
             protein-metal, and protein-DNA interactions
igrow.out.gz Input file in PDBsum directory giving protein-protein
             interactions
resdata.dat  Input file in PDBsum directory giving residue annotations


Compilation
-----------

cc -c common.c
cc -c readsum.c
cc -o readsum readsum.o common.o -lm -lz

For testing, compile with TEST flag:  cc -DTEST readsum.c


*/

/* Function calling-tree
   ---------------------

main
   -> s_append_grow_records
   -> s_get_grow_data
         -> string_truncate (in common.c)
         -> create_s_grow_record
   -> s_free_grow_data
   -> s_init_sort
   -> s_sort_grow_records
         -> s_grow_sort
         -> store_string (in common.c)
   -> s_get_protein_name
         -> string_truncate (in common.c)
   -> s_get_resdata
         -> find_pdbsum_dir (in common.c)
         -> s_get_equiv_chain
         -> string_truncate (in common.c)
         -> create_r_residue_record (in readpdb.c)
   -> s_transfer_resdata

*/


#include "common.h"
#include "readpdb.h"
#include "readsum.h"

/* Prototypes
   ---------- */

/* In common.c */
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
int string_chomp(char *string);
int string_truncate(char *string,int max_length);
void store_string(char **string_ptr,char *string);

/* In readpdb.c */
struct r_residue *create_r_residue_record(struct r_residue **fst_r_residue_ptr,
                                          struct r_residue **lst_r_residue_ptr,
                                          struct r_pdb *r_pdb_ptr,char athet,
                                          struct r_chain *r_chain_ptr,
                                          char *res_name,char *res_num,
                                          struct r_rasdef *r_rasdef_ptr);
struct r_residue *r_find_residue(struct r_residue *first_r_residue_ptr,
                                 char *res_name,char *res_num,char chain);


/***********************************************************************

s_append_grow_records  -  Append the second set of grow records onto the
                          end of the first

***********************************************************************/

void s_append_grow_records(struct s_grow **fst_s_grow_ptr,
			   struct s_grow *second_s_grow_ptr)
{
  struct s_grow *s_grow_ptr, *first_s_grow_ptr, *last_s_grow_ptr;

  /* Initialise variables */
  first_s_grow_ptr = *fst_s_grow_ptr;
  last_s_grow_ptr = NULL;
  
  /* Get the first grow */
  s_grow_ptr = first_s_grow_ptr;

  /* Loop through the grow records to find the last */
  while (s_grow_ptr != NULL)
    {
      /* Save the last grow record */
      last_s_grow_ptr = s_grow_ptr;

      /* Get the next grow record */
      s_grow_ptr = s_grow_ptr->next_s_grow_ptr;
    }

  /* If we have a valid last record, get it to point to the second set
     of records */
  if (last_s_grow_ptr != NULL)
    last_s_grow_ptr->next_s_grow_ptr = second_s_grow_ptr;

  /* Otherwise, make the second set the first */
  else
    first_s_grow_ptr = second_s_grow_ptr;

  /* Return the new first grow pointer */
  *fst_s_grow_ptr = first_s_grow_ptr;
}
/***********************************************************************

create_s_grow_record  -  Create and initialise a new s_grow record

***********************************************************************/

struct s_grow *create_s_grow_record(struct s_grow **fst_s_grow_ptr,
                                    struct s_grow **lst_s_grow_ptr)
{
  int i;

  struct s_grow *s_grow_ptr, *first_s_grow_ptr, *last_s_grow_ptr;

  /* Initialise s_grow pointers */
  s_grow_ptr = NULL;
  first_s_grow_ptr = *fst_s_grow_ptr;
  last_s_grow_ptr = *lst_s_grow_ptr;

  /* Allocate memory for structure to hold s_grow info */
  s_grow_ptr = (struct s_grow *) malloc(sizeof(struct s_grow));
  if (s_grow_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct s_grow\n");
      exit (1);
    }

  /* If this is the very first, save its pointer */
  if (first_s_grow_ptr == NULL)
    first_s_grow_ptr = s_grow_ptr;

  /* Add link from previous s_grow to the current one */
  if (last_s_grow_ptr != NULL)
    last_s_grow_ptr->next_s_grow_ptr = s_grow_ptr;
  last_s_grow_ptr = s_grow_ptr;

  /* Initialise the elements of the structure */
  for (i = 0; i < 2; i++)
    {
      s_grow_ptr->res_name[i][0] = '\0';
      s_grow_ptr->res_num[i][0] = '\0';
      s_grow_ptr->chain[i] = ' ';
    }
  s_grow_ptr->type = S_HBONDS;
  s_grow_ptr->distance = 0;
  s_grow_ptr->sort_string = &null;

  /* Initialise pointer to the next s_grow */
  s_grow_ptr->next_s_grow_ptr = NULL;

  /* Return current pointers */
  *fst_s_grow_ptr = first_s_grow_ptr;
  *lst_s_grow_ptr = last_s_grow_ptr;

  /* Return new s_grow pointer */
  return(s_grow_ptr);
}
/***********************************************************************

s_get_grow_data  -  Read in the data from grow.out or the igrow.out.gz
                    file to get all its interactions data

***********************************************************************/

int s_get_grow_data(char *pdb_code,char *pdbsum_dir,
                    struct s_grow **fst_s_grow_ptr,int type)
{
  char file_name[S_FILENAME_LEN], input_line[S_LINELEN + 1];
  char chain1, chain2, number_string[20];

  int err_no, len, line, natoms, ngrow, nresid;

  struct s_grow *first_s_grow_ptr, *last_s_grow_ptr, *s_grow_ptr;

#ifdef TEST
  FILE *file_ptr;
#else
  gzFile file_ptr;
#endif

  /* Initialise variables */
  line = 0;
  natoms = 0;
  ngrow = 0;
  nresid = 0;

  /* Initialise pointers */
  *fst_s_grow_ptr = first_s_grow_ptr = NULL;
  last_s_grow_ptr = NULL;
  s_grow_ptr = NULL;

  /* Form the file name of the relevant grow.out file */
  strcpy(file_name,pdbsum_dir);
  if (type == S_IGROW)
    strcat(file_name,"igrow.out.gz");
  else
    strcat(file_name,"grow.out");

#ifdef TEST
  /* Form name of input file */
  if (type == S_IGROW)
    strcpy(file_name,"igrow.out");

  /* Open the file, if it exists */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("\n*** Warning. Unable to open file [%s]\n",file_name);
      return(0);
    }

  /* Loop while reading in records from input file */
  while (fgets(input_line,S_LINELEN,file_ptr)!= NULL)
#else
  /* Open the file, if it exists */
  file_ptr = gzopen(file_name,"rb");
  if (file_ptr ==  NULL)
    {
      /* Check for insufficiency of memory */
      if (err_no)
        {
          printf("*** ERROR. Unable to open file [%s]\n",
                 file_name);
          return(ngrow);
        }

      /* Otherwise, error must be due to missing file */
      return(ngrow);
    }

  /* Loop while reading in records from the file */
  while (gzgets(file_ptr,input_line,S_LINELEN)!= NULL)
#endif
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,S_LINELEN);

      /* If this is not a comment line, process its data */
      if (input_line[0] != '#' && len > 65)
        {
          /* Get the two chains */
          chain1 = input_line[6];
          chain2 = input_line[54];

          /* Create a s_grow record */
          s_grow_ptr
            = create_s_grow_record(&first_s_grow_ptr,&last_s_grow_ptr);

          /* Increment count of s_grow records stored */
          ngrow++;

          /* Get the residue name, number and chain of first atom
             and residue */
          strncpy(s_grow_ptr->atom_name[0],input_line+24,4);
          s_grow_ptr->atom_name[0][4] = '\0';
          strncpy(s_grow_ptr->res_name[0],input_line+14,3);
          s_grow_ptr->res_name[0][3] = '\0';
          strncpy(s_grow_ptr->res_num[0],input_line+9,5);
          s_grow_ptr->res_num[0][5] = '\0';
          s_grow_ptr->chain[0] = input_line[6];

          /* Repeat for second atom */
          strncpy(s_grow_ptr->atom_name[1],input_line+33,4);
          s_grow_ptr->atom_name[1][4] = '\0';
          strncpy(s_grow_ptr->res_name[1],input_line+44,3);
          s_grow_ptr->res_name[1][3] = '\0';
          strncpy(s_grow_ptr->res_num[1],input_line+49,5);
          s_grow_ptr->res_num[1][5] = '\0';
          s_grow_ptr->chain[1] = input_line[54];

          /* Get the interaction type */
          if (input_line[0] == 'H')
            s_grow_ptr->type = S_HBONDS;
          else if (input_line[0] == 'N')
            s_grow_ptr->type = S_CONTACTS;
          else if (input_line[0] == 'S')
            s_grow_ptr->type = S_SALT_BRIDGES;
          else if (input_line[0] == 'D')
            s_grow_ptr->type = S_DISULPHIDES;
          else if (input_line[0] == 'M')
            s_grow_ptr->type = S_CONTACTS;

          /* Save the residue types */
          s_grow_ptr->res_type[0] = input_line[2];
          s_grow_ptr->res_type[1] = input_line[4];

          /* Get the interaction distance */
          strncpy(number_string,input_line+56,5);
          number_string[5] = '\0';
          s_grow_ptr->distance = atof(number_string);
        }
    }

  /* Return s_grow pointer */
  *fst_s_grow_ptr = first_s_grow_ptr;

  /* Close file */
#ifdef TEST
  fclose(file_ptr);
#else
  gzclose(file_ptr);
#endif

  /* Return number of records created */
  return(ngrow);
}
/***********************************************************************

s_free_grow_data  -  Free up memory used by grow data

***********************************************************************/

void s_free_grow_data(struct s_grow *first_s_grow_ptr)
{
  struct s_grow *s_grow_ptr, *next_s_grow_ptr;

  /* Get first s_grow */
  s_grow_ptr = first_s_grow_ptr;

  /* Loop through the s_grow records until done */
  while (s_grow_ptr != NULL)
    {
      /* Get pointer to the next s_grow record */
      next_s_grow_ptr = s_grow_ptr->next_s_grow_ptr;

      /* Free up memory taken by current s_grow */
      free(s_grow_ptr);

      /* Go to the next s_grow */
      s_grow_ptr = next_s_grow_ptr;
    }
}
/***********************************************************************

s_grow_sort  -  Shell-sort which sorts the grow records in order of
                chain and grow number

***********************************************************************/

void s_grow_sort(struct s_grow **s_grow_index_ptr,int nindex)
{
  int endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct s_grow *grow1_ptr, *grow2_ptr, *swap_grow_ptr;

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

              /* Get the two grow records */
              grow1_ptr = s_grow_index_ptr[indi];
              grow2_ptr = s_grow_index_ptr[indl];

              /* If in wrong order, the swap sort-pointers */
              if (strcmp(grow1_ptr->sort_string,
                         grow2_ptr->sort_string) > 0)
                {
                  swap_grow_ptr = s_grow_index_ptr[indi];
                  s_grow_index_ptr[indi] = s_grow_index_ptr[indl];
                  s_grow_index_ptr[indl] = swap_grow_ptr;

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

s_sort_grow_records  -  Re-order the grow records in order of
                        chain, residue number and interaction type

***********************************************************************/

void s_sort_grow_records(struct s_grow **fst_s_grow_ptr,int ngrow)
{
  int igrow;

  struct s_grow *s_grow_ptr, *first_s_grow_ptr, *prev_s_grow_ptr;

  struct s_grow **s_grow_index_ptr;

  /* Initialise variables */
  first_s_grow_ptr = *fst_s_grow_ptr;
  
  /* Create arrays of grow pointers and sort index */
  s_grow_index_ptr
    = (struct s_grow **) malloc(ngrow * sizeof(struct s_grow *));
  if (s_grow_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for grow-sort\n");
      exit(1);
    }
  
  /* Get the first grow */
  s_grow_ptr = first_s_grow_ptr;
  igrow = 0;

  /* Loop through the grows to place their pointers in the array */
  while (s_grow_ptr != NULL && igrow < ngrow)
    {
      /* Place grow pointer on the array */
      s_grow_index_ptr[igrow] = s_grow_ptr;

      /* Get the next grow */
      s_grow_ptr = s_grow_ptr->next_s_grow_ptr;
      igrow++;
    }
  ngrow = igrow;

  /* Sort the grow pointers into order of grow number and name */
  if (ngrow > 1)
    s_grow_sort(s_grow_index_ptr,ngrow);

  /* Reinitialise pointers */
  first_s_grow_ptr = prev_s_grow_ptr = NULL;

  /* After sort, re-order the grow pointers */
  for (igrow = 0; igrow < ngrow; igrow++)
    {
      /* Get this grow record */
      s_grow_ptr = s_grow_index_ptr[igrow];

      /* Set next pointer to NULL */
      s_grow_ptr->next_s_grow_ptr = NULL;

      /* If this is the first record, save its pointer */
      if (first_s_grow_ptr == NULL)
        first_s_grow_ptr = s_grow_ptr;
      else
        prev_s_grow_ptr->next_s_grow_ptr = s_grow_ptr;

      /* Save the current pointer */
      prev_s_grow_ptr = s_grow_ptr;
    }

  /* Return the new first grow pointer */
  *fst_s_grow_ptr = first_s_grow_ptr;
}
/***********************************************************************

s_init_sort  -  Initialise all the grow records' sort string to prepare
                for sorting

***********************************************************************/

void s_init_sort(struct s_grow *first_s_grow_ptr)
{
  char sort_string[12], type[2];

  struct s_grow *s_grow_ptr;

  /* Get the first grow record */
  s_grow_ptr = first_s_grow_ptr;

  /* Loop through the grows to place their pointers in the array */
  while (s_grow_ptr != NULL)
    {
      /* Form full sort string from chain, residue number and
         interaction type */
      sort_string[0] = s_grow_ptr->chain[0];
      sort_string[1] = '\0';
      strcat(sort_string,s_grow_ptr->res_num[0]);
      strcat(sort_string,s_grow_ptr->res_name[0]);
      sprintf(type,"%d",s_grow_ptr->type);
      strcat(sort_string,type);

      /* Store the sort string */
      store_string(&(s_grow_ptr->sort_string),sort_string);

      /* Get the next grow record */
      s_grow_ptr = s_grow_ptr->next_s_grow_ptr;
    }
}
/***********************************************************************

s_get_protein_name  -  Pick up the protein's name from the PDBsum 
                       database.dat file

***********************************************************************/

int s_get_protein_name(char *pdb_code,char *pdbsum_dir,
                       char *protein_name)
{
  char file_name[S_FILENAME_LEN], input_line[S_LINELEN + 1];
  char ch;

  char *directory, *string_ptr;

  int field_len[S_NFIELDS], i, len, line, type;
  int done, have_name, obsolete, wanted;

  FILE *file_ptr;

  /* Initialise variables */
  have_name = FALSE;
  protein_name[0] = '\0';
  line = 0;
  for (i = 0; i < S_NFIELDS; i++)
    field_len[i] = strlen(S_FIELD_NAME[i]);
  obsolete = FALSE;

  /* Form the file name of the database.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"database.dat");

  /* Open the PDBsum database.dat file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the file */
      while (fgets(input_line,S_LINELEN,file_ptr) != NULL)
        {
          /* Initialise field type */
          type = S_UNWANTED;

          /* Check whether this is one of the wanted fields */
          for (i = 0; i < S_NFIELDS && type == S_UNWANTED; i++)
            if (!strncmp(input_line,S_FIELD_NAME[i],field_len[i]))
              type = i;

          /* If field wanted, then process */
          if (type != S_UNWANTED)
            {
              /* Truncate the string */
              string_truncate(input_line,S_LINELEN);

              /* Write out according to field */
              switch(type)
                {
                  /* Protein name */
                case S_TITLE:
                case S_STRUCTURE:

                  /* Write out if we don't already have this protein's
                     name */
                  if (have_name == FALSE)
                    {
                      /* Store the protein name */
                      strcpy(protein_name,input_line+field_len[type]);

                      /* Set flag */
                      have_name = TRUE;
                    }
                  break;

                  /* Superseded entry */
                case S_OBSOLETE:

                  /* Set flag */
                  obsolete = TRUE;
                  break;

                  /* Ignore other fields */
                case S_DATE:
                case S_RELEASE_DATE:
                case S_XRAY:
                case S_RESOLUTION:
                case S_NMR:
                case S_NMODELS:
                case S_TH_MODEL:

                default :
                  break;
                }
            }
        }

      /* Close the file */
      fclose(file_ptr);
    }

  else
    printf("*** Warning. Unable to open %s\n",file_name);

  /* Return the obsolete flag */
  return(obsolete);
}
/***********************************************************************

s_get_equiv_chain  -  Read the database.dat file to see if current chain
                      is a copy of another

***********************************************************************/

char s_get_equiv_chain(char chain,char *pdbsum_dir)
{
  char unique_chain;
  char input_line[S_LINELEN + 1], file_name[S_FILENAME_LEN];

  char *string_ptr;

  int done;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  unique_chain = '\0';

  /* Form name of the database.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"database.dat");

  /* Open the database.dat file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find database.dat file, then return */
      return(chain);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,S_LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      string_truncate(input_line,S_LINELEN);

      /* Check for chain record */
      if (!strncmp(input_line,"CHAIN[",6))
        {
          /* Get the unique chain*/
          string_ptr = strchr(input_line,' ');
          if (string_ptr)
            unique_chain = string_ptr[1];
        }

      /* Check for copy record */
      else if (!strncmp(input_line,"COPY[",5))
        {
          /* If this is our chain, then return the representative chain */
          string_ptr = strchr(input_line,' ');
          if (string_ptr != NULL &&
              string_ptr[1] == chain &&
              unique_chain != '\0')
            {
              chain = unique_chain;
              done = TRUE;
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the equivalent chain */
  return(chain);
}
/***********************************************************************

s_get_resdata  -  Read in the residue numbers from the resdata.dat file

***********************************************************************/

int s_get_resdata(char *pdb_code,char chain,int pdb_type,
                  struct r_residue **fst_r_residue_ptr,
                  struct c_file_data file_name_ptr)
{
  char input_line[S_LINELEN + 1], file_name[S_FILENAME_LEN];
  char pdbsum_dir[S_FILENAME_LEN];
  char number_string[20], tmp_line[S_LINELEN + 1];
  char aa_code, athet, equiv_chain, inchain, res_name[4], res_num[6];
  char sst;

  char *string_ptr;

  static char subscript[] = { 'D', 'L', 'M' };

  int dir_type, i, ipos, len, nresid, number, type;
  int domain, lig_hb, lig_nb;
  int done, got_chain, reading_sequence, wanted;

  struct r_residue *r_residue_ptr, *first_r_residue_ptr, *last_r_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  athet = 'A';
  equiv_chain = chain;
  nresid = 0;
  done = FALSE;
  got_chain = FALSE;
  lig_hb = lig_nb = 0;
  wanted = FALSE;
  reading_sequence = FALSE;

  /* Initialise pointers */
  r_residue_ptr = first_r_residue_ptr = last_r_residue_ptr = NULL;
  *fst_r_residue_ptr = NULL;

  /* Find the PDBsum directory for this entry */
  dir_type = find_pdbsum_dir(pdb_code,file_name_ptr.pdbsum_template,
                             pdb_type,pdbsum_dir);

  /* Get the equivalent chain if different from the current one */
  if (chain != '\0')
    equiv_chain = s_get_equiv_chain(chain,pdbsum_dir);

  /* Form name of the resdata.dat file */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"resdata.dat");

  /* Open the resdata.dat file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find resdata.dat file, then return empty-handed */
      return(0);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,S_LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      // string_truncate(input_line,S_LINELEN);
      string_chomp(input_line);

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

              /* Reinitialise counts */
              lig_hb = lig_nb = 0;

              /* Get the domain number */
              strcpy(number_string,input_line+16);
              domain = atoi(number_string);

              /* If chain is wanted, process the record */
              if (inchain == equiv_chain || equiv_chain == '\0')
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

                  /* Create a residue record and store data just read in */
                  r_residue_ptr
                    = create_r_residue_record(&first_r_residue_ptr,
                                              &last_r_residue_ptr,NULL,
                                              athet,NULL,res_name,res_num,
                                              NULL);
                  r_residue_ptr->chain = inchain;

                  /* Get the corresponding secondary structure */
                  sst = input_line[14];

                  /* Save the secondary structure */
                  if (sst == 'H')
                    r_residue_ptr->sec_struc = R_HELIX;
                  else if (sst == 'E')
                    r_residue_ptr->sec_struc = R_STRAND;
                  else
                    r_residue_ptr->sec_struc = R_COIL;

                  /* Increment count of residues stored */
                  nresid++;
                  r_residue_ptr->seq_no = nresid;
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
                  r_residue_ptr->disulphide_id = atoi(number_string);
                }
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Save pointer to first residue and number of residues stored */
  *fst_r_residue_ptr = first_r_residue_ptr;

  /* Return number of residues read in */
  return(nresid);
}
/***********************************************************************

s_transfer_resdata  -  Transfer the secondary structure assignments read
                       in from the resdata.dat file to the corresponding
                       residues read in from the PDB file

***********************************************************************/

void s_transfer_resdata(struct r_residue *first_pdb_r_residue_ptr,
                        struct r_residue *first_r_residue_ptr)
{
  struct r_residue *r_residue_ptr, pdb_r_residue_ptr;
  struct r_residue *next_r_residue_ptr, *found_r_residue_ptr;

  /* Get pointer to the first residue from resdata.dat */
  r_residue_ptr = first_r_residue_ptr;

  /* Loop through all the resdata.dat residues  */
  while (r_residue_ptr != NULL)
    {
      /* Get pointer to the next residue record */
      next_r_residue_ptr = r_residue_ptr->next_r_residue_ptr;

      /* Find this residue in those read from the PDB */
      found_r_residue_ptr
        = r_find_residue(first_pdb_r_residue_ptr,
                         r_residue_ptr->res_name,
                         r_residue_ptr->res_num,r_residue_ptr->chain);

      /* If found, transfer the data across */
      if (found_r_residue_ptr != NULL)
        {
          /* Transfer the secondary structure assignment */
          found_r_residue_ptr->disulphide_id
            = r_residue_ptr->disulphide_id;
          found_r_residue_ptr->domain = r_residue_ptr->domain;
          found_r_residue_ptr->sec_struc = r_residue_ptr->sec_struc;
          found_r_residue_ptr->seq_no = r_residue_ptr->seq_no;
        }

      /* Free up memory taken by current residue */
      free(r_residue_ptr);

      /* Go to the next residue */
      r_residue_ptr = next_r_residue_ptr;
    }
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int s_main_test(int argc,char *argv[])
//int main(int argc,char *argv[])
{
  char pdb_code[PDB_CODELEN + 1], pdbdir_template[S_FILENAME_LEN];
  char pdbsumdir_template[S_FILENAME_LEN];

  int ngrow, ntokens, type;

  struct s_grow *first_s_grow_ptr;

  /* Initialise variables */
  ntokens = argc - 1;
  type = S_GROW;

  /* Initialise pointers */
  first_s_grow_ptr = NULL;

  /* Get the PDB code */
  if ((ntokens == 0 || ntokens > 2) ||
      (ntokens > 0 && strlen(argv[1]) > PDB_CODELEN))
    {
      printf("\n");
      printf("Correct usage is ...\n");
      printf("\n");
      printf("    readsum  pdb_code  [-i]\n");
      printf("\n");
      printf("where -i  indicates that igrow.out.gz file to be read\n");
      printf("\n");
      exit(1);
    }

  /* Get the PDB code from the first token */
  strcpy(pdb_code,argv[1]);
  printf("PDB code: %s\n",pdb_code);

  /* If have a second argument, check whether it is the -i flag */
  if (ntokens > 1 && !strcmp(argv[2],"-i"))
    type = S_IGROW;

  /* Use PDB location pattern to form full file-name */
  strcpy(pdbdir_template,"c:\\roman\\pdbsum\\pdb\\pdb[abcd].ent");
  strcpy(pdbsumdir_template,"c:\\roman\\pdbsum\\pdbsum\\[abcd]");

  /* Read in the data from the grow.out or igrow.out.gz file */
  ngrow
    = s_get_grow_data(pdb_code,pdbsumdir_template,&first_s_grow_ptr,
                      type);

  printf("Number of grow records created:   %8d\n",ngrow);

  exit(0);
}
