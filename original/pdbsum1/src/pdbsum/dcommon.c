/***********************************************************************

  dcommon.c - Common routines used by multiple DisaStr programs

***********************************************************************

Date:-         17 Nov 2017

Written by:-   Roman Laskowski


----------------------------------------------------------------------*/

/* Compilation
   -----------

cc -c dcommon.c

*/

/* Routines
   --------
   -> check_isoforms
         -> get_tokens
         -> idx_search
         -> read_in_sequence
         -> map_isoform
               -> isoform_search
         -> get_mapped_aa
         -> free_tokens (in common.c)
   -> compile_fasta_index
         -> create_idx_record
   -> compute_cadd_mark
   -> free_cadd_memory
   -> free_gnomad_memory
   -> get_isoforms
         -> string_chomp (in common.c)
         -> get_tokens (in common.c)
         -> store_string (in common.c)
         -> free_tokens (in common.c)
         -> append_string (in common.c)
         -> process_isoform
               -> store_string (in common.c)
               -> create_isoform_entry
         -> sort_isoforms
               -> isoform_sort
   -> read_isoforms_txt
         -> string_chomp (in common.c)
         -> get_tokens (in common.c)
         -> store_string (in common.c)
         -> free_tokens (in common.c)
         -> sort_isoforms
               -> isoform_sort
   -> get_label_letters
   -> get_tables
         -> read_in_table
               -> get_tokens (in common.c)
               -> free_tokens (in common.c)
   -> plot_label
   -> read_cadd_txt
         -> get_tokens (in common.c)
         -> create_cadd_entry
               -> store_string (in common.c)
         -> free_tokens (in common.c)
   -> read_gnomad_idx
         -> tokenize (in common.c)
         -> create_idx_record
   -> read_idx_file
         -> tokenize (in common.c)
         -> create_idx_record
   -> index_files
         -> index_file
               -> get_tokens (in common.c)
               -> free_tokens (in common.c)
         -> create_idx_index
   -> index_file
         -> tokenize (in common.c)
         -> create_idx_record
   -> get_colour_restype
   -> get_gnomad_offset
         -> tokenize (in common.c)
   -> read_gnomad_dat
         -> tokenize (in common.c)
         -> create_gnomad_record
               -> store_string (in common.c)
   -> create_array
   -> create_idx_record
   -> create_idx_index
         -> idx_sort
   -> exp3html
   -> idx_search
   -> calc_variant_scores
   -> get_data_dir
   -> get_domain_dir
   -> get_residue_name
         -> get_res_name (in common.c)

*/


#include "common.h"
#include "dcommon.h"
#include "readpdb.h"

/* Prototypes
   ---------- */

/* In common.c */
void append_string(char **string_ptr,char *string);
int apply_pixel(int bg_colour,float grey_shade,int *text_rgb);
void free_tokens(char **token,int ntokens);
float get_pixel_shade(char ch);
int get_res_name(char aa_code,char *res_name);
int get_tokens(char *line,char **token,int max_token,char token_sep);
void store_string(char **string_ptr,char *string);
int string_chomp(char *string);
int string_chop(char *string);
int string_truncate(char *string,int max_length);
int tokenize(char *string,char **token,int max_tokens,char token_sep);

/***********************************************************************

isoform_search  -  Perform a binary search for the first isoform record
                   for the given UniProt code

***********************************************************************/

struct isoform *isoform_search(struct isoform **isoform_index_ptr,
                               int nisoform,char *uniprot_acc,
                               int isoform)
{
  int bottom_point, ipoint, new_point, top_point;
  int last_lower, last_higher;
  int found;

  struct isoform *isoform_ptr, *found_isoform_ptr, *prev_isoform_ptr;

  /* Initialise variables */
  found = FALSE;
  found_isoform_ptr = NULL;
  if (nisoform == 0)
    return(found_isoform_ptr);
  last_lower = last_higher = -1;

  /* Make starting point half-way into the index */
  ipoint = nisoform / 2;
  bottom_point = 0;
  top_point = nisoform - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated GO pointer */
      isoform_ptr = isoform_index_ptr[ipoint];

      /* Check if index value is equal to the value we want */
      if (!strcmp(uniprot_acc,isoform_ptr->uniprot_acc) &&
          ((isoform == isoform_ptr->isoform &&
           isoform_ptr->order == 0) ||
           (isoform == 0 && isoform_ptr->count == 0)))
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_isoform_ptr = isoform_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (strcmp(uniprot_acc,isoform_ptr->uniprot_acc) < 0 ||
               (!strcmp(uniprot_acc,isoform_ptr->uniprot_acc) &&
                (isoform < isoform_ptr->isoform ||
                 (isoform == isoform_ptr->isoform &&
                  isoform_ptr->order > 0))))
        {
          /* If this is the first record in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, need to search earlier */
          else
            {
              /* Save the current position as the uppermost position
                 to search on */
              top_point = ipoint;

              /* Set pointer to halfway between the last search
                 low-search and here */
              ipoint = (ipoint + bottom_point) / 2;

              /* If we've already tried this one, then end here */
              if (ipoint == last_lower)
                found = TRUE;

              /* Save this as the last lower step made */
              last_lower = ipoint;
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
          if (ipoint >= nisoform - 1)
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

              /* If we've already tried this one, then end here */
              if (ipoint == last_higher)
                found = TRUE;

              /* Save this as the last higher step made */
              last_higher = ipoint;
            }
        }
    }

  /* Return the pointer */
  return(found_isoform_ptr);
}
/***********************************************************************

read_in_sequence  -  Read in the canonical protein sequence from the
                     FASTA file

***********************************************************************/

int read_in_sequence(struct idx *idx_ptr,char *seq_uniprot_acc,
                     char *sequence,char *protein_name,FILE *seqfile_ptr,
                     int gzipped
#ifdef TEST
#else
                   ,gzFile gzseqfile_ptr
#endif
                   )
{
  char input_line[D_LINELEN + 1];

  char *string_ptr;

  int le, line;
  int done, eof, error, ok, wanted;

  /* Initialise variables */
  eof = FALSE;
  ok = TRUE;
  protein_name[0] = '\0';
  sequence[0] = '\0';

  /* Find start-point in data file usinf fseek */
#ifdef TEST
#else
  if (gzipped == TRUE)
    error = gzseek(gzseqfile_ptr,idx_ptr->offset,SEEK_SET);
  else
#endif
    error = fseek(seqfile_ptr,idx_ptr->offset,SEEK_SET);
  
  /* If OK, then read in the data from this point onwards */
  if (error != TRUE)
    {
      /* Initialise variables */
      done = FALSE;
      line = 0;
      wanted = FALSE;

      /* Loop while reading in records from the input file */
      while (eof == FALSE && done == FALSE)
        {
          /* Read in the next record */
#ifdef TEST
#else
          if (gzipped == TRUE)
            {
              if (gzgets(gzseqfile_ptr,input_line,D_LINELEN) == NULL)
                eof = TRUE;
            }
          else
#endif
            if (fgets(input_line,D_LINELEN,seqfile_ptr) == NULL)
              eof = TRUE;

          /* If record read in OK, then process it */
          if (eof == FALSE)
            {
              /* If this is the key line, check UniProt code */
              if (!strncmp(input_line,">sp|",4))
                {
                  /* Check that it is the correct UniProt code */
                  if (!strncmp(input_line + 4,idx_ptr->key,
                               strlen(idx_ptr->key)))
                    {
                      /* Extract the protein name */
                      string_ptr = strchr(input_line,' ');
                      if (string_ptr != NULL)
                        {
                          strcpy(protein_name,string_ptr + 1);

                          /* Remove extraneous data */
                          string_ptr = strstr(protein_name," OS=");
                          if (string_ptr != NULL)
                            string_ptr[0] = '\0';
                        }
                      else
                        strcpy(protein_name,input_line + 4);
                          
                      /* Set flag that wanted */
                      wanted = TRUE;
                    }

                  /* Otherwise, see if we are done */
                  else if (wanted == TRUE)
                    {
                      done = TRUE;
                      wanted = FALSE;
                    }
                }

              /* Otherwise, if line is part of the sequence, append it */
              else if (wanted == TRUE)
                {
                  string_truncate(input_line,D_MAX_SEQLEN);
                  strcat(sequence,input_line);
                }
            }
        }
    }
  
  /* If we have found a sequence, save the UniProt code */
  if (sequence[0] != '\0')
    strcpy(seq_uniprot_acc,idx_ptr->key);

  /* Return that sequence read in OK */
  return(ok);
}
/***********************************************************************

map_isoform  -  Map the amino acids in the current isoform to the
                corresponding amino acids in the canonical sequence

***********************************************************************/

int map_isoform(char *uniprot_acc,int isoform,int seq_len,
                struct isoform **isoform_index_ptr,int nisoform,
                int *seqmap)
{
  static char last_uniprot_acc[UNIPROT_LEN];

  int ipos, nomap_len, nout, pos_start, res_start, seq_pos, shift;
  int done;

  static int first = TRUE;
  static int last_isoform = 0;
  static int map_len = 0;

  struct isoform *isoform_ptr, *first_isoform_ptr;
  static struct isoform *prev_isoform_ptr;

  /* Initialise variables */
  done = FALSE;
  seq_pos = 0;
  shift = 0;

  /* If first pass, initialise previous isoform details */
  if (first == TRUE)
    {
      /* Initialise details */
      last_uniprot_acc[0] = '\0';
      last_isoform = 0;
      prev_isoform_ptr = NULL;
      first = FALSE;
    }

  /* If this is the same isoform as before, then we already have the
     mapping, so can return */
  if (isoform == last_isoform && !strcmp(uniprot_acc,last_uniprot_acc))
    return(map_len);

  /* Otherwise, perform the search */
  else
    {
      /* Find the first isoform record for this UniProt code */
      first_isoform_ptr
        = isoform_search(isoform_index_ptr,nisoform,uniprot_acc,
                         isoform);

      /* Save the details in case required next time */
      prev_isoform_ptr = first_isoform_ptr;
      strcpy(last_uniprot_acc,uniprot_acc);
      last_isoform = isoform;
    }

  /* Initialise the map length */
  map_len = seq_len;
  done = FALSE;

  /* Get the first isoform record */
  isoform_ptr = first_isoform_ptr;

  /* Loop over the isoform records to compute mapping of current
     position */
  while (isoform_ptr != NULL && done == FALSE)
    {
      /* If this is our isoform, add the no-map length to the sequence
         length */
      if (isoform_ptr->isoform == isoform &&
          !strcmp(isoform_ptr->uniprot_acc,uniprot_acc))
        map_len = map_len + isoform_ptr->nomap_len;
        
      /* Otherwise, we are done */
      else
        done = TRUE;

      /* Get the next record */
      isoform_ptr = isoform_ptr->next_isoform_ptr;
    }

  /* Initialise the mapping */
  for (ipos = 0; ipos < map_len; ipos++)
    seqmap[ipos] = ipos + 1;

  /* If isoform not found, then can't map sequence positions to the
     canonical sequence */
  if (first_isoform_ptr == NULL)
    return(map_len);

  /* Get the first isoform record */
  isoform_ptr = first_isoform_ptr;
  done = FALSE;

  /* Loop over the isoform records to compute mapping of current
     position */
  while (isoform_ptr != NULL && done == FALSE)
    {
      /* If this is our isoform, then update the sequence mappings */
      if (isoform_ptr->isoform == isoform &&
          !strcmp(isoform_ptr->uniprot_acc,uniprot_acc))
        {
          /* Get the parameters */
          res_start = isoform_ptr->res_start;
          nomap_len = isoform_ptr->nomap_len;
          shift = isoform_ptr->shift;

          /* Get the mapping position of the residue start */
          pos_start = -1;
          for (ipos = 0; ipos < map_len && pos_start == -1; ipos++)
            {
              /* If this is the starting residue, save its position */
              if (seqmap[ipos] == res_start)
                pos_start = ipos;
            }

          /* If this isoform contains a sequence that can't be mapped onto
             the canonical sequence, set these residues mappings to -1 */
          if (nomap_len > 0 && pos_start > -1)
            {
              /* Initialise count of positions knocked out */
              nout = 0;
              
              /* Loop through the sequence to knock out the no-map
                 positions */
              for (ipos = pos_start; nout < nomap_len && ipos < map_len;
                   ipos++)
                {
                  /* Set to -1 and increment count */
                  seqmap[ipos] = -1;
                  nout++;
                }
            }

          /* If we have a shift, apply that */
          if (shift != 0)
            {
              /* Loop to apply shift to all subsequent residues */
              for (ipos = pos_start; ipos < map_len; ipos++)
                {
                  /* If a valid mapping, apply shift */
                  if (seqmap[ipos] > -1)
                    {
                      /* Apply shift */
                      seqmap[ipos] = seqmap[ipos] + shift;

                      /* If shift takes mapping off either end of the
                         sequence, set to -1 */
                      if (seqmap[ipos] < 1 || seqmap[ipos] > map_len)
                        seqmap[ipos] = -1;
                    }
                }
            }
        }

      /* Otherwise, we are done */
      else
        done = TRUE;

      /* Get the next record */
      isoform_ptr = isoform_ptr->next_isoform_ptr;
    }

  return(map_len);
}
/***********************************************************************

get_mapped_aa  -  Determine which position in the canonical sequence
                  the current residue maps onto

***********************************************************************/

int get_mapped_aa(char *sequence,int *seqmap,int map_len,int seq_len,
                  int seq_no,char aa_code)
{
  char aa_canonical;

  int ires;

  /* Get the mapping for this position */
  ires = -1;
  if (seq_no <= map_len)
    ires = seqmap[seq_no - 1];

  /* Check that the amino acid at the mapped position is what we're
     expecting */
  if (ires > 0 && ires <= seq_len)
    {
      /* Get the corresponding amino acid in the canonical sequence */
      aa_canonical = sequence[ires - 1];

      /* If this doesn't match what we're expecting, can't use this
         variant */
      if (aa_canonical != aa_code)
        ires = -1;
    }
  else
    ires = -1;

  /* Return the mapped sequence position */
  return(ires);
}
/***********************************************************************

idx_search  -  Perform a binary search to find the given key

***********************************************************************/

struct idx *idx_search(struct idx **idx_index_ptr,long nidx,char *key)
{
  int order;
  int found;

  long bottom_point, ipoint, jpoint, new_point, top_point;

  struct idx *idx_ptr, *found_idx_ptr, *prev_idx_ptr, *try_idx_ptr;

  /* Initialise variables */
  found = FALSE;
  found_idx_ptr = NULL;
  if (nidx == 0)
    return(found_idx_ptr);

  /* Make starting point half-way into the index */
  ipoint = nidx / 2;
  bottom_point = 0;
  top_point = nidx - 1;

  /* Loop until required index position found */
  while (found == FALSE)
    {
      /* Get the associated pointer */
      idx_ptr = idx_index_ptr[ipoint];

      /* Determine which order the records are in */
      if (!strcmp(key,idx_ptr->key) && idx_ptr->nentry == 0)
        order = D_EQUAL;
      else if (strcmp(key,idx_ptr->key) < 0 ||
               (!strcmp(key,idx_ptr->key) && idx_ptr->nentry > 0))
        order = D_LOWER;
      else
        order = D_HIGHER;

      /* Check if index value is equal to the value we want */
      if (order == D_EQUAL)
        {
          /* Set flag to indicate that match found and store pointer */
          found = TRUE;
          found_idx_ptr = idx_ptr;
        }

      /* If index value is higher than the one we want, then need
         to search lower */
      else if (order == D_LOWER)
        {
          /* If this is the first go in the sorted list, then can't
             go any lower */
          if (ipoint == 0)
            {
              /* Set flag that we are done */
              found = TRUE;
            }

          /* Otherwise, check whether the previous entry's key
             is below the one we're searching for */
          else
            {
              /* Get the pointer to the previous entry */
              prev_idx_ptr = idx_index_ptr[ipoint - 1];

              /* If value is below the search value, we've found the
                 position in the index we're after, but don't have an exact
                 match*/
              if (strcmp(prev_idx_ptr->key,key) < 0)
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
          if (ipoint >= nidx - 1)
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
  return(found_idx_ptr);
}
/***********************************************************************

check_isoforms  -  Perform a binary search for the first isoform record
                   for the given UniProt code

***********************************************************************/

int check_isoforms(struct isoform **isoform_index_ptr,int nisoform,
                   struct idx **key_idx_ptr,int nkeys,char **token,
                   char *uniprot_acc,char *seq_uniprot_acc,
                   char *sequence,char *protein_name,int *seqmap,
                   int *lst_isoform,FILE *seqfile_ptr,int gzipped
#ifdef TEST
#else
                   ,gzFile gzseqfile_ptr
#endif
                   )
{
  char aa[2], aa_change[4], aa_code, aa_mut, number_string[D_LINELEN];
  char codon_change[10], enst_id[20], note[100], strand_dir[100];

  static char prev_uniprot_acc[UNIPROT_LEN];

  char *string_ptr, *subtoken[D_MAX_TOKEN];

  int first_token, ipos, isoform, itoken, len, ntokens, seq_len, seq_no;
  int last_isoform, map_len, seq_pos;
  int done, start;

  static int first = TRUE;

  struct idx *idx_ptr;
  struct isoform *isoform_ptr, *found_isoform_ptr, *prev_isoform_ptr;

  static struct idx *prev_idx_ptr;

  /* Initialise variables */
  strcpy(aa_change,"-");
  strcpy(codon_change,"-");
  strcpy(strand_dir,"-");
  done = FALSE;
  first_token = 0;
  ipos = 0;
  isoform = -1;
  last_isoform = *lst_isoform;
  map_len = 0;
  seq_no = 0;
  uniprot_acc[0] = '\0';

  /* If the first pass, initialise previous UniProt entry */
  if (first == TRUE)
    {
      /* Initialise */
      prev_uniprot_acc[0] = '\0';
      prev_idx_ptr = NULL;
      first = FALSE;
    }

  /* Get the tokens from the all-transcripts field */
  ntokens = get_tokens(token[ALL_TRANSCRIPTS],subtoken,D_MAX_TOKEN,',');

  /* Loop over the tokens to find any alternative isoforms */
  for (itoken = 0; itoken < ntokens && done == FALSE; itoken++)
    {
     /* Initialise flag */
      start = FALSE;

      /* If this is the last token, then extract text */
      if (itoken == ntokens - 1)
        start = TRUE;

      /* If token has a slash, then at the start of the next ENST */
      else if ((string_ptr = strstr(subtoken[itoken]," / ")) != NULL)
        start = TRUE;
      
      /* If we have the start, retrieve the UniProt accession and 
         sequence position */
      if (start == TRUE)
        {
          /* Get the ENST id */
          string_ptr = strstr(subtoken[first_token]," / ");
          if (string_ptr != NULL)
            strcpy(enst_id,string_ptr + 3);
          else
            strcpy(enst_id,subtoken[first_token]);

          /* Get the UniProt isoform */
          strcpy(uniprot_acc,subtoken[first_token + SUB_UNIPROT_ACC]);

          /* Get the sequence position */
          if (subtoken[first_token + SUB_SEQ_NO][0] != '-')
            {
              /* Get the sequence position */
              seq_no = atoi(subtoken[first_token + SUB_SEQ_NO]);
            }
          else
            seq_no = 0;

          /* If we have a valid isoform and sequence number, extract
             the details */
          if ((string_ptr = strchr(uniprot_acc,'-')) != NULL &&
              seq_no > 0)
            {
              /* Extract the isoform number */
              isoform = atoi(string_ptr + 1);
              string_ptr[0] = '\0';

              /* Get the amino acid change */
              if (strlen(subtoken[first_token + SUB_AA_CHANGE]) == 3)
                strcpy(aa_change,subtoken[first_token + SUB_AA_CHANGE]);
              else if (strlen(subtoken[first_token + SUB_AA_CHANGE]) == 1)
                sprintf(aa_change,"%s/%s",
                        subtoken[first_token + SUB_AA_CHANGE],
                        subtoken[first_token + SUB_AA_CHANGE]);
              else
                strcpy(aa_change,"-");

              /* Save the codon change and strand direction */
              if (ntokens >= SUB_CODON_CHANGE)
                {
                  strcpy(codon_change,
                         subtoken[first_token + SUB_CODON_CHANGE]);
                  strcpy(strand_dir,
                         subtoken[first_token + SUB_STRAND_DIR]);

                  /* Truncate strand direction, if required */
                  string_ptr = strchr(strand_dir,' ');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                }

              /* Check we have a valid a.a. change */
              if (!strcmp(aa_change,"-") || aa_change[0] == '*')
                isoform = 0;

              /* If we have an isoform, check whether this UniProt code
                 is in our list of isoforms */
              if (isoform > 0)
                {
                  /* Get the original and variant residues */
                  aa_code = aa_change[0];
                  aa_mut = aa_change[2];

                  /* If the same UniProt code as last time, retrieve
                     the index */
                  if (!strcmp(uniprot_acc,prev_uniprot_acc))
                    idx_ptr = prev_idx_ptr;

                  /* Otherwise, perform the search */
                  else
                    {
                      /* Check if the UniProt code corresponds to a canonical
                         sequence */
                      idx_ptr
                        = idx_search(key_idx_ptr,nkeys,uniprot_acc);

                      /* If found, then store the details */
                      if (idx_ptr != NULL)
                        {
                          strcpy(prev_uniprot_acc,uniprot_acc);
                          prev_idx_ptr = idx_ptr;
                        }
                    }

                  /* If found, then retrieve the sequence */
                  if (idx_ptr != NULL)
                    {
                      /* If we don't yet have this protein's sequence,
                         read it in from the FASTA file */
                      if (sequence[0] == '\0' ||
                          strcmp(uniprot_acc,seq_uniprot_acc))
                        {
                          /* Read in the canonical protein sequence */
                          read_in_sequence(idx_ptr,seq_uniprot_acc,
                                           sequence,protein_name,seqfile_ptr,
                                           gzipped
#ifdef TEST
#else
                                           ,gzseqfile_ptr
#endif
                                           );

                          /* Re-initialise last isoform for which we have
                             a mapping */
                          last_isoform = 0;
                        }

                      /* Get the sequence length */
                      seq_len = strlen(sequence);

                      /* If we have a new isoform (or a new protein)
                         map the amino acids in the current isoform
                         to the corresponding positions in the
                         canonical sequence */
                      if (isoform != last_isoform)
                        {
                          /* Get the mapping */
                          map_len
                            = map_isoform(uniprot_acc,isoform,seq_len,
                                          isoform_index_ptr,nisoform,
                                          seqmap);

                          /* Save the isoform */
                          last_isoform = isoform;
                        }

                      /* Determine which position in the canonical sequence
                         the current residue maps onto */
                      if (map_len > 0)
                        seq_pos
                          = get_mapped_aa(sequence,seqmap,map_len,seq_len,
                                          seq_no,aa_code);
                      else
                        seq_pos = 0;

                      /* If this is a valid position, can use the mapping
                         to the canonical isoform */
                      if (seq_pos > 0)
                        {
                          /* Store the details for this isoform */
                          free(token[UNIPROT_ACCESSION]);
                          store_string(&(token[UNIPROT_ACCESSION]),
                                       uniprot_acc);
                          free(token[AA_CHANGE]);
                          store_string(&(token[AA_CHANGE]),aa_change);

                          /* Store the codon change */
                          if (strcmp(codon_change,"-"))
                            {
                              free(token[CODON_CHANGE]);
                              store_string(&(token[CODON_CHANGE]),
                                           codon_change);
                            }

                          /* Store the strand direction */
                          if (strcmp(strand_dir,"-"))
                            {
                              free(token[STRAND_DIR]);
                              store_string(&(token[STRAND_DIR]),strand_dir);

                              /* Update negative flag */
                              if (!strcmp(strand_dir,"positive") &&
                                  !strcmp(token[NEGATIVE],"TRUE"))
                                {
                                  free(token[NEGATIVE]);
                                  store_string(&(token[NEGATIVE]),"FALSE");
                                }

                              else if (!strcmp(strand_dir,"negative") &&
                                  !strcmp(token[NEGATIVE],"FALSE"))
                                {
                                  free(token[NEGATIVE]);
                                  store_string(&(token[NEGATIVE]),"TRUE");
                                }
                            }

                          /* Store the protein name */
                          free(token[PROTEIN_NAME]);
                          store_string(&(token[PROTEIN_NAME]),protein_name);

                          /* Store the amino acid code */
                          free(token[UNIPROT_AA]);
                          aa[0] = token[AA_CHANGE][0];
                          aa[1] = '\0';
                          store_string(&(token[UNIPROT_AA]),aa);

                          /* Store the mapped sequence number */
                          sprintf(number_string,"%d",seq_pos);
                          free(token[SEQ_NO]);
                          store_string(&(token[SEQ_NO]),number_string);

                          /* Store the gene name */
                          free(token[GENE]);
                          store_string(&(token[GENE]),
                                       subtoken[first_token + SUB_GENE]);

                          /* Change the transcript back */
                          free(token[TRANSCRIPT]);
                          store_string(&(token[TRANSCRIPT]),enst_id);

                          /* Append semi-colon */
                          append_string(&(token[TRANSCRIPT]),";");

                          /* Store the SIFTS and Polyphen scores */
                          free(token[SIFTS_SCORE]);
                          store_string(&(token[SIFTS_SCORE]),
                                       subtoken[first_token + SUB_SIFTS]);
                          free(token[POLYPHEN_SCORE]);
                          store_string(&(token[POLYPHEN_SCORE]),
                                       subtoken[first_token + SUB_POLYPHEN]);

                          /* Update the note */
                          sprintf(note,"Transcript mapped to Swissprot canonical "
                                  "sequence via isoform %s-%d",
                                  token[UNIPROT_ACCESSION],isoform);
                          free(token[NOTE]);
                          store_string(&(token[NOTE]),note);

                          /* Set flag that we're done */
                          done = TRUE;
                        }

                      /* If don't have a valid sequence position, set
                         isoform to zero */
                      else
                        isoform = 0;
                    }
                }
            }

          /* Set the first token to the current one */
          first_token = itoken;
        }
    }

  /* Free the tokens */
  free_tokens(subtoken,ntokens);

  /* Return last isoform stored */
  *lst_isoform = last_isoform;

  /* Return which isoform found, if any */
  return(isoform);
}
/***********************************************************************

create_idx_record  -  Create a new idx record

***********************************************************************/

struct idx *create_idx_record(struct idx **fst_idx_ptr,
                              struct idx **lst_idx_ptr,char *key,
                              long offset,int nentry)
{
  int itoken;

  struct idx *idx_ptr, *first_idx_ptr, *last_idx_ptr;

  /* Initialise idx pointers */
  idx_ptr = NULL;
  last_idx_ptr = *lst_idx_ptr;
  first_idx_ptr = *fst_idx_ptr;

  /* Create idx entry */
  idx_ptr = (struct idx*)malloc(sizeof(struct idx));
  if (idx_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct idx\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  store_string(&(idx_ptr->key),key);
  idx_ptr->offset = offset;
  idx_ptr->nentry = nentry;

  /* Initialise pointer to next idx record */
  idx_ptr->next_idx_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_idx_ptr == NULL)
    first_idx_ptr = idx_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_idx_ptr->next_idx_ptr = idx_ptr;
                      
  /* Save the current idx's pointer */
  last_idx_ptr = idx_ptr;

  /* Return current pointers */
  *fst_idx_ptr = first_idx_ptr;
  *lst_idx_ptr = last_idx_ptr;

  return(idx_ptr);
}
/***********************************************************************

compile_fasta_index  -  Read through the .fasta file to pick up the
                        offset positions of each key term

***********************************************************************/

long compile_fasta_index(char *file_name,struct idx **fst_idx_ptr,
                         char *key,char *current_uniprot_acc,int *gzip)
{
  char input_line[D_LINELEN + 1], uniprot_acc[UNIPROT_LEN];

  char *string_ptr;

  int ikey, len, lenkey;
  int eof, gzipped, wanted;

  long nkeys, nlines, offset;

  struct idx *idx_ptr, *first_idx_ptr, *last_idx_ptr;

  FILE *file_ptr;
#ifdef TEST
#else
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  eof = FALSE;
  *gzip = gzipped = FALSE;
  lenkey = strlen(key);
  nkeys = 0;
  nlines = 0;
  offset = 0;
  wanted = FALSE;

  /* Initialise pointers */
  idx_ptr = first_idx_ptr = last_idx_ptr = NULL;

  /* Determine if this is a gzipped file */
  if (strstr(file_name,".gz") != NULL)
    gzipped = TRUE;

#ifdef TEST
  if (gzipped == TRUE)
    {
      printf("*** ERROR. Cannot open gzipped file in TEST mode. Aborting\n");
      exit(1);
    }
#endif

  /* Open the file */
#ifdef TEST
#else
  if (gzipped == TRUE)
    gzfile_ptr = gzopen(file_name,"r");
  else
#endif
    file_ptr = fopen(file_name,"r");

  /* If file not found, then return */
  if ((gzipped == FALSE && file_ptr == NULL)
#ifdef TEST
#else
      || (gzipped == TRUE && gzfile_ptr == NULL)
#endif
      )   
    {
      printf("*** Warning. Unable to open file %s\n",file_name);
      return(nkeys);
    }

  /* Loop until we have finished reading the file */
  while (eof == FALSE)
    {
      /* Read in the next record */
#ifdef TEST
#else
      if (gzipped == TRUE)
        {
          if (gzgets(gzfile_ptr,input_line,D_LINELEN) == NULL)
            eof = TRUE;
        }
      else
#endif
        if (fgets(input_line,D_LINELEN,file_ptr) == NULL)
            eof = TRUE;

      /* If record read in OK, then process it */
      if (eof == FALSE)
        {
          /* If this is the key line, get the UniProt id */
          if (!strncmp(input_line,key,lenkey))
            {
              /* If this is a human sequence, store its offset */
              if (strstr(input_line,"_HUMAN ") != NULL)
                {
                  /* Extract the UniProt id */
                  string_ptr = strchr(input_line+lenkey,'|');
                  if (string_ptr != NULL)
                    {
                      /* Save the code */
                      string_ptr[0] = '\0';
                      strcpy(uniprot_acc,input_line+lenkey);

                      /* If entry is wanted, save offset */
                      if (current_uniprot_acc[0] == '\0' ||
                          !strcmp(current_uniprot_acc,uniprot_acc))
                        {
                          /* Increment count of index entries */
                          nkeys++;
                      
                          /* Create a new index record */
                          idx_ptr
                            = create_idx_record(&first_idx_ptr,
                                                &last_idx_ptr,
                                                uniprot_acc,offset,
                                                nkeys);

                          /* If only looking for a single UniProt entry,
                             can end read here */
                          if (!strcmp(current_uniprot_acc,uniprot_acc))
                            eof = TRUE;
                        }
                    }
                }
            }
        }

      /* Increment line count */
      nlines++;

      /* Get current file position */
#ifdef TEST
      offset = ftell(file_ptr);
#else
      offset = gztell(gzfile_ptr);
#endif
    }

  /* Close the input file */
#ifdef TEST
#else
  if (gzipped == TRUE)
    gzclose(gzfile_ptr);
  else
#endif
    fclose(file_ptr);

  /* Return the first index pointer */
  *fst_idx_ptr = first_idx_ptr;

  /* Return if file is gzipped */
  *gzip = gzipped;

  /* Return number of index entries stored */
  return(nkeys);
}
/***********************************************************************

compute_cadd_mark  -  Compute the CADD 'severity' marker

***********************************************************************/

int compute_cadd_mark(char *cadd_phred,char **cadd_mark)
{
  char number_string[20];
  
  char *string;
  
  int iband, mark;
  int have_cadd;

  float phred_score;

  /* Initialise variables */
  have_cadd = FALSE;
  store_string(&string,"-");
  *cadd_mark = string;

  /* Get the scores */
  if (strcmp(cadd_phred,"-"))
    phred_score = atof(cadd_phred);
  else
    return(have_cadd);

  /* Get the band the CADD scaores falls into */
  mark = 3;
  for (iband = 1; iband < NCADD_BANDS && mark == 3; iband++)
    if (phred_score < CADD_VALUE[iband])
      mark = iband - 1;

  /* If have a valid mark, save it */
  if (mark > -1)
    {
      /* Convert to string */
      sprintf(number_string,"%d",mark);

      /* Store the string */
      store_string(&string,number_string);

      /* Return the string */
      *cadd_mark = string;

      /* Set flag that we have CADD scores */
      have_cadd = TRUE;
    }

  /* Return flag */
  return(have_cadd);
}
/***********************************************************************

free_cadd_memory  -  Free up memory used by the CADD data

***********************************************************************/

void free_cadd_memory(struct cadd *first_cadd_ptr)
{
  struct cadd *cadd_ptr, *next_cadd_ptr;

  /* Get first cadd record */
  cadd_ptr = first_cadd_ptr;

  /* Loop through the cadd records until done */
  while (cadd_ptr != NULL)
    {
      /* Get pointer to the next cadd record */
      next_cadd_ptr = cadd_ptr->next_cadd_ptr;

      /* Free memory taken by individual fields */
      free(cadd_ptr->variant_id);

      /* Free up memory taken by current cadd */
      free(cadd_ptr);

      /* Go to the next cadd */
      cadd_ptr = next_cadd_ptr;
    }
}
/***********************************************************************

free_gnomad_memory  -  Free up memory used by the GNOMAD data

***********************************************************************/

void free_gnomad_memory(struct gnomad *first_gnomad_ptr)
{
  int ifield;

  struct gnomad *gnomad_ptr, *next_gnomad_ptr;

  /* Get first gnomad record */
  gnomad_ptr = first_gnomad_ptr;

  /* Loop through the gnomad records until done */
  while (gnomad_ptr != NULL)
    {
      /* Get pointer to the next gnomad record */
      next_gnomad_ptr = gnomad_ptr->next_gnomad_ptr;

      /* Free memory taken by individual fields */
      for (ifield = 0; ifield < NGNOMAD_FIELDS; ifield++)
        {
          if (gnomad_ptr->field[ifield] != &null)
            free(gnomad_ptr->field[ifield]);
        }

      /* Free up memory taken by current gnomad */
      free(gnomad_ptr);

      /* Go to the next gnomad */
      gnomad_ptr = next_gnomad_ptr;
    }
}
/***********************************************************************

create_isoform_entry  -  Create a isoform record

***********************************************************************/

struct isoform *create_isoform_entry(struct isoform **fst_isoform_ptr,
                                     struct isoform **lst_isoform_ptr,
                                     char *uniprot_acc,int isoform,
                                     int res_start,int shift,
                                     int nomap_len)
{
  char *string_ptr;

  struct isoform *isoform_ptr;
  struct isoform *first_isoform_ptr, *last_isoform_ptr;

  /* Initialise pointers */
  isoform_ptr = NULL;
  first_isoform_ptr = *fst_isoform_ptr;
  last_isoform_ptr = *lst_isoform_ptr;

  /* Create entry */
  isoform_ptr = (struct isoform *) malloc(sizeof(struct isoform));
  if (isoform_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct isoform\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  store_string(&(isoform_ptr->uniprot_acc),uniprot_acc);
  isoform_ptr->isoform = isoform;
  isoform_ptr->res_start = res_start;
  isoform_ptr->shift = shift;
  isoform_ptr->nomap_len = nomap_len;

  /* Initialise remaining fields */
  isoform_ptr->count = -1;
  isoform_ptr->order = -1;

  /* Initialise pointer to next entry */
  isoform_ptr->next_isoform_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_isoform_ptr == NULL)
    first_isoform_ptr = isoform_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_isoform_ptr->next_isoform_ptr = isoform_ptr;
                      
  /* Save the current residue's pointer */
  last_isoform_ptr = isoform_ptr;

  /* Return current pointers */
  *fst_isoform_ptr = first_isoform_ptr;
  *lst_isoform_ptr = last_isoform_ptr;

  return(isoform_ptr);
}
/***********************************************************************

process_isoform  -  Process the isoform details and create an isoform
                    record

***********************************************************************/

int process_isoform(char *uniprot_acc,char **vstring,
                    struct isoform **fst_isoform_ptr,
                    struct isoform **lst_isoform_ptr,int nisoform,
                    int res_start,int res_end,int *isoform_id,
                    char isoform_name[MAX_ISOFORM][ISOFORM_NAME_LEN + 1],
                    int niso)
{
  char number_string[20], tmp_string[D_MAX_SEQLEN];

  char *string_ptr, *tmp, *var_string;

  int i, ipos, isoform, isoform_store[MAX_ISOFORM], len, len1, len2;
  int istore, nblanks, next, nomap_len, nstored, number, shift, start;
  int done, open_bracket;

  struct isoform *isoform_ptr;
  struct isoform *first_isoform_ptr, *last_isoform_ptr;

  /* Initialise variables */
  isoform = -1;
  nomap_len = nstored = shift = 0;
  for (i = 0; i < MAX_ISOFORM; i++)
    isoform_store[i] = 0;
  open_bracket = FALSE;
  var_string = *vstring;

  /* Remove the note at the front and quote from the back */
  string_ptr = strchr(var_string,'"');
  if (string_ptr != NULL)
    {
      /* Copy string to temporary area */
      strcpy(tmp_string,string_ptr + 1);

      /* Remove closing quote */
      string_ptr = strchr(tmp_string,'"');
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Free up variants string */
      free(var_string);

      /* Copy edited version back */
      store_string(&var_string,tmp_string);
    }

  /* Initialise pointers */
  isoform_ptr = NULL;
  first_isoform_ptr = *fst_isoform_ptr;
  last_isoform_ptr = *lst_isoform_ptr;

  /* Extract the isoform number(s) */
  string_ptr = strstr(var_string," (in isoform ");
  if (string_ptr != NULL)
    {
      /* Save the string from the number onwards */
      store_string(&tmp,string_ptr);

      /* Terminate variants string before opening bracket */
      string_ptr[0] = '\0';

      /* Initialise flag */
      done = FALSE;
      start = 0;

      /* Check for " and " and replace by comma */
      string_ptr = strstr(tmp," and ");
      if (string_ptr != NULL)
        string_ptr[0] = ',';

      /* Loop until all isoforms have been extracted */
      while (done == FALSE)
        {
          /* Set flag */
          open_bracket = FALSE;

          /* Get the next isoform number */
          string_ptr = strstr(tmp + start,"isoform ");

          /* If we have one, then extract */
          if (string_ptr != NULL)
            {
              /* Get the number end */
              len = -1;
              for (ipos = 0; ipos < strlen(string_ptr + 8) && len == -1;
                   ipos++)
                {
                  /* If this is the end of the number, save its length */
                  if (string_ptr[ipos + 8] == ',' ||
                      (string_ptr[ipos + 8] == ')' &&
                       open_bracket == FALSE))
                    len = ipos;
                  else if (string_ptr[ipos + 8] == '(')
                    open_bracket = TRUE;
                }

              /* If we have the identifier end, can get its length */
              if (len > 0)
                {
                  /* Extract the id */
                  if (len <= ISOFORM_NAME_LEN)
                    {
                      /* Copy isoform name, skipping any spaces */
                      number_string[0] = '\0';
                      next = 0;
                      for (i = 0; i < len; i++)
                        {
                          if (string_ptr[i + 8] != ' ')
                            {
                              number_string[next] = string_ptr[i + 8];
                              next++;
                              number_string[next] = '\0';
                            }
                        }
                      //                      strncpy(number_string,string_ptr + 8,len);
                      //number_string[len] = '\0';

                      /* Determine which isoform this is */
                      number = -1;
                      for (i = 0; i < niso && number == -1; i++)
                        {
                          /* If correct one, save the number */
                          if (!strcmp(isoform_name[i],number_string))
                            number = i;
                        }

                      /* If we have a match, save the isoform id */
                      if (number > -1)
                        isoform = isoform_id[number];
                    }
                  else
                    isoform = -1;

                  /* Store the number */
                  if (nstored < MAX_ISOFORM && isoform > 0)
                    {
                      isoform_store[nstored] = isoform;
                      nstored++;
                    }

                  /* Shift the start position form next isoform */
                  start = string_ptr - tmp + 8 + len + 1;
                }

              /* Otherwise, end here */
              else
                done = TRUE;
            }

          /* If next isoform number not found, we are done */
          else
            done = TRUE;
        }
      
      /* Free up memory taken by temporary string */
      free(tmp);
    }

  /* If no isoform number identified, then return */
  if (nstored == 0)
    {
      /* Return the modified variant string */
      *vstring = var_string;
      return(nisoform);
    }

  /* If isoform has a missing set of residues, then calculate shift */
  if (!strcmp(var_string,"Missing"))
    {
      /* Get the shift to the residue numbering at this position */
      shift = res_end - res_start + 1;
    }

  /* For a substitution, calculate the shift */
  else if ((string_ptr = strstr(var_string," -> ")) != NULL)
    {
      /* Get the length of the original sequence substituted */
      len1 = string_ptr - var_string;
      nblanks = 0;

      /* Loop through the string to reduce count by any spaces found */
      for (ipos = 0; ipos < len1; ipos++)
        {
          if (var_string[ipos] == ' ')
            nblanks++;
        }
      len1 = len1 - nblanks;

      /* Get the length of the substituting sequence */
      len2 = strlen(string_ptr + 4);
      nblanks = 0;

      /* Loop through the string to reduce count by any spaces found */
      for (ipos = 0; ipos < len2; ipos++)
        {
          if (string_ptr[ipos + 4] == ' ')
            nblanks++;
        }
      len2 = len2 - nblanks;

      /* The latter length corresponds to the length that can't be
         mapped to the canonical isoform */
      nomap_len = len2;

      /* The length difference constitutes the shift in residue nunbering
         (after the no-map length) */
      shift = len1 - len2;
    }

  /* Loop over the stored isoforms to create a record for each */
  for (istore = 0; istore < nstored; istore++)
    {
      /* Get the stored isoform identifier */
      isoform = isoform_store[istore];

      /* Create an isoform record */
      isoform_ptr
        = create_isoform_entry(&first_isoform_ptr,&last_isoform_ptr,
                               uniprot_acc,isoform,res_start,shift,
                               nomap_len);

      /* Increment count */
      nisoform++;
    }

  /* Return current pointers */
  *fst_isoform_ptr = first_isoform_ptr;
  *lst_isoform_ptr = last_isoform_ptr;

  /* Return the modified variant string */
  *vstring = var_string;

  /* Return number of isoforms */
  return(nisoform);
}
/***********************************************************************

isoform_sort  -  Shell-sort which sorts the isoforms by UniProt id,
                 isoform no, and start residue

***********************************************************************/

void isoform_sort(long *index,struct isoform **isoform_index_ptr,
                  long nindex)
{
  int lognb2;
  int swap;
  long endloop, i, indi, indl, iswap, j, k, l, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct isoform *isoform1_ptr, *isoform2_ptr;

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

              /* Get the two isoform records */
              isoform1_ptr = isoform_index_ptr[index[indi]];
              isoform2_ptr = isoform_index_ptr[index[indl]];

              /* See if they are in the wrong order */
              swap = FALSE;
              if (strcmp(isoform1_ptr->uniprot_acc,
                         isoform2_ptr->uniprot_acc) > 0)
                swap = TRUE;
              else if (!strcmp(isoform1_ptr->uniprot_acc,
                               isoform2_ptr->uniprot_acc))
                {
                  if (isoform1_ptr->isoform > isoform2_ptr->isoform)
                    swap = TRUE;
                  else if (isoform1_ptr->isoform == isoform2_ptr->isoform)
                    {
                      /* Sort by start residue */
                      if (isoform1_ptr->res_start > 
                          isoform2_ptr->res_start)
                        swap = TRUE;
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

sort_isoforms  -  Sort the isoforms by their UniProt id

***********************************************************************/

long sort_isoforms(struct isoform **fst_isoform_ptr,long nisoform,
                   struct isoform ***isoform_idx_ptr)
{
  char last_uniprot_acc[UNIPROT_LEN];

  int count, last_isoform, nuniprot, order;

  long iorder, ipos;
  long *index;

  struct isoform *first_isoform_ptr, *last_isoform_ptr;
  struct isoform *isoform_ptr;
  struct isoform **isoform_index_ptr;

  /* Initialise variables */
  nuniprot = 0;

  /* Initialise pointers */
  first_isoform_ptr = *fst_isoform_ptr;
  last_isoform_ptr = NULL;

  /* Grab memory for temporary sort arrays */
  index = (long *) malloc(nisoform * sizeof(long));
  isoform_index_ptr
    = (struct isoform **) malloc(nisoform * sizeof(struct isoform *));

  /* Get pointer to first isoform record */
  isoform_ptr = first_isoform_ptr;
  ipos = 0;

  /* Loop through isoform records to initialise the sort index */
  while (isoform_ptr != NULL)
    {
      /* Initialise sort index */
      index[ipos] = ipos;
      isoform_index_ptr[ipos] = isoform_ptr;

      /* Increment isoform counter */
      ipos++;

      /* Get next isoform record */
      isoform_ptr = isoform_ptr->next_isoform_ptr;
    }

  /* Initialise isoform counter */
  nisoform = ipos;

  /* Sort the isoforms */
  isoform_sort(index,isoform_index_ptr,nisoform);

  /* Having sorted the isoforms, rearrange order of the isoform-record
     pointers in descending order of isoform */
  for (iorder = 0; iorder < nisoform; iorder++)
    {
      /* Get the array position of this isoform record */
      ipos = index[iorder];

      /* Get the current isoform record */
      isoform_ptr = isoform_index_ptr[ipos];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_isoform_ptr = isoform_ptr;

      /* Otherwise, get the previous isoform record to point to this one */
      else
        last_isoform_ptr->next_isoform_ptr = isoform_ptr;

      /* Blank out current record's next isoform pointer */
      isoform_ptr->next_isoform_ptr = NULL;

      /* Store current record as last encountered */
      last_isoform_ptr = isoform_ptr;
    }

  /* Free up memory used for temporary sort arrays */
  free(index);

  /* Initialise */
  last_uniprot_acc[0] = '\0';
  last_isoform = 0;
  count = order = 0;

  /* Get pointer to first isoform record */
  isoform_ptr = first_isoform_ptr;
  ipos = 0;

  /* Loop through isoform records to refill index */
  while (isoform_ptr != NULL)
    {
      /* Initialise sort index */
      isoform_index_ptr[ipos] = isoform_ptr;

      /* Increment code counter */
      ipos++;

      /* If a new UniProt id, increment count */
      if (strcmp(isoform_ptr->uniprot_acc,last_uniprot_acc))
        nuniprot++;

      /* Check if UniProt id or isoform have changed */
      if (strcmp(isoform_ptr->uniprot_acc,last_uniprot_acc) ||
          isoform_ptr->isoform != last_isoform)
        {
          /* If the UniProt id has changed, reset the count */
          if (strcmp(isoform_ptr->uniprot_acc,last_uniprot_acc))
            count = 0;

          /* Reset the order number */
          order = 0;

          /* Save the UniProt id and isoform */
          strcpy(last_uniprot_acc,isoform_ptr->uniprot_acc);
          last_isoform = isoform_ptr->isoform;
        }

      /* Save the count and order number, and increment */
      isoform_ptr->count = count;
      isoform_ptr->order = order;
      count++;
      order++;

      /* Get next code record */
      isoform_ptr = isoform_ptr->next_isoform_ptr;
    }

  /* Return current pointer */
  *fst_isoform_ptr = first_isoform_ptr;

  /* Return the isoform index */
  *isoform_idx_ptr = isoform_index_ptr;

  /* Show stats */
  printf("Number of UniProt sequences with isoforms:   %8d\n",nuniprot);
  fflush(stdout);
 
  /* Return number of sorted isoforms */
  return(nisoform);
}
/***********************************************************************

get_isoforms  -  Read in the UniProt sequence isoforms from the
                 uniprot.txt file

***********************************************************************/

int get_isoforms(char *file_name,char *current_uniprot_acc,
                 struct isoform ***isoform_idx_ptr)
{
  char input_line[D_LONG_LINELEN + 1], uniprot_acc[UNIPROT_LEN];
  char iso_name[ISOFORM_NAME_LEN + 1], number_string[D_LINELEN];
  char isoform_name[MAX_ISOFORM][ISOFORM_NAME_LEN + 1];

  char *string_ptr, *token[D_MAX_TOKEN], *var_string;

  int errno, ipos, len, nisoform, ntokens, res_end, res_start;
  int i, id, isoform, isoform_id[MAX_ISOFORM];
  int done, eof, gzipped, last_AC, have_ours, have_varseq, read_all;
  int ok, wanted;

  long line;

  struct isoform *isoform_ptr, *first_isoform_ptr, *last_isoform_ptr;

  struct isoform **isoform_index_ptr;

  FILE *file_ptr;
#ifdef TEST
#else
  gzFile gzfile_ptr;
#endif

  /* Initialise variables */
  done = FALSE;
  eof = FALSE;
  gzipped = FALSE;
  have_ours = FALSE;
  have_varseq = FALSE;
  isoform = 0;
  iso_name[0] = '\0';
  last_AC = FALSE;
  line = 0;
  nisoform = 0;
  read_all = FALSE;
  res_start = res_end = 0;
  uniprot_acc[0] = '\0';
  var_string = &null;
  wanted = FALSE;

  /* If all UniProt entries to be read in, set flag */
  if (current_uniprot_acc[0] == '\0')
    read_all = TRUE;

  /* Initialise pointers */
  first_isoform_ptr = last_isoform_ptr = NULL;
  isoform_index_ptr = NULL;

  /* Determine if this is a gzipped file */
  if (strstr(file_name,".gz") != NULL)
    gzipped = TRUE;

#ifdef TEST
  if (gzipped == TRUE)
    {
      printf("*** ERROR. Cannot open gzipped file in TEST mode. Aborting\n");
      exit(1);
    }
#endif

  /* Open the UniProt file */
#ifdef TEST
#else
  if (gzipped == TRUE)
    gzfile_ptr = gzopen(file_name,"r");
  else
#endif
    file_ptr = fopen(file_name,"r");

  /* If file not found, then abort */
  if ((gzipped == FALSE && file_ptr == NULL)
#ifdef TEST
#else
      || (gzipped == TRUE && gzfile_ptr == NULL)
#endif
      )   
    {
      /* Show error message and abort */
      printf("\n*** ERROR. Unable to open file [%s]\n",file_name);
      exit(1);
    }

  /* Loop until we have finished reading the file */
  while (eof == FALSE && done == FALSE)
    {
      /* Read in the next record */
#ifdef TEST
#else
      if (gzipped == TRUE)
        {
          if (gzgets(gzfile_ptr,input_line,D_LONG_LINELEN) == NULL)
            eof = TRUE;
        }
      else
#endif
        if (fgets(input_line,D_LONG_LINELEN,file_ptr) == NULL)
            eof = TRUE;

      /* If record read in OK, then process it */
      if (eof == FALSE)
        {
          /* If this is the AC line, get the UniProt id */
          if (!strncmp(input_line,"AC   ",5))
            {
              /* Only process if last record wasn't an AC one */
              if (last_AC == FALSE)
                {
                  /* Truncate the input line */
                  len = string_chomp(input_line);

                  /* Find semi-colon */
                  string_ptr = strchr(input_line + 5,';');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Extract the UniProt accession */
                  strcpy(uniprot_acc,input_line + 5);

                  /* Check if this UniProt entry is wanted */
                  if (read_all == TRUE ||
                      !strcmp(current_uniprot_acc,uniprot_acc))
                    wanted = TRUE;
                  else
                    wanted = FALSE;

                  /* If we have already processed our UniProt code,
                     then we are done */
                  if (have_ours == TRUE)
                    done = TRUE;
                }

              /* Re-initialise flags */
              last_AC = TRUE;
              have_varseq = FALSE;

              /* Initialise isoform ids */
              for (i = 0; i < MAX_ISOFORM; i++)
                {
                  isoform_id[i] = 0;
                  isoform_name[i][0] = '\0';
                }
              isoform = 0;
            }
          else
            last_AC = FALSE;

          /* Only process the line if it is wanted */
          if (wanted == TRUE)
            {
              /* If this is not the continuation of a VAR_SEQ entry, 
                 then process previous one */
              if (have_varseq == TRUE &&
                  strncmp(input_line,"FT           ",13))
                {
                  /* Save previous isoform */
                  nisoform
                    = process_isoform(uniprot_acc,&var_string,
                                      &first_isoform_ptr,&last_isoform_ptr,
                                      nisoform,res_start,res_end,
                                      isoform_id,isoform_name,isoform);
              
                  /* Initialise variables */
                  have_varseq = FALSE;
                  if (var_string != &null)
                    {
                      free(var_string);
                      var_string = &null;
                    }
                  iso_name[0] = '\0';

                  /* If this is our UniProt code, then set flag */
                  if (read_all == FALSE &&
                      !strcmp(current_uniprot_acc,uniprot_acc))
                    have_ours = TRUE;
                }

              /* If this is a new VAR_SEQ entry, grab details and initialise */
              if (!strncmp(input_line,"FT   VAR_SEQ ",13))
                {
                  /* Truncate the input line */
                  len = string_chomp(input_line);
                  
                  /* Get the tokens from the input line */
                  ntokens = get_tokens(input_line,token,D_MAX_TOKEN,' ');
              
                  /* If have enough tokens, then interpret and save */
                  if (ntokens == 3)
                    {
                      /* See if last token has a residue range */
                      string_ptr = strstr(token[2],"..");
                      if (string_ptr != NULL)
                        {
                          /* Get the second number */
                          strcpy(number_string,string_ptr + 2);
                          res_end = atoi(number_string);

                          /* Get the first number */
                          string_ptr[0] = '\0';
                          res_start = atoi(token[2]);
                        }

                      /* If only a single number, set as start and end */
                      else
                        res_start = res_end = atoi(token[2]);
                    }

                  /* Free up memory taken by tokens */
                  free_tokens(token,ntokens);

                  /* Set flag */
                  have_varseq = TRUE;
                }

              /* If this is a Name line, save the name */
              else if (!strncmp(input_line,"CC       Name=",14))
                {
                  /* Truncate the input line */
                  len = string_chomp(input_line);
                  
                  /* Chop off line at first space or semi-colon */
                  string_ptr = strchr(input_line + 14,' ');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  string_ptr = strchr(input_line + 14,';');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Check the length of the id */
                  len = strlen(input_line + 14);

                  /* Extract the isoform id */
                  if (len < ISOFORM_NAME_LEN)
                    {
                      /* Save the isoform name */
                      strcpy(iso_name,input_line + 14);
                    }
                }

              /* If this is an isoform id line, save the id */
              else if (!strncmp(input_line,"CC         IsoId=",17))
                {
                  /* Truncate the input line */
                  len = string_chomp(input_line);
                  
                  /* Chop off line at first space or semi-colon */
                  string_ptr = strchr(input_line + 17,' ');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';
                  string_ptr = strchr(input_line + 17,';');
                  if (string_ptr != NULL)
                    string_ptr[0] = '\0';

                  /* Find the start of the isoform id */
                  string_ptr = strchr(input_line,'-');
                  if (string_ptr != NULL)
                    {
                      /* Check the length of the id */
                      len = strlen(string_ptr + 1);

                      /* Extract the isoform id */
                      if (len < 10 && iso_name[0] != '\0')
                        {
                          strcpy(number_string,string_ptr + 1);
                          id = atoi(number_string);

                          /* If have a valid number and id, save */
                          if (isoform > -1 && isoform < MAX_ISOFORM &&
                              id > 0 && id < MAX_ISOFORM)
                            {
                              isoform_id[isoform] = id;
                              strcpy(isoform_name[isoform],iso_name);
                            }
                          isoform++;
                        }
                    }
                }

              /* Otherwise, if we have the continuation of a variants
                 sequence record, append */
              else if (have_varseq == TRUE)
                {
                  /* Truncate the input line */
                  len = string_chomp(input_line);

                  /* Get the tokens from the input line */
                  ntokens = get_tokens(input_line,token,D_MAX_TOKEN,' ');
              
                  /* Decide if this line is wanted */
                  ok = TRUE;
                  if (ntokens > 1 &&
                      (!strncmp(token[1],"/evidence=",10) ||
                       !strncmp(token[1],"/id=",4) ||
                       !strncmp(token[1],"ECO:",4)))
                    ok = FALSE;

                  /* If line wanted, append all the tokens */
                  if (ok == TRUE)
                    {
                      /* Get start position of the dara */
                      string_ptr = strstr(input_line,token[1]);

                      /* If OK, then add to variant sequence */
                      if (string_ptr != NULL)
                        {
                          /* Get start position */
                          ipos = string_ptr - input_line;

                          /* If the first, then store */
                          if (var_string == &null)
                            store_string(&var_string,input_line + ipos);
                          else
                            {
                              /* Append additional space */
                              append_string(&var_string," ");

                              /* Append the variation line */
                              append_string(&var_string,input_line + ipos);
                            }
                        }
                    }

                  /* Free up memory taken by tokens */
                  free_tokens(token,ntokens);
                }
            }
        }
    }

  /* If have a last VAR_SEQ entry, process it */
  if (have_varseq == TRUE)
    nisoform = process_isoform(uniprot_acc,&var_string,
                               &first_isoform_ptr,&last_isoform_ptr,
                               nisoform,res_start,res_end,isoform_id,
                               isoform_name,isoform);

  /* Close the input file */
#ifdef TEST
#else
  if (gzipped == TRUE)
    gzclose(gzfile_ptr);
  else
#endif
    fclose(file_ptr);

  /* Sort the isoform records */
  if (nisoform > 0)
    nisoform
      = sort_isoforms(&first_isoform_ptr,nisoform,&isoform_index_ptr);
  else
    nisoform = -1;

  /* Return the isoform index */
  *isoform_idx_ptr = isoform_index_ptr;

  /* Return number of records read in */
  return(nisoform);
}
/***********************************************************************

read_isoforms_txt  -  Read in the UniProt sequence isoforms from the
                      isoforms.txt file

***********************************************************************/

int read_isoforms_txt(struct isoform ***isoform_idx_ptr,
                      struct c_file_data file_name_ptr)
{
  char file_name[D_FILENAME_LEN], input_line[D_LINELEN + 1];

  char *string_ptr, *token[D_MAX_TOKEN];

  int len, line, nisoform, ntokens;

  struct isoform *isoform_ptr, *first_isoform_ptr, *last_isoform_ptr;

  struct isoform **isoform_index_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  nisoform = 0;

  /* Initialise pointers */
  first_isoform_ptr = last_isoform_ptr = NULL;
  *isoform_idx_ptr = isoform_index_ptr = NULL;

  /* Form the name of the isoforms.txt file */
  strcpy(file_name,file_name_ptr.other_dir[PDBSUM_REAL_DATA_DIR]);
  strcat(file_name,"/isoforms.txt");

  /* Open the isoforms.txt file */
  file_ptr = fopen(file_name,"r");
  if (file_ptr == NULL)
    {
      printf("*** Warning. File not found: %s\n",file_name);
      return(nisoform);
    }
  printf("Reading: %s\n",file_name);

  /* Loop through the input file */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL)
    {
      /* Truncate line and get its length */
      len = string_chomp(input_line);
      line++;

      /* Extract the space-separated tokens */
      ntokens = get_tokens(input_line,token,D_MAX_TOKEN,'\t');

      /* If have enough tokens on this line, create an isoform record */
      if (ntokens == NISOFORM_FIELDS)
        {
          /* Create an isoform record */
          isoform_ptr
            = create_isoform_entry(&first_isoform_ptr,&last_isoform_ptr,
                                   token[ISOFORM_UNIPROT_ACC],
                                   atoi(token[ISOFORM_ISOFORM]),
                                   atoi(token[ISOFORM_RES_START]),
                                   atoi(token[ISOFORM_SHIFT]),
                                   atoi(token[ISOFORM_NOMAP_LEN]));

          /* Increment count */
          nisoform++;
        }

      else
        printf("*** Warning. Wrong number of fields (%d) in isoforms.txt "
               "at line %d\n",ntokens,line);
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Sort the isoform records */
  if (nisoform > 0)
    nisoform
      = sort_isoforms(&first_isoform_ptr,nisoform,&isoform_index_ptr);
  else
    nisoform = -1;

  /* Return the isoform index */
  *isoform_idx_ptr = isoform_index_ptr;

  /* Return number of records read in */
  return(nisoform);
}
/***********************************************************************

get_label_letters  -  Get the letters making up this label and return its
                      width

***********************************************************************/

int get_label_letters(char *label,struct font font_ptr,
                      struct letter ***letter_arr_ptr,int letter_gap,
                      int *sty,int *edy,int *lab_height)
{
  char ch;

  char *string_ptr;
  
  int i, ichar, iletter, irow, starty, endy;
  int label_len, label_height, label_width, letter_height;

  int *used;

  struct letter *letter_ptr;

  struct letter **letter_array_ptr;

  /* Initialise variables */
  label_height = label_width = 0;
  letter_height = font_ptr.height;

  /* Get the label's length */
  label_len = strlen(label);

  /* Create an array to hold the letters in the name */
  letter_array_ptr
    = (struct letter **) malloc(label_len * sizeof(struct letter *));
  for (i = 0; i < label_len; i++)
    letter_array_ptr[i] = NULL;

  /* Initialise the label width */
  label_width = 0;

  /* Create temporary array to hold used rows */
  used = (int *) malloc(letter_height * sizeof(int));

  /* Initialise array indicating which label array rows are in use */
  for (irow = 0; irow < letter_height; irow++)
    used[irow] = FALSE;

  /* Loop over the characters making up the label */
  for (ichar = 0; ichar < label_len; ichar++)
    {
      /* Get this character */
      ch = label[ichar];

      /* If not a space, get the corresponding definition */
      if (ch != ' ')
        {
          /* Find the letter's font definition */
          string_ptr = strchr(font_ptr.character,ch);
          if (string_ptr != NULL)
            iletter = string_ptr - font_ptr.character;
          else
            iletter = CAPITAL_X;

          /* Get the corresponding character definition */
          letter_ptr = font_ptr.letter_index_ptr[iletter];
          letter_array_ptr[ichar] = letter_ptr;
 
          /* Add character width to label width */
          label_width 
            = label_width + letter_array_ptr[ichar]->width + letter_gap;

          /* Determine which rows this letter will use */
          for (irow = letter_ptr->miny; irow < letter_ptr->maxy + 1; irow++)
            used[irow] = TRUE;
        }
    }

  /* Add final gap to letter width */
  label_width = label_width + letter_gap;

  /* Determine range of rows in the label array that contain
     character parts, and consequently the label height */
  starty = 0;
  endy = letter_height - 1;
  for (i = 0; i < letter_height && used[i] == FALSE; i++)
    starty = i + 1;
  for (i = letter_height - 1; i > -1 && used[i] == FALSE; i--)
    endy = i - 1;
  label_height = endy - starty + 1;

  /* Return the letter array */
  *letter_arr_ptr = letter_array_ptr;

  /* Free memory */
  free(used);

  /* Return the computed parameters */
  *sty = starty;
  *edy = endy;
  *lab_height = label_height;
  
  /* Return the label width */
  return(label_width);
}
/***********************************************************************

create_cadd_entry  -  Create a record for this CADD variant

***********************************************************************/

struct cadd *create_cadd_entry(struct cadd **fst_cadd_ptr,
                               struct cadd **lst_cadd_ptr,
                               char *variant_id,char *cadd_raw,
                               char *cadd_phred)
{
  char *string_ptr;

  struct cadd *cadd_ptr;
  struct cadd *first_cadd_ptr, *last_cadd_ptr;

  /* Initialise pointers */
  cadd_ptr = NULL;
  first_cadd_ptr = *fst_cadd_ptr;
  last_cadd_ptr = *lst_cadd_ptr;

  /* Create entry */
  cadd_ptr = (struct cadd *) malloc(sizeof(struct cadd));
  if (cadd_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct cadd\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  store_string(&(cadd_ptr->variant_id),variant_id);

  /* Initialise remaining fields */
  cadd_ptr->cadd_raw = atof(cadd_raw);
  cadd_ptr->cadd_phred = atof(cadd_phred);

  /* Initialise pointer to next entry */
  cadd_ptr->next_cadd_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_cadd_ptr == NULL)
    first_cadd_ptr = cadd_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_cadd_ptr->next_cadd_ptr = cadd_ptr;
                      
  /* Save the current residue's pointer */
  last_cadd_ptr = cadd_ptr;

  /* Return current pointers */
  *fst_cadd_ptr = first_cadd_ptr;
  *lst_cadd_ptr = last_cadd_ptr;
  return(cadd_ptr);
}
/***********************************************************************

read_cadd_txt  -  Read in the CADD values from cadd.txt for the given
                  UniProt code

***********************************************************************/

int read_cadd_txt(struct cadd **fst_cadd_ptr,struct cadd **lst_cadd_ptr,
                  char *dir)
{
  char file_name[D_FILENAME_LEN + 1], input_line[D_LINELEN + 1];
  char command[D_FILENAME_LEN + 1];

  char *token[D_MAX_TOKEN];

  int len, line, ncadd, ntokens;

  struct cadd *cadd_ptr, *first_cadd_ptr, *last_cadd_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = ncadd = 0;

  /* Initialise pointers */
  cadd_ptr = first_cadd_ptr = last_cadd_ptr = NULL;
  *fst_cadd_ptr = *lst_cadd_ptr = NULL;

  /* Form name of the appropriate cadd_XXXXXX.txt file */
  strcpy(file_name,dir);
  strcat(file_name,"/cadd.txt");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",
             file_name);
      return(ncadd);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);

      /* Split the list up into its tokens */
      ntokens = get_tokens(input_line,token,D_MAX_TOKEN,'\t');

      /* If have enough tokens, save this record */
      if (ntokens > 2)
        {
          /* Create a cadd record */
          cadd_ptr
            = create_cadd_entry(&first_cadd_ptr,&last_cadd_ptr,
                                token[0],token[1],token[2]);

          /* Increment count of cadds */
          ncadd++;
        }

      /* Free up the memory taken by the tokens */
      free_tokens(token,ntokens);
    }

  /* Return current pointers */
  *fst_cadd_ptr = first_cadd_ptr;
  *lst_cadd_ptr = last_cadd_ptr;

  /* Close the input file */
  fclose(file_ptr);

  /* Return number of cadds read in */
  return(ncadd);
}
/***********************************************************************

read_in_table  -  Read in the table of variant counts from the
                  appropriate .tsv file

***********************************************************************/

int read_in_table(char *file_name,int *table,int *tot)
{
  char input_line[D_LINELEN + 1];

  char *string_ptr, *token[D_MAX_TOKEN];

  int iamino, ipos, itoken, jamino, len, line, ntokens, total;
  int ok;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  ok = FALSE;
  total = 0;

  /* Initialise amino acid positions */
  ipos = 0;
  
  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      /* If can't find file, then return empty-handed */
      printf("*** Warning. Unable to find input file: %s\n",
             file_name);
      return(FALSE);
    }

  /* Read in and ignore the first line */
  if (fgets(input_line,D_LINELEN,file_ptr) == NULL)
    return(FALSE);

  /* Loop while reading in records from the file */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_chomp(input_line);
      
      /* Split the line up into its tokens */
      ntokens = get_tokens(input_line,token,D_MAX_TOKEN,'\t');

      /* If have enough tokens, extract the counts */
      if (ntokens > NAMINO)
        {
          /* Increment line count */
          line++;

          /* Loop over the tokens, to store */
          for (itoken = 1; itoken < ntokens; itoken++)
            {
              /* Save in the table */
              table[ipos] = atoi(token[itoken]);

              /* Add to total */
              total = total + table[ipos];

              /* Increment subscript */
              ipos++;
            }
        }
      
      /* Free up the memory taken by the tokens */
      free_tokens(token,ntokens);
    }

  /* Close the input file */
  fclose(file_ptr);

  /* If have read in the right number of lines, set flag */
  if (line == NAMINO)
    ok = TRUE;

  /* Return this table's total */
  *tot = total;

  /* Return if all OK */
  return(ok);
}
/***********************************************************************

get_tables  -  Read in the table of variant counts from the appropriate
               .tsv file

***********************************************************************/

int get_tables(char *data_dir,int *base_changes,float *disease_propensity)
{
  char file_name[D_FILENAME_LEN + 1];

  static char *TSV_FILE[] = { "codon_changes", "counts", "variants" };

  int ipos, itable;
  int table[NTABLE_TYPES][NAMINO * NAMINO], total[NTABLE_TYPES];
  int ok;

  float expected;

  /* Initialise tables of counts */
  for (itable = 0; itable < NTABLE_TYPES; itable++)
    for (ipos = 0; ipos < NAMINO * NAMINO; ipos++)
      table[itable][ipos] = 0;

  /* Initialise final tables */
  for (ipos = 0; ipos < NAMINO * NAMINO; ipos++)
    disease_propensity[ipos] = base_changes[ipos] = 0;

  /* Read in the tables of counts from the .tsv files */
  for (itable = 0; itable < NTABLE_TYPES; itable++)
    {
      /* Form the name of the file to be read */
      strcpy(file_name,data_dir);
      strcat(file_name,"/");
      strcat(file_name,TSV_FILE[itable]);
      strcat(file_name,".tsv");

      /* Read in the counts from the table */
      ok = read_in_table(file_name,table[itable],&total[itable]);

      /* If table can't be read in, abort */
      if (ok == FALSE)
        exit(1);
    }

  /* Copy across the base changes data */
  for (ipos = 0; ipos < NAMINO * NAMINO; ipos++)
    base_changes[ipos] = table[BASE_CHANGES][ipos];

  /* Calculate the disease propensities */
  for (ipos = 0; ipos < NAMINO * NAMINO; ipos++)
    {
      /* Initialise to zero */
      disease_propensity[ipos] = 0;

      /* If this change is a single-base change, can compute the
         disease propensity */
      if (base_changes[ipos] == 1 && table[GNOMAD_VARIANTS][ipos] > 0)
        {
          /* Compute the expected count of disease variants, based on
             the natural variant counts */
          expected = table[GNOMAD_VARIANTS][ipos]
            * (float) total[UNIPROT_VARIANTS]
            / (float) total[GNOMAD_VARIANTS];

          /* Calculate the disease propensity for this type of a.a.
             change */
          disease_propensity[ipos]
            = (float) table[UNIPROT_VARIANTS][ipos] / expected;
        }
    }

  /* Return whether tables read in OK */
  return(ok);
}
/***********************************************************************

plot_label  -  Plot given label

***********************************************************************/

void plot_label(int *array,int nrows,int ncols,int xorigin,
                int yorigin,struct letter **letter_ptr,int lstart,
                int nletters,int *trgb,int letter_gap,int space_width,
                int *map,int flip)
{
  char ch;

  int ichar, icol, irow, x, xpos, xstart, ypos, ystart;
  int bg_colour, letter_height, letter_width;
  int first;

  long iloc, ipos;

  float grey_shade;

  /* Initialise variables */
  first = TRUE;
  letter_height = 0;
  letter_width = 0;
  xstart = xorigin;
  ystart = yorigin;

  /* Loop over the characters making up the label */
  for (ichar = 0; ichar < nletters; ichar++)
    {
      /* If not a space, then plot */
      if (letter_ptr[ichar] != NULL)
        {
          /* Save the letter height */
          letter_height = letter_ptr[ichar]->height;
          letter_width = letter_ptr[ichar]->width;

          /* Loop over the rows and columns of the letter */
          for (irow = 0; irow < letter_height; irow++)
            {
              for (icol = lstart; icol < letter_width; icol++)
                {
                  /* Get position in the letter array */
                  ipos = irow * letter_width + icol;

                  /* Get position in image array */
                  xpos = xstart + icol;

                  /* Flip y-value, if necessary */
                  if (flip == TRUE)
                    ypos = ystart - letter_height + irow;
                  else
                    ypos = ystart - irow;
                  iloc = ypos * ncols + xpos;

                  /* If within bounds, then get the pixel colour */
                  if (xpos >= 0 && xpos < ncols && ypos >= 0 && ypos < nrows)
                    {
                      /* Get this pixel */
                      ch = letter_ptr[ichar]->array[ipos];

                      /* If not empty, apply to image */
                      if (ch != ' ')
                        {
                          /* Get the grey-shade for this letter pixel */
                          grey_shade = get_pixel_shade(ch);

                          /* Get the current background colour at
                             this point */
                          bg_colour = array[iloc];

                          /* Apply pixel to array */
                          array[iloc]
                            = apply_pixel(bg_colour,grey_shade,trgb);

                          /* Adjust maximum and minimum y-values */
                          if (first == TRUE)
                            {
                              /* Store current point as max and min */
                              map[XMIN] = map[XMAX] = xpos;
                              map[YMIN] = map[YMAX] = ypos;

                              /* Unset flag */
                              first = FALSE;
                            }
                          else
                            {
                              /* Update max and min values */
                              if (xpos < map[XMIN])
                                map[XMIN] = xpos;
                              if (xpos > map[XMAX])
                                map[XMAX] = xpos;
                              if (ypos < map[YMIN])
                                map[YMIN] = ypos;
                              if (ypos > map[YMAX])
                                map[YMAX] = ypos;
                            }
                        }
                    }
                }
            }
        }
      else
        {
          /* Shift x-position */
          letter_width = 0;
          xpos = xstart + space_width - 1;
        }

      /* Increment x position in readiness for next character */
      xstart = xpos + letter_gap + 1;
    }
}
/***********************************************************************

create_array  -  Create an image array and initialise to white

***********************************************************************/

int *create_array(int nrows,int ncols,int init_val)
{
  int *array;

  long iloc, size;

  /* Get the image size */
  size = nrows * ncols;

  /* Create array to for the domain diagram */
  array = (int *) malloc(size * sizeof(int));
  if (array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for domain image "
             "array: %d x %d\n",ncols,nrows);
      exit(1);
    }

  /* Initialise array */
  for (iloc = 0; iloc < size; iloc++)
    array[iloc] = init_val;

  /* Return the array */
  return(array);
}
/***********************************************************************

read_gnomad_idx  -  Read gnomAD.idx index file to pick up the starting
                    offset for the corresponding gnomAD data file

***********************************************************************/

int read_gnomad_idx(struct idx **fst_idx_ptr,
                    struct c_file_data file_name_ptr)
{
  char file_name[D_FILENAME_LEN], input_line[D_LINELEN + 1];

  char *token[D_MAX_TOKEN];

  int nentry, nindex, nlines, ntokens;

  long offset;

  struct idx *idx_ptr, *first_idx_ptr, *last_idx_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nindex = 0;
  offset = -1;

  /* Initialise idx pointers */
  idx_ptr = *fst_idx_ptr = first_idx_ptr = last_idx_ptr = NULL;

  /* Form the name of the index file */
  strcpy(file_name,file_name_ptr.other_dir[GNOMAD_DIR]);
  strcat(file_name,"/gnomAD.idx");

  /* Open the index file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      printf("*** Warning. Unable to open data file [%s]\n",file_name);
      return(offset);
    }

  /* Read the file until we hit the right code */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL)
    {
      /* Split input line into tokens */
      ntokens = tokenize(input_line,token,D_MAX_TOKEN,' ');

      /* If have the right number of tokens, create an idx record */
      if (ntokens == 3)
        {
          /* Get the offset and number of entries */
          offset = atol(token[1]);
          nentry = atoi(token[2]);

          /* Create an idx record */
          idx_ptr = create_idx_record(&first_idx_ptr,&last_idx_ptr,
                                      token[0],offset,nentry);

          /* Increment count of records created */
          nindex++;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_idx_ptr = first_idx_ptr;

  /* Return the number of index entries */
  return(nindex);
}
/***********************************************************************

read_idx_file  -  Read .idx index file to pick up the starting offsets
                  for the corresponding .dat file

***********************************************************************/

int read_idx_file(char *file_name,struct idx **fst_idx_ptr)
{
  char input_line[D_LINELEN + 1];

  char *token[D_MAX_TOKEN];

  int nentry, nindex, nlines, ntokens;

  long offset;

  struct idx *idx_ptr, *first_idx_ptr, *last_idx_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nindex = 0;
  offset = -1;

  /* Initialise idx pointers */
  idx_ptr = *fst_idx_ptr = first_idx_ptr = last_idx_ptr = NULL;

  /* Open the index file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      printf("*** Warning. Unable to open index file [%s]\n",file_name);
      return(nindex);
    }

  /* Read the file until we hit the right code */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL)
    {
      /* Split input line into tokens */
      ntokens = tokenize(input_line,token,D_MAX_TOKEN,' ');

      /* If have the right number of tokens, create an idx record */
      if (ntokens >= 2)
        {
          /* Get the offset and number of entries */
          offset = atol(token[1]);
          if (ntokens > 2)
            nentry = atoi(token[2]);
          else
            nentry = 0;

          /* Create an idx record */
          idx_ptr = create_idx_record(&first_idx_ptr,&last_idx_ptr,
                                      token[0],offset,nentry);

          /* Increment count of records created */
          nindex++;
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_idx_ptr = first_idx_ptr;

  /* Return the number of index entries */
  return(nindex);
}
/***********************************************************************

index_file  -  Read given file and store the offsets for its records by
               its key field

***********************************************************************/

long index_file(char *file_name,struct idx **fst_idx_ptr,
                char token_sep,int min_tokens,int key_token,
                char end_char)
{
  char input_line[D_LINELEN + 1];

  char *token[D_MAX_TOKEN], *string_ptr;

  int nentry, ntokens;

  long nindex, nlines, offset;

  struct idx *idx_ptr, *first_idx_ptr, *last_idx_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nentry = 0;
  nindex = nlines = 0;
  offset = -1;

  /* Initialise idx pointers */
  idx_ptr = *fst_idx_ptr = first_idx_ptr = last_idx_ptr = NULL;

  /* Open the file to be indexed */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      printf("*** Warning. Unable to open index file [%s]\n",file_name);
      return(nindex);
    }

  /* Get starting file position */
  offset = ftell(file_ptr);

  /* Read the file until we hit the right code */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL)
    {
      /* Increment line count */
      nlines++;
      
      /* Split input line into tokens */
      ntokens = get_tokens(input_line,token,D_MAX_TOKEN,token_sep);

      /* If have the right number of tokens, create an idx record */
      if (ntokens >= min_tokens)
        {
          /* Truncate the key token if required */
          if (end_char != '\0')
            {
              /* Truncate at the given character */
              string_ptr = strchr(token[key_token],end_char);
              if (string_ptr != NULL)
                string_ptr[0] = '\0';
            }

          /* Create an idx record */
          idx_ptr = create_idx_record(&first_idx_ptr,&last_idx_ptr,
                                      token[key_token],offset,nentry);

          /* Increment count of records created */
          nindex++;
        }

      /* Free tokens */
      free_tokens(token,ntokens);

      /* Get current file position */
      offset = ftell(file_ptr);
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return current pointers */
  *fst_idx_ptr = first_idx_ptr;

  /* Return the number of index entries */
  return(nindex);
}
/***********************************************************************

get_gnomad_offset  -  Read archuacc.idx or archuid.idx index file to
                      pick up the starting offset for the corresponding
                      .dat file

***********************************************************************/

long get_gnomad_offset(char *uniprot_acc,int *nent,
                       struct c_file_data file_name_ptr)
{
  char file_name[D_FILENAME_LEN], input_line[D_LINELEN + 1];
  char number_string[20], search_string[UNIPROT_LEN + 2];

  char *token[D_MAX_TOKEN];

  int len, nentry, nlines, ntokens, slen;
  int done;

  long offset;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  nentry = 0;
  offset = -1;

  /* Form the search string */
  strcpy(search_string,uniprot_acc);
  strcat(search_string," ");
  slen = strlen(search_string);

  /* Form the name of the index file */
  strcpy(file_name,file_name_ptr.other_dir[GNOMAD_DIR]);
  strcat(file_name,"/gnomAD.idx");

  /* Open the index file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      printf("*** Warning. Unable to open data file [%s]\n",file_name);
      return(offset);
    }

  /* Read the file until we hit the right code */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* If this is our code, get the details */
      if (!strncmp(input_line,search_string,slen))
        {
          /* Split input line into tokens */
          ntokens = tokenize(input_line,token,D_MAX_TOKEN,' ');

          /* If have the right number of tokens, retrieve them */
          if (ntokens == 3)
            {
              /* Get the offset and number of entries */
              offset = atol(token[1]);
              nentry = atoi(token[2]);
            }

          /* Set flag that our search is over */
          done = TRUE;
        }

      /* If this code is larger than ours, then we are done */
      if (strcmp(input_line,search_string) > 0)
        done = TRUE;
    }

  /* Return the number of entries */
  *nent = nentry;

  /* Return the starting offset */
  return(offset);
}
/***********************************************************************

create_gnomad_record  -  Create a new gnomad record

***********************************************************************/

struct gnomad *create_gnomad_record(struct gnomad **fst_gnomad_ptr,
                                    struct gnomad **lst_gnomad_ptr,
                                    char *token[D_MAX_TOKEN],int ntokens)
{
  int itoken;

  struct gnomad *gnomad_ptr, *first_gnomad_ptr, *last_gnomad_ptr;

  /* Initialise gnomad pointers */
  gnomad_ptr = NULL;
  last_gnomad_ptr = *lst_gnomad_ptr;
  first_gnomad_ptr = *fst_gnomad_ptr;

  /* Create gnomad entry */
  gnomad_ptr = (struct gnomad*)malloc(sizeof(struct gnomad));
  if (gnomad_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct gnomad\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      if (strlen(token[itoken]) > 0)
        store_string(&(gnomad_ptr->field[itoken]),token[itoken]);
      else
        gnomad_ptr->field[itoken] = &null;
    }

  /* Initialise pointer to next gnomad record */
  gnomad_ptr->next_gnomad_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_gnomad_ptr == NULL)
    first_gnomad_ptr = gnomad_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_gnomad_ptr->next_gnomad_ptr = gnomad_ptr;
                      
  /* Save the current gnomad's pointer */
  last_gnomad_ptr = gnomad_ptr;

  /* Return current pointers */
  *fst_gnomad_ptr = first_gnomad_ptr;
  *lst_gnomad_ptr = last_gnomad_ptr;
  return(gnomad_ptr);
}
/***********************************************************************

read_gnomad_dat  -  Read in the natural variants for the given UniProt
                    code from the gnomAD data file

***********************************************************************/

int read_gnomad_dat(char *uniprot_acc,long offset,int nentry,
                    struct gnomad **fst_gnomad_ptr,
                    struct c_file_data file_name_ptr)
{
  char file_name[D_FILENAME_LEN], input_line[D_LINELEN + 1];
  char aa_string[4], number_string[20];

  char *token[D_MAX_TOKEN], *split[D_MAX_TOKEN], *string_ptr;

  int isplit, itoken, len, line, ngnomad, nread, ntokens;
  int done, error, wanted;
  
  static int first = TRUE;

  struct gnomad *gnomad_ptr, *first_gnomad_ptr, *last_gnomad_ptr;

  static FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  line = 0;
  ngnomad = 0;
  nread = 0;

  /* Initialise pointers */
  gnomad_ptr = first_gnomad_ptr = last_gnomad_ptr = NULL;

  /* If this is the first call, open the gnomAD file */
  if (first == TRUE)
    {
      /* Form the name of the input file */
      strcpy(file_name,file_name_ptr.other_dir[GNOMAD_FILE]);

      /* Open the file */
      if ((file_ptr = fopen(file_name,"r")) == NULL)
        {
          /* Otherwise, error must be due to missing file */
          printf("*** Warning. Unable to open data file [%s]\n",file_name);
          return(ngnomad);
        }

      /* Unset flag */
      first = FALSE;
    }

  /* Find start-point in data file using fseek */
  error = fseek(file_ptr,offset,SEEK_SET);

  /* Read the file until we have read in all the data */
  while (fgets(input_line,D_LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate line and get its length */
      len = string_chop(input_line);
      line++;

      /* Tokenize this line */
      ntokens = tokenize(input_line,token,D_MAX_TOKEN,'\t');

      /* If have a UniProt id, then see if it is our one */
      if (ntokens > GNOMAD_UNIPROT &&
          strstr(token[GNOMAD_UNIPROT],uniprot_acc) != NULL)
        {
          /* Determine whether this variant is wanted */
          wanted = TRUE;
          if (!strcmp(token[GNOMAD_PROTEIN_POSITION],"NA"))
            wanted = FALSE;
          else if (strlen(token[GNOMAD_AMINO_ACIDS]) != 3)
            wanted = FALSE;
          else
            {
              /* Get the amino acid string */
              strcpy(aa_string,token[GNOMAD_AMINO_ACIDS]);

              /* Check for unwanted cases */
              if (aa_string[2] == '*' || aa_string[2] == '-' ||
                  aa_string[2] == 'X' || aa_string[0] == '-' ||
                  aa_string[1] != '/')
                 wanted = FALSE;
            }

          /* If entry wanted, then store */
          if (wanted == TRUE)
            {
              /* Create a gnomad record to store the details */
              gnomad_ptr
                = create_gnomad_record(&first_gnomad_ptr,
                                       &last_gnomad_ptr,token,ntokens);

              /* Increase count of entries stored */
              ngnomad++;
            }

          /* Increase count of entries read in */
          nread++;

          /* If we have read in all the entries we were expecting, then
             we're done */
          if (nread >= nentry)
            done = TRUE;
        }
    }

  /* Return pointer */
  *fst_gnomad_ptr = first_gnomad_ptr;

  /* Return the number of entries actually used */
  return(ngnomad);
}
/***********************************************************************

idx_sort  -  Shell-sort the idx entries according to contact percentage

***********************************************************************/

void idx_sort(struct idx **sort_index_ptr,int nidxid)
{
  int match;
  int swap;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct idx *idx1_ptr, *idx2_ptr, *swap_idx_ptr;

  /*--Initialise values */
  n = nidxid;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

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
              idx1_ptr = sort_index_ptr[indi];
              idx2_ptr = sort_index_ptr[indl];

              /* See if pointers need to be swapped */
              swap = FALSE;
              if (strcmp(idx1_ptr->key,idx2_ptr->key) > 0)
                swap = TRUE;
              else if (!strcmp(idx1_ptr->key,idx2_ptr->key) &&
                       idx1_ptr->nentry > idx2_ptr->nentry)
                swap = TRUE;

              /* If in wrong order, then swap sort-pointers */
              if (swap == TRUE)
                {
                  swap_idx_ptr = sort_index_ptr[indi];
                  sort_index_ptr[indi] = sort_index_ptr[indl];
                  sort_index_ptr[indl] = swap_idx_ptr;

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

create_idx_index  -  Create index for the idx records

***********************************************************************/

long create_idx_index(struct idx *first_idx_ptr,struct idx ***idx_ind_ptr,
                      long nindex)
{
  char last_key[D_LINELEN];
  
  long ientry, nentry;

  struct idx *idx_ptr;

  struct idx **idx_index_ptr;

  /* Create the index array as one large array */
  idx_index_ptr = (struct idx **) malloc(nindex * sizeof(struct idx *));

  /* Check that sufficient memory was allocated */
  if (idx_index_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for index array\n");
      exit(1);
    }

  /* Get pointer to the first record */
  idx_ptr = first_idx_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the sorted index
     array */
  while (idx_ptr != NULL)
    {
      /* Store current pointer at end of the index array */
      idx_index_ptr[nentry] = idx_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      idx_ptr = idx_ptr->next_idx_ptr;
    }

  /* Save number actually stored */
  nindex = nentry;
  nentry = 0;
  last_key[0] = '\0';

  /* Sort the idx ids */
  if (nindex > 1)
    idx_sort(idx_index_ptr,nindex);

  /* Pass through the index entries to number any duplicates */
  for (ientry = 0; ientry < nindex; ientry++)
    {
      /* Get this entry */
      idx_ptr = idx_index_ptr[ientry];

      /* If this is the same key value as before, increment entry count */
      if (!strcmp(idx_ptr->key,last_key))
        nentry++;

      /* Otherwise, re-initialise count */
      else
        {
          /* Reinitialise */
          nentry = 0;

          /* Save the key */
          strcpy(last_key,idx_ptr->key);
        }

      /* Save the entry number */
      idx_ptr->nentry = nentry;
    }

  /* Return the array pointer */
  *idx_ind_ptr = idx_index_ptr;

  /* Return number of idxs */
  return(nindex);
}
/***********************************************************************

index_files  -  Read through the UniProt and Ensembl data files to index
                them

***********************************************************************/

int index_files(struct idx ***idx_idx_ptr,long *nindex,
                struct c_file_data file_name_ptr)
{
  char file_name[D_FILENAME_LEN + 1];

  int ifile;
  
  struct idx *idx_ptr, *first_idx_ptr, *last_idx_ptr;

  struct idx **idx_index_ptr;

  /* Loop over the files to be indexed */
  for (ifile = 0; ifile < NINDEX_FILES; ifile++)
    {
      /* Initialise variables */
      nindex[ifile] = 0;

      /* Initialise pointers */
      idx_ptr = first_idx_ptr = last_idx_ptr = NULL;

      /* Form the name of the file to be indexed */
      strcpy(file_name,file_name_ptr.other_dir[COORUN_DIR]);
      strcat(file_name,"/");
      strcat(file_name,INDEX_FILE_NAME[ifile]);

      /* Call the indexing function */
      nindex[ifile]
        = index_file(file_name,&first_idx_ptr,INDEX_TOKEN_SEP[ifile],
                     INDEX_MIN_TOKENS[ifile],INDEX_KEY_TOKEN[ifile],
                     INDEX_END_CHAR[ifile]);

      /* Create a sorted index to allow binary search by the key field */
      nindex[ifile]
        = create_idx_index(first_idx_ptr,&idx_index_ptr,nindex[ifile]);

      /* Save the index */
      idx_idx_ptr[ifile] = idx_index_ptr;
    }
}
/***********************************************************************

exp2html  -  Convert exponential-format number to HTML notation

***********************************************************************/

void exp2html(char *number_string,char *html_string)
{
  char tmp[50];

  char *string_ptr;

  int exponent;

  double number;

  /* Get the number from the string */
  number = atof(number_string);

  /* Convert back into a string */
  sprintf(number_string,"%g",number);
  
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
}
/***********************************************************************

calc_variant_scores  -  Compute this variant's deleteriousness scores

***********************************************************************/

int calc_variant_scores(char *gnomad_sift,int *sft_band,float *sft_score,
                        char *gnomad_polyphen,int *polyph_band,
                        float *polyph_score)
{
  char tmp[D_LINELEN + 1];

  char *string_ptr;

  int iband;
  int deleterious, sift_band, polyphen_band;

  float sift_score, polyphen_score;

  /* Initialise variables */
  deleterious = 0;
  sift_band = polyphen_band = -1;
  sift_score = polyphen_score = 0;
  
  /* Extract the SIFT band for this variant */
  strcpy(tmp,gnomad_sift);
  sift_band = SIFT_UNKNOWN;
  if (!strcmp(gnomad_sift,"-"))
    sift_band = SIFT_NONE;
  else
    for (iband = 1; iband < NSIFT_BANDS; iband++)
      if (strstr(tmp,SIFT_BAND[iband]) != NULL)
        sift_band = iband;

  /* If just the score, extract it and calculate corresponding band */
  if (sift_band == SIFT_UNKNOWN)
    {
      /* Extract the score */
      sift_score = atof(gnomad_sift);

      /* See which band this falls in */
      for (iband = 1; iband < NSIFT_BANDS; iband++)
        if (sift_score >= SIFT_SCORE[iband])
          sift_band = iband;
    }

  /* Otherwise, separate consequence from the SIFT score */
  else if (sift_band != SIFT_NONE)
    {
      string_ptr = strchr(tmp,')');
      if (string_ptr != NULL)
        string_ptr[0] = '\0';
      string_ptr = strchr(tmp,'(');
      if (string_ptr != NULL)
        sift_score = atof(string_ptr + 1);
    }

  /* Repeat for POLYPHEN band */
  strcpy(tmp,gnomad_polyphen);
  polyphen_band = POLYPHEN_UNKNOWN;
  if (!strcmp(gnomad_polyphen,"-"))
    polyphen_band = POLYPHEN_NONE;
  else
    for (iband = 1; iband < NPOLYPHEN_BANDS; iband++)
      if (strstr(tmp,POLYPHEN_BAND[iband]) != NULL)
        polyphen_band = iband;

  /* If just the score, extract it and calculate corresponding band */
  if (polyphen_band == POLYPHEN_UNKNOWN)
    {
      /* Extract the score */
      polyphen_score = atof(gnomad_polyphen);

      /* See which band this falls in */
      for (iband = 1; iband < NPOLYPHEN_BANDS; iband++)
        if (polyphen_score <= POLYPHEN_VALUE[iband])
          polyphen_band = iband;
    }

  /* Otherwise, separate consequence from the POLYPHEN score */
  else if (polyphen_band != POLYPHEN_NONE)
    {
      /* Extract POLYPHEN score */
      string_ptr = strchr(tmp,')');
      if (string_ptr != NULL)
        string_ptr[0] = '\0';
      string_ptr = strchr(tmp,'(');
      if (string_ptr != NULL)
        polyphen_score = atof(string_ptr + 1);
    }

  /* Calculate a 'deleteriousness score' by combining the two */
  deleterious = (sift_band + 1) * (polyphen_band + 1);

  /* Return the various scores */
  *sft_band = sift_band;
  *sft_score = sift_score;
  *polyph_band = polyphen_band;
  *polyph_score = polyphen_score;
  
  /* Return the deleteriousness score */
  return(deleterious);
}
/***********************************************************************

get_colour_restype  -  Get colour from associated residue type

***********************************************************************/

void get_colour_restype(char *c_name,char aa_code)
{
  char *string_ptr;
  
  int iamino, itype, res_type;
  int done;

  /* Check for special characters */
  done = FALSE;
  for (itype = 0; itype < NMUT_TYPES && done == FALSE; itype++)
    {
      /* If this is the right one, colour red */
      if (aa_code == MUT_CODE[itype])
        {
          /* Set colour as red, and flag as done */
          strcpy(c_name,"red");
          done = TRUE;
        }
    }

  /* If we have our colour, then return */
  if (done == TRUE)
    return;

  /* Identify the amino acid */
  iamino = -1;
  string_ptr = strchr(AMINO,aa_code);
  if (string_ptr != NULL)
    iamino = string_ptr - AMINO;

  /* Get the corresponding residue type */
  res_type = RES_TYPE[iamino + 1];

  /* Get the corresponding colour */
  strcpy(c_name,RES_COLOUR_NAME[res_type]);
}
/***********************************************************************

get_data_dir  -  Get the directory where the data are

***********************************************************************/

void get_data_dir(char **ddir,char *root_dir,char *id)
{
  char file_name[C_FILENAME_LEN], parent_dir[C_FILENAME_LEN];
  char command[C_FILENAME_LEN], directory[C_FILENAME_LEN];

  char *data_dir;

  int len, loop;
  int error;

  struct stat s;

  /* Get the length of the UniProt account */
  len = strlen(id);

  /* Form the name of the parent directory */
  strcpy(file_name,root_dir);
  strcat(file_name,"/");
  strcat(file_name,id + len - 2);
  strcat(file_name,"/");
  strcpy(parent_dir,file_name);
  strcat(file_name,id);
  strcat(file_name,"/");
  store_string(&data_dir,file_name);

  /* On the first loop, check the parent directory */
  strcpy(directory,parent_dir);

  /* Loop twice to see if the directory, and its parent, need to be
     created */
  for (loop = 0; loop < 2; loop++)
    {
      /* Check if the directory exists */
      error = stat(directory,&s);

      /* If the directory does not exist, then create it */
      if(error == -1)
        {
          /* Create directory */
          strcpy(command,"mkdir ");
          strcat(command,directory);
#ifdef TEST
          error = FALSE;
#else
          error = system(command);
#endif
          if (error == TRUE)
            printf("*** Warning. Unable to create directory: %s\n",
                   directory);

          /* Change directory premission */
          strcpy(command,"chmod 777 ");
          strcat(command,directory);
#ifdef TEST
          error = FALSE;
#else
          system(command);
#endif
        }

      /* On the second loop, check the data directory */
      strcpy(directory,data_dir);
    }

  /* Return the data directory */
  *ddir = data_dir;
}
/***********************************************************************

get_domain_dir  -  Get the directory where the domain data are

***********************************************************************/

void get_domain_dir(char **ddir,char *root_dir,char *id)
{
  char command[C_FILENAME_LEN], directory[C_FILENAME_LEN];
  char parent_dir[4][C_FILENAME_LEN];
  char domain_type[10];

  char *data_dir, *string_ptr, *token[D_MAX_TOKEN];

  int have_file, itoken, len, loop, nparents, ntokens;
  int error;

  /* Initialise variables */
  for (itoken = 0; itoken < 4; itoken++)
    parent_dir[itoken][0] = '\0';
  directory[0] = '\0';
  
  /* Get the length of the UniProt account */
  len = strlen(id);

  /* Get the domain type */
  if (!strncmp(id,"CATH",4))
    {
      /* For CATH domain, have 4 parent directories */
      nparents = 4;

      /* Extract the CATH levels making up this code */
      ntokens = get_tokens(id,token,D_MAX_TOKEN,'.');

      /* If have enough tokens, form all the directory names, creating
         the directories where necessary */
      if (ntokens > 3)
        {
          /* Loop over the tokens */
          for (itoken = 0; itoken < 4; itoken++)
            {
              /* For first token, strip out the CATH_ prefix and form
                 directory name */
              if (itoken == 0)
                {
                  /* Form the name of the directory */
                  strcpy(parent_dir[0],root_dir);
                  strcat(parent_dir[0],"/CATH/");
                  strcat(parent_dir[0],token[0] + 5);
                }

              /* Otherwise, append level name to previous parent
                 directory */
              else
                {
                  /* If this is the last token, strip off the FunFam id */
                  if (itoken == 3)
                    {
                      /* Get start of FunFam id and remove */
                      string_ptr = strchr(token[itoken],'-');
                      if (string_ptr != NULL)
                        string_ptr[0] = '\0';
                    }
                  
                  /* Form this directory name */
                  strcpy(parent_dir[itoken],parent_dir[itoken - 1]);
                  strcat(parent_dir[itoken],"/");
                  strcat(parent_dir[itoken],token[itoken]);
                }
            }
        }

      /* Free the memory taken up by the tokens */
      free_tokens(token,ntokens);

      /* Form the last parent, and this domain's directory */
      strcpy(directory,parent_dir[3]);
      strcat(directory,"/");
      strcat(directory,id);
    }

  /* Otherwise, have a simpler directory structure */
  else
    {
      /* Only one parent directory for Pfam domains */
      nparents = 1;

      /* Form the name of the parent directory */
      strcpy(parent_dir[0],root_dir);
      strcat(parent_dir[0],"/Pfam/");
      strcat(parent_dir[0],id + len - 2);

      /* Form the name of this domain's directory */
      strcpy(directory,parent_dir[0]);
      strcat(directory,"/");
      strcat(directory,id);
    }

  /* Save the directory name */
  store_string(&data_dir,directory);

  /* Loop to check which directories exist, and which need to be created */
  for (loop = 0; loop < nparents + 1; loop++)
    {
      /* On last loop, check the directory itself */
      if (loop == nparents)
        strcpy(directory,data_dir);
      else
        strcpy(directory,parent_dir[loop]);

      /* See if the directory exists */
      strcpy(command,"ls -1d ");
      strcat(command,directory);
      strcat(command," >& /dev/null");
#ifdef TEST
      error = FALSE;
#else
      error = system(command);
#endif

      /* If it doesn't exist, create it */
      if (error != FALSE)
        {
          /* Create directory */
          strcpy(command,"mkdir ");
          strcat(command,directory);
#ifdef TEST
          error = FALSE;
#else
          error = system(command);
#endif
          if (error != FALSE)
            printf("*** Warning. Unable to create directory: %s\n",
                   directory);

          /* Change directory premission */
          strcpy(command,"chmod 755 ");
          strcat(command,directory);
#ifdef TEST
          error = FALSE;
#else
          system(command);
#endif
        }
    }

  /* Return the data directory */
  *ddir = data_dir;
}
/***********************************************************************

get_residue_name  -  Get the three-letter residue name from the
                     single-letter code

***********************************************************************/

void get_residue_name(char aa_code,char *res_name)
{
  int type;

  /* Initialise variables */
  res_name[0] = '\0';

  /* Check for one of the special codes */
  for (type = 0; type < NMUT_TYPES && res_name[0] == '\0'; type++)
    {
      if (aa_code == MUT_CODE[type])
        strcpy(res_name,MUT_TYPE_DESC_SHORT[type]);
    }

  /* If not found, assume standard amino acid */
  if (res_name[0] == '\0')
    get_res_name(aa_code,res_name);
}


