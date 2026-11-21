/***********************************************************************

  align.c - Sequence alignment routines

***********************************************************************

Date:-         20 feb 2015

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Compilation
   -----------

cc -c align.c

*/

/* Routines
   --------

   -> align_sequences
         -> create_alignment_arrays
         -> fill_array
         -> nw_sweep
         -> path_len
         -> perform_alignment

*/

#include "common.h"

/***********************************************************************

create_alignment_arrays  -  Create the arrays for the sequence alignment

***********************************************************************/

void create_alignment_arrays(int *seq_len,int **arr,int **tmparr,
                             int **pth)
{
  int iresid, iseq, max_path, narray;

  int *array, *path, *tmp_array;

  /* Initialise variables */
  narray = seq_len[0] * seq_len[1];
  max_path = seq_len[0] + seq_len[1];

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

  /* Create working array for use in Needleman-Wunsch process */
  tmp_array = (int *) malloc(narray * sizeof(int));
  if (tmp_array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for temporary array\n");
      exit(1);
    }

  /* Return the array pointers */
  *arr = array;
  *pth = path;
  *tmparr = tmp_array;
}
/***********************************************************************

fill_array  -  Initialise array with standard residue similarities

***********************************************************************/

void fill_array(int *seq_len,char *sequence1,char *sequence2,
                int *array,int *tmp_array)
{
#define MIN_VALUE    0
#define NAMINO      20

  char *string_ptr;

  int ipoint, ires, jres, res_type1, res_type2;

  static char amino[] = "XARNDCQEGHILKMFPSTWYV";

  static int matrix[NAMINO + 1][NAMINO + 1] =
    {
  /*  X  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* X */
      0, 2,-2, 0, 0,-2, 0, 0, 1,-1,-1,-2,-1,-1,-4, 1, 1, 1,-6,-3, 0, /* A */
      0,-2, 6, 0,-1,-4, 1,-1,-3, 2,-2,-3, 3, 0,-4, 0, 0,-1, 2,-4,-2, /* R */
      0, 0, 0, 2, 2,-4, 1, 1, 0, 2,-2,-3, 1,-2,-4,-1, 1, 0,-4,-2,-2, /* N */
      0, 0,-1, 2, 4,-5, 2, 3, 1, 1,-2,-4, 0,-3,-6,-1, 0, 0,-7,-4,-2, /* D */
      0,-2,-4,-4,-5,12,-5,-5,-3,-3,-2,-6,-5,-5,-4,-3, 0,-2,-8, 0,-2, /* C */
      0, 0, 1, 1, 2,-5, 4, 2,-1, 3,-2,-2, 1,-1,-5, 0,-1,-1,-5,-4,-2, /* Q */
      0, 0,-1, 1, 3,-5, 2, 4, 0, 1,-2,-3, 0,-2,-5,-1, 0, 0,-7,-4,-2, /* E */
      0, 1,-3, 0, 1,-3,-1, 0, 5,-2,-3,-4,-2,-3,-5,-1, 1, 0,-7,-5,-1, /* G */
      0,-1, 2, 2, 1,-3, 3, 1,-2, 6,-2,-2, 0,-2,-2, 0,-1,-1,-3, 0,-2, /* H */
      0,-1,-2,-2,-2,-2,-2,-2,-3,-2, 5, 2,-2, 2, 1,-2,-1, 0,-5,-1, 4, /* I */
      0,-2,-3,-3,-4,-6,-2,-3,-4,-2, 2, 6,-3, 4, 2,-3,-3,-2,-2,-1, 2, /* L */
      0,-1, 3, 1, 0,-5, 1, 0,-2, 0,-2,-3, 5, 0,-5,-1, 0, 0,-3,-4,-2, /* K */
      0,-1, 0,-2,-3,-5,-1,-2,-3,-2, 2, 4, 0, 6, 0,-2,-2,-1,-4,-2, 2, /* M */
      0,-4,-4,-4,-6,-4,-5,-5,-5,-2, 1, 2,-5, 0, 9,-5,-3,-3, 0, 7,-1, /* F */
      0, 1, 0,-1,-1,-3, 0,-1,-1, 0,-2,-3,-1,-2,-5, 6, 1, 0,-6,-5,-1, /* P */
      0, 1, 0, 1, 0, 0,-1, 0, 1,-1,-1,-3, 0,-2,-3, 1, 2, 1,-2,-3,-1, /* S */
      0, 1,-1, 0, 0,-2,-1, 0, 0,-1, 0,-2, 0,-1,-3, 0, 1, 3,-5,-3, 0, /* T */
      0,-6, 2,-4,-7,-8,-5,-7,-7,-3,-5,-2,-3,-4, 0,-6,-2,-5,17, 0,-6, /* W */
      0,-3,-4,-2,-4, 0,-4,-4,-5, 0,-1,-1,-4,-2, 7,-5,-3,-3, 0,10,-2, /* Y */
      0, 0,-2,-2,-2,-2,-2,-2,-1,-2, 4, 2,-2, 2,-1,-1,-1, 0,-6,-2, 4  /* V */
    };


  /* Loop over the residues in first sequence */
  for (ires = 0; ires < seq_len[0]; ires++)
    {
      /* Get residue type */
      string_ptr = strchr(amino,sequence1[ires]);
      if (string_ptr != NULL)
        res_type1 = string_ptr - amino;
      else
        res_type1 = 0;

      /* Loop over the residues in second sequence */
      for (jres = 0; jres < seq_len[1]; jres++)
        {
          /* Get residue type */
          string_ptr = strchr(amino,sequence2[jres]);
          if (string_ptr != NULL)
            res_type2 = string_ptr - amino;
          else
            res_type2 = 0;

          /* Get the array location of this residue pair */
          ipoint = ires * seq_len[1] + jres;

          /* Get the similarity between these two residues */
          array[ipoint] = matrix[res_type1][res_type2] - MIN_VALUE;
          tmp_array[ipoint] = array[ipoint];
        }
    }
}
/***********************************************************************

nw_sweep  -  Sweep through the alignment matrix, updating the cells
             according to the best path from any cell on the previous
             row or column - slow version

***********************************************************************/

void nw_sweep(int *array,int seq_len[2],int diff)
{
#define GAP_PENALTY           5

  int ipoint, ires, jend, jpoint, jres, jstart, kres, lres;
  int max_score, score;

  float gradient;

  /* Calculate gradient of diagonal */
  gradient = (float) seq_len[1] / (float) seq_len[0];

  /* Loop backwards through the array, starting at the bottom right
     cell and working leftwards and upwards */
  for (ires = seq_len[0] - 2; ires > -1; ires--)
    { 
      /* Calculate the range of j-values so that we only traverse
         close to the diagonal */
      jstart = gradient * ires + diff;
      if (jstart > seq_len[1] - 2)
        jstart = seq_len[1] - 2;
      if (jstart < 0)
        jstart = 0;
      jend = gradient * ires - diff;
      if (jend < -1)
        jend = -1;
      if (jend > seq_len[1] - 2)
        jend = seq_len[1] - 2;
      if (jstart < jend)
        jstart = jend;

      /* Loop over the given range of j-values */
      for (jres = jstart; jres > jend; jres--)
        {
          /* Get the array location of this residue pair */
          ipoint = ires * seq_len[1] + jres;

          /* Initialise maximum score */
          max_score = -100;

          /* Loop over the cells to the right in the row below the
             current one */
          for (kres = ires + 1, lres = jres + 1;
               kres < seq_len[0] && lres < seq_len[1]; lres++)
            {
              /* Get this cell's array location */
              jpoint = kres * seq_len[1] + lres;

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
               kres < seq_len[0] && lres < seq_len[1]; kres++)
            {
              /* Get this cell's array location */
              jpoint = kres * seq_len[1] + lres;

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

int find_best_path(int *array,int seq_len[2],int *path,int *alnlen)
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
  max_path = seq_len[0] + seq_len[1];

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
      for (jres = start_jres; jres < seq_len[1]; jres++)
        {
          /* Get the array location of this cell */
          ipoint = ires * seq_len[1] + jres;
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
      for (ires = start_ires; ires < seq_len[0]; ires++)
        {
          /* Get the array location of this cell */
          ipoint = ires * seq_len[1] + jres;
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
          ipoint = max_pos[0] * seq_len[1] + max_pos[1];

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
      if (start_ires > seq_len[0] - 1 || start_jres > seq_len[1] - 1)
        {
          /* If haven't finished on the very last cell of the array,
             make path end one cell beyond it */
          if (start_ires != seq_len[0] || start_jres != seq_len[1])
            {
              /* Check for gap in first sequence */
              if (seq_len[0] - last_pos[0] > 1)
                {
                  /* Add insertion in first sequence to overall alignment
                     length */
                  align_len = align_len + seq_len[0] - last_pos[0] - 1;
                }

              /* Check for gap in second sequence */
              else if (seq_len[1] - last_pos[1] > 1)
                {
                  /* Add insertion in second sequence to overall alignment
                     length */
                  align_len = align_len + seq_len[1] - last_pos[1] - 1;
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

perform_alignment  -  Get the sequence alignment and calculate the sequence
                      similarity

***********************************************************************/

float perform_alignment(int *array,int *path,int path_len,char *sequence1,
                        char *sequence2,int *seq_len,int align_len,
                        char **alignment,int *nident)
{
  char ch, ch1, ch2;

  int aline, last_ires, last_jres, ipoint, ipos, ires, jres, kres, lpos;
  int align_pos, nidentical, noverlap, len_short, start_pos;
  int done, overlap;

  double percent_id, percent_ol;

  /* Initialise variables */
  align_pos = 0;
  last_ires = -1;
  last_jres = -1;
  overlap = FALSE;
  noverlap = 0;
  nidentical = 0;

  /* Get the length of the shortest sequence */
  len_short = seq_len[0];
  if (seq_len[1] < len_short)
    len_short = seq_len[1];
  if (len_short == 0)
    len_short = 1;

  /* Loop over the points along the path */
  for (ipos = 0; ipos < path_len; ipos++)
    {
      /* Get this point in the path */
      ipoint = path[ipos];

      /* Extract row and column number of this cell */
      if (ipoint == -1)
        {
          ires = seq_len[0];
          jres = seq_len[1];
          overlap = FALSE;
        }
      else
        {
          ires = ipoint / seq_len[1];
          jres = ipoint - ires * seq_len[1];
        }

      /* Check for gap in first sequence */
      if (ires - last_ires > 1)
        {
          /* Write out insertion from first sequence */
          for (kres = last_ires + 1; kres < ires; kres++)
            {
              /* Retrieve the relevant residue */
              ch = sequence1[kres];

              /* Store in alignment array */
              if (align_pos < align_len)
                {
                  /* Store current alignment position */
                  alignment[0][align_pos]= ch;
                  alignment[1][align_pos] = '-';
                  align_pos++;
                }

              /* Increment overlap length */
              if (overlap == TRUE)
                noverlap++;
            }
        }

      /* Check for gap in second sequence */
      else if (jres - last_jres > 1)
        {
          /* Write out insertion from second sequence */
          for (kres = last_jres + 1; kres < jres; kres++)
            {
              /* Retrieve the relevant residue */
              ch = sequence2[kres];

              /* Store in alignment array */
              if (align_pos < align_len)
                {
                  /* Store current alignment position */
                  alignment[1][align_pos] = ch;
                  alignment[0][align_pos] = '-';
                  align_pos++;
                }

              /* Increment overlap length */
              if (overlap == TRUE)
                noverlap++;
            }
        }

      /* If valid, write out the current alignment position */
      if (ires < seq_len[0] && jres < seq_len[1])
        {
          /* Retrieve the aligned residue from the first sequence */
          ch1 = sequence1[ires];

          /* Retrieve the aligned residue from the second sequence */
          ch2 = sequence2[jres];

          /* Store in alignment array */
          if (align_pos < align_len)
            {
              /* Store the a.a. codes of each sequence */
              alignment[0][align_pos] = ch1;
              alignment[1][align_pos] = ch2;
              align_pos++;
            }

          /* If residues identical, increment count of sequence
             identities */
          if (ch1 == ch2)
            nidentical++;

          /* Increment overlap length */
          overlap = TRUE;
          noverlap++;
        }

      /* Store current position */
      last_ires = ires;
      last_jres = jres;
    }

  /* Terminate the alignment strings */
  alignment[0][align_len] = '\0';
  alignment[1][align_len] = '\0';

  /* Calculate percentage identity */
  percent_id = 100.0 * (float) nidentical / (float) len_short;
  if (noverlap > 0)
    percent_ol = 100.0 * (float) nidentical / (float) noverlap;
  else
    percent_ol = 0.0;

  /* Return number of identical residues */
  *nident = nidentical;

  /* Return the percentage similarity */
  return(percent_id);
}
/***********************************************************************

align_sequences  -  Align the query sequence with the current PDB
                    sequence, heavily weighting the equivalent
                    ligand-binding residues

***********************************************************************/

float align_sequences(char *sequence1,char *sequence2,
                      char *alignment[2],int *aln_len,int *nident)
{
  int align_len, i, ipos, istart, nidentical, seqlen1, lenkey;
  int diff, dir_type, max_len, path_len, seqlen2, seq_len[2];
  int found;

  int *array, *path, *tmp_array;

  float similarity;

  /* Initialise variables */
  seqlen1 = strlen(sequence1);
  seqlen2 = strlen(sequence2);
  similarity = 0.0;

  /* Get the difference between the two sequence lengths */
  diff = abs(seqlen1 - seqlen2);
  max_len = seq_len[0] = seqlen1;
  seq_len[1] = seqlen2;
  if (seq_len[1] > max_len)
    max_len = seq_len[1];
  diff = max_len / 2;

  /* Create the alignment arrays */
  create_alignment_arrays(seq_len,&array,&tmp_array,&path);

  /* Initialise array with standard residue similarities */
  fill_array(seq_len,sequence1,sequence2,array,tmp_array);

  /* Apply Needleman and Wunsch algorithm to adjust scores matrix */
  nw_sweep(array,seq_len,diff);

  /* Find best route through matrix to give best alignment of the two
     sequences */
  path_len = find_best_path(array,seq_len,path,&align_len);

  /* Create alignment arrays */
  for (i = 0; i < 2; i++)
    {
      alignment[i] = (char *) malloc((align_len + 1) * sizeof(char));
      if (alignment[i]  == NULL)
        {
          printf("*** ERROR. Unable to allocate memory for alignment"
                 " array\n");
          exit(1);
        }
    }

  /* Calculate the sequence similarity */
  similarity
    = perform_alignment(tmp_array,path,path_len,sequence1,
                        sequence2,seq_len,align_len,alignment,
                        &nidentical);

  /* Free up memory taken by the alignment arrays */
  free(array);
  free(path);
  free(tmp_array);

  /* Return the alignment length and number of identical residues */
  *aln_len = align_len;
  *nident = nidentical;

  /* Return the sequence similarity */
  return(similarity);
}
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int test_main(int argc,char *argv[])
{
  static char *sequence1
    = "SNISRQAYADMFGPTVGDKVRLADTELWIEVEDDLTTYGEEVKFGGGKVIRDGMGQGQMLAADCVDLVLTNALIVDHWGIVKADIGVKDGRIFAIGKAGNPDIQPNVTIPIGAATEVIAAEGKIVTAGGIDTHIHWICPQQAEEALVSGVTTMVGGGTGPAAGTHATTCTPGPWYISRMLQAADSLPVNIGLLGKGNVSQPDALREQVAAGVIGLXIHEDWGATPAAIDCALTVADEMDIQVALHSDTLNESGFVEDTLAAIGGRTIHTFHTEGAGGGHAPDIITACAHPNILPSSTNPTLPYTLNTIDEHLDMLMVYHHLDPDIAEDVAFAESRIRRETIAAEDVLHDLGAFSLTSSDSQAMGRVGEVILRTWQVAHRMKVQRGALAEETGDNDNFRVKRYIAKYTINPALTHGIAHEVGSIEVGKLADLVVWSPAFFGVKPATVIKGGMIAIAPMGDINASIPTPQPVHYRPMFGALGSARHHCRLTFLSQAAAANGVAERLNLRSAIAVVKGCRTVQKADMVHNSLQPNITVDAQTYEVRVDGELITSEPADVLPMAQRYFLF";
  static char *sequence2
    = "SNISRQAYADMFGPTVGDKVRLADTELWIEVEDDLTTYGEEVKFGGGKVIRDGMGQGQMLAADCVDLVLTNALIVDHWGIVKADIGVKDGRIFAIGKAGNPDIQPNVTIPIGAATEVIAAEGKIVTASIEVGKLADLVVWSPAFFGVKPATVIKGGMIAIAPMGDINASIPTPQPVHYRPMFGAL";

  char *alignment[2];

  int align_len, nidentical;

  float identity;

  /* Call the alignment routines */
  identity = align_sequences(sequence1,sequence2,alignment,&align_len,
                             &nidentical);

  /* Show the percentage identity */
  printf("Number of identical residues = %8d\n",nidentical);
  printf("Percentage identity          = %8.3f\n",identity);
  printf("Alignment length             = %8d\n",align_len);
  printf("%s\n",alignment[0]);
  printf("%s\n",alignment[1]);



  return(0);
}
