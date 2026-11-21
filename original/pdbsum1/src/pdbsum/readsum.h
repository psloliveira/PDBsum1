/***********************************************************************

  readsum.h - Include file for readsum.c

***********************************************************************/

/* Array parameters */
#define S_COLNAM_LEN        30
#define S_FILENAME_LEN    1024
#define S_LINELEN          512
#define S_MAX_TOKEN         10
#define S_STRING_LEN       (10 * S_LINELEN)
#define S_TOKEN_LEN         20

/* Template characters (allowing up to 6-character PDB code) */
static char S_TEMPLATE_CHARS[] = "abcdef";

/* Input file types */
#define S_GROW                0
#define S_IGROW               1

/* Interaction types */
#define S_HBONDS              0
#define S_CONTACTS            1
#define S_SALT_BRIDGES        2
#define S_DISULPHIDES         3
#define S_COVALENT            4

#define S_NTYPES              5

static char S_INT_TYPE[S_NTYPES] = { 'H', 'N', 'S', 'D', 'C' };

/* Record types defined in the rasmol.dfn file */
#define S_UNKNOWN            -1
#define S_CALPHA_ONLY         0
#define S_LIGAND              1
#define S_METAL               2
#define S_NUCLEIC_ACID        3
#define S_PROTEIN             4
#define S_LIGAND_RANGE        5

/* Field types */
#define S_UNWANTED           -1
#define S_TITLE               0
#define S_STRUCTURE           1
#define S_OBSOLETE            2
#define S_DATE                3
#define S_RELEASE_DATE        4
#define S_XRAY                5
#define S_RESOLUTION          6
#define S_NMR                 7
#define S_NMODELS             8
#define S_TH_MODEL            9

#define S_NFIELDS            10

static char *S_FIELD_NAME[S_NFIELDS] = {
  "TITLE ", "STRUCTURE ", "OBSOLETE TRUE", "DATE ", "RELEASE_DATE ",
  "XRAY TRUE", "RESOLUTION ", "NMR TRUE", "NMODELS ",
  "THEORETICAL_MODEL TRUE",
};

/* Structure for each grow data item
   --------------------------------- */
struct s_grow
{
  int                type;
  char               atom_name[2][5];
  char               res_name[2][4];
  char               res_num[2][6];
  char               chain[2];
  char               res_type[2];
  float              distance;
  char               *sort_string;
  struct s_grow      *next_s_grow_ptr;
};

