/***********************************************************************

  fixpdb.h - Include file for fixpdb.c

***********************************************************************/

/* Array parameters */
#define FILENAME_LEN  512
#define LINELEN      1000

/* Chain types */
#define UNKNOWN              'U'
#define PROTEIN              'P'
#define NUCLEIC_ACID         'N'
#define DNA                  'D'
#define RNA                  'R'

#define NAMINO               20

static char *BASE_NAMES
= "ADE|CYT|GUA|THY|URA|A  |C  |G  |T  |U  |  A|  C|  G|  T|  U|"
  " DA| DC| DG| DT|  U";

static char *DNA_NAME[5] = { " DA", " DC", " DG", " DT", " DU" };
static char *RNA_NAME[5] = { "  A", "  C", "  G", "  T", "  U" };

#define NWATER_NAMES         5

static char *WATER_NAME[NWATER_NAMES]
= { "WAT", "0W0", "OWO", "SOL", "H2O" };

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    *file_name;
  int                     run_64;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  struct chain       *next_chain_ptr;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  char               aa_code;
  char               res_name[4];
  char               res_num[6];
  struct atom        *first_atom_ptr;
  struct chain       *chain_ptr;
  struct residue     *next_residue_ptr;
};

/* Structure for each atom record
   ------------------------------ */
struct atom
{
  int                     atom_number;
  char                    atom_name[5];
  double                  coords[3];
  double                  bvalue;
  double                  occupancy;
  struct atom             *next_atom_ptr;
};

