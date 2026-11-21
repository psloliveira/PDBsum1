/***********************************************************************

  showmain.h - Include file for showmain.c

***********************************************************************/

/* Flag-values */
#define START              0
#define END                1

/* Array parameters */
#define CODE_LEN          15
#define EC_LEN            15
#define FILENAME_LEN     512
#define LINELEN       100240
#define MAX_HTML_LINES  8000
#define MAX_LEVEL          5
#define MAX_TAG_LEN       20
#define MAXTEMPLATE_LEN   80
#define MAX_TOKEN         80
#define NAMLEN           100
#define NCOUNTERS         10
#define NPAR              10
#define RES_PER_LINE      50
#define TMP_STRING_LEN  (LINELEN)
#define TOKEN_LEN         30

FILE *dump_file_ptr;

/* Dummy loop array */
static int LOOP[MAX_LEVEL] = { -1, -1, -1, -1, -1 };

/* Number of Contents entries */
#define NCONTENT_HEADERS          6

/* Contents keywords in the database.dat file */
static char *CONTENTS_HEADER[NCONTENT_HEADERS] =
  {
    "PDB_CODE",
    "PROTEIN",
    "DNA",
    "LIGAND",
    "METAL",
    "WATER"
  };

/* Number of headers relating to tabs at top of page */
#define NTAB_HEADERS              4

/* Tab keywords in the database.dat and ligands.dat files */
static char *TAB_HEADER[NTAB_HEADERS] =
  {
    "PROTEIN",
    "DNA",
    "NLIGAND_TYPES",
    "LINKS"
  };

/* Number of structure entries */
#define NSTRUCTURE_HEADERS        6

/* Structure keywords in the database.dat file */
static char *STRUCTURE_HEADER[NSTRUCTURE_HEADERS] =
  {
    "WIRPLOT",
    "PROMOTIF",
    "CSA",
    "SITE",
    "VARIANTS",
    "PROSITE"
  };

/* Number of annotation tags used in the template file */
#define NTAGS                   4

/* Tag types */
#define DATA_TAG                0
#define LNK_TAG                 1

/* Tags used in template file */
static char *TAG_END = "%>";
static char *TAG[NTAGS][2] = {
  { "<NO%>", "</NO%>" },
  { "<FLD%>", "</FLD%>" },
  { "<LNK%>", "</LNK%>" },
  { "<EXT%>", "</EXT%>" }
};

/* States used in formatting the formula */
#define BLANK             -1
#define ELEMENT            0
#define NUMBER             1
#define CHARGE             2
#define UPPER_CASE         3
#define LOWER_CASE         4

/* Special files */
#define STANDARD_PDB          0
#define UPLOAD_PDB            1
#define MCSG_PDB              2
#define PROFUNC               3
#define EBI_COLOUR            4
#define DRAT_DB               5
#define SAS                   6
#define DRUGPORT              7
#define LIGSEARCH             8
#define LIGPLUS               9
#define VARSITE              10
#define ALPHA_FOLD           11

#define NSPECIAL             12

/* PDBsum bar colours */
static char *BAR_COLOUR[NSPECIAL]
= { "#207a7a", "#ff0000", "#a7c6e5", "#ffc125", "#207a7a", "#207a7a",
"#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff", "#0000ff" };

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  int                     version;
  char                    pdb_code[CODE_LEN];
  char                    chain;
  char                    template_name[MAXTEMPLATE_LEN + 1];
  char                    template_dir[FILENAME_LEN];
  char                    template_real_dir[FILENAME_LEN];
  char                    database_dat[FILENAME_LEN];
  int                     enzymes;
  int                     drat;
  int                     drugport;
  int                     profunc;
  int                     ligsearch;
  int                     ligplus;
  int                     varsite;
  int                     sas;
  int                     highlights;
  int                     bottom_level;
  int                     nlevels;
  int                     instance[MAX_LEVEL];
  int                     ec_number[MAX_LEVEL];
  char                    ec_code[TOKEN_LEN];
  char                    output_file[FILENAME_LEN];
  char                    option[NAMLEN];
  char                    include_file[FILENAME_LEN];
  char                    dir_name[FILENAME_LEN];
  char                    subdir[FILENAME_LEN];
  char                    subdir_real[FILENAME_LEN];
  char                    subdir_dir[FILENAME_LEN];
  char                    work_dir[FILENAME_LEN];
  char                    parameter[NPAR][FILENAME_LEN];
  int                     fullname;
  int                     call_no;
  long                    start_no;
  int                     have_error;
  char                    *user_id;
  char                    security_code[20];
  int                     pdb_type;
  int                     relinks;
  int                     fake_entry;
  char                    score[20];
  int                     from_ebi;
  char                    drug_id[20];
  char                    uniprot_acc[20];
  char                    *ligand_ranges;
  char                    *link_path;
  int                     proprietary;
  int                     list_fields;
};

/* Structure for PDB entry data
   ---------------------------- */
struct field
{
  char                    *field_id;
  int                     read;
  int                     nlevel;
  char                    *value;
  struct field            *left_field_ptr;
  struct field            *right_field_ptr;
};

/* Structure for each iname record
   ------------------------------- */
struct iname
{
  char                   *image_name;
  struct iname           *next_iname_ptr;
};

