/***********************************************************************

  genpymol.h - Include file for genpymol.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN          12
#define COLNO_LEN           15
#define FILENAME_LEN       512
#define LINELEN          50000
#define LONG_LINELEN     40960
#define RASMOL_LINELEN     256
#define MAX_CHAINS         100
#define MAX_DOM           1000
#define MAX_LIGRANGE       220
#define MAX_TOKEN          100
#define MAX_TOKENS       10000
#define NAMINO              20
#define NAMLEN             512
#define OBJ_NAME_LEN        41
#define STRING_LEN         512
#define TOKEN_LEN           20

/* Maximum and minimum subscripts */
#define MIN                   0
#define MAX                   1

/* Additional object types */
#define DISULPH               8
#define ACTSITE               9

/* Object types */
#define CHAIN                 0
#define LIGMET                1                  

/* Transparency levels */
#define NONE                  0
#define PARTIAL               1                  
#define FULL                  2                  

/* Pfam domain types */
#define UNKNOWN              -1
#define PFAMA_DOMAIN          0
#define PFAMB_DOMAIN          1

#define NDOMAIN_TYPES         2

/* Degrees of transparency in PyMol */
static const float PYMOL_TRANSPARENCY[3] = { 0.0, 0.3, 0.6 };

/* Output types */
#define PYMOL_SCRIPT          1
#define RASMOL_SCRIPT         2
#define JMOL_SCRIPT           3
#define JSMOL                 4
#define ATOM_COORDS           5

/* Data in pdb.align file */
#define ALIGN_SEQ_NO         0
#define ALIGN_AA_CODE        1
#define ALIGN_PDB_CODE       2
#define ALIGN_CHAIN          3
#define ALIGN_ENTRY_NO       4
#define ALIGN_DOMAIN         5
#define ALIGN_ISEQ           6
#define ALIGN_PDB_AA         7
#define ALIGN_RES_NAME       8
#define ALIGN_RES_NUM        9
#define ALIGN_START         10
#define ALIGN_END           11

#define NALIGN_FIELDS       12

/* Cut-off distance for structure surrounding chosen residue */
#define CUTOFF_DIST   ((float) 15.0)
#define CUTOFF_DIST_SQRD   (CUTOFF_DIST * CUTOFF_DIST)

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  char                    uniprot_acc[20];
  char                    chain;
  int                     pdb_type;
  int                     disastr;
  int                     user_id;
  int                     output_type;
  int                     labels;
  int                     neighbours;
  int                     start;
  int                     end;
  int                     variants;
  int                     run_64;
  int                     covid19;
  int                     start_colour;
  int                     colour_by_domain;
  double                  cutoff_dist;
  int                     naddmut;
  int                     res_cons;
  char                    *pml_file;
  struct addmut           *first_addmut_ptr;
};

/* Structure for each addmut data item
   ----------------------------------- */
struct addmut
{
  int                seq_no;
  char               aa_code;
  char               aa_mut;
  char               *added_variant;
  char               *aa_extra;
  char               *mut_extra;
  int                mut_length;
  int                type;
  struct addmut      *next_addmut_ptr;
};

/* Structure for each residue data item
   ------------------------------------ */
struct residue
{
  int                seq_no;
  char               aa_code;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  char               pdb_res_name[4];
  char               pdb_res_num[6];
  char               pdb_code[5];
  char               pdb_chain;
  char               sec_struc;
  float              conservation;
  char               last_atom_name[5];
  struct r_residue   *r_residue_ptr;
  struct alist       *first_alist_ptr;
  struct alist       *last_alist_ptr;
  struct residue     *next_residue_ptr;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  char               copy_of;
  struct chain       *next_chain_ptr;
};

/* Structure for annot item
   ------------------------- */
struct annot
{
  int                     annot_no;
  int                     start_res;
  int                     end_res;
  int                     length;
  int                     pre_start;
  int                     post_end;
  char                    *sequence;
  int                     domain;
  int                     deleted;
  int                     nrefs;
  int                     res_match;
  int                     serial_no;
  struct aamut            *aamut_ptr;
  struct desc             *desc_ptr;
  struct annot            *next_annot_ptr;
};

/* Structure for each aamut data item
   ---------------------------------- */
struct aamut
{
  int                seq_no;
  char               aa_code;
  char               aa_mut;
  int                ndiseases;
  char               *aa_extra;
  char               *mut_extra;
  int                mut_length;
  int                nwritten;
  int                var_type;
  struct annot       *annot_ptr;
  struct disease     *disease_ptr;
  struct aamut       *next_aamut_ptr;
};

/* Structure for annot name
   ------------------------- */
struct desc
{
  char                    *annot_id;
  char                    *name;
  char                    *description;
  int                     type;
  int                     colour_num;
  int                     nannot;
  struct desc             *next_desc_ptr;
};

/* Structure for storing disease details
   ------------------------------------- */
struct disease
{
  char               *id;
  int                type;
  char               *name;
  char               *description;
  char               *omim_id;
  int                nrefs;
  int                ncodes;
  int                nmut;
  int                is_note;
  struct disease     *next_disease_ptr;
};

/* Structure for list of residue annotations
   ----------------------------------------- */
struct alist
{
  struct annot            *annot_ptr;
  struct alist            *next_alist_ptr;
};

/* Structure for each pymol object item
   ------------------------------------ */
struct pobject
{
  char                    *name;
  char                    *display_type;
  int                     type;
  char                    chain_id;
  int                     ligmet_type;
  int                     got_sites;
  int                     got_disulph;
  struct pobject          *next_pobject_ptr;
};

/* Structure for domain item
   ------------------------- */
struct domain
{
  char                    *pfam_id;
  char                    *name;
  char                    *description;
  int                     type;
  int                     colour;
  int                     start;
  int                     end;
  struct r_residue        *start_r_residue_ptr;
  struct r_residue        *end_r_residue_ptr;
  struct domain           *next_domain_ptr;
};

