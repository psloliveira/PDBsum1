/***********************************************************************

  resdata.h - Include file for resdata.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN          11
#define FILENAME_LEN       512
#define LINELEN           1024
#define MAXDOMCHAIN        200
#define MAX_INT_CHAINS      10
#define MAX_PATTERNS         5
#define MAX_SITES            5
#define MAX_TOKEN          100
#define PROSITE_PATLEN      15
#define STRING_LEN        1000
#define TOKEN_LEN          100

/* Record types defined in the rasmol.dfn file */
#define DELETED              -9
#define OTHER                -1
#define CA_ONLY               0
#define LIGAND                1
#define METAL                 2
#define NUCLEIC_ACID          3
#define PROTEIN               4
#define LIGAND_RANGE          5

/* Contact types */
#define CONTACT_TYPES         3
#define NOT_WANTED           -9
#define TO_DNA                0
#define TO_LIGAND             1
#define TO_METAL              2

/* Secondary structure types */
#define COIL                  0
#define HELIX                 1
#define STRAND                2

/* Secondary structure symbols */
static char SEC_STRUC[] = ".HE";

/* States during reading through the domain line */
#define WOLF_CODE         0
#define NDOMAINS          1
#define NFRAGS            2
#define NSEGMENTS         3
#define CHAIN1            4
#define RES1              5
#define HYPHEN1           6
#define CHAIN2            7
#define RES2              8
#define HYPHEN2           9
#define DONE_SEG         10

/* SITE records */
#define SITE_DESCRIPTION  0
#define RESIDUE_DETAILS   1

FILE *dump_file_ptr;

/* Special files */
#define STANDARD_PDB          0
#define UPLOADED_PDB          1
#define MCSG_PDB              2

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  int                     pdb_type;
  int                     pdbsum1;
  int                     run_64;
};

/* Structure for each RasMol definition entry
   ------------------------------------------ */
struct rasdef
{
  int                type;
  int                number;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                count;
  char               colour[COLNAM_LEN];
  char               *rasmol_ligand_def;
  char               ligand_from[7];
  char               ligand_to[7];
  char               interacts_with[MAX_INT_CHAINS];
  int                nhets;
  struct het         *first_het_ptr;
  struct het         *last_het_ptr;
  struct rasdef      *range_rasdef_ptr;
  struct rasdef      *next_rasdef_ptr;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  char               copy_of;
  int                nresid;
  struct residue     *first_residue_ptr;
  struct chain       *next_chain_ptr;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  char               aa_code;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  char               sec_struc;
  int                beta_turn;
  int                gamma_turn;
  int                disulphide_id;
  int                nhbonds[CONTACT_TYPES];
  int                ncontacts[CONTACT_TYPES];
  int                domain;
  int                nsites;
  int                npatterns;
  int                pattern_pos[MAX_PATTERNS];
  struct atom        *first_atom_ptr;
  struct chain       *chain_ptr;
  struct rlink       *first_rlink_ptr;
  struct site        *site_ptr[MAX_SITES];
  struct pattern     *pattern_ptr[MAX_PATTERNS];
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

/* Structure for each active site residue
   -------------------------------------- */
struct site
{
  char               *name;
  char               *desc;
  int                type;
  int                number;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  struct site        *next_site_ptr;
};

/* Structure for each MOTIF data item
   ---------------------------------- */
struct motif
{
  int                data_type;
  char               sec_struc;
  char               res_num1[6];
  char               chain1;
  char               res_num2[6];
  char               chain2;
  struct motif       *next_motif_ptr;
};

/* Structure for each PROSITE pattern in a given PDB file
   ------------------------------------------------------ */
struct pattern
{
  int                     pattern_no;
  char                    prosite_id[PROSITE_PATLEN];
  char                    chain;
  int                     occurrence;
  char                    first_res_num[6];
  char                    last_res_num[6];
  int                     position;
  int                     true_false;
  char                    *residue_string;
  char                    *colours_string;
  struct pattern          *next_pattern_ptr;
};

/* Structure for each GROW data item
   ------------------------------------- */
struct grow
{
  char                    data_type;
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  char                    other_molecule_type;
  struct rasdef           *rasdef_ptr;
  struct grow             *next_grow_ptr;
};

/* Structure for each het data item
   ------------------------------------- */
struct het
{
  char                    name[4];
  char                    *desc;
  struct het              *next_het_ptr;
};

/* Structure for each rlink data item
   ------------------------------------- */
struct rlink
{
  struct rasdef           *rasdef_ptr;
  int                     nhbonds;
  int                     ncontacts;
  struct rlink            *next_rlink_ptr;
};

