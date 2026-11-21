/***********************************************************************

  readpdb.h - Include file for readpdb.c

***********************************************************************/

/* Array parameters */
#define R_COLNAM_LEN         30
#define R_FILENAME_LEN     1024
#define R_LINELEN          5120
#define R_LONG_LINELEN    40960
#define R_MAX_AUTHORS       100
#define R_MAX_BVALUE   ( 999.99)
#define R_MAX_CONECT          5
#define R_MAX_LEVEL           4
#define R_MAX_MOLECS       2000
#define R_MAX_TOKEN        1000
#define R_PUBMED_LEN         24
#define R_STRING_LEN        (10 * R_LINELEN)
#define R_TOKEN_LEN         200

/* Record types defined in the rasmol.dfn file */
#define R_UNKNOWN            -1
#define R_CALPHA_ONLY         0
#define R_LIGAND              1
#define R_METAL               2
#define R_NUCLEIC_ACID        3
#define R_PROTEIN             4
#define R_LIGAND_RANGE        5
#define R_WATER               6
#define R_GAP                 7

#define R_NRESID_TYPES        8

static char *R_RESTYPE_NAME[R_NRESID_TYPES] = {
  "R_CALPHA_ONLY ", "R_LIGAND      ", "R_METAL       ",
  "R_NUCLEIC_ACID", "R_PROTEIN     ", "R_LIGAND_RANGE",
  "R_WATER       ", "R_GAP         "
};

/* PDB key words */
#define NONE                  0
#define TITLE                 1
#define COMPOUND              2
#define SOURCE                3
#define KEY_WORDS             4
#define AUTHOR                5
#define JOURNAL               6
#define REFERENCE             7
#define MODEL                 8
#define ENDMDL                9
#define ATOM                 10
#define HETATM               11
#define TER                  12
#define CONECT               13

#define NKEY_WORDS           14

static char *key[NKEY_WORDS] = {
  "NONE  ", "TITLE ", "COMPND", "SOURCE", "KEYWDS", "AUTHOR", "JRNL  ",
  "REMARK   1 ", "MODEL ", "ENDMDL", "ATOM  ", "HETATM", "TER   ",
  "CONECT"
};

/* PDB data required from PDB file */
#define COORDS_ONLY           0
#define HEADERS_ONLY          1
#define ALL_DATA              2

/* Flags relating to specific keywords */
#define MOLECULE_NAME         1
#define SYNONYMS              2
#define SCIENTIFIC_NAME       3
#define COMMON_NAME           4
#define TAXONOMY_ID           5

#define NCONECT_PER_LINE      5

/* Reference types */
#define KEY_REF               0
#define SECONDARY_REF         1
#define ADDITIONAL_REF        2
#define CITING_REF            3

#define R_NREF_TYPES          4

static char *R_REFERENCE_TYPE[R_NREF_TYPES] = {
  "KEY_REF", "SECONDARY_REF", "ADDITIONAL_REF", "CITING_REF"
};

/* Record types */
#define R_ATOM_REC            0
#define R_HETATM_REC          1
#define R_NOT_WANTED          9

/* Chain types */
#define R_CHAIN_UNKNOWN              -1
#define R_CHAIN_PROTEIN               0
#define R_CHAIN_DNA                   1
#define R_CHAIN_LIGAND                2
#define R_CHAIN_CALPHA_ONLY           3
#define R_CHAIN_METAL                 4
#define R_CHAIN_ATTACHMENT            5
#define R_CHAIN_PP_INTERFACE          6
#define R_CHAIN_PD_INTERFACE          7

#define R_NCHAIN_TYPES          8

static char *R_CHAIN_TYPE[R_NCHAIN_TYPES] = {
  "PROTEIN", "DNA/RNA", "LIGAND", "CA_ONLY", "METAL", "ATTACHMENT",
  "PP_INTERFACE", "PD_INTERFACE"
};

/* Ligand parameters */
#define R_LIGAND_LENGTH      10

/* Residue types */
#define R_TYPE_UNKNOWN       -1
#define R_TYPE_STANDARD       0
#define R_TYPE_NON_STANDARD   1
#define R_TYPE_NUCLEIC_ACID   2
#define R_TYPE_METAL_ION      3

/* Secondary structure types */
#define R_OTHER              -1
#define R_COIL                0
#define R_HELIX               1
#define R_STRAND              2

#define R_NSST                3

static char *R_SST_NAME[R_NSST] = { "Coil", "Helix", "Strand" };

static char *R_SST = "-HE";

enum LineType {
    ATOM_RECORD, CONECT_RECORD, ENDMDL_RECORD, HEADER_RECORD, HETATM_RECORD,
    HETNAM_RECORD, MODEL_RECORD, TER_RECORD, TITLE_RECORD
};

/* Metals */
#define R_NMETALS           100

/* Note "AC  " changed to "ac  " below as carbon atoms in certain
   molecules labelled AC1, AC1, etc (eg in FAD). Similarly NP, NO
   and PO */
static char *R_METAL_NAME[R_NMETALS] = {
  "HE  ", "LI  ", "BE  ", " B  ", " F  ",
  "NE  ", "NA  ", "MG  ", "AL  ", "SI  ",
  "CL  ", "AR  ", " K  ", "CA  ", "SC  ",
  "TI  ", " V  ", "CR  ", "MN  ", "FE  ",
  "CO  ", "NI  ", "CU  ", "ZN  ", "GA  ",
  "GE  ", "AS  ", "SE  ", "BR  ", "KR  ",
  "RB  ", "SR  ", " Y  ", "ZR  ", "NB  ",
  "MO  ", "TC  ", "RU  ", "RH  ", "PD  ",
  "AG  ", "CD  ", "IN  ", "SN  ", "SB  ",
  "TE  ", " I  ", "XE  ", "CS  ", "BA  ",
  "LA  ", "CE  ", "PR  ", "ND  ", "PM  ",
  "SM  ", "EU  ", "GD  ", "TB  ", "DY  ",
  "HO  ", "ER  ", "TM  ", "YB  ", "LU  ",
  "HF  ", "TA  ", " W  ", "RE  ", "OS  ",
  "IR  ", "PT  ", "AU  ", "HG  ", "TL  ",
  "PB  ", "BI  ", "po  ", "AT  ", "RN  ",
  "FR  ", "RA  ", "ac  ", "TH  ", "PA  ",
  " U  ", "np  ", "PU  ", "AM  ", "CM  ",
  "BK  ", "CF  ", "ES  ", "FM  ", "MD  ",
  "no  ", "LR  ", "....", "....", "...."
};

/* Predefine the commonly occuring metals */
#define R_NCOMMON_METALS        5
static char *R_COMMON_METAL[R_NCOMMON_METALS] = {
  "CL", "ZN", "CU", "CA", "FE"
};

/* Metal colours */
#define R_NMETAL_COLOURS      5
static char *R_METAL_COLOUR[R_NMETAL_COLOURS] = {
  "green", "cyan", "orange", "pink", "brown"
};

/* Standard atom types and radii */
#define R_NATOM_TYPES         6
static char *R_ELEMENT[R_NATOM_TYPES] = {
  "XX", " N", " C", " O", " S", " P"
};

/* Distance parameters */
#define R_COVALENT_DIST       2.0
#define R_COVALENT_DIST_SQRD (R_COVALENT_DIST * R_COVALENT_DIST)
#define R_MAX_CONECT_LEN      2.5
#define R_MAX_CONECT_LEN_SQRD (R_MAX_CONECT_LEN * R_MAX_CONECT_LEN)
#define R_MAX_METAL_LEN       3.0
#define R_MAX_METAL_LEN_SQRD  (R_MAX_METAL_LEN * R_MAX_METAL_LEN)
#define R_TOO_FAR_DIST        5.0

/* Max atom radius */
#define R_MAX_ATOM_RADIUS     2.0

#define NAMINO              20

/* Residue names */
static char *AA_NAME[] = {
  "XXX", "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
  "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
  "SER", "THR", "VAL", "TRP", "TYR" };

static char *RES_COLOUR_NAME[] = { 
  "white", "sky_blue", "red", "green", "grey", "purple", "brown", "yellow" };
  
static char *RES_COLOUR_DESC[] = { 
  "unknown", "positive", "negative", "neutral", "aliphatic", "aromatic",
  "Pro&Gly", "cysteine" };
  
static int RES_TYPE[NAMINO + 1]
= { 0, 4, 7, 2, 2, 5, 6, 1, 4, 1, 4, 4, 3, 6, 3, 1, 3, 3, 4, 5, 5 };
//  X  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
  
/* Fauchere and Pliska hydrophobicity scale */
static double HYDROPHOBICITY_SCALE[NAMINO] =
{
  0.31,  /* A */
  1.54,  /* C */
  -0.77, /* D */
  -0.64, /* E */
  1.79,  /* F */
  0.00,  /* G */
  0.13,  /* H */
  1.80,  /* I */
  -0.99, /* K */
  1.70,  /* L */
  1.23,  /* M */
  -0.60, /* N */
  0.72,  /* P */
  -0.22, /* Q */
  -1.01, /* R */
  -0.04, /* S */
  0.26,  /* T */
  1.22,  /* V */
  2.25,  /* W */
  0.96,  /* Y */
};

/* Amino acids ordered by increasing hydrophobicity scale */
static char AMINO_INCREASING[]={"RKDENQSGHTAPYVMCLFIW"};

/* Amino acids ordered by decreasing hydrophobicity scale */
static char AMINO_DECREASING[]={"WIFLCMVYPATHGSQNEDKR"};

/* Numbers of atoms in each residue */
static int AA_SIZE[NAMINO] =
{
  5,  /* A */
  6,  /* C */
  8,  /* D */
  9,  /* E */
 11,  /* F */
  4,  /* G */
 10,  /* H */
  8,  /* I */
  9,  /* K */
  8,  /* L */
  8,  /* M */
  8,  /* N */
  7,  /* P */
  9,  /* Q */
 11,  /* R */
  6,  /* S */
  7,  /* T */
  7,  /* V */
 14,  /* W */
 12,  /* Y */
};

/* Conservation colours */
#define NCONS_COLOURS        10

static char *CONSERVATION_COLOUR[NCONS_COLOURS]
= {
  "black", "blue", "purple", "sky_blue", "cyan",
  "green", "yellow", "orange", "pink", "red"
};


/* Structure for each PDB entry
   ---------------------------- */
struct r_pdb
{
  char                *pdb_code;
  char                *protein_name;
  int                 ninfo;
  int                 nmodels;
  int                 natoms;
  int                 nresid;
  int                 nchains;
  int                 nmolecs;
  int                 nligands;
  int                 nmetals;
  struct r_rasdef     *first_r_rasdef_ptr;
  struct r_info       *first_r_info_ptr;
  struct r_ref        *first_r_ref_ptr;
  struct r_chain      *first_r_chain_ptr;
  struct r_molec      *first_r_molec_ptr;
  struct r_ligand     *first_r_ligand_ptr;
  struct r_residue    *first_r_residue_ptr;
  struct r_atom       *first_r_atom_ptr;
  struct r_info       *last_r_info_ptr;
  struct r_ref        *last_r_ref_ptr;
  struct r_hetdesc    *first_r_hetdesc_ptr;
  struct r_chain      *last_r_chain_ptr;
  struct r_residue    *last_r_residue_ptr;
  struct r_atom       *last_r_atom_ptr;
};

/* Structure for each header info item
   ----------------------------------- */
struct r_info
{
  int                     nlevel;
  int                     instance[R_MAX_LEVEL];
  char                    *name;
  char                    *value;
  struct r_info           *next_r_info_ptr;
};

/* Structure for jnl reference item
   -------------------------------- */
struct r_ref
{
  int                     type;
  char                    number[5];
  int                     nauthors;
  char                    *author[R_MAX_AUTHORS];
  int                     to_be_published;
  int                     new;
  char                    *title;
  char                    *journal;
  int                     year;
  int                     volume;
  char                    *start_page;
  char                    *end_page;
  char                    pubmed_id[R_PUBMED_LEN];
  char                    *doi_number;
  struct r_ref            *next_r_ref_ptr;
};

/* Structure for each chain
   ------------------------ */
struct r_chain
{
  char               chain_id;
  int                type;
  int                nresid;
  int                nhelices;
  int                nstrands;
  int                nstand;
  int                nbase;
  int                natom_recs;
  int                nhetatm_recs;
  int                dominant_sst;
  int                symm_chain;
  struct r_chain     *link_r_chain_ptr;
  struct r_pdb       *r_pdb_ptr;
  struct r_residue   *first_r_residue_ptr;
  struct r_chain     *next_r_chain_ptr;
};

/* Structure for coordinate limits
   ------------------------------- */
struct r_limits
{
  int                     first;
  double                  valmin[3];
  double                  valmax[3];
  double                  c_of_g[3];
  int                     npoints;
};

/* Structure for each molecule
   --------------------------- */
struct r_molec
{
  int                type;
  int                nresid;
  struct r_limits    r_limits_ptr;
  struct r_pdb       *r_pdb_ptr;
  struct r_residue   *first_r_residue_ptr;
  struct r_molec     *next_r_molec_ptr;
};

/* Structure for each residue
   -------------------------- */
struct r_residue
{
  char               athet;
  int                type;
  char               aa_code;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                seq_no;
  int                ref_resid;
  double             coords_min[3];
  double             coords_max[3];
  double             centroid[3];
  int                natoms;
  int                domain;
  int                sec_struc;
  int                disulphide_id;
  char               *site;
  int                nmain_chain;
  int                have_calpha;
  int                connected;
  float              accessibility;
  int                wanted;
  struct r_limits    r_limits_ptr;
  struct r_chain     *r_chain_ptr;
  struct r_pdb       *r_pdb_ptr;
  struct r_rasdef    *r_rasdef_ptr;
  struct r_ligand    *r_ligand_ptr;
  struct r_hetdesc   *r_hetdesc_ptr;
  struct r_atom      *first_r_atom_ptr;
  struct r_rbond     *first_r_rbond_ptr;
  struct r_residue   *next_r_residue_ptr;
};

/* Structure for each atom
   ----------------------- */
struct r_atom
{
  int                atom_number;
  char               atom_name[5];
  char               element[3];
  double             coords[3];
  double             orig_coords[3];
  float              bvalue;
  float              occupancy;
  int                type;
  int                rec_type;
  struct r_conect    *first_r_conect_ptr;
  struct r_conect    *last_r_conect_ptr;
  struct r_residue   *r_residue_ptr;
  struct r_atom      *next_r_atom_ptr;
};

/* Structure for each atom CONECT record
   ------------------------------------- */
struct r_conect
{
  struct r_atom           *to_r_atom_ptr;
  struct r_conect         *next_r_conect_ptr;
};

/* Structure for each covalent bond between a pair of residues
   ----------------------------------------------------------- */
struct r_rbond
{
  struct r_residue        *to_r_residue_ptr;
  struct r_rbond          *next_r_rbond_ptr;
};

/* Structure for each RasMol definition entry
   ------------------------------------------ */
struct r_rasdef
{
  int                type;
  int                number;
  int                copy;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  char               colour[R_COLNAM_LEN];
  char               *ligplot_id;
  char               *ligand_name;
  char               *ligand_range;
  char               ligand_from[7];
  char               ligand_to[7];
  int                nresid;
  int                ncopy;
  int                wanted;
  struct r_residue   *first_r_residue_ptr;
  struct r_rasdef    *next_r_rasdef_ptr;
};

/* Structure for each ligand data item
   ---------------------------------- */
struct r_ligand
{
  int                     number;
  char                    *residue_sequence;
  char                    *residue_list;
  char                    *residue_range;
  char                    first_residue[10];
  char                    last_residue[10];
  int                     metal;
  int                     metal_colour;
  char                    metal_name[3];
  int                     type;
  int                     count;
  int                     copy;
  int                     unique;
  int                     unique_count;
  int                     ncopy;
  int                     duplicate_count;
  int                     nduplicates;
  int                     nresid;
  int                     natoms;
  int                     sort_score;
  double                  coord_min[3];
  double                  coord_max[3];
  int                     got_bonds;
  struct r_rasdef         *r_rasdef_ptr;
  struct r_residue        *first_r_residue_ptr;
  struct r_list           *first_r_list_ptr;
  struct r_ligand         *duplicate_r_ligand_ptr;
  struct r_ligand         *next_r_ligand_ptr;
};

/* Structure for each list data item
   --------------------------------- */
struct r_list
{
  struct r_residue        *r_residue_ptr;
  struct r_list           *next_r_list_ptr;
};

/* Structure for each hetdesc data item
   ------------------------------------ */
struct r_hetdesc
{
  char                    res_name[4];
  char                    *details[3];
  int                     count;
  struct r_hetdesc        *next_r_hetdesc_ptr;
};

