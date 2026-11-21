/***********************************************************************

  wiring.h - Include file for wiring.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN          12
#define FILENAME_LEN       512
#define FRAG_LEN            60
#define LINELEN         100000
#define MAX_ALLSITES       500
#define MAX_CHANGE          10
#define MAX_LABEL_WIDTH    200
#define MAX_LEVEL          100
#define MAX_PATTERNS         5
#define MAX_SEQLEN        5000
#define MAX_LIGANDS          5
#define MAX_METALS           5
#define MAX_SITES            5
#define MAX_VARIANTS        20
#define MAX_TOKEN          100
#define NAMINO              20
#define NATOMS               5
#define PFAM_LEN            11
#define TOKEN_LEN          100

/* Ligand record types */
#define LIGAND                0
#define METAL                 1

/* Coordinate names */
#define X                     0
#define Y                     1
#define Z                     2

/* Min/max names */
#define MIN                   0
#define MAX                   1

/* Pfam domain types */
#define UNKNOWN              -1
#define PFAMA                 0
#define PFAMB                 1

#define NDOMAIN_TYPES         2

/* Image dimension subscripts */
#define COLUMN                0
#define ROW                   1

/* Fonts read in from .ppm files */
#define NFONTS                3

#define FONT_TEXT             0
#define FONT_NUM              1
#define GREEK                 2

/* Characters defined in wirletters.ppm file */
#define CAPITAL_I             9
#define CAPITAL_X            23
#define MAX_LETTER_HEIGHT    60

/* Gap between letters */
#define LETTER_WIDTH          1
#define LETTER_WIDTH_LARGE    4
#define LABEL_SHIFT           1
#define LABEL_SHIFT_LARGE     3

/* Residues per line parameters */
#define MIN_RESID_PER_LINE   20
#define NRESID_PER_LINE      60
#define SAS_NRESID_PER_LINE 100

/* Global parameter defining plot parameters */
static int I;

/* Plot size options */
#define STANDARD_WIRING       0
#define LARGE_WIRING          1

#define NSIZE                 2

#define STANDARD_PLOT         0
#define GO_LARGE              1

/* Plot size parameters */
#define OVER_CONTACTS         0
#define UNDER_CONTACTS        1
#define SITE_HEIGHT           7
#define ERROR_MARGIN         50

static const int COIL_DROP[NSIZE] = { 4, 8 };
static const int HELIX_DROP[NSIZE] = { 8, 25 };
static const int CONTACT_HEIGHT[NSIZE] = { 3, 9 };
static const int CONTACT_WIDTH[NSIZE] = { 3, 9 };
static const int CONTACT_START1[NSIZE] = { 2, 15 };
static const int CONTACT_START2[NSIZE] = { 0, 7 };
static const int CONTACT_START3[NSIZE] = { 0, 5 };
static const int CONTACT_GAP2[NSIZE] = { 1, 5 };
static const int CONTACT_GAP3[NSIZE] = { -1, 2 };
static const int CONTACT_YSHIFT3[NSIZE] = { 2, 0 };

static const int LABEL_MARGIN[NSIZE] = { 4, 10 };
static const int LINE_GAP[NSIZE] = { 5, 25 };
static const int LINE_WIDTH[NSIZE] = { 1, 3 };
static const int LR_MARGIN[NSIZE] = { 10, 50 };
static const int PROSITE_HEIGHT[NSIZE] = { 3, 15 };
static const int PROSITE_WIDTH[NSIZE]  = { 5, 25 };
static const int OVER_LETTER_MARGIN[NSIZE] = { 2, 10 };
static const int RESIDUE_WIDTH[NSIZE] = { 7, 35 };
static const int SEC_STRUC_HEIGHT[NSIZE] = { 20, 100 };
static const int SITE_WIDTH[NSIZE] = { 9, 45 };
static const int TICK_HEIGHT[NSIZE]  = { 2, 10 };
static const int UNDER_LETTER_MARGIN[NSIZE] = { 2, 10 };

/* Helix splice parameters */
static const int FIRST_SPLICE[NSIZE][8]
  = {  9, 20, 26, 42, 48, 63,  84, 69,
      50,104,161,214,270,323, 423,378 };
static const int END_SPLICE1[NSIZE][8]
  = { 10,  0, 15,  0, 15,  0,   0, 15,
      45,  0, 45,  0, 45,  0,   0, 45 };
static const int END_SPLICE2[NSIZE][8]
  = {  0, 10,  0, 10,  0, 10,   0, 10,
       0, 51,  0, 51,  0, 51,   0, 51 };
static const int TURN_START[NSIZE] = { 48, 270 };
static const int FULL_TURN[NSIZE] = { 22, 109 };
static const int PART_TURN[NSIZE] = { 15, 51 };
static const int END2_START[NSIZE] = { 84, 423 };

/* Strand splice parameters */
static const int ARROW_START[NSIZE] = { 14, 72 };
static const int ARROW_MIDDLE[NSIZE] = { 21, 105 };
static const int STRAND_MIDDLE[NSIZE] = { 7, 34 };

/* Disulphide splice parameters */
static const int LEVEL_GAP[NSIZE] = { 4, 30 };
static const int LEVEL_HEIGHT[NSIZE] = { 5, 13 };
static const int DISARROW_HEIGHT[NSIZE] = { 4, 12 };
static const int DISULPHIDE_COLWIDTH[NSIZE] = { 3, 6 };
static const int DISULPHIDE_HEIGHT[NSIZE] = { 5, 13 };
static const int DISULPHIDE_SPAN[NSIZE] = { 7, 35 };
static const int DISULPHIDE_WIDTH[NSIZE] = { 3, 9 };
static const int EXTERNAL_START[NSIZE] = { 7, 35 };
static const int EXTERNAL_WIDTH[NSIZE] = { 7, 16 };
static const int DIMARKER_LETTER_WIDTH[NSIZE] = { 3, 20 };

/* Catalytic residue border size */
static const int CATRES_BORDER[NSIZE] = { 1, 3 };

/* Old version
static const int FIRST_SPLICE[8]
  = { 12, 26, 34, 54, 60, 82, 110, 90 };
static const int END_SPLICE1[8]
  = { 14,  0, 20,  0, 20,  0,   0, 20 };
static const int END_SPLICE2[8]
  = {  0, 12,  0, 12,  0, 12,   0, 12 };
static const int TURN_START = 60;
static const int FULL_TURN = 30;
static const int PART_TURN = 22;
static const int END2_START = 110;
*/

/* Disulphide ends */
#define N_A                   0
#define START                 1
#define END                   2

/* Disulphide levels */
#define HELIX_LEVEL           1
#define STRAND_LEVEL          1

/* Disuphide marker heights */
#define OVER_HELIX            0
#define OVER_STRAND           0

/* Contact types */
#define CONTACT_TYPES         3
#define NOT_WANTED           -9
#define TO_DNA                0
#define TO_LIGAND             1
#define TO_METAL              2

/* Contact colours */
static const int CONTACT_COLOUR[CONTACT_TYPES] = { GREEN, RED, BLUE };

/* Annotation types */
#define CATRES_RECORD         0
#define PROSITE_PATTERN       1
#define PDB_SITE              2
#define LIGAND_CONTACT        3
#define METAL_CONTACT         4
#define DNA_CONTACT           5

#define NANNOT                6

static char *ANNOT_TYPE[NANNOT] = { 
  "CATRES", "PROSITE", "SITE", "LIGAND", "METAL", "DNA"
};

/* Data items to be plotted on wiring diagram */
#define UNDEFINED            -1
#define DOMAINS_CATH          0
#define SITE_RECORDS          1
#define DNA_CONTACTS          2
#define LIGAND_CONTACTS       3
#define METAL_CONTACTS        4
#define DISULPH               5
#define BTURNS                6
#define GTURNS                7
#define PROSITE_PATTERNS      8
#define CSA_RESIDUES          9
#define BHAIRPINS            10
#define NON_STANDARD_AA      11
#define EXTERN_DISULPH       12
#define CALPHA_ONLY          13
#define ALPHA_HELIX          14
#define BETA_STRAND          15
#define RESIDUE_CODE         16
#define HUMAN_VARIANTS       17
#define PFAM_DOMAINS         18

#define NITEMS               19

static const char *ITEM_NAME[NITEMS] = {
  "CATH_DOMAINS", "SITE_RECORDS", "DNA_CONTACTS", "LIGAND_CONTACTS",
  "METAL_CONTACTS", "DISULPH", "BTURNS", "GTURNS", "PROSITE_PATTERNS",
  "CSA_RESIDUES", "BHAIRPINS", "NON_STANDARD_AA", "EXTERN_DISULPH",
  "CALPHA_ONLY", "ALPHA_HELIX", "BETA_STRAND", "RESIDUE_CODE",
  "HUMAN_VARIANTS", "PFAM_DOMAINS"
};

#define NONE                  0
#define MOTIF                 1
#define CONTACT               2

static const int ITEM_TYPE[NITEMS] = {
  NONE, NONE, CONTACT, CONTACT,
  CONTACT, MOTIF, MOTIF, MOTIF, NONE,
  NONE, MOTIF, NONE, NONE,
  NONE, NONE, NONE, NONE, NONE
};

/* PROMOTIF data */
#define NMOTIFS               6

static const char *MOTIF_NAME[NMOTIFS] = {
  "HELICES", "STRANDS", "BETA_TURNS", "GAMMA_TURNS", "DISULPHIDES",
  "BETA_HAIRPINS"
};

/* Annotation description types */
#define LIGAND_DESC           0
#define METAL_DESC            1
#define SITE_DESC             2
#define VARIANT_DESC          3

#define NDESC_TYPES           4

static const char *DESC_NAME[NDESC_TYPES] = {
  "LIGAND", "METAL", "SITE", "VARIANT"
};

/* PROMOTIF types */
#define HELICES               0
#define STRANDS               1
#define BETA_TURNS            2
#define GAMMA_TURNS           3
#define DISULPHIDES           4
#define BETA_HAIRPINS         5

/* Helix parameters for generating template images */
#define MIN_HELIX            14
#define MAX_HELIX            14
#define DUMMY_X      ((float) -5.0)
#define DUMMY_Y      ((float)  0.0)
#define DUMMY_Z      ((float)  0.0)

/* Secondary structure origin */
#define ACTUAL            0
#define FROM_HOMOLOGUE    1
#define PREDICTED         2

/* Number of secondary structure colour types */
#define NSTRUC_TYPES      3

/* Corresponding chain numbers for appropriate colours */
static const int STRUC_CHAIN[NSTRUC_TYPES] = { 0, 1, 5 };

/* Main-chain atom names and corresponding bond lengths */
static const char *ATOM_NAME[] = { " N  ", " CA ", " C  " };
static const float BOND_LENGTH[3] = { 1.47, 1.53, 1.32 };

/* Special files */
#define NO_PDBSUM            -1
#define STANDARD_PDB          0
#define UPLOADED_PDB          1
#define MCSG_PDB              2

/* Secondary structure types */
#define GAP_REGION           -9
#define COIL                  0
#define HELIX                 1
#define STRAND                2
#define SEQUENCE_BREAK        9

/* Parameters for determining which is the dominant secondary structure */
#define DOMIN_HELIX       ( 45.0 )
#define DOMIN_STRAND          6

/* Secondary structure type colours */
static const int SST_TYPE[3] = { HELIX, STRAND, COIL };
static const char *SST_NAME[3] = { "COIL", "HELIX", "STRAND" };

/* Wiring diagram parameters */

/* PfamB template repeat width */
#define PFAM_TEMP_WIDTH       12
#define PFAM_TEMP_WIDTH_LARGE 62

static const int PFAM_TEMP_COL[3] = { 139, 87, 187 };

#define PFAM_MARGIN            6
#define PFAM_MARGIN_LARGE     30

#define CATH_MARGIN            5
#define CATH_MARGIN_LARGE     20

static const int CATH_TEMP_COL[3] = { 139, 87, 187 };

/* Cutoff for applying greying factor for PfamB domain colouring */
#define GREYING_CUTOFF     10.0

/* Image smoothing parameters */
#define NSMOOTH               1

/* Domain types */
#define NO_DOMAIN            -1
#define CATH                  0
#define PFAM                  1

/* States during reading through the CATH domain line */
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

/* Conservation colours */
#define NCONS_COLOURS        10

static char *CONSERVATION_COLOUR[NCONS_COLOURS]
= {
  "black", "blue", "purple", "sky_blue", "cyan",
  "green", "yellow", "orange", "pink", "red"
};

/* Definitions of SITE markers for residues belonging to 1, 2 or 3 sites */
static const char *SITE_DEFN[3][SITE_HEIGHT]
= {
  "000000000",
  " 0111110 ",
  " 0111110 ",
  "  01110  ",
  "   010   ",
  "   010   ",
  "    0    ",
  "000000000",
  " 0112220 ",
  " 0111220 ",
  "  01220  ",
  "   010   ",
  "   020   ",
  "    0    ",
  "000000000",
  " 0111220 ",
  " 0112220 ",
  "  03330  ",
  "   030   ",
  "   030   ",
  "    0    ",
};

#define NSITE_COLOURS         8

static char *SITE_COLOUR[NSITE_COLOURS]
= {
  "green", "blue", "red", "yellow", "purple", "brown", "pink", "cyan"
};

#define NPROSITE_COLOURS      9

static char *PROSITE_COLOUR[NPROSITE_COLOURS]
= {
  "red", "redorange", "orange", "yellow", "green", "greenblue",
  "purple", "sky_blue", "blue"
};

#define NDOMAIN_COLOURS       8

static char *DOMAIN_COLOUR[NDOMAIN_COLOURS]
= {
  "red", "blue", "green", "purple", "brown", "orange", "cyan", "pink"
};

/* Image parameters */
#define COLOUR_MAP_DIFF    ( 50.0 )
#define IMAGE_WIDTH         250
#define IMAGE_ALLOW          20
#define MIN_BLANKS            5
#define MIN_BLANKS_LARGE     25

/* Max-min labels */
#define XMIN                 0
#define XMAX                 1
#define YMIN                 2
#define YMAX                 3

static char *MAP_TYPE[4] = { "XMIN", "XMAX", "YMIN", "YMAX" };

/* Domain label parameters */
#define SPACE_WIDTH          3

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  char                    chain;
  int                     generate;
  int                     residue_conservation;
  int                     standard;
  int                     from_sas;
  int                     nresid_per_line;
  int                     go_large;
  int                     chain_colour;
  int                     user_id;
  int                     pdb_type;
  int                     pfam_domains;
  int                     motifs;
  int                     run_64;
};

/* Structure for each atom record written out
   ------------------------------------------ */
struct atom
{
  int                atom_number;
  char               atom_name[5];
  double             coords[3];
  float              bvalue;
  struct residue     *residue_ptr;
  struct atom        *next_atom_ptr;
};

/* Structure for each annotation description
   ----------------------------------------- */
struct desc
{
  int                type;
  int                number;
  char               *description;
  struct desc        *next_desc_ptr;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  char               aa_code;
  int                seq_no;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  char               sec_struc;
  int                sec_struc_origin;
  char               *sst_id;
  int                key;
  int                cat_res;
  int                beta_turn;
  int                gamma_turn;
  int                disulphide_id;
  int                disulphide_end;
  int                disulphide_level;
  int                disulphide_span;
  char               *link;
  int                hpstart;
  int                hpend;
  int                hp_drop;
  int                nhbonds[CONTACT_TYPES];
  int                ncontacts[CONTACT_TYPES];
  int                have_contacts;
  int                domain;
  int                nsites;
  int                site[MAX_SITES];
  int                have_sites;
  int                nligands;
  int                ligand_no[MAX_LIGANDS];
  int                nmetals;
  int                metal_no[MAX_METALS];
  int                npatterns;
  int                pattern_id[MAX_PATTERNS];
  int                pattern_colour[MAX_PATTERNS];
  int                nvariants;
  int                variant[MAX_VARIANTS];
  int                pdrop;
  int                similarity;
  float              pred_confidence;
  int                conserv;
  int                natoms;
  int                from_annot[NITEMS];
  char               aa_uniprot;
  char               range;
  struct annot       *first_annot_ptr;
  struct atom        *first_atom_ptr;
  struct domlist     *first_domlist_ptr;
  struct domlist     *last_domlist_ptr;
  struct residue     *next_residue_ptr;
};

/* Structure for domain item
   ------------------------- */
struct domain
{
  int                     type;
  int                     domain_no;
  int                     start_res;
  char                    start_res_name[4];
  char                    start_res_num[6];
  char                    start_chain;
  int                     end_res;
  char                    end_res_name[4];
  char                    end_res_num[6];
  char                    end_chain;
  int                     colour;
  int                     map[4];
  struct cath             *cath_ptr;
  struct pfam             *pfam_ptr;
  struct domain           *next_domain_ptr;
};

/* Structure for residue's domain list
   ----------------------------------- */
struct domlist
{
  struct domain           *domain_ptr;
  struct domlist          *next_domlist_ptr;
};

/* Structure for CATH domain
   ------------------------- */
struct cath
{
  char                    *cath_code;
  char                    *class;
  char                    *architecture;
  struct cath             *next_cath_ptr;
};

/* Structure for Pfam domain
   ------------------------- */
struct pfam
{
  char                    *pfam_id;
  int                     pfam_no;
  int                     colour;
  char                    *name;
  char                    *description;
  int                     type;
  int                     ncopy;
  struct domain           *first_domain_ptr;
  struct pfam             *next_pfam_ptr;
};

/* Structure for each PROMOTIF data item
   ------------------------------------- */
struct promotif
{
  int                motif;
  int                identifier[2];
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  struct promotif    *next_promotif_ptr;
};

/* Structure for each colour mapping
   --------------------------------- */
struct colourmap
{
  int                colour;
  int                other_colour[3];
  int                missing;
  struct colourmap   *use_colourmap_ptr;
  struct colourmap   *next_colourmap_ptr;
};

/* Structure for image lines
   ------------------------- */
struct imgline
{
  int                row;
  int                start_col;
  int                end_col;
  int                gap;
  int                blank;
  int                deleted;
};

/* Structure for each annotation entry
   ----------------------------------- */
struct annot
{
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                cat_res;
  int                nhbonds[CONTACT_TYPES];
  int                ncontacts[CONTACT_TYPES];
  int                domain;
  int                nsites;
  int                site[MAX_SITES];
  int                npatterns;
  int                pattern_id[MAX_PATTERNS];
  int                pattern_colour[MAX_PATTERNS];
  struct annot       *next_annot_ptr;
};

/* Structure for each image map entry
   ---------------------------------- */
struct map
{
  int                type;
  char               *description;
  int                limits[4];
  int                deleted;
  struct residue     *residue_ptr;
  struct map         *next_map_ptr;
};

/* Structure for each label entry
   ------------------------------ */
struct label
{
  int                xpos;
  int                ypos;
  int                len;
  int                letter_gap;
  int                space_width;
  int                backfill;
  int                colour;
  struct letter      **letter_array_ptr;
  struct label       *next_label_ptr;
};

