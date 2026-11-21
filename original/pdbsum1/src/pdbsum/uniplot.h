/***********************************************************************

  uniplot.h - Include file for uniplot.c

***********************************************************************/

/* Testing without generating wiring diagrams */
/* debug */
#define SKIP              FALSE
/* debug */

/* Array parameters */
#define CODELEN             20
#define COLNAM_LEN          12
#define FILENAME_LEN       512
#define FRAG_LEN            60
#define LETTER_ROWS          9
#define LINELEN         100000
#define LONG_LINELEN    500000
#define MAX_CODES          100
#define MAX_DOMAINS       1000
#define MAX_LABEL_WIDTH    200
#define MAX_RESNAMES      1000
#define MAX_PATTERNS         5
#define MAX_SEQLEN      100000
#define MAX_SITES            5
#define MAX_TOKEN       200000
#define MAX_TOKENS         100
#define NAME_LEN           512
#define NAMINO              20
#define NLETTERS            80
#define PFAM_LEN            20
#define RESNAMES_LEN    (11 * MAX_RESNAMES)
#define SPACE_WIDTH          3
#define UNIPROT_LEN         20
#define TOKEN_LEN          100

/* Record types defined in the rasmol.dfn file */
#define OTHER                -1
#define CA_ONLY               0
#define LIGAND                1
#define METAL                 2
#define NUCLEIC_ACID          3
#define PROTEIN               4
#define LIGAND_RANGE          5

/* Max-min labels */
#define XMIN                 0
#define XMAX                 1
#define YMIN                 2
#define YMAX                 3

static char *MAP_TYPE[4] = { "XMIN", "XMAX", "YMIN", "YMAX" };

/* Update types */
#define STANDARD              0
#define WEEKLY_UPDATE         1

/* PDB file types */
#define STANDARD_PDB          0
#define UPLOAD_PDB            1
#define MCSG_PDB              2

/* Contact types */
#define CONTACT_TYPES         3
#define NOT_WANTED           -9
#define TO_DNA                0
#define TO_LIGAND             1
#define TO_METAL              2

/* Secondary structure types */
#define GAP_REGION           -9
#define COIL                  0
#define HELIX                 1
#define STRAND                2
#define SEQUENCE_BREAK        9

/* Wiring diagram parameters */
#define CATH_HEIGHT           5
#define SEQ_HEIGHT           11
#define WIR_HEIGHT           11
#define WIR_GAP               4
#define PIC_HEIGHT  (CATH_HEIGHT + SEQ_HEIGHT + WIR_HEIGHT + WIR_GAP)
#define BREAK_HEIGHT          7
#define COIL_HEIGHT           3
#define HELIX_HEIGHT          9
#define STRAND_HEIGHT        11
#define TICK_HEIGHT           5
#define OLD_HELIX_HEIGHT     11
#define OLD_STRAND_HEIGHT    11
#define PFAMA_HEIGHT         11
#define PFAMB_HEIGHT  ( PFAMA_HEIGHT )

/* Domain types */
#define CATH                  0
#define PFAM                  1

/* Upper CATH domain colour */
#define UPPER_DOMAIN_COL 7

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

/* Pfam domain types */
#define UNKNOWN              -1
#define PFAMA                 0
#define PFAMB                 1
#define COVID19               2

#define NDOMAIN_TYPES         3

/* PDB residue number components */
#define HYPHENS               0
#define NUMBERS               1

/* Sort orders */
#define LIG_PFAM              0
#define PFAM_LIG              1

/* Pfam-ligand domain types */
#define SINGLE                0
#define PART_OF_MULTIPLE      1
#define MULTIPLE              2

static char DOMAIN_TYPE[] = "SPM";

/* Entry types */
#define NO_COVID              0
#define COVID_ONLY            1
#define MIXED                 2

/* Commands set according to operating system */
char CP_COMMAND[6];
char MV_COMMAND[6];
char RM_COMMAND[6];
char SLEEP_COMMAND[30];

/* Minimum ligand-domain contacts for multi-domain interactions */
#define MIN_CONTACTS          0

static int DOMAIN_SHADE[PFAMA_HEIGHT]
= { 50, 75, 100, -1, 100, 100, 90, 75, 60, 40, 20 };

static float ADJUST_COL[10]
= { 0, 0.2, 0.4, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8 };

static char *COIL_DEFN[COIL_HEIGHT]
= {
  "8",
  "5",
  "2"
};
#define LEN_COIL_DEFN        1

static char *HELIX_DEFN[HELIX_HEIGHT]
= {
  "9 3",
  "85 ",
  "58 ",
  " 57",
  " 55",
  " 55",
  " 22",
  "  5",
  "3 8"
};
#define LEN_HELIX_DEFN        3

static char *OLD_HELIX_DEFN[OLD_HELIX_HEIGHT]
= {
  " @    ",
  " *    ",
  " *    ",
  "@*@ @ ",
  "*@* * ",
  "* * **",
  "  *@**",
  "  @*@@",
  "   *  ",
  "   *  ",
  "   @  ",
};
#define LEN_OLD_HELIX_DEFN    6

static char *STRAND_DEFN[STRAND_HEIGHT]
= {
  " 9  ",
  " 91 ",
  "999 ",
  "8881",
  "7777",
  "5555",
  "5555",
  "2252",
  "777 ",
  " 72 ",
  " 8  ",
};
#define LEN_STRAND_DEFN       4

static char *OLD_STRAND_DEFN[OLD_STRAND_HEIGHT]
= {
  "  *  ",
  "  *@ ",
  "**** ",
  "*  *@",
  "*   *",
  "*   *",
  "*   *",
  "*  *@",
  "**** ",
  "  *@ ",
  "  *  "
};
#define LEN_OLD_STRAND_DEFN       5

static char *BREAK_DEFN[BREAK_HEIGHT]
= {
  " ",
  " ",
  "#",
  "#",
  "#",
  " ",
  " "
};
#define LEN_BREAK_DEFN        1

static char *PFAMA_DEFN[PFAMA_HEIGHT]
= {
  "*",
  "*",
  "*",
  "*",
  "*",
  "*",
  "*",
  "*",
  "*",
  "*",
  "*",
};
#define LEN_PFAMA_DEFN        1

static char *PFAMB_DEFN[PFAMB_HEIGHT]
= {
  "*:::**",
  ":::***",
  "::***:",
  ":***::",
  "***:::",
  "**:::*",
  "*:::**",
  ":::***",
  "::***:",
  ":***::",
  "***:::"
};
#define LEN_PFAMB_DEFN        6

static char *CATH_DEFN[CATH_HEIGHT]
= {
  "  * *  ",
  " ** ** ",
  "*******",
  " ** ** ",
  "  * *  "
};
#define LEN_CATH_DEFN         7

/* Sequence subscripts */
#define NO_SEQUENCE          -1
#define PDB_SEQ               0
#define UNIPROT_SEQ           1
#define SIMILARITY_MARKERS    2

/* Audit parameters */
#define NAUDIT           100000
#define NAUDIT_PDB         1000

/* Image parameters */
#define IMAGE_WIDTH         250
#define IMAGE_ALLOW          20

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  int                     pdb_type;
  int                     pdbsum_type;
  int                     run_type;
  int                     covid19;
  int                     wget;
  int                     run_64;
};

/* Structure for each PDB entry
   ---------------------------- */
struct pdb
{
  char               pdb_code[5];
  int                nchains;
  int                nunique;
  int                have_cath;
  int                have_pfam;
  int                coloured;
  int                done;
  struct chain       *first_chain_ptr;
  struct pdb         *next_pdb_ptr;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char               chain_id;
  char               copy_of;
  int                ncodes;
  int                nimages[MAX_CODES];
  int                chain_no;
  int                plot_no[MAX_CODES];
  int                nresid;
  int                image_height;
  int                image_width;
  int                last_width;
  int                ncath_dom;
  int                nligs;
  int                nligands;
  char               *alignment[2];
  int                align_len;
  int                nidentical;
  float              identity;
  struct pdb         *pdb_ptr;
  struct code        *code_ptr[MAX_CODES];
  struct domain      *first_cath_domain_ptr;
  struct smatch      *smatch_ptr[MAX_CODES];
  struct ligand      *first_ligand_ptr;
  struct ligand      *last_ligand_ptr;
  struct liglist     *first_liglist_ptr;
  struct liglist     *last_liglist_ptr;
  struct chain       *next_chain_ptr;
};

/* Structure for sequence match item
   --------------------------------- */
struct smatch
{
  int                align_len;
  int                nmismatch;
  char               *uniprot_seq;
  char               *pdb_seq;
  char               *pdb_hyphens;
  char               *pdb_resno;
  char               *diffs;
  struct code        *code_ptr;
  struct domain      *first_domain_ptr;
};

/* Structure for each code data item
   ---------------------------------- */
struct code
{
  char               uniprot_acc[UNIPROT_LEN];
  char               uniprot_id[UNIPROT_LEN];
  char               *name;
  int                order;
  int                seq_len;
  int                pfam_seq_len;
  int                have_seq;
  int                covid_protein;
  int                done;
  struct clist       *first_clist_ptr;
  struct clist       *last_clist_ptr;
  struct domain      *first_domain_ptr;
  struct code        *next_code_ptr;
};

/* Structure for clist item
   ------------------------ */
struct clist
{
  struct chain            *chain_ptr;
  int                     icode;
  struct clist            *next_clist_ptr;
};

/* Structure for domain item
   ------------------------- */
struct domain
{
  char                    *id;
  int                     type;
  int                     domain_no;
  int                     start_res;
  int                     end_res;
  int                     colour;
  int                     map[4];
  struct cath             *cath_ptr;
  struct pfam             *pfam_ptr;
  struct domain           *next_domain_ptr;
};

/* Structure for Pfam domain
   ------------------------- */
struct pfam
{
  char                    *pfam_id;
  char                    *name;
  char                    *description;
  int                     type;
  struct pfam             *next_pfam_ptr;
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
  int                cat_res;
  int                beta_turn;
  int                gamma_turn;
  int                disulphide_id;
  int                nhbonds[CONTACT_TYPES];
  int                ncontacts[CONTACT_TYPES];
  int                domain;
  int                nsites;
  char               *site_name[MAX_SITES];
  int                npatterns;
  int                pattern_id[MAX_PATTERNS];
  int                pattern_colour[MAX_PATTERNS];
  int                similarity;
  float              pred_confidence;
  int                pfam_no;
  int                nligs;
  int                range;
  struct chain       *chain_ptr;
  struct pfam        *pfam_ptr;
  struct liglist     *first_liglist_ptr;
  struct liglist     *last_liglist_ptr;
  struct residue     *next_residue_ptr;
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

/* Structure for pfamlig item
   ------------------------ */
struct pfamlig
{
  char               *ligand_name;
  char               *domain_list;
  char               ligand_id[7];
  char               pdb_code[5];
  char               chain;
  int                seq_no;
  int                max_contacts;
  int                ncontacts;
  int                domain_no;
  int                done;
  int                type;
  struct code        *code_ptr;
  struct pfam        *pfam_ptr;
  struct plist       *first_plist_ptr;
  struct pfamlig     *next_pfamlig_ptr;
};

/* Structure for plist item
   ------------------------ */
struct plist
{
  int                domain_no;
  int                ncontacts;
  struct pfam        *pfam_ptr;
  struct plist       *next_plist_ptr;
};

/* Structure for liglist item
   -------------------------- */
struct liglist
{
  struct r_rasdef    *r_rasdef_ptr;
  struct liglist     *next_liglist_ptr;
};

/* Structure for ligand
   -------------------- */
struct ligand
{
  char                    *ligand_name;
  char                    *description;
  int                     type;
  int                     number;
  int                     nhets;
  struct het              *first_het_ptr;
  struct ligand           *next_ligand_ptr;
};

/* Structure for each het data item
   ------------------------------------- */
struct het
{
  char                    name[4];
  char                    *desc;
  struct het              *next_het_ptr;
};

/* Structure for Het Group -> drug id record
   ------------------------------------------ */
struct hetdrug
{
  char              het_name[4];
  char              *generic_name;
  char              *drug_id;
  char              *drug_desc;
  struct hetdrug    *next_hetdrug_ptr;
};

