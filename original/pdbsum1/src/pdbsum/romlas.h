/***********************************************************************

  romlas.h - Include file for romlas.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN          30
#define COLNO_LEN           14
#define FILENAME_LEN       512
#define LABEL_LEN           15
#define LINELEN         100000
#define MAX_ATOMNO       99999
#define MAX_LOCAL_PDB       27
#define MAX_LOOPS           10
#define MAXPDB            1000
#define MAX_OBJ_COL         14
#define MAX_CHAINS         500
#define MAXCONECT            5
#define MAXDOMCHAIN        200
#define MAX_FIT_ATOMS     5000
#define MAX_FIT_LIST      1000
#define MAX_TOKEN         5000
#define NAME_LEN            80
#define OBJ_NAME_LEN        41
#define PROSITE_PATLEN      15
#define STRING_LEN        1000
#define TOKEN_LEN           20

/* Display options */
#define NONE              0
#define ALL_ATOMS         1
#define MAIN_CHAIN        2
#define C_ALPHA           3
#define CARTOON           4
#define LIGAND_ONLY       5
#define INT_RESIDUES      6
#define LIGAND_AND_RES    7
#define BINDING_SITE      8
#define DOT_SURFACE       9
#define LIGAND_AND_BOX   10

/* Entry types */
#define PDB               0
#define LIGPLOT           1
#define MASK              2
#define SITE_RES          3
#define M_C               4
#define C_A               5
#define CART              6
#define ALL_AT            7
#define DOTS              8
#define XSITE             9
#define LIG_BX           10

/* States during reading through the domain line */
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

/* PROMOTIF motif types */
#define BETA_TURNS        0
#define GAMMA_TURNS       1
#define BETA_BULGES       2
#define BETA_HAIRPINS     3
#define HELICES           4
#define BAB_MOTIFS        5
#define PSI_LOOPS         6
#define SHEETS            7
#define HELIX_INTERACTIONS 8

#define NMOTIFS           9

/* Record types defined in the rasmol.dfn file */
#define UNKNOWN           -1
#define CA_ONLY            0
#define LIGAND             1
#define METAL              2
#define NUCLEIC_ACID       3
#define PROTEIN            4
#define LIGAND_RANGE       5

/* Start and end types */
#define START             0
#define END               1

/* Residues in named list */
#define IN_LIST           0
#define NOT_IN_LIST      -1

/* Global variables
   ---------------- */
int verbose;
FILE *dump_file_ptr;

/* Special files */
#define STANDARD_PDB          0
#define UPLOAD_PDB            1
#define MCSG_PDB              2

/* Factor for rmsd cut-off used during structure fitting */
#define MIN_STORED            3
#define RMSD_FACTOR           2

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  int                     any_het;
  int                     binding_sites;
  int                     rasmol_script;
  int                     write_to_file;
  int                     cath_alignment;
  int                     from_fasta;
  int                     prosite_pattern;
  int                     promotif;
  int                     pdb_site;
  char                    site_id[5];
  int                     motif_type;
  int                     motif_no;
  int                     csa;
  char                    *residue_list;
  char                    rlist_colour[COLNAM_LEN + 1];
  int                     colour_by_pattern;
  int                     display_ligand_only;
  int                     display_ligand_and_residues;
  int                     display_backbone_only;
  int                     display_cartoon;
  int                     display_mainchain_only;
  int                     display_all_atoms;
  int                     display_binding_site_waters;
  int                     all_chains;
  int                     whole;
  int                     metal;
  char                    chains[MAX_CHAINS + 1];
  int                     nchains;
  int                     single_chain;
  int                     single_domain;
  int                     whole_protein;
  int                     include_hets;
  int                     add_sidechains;
  int                     fit_to_backbone;
  int                     fit_to_mainchain;
  int                     label_prot_residues;
  int                     label_ligand;
  int                     label_atoms;
  int                     colour_by_structure;
  int                     sas_colour;
  int                     side_sas;
  int                     nsite_surface_files;
  int                     xsite_spheres;
  int                     profunc;
  int                     reverse_search;
  int                     upload;
  int                     use_wolf;
  int                     split_pdb;
  int                     all_entries;
  int                     fit_out;
  int                     atnum_reset;
  double                  box_cutoff;
  char                    work_dir[FILENAME_LEN];
  struct reslist          *first_reslist_ptr;
  int                     pdb_type;
  int                     jmol;
  int                     pymol;
  char                    *pdb_file;
  int                     run_64;
};

/* Structure for each RasMol definition entry
   ------------------------------------------ */
struct rasdef
{
  int                type;
  int                deleted;
  int                number;
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                count;
  char               colour[COLNAM_LEN];
  char               *rasmol_ligand_def;
  char               ligand_from[7];
  char               ligand_to[7];
  struct atom        *first_atom_ptr;
  struct atom        *last_atom_ptr;
  struct rasdef      *next_rasdef_ptr;
};

/* Structure for residue list
   -------------------------- */
struct reslist
{
  int                     number;
  char                    chain;
  char                    res_num1[6];
  char                    res_num2[6];
  struct residue          *residue1_ptr;
  struct residue          *residue2_ptr;
  struct reslist          *next_reslist_ptr;
};

/* Structure for each sequence in the alignment
   -------------------------------------------- */
struct sequence
{
  char                    pdb_code[5];
  char                    chain;
  int                     domain_no;
  char                    ligand_id[10];
  char                    reslist_name[FILENAME_LEN];
  int                     display_option;
  int                     alignment_column;
  int                     grey;
  char                    colour_name[COLNAM_LEN];
  char                    colour_numbers[COLNO_LEN];
  char                    label[LABEL_LEN];
  int                     nconect;
  struct list             *first_list_ptr;
  struct list             *last_list_ptr;
  struct atom             *first_atom_ptr;
  struct conect           *first_conect_ptr;
  struct hbond            *first_hbond_ptr;
  struct residue          *first_residue_ptr;
  struct residue          *last_residue_ptr;
  struct fitlist          *fitlist_ptr;
  struct fitlist          *first_fitlist_ptr;
  struct fitlist          *last_fitlist_ptr; 
 struct prosite_pattern   *first_prosite_pattern_ptr;
};
struct sequence *first_sequence_ptr, *stored_sequence_ptr[MAXPDB];

/* Structure for coordinate limits
   ------------------------------- */
struct limits
{
  int                     first;
  double                  min_x;
  double                  max_x;
  double                  min_y;
  double                  max_y;
  double                  min_z;
  double                  max_z;
};

/* Structure for each residue required, as read from a list
   -------------------------------------------------------- */
struct list
{
  int                     blank;
  char                    aa_code;
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  int                     colour;
  int                     is_ligand;
  int                     interacting_residue;
  int                     fit;
  struct list             *reference_list_ptr;
  struct residue          *residue_ptr;
  struct list             *next_list_ptr;
};

/* Structure for each residue to be written out
   -------------------------------------------- */
struct residue
{
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  int                     natoms;
  int                     colour;
  int                     type;
  int                     is_ligand;
  int                     in_list;
  int                     interacting_residue;
  int                     from_alignment;
  int                     first_atom_no;
  int                     last_atom_no;
  int                     position;
  int                     check_distance;
  double                  rmsd;
  int                     fit;
  int                     fit_used;
  struct atom             *first_atom_ptr;
  struct atom             *N_atom_ptr;
  struct atom             *CA_atom_ptr;
  struct atom             *C_atom_ptr;
  struct atom             *O_atom_ptr;
  struct rasdef           *rasdef_ptr;
  struct residue          *reference_residue_ptr;
  struct residue          *next_residue_ptr;
};

/* Structure for each atom record to be written out
   ------------------------------------------------ */
struct atom
{
  int                     atom_number;
  int                     out_number;
  int                     display_atom;
  int                     fit_atom;
  char                    at_het;
  char                    atom_name[5];
  double                  x, y, z;
  char                    occup_bvalue[13];
  int                     written;
  struct residue          *residue_ptr;
  struct atom             *next_atom_ptr;
};


/* Structure for each H-bond read in from the grow.out file
   -------------------------------------------------------- */
struct hbond
{
  struct atom             *acceptor_atom_ptr;
  struct atom             *donor_atom_ptr;
  double                  length;
  struct hbond            *next_hbond_ptr;
};


/* Structure for each of the binding-site mask files
   ------------------------------------------------- */
struct surface
{
  char                    site_surface_code[7];
  char                    xsite_probe[3];
  int                     surface_number;
  double                  grid_spacing;
  char                    colour_name[COLNAM_LEN];
  char                    colour_numbers[COLNO_LEN];
  struct vertex           *first_vertex_ptr;
  struct surface          *next_surface_ptr;
};
struct surface *first_surface_ptr;


/* Structure for each of the binding-site mask vertices
   ---------------------------------------------------- */
struct vertex
{
  int                     atom_number;
  double                  x, y, z;
  struct vertex           *next_vertex_ptr;
};


/* Structure for each CONECT record to be written out
   -------------------------------------------------- */
struct conect
{
  struct atom             *atom_ptr[MAXCONECT];
  struct conect           *next_conect_ptr;
};
struct conect *first_conect_ptr, *last_conect_ptr;


/* Structure for each list of fitted residues
   ------------------------------------------ */
struct fitlist
{
  char                    range[10];
  int                     done;
  struct fitlist          *next_fitlist_ptr;
};

/* Structure for each PROSITE pattern to be displayed
   -------------------------------------------------- */
struct match_pattern
{
  char                    prosite_id[PROSITE_PATLEN];
  char                    chain;
  int                     occurrence;
  struct match_pattern    *next_match_pattern_ptr;
};
struct match_pattern *first_match_pattern_ptr;

/* Structure for each PROSITE pattern in a given PDB file
   ------------------------------------------------------ */
struct prosite_pattern
{
  char                    prosite_id[PROSITE_PATLEN];
  char                    chain;
  int                     occurrence;
  char                    first_res_num[6];
  char                    last_res_num[6];
  int                     position;
  char                    residue_string[STRING_LEN];
  char                    colours_string[STRING_LEN];
  struct prosite_pattern  *next_prosite_pattern_ptr;
};
struct prosite_pattern *first_prosite_pattern_ptr;

/* Structure for residue range and how it should be displayed
   ---------------------------------------------------------- */
struct range
{
  char                    range_name[NAME_LEN + 1];
  char                    name[NAME_LEN + 1];
  char                    chain;
  char                    res_name[2][4];
  char                    res_num[2][6];
  char                    colour[COLNAM_LEN + 1];
  char                    label[10];
  int                     sidechain;
  float                   wireframe;
  int                     hbonds;
  struct range            *next_range_ptr;
};


/*************************************************************************

Include variables for Andrew Martin's least-squares fitting routines,
matfit and qikfit

*************************************************************************/

#define SMALL  1.0e-20     /* Convergence cutoffs                       */
#define SMALSN 1.0e-10

typedef double REAL;
typedef struct
{  REAL x, y, z;
}  VEC3F;

typedef VEC3F COOR;

#define ABS(x)   (((x)<0)   ? (-(x)) : (x))

/* Prototypes for functions defined in fit.c */
/* BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, 
            BOOL column);

static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column);  */


