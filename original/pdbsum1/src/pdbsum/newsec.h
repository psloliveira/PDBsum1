/***********************************************************************

  newsec.h - Include file for newsec.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN          20
#define FILENAME_LEN       512
#define LINELEN           1024
#define MAXCONECT            5
#define MAX_BIOMT          500
#define MAX_CHAINS         100
#define MAX_COPIES         100
#define MAX_DESC_LEN      1024
#define MAX_INTERFACES    1000
#define MAX_NUCPLOT_DNA_LEN 1000
#define MAX_LEN             16
#define MAX_LINES           10
#define MAX_OUTLEN          80
#define MAX_SEQLEN       10000
#define MAX_SEQ_HTML_LINE  600
#define MAX_SEQ_LINE       600
#define MAX_SURF_ATOMS   20000
#define MAX_SURF_CHAINS     10
#define MAX_SYMM_CHAINS   1000
#define MAX_SYMOP            6
#define MAX_TOKEN           10
#define NAMINO              20
#define NAMLEN             512
#define NSTRAND_VERTICES    14
#define TOKEN_LEN           20

/* Sequence comparison parameters */
#define MIN_SEQLEN          30
#define PERCEN_SEQLEN_DIFF ( 5.0)
#define TUPLEN               5

/* Raster3D flag */
#define RUN_RASTER3D      FALSE
#define MAX_PYMOL_ATOMS  100000

/* Image flags */
#define GIFS                 0
#define NO_GIFS              1

/* Flag indicating whether image files are GIF or JPEG format */
#define GIF                  0
#define JPEG                 1
#define IMAGE_FORMAT         JPEG

/* Commands set according to operating system */
char CP_COMMAND[6];
char MV_COMMAND[6];
char RM_COMMAND[6];
char SLEEP_COMMAND[30];

/* Ligand parameters */
#define LIGAND_LENGTH       10

/* Residue charges */
#define NOT_WANTED            0
#define POSITIVE              1                  
#define NEGATIVE              2                  

/* Distance parameters */
#define COVALENT_DIST      2.0
#define COVALENT_DIST_SQRD (COVALENT_DIST * COVALENT_DIST)
#define MAX_CONECT_LEN      2.5
#define MAX_CONECT_LEN_SQRD (MAX_CONECT_LEN * MAX_CONECT_LEN)
#define MAX_METAL_LEN       3.0
#define MAX_METAL_LEN_SQRD  (MAX_METAL_LEN * MAX_METAL_LEN)
#define SALT_BRIDGE_DIST    4.0
#define SBRIDGE_DIST_SQRD   (SALT_BRIDGE_DIST * SALT_BRIDGE_DIST)
#define TOO_FAR_DIST        5.0
#define TOO_FAR_CONTACT    10.0

/* Raster3D parameters */
#define BACKBONE_RADIUS     0.30
#define BACK_WALL_GAP       0.60
#define BOND_RADIUS         0.20
#define CA_CA_CUTOFF        6.00
#define CA_CA_CUTOFF_SQRD  (CA_CA_CUTOFF * CA_CA_CUTOFF)
#define HELIX_RADIUS        1.80
#define LIGAND_ATOM_PERCEN  0.25
#define LIGAND_BOND_RADIUS  0.20
#define MIN_ATOM_SIZE       0.015
#define N_HEIGHT            0.75
#define N_RADIUS            0.10
#define N_WIDTH             0.50
#define P_P_CUTOFF          8.00
#define P_P_CUTOFF_SQRD    (P_P_CUTOFF * P_P_CUTOFF)
#define RASTER3D_MIN       -0.47
#define RASTER3D_MAX        0.47
#define SPACEFILL_PERCEN    1.00
#define SSBOND_RADIUS       0.10
#define STRAND_HEAD_LEN     2.50
#define STRAND_HEAD_WIDTH   3.50
#define STRAND_THICKNESS    0.80
#define STRAND_WIDTH        1.80

/* Record types */
#define ATOM_REC              0
#define HETATM_REC            1

/* Link types */
#define ATOM_LINK             0
#define RESIDUE_LINK          1
#define CHAIN_LINK            2

/* Residue types */
#define UNKNOWN              -1
#define STANDARD              0
#define NON_STANDARD          1
#define NUCLEIC_ACID          2
#define METAL_ION             3

/* Chain types */
#define PROTEIN               0
#define DNA                   1
#define LIGAND                2
#define CALPHA_ONLY           3
#define METAL                 4
#define ATTACHMENT            5
#define PP_INTERFACE          6
#define PD_INTERFACE          7

#define NCHAIN_TYPES          8

static char *CHAIN_TYPE[NCHAIN_TYPES] = {
  "PROTEIN", "DNA/RNA", "LIGAND", "CA_ONLY", "METAL", "ATTACHMENT",
  "PP_INTERFACE", "PD_INTERFACE"
};

/* Quaternary structure types */
#define NOT_GIVEN             0

#define NUNIT_TYPES          12

static char *UNIT_TYPE[NUNIT_TYPES + 3] = {
  "NOT GIVEN", "MONOMER", "DIMER", "TRIMER", "TETRAMER", "PENTAMER",
  "HEXAMER", "HEPTAMER", "OCTAMER", "NONAMER", "DECAMER", "UNDECAMER",
  "DODECAMER", "COMPLEX", "VIRUS-COMPLEX"
};

/* Quaternary structure types */
#define COMPLEX         (NUNIT_TYPES + 1)
#define VIRUS_COMPLEX   (NUNIT_TYPES + 2)

#define NEXTRA_TYPES          7

static char *EXTRA_UNIT_TYPE[NEXTRA_TYPES] = {
  "TRIDECAMER", "TETRADECAMER", "PENTADECAMER", "HEXADECAMER",
  "HEPTADECAMER", "OCTADECAMER", "NONADECAMER"
};

static int EXTRA_NUMBER[NEXTRA_TYPES] = {
  30, 40, 50, 60, 70, 80, 90
};

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

/* Degrees of transparency in PyMol */
static const float PYMOL_TRANSPARENCY[3] = { 0.0, 0.3, 0.6 };

/* Interface chain equivalences */
#define NO                    0
#define DIRECT                1                  
#define ACROSS                2                  

/* Secondary structure types */
#define OTHER                -1
#define COIL                  0
#define HELIX                 1
#define STRAND                2

/* Parameters determining secondary structure dominance in a chain */
#define MIN_HELIX         ( 45.0 )
#define MIN_STRAND            6

/* Hetgroup types */
#define NHETYPE               3

#define HET_GROUP             0
#define AMINO_ACID            1
#define NON_STANDARD_AA       2
#define RESIDUE_ATTACHMENT    3

/* Hetdesc types */
#define DESCRIPTION           0
#define SYNONYM               1
#define FORMULA               2

/* Interaction types */
#define HBONDS                0
#define CONTACTS              1
#define SALT_BRIDGES          2
#define DISULPHIDES           3

#define NTYPES                4

/* Max sulphur-sulphur separation in disulphide */
#define MAX_SSDIST        ( 2.5 )
#define MAX_SSDIST2       ( MAX_SSDIST * MAX_SSDIST )

/* Interface types */
#define PROTEIN_PROTEIN       0
#define PROTEIN_DNA           1

/* Transformation types */
#define NO_MATRIX             0
#define SYMMETRY              1
#define BIOMOLECULE           2

/* Standard atom types and radii */
#define NATOM_TYPES           6
static char *ELEMENT[] = {
  "XX", " N", " C", " O", " S", " P"
};

static double VDW_RADIUS[NATOM_TYPES] = {
  1.80,  /*  0. Unknown */
  1.65,  /*  1. Nitrogen */
  1.87,  /*  2. Carbon */
  1.40,  /*  3. Oxygen */
  1.85,  /*  4. Sulphur */
  1.90,  /*  5. Phosphorus */
};

#define METAL_RADIUS        1.50

/* Corresponding atom colours */
static char *ATOM_COLOUR[NATOM_TYPES] = {
  "magenta", "blue", "black", "red", "yellow", "purple"
};

/* Metals */
#define NMETALS             100

/* Note "AC  " changed to "ac  " below as carbon atoms in certain
   molecules labelled AC1, AC1, etc (eg in FAD). Similarly NP, NO
   and PO */
static char *METAL_NAME[NMETALS] = {
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
#define NCOMMON_METALS        5
static char *COMMON_METAL[NCOMMON_METALS] = {
  "CL", "ZN", "CU", "CA", "FE"
};

#define LIGHT_GREY_RED      0.8
#define LIGHT_GREY_GREEN    0.8
#define LIGHT_GREY_BLUE     0.8

/* Metal colours */
#define NMETAL_COLOURS        5
static char *METAL_COLOUR[NMETAL_COLOURS] = {
  "green", "cyan", "orange", "pink", "brown"
};

/* SURFNET colours */
#define NSURFNET_COLS        10
static char *SURFNET_COLOUR[NSURFNET_COLS] = {
  "yellow", "cyan", "lime green", "green", "blue",
  "pink", "brown", "red", "purple", "black"
};

/* Special files */
#define STANDARD_PDB          0
#define UPLOADED_PDB          1
#define MCSG_PDB              2

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  char                    file_name[FILENAME_LEN + 1];
  int                     no_gifs;
  int                     uploaded_file;
  int                     pdb_type;
  int                     calc_interfaces;
  int                     gen_pymol;
  int                     molecule_align;
  int                     biomol_only;
  int                     run_raster3d;
  int                     ligand_length;
  int                     run_64;
};

/* Structure for info extracted from the PDB file
   ---------------------------------------------- */
struct data
{
  char                    pdb_code[5];
  double                  resolution;
  double                  rfactor;
  double                  rfree;
  int                     nmr;
  int                     theoretical_model;
  int                     ensemble;
  int                     nmodels;
  int                     nwaters;
  int                     nhydrogens;
  int                     nchains;
  int                     nprotein;
  int                     nca_only;
  int                     ndna;
  int                     nligands;
  int                     nmetals;
  int                     nattach;
  int                     nunknown;
  int                     nresidues;
  int                     natom_records;
  int                     nhetatm_records;
};

/* Structure for info about biomolecules and how to generate them
   -------------------------------------------------------------- */
struct biomol
{
  int                     first_matrix[MAX_BIOMT];
  int                     type[MAX_BIOMT];
  int                     is_biomolecule;
  int                     ntransform;
  int                     nsymmetry;
  int                     nbiomolecule;
  int                     nmolecules;
  float                   biomol_transform[MAX_BIOMT][3][3];
  int                     nmatrix;
  int                     nidentity;
  int                     new_style;
  int                     skip;
  int                     split_asu;
  char                    prefix[10];
  int                     quaternary_structure;
  char                    qstruc_desc[30];
  int                     nsymop;
  int                     in_pqs;
  int                     use_pqs;
  int                     pqs_skip;
  int                     pqs_split_asu;
  char                    pqs_prefix[10];
  int                     pqs_quaternary_structure;
  char                    pqs_qstruc_desc[30];
  int                     have_author_info;
  int                     author_skip;
  char                    author_prefix[10];
  int                     author_quaternary_structure;
  char                    author_qstruc_desc[30];
  char                    *tchains[MAX_BIOMT];
  int                     have_biomol_image;
};

/* Structure for each chain data item
   ---------------------------------- */
struct chain
{
  char                    chain_id;
  int                     offset;
  char                    name[LINELEN + 1];
  char                    *sequence;
  int                     type;
  int                     nresid;
  int                     nstandard;
  int                     natoms;
  int                     first_atom_number;
  int                     last_atom_number;
  int                     nhelices;
  int                     nstrands;
  int                     nstand;
  int                     nbase;
  char                    copy_of;
  int                     ncopies;
  int                     same_length;
  int                     natom_recs;
  int                     nhetatm_recs;
  double                  coord_min[3];
  double                  coord_max[3];
  int                     got_bonds;
  int                     got_sites;
  double                  asa;
  int                     nres_sst[3];
  int                     dominant_sst;
  int                     symm_chain;
  int                     used;
  struct chain            *link_chain_ptr;
  struct residue          *first_residue_ptr;
  struct chain            *next_chain_ptr;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  int                     number;
  int                     type;
  int                     ref;
  char                    aa_code;
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  int                     sec_struc;
  int                     natoms;
  int                     nmain_chain;
  int                     have_calpha;
  int                     connected;
  double                  coord_min[3];
  double                  coord_max[3];
  double                  ave_iface_coords[3];
  int                     niface_atoms;
  double                  centroid[3];
  int                     ncentroid_atoms;
  int                     got_bonds;
  char                    *site;
  int                     disulphide;
  int                     range;
  struct atom             *first_atom_ptr;
  struct chain            *chain_ptr;
  struct rbond            *first_rbond_ptr;
  struct hetdesc          *hetdesc_ptr;
  struct ligand           *ligand_ptr;
  struct residue          *link_residue_ptr;
  struct residue          *attach_residue_ptr;
  struct residue          *next_residue_ptr;
};

/* Structure for each atom record to be written out
   ------------------------------------------------ */
struct atom
{
  int                     atom_number;
  char                    atom_name[5];
  int                     rec_type;
  int                     type;
  char                    element[3];
  double                  coords[3];
  double                  original_coords[3];
  double                  view_coords[3];
  double                  bvalue;
  double                  occupancy;
  char                    alt;
  struct residue          *residue_ptr;
  struct conect           *first_conect_ptr;
  struct conect           *last_conect_ptr;
  struct atom             *link_atom_ptr;
  struct atom             *next_atom_ptr;
};

/* Structure for each atom CONECT record
   ------------------------------------- */
struct conect
{
  int                     from_conect;
  struct atom             *to_atom_ptr;
  struct conect           *next_conect_ptr;
};

/* Structure for each covalent bond between a pair of residues
   ----------------------------------------------------------- */
struct rbond
{
  struct residue          *to_residue_ptr;
  struct rbond            *next_rbond_ptr;
};

/* Structure for each disulphide bond between a pair of residues
   ------------------------------------------------------------- */
struct ssbond
{
  struct residue          *residue1_ptr;
  struct residue          *residue2_ptr;
  struct atom             *atom1_ptr;
  struct atom             *atom2_ptr;
  struct ssbond           *next_ssbond_ptr;
};

/* Structure for each active site residue
   -------------------------------------- */
struct site
{
  char                    *name;
  int                     number;
  char                    res_name[4];
  char                    res_num[6];
  char                    chain;
  struct site             *next_site_ptr;
};

/* Structure for each list data item
   --------------------------------- */
struct list
{
  struct residue          *residue_ptr;
  struct list             *next_list_ptr;
};

/* Structure for each ligand data item
   ---------------------------------- */
struct ligand
{
  int                     number;
  char                    *residue_sequence;
  char                    *residue_list;
  char                    *residue_range;
  int                     metal;
  int                     metal_colour;
  char                    metal_name[3];
  int                     type;
  int                     count;
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
  struct residue          *first_residue_ptr;
  struct list             *first_list_ptr;
  struct ligand           *duplicate_ligand_ptr;
  struct ligand           *next_ligand_ptr;
};

/* Structure for each hbadd data item
   ---------------------------------- */
struct hbadd
{
  struct residue          *residue_ptr;
  int                     nmissed;
  struct hbadd            *next_hbadd_ptr;
};

/* Structure for each hetgroup data item
   ------------------------------------- */
struct hetgroup
{
  int                     type;
  char                    name[4];
  int                     count;
  struct hetgroup         *next_hetgroup_ptr;
};

/* Structure for each hetdesc data item
   ------------------------------------ */
struct hetdesc
{
  char                    res_name[4];
  char                    *details[3];
  int                     count;
  struct hetdesc          *next_hetdesc_ptr;
};

/* Structure for each interface record
   ----------------------------------- */
struct interface
{
  int                     type;
  char                    chains[3];
  int                     nlinks;
  int                     unique;
  int                     ncopies;
  int                     nresid[2];
  int                     count[NTYPES];
  double                  chain_area[2];
  double                  interface_area[2];
  struct respair          *first_respair_ptr;
  struct respair          *last_respair_ptr;
  struct interface        *copyof_interface_ptr;
  struct interface        *next_interface_ptr;
};

/* Structure for each residue pair record
   -------------------------------------- */
struct respair
{
  struct residue          *residue_ptr[2];
  double                  ave_coords[2][3];
  int                     natoms[2];
  int                     count[NTYPES];
  struct respair          *next_respair_ptr;
};

/* Structure for each pymol data item
   ---------------------------------- */
struct pymol
{
  char                    *gif_name;
  char                    *png_name;
  int                     crop;
  struct pymol            *next_pymol_ptr;
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



