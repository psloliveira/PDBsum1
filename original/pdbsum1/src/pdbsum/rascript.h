/***********************************************************************

  rascript.h - Include file for rascript.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN        20
#define COLNO_LEN         14
#define DISPLAY_MODE_LEN  21
#define FILENAME_LEN     512
#define LINELEN         1100
#define MAX_ATOM_COL      10
#define MAXATOMS       99999
#define MAX_CHAINS        50
#define MAXCONECT          5
#define MAXDOMCHAIN      200
#define MAXLIGANDS       100
#define MAX_LOCAL_PDB     27
#define MAX_METAL_COL      5
#define MAX_METAL_TYPES   50
#define MAX_RANGE_LEN     41
#define MAX_OBJ_COL        8
#define MAXSITE           50
#define MAX_VDW_RADII      6
#define OBJ_NAME_LEN      41
#define TOKEN_LEN       20

/* Record types defined in the rasmol.dfn file */
#define UNKNOWN           -1
#define CA_ONLY            0
#define LIGAND             1
#define METAL              2
#define NUCLEIC_ACID       3
#define PROTEIN            4
#define LIGAND_RANGE       5

/* Upper domain colour */
#define UPPER_DOMAIN_COL 7

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

/* Output options */
#define OUT_RASMOL        1
#define OUT_MOLSCRIPT     2
#define OUT_RASTER3D      3

/* Colour options */
#define COLOUR_BY_CHAIN   1
#define COLOUR_BY_SST     2

/* Secondary structures */
#define COIL              0
#define HELIX             1
#define STRAND            2
#define TURN              3
#define DHELIX            4

/* Special files */
#define STANDARD_PDB          0
#define UPLOADED_PDB          1
#define MCSG_PDB              2

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  int                colour_by_domain;
  int                all_chains;
  int                define_dots;
  int                define_ligands;
  int                define_mask;
  int                define_metals;
  int                define_nucleic;
  int                define_protein;
  int                define_sites;
  int                define_surfaces;
  int                show_disulphides;
  int                show_dots;
  int                show_ligands;
  int                show_mask;
  int                show_metals;
  int                show_nucleic;
  int                show_protein;
  int                show_sites;
  int                show_surfaces;
  int                display_option;
  int                write_to_file;
  int                write_gif_only;
  int                change_MSE;
  int                coords_only;
  int                output_option;
  char               background_colour[COLNAM_LEN];
  char               dhelix_colour[COLNAM_LEN];
  char               helix_colour[COLNAM_LEN];
  char               strand_colour[COLNAM_LEN];
  char               coil_colour[COLNAM_LEN];
  float              rotation_angle[3];
  int                revrot;
  int                have_transformations;
  int                colour_option;
  char               ligand_colour[COLNAM_LEN];
  float              metatom_radius;
  float              ligatom_radius;
  float              ligbond_radius;
  int                pdb_type;
  int                jmol;
  int                run_64;
};

/* Structure for each residue to be written out
   -------------------------------------------- */
struct residue
{
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                type;
  int                sec_struc;
  int                natoms;
  struct rasdef      *rasdef_ptr;
  struct atom        *first_atom_ptr;
  struct residue     *next_residue_ptr;
};

/* Structure for coordinate limits
   ------------------------------- */
struct limits
{
  int                     first;
  double                  valmin[3];
  double                  valmax[3];
  double                  c_of_g[3];
  int                     npoints;
};

/* Structure for each atom record written out
   ------------------------------------------ */
struct atom
{
  int                atom_number;
  char               atom_name[5];
  double             x;
  double             y;
  double             z;
  double             bvalue;
  double             occupancy;
  char               element[3];
  struct residue     *residue_ptr;
  struct atom        *next_atom_ptr;
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
  struct limits      atom_limits;
  struct atom        *first_atom_ptr;
  struct atom        *last_atom_ptr;
  struct rasdef      *next_rasdef_ptr;
};

/* Structure for each PROMOTIF data item
   ------------------------------------- */
struct promotif
{
  int                data_type;
  char               sec_struc;
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  int                class;
  struct promotif    *next_promotif_ptr;
};

/* Structure for each GROW data item
   ------------------------------------- */
struct grow
{
  char               data_type;
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  struct grow        *next_grow_ptr;
};

/* Structure for each CONECT record to be written out
   -------------------------------------------------- */
struct conect
{
  int                     atom_number[MAXCONECT];
  int                     nconect;
  struct conect           *next_conect_ptr;
};

