/***********************************************************************

  ppcont.h - Include file for ppcont.c

***********************************************************************/

/* Array parameters */
#define COLNAMLEN           12
#define FILENAME_LEN       512
#define LINELEN           5000
#define MAX_CHAIN_COL        8
#define MAX_TOKEN           40
#define NAMLEN             100
#define TOKEN_LEN           30

/* PostScript parameters */
#define BBOXX1                ((float)  10.0)
#define BBOXX2                ((float) 580.0)
#define BBOXY1                ((float)  20.0)
#define BBOXY2                ((float) 790.0)
#define BORDER_MARGIN         ((float)   0.93)
#define CHAINID1_HEIGHT       ((float)  30.0)
#define CHAINID2_HEIGHT       ((float)  10.0)
#define CHAR_ASPECT           ((float)   0.484)
#define CIRCLE_BORDER_WIDTH   ((float)   1.0)
#define COIWID                ((float)   0.25)
#define DEFAULT_WIDTH         ((float) 260.0)
#define HEADTEXT_HEIGHT       ((float)  15.0)
#define INNER_FRACTION        ((float)   0.80)
#define INTERCIRCLE_GAP       ((float)  18.0)
#define LABEL_HEIGHT          ((float)   8.0)
#define LINE_SEP              ((float)   4.0)
#define LINE_WIDTH            ((float)   2.0)
#define MAX_RADIUS            ((float)  70.0)
#define MIN_HEAD_HEIGHT       ((float)  20.0)
#define MIN_SCALED_SIZE       ((float)  16.0)
#define MIN_TEXT_SIZE         ((float)   8.0)
#define NPPCOLOURS                      16
#define PS_STRING_LENGTH               132
#define PICTURE_HEIGHT        ((float) 700.0)
#define PICTURE_WIDTH         ((float) 500.0)
#define TEXT_HEIGHT           ((float)  25.0)
#define WEDGE_LINE_WIDTH      ((float)   1.0)
#define XMARGIN               ((float)  30.0)
#define XWIDTH                ((float) 570.0)
#define YGAP                  ((float)  20.0)
#define YHEIGHT               ((float) 800.0)
#define YHEIGHT_LAND          ((float) 460.0)
#define YMARGIN               ((float)  50.0)

/* Angle constants */
#define FOUR_PI    ((float)  (4.0 * PI))

/* Interaction types */
#define HBONDS                0
#define CONTACTS              1
#define SALT_BRIDGES          2
#define DISULPHIDES           3

#define NTYPES                4

/* Interface angles */
#define CENTRE                0
#define START                 1
#define END                   2

/* Special files */
#define STANDARD_PDB          0
#define UPLOAD_PDB            1
#define MCSG_PDB              2

/* Random number generator parameters */
/* debug */
#define RANDOM_START       FALSE
/* debug */

/* Energy minimization parameters */
#define DELTA     ((float) 10.0)
#define MAXSTEPS           100

/* Energy minimization temrs */
#define NTERMS             4
#define CLASH_ENERGY       0
#define OVERLAP_ENERGY     1
#define SEPARATION_ENERGY  2
#define WEDGE_ENERGY       3

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char                    pdb_code[5];
  char                    chains[3];
  int                     pdb_type;
  char                    work_dir[FILENAME_LEN + 1];
  int                     all_chains;
  int                     read;
  int                     save;
  int                     run_64;
};

/* Structure for input PostScript parameters
   ----------------------------------------- */
struct psparams
{
  int                landscape;
  int                splitps;
  int                in_colour;
  char               program_name[FILENAME_LEN];
};

/* Structure for PostScript page info
   ---------------------------------- */
struct pspage
{
  char               ps_name[FILENAME_LEN];
  char               gif_name[FILENAME_LEN];
  FILE               *ps_file;
  int                page;
  float              Page_Max_x, Page_Min_x, Page_Max_y, Page_Min_y;
  float              origin[2];
  float              scale_factor;
  float              x, y;
  float              lim[4];
};

/* Structure for each interface record
   ----------------------------------- */
struct interface
{
  int                     nresid[2];
  double                  chain_area[2];
  double                  interface_area[2];
  int                     count[NTYPES];
  int                     ntypes;
  float                   angle;
  struct chain            *chain_ptr;
  struct interface        *next_interface_ptr;
};

/* Structure for each chain record
   ------------------------------- */
struct chain
{
  char                    chain;
  int                     number;
  int                     chain_no;
  double                  chain_area;
  double                  total_interface_area;
  int                     ninterfaces;
  int                     cluster;
  int                     done;
  struct interface        *first_interface_ptr;
  struct interface        *last_interface_ptr;
  struct wedge            *first_wedge_ptr;
  struct wedge            *last_wedge_ptr;
  struct chain            *next_chain_ptr;
  struct chain            *next_cluster_chain_ptr;
};

/* Structure for each wedge record
   ------------------------------- */
struct wedge
{
  float                   angle[3];
  int                     first_chain;
  struct interface        *interface_ptr;
  struct chain            *chain_ptr;
  struct chain            *other_chain_ptr;
  struct wedge            *next_wedge_ptr;
};

/* Structure for each cluster record
   --------------------------------- */
struct cluster
{
  int                     number;
  int                     nchains;
  struct chain            *first_chain_ptr;
  struct cluster          *next_cluster_ptr;
};

