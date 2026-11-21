/***********************************************************************

  topoldia.h - Include file for topoldia.c

***********************************************************************/

/* Array parameters */
#define COLNAM_LEN        30
#define COLNAMLEN         12
#define COLNO_LEN         14
#define FILENAME_LEN     512
#define LINELEN         1100
#define MAX_DOMCOL         7
#define MAX_SEG            4
#define MAX_TOKEN         20
#define TOKEN_LEN         20

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

/* Motifs */
#define HELICES            0
#define HELIX_INTERACTIONS 1
#define STRANDS            2
#define BETA_TURNS         3
#define GAMMA_TURNS        4
#define DISULPHIDES        5

#define NMOTIFS            6

static const char *MOTIF_NAME[NMOTIFS] = {
  "HELICES", "HELIX_INTERACTIONS", "STRANDS", "BETA_TURNS", "GAMMA_TURNS",
  "DISULPHIDES"
};

/* Secondary structures */
#define NTERMINUS         0
#define COIL              1
#define HELIX             2
#define STRAND            3
#define TURN              4
#define DHELIX            5
#define CTERMINUS         6
#define OVERHANG1         7
#define OVERHANG2         8

#define NSSE_TYPES        9

static const char *SSE_TYPE[NSSE_TYPES] =
  { "NTERMINUS",
    "COIL     ",
    "HELIX    ",
    "STRAND   ",
    "TURN     ",
    "DHELIX   ",
    "CTERMINUS",
    "OVERHANG1",
    "OVERHANG2" };

/* Special files */
#define STANDARD_PDB          0
#define UPLOADED_PDB          1
#define MCSG_PDB              2

/* Conversion factors from HERA layout rows and columns to x-, y-coords */
#define CONVERTX              4
#define CONVERTY              2

/* Arrow orientations and directions */
#define HORIZONTAL            0
#define VERTICAL              1
#define FORWARD               0
#define BACKWARD              1

/* Curve-coil limit */
#define MIN_CURVE_DIST       10

/* Scores for picture cross-overs */
#define SSE_SCORE           100
#define LINE_SCORE            1
#define MAX_SCORE             3

/* Identifiers for path-search arrays */
#define STEP_DIR   0
#define STEP_SIZE  1

/* Axis identifiers */
#define X_AXIS    -1
#define Y_AXIS     1

/* Path-search step directions and axes */
static const int DIRECTION[4][2] = { 0, 1, 0, -1, 1, 0, -1, 0 };
static const int DIR_AXIS[4] = { X_AXIS, X_AXIS, Y_AXIS, Y_AXIS };

/* Domain colours */
static const char *SSE_COLOUR[MAX_DOMCOL + 1][2] =
  { "Red", "Gold",
    "Red", "HotPink",
    "Blue", "DeepSkyBlue",
    "Green", "LightGreen",
    "Purple", "SlateBlue",
    "SaddleBrown", "Orange",
    "Cyan", "SkyBlue",
    "Magenta", "MediumOrchid" };

/* PostScript parameters */
#define ARRAY_MARGIN                    10
#define ARROW_CUT             ((float)   1.0)
#define ARROW_SIZE            ((float)   2.0)
#define BBOXX1                ((float)  10.0)
#define BBOXX2                ((float) 580.0)
#define BBOXY1                ((float)  20.0)
#define BBOXY2                ((float) 790.0)
#define BORDER_MARGIN         ((float)   0.93)
#define CHAR_ASPECT           ((float)   0.484)
#define COIL_ARROW_WIDTH      ((float)   0.5)
#define COIL_ARROW_HEIGHT     ((float)   0.5)
#define COIL_OVERHANG         ((float)   1.0)
#define COIL_WIDTH            ((float)   0.2)
#define DEFAULT_WIDTH         ((float)  40.0)
#define HELIX_ELLIPSE         ((float)   2.0)
#define HELIX_WIDTH           ((float)   2.8)
#define LINE_WIDTH            ((float)   0.1)
#define PICTURE_HEIGHT        ((float) 700.0)
#define PICTURE_WIDTH         ((float) 500.0)
#define RES_NUMBERS           ((float)   2.0)
#define SHADOW_WIDTH          ((float)   0.6)
#define SSE_WIDTH                        3
#define TERMINUS_HEIGHT                  2
#define TERMINUS_WIDTH        ((float)   2.0)
#define TERMINUS_TEXT         ((float)   3.0)
#define TEXT_SIZE_MIN         ((float)   0.2)
#define XMARGIN               ((float)  30.0)
#define XWIDTH                ((float) 570.0)
#define YHEIGHT               ((float) 800.0)
#define YHEIGHT_LAND          ((float) 460.0)
#define YMARGIN               ((float)  50.0)

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  char               pdb_code[5];
  int                pdb_type;
  char               chain;
  int                domain;
  int                show_residue_nos;
  int                show_sse_labels;
  int                use_data;
  char               *plt_file;
  int                save_data;
  int                local;
  int                new;
  int                run_64;
};

/* Structure for each PROMOTIF data item
   ------------------------------------- */
struct promotif
{
  int                motif;
  char               sec_struc;
  int                identifier[2];
  char               res_name1[4];
  char               res_num1[6];
  char               chain1;
  char               res_name2[4];
  char               res_num2[6];
  char               chain2;
  struct promotif    *next_promotif_ptr;
};

/* Structure for each residue
   -------------------------- */
struct residue
{
  char               res_name[4];
  char               res_num[6];
  char               chain;
  int                sec_struc;
  struct residue     *next_residue_ptr;
};

/* Structure for each helix-helix interaction record
   ------------------------------------------------- */
struct helint
{
  float              distance;
  char               int_type1;
  char               int_type2;
  float              omega;
  int                nint;
  struct promotif    *promotif1_ptr;
  struct promotif    *promotif2_ptr;
  struct helint      *next_helint_ptr;
};

/* Structure for each SSE data item
   ------------------------------------- */
struct sse
{
  int                type;
  int                first_pos;
  int                last_pos;
  char               first_resnum[6];
  char               last_resnum[6];
  int                orientation;
  int                direction;
  int                xorig[2];
  int                yorig[2];
  int                xpos[2];
  int                ypos[2];
  int                xmin;
  int                xmax;
  int                ymin;
  int                ymax;
  int                endseg;
  int                sort;
  char               *sst_id;
  struct sse         *next_sse_ptr;
};

/* Structure for coordinate limits
   ------------------------------- */
struct limits
{
  int                     first;
  double                  min_x;
  double                  max_x;
  double                  min_y;
  double                  max_y;
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
  FILE               *tdf_file;
  int                page;
  float              Page_Max_x, Page_Min_x, Page_Max_y, Page_Min_y;
  float              origin[2];
  float              scale_factor;
  float              x, y;
  float              lim[4];
};

