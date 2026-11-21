/***********************************************************************

  dimps.h - Include file for dimps.c

***********************************************************************/

/* PostScript parameters */
#define ADD_CHAIN                FALSE
#define BBOXX1          ((float)  10.0)
#define BBOXX2          ((float) 580.0)
#define BBOXY1          ((float)  20.0)
#define BBOXY2          ((float) 790.0)
#define BORDER_MARGIN   ((float)   0.93)
#define CHAR_ASPECT     ((float)   0.484)
#define COIWID          ((float)   0.25)
#define DEFAULT_LINE_SEP ((float) 10.0)
#define ELLIPSE_ASPECT  ((float)   2.0)
#define ELLIPSE_SIZE    ((float)  30.0)
#define EXTRA_FACTOR    ((float)   1.2)
#define HEADTEXT_HEIGHT ((float)  30.0)
#define INTERELLIPSE_GAP ((float) 10.0)
#define INTERFACE_GAP   ((float)  80.0)
#define INTERLINE_WIDTH ((float)  1.5)
#define LABEL_GAP       ((float)  16.0)
#define MIN_HEAD_HEIGHT ((float)  20.0)
#define MIN_SCALED_SIZE ((float)  16.0)
/* v.5.3.4--> */
// #define NCOLOURS                  16
/* <--v.5.3.4 */
#define PS_STRING_LENGTH         132
#define PICTURE_HEIGHT  ((float) 700.0)
#define PICTURE_WIDTH   ((float) 500.0)
#define RES_BORDER_WIDTH ((float)  2.0)
#define STRIPE_GAP       ((float)  8.0)
#define STRIPE_WIDTH     ((float)  7.5)
#define TEXT_HEIGHT     ((float)  25.0)
#define XMARGIN         ((float)  30.0)
#define XWIDTH          ((float) 470.0)
#define YMARGIN         ((float)  50.0)
#define YHEIGHT         ((float) 800.0)
#define YHEIGHT_LAND    ((float) 460.0)

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

