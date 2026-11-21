/**************   P O S T S C R I P T   R O U T I N E S   *************/

#include "common.h"
#include "dimhtml.h"
#include "dimps.h"

/* Prototypes
   ---------- */

/* In common.c */
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
int string_truncate(char *string,int max_length);
void to_upper(char *search_string,int length);


/***********************************************************************

adjust_for_landscape  -  Perform boundary adjustments for printing
                         in Landscape, as opposed to Portrait,
                         orientation

***********************************************************************/

void adjust_for_landscape(struct pspage *pspage_ptr)
{
  /* Make necessary adjustments */
  pspage_ptr->Page_Min_x = BBOXY1 + XMARGIN;
  pspage_ptr->Page_Max_x = pspage_ptr->Page_Min_x + YHEIGHT;
  pspage_ptr->Page_Min_y = BBOXX1 + YMARGIN;
  pspage_ptr->Page_Max_y = pspage_ptr->Page_Min_y + YHEIGHT_LAND;
}
/***********************************************************************

initialise_psparams  -  Initialise PostScript parameters

***********************************************************************/

void initialise_psparams(struct psparams *psparams_ptr,
                         struct pspage *pspage_ptr)
{
  /* Initialise PostScript parameters */
  psparams_ptr->landscape = FALSE;
  psparams_ptr->splitps = FALSE;
  psparams_ptr->in_colour = TRUE;
  strcpy(psparams_ptr->program_name,"dimhtml");

  /* Initialise PostScript plot variables */
  pspage_ptr->ps_name[0] = '\0';
  pspage_ptr->ps_file = NULL;
  pspage_ptr->Page_Min_x = BBOXX1 + XMARGIN;
  pspage_ptr->Page_Max_x = pspage_ptr->Page_Min_x + XWIDTH;
  pspage_ptr->Page_Min_y = BBOXY1 + YMARGIN;
  pspage_ptr->Page_Max_y = pspage_ptr->Page_Min_y + YHEIGHT;
}
/***********************************************************************

init_pspage  -  Initialise PostScript page variables

***********************************************************************/

void init_pspage(struct psparams psparams_ptr,struct pspage *pspage_ptr)
{
  /* If PostScript plot to be in landscape mode, adjust relevant
     page parameters */
  if (psparams_ptr.landscape == TRUE)
    adjust_for_landscape(pspage_ptr);

  /* Initialise page variables */
  pspage_ptr->page = 0;

  /* Initialise starting y-position */
  pspage_ptr->y = pspage_ptr->Page_Max_y;

  /* Initialise page plot limits */
  pspage_ptr->lim[0] = XMARGIN + XWIDTH;
  pspage_ptr->lim[1] = 0.0;
  pspage_ptr->lim[2] = YMARGIN + YHEIGHT;
  pspage_ptr->lim[3] = 0.0;
}
/***********************************************************************

psatxt_   Write out text in Arial-Bold to PostScript file

***********************************************************************/

void psatxt_(FILE *ps_file,float x,float y,float size,char *text)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %8.3f Aprint\n",text,size);
}
/***********************************************************************

psbbox_ - Function to write out bounded box to Postscript file

**********************************************************************/

void psbbox_(FILE *ps_file,float x1,float y1,float x2,float y2,float x3,
             float y3,float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pline4\n",
          x1,y1,x2,y2,x3,y3,x4,y4);
}
/***********************************************************************

pscale_   Scale by given dimensions

***********************************************************************/

void pscale_(FILE *ps_file,float x_scale,float y_scale)
{
  fprintf(ps_file,"%7.2f %7.2f scale\n", x_scale, y_scale);
}
/***********************************************************************

psccol_  -  Define colour for filled circle

***********************************************************************/

void psccol_(FILE *ps_file,char *text_string)
{
    fprintf(ps_file,"/Circol { %s } def\n",text_string);
}
/***********************************************************************

pscirc_  -  Write out a coloured circle to PostScript file

***********************************************************************/

void pscirc_(FILE *ps_file,float x,float y,float radius)
{
    fprintf(ps_file,"%7.2f %7.2f %7.2f Circle\n",x,y,radius);
}
/***********************************************************************

pscolb_  -  Write background colour level (for text, etc)

***********************************************************************/

void pscolb_(FILE *ps_file,char *text_string)
{
    fprintf(ps_file,"%s\n",text_string);
}
/***********************************************************************

pscomm_   Write out a comment to the PostScript file

***********************************************************************/

void pscomm_(FILE *ps_file,char *text_string)
{
  fprintf(ps_file,"\n");
  fprintf(ps_file,"%% %s\n",text_string);
  fprintf(ps_file,"\n");
}
/***********************************************************************

pscrgb_  -  Write background colour level (for text, etc) using given
            RGB values

***********************************************************************/

void pscrgb_(FILE *ps_file,float red,float green,float blue)
{
  fprintf(ps_file,"%8.3f %8.3f %8.3f setrgbcolor\n",red,green,blue);
}
/***********************************************************************

psctxt_   Write out centred text to PostScript file

***********************************************************************/

void psctxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f Center\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f Print\n",text_string,size);
}
/***********************************************************************

pscxtx_   Write out text (centred on x only) to PostScript file

***********************************************************************/

void pscxtx_(FILE *ps_file,float x,float y,int size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %d Centerx\n",text_string,size);
  fprintf(ps_file,"(%s) %d Print\n",text_string,size);
}
/***********************************************************************

psdash_  -   Set line to dashed, dotted or continuous

***********************************************************************/

void psdash_(FILE *ps_file,int on_off)
{
  if (on_off == 0)
    fprintf(ps_file,"[]  0 setdash\n");
  else
    fprintf(ps_file,"[ %2d %2d ] 0 setdash\n",on_off,on_off);
}
/* v.5.3.8--> */
/***********************************************************************

psdots_  -   Set dotted line

***********************************************************************/

void psdots_(FILE *ps_file,int on,int off)
{
  fprintf(ps_file,"[ %2d %2d ] 0 setdash\n",on,off);
}
/* <--v.5.3.8 */
/***********************************************************************

pselip_  -  Write out a coloured ellipse to PostScript file

***********************************************************************/

void pselip_(FILE *ps_file,float x,float y,float radius,float elong_x,
             float elong_y)
{
    fprintf(ps_file,"0.0 0.0 %7.2f %7.2f %7.2f %7.2f %7.2f Ellipse\n",
            radius,elong_x,elong_y,x,y);
}
/* v.5.3.8--> */
/***********************************************************************

psoeli_  -  Write out an open ellipse to PostScript file

***********************************************************************/

void psoeli_(FILE *ps_file,float x,float y,float width,float ratio,
             float angle)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n",x,y);
  fprintf(ps_file,"%7.2f RotAngle\n",-angle);
  fprintf(ps_file,"0.0 0.0 %7.2f 1.0 %7.2f 0 0 Oellipse\n",width,ratio);
  fprintf(ps_file,"grestore\n");
}
/* <--v.5.3.8 */
/***********************************************************************

psltxt_ - Function to write out left aligned text to PostScript file

***********************************************************************/
void psltxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f  Lalign\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f  Bprint\n",text_string,size);
}
/***********************************************************************

psrtxt_ - Function to write out right aligned text to PostScript file

***********************************************************************/
void psrtxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f  Ralign\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f  Bprint\n",text_string,size);
}
/***********************************************************************

pshade_  -  Write colour/shade out to PostScript file

***********************************************************************/

void pshade_(FILE *ps_file,char *colour)
{
    fprintf(ps_file,"G %s\n",colour);
}
/***********************************************************************

pshcol_  -  Write colour/shade as RGB values out to PostScript file

***********************************************************************/

void pshcol_(FILE *ps_file,float rgb[3])
{
    fprintf(ps_file,"G %8.3f %8.3f %8.3f setrgbcolor\n",rgb[0],rgb[1],
            rgb[2]);
}
/***********************************************************************

psline_  -   Write line out to PostScript file

***********************************************************************/

void psline_(FILE *ps_file,float x1,float y1,float x2,float y2)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f L\n",x1, y1, x2, y2);
}
/***********************************************************************

pslwid_  -  Write line-width out to PostScript file

***********************************************************************/

void pslwid_(FILE *ps_file,float width)
{
    fprintf(ps_file,"%8.3f W\n",width);
}
/***********************************************************************

psrest_  -  Restore graphics state

***********************************************************************/

void psrest_(FILE *ps_file)
{
    fprintf(ps_file,"grestore\n");
}
/***********************************************************************

psrgbc_  -  Define RGB colour for filled circle

***********************************************************************/

void psrgbc_(FILE *ps_file,float red,float green,float blue)
{
  fprintf(ps_file,"/Circol { %8.3f %8.3f %8.3f setrgbcolor } def\n",
          red,green,blue);
}
/***********************************************************************

psrot_  -  Rotate graphics page through 90 degrees

***********************************************************************/

void psrot_(FILE *ps_file,float new_origin_x, float new_origin_y)
{
  fprintf(ps_file," %7.2f %7.2f moveto Rot90\n",new_origin_x,
          new_origin_y);
}
/***********************************************************************

pssave_  -  Save graphics state

***********************************************************************/

void pssave_(FILE *ps_file)
{
    fprintf(ps_file,"gsave\n");
}
/***********************************************************************

pstext_   Write out text to PostScript file

***********************************************************************/

void pstext_(FILE *ps_file,float x,float y,float size,char *text)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %8.3f Print\n",text,size);
}
/***********************************************************************

pstran_   Translate to given (x,y) position

***********************************************************************/

void pstran_(FILE *ps_file,float x,float y)
{
  fprintf(ps_file,"%7.2f %7.2f translate\n", x, y);
}
/***********************************************************************

pstxtg_  -  Write out Greek text to PostScript file

***********************************************************************/

void pstxtg_(FILE *ps_file,float x,float y,float size,char *text)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %8.3f Gprint\n",text,size);
}
/***********************************************************************

psubox_  -  Write out unbounded box to PostScript file

***********************************************************************/

void psubox_(FILE *ps_file,float x1,float y1,float x2,float y2,float x3,
             float y3,float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl4\n",
          x1, y1, x2, y2, x3, y3, x4, y4);
}
/***********************************************************************

psucir_  -  Write out an unbounded coloured circle to PostScript file

***********************************************************************/

void psucir_(FILE *ps_file,float x,float y,float radius)
{
    fprintf(ps_file,"%7.2f %7.2f %7.2f Ucircle\n",x,y,radius);
}
/***********************************************************************

psutri_  -  Write out unbounded triangle to PostScript file

***********************************************************************/

void psutri_(FILE *ps_file,float x1,float y1,float x2,float y2,float x3,
             float y3)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl3\n",
          x1, y1, x2, y2, x3, y3);
}
/***********************************************************************

get_col_rgb  -  Get the colour name for the given colour number and
                   the corresponding RGB values

***********************************************************************/

void get_col_rgb(int colour,char *c_name,float *rgb)
{
/* v.5.3.4--> 
  static char *colour_list[] = {
    "black", "red", "green", "blue", "yellow", "purple", "brown",
    "sky_blue", "orange", "cyan", "grey", "white", "magenta", "pink",
    "dark_grey", "redorange", "greenblue"
  };

  static float rgb_triplet[NCOLOURS + 1][3] = {
    (float) 0.000,  (float) 0.000,  (float) 0.000,
    (float) 1.000,  (float) 0.000,  (float) 0.000,
    (float) 0.100,  (float) 0.500,  (float) 0.000,
    (float) 0.000,  (float) 0.000,  (float) 1.000,
    (float) 1.000,  (float) 1.000,  (float) 0.000,
    (float) 0.800,  (float) 0.400,  (float) 1.000,
    (float) 0.800,  (float) 0.600,  (float) 0.000,
    (float) 0.000,  (float) 0.600,  (float) 1.000,
    (float) 1.000,  (float) 0.600,  (float) 0.000,
    (float) 0.200,  (float) 0.800,  (float) 1.000,
    (float) 0.500,  (float) 0.500,  (float) 0.500,
    (float) 1.000,  (float) 1.000,  (float) 1.000,
    (float) 0.938,  (float) 0.000,  (float) 0.500,
    (float) 1.000,  (float) 0.412,  (float) 0.706,
    (float) 0.200,  (float) 0.200,  (float) 0.200,
    (float) 1.000,  (float) 0.400,  (float) 0.000,
    (float) 0.400,  (float) 0.600,  (float) 1.000
  };
 <--v.5.3.4 */

  /* Get the name and RGB triplet for this colour */
/* v.5.3.4--> */
//  strcpy(c_name,colour_list[colour]);
//  rgb[0] = rgb_triplet[colour][0];
//  rgb[1] = rgb_triplet[colour][1];
//  rgb[2] = rgb_triplet[colour][2];
  strcpy(c_name,COLOUR_LIST[colour]);
  rgb[0] = RGB_TRIPLET[colour][0];
  rgb[1] = RGB_TRIPLET[colour][1];
  rgb[2] = RGB_TRIPLET[colour][2];
/* <--v.5.3.4 */
}
/***********************************************************************

psopen_  -  Open PostScript file and write out header records

***********************************************************************/

FILE *psopen_(char *filename,char *program_name)
{
  char  c_name[COLNAMLEN], colno[3], exclamation, percent;
  float grey_shade, rgb[3];
  int   icolour, igrey;
/* v.5.3.8--> */
  int i;
/* <--v.5.3.8 */

  FILE *ps_file;

  /* Initialise variables */
  percent = '%';
  exclamation = '!';

  /* Open the PostScript file */
  if ((ps_file = fopen(filename,"w")) == NULL)
    {
      printf("\n*** Unable to open PostScript file: %s\n",filename);
      exit(1);
    }

  /* Write out headings to PostScript file*/
  fprintf(ps_file,"%c%cPS-Adobe-3.0\n",percent,exclamation);
  fprintf(ps_file,"%%%%Creator: %s\n",program_name);
  fprintf(ps_file,"%%%%DocumentNeededResources: font Times-Roman Symbol\n");
  fprintf(ps_file,"%%%%BoundingBox: (atend)\n");
  // fprintf(ps_file,"%%%%Pages: 1\n");
  fprintf(ps_file,"%%%%EndComments\n");
  fprintf(ps_file,"%%%%BeginProlog\n");
  fprintf(ps_file,"/L { moveto lineto stroke } bind def\n");
  fprintf(ps_file,"/G { gsave } bind def\n");
  fprintf(ps_file,"/W { setlinewidth } bind def\n");
  fprintf(ps_file,"/D { setdash } bind def\n");
  fprintf(ps_file,"/Col { setrgbcolor } bind def\n");
  fprintf(ps_file,"/Zero_linewidth { 0.0 } def\n");
  fprintf(ps_file,"/Sphcol { 1 setgray } def\n");
  fprintf(ps_file,"/Sphere {");
  fprintf(ps_file,"  newpath 3 copy 0 360 arc gsave Sphcol fill 0\n");
  fprintf(ps_file,"  setgray 0.5 setlinewidth\n");
  fprintf(ps_file,"  3 copy 0.94 mul 260 350 arc stroke 3 copy 0.87");
  fprintf(ps_file," mul 275 335 arc stroke\n");
  fprintf(ps_file,"  3 copy 0.79 mul 295 315 arc stroke 3 copy 0.8");
  fprintf(ps_file," mul 115 135 arc\n");
  fprintf(ps_file,"  3 copy 0.6 mul 135 115 arcn closepath gsave 1");
  fprintf(ps_file," setgray fill grestore stroke\n");
  fprintf(ps_file,"  3 copy 0.7 mul 115 135 arc stroke 3 copy 0.6");
  fprintf(ps_file," mul 124.9 125 arc\n");
  fprintf(ps_file,"  0.8 mul 125 125.1 arc stroke grestore stroke");
  fprintf(ps_file," } bind def\n");

  /* Write out the different colours */
  for (icolour = 0; icolour < NCOLOURS; icolour++)
    {
      /* Get the colour name and RGB triplet for this colour */
      get_col_rgb(icolour,c_name,rgb);

      /* Write out the colour definition */
      fprintf(ps_file,"/%s    { %5.3f %5.3f %5.3f  setrgbcolor } def\n",
              c_name,rgb[0],rgb[1],rgb[2]);
    }

  /* Loop to write out the different grey-shades */
  for (igrey = 0; igrey < 11; igrey++)
  {
      grey_shade = (float)igrey / (float) 10.0;
      sprintf(colno,"%2d",igrey);
      colno[2] = '\0';

      if (colno[0] == ' ')
          colno[0] = '0';
      fprintf(ps_file,"/Grey%s { %6.2f setgray } def\n",colno,grey_shade);
  }

  /* Print the remainder of the definitions */
  fprintf(ps_file,"/Poly3 { moveto lineto lineto fill grestore } bind ");
  fprintf(ps_file,"def\n");
  fprintf(ps_file,"/Pl3 { 6 copy Poly3 moveto moveto moveto closepath ");
  fprintf(ps_file,"stroke } bind def\n");
  fprintf(ps_file,"/Pline3 { 6 copy Poly3 moveto lineto lineto closepa");
  fprintf(ps_file,"th stroke } bind def\n");
  fprintf(ps_file,"/Poly4 { moveto lineto lineto lineto fill grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Pl4 { 8 copy Poly4 moveto moveto moveto moveto ");
  fprintf(ps_file,"closepath stroke } bind def\n");
  fprintf(ps_file,"/Pline4 { 8 copy Poly4 moveto lineto lineto lineto");
  fprintf(ps_file," closepath stroke } bind def\n");
  fprintf(ps_file,"/Circol { 1 setgray } def\n");
  fprintf(ps_file,"/Circle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore stroke grestore } bind def\n");
  fprintf(ps_file,"/Ucircle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore grestore } bind def\n");
  fprintf(ps_file,"/Arc { newpath arc stroke newpath } bind def\n");
  fprintf(ps_file,"/Ellipse { gsave translate scale Circle grestore } ");
  fprintf(ps_file,"bind def\n");
/* v.5.3.8--> */
  fprintf(ps_file,"/Ocircle { gsave newpath 0 360 arc stroke grestore } ");
  fprintf(ps_file,"bind def\n");

  /* Loop over the ellipse colours used for highlighting variants */
  for (i = 0; i < NVAR_TYPES + 2; i++)
    fprintf(ps_file,"/Ellipse_col_%s { %s setrgbcolor } def\n",
            ELLIPSE_COL_DEF[i],ELLIPSE_COLOUR[i]);
  fprintf(ps_file,"/Ellipse { gsave translate scale Circle grestore } ");
  fprintf(ps_file,"bind def\n");
  fprintf(ps_file,"/Oellipse { gsave translate scale Ocircle grestore } ");
  fprintf(ps_file,"bind def\n");
/* <--v.5.3.8 */
  fprintf(ps_file,"/Print { /Times-Roman findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Bprint { /Times-Bold findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Aprint { /Arial-Bold findfont exch scalefont setfont");
  fprintf(ps_file," show } bind def\n");
  fprintf(ps_file,"/Gprint { /Symbol findfont exch scalefont setfont show");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Lalign {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop pop -3 div 0.0 exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Ralign {\n");
  fprintf(ps_file,"  dup /Times-Bold findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -1 mul exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Center {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch -3 div rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/Centerx {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch 0 mul rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/CenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch 3 div exch rmoveto\n");
  fprintf(ps_file," } bind def\n");
  fprintf(ps_file,"/UncenterRot90 {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth } bind def\n");
  fprintf(ps_file,"/Rot90 { gsave currentpoint translate 90 rotate }");
  fprintf(ps_file," bind def\n");
/* v.5.3.8--> */
  fprintf(ps_file,"/RotAngle { gsave currentpoint translate rotate } "
          "bind def\n");
/* <--v.5.3.8 */
  fprintf(ps_file,"%%%%EndProlog\n");
  fprintf(ps_file,"%%%%BeginSetup\n");
  fprintf(ps_file,"1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray");
  fprintf(ps_file," [ ] 0 setdash newpath\n");
  fprintf(ps_file,"%%%%EndSetup\n");

  /* Return file pointer */
  return(ps_file);
}
/***********************************************************************

pspage_  -  Write out start of new PostScript page

***********************************************************************/

void pspage_(FILE *ps_file,int page,int in_colour,int landscape,
             char *program_name)
{
  float xx1, xx2, xy1, xy2;

  /* Define border round picture */
  xx1 = (float)BBOXX1;
  xx2 = (float)BBOXX2;
  xy1 = (float)BBOXY1;
  xy2 = (float)BBOXY2;

  /* No border round picture */
  xx1 = -1.0;
  xx2 = 650.0;
  xy1 = -1.0;
  xy2 = 951.0;

  /* Print the page-opening commands */
  fprintf(ps_file,"%%%%Page: %d %d \n",page,page);
  fprintf(ps_file,"/%sSave save def\n",program_name);
  //fprintf(ps_file,"%6.2f%6.2f moveto %7.2f%7.2f lineto%7.2f%7.2f lineto\n",
  //        xx1,xy1,xx2,xy1,xx2,xy2);
  //fprintf(ps_file,"%7.2f%7.2f lineto closepath\n",xx1,xy2);
  //fprintf(ps_file,"gsave 1.0000 setgray fill grestore\n");
  //fprintf(ps_file,"stroke gsave\n");
  fprintf(ps_file,"gsave\n");

  /*Draw in the background colour*/
  //  if (in_colour == TRUE)
  //  {
  //    pscomm_(ps_file,"Background");
  //    pshade_(ps_file,"white");
  //    psubox_(ps_file,BBOXX1, BBOXY1, BBOXX2, BBOXY1, BBOXX2, BBOXY2, BBOXX1, BBOXY2);
  //  }

  /* If picture to be plotted in landscape mode, then rotate page */
  if (landscape == TRUE)
    psrot_(ps_file,BBOXX2 + BBOXX1,0.0);
}
/***********************************************************************

psendp_  -  Close off the current PostScript page

***********************************************************************/

void psendp_(FILE *ps_file,float lim[4],char *program_name)
{
  /* Write out the clip command (only works for Portrait mode) */
  fprintf(ps_file,"%%convert -clip %dx%d+%d-%d\n",
          (int) (lim[1] - lim[0] + 10),(int) (lim[3] - lim[2] + 10),
          (int) lim[0],(int) lim[2]);

  fprintf(ps_file,"\n\n");
  fprintf(ps_file,"%sSave restore\n",program_name);
  fprintf(ps_file,"showpage\n");

  /* Re-initialise limits */
  lim[1] = 0.0;
  lim[2] = YMARGIN + YHEIGHT;
}
/***********************************************************************

psclos_   Write final lines to PostScript file and close

***********************************************************************/

void psclos_(FILE *ps_file)
{
  /* Write out closing lines to PostScript file*/
  fprintf(ps_file,"%%%%Trailer\n");
  fprintf(ps_file,"%%%%BoundingBox: %d %d %d %d\n",(int) BBOXX1 - 1,
          (int) BBOXY1 - 1, (int) BBOXX2 + 1, (int) BBOXY2 + 1);
  fprintf(ps_file,"%%%%EOF");
  fclose(ps_file);
}
/***********************************************************************

start_new_page  -  Start a new PostScript page

***********************************************************************/

void start_new_page(struct psparams psparams_ptr,
                    struct pspage *pspage_ptr)
{
  char filename[FILENAME_LEN];

  /* If have a previous page, then close off */
  if (pspage_ptr->page > 0)
    psendp_(pspage_ptr->ps_file,pspage_ptr->lim,
            psparams_ptr.program_name);

  /* If writing each page to a separate PostScript file, then
     close previous file and form name of next one */
  if (psparams_ptr.splitps == TRUE)
    {
      /* Close off previous file, if there is one */
      if (pspage_ptr->page > 0)
        psclos_(pspage_ptr->ps_file);

      /* Form name of next file */
      sprintf(filename,"%s%4.4d.ps",pspage_ptr->ps_name,
              pspage_ptr->page + 1);
    }
  else
    strcpy(filename,pspage_ptr->ps_name);

  /* If this is the first page, or printing every page to a new
     file, then open the PostScript file */
  if (pspage_ptr->page == 0 || psparams_ptr.splitps == TRUE)
    {
      /* Open the file */
      pspage_ptr->ps_file = psopen_(filename,psparams_ptr.program_name);
    }

  /* Increment page count */
  pspage_ptr->page++;

  /* Initialise the new PostScript page */
  pspage_(pspage_ptr->ps_file,pspage_ptr->page,psparams_ptr.in_colour,
          psparams_ptr.landscape,psparams_ptr.program_name);

  /* Write out the headings for the plot */
          
  /* Initialise all the variables for the current page */
  pspage_ptr->y = pspage_ptr->Page_Max_y;
}
/***********************************************************************

get_chain_col  -  Determine the colour for the given chain

***********************************************************************/

int get_chain_col(char chain,char *colour_name)
{
  int colour, ichain, nchains;

  static char chain_list[]
    = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890!#_-=:<>|"
    "abcdefghijklmnopqrstuvwxyz";
  static char *col_name[] = {
    "purple", "red", "brown", "pink", "blue",
    "green", "cyan", "yellow"
  };

  /* Initialise variable */
  colour = 0;
  colour_name[0] = '\0';

  /* Get number of chain types */
  nchains = strlen(chain_list);

  /* Get the chain's position in the list */
  for (ichain = 0 ; ichain < nchains && colour == 0; ichain++)
    {
      if (chain == chain_list[ichain])
        {
          /* Store the colour number */
          colour = ichain;

          /* Get the corresponding colour name */
          colour = colour % MAX_CHAIN_COL;
          strcpy(colour_name,col_name[colour]);
        }
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

get_domain_colour  -  Get the colour of the given residue, according to
                      its CATH domain assignment

***********************************************************************/

void get_domain_colour(int domain,char *colour_name)
{
  int c_numb;

  static char *col_name[] = { 
    "black", "red", "blue", "green", "purple",
    "orange", "cyan", "pink" };

  /* Initialise colour */
  c_numb = 0;

  /* Find the colour corresponding to the given domain */
  if (domain < 0)
    c_numb = 0;
  else if (domain > UPPER_DOMAIN_COL)
    c_numb = domain % UPPER_DOMAIN_COL + 1;
  else
    c_numb = domain;

  /* Get the corresponding colour name */
  strcpy(colour_name,col_name[c_numb]);
}
/***********************************************************************

get_angle  -  Function to calculate the angle between line from origin
              to supplied point (x,y), and the x-axis. The value
              returned is between 0 and 2 pi radians.

***********************************************************************/

float get_angle(float x, float y)
{
  float ang;

  /* if x is very close to zero, then return pi/2 or 3pi/2, depending
     on y-value */
  if (fabs(x) < 0.0000001)
    {
      if (y < 0.0)
        ang = (float) (3.0 * PI / 2.0);
      else
        ang =       (float) (PI / 2.0);
    }

  /* If y is very close to zero, then return 0 or pi, depending on
     x-value */
  else if (fabs(y) < 0.0000001)
    {
      if (x < 0.0)
        ang = PI;
      else
        ang = 0.0;
    }

  /* Otherwise, calculate angle from arctan y/x */
  else
    {
      ang = (float) atan((double) (y / x));
      if (ang > 0.0)
        {
          if (y < 0.0) 
            ang = ang + PI;
        }
      else
        {
          if (y < 0.0)
            ang = (float) (2.0 * PI + ang);
          else
            ang = PI + ang;
        }
    }
  return ang;
}
/***********************************************************************

draw_interactions  -  Draw the lines representing the interactions
                      between the given two residues

***********************************************************************/

void draw_interactions(float scale_factor,float xleft,float yleft,
                       float xright,float yright,int count[2],
                       int disulphide,int salt_bridge,float line_sep,
                       float stripe_width,char *comment,
                       struct pspage *pspage_ptr)
{
  static char *colour[] = { "sky_blue", "orange", "red", "yellow" };

  int flip, istripe, line, nlines, nstripes;
  int col;

  static int last_col = -1;

  float dx, dy, length, stripe_sep, theta, width, x, x1, x2, y, y1, y2;
  float xstep, ystep;

  /* Initialise variables */
  pscomm_(pspage_ptr->ps_file,comment);

  /* Calculate angle made by the line */
  x = xright - xleft;
  y = yright - yleft;
  theta = get_angle(x,y);

  /* Plot stripes indicating non-bonded contacts */
  if (count[CONTACTS] > 0)
    {
      /* Set colour */
      col = CONTACTS;

      /* If colour has changed, define line colour */
      if (col != last_col)
        {
          /* Define line colour */
          pscolb_(pspage_ptr->ps_file,colour[col]);

          /* Save current colour */
          last_col = col;
        }

      /* Determine width of the stripes */
      width = stripe_width * count[CONTACTS];

      /* Calculate adjustment to each stripe point */
      dx = width * (float) sin(theta) / (float) 2.0;
      dy = width * (float) cos(theta) / (float) 2.0;

      /* Get the length of the line */
      length = (xright - xleft) * (xright - xleft)
       + (yright - yleft) * (yright - yleft);
      length = (float) sqrt(length);

      /* Determine the number of stripe points to be plotted */
      stripe_sep = scale_factor * STRIPE_GAP;
      nstripes = (int) (length / stripe_sep) + 1;

      /* Calculate steps between stripe centre-points */
      xstep = (xright - xleft) / (nstripes - 1);
      ystep = (yright - yleft) / (nstripes - 1);

      /* Draw the stripes */
      for (istripe = 0; istripe < nstripes; istripe++)
       {
         /* Calculate start- and end-points of this stripe */
         x1 = xleft + istripe * xstep + dx;
         x2 = xleft + istripe * xstep - dx;
         y1 = yleft + istripe * ystep - dy;
         y2 = yleft + istripe * ystep + dy;

         /* Draw the stripe */
         psline_(pspage_ptr->ps_file,x1,y1,x2,y2);
       }
    }

  /* If have disulphide bond or salt bridge, draw extra lines */
  nlines = count[HBONDS];
  if (disulphide == TRUE)
    nlines++;
  if (salt_bridge == TRUE)
    nlines++;

  /* Initialise the flip factor (zero if an odd number of lines and
     1 if an even number) */
  flip = 0;
  if (nlines % 2 == 0)
    flip = 1;

  /* Loop over the H-bonds to be plotted */
  col = HBONDS;
  for (line = 0; line < nlines; line++)
    {
      /* If first line, check whether it is a special one */
      if (line == 0 && disulphide == TRUE)
        col = DISULPHIDES;
      else if (line == 0 && salt_bridge == TRUE)
        col = SALT_BRIDGES;
      else
        col = HBONDS;

      /* If colour has changed, define line colour */
      if (col != last_col)
        {
          /* Define line colour */
          pscolb_(pspage_ptr->ps_file,colour[col]);

          /* Save current colour */
          last_col = col;
        }

      /* Calculate coords of this line */
      dx = flip * line_sep * (float) sin(theta) / 2;
      dy = flip * line_sep * (float) cos(theta) / 2;
      x1 = xleft + dx;
      x2 = xright + dx;
      y1 = yleft + dy;
      y2 = yright + dy;

      /* Draw line joining the two residues */
      psline_(pspage_ptr->ps_file,x1,y1,x2,y2);

      /* Adjust the flip flag */
      if (flip <= 0)
        flip = abs(flip) + 2;
      else
        flip = - flip;
    }
}
/***********************************************************************

plot_interactions  -  Plot the interactions between residues either side
                      of the interface

***********************************************************************/

void plot_interactions(struct residue **alignment,int align_len,
                       int total_len,float scale_factor,float origin[2],
                       float line_sep,float stripe_width,
                       struct pspage *pspage_ptr)
{
  char comment[80];

  int ipos, jpos;
  int found;

  float line_width, radius, xleft, xright, yleft, yright;

  struct link *link_ptr;
  struct residue *residue2_ptr;
  struct residue *residue_ptr,  *other_residue_ptr;

  /* Initialise variables */
  pscomm_(pspage_ptr->ps_file,"Interactions");

  /* Define colours and line widths */
  line_width = scale_factor * INTERLINE_WIDTH;
/* v.5.4.4 --> 
  if (line_width < 0.75)
    line_width = (float) 0.75;
 <-- v.5.4.4 */
  pscolb_(pspage_ptr->ps_file,"black");
  pslwid_(pspage_ptr->ps_file,line_width);

  /* Define starting location */
  xleft = origin[0];
  xright = xleft + scale_factor
    * (2 * ELLIPSE_SIZE * ELLIPSE_ASPECT + INTERFACE_GAP);
  yleft = origin[1];
  radius = scale_factor * ELLIPSE_SIZE;

  /* Loop through all the positions in the alignment to fill any blanks */
  for (ipos = 0; ipos < total_len; ipos++)
    {
      /* Get the first residue pointer */
      residue_ptr = alignment[ipos];

      /* If present, plot all the interactions this residue makes */
      if (residue_ptr != NULL)
        {
          /* Get the first of this residue's interactions */
          link_ptr = residue_ptr->first_link_ptr;

          /* Loop through all the links to see if first residue
             interacts with the second */
          while (link_ptr != NULL)
            {
              /* Get other residue */
              other_residue_ptr = link_ptr->residue_ptr[1];

              /* Initialise flag and other residue's y-value */
              found = FALSE;
              yright = origin[1];
              
              /* Loop over all positions in the alignment until
                 current interacting residue found */
              for (jpos = 0; jpos < total_len && found == FALSE; jpos++)
                {
                  /* Get the other residue pointers at this point in
                     the alignment */
                  residue2_ptr = alignment[jpos + align_len];

                  /* If this is the residue we're after, can plot the
                     interaction(s) */
                  if (residue2_ptr == other_residue_ptr)
                    {
                      /* Form comment to be written to PostScript file */
                      sprintf(comment,"%s %s %c <-> %s %s %c [%d %d",
                              residue_ptr->res_name,
                              residue_ptr->res_num,
                              residue_ptr->chain,
                              residue2_ptr->res_name,
                              residue2_ptr->res_num,
                              residue2_ptr->chain,
                              link_ptr->count[HBONDS],
                              link_ptr->count[CONTACTS]);
                      if (link_ptr->disulphide == TRUE)
                        strcat(comment," d");
                      if (link_ptr->salt_bridge == TRUE)
                        strcat(comment," s");
                      strcat(comment,"]");

                      /* Draw the interactions */
                      draw_interactions(scale_factor,xleft,yleft,xright,
                                        yright,link_ptr->count,
                                        link_ptr->disulphide,
                                        link_ptr->salt_bridge,line_sep,
                                        stripe_width,comment,pspage_ptr);

                      /* Set flag that residue found and that we can
                         move on */
                      found = TRUE;
                    }

                  /* Decrement y-position of right-hand residue */
                  yright = yright
                    - scale_factor * (2 * ELLIPSE_SIZE + INTERELLIPSE_GAP);
                }

              /* Get the next link */
              link_ptr = link_ptr->nextres_link_ptr[0];
            }
        }

      /* Decrement y-value */
      yleft = yleft - scale_factor * (2 * ELLIPSE_SIZE + INTERELLIPSE_GAP);
    }
}
/* v.5.3.8--> */
/***********************************************************************

find_highres_entry  -  Find the given record

***********************************************************************/

struct highres *find_highres_entry(struct highres *first_highres_ptr,
                                   char *res_name,char *res_num,
                                   char chain)
{
  struct highres *highres_ptr, *found_highres_ptr;

  /* Initialise variables */
  found_highres_ptr = NULL;

  /* Get pointer to the first record */
  highres_ptr = first_highres_ptr;

  /* Loop over all the records to find the one required */
  while (highres_ptr != NULL && found_highres_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (chain == highres_ptr->chain &&
          !strcmp(res_num,highres_ptr->res_num) &&
          !strcmp(res_name,highres_ptr->res_name))
        found_highres_ptr = highres_ptr;

      /* Get pointer to the next domain */
      highres_ptr = highres_ptr->next_highres_ptr;
    }

  /* Return the matching record pointer */
  return(found_highres_ptr);
}
/***********************************************************************

create_highres_entry  -  Create a record holding residue to be
                         highlighted

***********************************************************************/

struct highres *create_highres_entry(struct highres **fst_highres_ptr,
                                     struct highres **lst_highres_ptr,
                                     int type,int ref,char *res_name,
                                     char *res_num,char chain)
{
  struct highres *highres_ptr;
  struct highres *first_highres_ptr, *last_highres_ptr;

  /* Initialise pfam pointers */
  highres_ptr = NULL;
  first_highres_ptr = *fst_highres_ptr;
  last_highres_ptr = *lst_highres_ptr;

  /* Create pfam entry */
  highres_ptr = (struct highres *) malloc(sizeof(struct highres));
  if (highres_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct highres\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the known data */
  highres_ptr->type = type;
  highres_ptr->ref = ref;
  highres_ptr->chain = chain;
  strcpy(highres_ptr->res_name,res_name);
  strcpy(highres_ptr->res_num,res_num);

  /* Initialise pointer to next entry */
  highres_ptr->next_highres_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_highres_ptr == NULL)
    first_highres_ptr = highres_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_highres_ptr->next_highres_ptr = highres_ptr;
                      
  /* Save the current residue's pointer */
  last_highres_ptr = highres_ptr;

  /* Return current pointers */
  *fst_highres_ptr = first_highres_ptr;
  *lst_highres_ptr = last_highres_ptr;
  return(highres_ptr);
}
/***********************************************************************

format_res_number  -  Interpret the residue-number from the input
                      parameter for the residue to be highlighted

***********************************************************************/

void format_res_number(char *tmpstr,char *res_num)
{
  char ch, insertion_code, new_number[6], res_number[6];
  int cpos, idigit, ipos, ndigits;
  int nchars;
  int at_end, got_number;

  /* Initialise variables */
  at_end = FALSE;
  got_number = FALSE;
  nchars = 0;
  ndigits = 0;
  strcpy(new_number,"     ");
  strcpy(res_number,"     ");

  /* Initialise position of insertion code */
  cpos = 0;
  insertion_code = ' ';

  /* Loop through the characters searching for end of number or
     insertion-code */
  for (ipos = 0; ipos < 5 && at_end == FALSE; ipos++)
    {
      /* Get this character */
      ch = tmpstr[ipos];

      /* If this is the end of the string, store end-position */
      if (ch == '\0')
        {
          cpos = ipos;
          at_end = TRUE;
        }

      /* If this is a digit, then have part of a residue number */
      else if ((ch >= '0' && ch <= '9') ||
               ch == '-')
        {
          res_number[ndigits] = ch;
          ndigits++;
          got_number = TRUE;
        }

      /* If this is a character it might be the insertion code */
      else if ((ch >= 'A' && ch <= 'Z') ||
               ch == '\'')
        {
          /* Store position of insertion code */
          cpos = ipos;
          insertion_code = ch;

          /* Increment count of characters */
          nchars++;
        }

      /* If this is a space, have reached the end of the number */
      else if (ch == ' ' || ch == '\0')
        {
          if (got_number == TRUE)
            cpos = ipos;
        }

      /* Otherwise, have an invalid character in the residue
         number */
      else
        cpos = -1;
    }

  /* If have too many digits, then last one is an insertion code */
  if (ndigits > 4)
    {
      insertion_code = res_number[4];
      ndigits--;
    }

  /* Reformat the number */
  if (cpos == 0)
    cpos = 4;

  /* Now re-form the number such the first four characters
     contain the digits and the 5th contains the insertion
     code (if any) */
  for (idigit = 0; idigit < ndigits; idigit++)
    {
      ipos = idigit + (4 - ndigits);
      new_number[ipos] = res_number[idigit];
    }

  /* Add the insertion code */
  new_number[4] = insertion_code;
  new_number[5] = '\0';

  /* Transfer the new number across */
  strncpy(res_num,new_number,5);
  res_num[5] = '\0';
}
/***********************************************************************

extract_residue_name  -  Extract resisdue name, number and chain id

***********************************************************************/

void extract_residue_name(char *highlight_res,char *res_name,
                          char *res_num,char *chn)
{
  char chain, number_string[20], tmp_string[20];

  char *string_ptr;

  int name_end;

  /* Initialise variables */
  res_name[0] = res_num[0] = '\0';
  chain = ' ';
  
  /* Check for square brackets defining residue name */
  if (highlight_res[0] == '[')
    {
      /* Get the name bounded by square brackets */
      strcpy(tmp_string,highlight_res + 1);

      /* Check for close-bracket */
      string_ptr = strstr(tmp_string,"]");
      if (string_ptr != NULL)
        {
          string_ptr[0] = '\0';
          name_end = (string_ptr - tmp_string) + 1;
        }
      else
        {
          tmp_string[3] = '\0';
          name_end = 3;
        }

      /* Store the residue name */
      strcpy(res_name,tmp_string);
    }

  /* Otherwise, take residue name to be first 3 characters */
  else
    {
      /* Transfer name */
      strncpy(res_name,highlight_res,3);
      res_name[3] = '\0';
      name_end = 2;
    }

  /* Convert residue name to upper case */
  to_upper(res_name,3);

  /* Get the residue number */
  strcpy(number_string,highlight_res + name_end + 1);
  string_ptr = strstr(number_string,":");
  if (string_ptr != NULL)
    {
      /* Extract the chain identifier after the colon */
      chain = string_ptr[1];
      string_ptr[0] = '\0';
    }

  /* Convert the residue number into the standard fixed format */
  format_res_number(number_string,res_num);

  /* Return the chain id */
  *chn = chain;
}
/***********************************************************************

get_other_chain  -  Read in the chain equivalences from the chains.dat
                   file and return the other chain in the interface if
                   it is equivalent to the one with the highlighted
                   residue

***********************************************************************/

char get_other_chain(char chain,char *chains,char *pdbsum_dir)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char chain_id, copy_of, other_chain;

  int len, nlines;
  int done, our_interface;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  nlines = 0;
  other_chain = '\0';

  /* Form the name of the chains.dat file in the PDBsum directory */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"chains.dat");

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL &&
             done == FALSE)
        {
          /* Increment count of lines read in */
          nlines++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If line is blank, then take to be chain blank */
          if (len == 0)
            chain_id = ' ';
          else
            chain_id = input_line[0];

          /* If this is a copy of another chain, store that */
          if (len > 2)
            {
              /* Get the chain this is a copy of */
              copy_of = input_line[2];

              /* See if the two chains correspond to the two chains
                 in the protein-protein interface we are plotting */
              our_interface = FALSE;
              if ((chain_id == chains[0] && copy_of == chains[1]) ||
                  (chain_id == chains[1] && copy_of == chains[0]))
                {
                  /* Save the other chain */
                  if (chain == chain_id)
                    other_chain = copy_of;
                  else
                    other_chain = chain_id;

                  /* Set flag that we're done */
                  done = TRUE;
                }
            }

          /* Otherwise, chain is unique */
          else
            {
              copy_of = '\0';
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Show warning tnat chain.dat file not found */
  else
    printf("*** Warning. File not found: %s\n",file_name);

  /* Return the other chain */
  return(other_chain);
}
/***********************************************************************

get_variants  -  Read in the DisaStr veriants on the current LIGPLOT
                 diagram

***********************************************************************/

int get_variants(struct highres **fst_highres_ptr,
                 struct highres **lst_highres_ptr,int nhigh,
                 char other_chain,struct parameters params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chn, chain, res_name[4], res_num[6];
  char number_string[20], tmp[LINELEN + 1];

  char *string_ptr;

  int len, line, loop, nloops, type, var_type;
  int ref, wanted;

  struct highres *highres_ptr, *first_highres_ptr, *last_highres_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  nloops = 1;
  if (other_chain != '\0')
    nloops = 2;
  ref = FALSE;
  var_type = NO_VARIANT;
  type = NOT_WANTED;
  wanted = FALSE;

  /* Initialise pointers */
  first_highres_ptr = *fst_highres_ptr;
  last_highres_ptr = *lst_highres_ptr;
  highres_ptr = NULL;

  /* Form name of input annotations.dat file */
  strcpy(file_name,params.disastr_dir);
  strcat(file_name,"/wiring/annotations.dat");

  /* If file not found, abort */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. File not found: %s\n",file_name);
      return(nhigh);
    }

  /* Loop while reading in records from input file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the input line */
      len = string_truncate(input_line,LINELEN);

      /* Get the position of the space */
      string_ptr = strchr(input_line,' ');

      /* If this is the start of a new residue, then re-initialise */
      if (!strncmp(input_line,"KEY: ",5))
        {
          highres_ptr = NULL;
          var_type = NO_VARIANT;
          type = NOT_WANTED;
          wanted = FALSE;
        }

      /* If variant type, save the type */
      else if (!strncmp(input_line,"INT_VAR_TYPE ",13))
        {
          /* Get the variant type */
          if (string_ptr != NULL)
            {
              /* Get the type */
              strcpy(number_string,string_ptr + 1);
              var_type = atoi(number_string);
              if (var_type < NO_VARIANT)
                var_type = NO_VARIANT;
            }
        }

      /* If interaction type, save the type */
      else if (!strncmp(input_line,"INT_TYPE[",9))
        {
          /* Get the interaction type */
          if (string_ptr != NULL)
            {
              /* Get the type */
              strcpy(number_string,string_ptr + 1);
              type = atoi(number_string);
            }
        }

      /* If this is the start of a new residue, then re-initialise */
      else if (!strncmp(input_line,"INT_PDB_CODE[",13))
        {
          /* Reinitialise flag */
          wanted = FALSE;

          /* Check if this is the right PDB code */
          if (string_ptr != NULL && !strcmp(string_ptr + 1,
                                            params.pdb_code))
            wanted = TRUE;

          /* If the wrong type of interaction, then we don't want it */
          if (type != TO_PROTEIN || var_type == NO_VARIANT)
            wanted = FALSE;
        }

      /* Only consider line if we have the right type of interaction */
      if (wanted == TRUE && string_ptr != NULL)
        {
          /* If this is the residue name, then save */
          if (!strncmp(input_line,"INT_RES_NAME[",13))
            strcpy(res_name,string_ptr + 1);

          /* If this is the residue number, then save */
          else if (!strncmp(input_line,"INT_RES_NUM[",12))
            {
              /* Get the raw number */
              strcpy(tmp,string_ptr + 1);

              /* Format the number */
              format_res_number(tmp,res_num);
            }
              
          /* If this is the chain, then save */
          else if (!strncmp(input_line,"INT_CHAIN[",10))
            chain = string_ptr[1];

          /* If this is the DIMPLOT identifer, save the residue */
          else if (!strncmp(input_line,"INT_DIMPLOT[",12))
            {
              /* Initialise chain */
              chn = chain;

              /* Loop twice if we have a homodimer */
              for (loop = 0; loop < nloops; loop++)
                {
                  /* See if we already have this residue */
                  highres_ptr
                    = find_highres_entry(first_highres_ptr,res_name,
                                         res_num,chn);

                  /* If not the reference residue, then create a record
                     for it */
                  if (highres_ptr == NULL)
                    {
                      /* Create record */
                      highres_ptr
                        = create_highres_entry(&first_highres_ptr,
                                               &last_highres_ptr,var_type,
                                               ref,res_name,res_num,chn);

                      /* Increment count of variants stored */
                      nhigh++;

                      /* Set chain to the other member of the homodimer */
                      chn = other_chain;
                    }

                  /* If this is the reference residue, then save the type */
                  else if (highres_ptr->ref == TRUE)
                    highres_ptr->type = type;
                }
            }
        }
    }

  /* Close the input file */
  fclose(file_ptr);

  /* Return the pointers */
  *fst_highres_ptr = first_highres_ptr;
  *lst_highres_ptr = last_highres_ptr;

  /* Return the number of highlighted residues */
  return(nhigh);
}
/***********************************************************************

store_variants  -  If called from DisaStr, pick up all the disease-
                   associated variants for the given UniProt entry

***********************************************************************/

struct highres *store_variants(struct highres **fst_highres_ptr,
                               char *pdbsum_dir,struct parameters params)
{
  char chain, other_chain, res_name[4], res_num[6];

  int dir_type, nhigh, pdb_type, var_type;
  int ref;

  struct highres *highres_ptr, *first_highres_ptr, *last_highres_ptr;
  struct highres *our_highres_ptr;
    
  /* Initialise variables */
  nhigh = 0;
  pdb_type = STANDARD_PDB;
  ref = TRUE;
  var_type = NO_VARIANT;

  /* Initialise pointers */
  our_highres_ptr = highres_ptr = first_highres_ptr
    = last_highres_ptr = NULL;

  /* If we have a residue to highlight, then process variants */
  if (params.highlight_res[0] != '\0')
    {
      /* Unpack the name of the highlighted residue */
      extract_residue_name(params.highlight_res,res_name,res_num,
                           &chain);

      /* Pick up the chains from the chains.dat file, if present */
      other_chain = get_other_chain(chain,params.chain,pdbsum_dir);

      /* Create a highlight record for this residue */
      highres_ptr
        = create_highres_entry(&first_highres_ptr,&last_highres_ptr,
                               var_type,ref,res_name,res_num,chain);
      nhigh++;

      /* Save the pointer to our highres record */
      our_highres_ptr = highres_ptr;

      /* Unset flag */
      ref = FALSE;

      /* If we have a homodimer interface, create a highlight record
         for the same residue on the other chain */
      if (other_chain != '\0')
        {
          /* Create a highlight record for this residue */
          highres_ptr
            = create_highres_entry(&first_highres_ptr,&last_highres_ptr,
                                   var_type,ref,res_name,res_num,
                                   other_chain);
          nhigh++;
        }

      /* Get all the other residues to be highlighted */
      nhigh = get_variants(&first_highres_ptr,&last_highres_ptr,nhigh,
                           other_chain,params);
    }

  /* Return current pointers */
  *fst_highres_ptr = first_highres_ptr;

  /* Return pointer to our residue */
  return(our_highres_ptr);
}  
/* <--v.5.3.8 */
/***********************************************************************

convert_residue_name  -  Convert residue name and number as it is to
                         appear on the output

***********************************************************************/

int convert_residue_name(char *res_name,char *res_num,char chain,
                         char *res_name_string,int show_chains)
{
  char converted_res_name[4], chain_string[5];
  char upper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char lower[] = "abcdefghijklmnopqrstuvwxyz";
  char digit[] = "0123456789";

  int i, j, len, found, number;

  /* Check whether any digit is a number */
  number = FALSE;
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 10 && number == FALSE; j++)
        {
          if (res_name[i] == digit[j])
            number = TRUE;
        }
    }

  /* If residue name is a water, or last character is a number,
     or first character is a space, then keep in upper case */
  if (!strncmp(res_name,"HOH",3) || !strncmp(res_name,"WAT",3) ||
      res_name[0] == ' ' || number == TRUE)
    {
      strncpy(converted_res_name,res_name,3);
      converted_res_name[3] = '\0';
    }
  
  /* Otherwise, convert 2nd and 3rd letters of residue-name to lower-case */
  else
    {
      converted_res_name[0] = res_name[0];
      for (i = 1; i < 3; i++)
        {
          converted_res_name[i] = res_name[i];
          found = FALSE;
          for (j = 0; j < 26 && found == FALSE; j++)
            {
              if (res_name[i] == upper[j])
                {
                  found = TRUE;
                  converted_res_name[i] = lower[j];
                }
            }
        }
      converted_res_name[3] = '\0';
    }

  /* Transfer new residue name into string holding final residue
     details */
  strcpy(res_name_string,converted_res_name);

  /* Add the residue number */
  if (res_num[2] == ' ')
    strcat(res_name_string,res_num+3);
  else if (res_num[1] == ' ')
    strcat(res_name_string,res_num+2);
  else if (res_num[0] == ' ')
    strcat(res_name_string,res_num+1);
  else
    strcat(res_name_string,res_num);

  /* Check for blank insertion code */
  len = strlen(res_name_string);
  if (res_name_string[len - 1] == ' ')
    res_name_string[len - 1] = '\0';

  /* Add the chain-id, if there is one */
  if (chain != ' ' && show_chains == TRUE)
    {
      sprintf(chain_string,"(%c)",chain);
      strcat(res_name_string,chain_string);
    }

  /* Return the string length */
  len = strlen(res_name_string);
  return(len);
}
/* v.5.3.8--> */
/***********************************************************************

trim_plot  -  Trim the more distant residues from the variant residue

***********************************************************************/

int trim_plot(struct residue **alignment,int align_len,int total_len,
              struct highres *our_highres_ptr)
{
  int end_pos, hpos, ipos, ires, jpos, start_pos;

  struct residue *residue_ptr[2];

  /* Initialise variables */
  hpos = -1;
  jpos = 0;

  /* Loop through all the positions in the alignment to identify the
     highlighted residue */
  for (ipos = 0; ipos < total_len && hpos == -1; ipos++)
    {
      /* Get the corresponding residue pointers at this point in
         the alignment */
      residue_ptr[0] = alignment[ipos];
      residue_ptr[1] = alignment[ipos + align_len];

      /* If either residue is our residue, save the position at which
         it occurs */
      if ((residue_ptr[0] != NULL &&
           residue_ptr[0]->chain == our_highres_ptr->chain &&
           !strcmp(residue_ptr[0]->res_num,our_highres_ptr->res_num) &&
           !strcmp(residue_ptr[0]->res_name,our_highres_ptr->res_name)) ||
          (residue_ptr[1] != NULL &&
           residue_ptr[1]->chain == our_highres_ptr->chain &&
           !strcmp(residue_ptr[1]->res_num,our_highres_ptr->res_num) &&
           !strcmp(residue_ptr[1]->res_name,our_highres_ptr->res_name)))
        hpos = ipos;
    }

  /* If not found, set the position at halfway along */
  if (hpos == -1)
    hpos = total_len / 2;

  /* Determine the trim points */
  start_pos = hpos - MAX_ALIGN_LEN / 2;
  if (start_pos < 0)
    start_pos = 0;
  end_pos = start_pos + MAX_ALIGN_LEN - 1;
  if (end_pos > total_len)
    end_pos = total_len;

  /* Loop through all the alignment positions again to delete the
     terminal residues */
  for (ipos = 0; ipos < total_len; ipos++)
    {
      /* Get the corresponding residue pointers at this point in
         the alignment */
      residue_ptr[0] = alignment[ipos];
      residue_ptr[1] = alignment[ipos + align_len];

      /* If position is inside the wanted range, place in next available
         slot */
      if (ipos >= start_pos && ipos <= end_pos)
        {
          /* Save in next available slot */
          alignment[jpos] = alignment[ipos];
          alignment[jpos + align_len] = alignment[ipos + align_len];

          /* Increment location of next free spot */
          jpos++;
        }

      /* Otherwise, mark the residues as deleted */
      else
        {
          /* Delete residues */
          if (residue_ptr[0] != NULL)
            residue_ptr[0]->deleted = TRUE;
          if (residue_ptr[1] != NULL)
            residue_ptr[1]->deleted = TRUE;
        }
    }

  /* Save the new total length */
  total_len = jpos;

  /* Return the new total length */
  return(total_len);
}
/***********************************************************************

plot_ellipse  -  Plot the ellipse around the highlighted residue

***********************************************************************/

void plot_ellipse(FILE *ps_file,float centre_x,float centre_y,
                  float scale_factor,float ellipse_radius,
                  struct highres *highres_ptr)
{
  char colour[60], comment[60];

  int type;

  float angle, line_width, ratio, width;

  /* Initialise variables */
  angle = 0;
  type = highres_ptr->type;

  /* Print comment */
  if (highres_ptr->ref == TRUE)
    {
      /* Define the type as the added residue */
      type = ADDED_VAR;

      /* Show heading */
      sprintf(comment,"Highlighted residue: %s%s%c  Type: %s",
              highres_ptr->res_name,highres_ptr->res_num,
              highres_ptr->chain,ELLIPSE_COL_DEF[type + 2]);
    }
  else
    sprintf(comment,"Variant residue: %s%s%c  Type: %s",
            highres_ptr->res_name,highres_ptr->res_num,
            highres_ptr->chain,ELLIPSE_COL_DEF[type + 2]);
  pscomm_(ps_file,comment);
  pssave_(ps_file);

  /* Define ellipse's line width */
  line_width = scale_factor * ELLIPSE_LINEWIDTH;
/* v.5.4.4 --> 
  if (line_width < 3.0)
    line_width = (float) 3.0;
 <-- v.5.4.4 */
  pslwid_(ps_file,line_width);

  /* Set the line colour */
  sprintf(colour,"Ellipse_col_%s",ELLIPSE_COL_DEF[type + 2]);
  pscolb_(ps_file,colour);

  /* If not the reference residue, make the ellipse dashed */
  if (highres_ptr->ref == FALSE)
    psdots_(ps_file,1,8);

  /* Calculate the height to width ratio */
  ratio = 1.0 / ELLIPSE_ASPECT;

  /* Convert the ellipse width to PostScript equivalent */
  width = 2 * ELLIPSE_FACTOR * ellipse_radius;

  /* Draw the ellipse */
  psoeli_(ps_file,centre_x,centre_y,width,ratio,angle);
  
  /* If not the reference residue, unset dashed lines */
  if (highres_ptr->ref == FALSE)
    psdash_(ps_file,0);

  /* Reset graphics state */
  psrest_(ps_file);
  pscomm_(ps_file,"---");
}
/***********************************************************************

plot_variant  -  If this is a variant residue, plot an ellipse around it

***********************************************************************/

void plot_variant(FILE *ps_file,struct highres *first_highres_ptr,
                  char *res_name,char *res_num,char chain,
                  float centre_x,float centre_y,float ellipse_radius,
                  float scale_factor)
{
  int var_type;

  struct highres *highres_ptr;

  /* See if this residue is to be highlighted */
  highres_ptr
    = find_highres_entry(first_highres_ptr,res_name,res_num,chain);

  /* If it is a variant residue, work out parameters of ellipse
     to be plotted */
  if (highres_ptr != NULL &&
      (highres_ptr->type != NO_VARIANT || highres_ptr->ref == TRUE))
    {
      /* Get its variant type */
      var_type = highres_ptr->type;

      /* Plot the ellipse around the highlighted residue */
      plot_ellipse(ps_file,centre_x,centre_y,scale_factor,ellipse_radius,
                   highres_ptr);
    }
}
/* <--v.5.3.8 */
/***********************************************************************

plot_ellipses  -  Plot the residues as ellipses

***********************************************************************/

float plot_ellipses(struct residue **alignment,int align_len,int total_len,
                    float scale_factor,float origin[2],int labels_outside,
/* v.5.3.8--> */
//                    int show_chains,struct pspage *pspage_ptr)
                    int show_chains,struct highres *first_highres_ptr,
                    struct pspage *pspage_ptr)
/* <--v.5.3.8 */
{
  char colour[COLNAMLEN + 1], last_colour[COLNAMLEN + 1];
  char residue_name[20];

  int ipos, ires, len, nextra;
  int labels_inside;

  float label_gap, line_width, radius, text_size, x, xlabel, y;
  float yextra, ygap;

  struct residue *residue_ptr[2];

  /* Initialise variables */
  last_colour[0] = '\0';
  label_gap = scale_factor * LABEL_GAP;
  labels_inside = TRUE;
  nextra = 0;

  /* Define the y separation between the ellipses */
  ygap = scale_factor * (2 * ELLIPSE_SIZE + INTERELLIPSE_GAP);
  yextra = ygap * EXTRA_FACTOR;

  /* Define colours and line widths */
  line_width = scale_factor * RES_BORDER_WIDTH;
/* v.5.4.4 --> 
  if (line_width < 1.0)
    line_width = (float) 1.0;
 <-- v.5.4.4 */
  pscomm_(pspage_ptr->ps_file,"Residues");
  pscolb_(pspage_ptr->ps_file,"black");
  pslwid_(pspage_ptr->ps_file,line_width);

  /* Calculate text size, adjusting if too small */
  text_size = scale_factor * TEXT_HEIGHT;
  if (text_size < MIN_SCALED_SIZE || labels_outside == TRUE)
    {
      /* Fix text size to minimum and set flag to print outside residue
         ellipses */
/* v.5.4.4 --> */
//      text_size = MIN_SCALED_SIZE;
      text_size = 2.0 * text_size;
/* <-- v.5.4.4 */
      labels_inside = FALSE;
    }

  /* Define starting location */
  x = origin[0];
  y = origin[1];
  radius = scale_factor * ELLIPSE_SIZE;

  /* Loop through all the positions in the alignment to fill any blanks */
  for (ipos = 0; ipos < total_len; ipos++)
    {
      /* Get the corresponding residue pointers at this point in
         the alignment */
      residue_ptr[0] = alignment[ipos];
      residue_ptr[1] = alignment[ipos + align_len];

      /* If either residue is involved in H-bonds, need to allow an
         extra spacing for DIMPLOT */
      if ((residue_ptr[0] != NULL && residue_ptr[0]->nhbonds > 0) ||
          (residue_ptr[1] != NULL && residue_ptr[1]->nhbonds > 0))
        nextra++;

      /* Loop over the two residues at this position */
      for (ires = 0; ires < 2; ires++)
        {
          /* If have residue at this position, plot it */
          if (residue_ptr[ires] != NULL)
            {
              /* Write out residue name to PostScript file */

              /* Get this residue's colour */
              strcpy(colour,residue_ptr[ires]->colour);

              /* Define ellipse colour */
              if (strcmp(colour,last_colour))
                psccol_(pspage_ptr->ps_file,colour);

              /* Save current colour */
              strcpy(last_colour,colour);

              /* Plot the ellipse */
              pselip_(pspage_ptr->ps_file,x,y,radius,ELLIPSE_ASPECT,1.0);

              /* Save this residue's y-coordinate for DIMPLOT restraints */
              residue_ptr[ires]->ycoord = y - nextra * yextra;

              /* Format residue name */
              len
                = convert_residue_name(residue_ptr[ires]->res_name,
                                       residue_ptr[ires]->res_num,
                                       residue_ptr[ires]->chain,
                                       residue_name,show_chains);

              /* Print residue name */
              if (labels_inside == TRUE)
                psctxt_(pspage_ptr->ps_file,x,y,text_size,residue_name);

              /* Print residue name to left or right of residue */
              else
                {
                  /* Print label to left */
                  if (ires == 0)
                    {
                      /* Calculate label's x-coord */
                      xlabel = x - radius * ELLIPSE_ASPECT - label_gap;

                      /* Print residue name */
                      psrtxt_(pspage_ptr->ps_file,xlabel,y,text_size,
                              residue_name);
                    }

                  /* Print label to right */
                  else
                    {
                      /* Calculate label's x-coord */
                      xlabel = x + radius * ELLIPSE_ASPECT + label_gap;

                      /* Print residue name */
                      psltxt_(pspage_ptr->ps_file,xlabel,y,text_size,
                              residue_name);
                    }
                }

/* v.5.3.8--> */
              /* If we have variants to be highlighted, see if this
                 residue is one of them */
              if (first_highres_ptr != NULL)
                plot_variant(pspage_ptr->ps_file,first_highres_ptr,
                             residue_ptr[ires]->res_name,
                             residue_ptr[ires]->res_num,
                             residue_ptr[ires]->chain,x,y,radius,
                             scale_factor);
/* <--v.5.3.8 */
            }

          /* Increment x position */
          x = x + scale_factor
            * (2 * ELLIPSE_SIZE * ELLIPSE_ASPECT + INTERFACE_GAP);
        }

      /* Decrement y-value and reset x-value */
      x = origin[0];
      y = y - ygap;
    }

  /* Return the y separation */
  return(ygap);
}
/***********************************************************************

close_page  -  Close current PostScript file if have an open one

***********************************************************************/

void close_page(struct psparams psparams_ptr,struct pspage *pspage_ptr)
{
  /* Close off previous file, if there is one */
  if (pspage_ptr->page > 0)
    {
      psendp_(pspage_ptr->ps_file,pspage_ptr->lim,
            psparams_ptr.program_name);
      psclos_(pspage_ptr->ps_file);
    }
}
/***********************************************************************

plot_diagram  -  Plot the residue interactions across the interface

***********************************************************************/

float plot_diagram(struct residue **alignment,int align_len,
                   int total_len,int *max_count,
                   struct psparams psparams_ptr,
                   struct pspage *pspage_ptr,char *pdbsum_dir,
                   struct parameters params)
{
  char name[20], colour[COLNAMLEN + 1];

  float frame_aspect, height, picture_aspect, scale_factor, width;
  float line_sep, origin[2], stripe_width, text_size, x, y, ygap;

/* v.5.3.8--> */
  struct highres *highres_ptr, *our_highres_ptr, *first_highres_ptr;
/* <--v.5.3.8 */

  /* Initialise variables */
  our_highres_ptr = first_highres_ptr = NULL;

/* v.5.3.8--> */
  /* If we are plotting variant residues, pick up all the variants to
     be highlighted */
  if (params.highlight_res[0] != '\0')
    {
      /* If we have a DisaStr residue to be highlighted, get all the variants
         for this UniProt entry */
      our_highres_ptr
        = store_variants(&first_highres_ptr,pdbsum_dir,params);

      /* If alignment too long to show clearly, trim the more distant
         residues */
      if (total_len > MAX_ALIGN_LEN)
        {
          /* Trim the plot by removing the outs sets of residues */
          total_len
            = trim_plot(alignment,align_len,total_len,our_highres_ptr);
        }
    }
/* <--v.5.3.8 */

  /* Form name of output PostScript file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(pspage_ptr->ps_name,params.work_dir);
      strcat(pspage_ptr->ps_name,DIR_SEP);
    }
  strcat(pspage_ptr->ps_name,"dimhtml.ps");

  /* Determine extent of picture */
  width = 4 * ELLIPSE_ASPECT * ELLIPSE_SIZE + INTERFACE_GAP;
  height = 2 * total_len * ELLIPSE_SIZE + (total_len - 1) * INTERELLIPSE_GAP
    + (float) 1.5 * HEADTEXT_HEIGHT;

  /* Calculate aspect ratio for picture and for frame within picture to
     be placed */
  picture_aspect = height / width;
  frame_aspect = PICTURE_HEIGHT / PICTURE_WIDTH;

  /* Determine optimal scaling factor for the picture */
  if (picture_aspect > frame_aspect)
    {
      /* Bounded by y: calculate the scaling factor */
      scale_factor = PICTURE_HEIGHT / height;

      /* Determine origin of top left ellipse */
      origin[0] = (PICTURE_WIDTH - scale_factor * width) / 2
        + scale_factor * ELLIPSE_SIZE * ELLIPSE_ASPECT;
      origin[1] = PICTURE_HEIGHT
        - scale_factor * ((float) 1.5 * HEADTEXT_HEIGHT + ELLIPSE_SIZE);
    }
  else
    {
      /* Bounded by x: calculate the scaling factor */
      scale_factor = PICTURE_WIDTH / width;

      /* Determine origin of bottom left ellipse */
      origin[0] = scale_factor * ELLIPSE_ASPECT * ELLIPSE_SIZE;
      origin[1] = PICTURE_HEIGHT
        - (PICTURE_HEIGHT - scale_factor * height) / 2
        - scale_factor * ((float) 1.5 * HEADTEXT_HEIGHT + ELLIPSE_SIZE);
    }

  /* Determine the picture origin */
  origin[0] = origin[0] + XMARGIN;
  origin[1] = origin[1] + YMARGIN;

  /* Determine the separation between H-bond lines between any
     residue pair */
  line_sep = DEFAULT_LINE_SEP;
  if (line_sep * (max_count[HBONDS] - 1) > ELLIPSE_SIZE)
    line_sep = ELLIPSE_SIZE / (max_count[HBONDS] - 1);
  line_sep = line_sep * scale_factor;

  /* Determine the width of the stripes representing non-bonded
     contacts  */
  stripe_width = STRIPE_WIDTH;
  if (stripe_width * (max_count[CONTACTS]) > ELLIPSE_SIZE)
    stripe_width = ELLIPSE_SIZE / (max_count[CONTACTS]);
  stripe_width = stripe_width * scale_factor;

  /* Initialise PostScript page variables */
  init_pspage(psparams_ptr,pspage_ptr);

  /* Start the page holding the plot */
  start_new_page(psparams_ptr,pspage_ptr);

  /* Write out chain headings */
  text_size = scale_factor * HEADTEXT_HEIGHT;
  if (text_size < MIN_HEAD_HEIGHT)
    text_size = MIN_HEAD_HEIGHT;
  x = origin[0] + scale_factor * ELLIPSE_ASPECT * ELLIPSE_SIZE;
  y = origin[1] + scale_factor * ELLIPSE_SIZE + text_size;

  /* Get colour of first chain */
  if (params.chains == TRUE)
    get_chain_col(params.chain[0],colour);
  else
    get_domain_colour(params.domain[0],colour);
  pscolb_(pspage_ptr->ps_file,colour);

  /* Write out heading */
  if (params.chains == TRUE)
    sprintf(name,"Chain %c",params.chain[0]);
  else
    sprintf(name,"Domain %d",params.domain[0]);
  psrtxt_(pspage_ptr->ps_file,x,y,text_size,name);
  x = x + scale_factor * INTERFACE_GAP;

  /* Get colour of second chain */
  if (params.chains == TRUE)
    get_chain_col(params.chain[1],colour);
  else
    get_domain_colour(params.domain[1],colour);
  pscolb_(pspage_ptr->ps_file,colour);

  /* Write out heading */
  if (params.chains == TRUE)
    sprintf(name,"Chain %c",params.chain[1]);
  else
    sprintf(name,"Domain %d",params.domain[1]);
  psltxt_(pspage_ptr->ps_file,x,y,text_size,name);
  pscolb_(pspage_ptr->ps_file,"black");

  /* Draw lines representing the interactions across the interface */
  plot_interactions(alignment,align_len,total_len,scale_factor,origin,
                    line_sep,stripe_width,pspage_ptr);

/* <--v.5.3.8 */

  /* Draw ellipses */
  ygap
    = plot_ellipses(alignment,align_len,total_len,scale_factor,origin,
/* v.5.3.8--> */
//                    params.labels_outside,params.show_chains,pspage_ptr);
                    params.labels_outside,params.show_chains,
                    first_highres_ptr,pspage_ptr);
/* <--v.5.3.8 */

  /* Close the current page */
  close_page(psparams_ptr,pspage_ptr);

  /* Return the y-gap separation of ellipses on the plot */
  return(ygap);
}
