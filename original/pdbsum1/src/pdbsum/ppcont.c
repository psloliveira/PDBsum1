/***********************************************************************

  ppcont.c - Program to plot schematic diagram of protein-protein
             interface as two touching spheres

************************************************************************

Date:-          6 Dec 2004.
Dir update:-   17 Oct 2012

Written by:-   Roman Laskowski


----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name       Description
-----------       -----------
interfaces.dat    Input file containing protein-protein interface data
ppcont.ps         Output PostScript file

*/

/* Compilation
   -----------

cc -c ppcont.c
cc -c common.c
cc -c reapar.c
cc -o ppcont ppcont.o common.o reapar.o -lm -lz

Optimization flags:-  

cc -O2 -funroll-loops -c ppcont.c
cc -O2 -funroll-loops -c common.c
cc -O2 -funroll-loops -c reapar.c
cc -o ppcont ppcont.o common.o reapar.o -lm -lz

*/

/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> create_field_record
               -> store_string (in common.c)
   -> find_pdbsum_dir (in common.c)
   -> read_database_file
         -> string_truncate (in common.c)
   -> plot_diagram
         -> init_pspage
               -> adjust_for_landscape
         -> start_new_page
               -> psendp_
               -> psclos_
               -> pspage_
                     -> psrot_
               -> psopen_
                     ->  get_ppcolour_rgb
         -> get_ppchain_colour (in common.c)
         -> pscolb_
         -> psrtxt_
         -> psltxt_
         -> plot_circles
               -> pscomm_
               -> pscolb_
               -> pslwid_
               -> get_ppchain_colour (in common.c)
               -> psccol_
               -> pscirc_
               -> pspcir_
         -> close_page
   -> read_chains
         -> find_chain
         -> create_chain_record
         -> create_interface_record
   -> count_clusters
         -> create_cluster_record
   -> arrange_clusters
         -> compute_interaction_wedges
               -> create_wedge_record
         -> calc_global_radius
               -> calc_circle_coords
                     -> calc_separation
         -> place_chains
               -> get_random_number
               -> sort_chains
               -> calc_energy
                     -> chain_clash_energy
                     -> overlap_energy
                     -> calc_interaction_angle
                           -> get_angle
                     -> wedge_clash_energy
         -> energy_minimize
               -> get_random_number
               -> calc_energy (as above)
               -> calc_force_vector
                     -> calc_energy (as above)
               -> shift_coords
   -> plot_allchains_diagram
         -> write_coords
               -> compact_number
         -> init_pspage
               -> adjust_for_landscape
         -> start_new_page
               -> psendp_
               -> psclos_
               -> pspage_
                     -> psrot_
               -> psopen_
         -> plot_all_interactions
               -> pdlwid_
               -> get_angle
               -> pscolb_
               -> psline_
         -> plot_all_circles
               -> pslwid_
               -> pscolb_
               -> get_ppchain_colour (in common.c)
               -> psccol_
               -> pscirc_
         -> plot_wedges
               -> pslwid_
               -> get_ppchain_colour (in common.c)
               -> pscolb_
               -> psccol_
               -> pspcir_
               -> pscirc_
               -> psctxt_
         -> close_page

*/

#include "common.h"
#include "ppcont.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
float get_random_number(int *iseed,int random_start);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int ntoken,
                           struct parameters *params)
{
  int itoken;

  /* Initialise parameters */
  params->pdb_code[0] = '\0';
  params->chains[0] = '\0';
  params->work_dir[0] = '\0';
  params->pdb_type = STANDARD_PDB;
  params->all_chains = FALSE;
  params->save = FALSE;
  params->read = FALSE;
  params->run_64 = FALSE;

  /* If nothing entered on the command line, show options available */
  if (ntoken < 2 || (ntoken == 1 && (!strcmp(strptrs[1],"-help") ||
                                     !strcmp(strptrs[1],"-h"))))
    {
      printf("\n");
      printf("Correct usage is ...\n"
             "\n"
             "  ppcont  pdb_codeXY  [-dir work_dir]  [UPLOAD | MCSG]  [-all]"
             "  [-save]  [-read]  [-64]\n"
             "          [-pdbsum1]\n"
             "\n"
             "where\n"
             "   * pdb_code      = PDB entry to be displayed\n"
             "   * XY            = Chains defining interface\n"
             "   * UPLOAD        = Uploaded PDB file for PDBsum generation\n"
             "   * MCSG          = MCSG structure\n"
             "   * -dir work_dir = Directory containing input and "
             "output files\n"
             "   * -all          = Plot all chains and their partners\n"
             "   * -read         = Read in the all-chains coords from "
             "ppcont.dat file\n"
             "   * -save         = Save the all-chains coords to the "
             "ppcont.dat file\n"
             "   * -pdbsum1      = PDBsum1 structure\n"
             "   * -64           = run on 64-bit processor\n"
             "\n"
             "\n");
      exit(1);
    }

  /* Take first parameter to be the PDB code and chain identifiers */
  strncpy(params->pdb_code,strptrs[1],4);
  params->pdb_code[4] = '\0';
  if (strlen(strptrs[1]) > 5)
    {
      params->chains[0] = strptrs[1][4];
      params->chains[1] = strptrs[1][5];
    }

  /* Loop through any additional command arguments */
  for (itoken = 2; itoken < ntoken + 1; itoken++)
    {
      /* Check for the -dir option giving the name of the working
         directory */
      if (!strcmp(strptrs[itoken],"-dir"))
        {
          /* Get the directory name from the next token, if there is one */
          if (itoken < ntoken)
            {
              strcpy(params->work_dir,strptrs[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the UPLOAD option */
      else if (!strcmp(strptrs[itoken],"UPLOAD"))
        params->pdb_type = UPLOAD_PDB;

      /* Check for the MCSG option */
      else if (!strcmp(strptrs[itoken],"MCSG"))
        params->pdb_type = MCSG_PDB;

      /* Check for the all option to plot all chains */
      else if (!strcmp(strptrs[itoken],"-all"))
        params->all_chains = TRUE;

      /* Check for option to save the all-chains coords */
      else if (!strcmp(strptrs[itoken],"-save"))
        params->save = TRUE;

      /* Check for the read option to pick up the all-chains coords */
      else if (!strcmp(strptrs[itoken],"-read"))
        params->read = TRUE;

      /* Check for the PDBsum1 option */
      else if (!strcmp(strptrs[itoken],"-pdbsum1"))
        params->pdb_type = PDBSUM1;

      /* Check for the 64-bit option */
      else if (!strcmp(strptrs[itoken],"-64"))
        params->run_64 = TRUE;
    }
}
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
  strcpy(psparams_ptr->program_name,"ppcont");

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

read_database_file  -  Read in the required data for the plot from the
                       interfaces.dat file

***********************************************************************/

void read_database_file(char *pdbsum_dir,char *fname,
                        struct interface *interface_ptr,
                        struct parameters params)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char chains[2], number_string[20];

  char *string_ptr;

  int i, ichain, len;
  int done, found;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  found = FALSE;
  interface_ptr->nresid[0] = 0;
  interface_ptr->nresid[1] = 0;
  interface_ptr->chain_area[0] = 0.0;
  interface_ptr->chain_area[1] = 0.0;
  interface_ptr->interface_area[0] = 0.0;
  interface_ptr->interface_area[1] = 0.0;
  interface_ptr->angle = 0.0;
  interface_ptr->ntypes = 0;
  for (i = 0; i < NTYPES; i++)
    interface_ptr->count[i] = 0;

  /* Form name of the interfaces.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,fname);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open interfaces.dat file: %s\n",file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* If this is a chain identifier, store */
      if (!strncmp(input_line,"PPIFACE_CHAIN",13))
        {
          /* Get which chain this is */
          number_string[0] = input_line[13];
          number_string[1] = '\0';
          ichain = atoi(number_string) - 1;

          /* Get the chain */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            chains[ichain] = string_ptr[1];
          else
            chains[ichain] = '\0';

          /* If this is the second chain, check whether this is the
             interface we're after */
          if (ichain == 1)
            {
              /* If we've already found our interface, we're done */
              if (found == TRUE)
                done = TRUE;

              /* Check for match */
              else if ((params.chains[0] == chains[0] &&
                   params.chains[1] == chains[1]) ||
                  (params.chains[0] == chains[1] &&
                   params.chains[1] == chains[0]))
                {
                  /* Set flag */
                  found = TRUE;
                }

              /* Otherwise, not a match */
              else
                found = FALSE;
            }
        }

      /* If this is a the number of residues, extract */
      else if (found == TRUE &&
               !strncmp(input_line,"PPIFACE_NRESID",14))
        {
          /* Get which chain this is */
          string_ptr = strstr(input_line,"[");
          if (string_ptr != NULL)
            number_string[0] = string_ptr[-1];
          else
            number_string[0] = '0';
          number_string[1] = '\0';
          ichain = atoi(number_string) - 1;

          /* Get the number of residues */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            strcpy(number_string,string_ptr+1);
          else
            number_string[0] = '\0';

          /* Store appropriate number */
          interface_ptr->nresid[ichain] = atoi(number_string);
        }

      /* If this is a chain area, extract */
      else if (found == TRUE &&
               (!strncmp(input_line,"PPIFACE_AREA",12) ||
                !strncmp(input_line,"PPIFACE_BURIED",14)))
        {
          /* Get which chain this is */
          string_ptr = strstr(input_line,"[");
          if (string_ptr != NULL)
            number_string[0] = string_ptr[-1];
          else
            number_string[0] = '0';
          number_string[1] = '\0';
          ichain = atoi(number_string) - 1;

          /* Get the area */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            strcpy(number_string,string_ptr+1);
          else
            number_string[0] = '\0';

          /* Store appropriate area */
          if (!strncmp(input_line,"PPIFACE_AREA",12))
            interface_ptr->chain_area[ichain] = atof(number_string);
          else
            interface_ptr->interface_area[ichain] = atof(number_string);
        }

      /* If this is a count of interactions, extract */
      else if (found == TRUE &&
               (!strncmp(input_line,"PPIFACE_NHBONDS",15) ||
                !strncmp(input_line,"PPIFACE_NCONTACTS",17) ||
                !strncmp(input_line,"PPIFACE_SALT_BRIDGES",20) ||
                !strncmp(input_line,"PPIFACE_DISULPHIDES",19)))
        {
          /* Get the count */
          string_ptr = strstr(input_line," ");
          if (string_ptr != NULL)
            strcpy(number_string,string_ptr+1);
          else
            number_string[0] = '\0';

          /* Store appropriate count */
          if (strstr(input_line,"HBONDS"))
            interface_ptr->count[HBONDS] = atoi(number_string);
          else if (strstr(input_line,"CONTACTS"))
            interface_ptr->count[CONTACTS] = atoi(number_string);
          else if (strstr(input_line,"SALT_BRIDGES"))
            interface_ptr->count[SALT_BRIDGES] = atoi(number_string);
          else if (strstr(input_line,"DISULPHIDES"))
            interface_ptr->count[DISULPHIDES] = atoi(number_string);
        }
    }

  /* Close the file */
  fclose(file_ptr);
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

psbtxt_   Write out bold centred text to PostScript file

***********************************************************************/

void psbtxt_(FILE *ps_file,float x,float y,float size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %7.2f Center\n",text_string,size);
  fprintf(ps_file,"(%s) %7.2f Bprint\n",text_string,size);
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
/***********************************************************************

pselip_  -  Write out a coloured ellipse to PostScript file

***********************************************************************/

void pselip_(FILE *ps_file,float x,float y,float radius,float elong_x,
             float elong_y)
{
    fprintf(ps_file,"0.0 0.0 %7.2f %7.2f %7.2f %7.2f %7.2f Ellipse\n",
            radius,elong_x,elong_y,x,y);
}
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

pspcir_  -  Write out a shaded/coloured part-circle to PostScript file
            Start- and end-angles supplied in degrees, the shading going
            anti-clockwise from angle_start to angle_end, the zero point
            being at 3 o'clock

***********************************************************************/

void pspcir_(FILE *ps_file,float x,float y,float radius,float angle_start,
             float angle_end)
{
    fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pcircle\n",
            x,y,radius,angle_start,angle_end,x,y);
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

get_ppcolour_rgb  -  Get the colour name for the given colour number and
                     the corresponding RGB values

***********************************************************************/

void get_ppcolour_rgb(int colour,char *c_name,float *rgb)
{
#define NPPCOLOURS 16
  
  static char *colour_list[] = {
    "black", "red", "green", "blue", "yellow", "purple", "brown",
    "sky_blue", "orange", "cyan", "grey", "white", "magenta", "pink",
    "dark_grey", "redorange", "greenblue"
  };

  static float rgb_triplet[NPPCOLOURS + 1][3] = {
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

  /* Get the name and RGB triplet for this colour */
  strcpy(c_name,colour_list[colour]);
  rgb[0] = rgb_triplet[colour][0];
  rgb[1] = rgb_triplet[colour][1];
  rgb[2] = rgb_triplet[colour][2];
}
/***********************************************************************

psopen_  -  Open PostScript file and write out header records

***********************************************************************/

FILE *psopen_(char *filename,char *program_name)
{
  char  c_name[COLNAMLEN], colno[3], exclamation, percent;
  float grey_shade, rgb[3];
  int   icolour, igrey;

  FILE *ps_file;

  /* Initialise variables */
  percent = '%';
  exclamation = '!';

  /* Open the PostScript file */
  ps_file = open_output_file(filename,TRUE);

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
  for (icolour = 0; icolour < NPPCOLOURS + 1; icolour++)
    {
      /* Get the colour name and RGB triplet for this colour */
      get_ppcolour_rgb(icolour,c_name,rgb);

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
  fprintf(ps_file,"/Ocircle { gsave newpath 0 360 arc stroke grestore } ");
  fprintf(ps_file," bind def\n");
  fprintf(ps_file,"/PC { moveto arc gsave Circol fill grestore } bind def\n");
  fprintf(ps_file,"/Pcircle { 7 copy PC moveto arc closepath stroke } bind ");
  fprintf(ps_file,"def\n");
  fprintf(ps_file,"/Ucircle { gsave newpath 0 360 arc gsave Circol fill ");
  fprintf(ps_file,"grestore grestore } bind def\n");
  fprintf(ps_file,"/Arc { newpath arc stroke newpath } bind def\n");
  fprintf(ps_file,"/Ellipse { gsave translate scale Circle grestore } ");
  fprintf(ps_file,"bind def\n");
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
  //fprintf(ps_file,"stroke\n");
  fprintf(ps_file,"gsave\n");

  /*Draw in the background colour*/
  //if (in_colour == TRUE)
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

get_ppchain_colour  -  Determine the colour for the given chain

***********************************************************************/

int get_ppchain_colour(char chain,char *colour_name)
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

plot_circles  -  Plot the two chains as circles

***********************************************************************/

void plot_circles(struct interface interface_ptr,char *chains,
                  float scale_factor,float *origin,float *radius,
                  float *angle,float line_sep,int nlines,
                  struct pspage *pspage_ptr)
{
  char colour[COLNAMLEN + 1], number_string[20];
  char other_colour[COLNAMLEN + 1];

  static char *line_colour[] = { "sky_blue", "orange", "red", "yellow" };

  int flip, i, ichain, iside, nchars;

  float angle_start, angle_end, line_width, r, x, x1, x2, y, y1;
  float text_size, text_width;

  /* Define starting location */
  x = origin[0];
  y = origin[1];

  /* Define text size for interaction counts */
  text_size = scale_factor * LABEL_HEIGHT;
  if (text_size < MIN_TEXT_SIZE)
    text_size = MIN_TEXT_SIZE;

  /* Initialise the flip factor (zero if an odd number of lines and
     1 if an even number) */
  flip = 0;
  if (nlines % 2 == 0)
    flip = 1;

  /* If have any interactions to plot, define line-width */
  if (nlines > 0)
    {
      pscomm_(pspage_ptr->ps_file,"Interactions");
      line_width = scale_factor * LINE_WIDTH;
      if (line_width < 1.0)
        line_width = (float) 1.0;
      pslwid_(pspage_ptr->ps_file,line_width);
    }

  /* Loop over the lines showing interaction counts */
  for (i = 0; i < NTYPES; i++)
    {
      if (interface_ptr.count[i] > 0)
        {
          /* Get the number of characters defining the count and hence
             text width */
          if (interface_ptr.count[i] < 10)
            nchars = 1;
          else if (interface_ptr.count[i] < 100)
            nchars = 2;
          else
            nchars = 3;
          text_width = (nchars + 1) * text_size * CHAR_ASPECT;

          /* Define this interaction's line colour */
          pscolb_(pspage_ptr->ps_file,line_colour[i]);

          /* Define start and end-coords of this line */
          x1 = x - scale_factor
            * (radius[0] + INTERCIRCLE_GAP / (float) 2.0);
          x2 = x - text_width / (float) 2.0;
          y1 = y - flip * line_sep;

          /* Draw the first half of the line */
          psline_(pspage_ptr->ps_file,x1,y1,x2,y1);

          /* Define start and end-coords of this line */
          x1 = x + text_width / (float) 2.0;
          x2 = x + scale_factor
            * (radius[1] + INTERCIRCLE_GAP / (float) 2.0);

          /* Draw the second half of the line */
          psline_(pspage_ptr->ps_file,x1,y1,x2,y1);

          /* Show the count of interactions above the line */
          sprintf(number_string,"%d",interface_ptr.count[i]);
          psbtxt_(pspage_ptr->ps_file,x,y1,text_size,number_string);

          /* Adjust the flip flag */
          if (flip <= 0)
            flip = abs(flip) + 2;
          else
            flip = - flip;
        }
    }

  /* Define colours and line widths */
  line_width = scale_factor * CIRCLE_BORDER_WIDTH;
  if (line_width < 1.0)
    line_width = (float) 1.0;
  pscomm_(pspage_ptr->ps_file,"Circles");
  pslwid_(pspage_ptr->ps_file,line_width);

  /* Loop over the two chains to plot */
  for (ichain = 0; ichain < 2; ichain++)
    {
      /* Set line colour for circle border */
      pscolb_(pspage_ptr->ps_file,"black");

      /* Get which side of the central line the x-coord lies on */
      iside = 2 * ichain - 1;

      /* Calculate coords of centre of circle */
      x = origin[0]
        + iside * scale_factor
        * (radius[ichain] + INTERCIRCLE_GAP / (float) 2.0);
      y = origin[1];

      /* Calculate radius of this circle */
      r = scale_factor * radius[ichain];

      /* Get this chain's colour */
      get_ppchain_colour(chains[ichain],colour);

      /* Define circle colour */
      psccol_(pspage_ptr->ps_file,colour);

      /* Plot the circle */
      pscirc_(pspage_ptr->ps_file,x,y,r);

      /* Calculate start and end angles for the wedge showing the
         extent of the interaction */
      angle_start = - angle[ichain] / (float) 2.0 + ichain * (float) 180.0;
      angle_end = angle[ichain] / (float) 2.0 + ichain * (float) 180.0;

      /* Get the other chain's colour */
      get_ppchain_colour(chains[1 - ichain],other_colour);

      /* Define wedge colour */
      psccol_(pspage_ptr->ps_file,other_colour);

      /* Plot the wedge representing the contact area */
      pspcir_(pspage_ptr->ps_file,x,y,r,angle_start,angle_end);

      /* Plot inner line of wedge */
      pspcir_(pspage_ptr->ps_file,x,y,r * INNER_FRACTION + line_width,
              angle_start,angle_end);

      /* Redefine circle colour */
      pscolb_(pspage_ptr->ps_file,colour);
      psccol_(pspage_ptr->ps_file,colour);

      /* Plot the inner circle to overwite wedge */
      pscirc_(pspage_ptr->ps_file,x,y,r * INNER_FRACTION);

      /* Show number of residues involved in interface */
      sprintf(number_string,"%dres",interface_ptr.nresid[ichain]);

      /* Write out number of residues */
      pscolb_(pspage_ptr->ps_file,"black");
      if (ichain == 0)
        {
          x1 = x + r * INNER_FRACTION * (float) 0.95;
          psrtxt_(pspage_ptr->ps_file,x1,y,text_size,number_string);
        }
      else
        {
          x1 = x - r * INNER_FRACTION * (float) 0.95;
          psltxt_(pspage_ptr->ps_file,x1,y,text_size,number_string);
        }
    }
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

plot_diagram  -  Plot the circle diagrams representing summary interface
                 statistics

***********************************************************************/

void plot_diagram(struct interface interface_ptr,
                  struct psparams psparams_ptr,
                  struct pspage *pspage_ptr,struct parameters params)
{
  char chain_name[20], colour[COLNAMLEN + 1];

  int i, ichain, nlines;

  float height, max_radius, min_radius, scale_factor, width;
  float angle[2], origin[2], radius[2], text_size, x, y;
  float line_sep;

  /* Initialise variables */
  max_radius = 0.0;
  min_radius = 0.0;

  /* Count number of lines to be drawn between the two circles */
  nlines = 0;
  for (i = 0; i < NTYPES; i++)
    if (interface_ptr.count[i] > 0)
      nlines++;

  /* Form name of output PostScript file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(pspage_ptr->ps_name,params.work_dir);
      strcat(pspage_ptr->ps_name,DIR_SEP);
    }
  strcat(pspage_ptr->ps_name,"ppcont.ps");

  /* Determine parameters for circle representations of protein chains */
  for (ichain = 0; ichain < 2; ichain++)
    {
      /* If have a chain surface area, calculate parameters */
      if (interface_ptr.chain_area[ichain] > 0.0)
        {
          /* Convert chain's total surface area into an equivalent radius */
          radius[ichain]
            = (float) sqrt(interface_ptr.chain_area[ichain] / FOUR_PI);

          /* Clip if above maximum allowed radius */
          if (radius[ichain] > MAX_RADIUS)
            radius[ichain] = MAX_RADIUS;

          /* Save highest radius */
          if (radius[ichain] > max_radius)
            max_radius = radius[ichain];

          /* Determine angle for "interaction wedge" */
          if (interface_ptr.interface_area[ichain] > 0.0)
            angle[ichain]
              = (float) (RADDEG * TWO_PI * interface_ptr.interface_area[ichain]
                         / interface_ptr.chain_area[ichain]);
          else
            angle[ichain] = 0.0;

          /* Clip is above maximum */
          if (angle[ichain] > 180.0)
            angle[ichain] = (float) 180.0;
        }

      /* Otherwise, set all to zero */
      else
        {
          radius[ichain] = 0.0;
          angle[ichain] = 0.0;
        }
    }

  /* Get the radius of the smallest circle */
  min_radius = radius[0];
  if (radius[1] < min_radius)
    min_radius = radius[1];

  /* Determine extent of picture */
  width = 4 * MAX_RADIUS + INTERCIRCLE_GAP;
  height = MAX_RADIUS + (float) 1.5 * HEADTEXT_HEIGHT;

  /* Determine the scaling factor */
  scale_factor = PICTURE_WIDTH / width;

  /* Determine the separation between interaction lines joining
     the two circles */
  line_sep = LINE_SEP;
  if (line_sep * (nlines - 1) > min_radius)
    line_sep = min_radius / (nlines - 1);
  line_sep = line_sep * scale_factor;

  /* Determine origin of central point */
  origin[0] = scale_factor * width / (float) 2.0;
  origin[1] = PICTURE_HEIGHT / (float) 2.0;

  /* Determine the picture origin */
  origin[0] = origin[0] + XMARGIN;
  origin[1] = origin[1] + YMARGIN;

  /* Initialise PostScript page variables */
  init_pspage(psparams_ptr,pspage_ptr);

  /* Start the page holding the plot */
  start_new_page(psparams_ptr,pspage_ptr);

  /* Write out chain headings */
  text_size = scale_factor * HEADTEXT_HEIGHT;
  if (text_size < MIN_HEAD_HEIGHT)
    text_size = MIN_HEAD_HEIGHT;
  x = origin[0] - scale_factor * INTERCIRCLE_GAP / (float) 2.0;
  y = origin[1] + scale_factor * max_radius + text_size;

  /* Get colour of first chain */
  get_ppchain_colour(params.chains[0],colour);
  pscolb_(pspage_ptr->ps_file,colour);

  /* Write out heading */
  sprintf(chain_name,"Chain %c",params.chains[0]);
  psrtxt_(pspage_ptr->ps_file,x,y,text_size,chain_name);
  x = x + scale_factor * INTERCIRCLE_GAP;

  /* Get colour of second chain */
  get_ppchain_colour(params.chains[1],colour);
  pscolb_(pspage_ptr->ps_file,colour);

  /* Write out heading */
  sprintf(chain_name,"Chain %c",params.chains[1]);
  psltxt_(pspage_ptr->ps_file,x,y,text_size,chain_name);
  pscolb_(pspage_ptr->ps_file,"black");

  /* Draw the two circles representing the interacting chains */
  plot_circles(interface_ptr,params.chains,scale_factor,origin,
               radius,angle,line_sep,nlines,pspage_ptr);

  /* Close the current page */
  close_page(psparams_ptr,pspage_ptr);
}
/***********************************************************************

find_chain  -  Find the given chain

***********************************************************************/

struct chain *find_chain(struct chain *first_chain_ptr,char chain)
{
  struct chain *chain_ptr, *found_chain_ptr;

  /* Initialise variables */
  found_chain_ptr = NULL;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains to find the one just read in */
  while (chain_ptr != NULL && found_chain_ptr == NULL)
    {
      /* If have a match, store the pointer */
      if (chain == chain_ptr->chain)
        found_chain_ptr = chain_ptr;

      /* Get pointer to the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }

  /* Return the matching chain pointer */
  return(found_chain_ptr);
}
/***********************************************************************

create_chain_record  -  Create and initialise a new chain record

***********************************************************************/

struct chain *create_chain_record(struct chain **fst_chain_ptr,
                                  struct chain **lst_chain_ptr,
                                  char chain,int number)
{
  struct chain *chain_ptr, *first_chain_ptr, *last_chain_ptr;

  /* Initialise chain pointers */
  chain_ptr = NULL;
  first_chain_ptr = *fst_chain_ptr;
  last_chain_ptr = *lst_chain_ptr;

  /* Allocate memory for structure to hold chain info */
  chain_ptr = (struct chain *) malloc(sizeof(struct chain));
  if (chain_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct chain\n");
      exit (1);
    }

  /* If this is the very first chain, then store pointer */
  if (first_chain_ptr == NULL)
    first_chain_ptr = chain_ptr;

  /* Add link from previous chain to the current one */
  if (last_chain_ptr != NULL)
    last_chain_ptr->next_chain_ptr = chain_ptr;
  last_chain_ptr = chain_ptr;

  /* Store the chain details */
  chain_ptr->chain = chain;
  chain_ptr->number = number;
  chain_ptr->chain_no = number;
  chain_ptr->cluster = 0;
  chain_ptr->done = FALSE;
  chain_ptr->chain_area = 0.0;
  chain_ptr->total_interface_area = 0.0;
  chain_ptr->first_interface_ptr = NULL;
  chain_ptr->last_interface_ptr = NULL;
  chain_ptr->first_wedge_ptr = NULL;
  chain_ptr->last_wedge_ptr = NULL;
  chain_ptr->ninterfaces = 0;
  chain_ptr->next_chain_ptr = NULL;
  chain_ptr->next_cluster_chain_ptr = NULL;

  /* Return current pointers */
  *fst_chain_ptr = first_chain_ptr;
  *lst_chain_ptr = last_chain_ptr;

  /* Return new chain pointer */
  return(chain_ptr);
}
/***********************************************************************

create_interface_record  -  Create and initialise a new interface record

***********************************************************************/

struct interface *create_interface_record(struct interface **fst_interface_ptr,
                                          struct interface **lst_interface_ptr,
                                          struct chain *other_chain_ptr)
{
  int i;

  struct interface *interface_ptr, *first_interface_ptr, *last_interface_ptr;

  /* Initialise interface pointers */
  interface_ptr = NULL;
  first_interface_ptr = *fst_interface_ptr;
  last_interface_ptr = *lst_interface_ptr;

  /* Allocate memory for structure to hold interface info */
  interface_ptr = (struct interface *) malloc(sizeof(struct interface));
  if (interface_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct interface\n");
      exit (1);
    }

  /* If this is the very first interface, then store pointer */
  if (first_interface_ptr == NULL)
    first_interface_ptr = interface_ptr;

  /* Add link from previous interface to the current one */
  if (last_interface_ptr != NULL)
    last_interface_ptr->next_interface_ptr = interface_ptr;
  last_interface_ptr = interface_ptr;

  /* Store the interface details */
  for (i = 0; i < 2; i++)
    {
      interface_ptr->nresid[i] = 0;
      interface_ptr->chain_area[i] = 0.0;
      interface_ptr->interface_area[i] = 0.0;
    }
  for (i = 0; i < NTYPES; i++)
    interface_ptr->count[i] = 0;
  interface_ptr->angle = 0.0;
  interface_ptr->ntypes = 0;
  interface_ptr->chain_ptr = other_chain_ptr;
  interface_ptr->next_interface_ptr = NULL;

  /* Return current pointers */
  *fst_interface_ptr = first_interface_ptr;
  *lst_interface_ptr = last_interface_ptr;

  /* Return new interface pointer */
  return(interface_ptr);
}
/***********************************************************************

read_chains  -  Read in the all the chains and their interface data from
                the interfaces.dat file

***********************************************************************/

int read_chains(char *pdbsum_dir,char *fname,
                struct chain **fst_chain_ptr,struct parameters params)
{
#define NKEY  8

  char input_line[LINELEN + 1], file_name[FILENAME_LEN];
  char tmp_line[LINELEN + 1];
  char value[LINELEN + 1];
  char chains[2], number_string[20];

  char *string_ptr;

  static char *keyword[] = {
    "PPIFACE_CHAIN", "PPIFACE_NRESID", "PPIFACE_AREA", "PPIFACE_BURIED",
    "PPIFACE_NHBONDS", "PPIFACE_NCONTACTS", "PPIFACE_SALT_BRIDGES",
    "PPIFACE_DISULPHIDES",
  };

  int ichain, ikey, key, klen, len, nchains;

  static int keylen[NKEY] = { 13, 14, 12, 14, 15, 17, 20, 19 };

  struct chain *chain_ptr[2], *first_chain_ptr, *last_chain_ptr;
  struct interface *interface_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;

  /* Initialise chain pointers */
  chain_ptr[0] = NULL;
  chain_ptr[1] = NULL;
  first_chain_ptr = NULL;
  last_chain_ptr = NULL;

  /* Form name of the interfaces.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"newsec/");
  strcat(file_name,fname);

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) ==  NULL)
    {
      printf("*** ERROR. Unable to open interfaces.dat file: %s\n",file_name);
      exit(-1);
    }

  /* Loop while reading in records from the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Check for a "copy" interface */
      if (!strncmp(input_line,"COPY_",5))
        {
          /* Remove the "COPY_" part of the keyword */
          strcpy(tmp_line,input_line+5);
          strcpy(input_line,tmp_line);
          len = len -5;
        }

      /* Initisliae flag */
      key = -1;

      /* Check whether have one of the required key words */
      for (ikey = 0; ikey < NKEY && key == -1; ikey++)
        {
          /* Check whether this is the required keyword */
          if (!strncmp(input_line,keyword[ikey],keylen[ikey]))
            key = ikey;
        }

      /* If have one of the keywords applying to a specific chain, then
         extract which chain it is */
      if (key > -1 && key < 4)
        {
          /* Get legnth of current keyword */
          klen = strlen(keyword[key]);

          /* Get which chain this is */
          number_string[0] = input_line[klen];
          number_string[1] = '\0';
          ichain = atoi(number_string) - 1;
        }
      else
        ichain = -1;

      /* Extract the field value */
      string_ptr = strstr(input_line," ");
      if (string_ptr != NULL)
        strcpy(value,string_ptr+1);
      else
        value[0] = '\0';

      /* If this is a chain identifier, find/create chain record */
      if (key == 0)
        {
          /* Initialise interface pointer */
          interface_ptr = NULL;

          /* Store the chain */
          chains[ichain] = value[0];

          /* Check whether we already have this chain */
          chain_ptr[ichain] = find_chain(first_chain_ptr,chains[ichain]);

          /* If don't have chain, create chain record */
          if (chain_ptr[ichain] == NULL)
            {
              /* Create chain record */
              chain_ptr[ichain]
                = create_chain_record(&first_chain_ptr,&last_chain_ptr,
                                      chains[ichain],nchains);

              /* Increment count of new chains */
              nchains++;
            }

          /* If this is the second chain, create a link record representing
             the interface between the two chains */
          if (ichain == 1 && chain_ptr[0] != NULL && chain_ptr[1] != NULL)
            {
              /* Create the interface record for first chain */
              interface_ptr
                = create_interface_record(&(chain_ptr[0]->first_interface_ptr),
                                          &(chain_ptr[0]->last_interface_ptr),
                                          chain_ptr[1]);

              /* Increment each chains count of interfaces */
              chain_ptr[0]->ninterfaces++;
              chain_ptr[1]->ninterfaces++;
            }
        }

      /* If have wanted keyword and valid interface record store details */
      else if (key > -1 && interface_ptr != NULL)
        {
          /* Number of interface residues */
          if (key == 1 && ichain > -1)
            interface_ptr->nresid[ichain] = atoi(value);

          /* Chain area */
          else if (key == 2 && ichain > -1 && chain_ptr[ichain] != NULL)
            {
              /* Store chain's area */
              interface_ptr->chain_area[ichain] = atof(value);
              chain_ptr[ichain]->chain_area = atof(value);
            }

          /* Interface area */
          else if (key == 3 && ichain > -1 && chain_ptr[ichain] != NULL)
            {
              /* Store interface area */
              interface_ptr->interface_area[ichain] = atof(value);

              /* Add to this chain's overall interface area */
              chain_ptr[ichain]->total_interface_area
                = chain_ptr[ichain]->total_interface_area + atof(value);
            }

          /* Number of hydrogen bonds */
          else if (key == 4)
            {
              /* Store number of hydrogen bonds */
              interface_ptr->count[HBONDS] = atoi(value);
              if (interface_ptr->count[HBONDS] > 0)
                interface_ptr->ntypes++;
            }

          /* Number of non-bonded contacts */
          else if (key == 5)
            {
              /* Store number of non-bonded contacts */
              interface_ptr->count[CONTACTS] = atoi(value);
              if (interface_ptr->count[CONTACTS] > 0)
                interface_ptr->ntypes++;
            }

          /* Number of salt bridges */
          else if (key == 6)
            {
              /* Store number of salt bridges */
              interface_ptr->count[SALT_BRIDGES] = atoi(value);
              if (interface_ptr->count[SALT_BRIDGES] > 0)
                interface_ptr->ntypes++;
            }

          /* Number of disulphide bonds */
          else if (key == 7)
            {
              /* Store number of disulphide bonds */
              interface_ptr->count[DISULPHIDES] = atoi(value);
              if (interface_ptr->count[DISULPHIDES] > 0)
                interface_ptr->ntypes++;
            }
        }
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return pointer to first chain record */
  *fst_chain_ptr = first_chain_ptr;

  /* Return number of chains read in */
  return(nchains);
}
/***********************************************************************

create_cluster_record  -  Create and initialise a new cluster record

***********************************************************************/

struct cluster *create_cluster_record(struct cluster **fst_cluster_ptr,
                                      struct cluster **lst_cluster_ptr,
                                      int number)
{
  struct cluster *cluster_ptr, *first_cluster_ptr, *last_cluster_ptr;

  /* Initialise cluster pointers */
  cluster_ptr = NULL;
  first_cluster_ptr = *fst_cluster_ptr;
  last_cluster_ptr = *lst_cluster_ptr;

  /* Allocate memory for structure to hold cluster info */
  cluster_ptr = (struct cluster *) malloc(sizeof(struct cluster));
  if (cluster_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct cluster\n");
      exit (1);
    }

  /* If this is the very first cluster, then store pointer */
  if (first_cluster_ptr == NULL)
    first_cluster_ptr = cluster_ptr;

  /* Add link from previous cluster to the current one */
  if (last_cluster_ptr != NULL)
    last_cluster_ptr->next_cluster_ptr = cluster_ptr;
  last_cluster_ptr = cluster_ptr;

  /* Store the cluster details */
  cluster_ptr->number = number;
  cluster_ptr->nchains = 0;
  cluster_ptr->first_chain_ptr = NULL;
  cluster_ptr->next_cluster_ptr = NULL;

  /* Return current pointers */
  *fst_cluster_ptr = first_cluster_ptr;
  *lst_cluster_ptr = last_cluster_ptr;

  /* Return new cluster pointer */
  return(cluster_ptr);
}
/***********************************************************************

count_clusters  -  Count and store clusters of interacting chains

***********************************************************************/

int count_clusters(struct chain *first_chain_ptr,int nchains,
                   struct cluster **fst_cluster_ptr)
{
  int nassigned, nclusters;
  int assigned, done, new;

  struct chain *chain_ptr, *last_chain_ptr, *other_chain_ptr;
  struct cluster *cluster_ptr, *first_cluster_ptr, *last_cluster_ptr;
  struct interface *interface_ptr;

  /* Initialise variables */
  nassigned = 0;
  nclusters = 0;

  /* Initialise pointers */
  cluster_ptr = NULL;
  first_cluster_ptr = NULL;
  last_cluster_ptr = NULL;

  /* Loop until all chains have been assigned to clusters */
  while (nassigned < nchains)
    {
      /* Increment cluster number */
      nclusters++;

      /* Create a new cluster record */
      cluster_ptr
        = create_cluster_record(&first_cluster_ptr,&last_cluster_ptr,
                                nclusters);

      /* Get pointer to the first chain */
      chain_ptr = first_chain_ptr;
      last_chain_ptr = NULL;
      done = FALSE;

      /* Loop over all the chains until find one not yet assigned to
         a cluster */
      while (chain_ptr != NULL && done == FALSE)
        {
          /* If chain not assigned to a cluster, assign it to current
             cluster */
          if (chain_ptr->cluster == 0)
            {
              /* Assign chain and set flag */
              chain_ptr->cluster = nclusters;
              done = TRUE;

              /* Mark this as being the first and last chain of the
                 current cluster */
              cluster_ptr->first_chain_ptr = chain_ptr;
              last_chain_ptr = chain_ptr;
 
              /* Renumber this chain */
              chain_ptr->number = cluster_ptr->nchains;
              chain_ptr->chain_no = nassigned;

             /* Increment count of chains in this cluster */
              cluster_ptr->nchains++;
              nassigned++;
            }

          /* Get the next chain */
          chain_ptr = chain_ptr->next_chain_ptr;
        }

      /* Initialise flag */
      new = TRUE;

      /* Loop while new chains are being assigned to the current cluster */
      while (new == TRUE)
        {
          /* Initialise flag */
          new = FALSE;

          /* Get pointer to the first chain of this cluster */
          chain_ptr = cluster_ptr->first_chain_ptr;

          /* Loop over all the chains to find those belonging to the
             current cluster which haven't yet has their interaction
             partners added */
          while (chain_ptr != NULL)
            {
              /* If chain belongs to current cluster, add its partners
                 to the cluster */
              if (chain_ptr->cluster == nclusters &&
                  chain_ptr->done == FALSE)
                {
                  /* Get the first of this chain's interactions */
                  interface_ptr = chain_ptr->first_interface_ptr;

                  /* Loop over the interfaces */
                  while (interface_ptr != NULL)
                    {
                      /* Get the other chain */
                      other_chain_ptr = interface_ptr->chain_ptr;

                      /* If chain not already assigned, add it to the 
                         current cluster */
                      if (other_chain_ptr->cluster == 0)
                        {
                          /* Assign chain to cluster */
                          other_chain_ptr->cluster = nclusters;

                          /* Add into linked list for this cluster */
                          last_chain_ptr->next_cluster_chain_ptr
                            = other_chain_ptr;
                          last_chain_ptr = other_chain_ptr;

                          /* Renumber this chain */
                          other_chain_ptr->number = cluster_ptr->nchains;
                          other_chain_ptr->chain_no = nassigned;

                          /* Increment count of chains in this cluster */
                          cluster_ptr->nchains++;
                          nassigned++;

                          /* Set flag that new chain added */
                          new = TRUE;
                        }

                      /* Get the next interface */
                      interface_ptr = interface_ptr->next_interface_ptr;
                    }

                  /* Mark this chain as done */
                  chain_ptr->done = TRUE;
                }

              /* Otherwise, if chain not assigned to any cluster, check
                 to see if any of its partners are in the current
                 cluster */
              else if (chain_ptr->cluster == 0)
                {
                  /* Initialie flag */
                  assigned = FALSE;

                  /* Get the first of this chain's interactions */
                  interface_ptr = chain_ptr->first_interface_ptr;

                  /* Loop over the interfaces */
                  while (interface_ptr != NULL && assigned == FALSE)
                    {
                      /* Get the other chain */
                      other_chain_ptr = interface_ptr->chain_ptr;

                      /* If other chain belongs to the current cluster
                         need to add current chain to the cluster */
                      if (other_chain_ptr->cluster == nclusters)
                        {
                          /* Assign current chain to cluster */
                          chain_ptr->cluster = nclusters;

                          /* Add into linked list for this cluster */
                          last_chain_ptr->next_cluster_chain_ptr
                            = chain_ptr;
                          last_chain_ptr = chain_ptr;

                          /* Renumber this chain */
                          chain_ptr->number = cluster_ptr->nchains;
                          chain_ptr->chain_no = nassigned;

                          /* Increment count of chains in this cluster */
                          cluster_ptr->nchains++;
                          nassigned++;

                          /* Set flags */
                          new = TRUE;
                          assigned = TRUE;
                        }

                      /* Get the next interface */
                      interface_ptr = interface_ptr->next_interface_ptr;
                    }
                }

              /* Get the next chain */
              chain_ptr = chain_ptr->next_chain_ptr;
            }
        }
    }

  /* Return pointer to first cluster record */
  *fst_cluster_ptr = first_cluster_ptr;
  printf("Number of clusters = %d\n",nclusters);

  /* Return number of clusters */
  return(nclusters);
}
/***********************************************************************

calc_separation  -  Calculate separation between circles for current
                    global radius

***********************************************************************/

float calc_separation(float *radius,float global_radius,int nchains)
{
  int ichain;

  float separation;
  float chord, chord_sum, circumference, sinang, theta;

  /* Initialise variables */
  separation = 0.0;

  /* Get the circumference corresponding to the current global radius */
  circumference = TWO_PI * global_radius;

  /* Initialise sum of chords subtended by the chain circles */
  chord_sum = 0.0;

  /* Loop through all chains to get the total chord that their circles
     subtend */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Get the length of the chord subtended by the circle
         representing the current chain */
      sinang = radius[ichain] / (2 * global_radius);

      /* Get the corresponding angle */
      theta = 2 * (float) asin((double) sinang);

      /* Calculate the corresponding chord length */
      chord = circumference * theta / TWO_PI;

      /* Add to accumuloated chord sum */
      chord_sum = chord_sum + 2 * chord;
    }

  /* Calculate adjusted minimum separation from the chord sums */
  separation = (circumference - chord_sum) / nchains;

  /* Return current separation */
  return(separation);
}
/***********************************************************************

calc_circle_coords  -  Calculate coordinates of circles on the global
                       circle

***********************************************************************/

float calc_circle_coords(float *radius,float *circle_coords,int nchains,
                         float global_radius)
{
#define MAX_LOOPS  10

  int ichain, jchain, loop;
  int done;

  float chord_sum, circumference, separation;
  float dtheta, last_dtheta, sinang, theta;
  float dist, dist_sqrd, gap_angle, min_sep, x, x1, x2, y, y1, y2;

  /* Initialise variables */
  done = FALSE;
  loop = 0;

  /* Loop until have a reasonable placement of the circles on the global
     circle */
  while (done == FALSE && loop < MAX_LOOPS)
    {
      /* Initialise variables */
      dtheta = 0.0;
      last_dtheta = 0.0;
      //theta = - PI_BY_2;
      theta = - PI;

      /* Calculate current separation */
      separation = calc_separation(radius,global_radius,nchains);

      /* Get the current circumference */
      circumference = TWO_PI * global_radius;

      /* Get the angle corresponding to the current separation */
      gap_angle = TWO_PI * separation / circumference;

      /* Initialise chord sum */
      chord_sum = 0.0;

      /* Loop over the circles representing each chain to place on 
         circumference of global circle */
      for (ichain = 0; ichain < nchains; ichain++)
        {
          /* Calculate angle corresponding to current circle radius */
          sinang = radius[ichain] / (2 * global_radius);
          dtheta = 2 * (float) asin((double) sinang);

          /* Calculate the angular location of this circle */
          if (ichain > 0)
            {
              /* Add previous angle and one gap to last angle
                 position */
              theta = theta + last_dtheta + gap_angle;

              /* Add current angle */
              theta = theta + dtheta;
            }

          /* Save current angle */
          last_dtheta = dtheta;

          /* Convert angular position to x- and y-coords */
          x = global_radius * (float) cos(theta);
          y = - global_radius * (float) sin(theta);

          /* Store the coordinates of this circle */
          circle_coords[2 * ichain] = x;
          circle_coords[2 * ichain + 1] = y;
        }

      /* Initialise minimum separation */
      min_sep = global_radius;

      /* Loop over all circles again to get the minimum separation */
      for (ichain = 0; ichain < nchains; ichain++)
        {
          /* Extract this circle's coords */
          x1 = circle_coords[2 * ichain];
          y1 = circle_coords[2 * ichain + 1];

          /* Loop over other circles to calculate separations */
          for (jchain = ichain + 1; jchain < nchains; jchain++)
            {
              /* Extract this circle's coords */
              x2 = circle_coords[2 * jchain];
              y2 = circle_coords[2 * jchain + 1];

              /* Calculate squared distance between the centres of these
                 two circles */
              dist_sqrd = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

              /* Calculate distance between circumferences */
              dist = (float) sqrt(dist_sqrd)
                - radius[ichain] - radius[jchain];

              /* If this is the smallest so far, store */
              if (dist < min_sep)
                min_sep = dist;
            }
        }

      /* If separation is above minimum we are done */
      if (min_sep >= INTERCIRCLE_GAP)
        done = TRUE;

      /* Otherwise, increase global radius size and recompute */
      else
        global_radius
          = global_radius + INTERCIRCLE_GAP - min_sep + (float) 1.0;

      /* Increment loop counter */
      loop++;
    }

  /* Return the global radius */
  return(global_radius);
}
/***********************************************************************

calc_global_radius  -  Calculate overall radius of circle on which all
                       the chain-circles are to be placed

***********************************************************************/

float calc_global_radius(struct chain *first_chain_ptr,int nchains,
                         float *radius,float *mx_radius,
                         float *circle_coords)
{
  int ichain;

  float global_radius, max_radius[2], separation, tot_radius;

  struct chain *chain_ptr;

  /* Initialise variables */
  global_radius = 0.0;
  max_radius[0] = 0.0;
  max_radius[1] = 0.0;
  tot_radius = 0.0;
  separation = INTERCIRCLE_GAP;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;
  ichain = 0;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get this chain's radius */
      radius[ichain] = (float) sqrt(chain_ptr->chain_area / FOUR_PI);

      /* If largest so far, store */
      if (radius[ichain] > max_radius[0])
        {
          max_radius[1] = max_radius[0];
          max_radius[0] = radius[ichain];
        }

      /* If second-largest so far, store */
      else if (radius[ichain] > max_radius[1])
        {
          max_radius[1] = radius[ichain];
        }

      /* Add to total area */
      tot_radius = tot_radius + radius[ichain];

      /* Get the next chain */
      chain_ptr = chain_ptr->next_cluster_chain_ptr;

      /* Increment chain counter */
      ichain++;
    }

  /* Estimate lower bound for the global radius */
  global_radius = (tot_radius + nchains * separation) / TWO_PI;

  /* If this is less than the sum of the two largest radii plus the
     minimum separation, adjust */
  if (global_radius < max_radius[0] + max_radius[1] + separation)
    global_radius = max_radius[0] + max_radius[1] + separation;

  /* Calculate optimum coordinates of circles on the global circle */
  global_radius
    = calc_circle_coords(radius,circle_coords,nchains,global_radius);

  /* Return maximum radius */
  *mx_radius = max_radius[0];

  /* Return the calculated global radius */
  return(global_radius);
}
/***********************************************************************

create_wedge_record  -  Create and initialise a new wedge record

***********************************************************************/

struct wedge *create_wedge_record(struct wedge **fst_wedge_ptr,
                                  struct wedge **lst_wedge_ptr,
                                  struct interface *interface_ptr,
                                  int first_chain,struct chain *chain_ptr,
                                  struct chain *other_chain_ptr,
                                  double chain_area,
                                  double interface_area)
{
  int i;

  float angle;

  struct wedge *wedge_ptr, *first_wedge_ptr, *last_wedge_ptr;

  /* Initialise wedge pointers */
  wedge_ptr = NULL;
  first_wedge_ptr = *fst_wedge_ptr;
  last_wedge_ptr = *lst_wedge_ptr;

  /* Allocate memory for structure to hold wedge info */
  wedge_ptr = (struct wedge *) malloc(sizeof(struct wedge));
  if (wedge_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct wedge\n");
      exit (1);
    }

  /* If this is the very first wedge, then store pointer */
  if (first_wedge_ptr == NULL)
    first_wedge_ptr = wedge_ptr;

  /* Add link from previous wedge to the current one */
  if (last_wedge_ptr != NULL)
    last_wedge_ptr->next_wedge_ptr = wedge_ptr;
  last_wedge_ptr = wedge_ptr;

  /* Store the wedge details */
  for (i = 0; i < 3; i++)
    wedge_ptr->angle[i] = 0.0;
  wedge_ptr->angle[CENTRE] = 0.0;
  wedge_ptr->interface_ptr = interface_ptr;
  wedge_ptr->first_chain = first_chain;
  wedge_ptr->chain_ptr = chain_ptr;
  wedge_ptr->other_chain_ptr = other_chain_ptr;
  wedge_ptr->next_wedge_ptr = NULL;

  /* Calculate initial start and end angles for this wedge */
  if (chain_area > 0.0)
    angle = (float) (RADDEG * TWO_PI * interface_area / chain_area);
  else
    angle = 0.0;
  wedge_ptr->angle[START] = - angle / 2;
  wedge_ptr->angle[END] = angle / 2;

  /* Return current pointers */
  *fst_wedge_ptr = first_wedge_ptr;
  *lst_wedge_ptr = last_wedge_ptr;

  /* Return new wedge pointer */
  return(wedge_ptr);
}
/***********************************************************************

compute_interaction_wedges  -  Calculate start- and end-angles for each
                               interaction wedge
  
***********************************************************************/

void compute_interaction_wedges(struct chain *first_chain_ptr)
{
  int ichain;

  struct chain *chain_ptr, *other_chain_ptr;
  struct interface *interface_ptr;
  struct wedge *wedge_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;
  ichain = 0;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get the first of this chain's interactions */
      interface_ptr = chain_ptr->first_interface_ptr;

      /* Loop over the interfaces */
      while (interface_ptr != NULL)
        {
          /* Get the other chain involved in this interface */
          other_chain_ptr = interface_ptr->chain_ptr ;

          /* Create a new wedge record */
          wedge_ptr
            = create_wedge_record(&(chain_ptr->first_wedge_ptr),
                                  &(chain_ptr->last_wedge_ptr),
                                  interface_ptr,TRUE,
                                  chain_ptr,other_chain_ptr,
                                  interface_ptr->chain_area[0],
                                  interface_ptr->interface_area[0]);

          /* Repeat for the other chain */
          wedge_ptr
            = create_wedge_record(&(other_chain_ptr->first_wedge_ptr),
                                  &(other_chain_ptr->last_wedge_ptr),
                                  interface_ptr,FALSE,
                                  other_chain_ptr,chain_ptr,
                                  interface_ptr->chain_area[1],
                                  interface_ptr->interface_area[1]);

          /* Get the next interface */
          interface_ptr = interface_ptr->next_interface_ptr;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;

      /* Increment chain counter */
      ichain++;
    }
}
/***********************************************************************

chain_clash_energy  -  Calculate contribution to "energy" of any chain
                       overlaps

***********************************************************************/

float chain_clash_energy(int nchains,float *radius,float *circle_coords)
{
  int ichain, jchain;

  float energy, clash_energy;
  float dist, dist_sqrd, r, x1, x2, y1, y2;

  /* Initialise energy */
  energy = 0.0;

  /* Loop over all the chains */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Extract this circle's coords */
      x1 = circle_coords[2 * ichain];
      y1 = circle_coords[2 * ichain + 1];

      /* Loop over all other chains to check for clashes */
      for (jchain = ichain + 1; jchain < nchains; jchain++)
        {
          /* Extract this circle's coords */
          x2 = circle_coords[2 * jchain];
          y2 = circle_coords[2 * jchain + 1];

          /* Calculate squared distance between the centres of these
             two circles */
          dist_sqrd = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

          /* Calculate distance between circumferences */
          dist = (float) sqrt(dist_sqrd)
            - radius[ichain] - radius[jchain];

          /* Define push radius */
          r = INTERCIRCLE_GAP + radius[ichain] + radius[jchain]; 

          /* If distance is within the inter-circle spacing, penalise
             by adding extra energy term */
          clash_energy = 0.0;
          if (dist < r && dist >= INTERCIRCLE_GAP)
            clash_energy = - (float) fabs(r - dist);
          else if (dist < INTERCIRCLE_GAP)
            clash_energy = - (r - dist) * (r - dist);
          if (dist < 0.0)
            clash_energy = - clash_energy * clash_energy;

          /* Add energy to ovewrall overlap energy */
          energy = energy + clash_energy;
        }
    }

  /* Return computed energy */
  return(energy);
}
/***********************************************************************

overlap_energy  -  Calculate contribution to "energy" of any interaction
                   lines overlapping any circles

***********************************************************************/

float overlap_energy(struct chain *first_chain_ptr,int nchains,
                     float *radius,float *circle_coords,float x1,
                     float y1,float x2,float y2,int ichain,int jchain)
{
  int kchain;

  float energy, olap_energy;
  float dist, minx, miny, maxx, maxy, r, x, y;
  float modv, v1, v2, w1, w2;

  struct chain *chain_ptr;

  /* Initialise energy */
  energy = 0.0;

  /* Get maximum and minimum coords defining the interaction line */
  minx = x1;
  maxx = x1;
  miny = y1;
  maxy = y1;
  if (x2 < minx)
    minx = x2;
  if (x2 > maxx)
    maxx = x2;
  if (y2 < miny)
    miny = y2;
  if (y2 > maxy)
    maxy = y2;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;
  kchain = 0;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Only consider if not one of the two involved in the given
         interface */
      if (kchain != ichain && kchain != jchain)
        {
          /* Extract this circle's coords */
          x = circle_coords[2 * kchain];
          y = circle_coords[2 * kchain + 1];

          /* Get effective radius */
          r = radius[kchain] + INTERCIRCLE_GAP;

          /* Check whether this circle is in range */
          if (x - r < maxx && x + r > minx &&
              y - r < maxy && y + r > miny)
            {
              /* Calculate normalised vector along interaction line */
              v1 = x2 - x1;
              v2 = y2 - y1;
              modv = (float) sqrt(v1 * v1 + v2 * v2);
              if (modv > 0.0)
                {
                  v1 = v1 / modv;
                  v2 = v2 / modv;
                }

              /* Get vector from circle to one end to the interaction
                 line */
              w1 = x - x1;
              w2 = y - y1;

              /* Calculate cross product between the two vectors to get
                 the distance from point to line */
              dist = (float) fabs(v1 * w2 - v2 * w1) - radius[kchain];

              /* If distance is within the inter-circle spacing,
                 penalise by adding extra energy term */
              olap_energy = 0.0;
              if (dist < r && dist >= INTERCIRCLE_GAP)
                olap_energy = - (float) fabs(r - dist);
              else if (dist < INTERCIRCLE_GAP)
                olap_energy = - (r - dist) * (r - dist);
              if (dist < 0.0)
                olap_energy = - olap_energy * olap_energy;

              /* Add energy to ovewrall overlap energy */
              energy = energy + olap_energy;
            }
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_cluster_chain_ptr;

      /* Increment chain counter */
      kchain++;
    }

  /* Return computed energy */
  return(energy);
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

calc_interaction_angle  -  Calculate angle made by given interaction
                           line

***********************************************************************/

float calc_interaction_angle(float x1,float y1,float x2,float y2)
{
  float angle, x, xleft, xright, y, yleft, yright;

  /* Initialise variables */
  angle = 0.0;

  /* Determine whyich point is on the left and which on the right */
  if (x1 < x2)
    {
      xleft = x1;
      yleft = y1;
      xright = x2;
      yright = y2;
    }
  else
    {
      xleft = x2;
      yleft = y2;
      xright = x1;
      yright = y1;
    }

  /* Calculate angle made by the line */
  x = xright - xleft;
  y = yright - yleft;
  angle = get_angle(x,y);

  /* Adjust, if necessary */
  if (x1 >= x2)
    angle = angle - PI;

  /* Return calculated angle */
  return(angle);
}
/***********************************************************************

wedge_clash_energy  -  Calculate contribution to "energy" from any
                       clashing wedges

***********************************************************************/

float wedge_clash_energy(struct chain *first_chain_ptr)
{
  int end_angle, other_end_angle, other_start_angle, start_angle;
  int add, i, overlap;

  float angle;
  float energy, wedge_energy;

  struct chain *chain_ptr;
  struct interface *interface_ptr;
  struct wedge *wedge_ptr, *other_wedge_ptr;

  /* Initialise energy */
  energy = 0.0;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get pointer to first of this circle's interface wedges */
      wedge_ptr = chain_ptr->first_wedge_ptr;

      /* Loop over the wedges to plot */
      while (wedge_ptr != NULL)
        {
          /* Get the interface that this wedge corresponds to */
          interface_ptr = wedge_ptr->interface_ptr;

          /* Get interaction angle */
          if (wedge_ptr->first_chain == TRUE)
            angle = RADDEG * interface_ptr->angle;
          else
            angle = RADDEG * (interface_ptr->angle - PI);

          /* Calculate start and end angles for this wedge */
          start_angle = (int) (wedge_ptr->angle[START] + angle);
          end_angle = (int) (wedge_ptr->angle[END] + angle);

          /* Put into the range 0-360 */
          start_angle = start_angle % 360;
          end_angle = end_angle % 360;
          if (start_angle < 0)
            {
              start_angle = start_angle + 360;
              end_angle = end_angle + 360;
            }

          /* Get the next wedge */
          other_wedge_ptr = wedge_ptr->next_wedge_ptr;

          /* Loop over the wedges to plot */
          while (other_wedge_ptr != NULL)
            {
              /* Get the interface that this wedge corresponds to */
              interface_ptr = other_wedge_ptr->interface_ptr;

              /* Get interaction angle */
              if (other_wedge_ptr->first_chain == TRUE)
                angle = RADDEG * interface_ptr->angle;
              else
                angle = RADDEG * (interface_ptr->angle - PI);

              /* Calculate start and end angles for this wedge */
              other_start_angle
                = (int) (other_wedge_ptr->angle[START] + angle);
              other_end_angle
                = (int) (other_wedge_ptr->angle[END] + angle);

              /* Put into the range 0-360 */
              other_start_angle = other_start_angle % 360;
              other_end_angle = other_end_angle % 360;
              if (other_start_angle < 0)
                {
                  other_start_angle = other_start_angle + 360;
                  other_end_angle = other_end_angle + 360;
                }

              /* Initialise for overlap check */
              add = -360;

              /* Loop over the three overlap check */
              for (i = 0; i < 3; i++)
                {
                  /* Check whether the two wedges overlap */
                  if (start_angle + add > other_start_angle &&
                      start_angle + add < other_end_angle)
                    {
                      /* Calculate overlap */
                      overlap = other_end_angle - (start_angle + add);

                      /* Convert to an overlap energy */
                      wedge_energy = (float) (-10 * overlap);

                      /* Add to overall energy */
                      energy = energy + wedge_energy;
                    }

                  /* Check for overlap in other direction */
                  if (end_angle + add > other_start_angle &&
                      end_angle + add < other_end_angle)
                    {
                      /* Calculate overlap */
                      overlap = end_angle + add - other_start_angle;

                      /* Convert to an overlap energy */
                      wedge_energy = (float) (-10 * overlap);

                      /* Add to overall energy */
                      energy = energy + wedge_energy;
                    }

                  /* Increment the added angle */
                  add = add + 360;
                }

              /* Get the next wedge */
              other_wedge_ptr = other_wedge_ptr->next_wedge_ptr;
            }

          /* Get the next wedge */
          wedge_ptr = wedge_ptr->next_wedge_ptr;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_cluster_chain_ptr;
    }


  /* Return computed energy */
  return(energy);
 }
/***********************************************************************

calc_energy  -  Calculate "energy" of current arrangement of circles

***********************************************************************/

float calc_energy(struct chain *first_chain_ptr,int nchains,
                  float *radius,float *circle_coords,int *include)
{
  int ichain, jchain;

  float energy, olap_energy, sep_energy, wedge_energy;
  float dist, dist_sqrd, x1, x2, y1, y2;

  struct chain *chain_ptr, *other_chain_ptr;
  struct interface *interface_ptr;

  /* Initialise energy */
  energy = 0.0;

  /* Calculate the chain-clash energy */
  if (include[CLASH_ENERGY] == TRUE)
    energy = chain_clash_energy(nchains,radius,circle_coords);

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get the chain number */
      ichain = chain_ptr->number;

      /* Extract this circle's coords */
      x1 = circle_coords[2 * ichain];
      y1 = circle_coords[2 * ichain + 1];

      /* Get the first of this chain's interactions */
      interface_ptr = chain_ptr->first_interface_ptr;

      /* Loop over the interfaces */
      while (interface_ptr != NULL)
        {
          /* Get the other chain involved in this interface */
          other_chain_ptr = interface_ptr->chain_ptr;

          /* Get this chain */
          jchain = other_chain_ptr->number;

          /* Extract this circle's coords */
          x2 = circle_coords[2 * jchain];
          y2 = circle_coords[2 * jchain + 1];

          /* Calculate term for obtaining target distance between
             interacting chains */
          if (include[SEPARATION_ENERGY] == TRUE)
            {
              /* Calculate squared distance between the centres of these
                 two circles */
              dist_sqrd = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

              /* Calculate distance between circumferences */
              dist
                = (float) sqrt(dist_sqrd) - radius[ichain] - radius[jchain];

              /* Calculate energy from this separation */
              sep_energy = (float) fabs (dist - INTERCIRCLE_GAP);
              sep_energy = - sep_energy * sep_energy;
              if (dist < INTERCIRCLE_GAP)
                sep_energy = - sep_energy * sep_energy;

              /* Add to total energy */
              energy = energy + sep_energy;
            }

          /* Calculate overlap energy between the lines representing
             this interface and any other circles */
          if (include[OVERLAP_ENERGY] == TRUE)
            {
              olap_energy
                = overlap_energy(first_chain_ptr,nchains,radius,circle_coords,
                                 x1,y1,x2,y2,ichain,jchain);

              /* Add to total energy */
              energy = energy + olap_energy;
            }

          /* Calculate the angle the interaction line makes */
          interface_ptr->angle = calc_interaction_angle(x1,y1,x2,y2);

          /* Get the next interface */
          interface_ptr = interface_ptr->next_interface_ptr;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_cluster_chain_ptr;
    }

  /* Calculate contribution from any wedge clashes */
  if (include[WEDGE_ENERGY] == TRUE)
    {
      wedge_energy = wedge_clash_energy(first_chain_ptr);

      /* Add wedge clash energy to total energy */
      energy = energy + wedge_energy;
    }

  /* Return computed energy */
  return(energy);
}
/***********************************************************************

sort_chains  -  Bubble-sort the chains by their number of interfaces

***********************************************************************/

void sort_chains(int *ninterfaces,int *index,int nchains)
{
  int ichain, ipos, jchain;
  int swapped;

  /* Initialise variables */
  swapped = TRUE;

  /* Loop until chains in required order */
  while (swapped == TRUE)
    {
      /* Initialise swapped flag */
      swapped = FALSE;

      /* Loop over all the chains */
      for (ipos = 0; ipos < nchains - 1; ipos++)
        {
          /* Get the corresponding chain */
          ichain = index[ipos];

          /* Get the next chain in the current list */
          jchain = index[ipos + 1];

          /* Determine whether the chains are in the correct order */
          if (ninterfaces[jchain] > ninterfaces[ichain])
            {
              /* Chains in the wrong order, so need to swap */
              index[ipos + 1] = ichain;
              index[ipos] = jchain;

              /* Set flag that swap made */
              swapped = TRUE;
            }
        }
    }
}
/***********************************************************************

place_in_grid  -  Place current chain in the grid, biasing its location
                  to grid-points neighbouring chains it interacts with

***********************************************************************/

int place_in_grid(int ichain,int *grid,int *map,int *interface_list,
                  int *first_interface,int *ninterfaces,
                  int *top_location,int *is,int ngrid,int nside)
{
  int i, ipos, j, jchain, jpos, k, l, iseed, nchains, score, top_score;
  int kstart, kend, lstart, lend, ntop;
  int edge1, edge2, interacts;

  /* Initialise variables */
  iseed = *is;
  ntop = 0;
  top_score = -1;

  /* Initialise the map */
  for (ipos = 0; ipos < ngrid; ipos++)
    {
      map[ipos] = 0;
    }

  /* Use chains currently placed in the grid to prepare a map of 
     favourable locations for the current chain */
  for (ipos = 0; ipos < ngrid; ipos++)
    {
      /* If current grid location is occupied, get the chain */
      if (grid[ipos] > -1)
        {
          /* Get the current chain identifier */
          jchain = grid[ipos];

          /* Make the corresponding point in the map inaccessible */
          map[ipos] = -1;

          /* Initialise flag */
          interacts = FALSE;

          /* Get the first of the chains that this one interacts with
             and their number */
          jpos = first_interface[jchain];
          nchains = ninterfaces[jchain];

          /* Check whether the current chain interacts with this one */
          for (i = 0; i < nchains && interacts == FALSE; i++)
            {
              /* If have the right chain, set flag */
              if (interface_list[jpos] == ichain)
                interacts = TRUE;

              /* Increment chain position */
              jpos++;
            }

          /* If our current chain interacts with this chain, then all
             neighbouring points that are blank are favourable */
          if (interacts == TRUE)
            {
              /* Get the corresponding subscripts of this grid point */
              i = ipos / nside;
              j = ipos % nside;

              /* Get start and end loop subscripts */
              kstart = i - 1;
              if (kstart < 0)
                kstart = 0;
              kend = i + 1;
              if (kend > nside - 1)
                kend = nside - 1;
              lstart = j - 1;
              if (lstart < 0)
                lstart = 0;
              lend = j + 1;
              if (lend > nside - 1)
                lend = nside - 1;

              /* Loop over all surrounding grid points to see which are
                 available */
              for (k = kstart; k < kend + 1; k++)
                {
                  /* Check whether this is an edge grid point */
                  if (k == i - 1 || k == i + 1)
                    edge1 = TRUE;
                  else
                    edge1 = FALSE;

                  /* Loop over columns */
                  for (l = lstart; l < lend + 1; l++)
                    {
                      /* Check whether this is an edge grid point */
                      if (l == j - 1 || l == j + 1)
                        edge2 = TRUE;
                      else
                        edge2 = FALSE;

                      /* Get this grid point's subscript */
                      jpos = k * nside + l;

                      /* If this is a valid location, update map score */
                      if (map[jpos] > -1)
                        {
                          /* Determine update score */
                          if (edge1 == TRUE && edge2 == TRUE)
                            score = 1;
                          else
                            score = 2;

                          /* Update this map position */
                          map[jpos] = map[jpos] + score;
                        }
                    }
                }
            }
        }
    }

  /* Determine the top score in the map */
  for (ipos = 0; ipos < ngrid; ipos++)
    if (map[ipos] > top_score)
      top_score = map[ipos];

  /* Determine locations of the highest-scoring grid positions */
  for (ipos = 0; ipos < ngrid; ipos++)
    {
      /* If this is one of the top-scoring positions, store its location */
      if (map[ipos] == top_score)
        {
          /* Store location */
          top_location[ntop] = ipos;

          /* Increment count of top locations */
          ntop++;
        }
    }

  /* Return current random number seed */
  *is = iseed;

  /* Return number of top locations */
  return(ntop);
}
/***********************************************************************

place_chains  -  Place circles representing chains to give a clear
                 representation of which w=chains interact with one
                 another

***********************************************************************/

int place_chains(struct chain *first_chain_ptr,int nchains,
                 float *radius,float *wd,float *ht,float *circle_coords)
{
  int i, ichain, iorder, ipos, iseed, j, jchain, loop;
  int interface, ngrid, nloops, nside, ntop;
  int abort, done, have_placement, placed;

  static int include[NTERMS] = { TRUE, TRUE, TRUE, TRUE };

  int *grid, *best_grid, *first_interface, *index, *interface_list;
  int *map, *ninterfaces, *top_location;

  float best_energy, energy, grid_spacing, height, max_radius, width;
  float minx, miny, maxx, maxy, midx, midy, x, y;

  struct chain *chain_ptr, *other_chain_ptr;
  struct wedge *wedge_ptr;

  /* Initialise variables */
  best_energy = - (float) HUGE_VAL;
  done = FALSE;
  have_placement = FALSE;
  height = 0.0;
  interface = 0;
  loop = 0;
  max_radius = 0.0;
  nloops = 2 * nchains;
  nloops = 500;
  width = 0.0;
  grid_spacing = 0.0;

  /* Initialise random number generator */
  iseed = -1;
  get_random_number(&iseed,RANDOM_START);

  /* Allocate memory to store the number of interfaces each chain has and
     an index giving their order sorted by this number */
  index = (int *) malloc(nchains * sizeof(int));
  ninterfaces = (int *) malloc(nchains * sizeof(int));
  first_interface = (int *) malloc(nchains * sizeof(int));
  interface_list = (int *) malloc(nchains * nchains * sizeof(int));

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains to get all the radii */
  while (chain_ptr != NULL)
    {
      /* Get the chain number */
      ichain = chain_ptr->number;
/* debug */
      printf("Chain %d  Area = %8.3f\n",ichain,chain_ptr->chain_area);
  fflush(stdout);
/* debug */

      /* Get this chain's radius */
      radius[ichain] = (float) sqrt(chain_ptr->chain_area / FOUR_PI);

      /* If largest so far, store */
      if (radius[ichain] > max_radius)
        max_radius = radius[ichain];

      /* Store the chain's number of interfaces and its current order */
      ninterfaces[ichain] = 0;
      index[ichain] = ichain;

      /* Store location of first interface chain in the interfaces list */
      first_interface[ichain] = interface;

      /* Get the first of this chain's interactions */
      wedge_ptr = chain_ptr->first_wedge_ptr;

      /* Loop over the interfaces */
      while (wedge_ptr != NULL)
        {
          /* Get the other chain involved in the interface */
          if (chain_ptr == wedge_ptr->chain_ptr)
            other_chain_ptr = wedge_ptr->other_chain_ptr;
          else
            other_chain_ptr = wedge_ptr->chain_ptr;
          jchain = other_chain_ptr->number;

          /* Store this chain identifier */
          interface_list[interface] = jchain;
          interface++;

          /* Increment count of this chain's interfaces */
          ninterfaces[ichain]++;

          /* Get the next interface */
          wedge_ptr = wedge_ptr->next_wedge_ptr;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_cluster_chain_ptr;
    }

  /* Sort the chains according to the number of interfaces they have */
  if (nchains > 1)
    sort_chains(ninterfaces,index,nchains);

  /* Calculate the gap separation for the grid onto which the circles
     will be placed */
  grid_spacing = 2 * max_radius + INTERCIRCLE_GAP;

  /* Determine size of grid onto which circles will be placed */
  nside = (nchains + 2) / 2;
  ngrid = nside * nside;

  /* Allocate memory for the grid and best grid */
  grid = (int *) malloc(ngrid * sizeof(int));
  best_grid = (int *) malloc(ngrid * sizeof(int));
  map = (int *) malloc(ngrid * sizeof(int));
  top_location = (int *) malloc(ngrid * sizeof(int));

  /* Initialise best grid */
  for (i = 0; i < ngrid; i++)
    best_grid[i] = -1;

  /* Loop until have a reasonable solution */
  while (done == FALSE && loop < nloops)
    {
      /* Initialise grid */
      for (i = 0; i < ngrid; i++)
        grid[i] = -1;
      abort = FALSE;

      /* Loop through all the chains to place them on the grid */
      for (iorder = 0; iorder < nchains && abort == FALSE; iorder++)
        {
          /* Get this chain */
          ichain = index[iorder];

          /* Place it in the grid */
          ntop
            = place_in_grid(ichain,grid,map,interface_list,first_interface,
                            ninterfaces,top_location,&iseed,ngrid,nside);

          /* Initialise flag */
          placed = FALSE;

          /* If have no top locations where chain can be placed, have
             run into a dead-end and cannot complete this configuration */
          if (ntop == 0)
            abort = TRUE;

          /* Loop until this circle has been placed */
          while (abort == FALSE && placed == FALSE)
            {
              /* Generate random grid location */
              ipos = (int) (ntop * get_random_number(&iseed,RANDOM_START));

              /* If this place is empty, place the circle here */
              if (ipos < ntop)
                {
                  /* Get the grid location */
                  ipos = top_location[ipos];

                  /* Store chain number */
                  grid[ipos] = ichain;

                  /* Set flag that this chain has been placed */
                  placed = TRUE;

                  /* Get the corresponding subscripts */
                  i = ipos / nside;
                  j = ipos % nside;

                  /* Calculate this circle's coords */
                  circle_coords[2 * ichain] = i * grid_spacing;
                  circle_coords[2 * ichain + 1] = j * grid_spacing;
                }
            }
        }

      /* Calculate the "energy" of the current configuration */
      if (abort == FALSE)
        energy = calc_energy(first_chain_ptr,nchains,radius,circle_coords,
                             include);
      else
        energy = - (float) HUGE_VAL;
        
      /* If this is the best energy so far, store the grid as the best
         so far */
      if (abort == FALSE && (have_placement == FALSE || energy > best_energy))
        {
          /* Transfer the grid */
          for (i = 0; i < ngrid; i++)
            best_grid[i] = grid[i];

          /* Store the best energy */
          best_energy = energy;

          /* Set flag that we have a reasonabloe placement */
          have_placement = TRUE;
        }

      /* Increment loop counter */
      loop++;
    }

  /* Initialise maximum and minimum coords */
  minx = miny = nside * grid_spacing;
  maxx = maxy = 0.0;

  /* Calculate the coordinates corresponding to the best grid */
  for (ipos = 0; ipos < ngrid; ipos++)
    {
      /* If have a chain at this position, compute its coords */
      if (best_grid[ipos] > -1)
        {
          /* Get the chain identifier */
          ichain = best_grid[ipos];

          /* Get the corresponding subscripts */
          i = ipos / nside;
          j = ipos % nside;

          /* Calculate this circle's coords */
          x = i * grid_spacing;
          y = j * grid_spacing;
          circle_coords[2 * ichain] = x;
          circle_coords[2 * ichain + 1] = y;

          /* Adjust maximum and minimum coords, if necessary */
          if (x - radius[ichain] < minx)
            minx = x - radius[ichain];
          if (x + radius[ichain] > maxx)
            maxx = x + radius[ichain];
          if (y - radius[ichain] < miny)
            miny = y - radius[ichain];
          if (y + radius[ichain] > maxy)
            maxy = y + radius[ichain];
        }
    }

  /* Calculate width and height of the current picture */
  width = maxx - minx;
  height = maxy - miny;

  /* Calculate mid-point of picture */
  midx = (maxx + minx) / 2;
  midy = (maxy + miny) / 2;

  /* Adjust the circle coords to place the origin in the centre of
     the picture */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Adjust x and y coords */
      circle_coords[2 * ichain] = circle_coords[2 * ichain] - midx;
      circle_coords[2 * ichain + 1]
        = circle_coords[2 * ichain + 1] - midy;
    }

  /* Return height and width of the picture */
  *ht = height;
  *wd = width;

  /* Free up memory */
  free(grid);
  free(best_grid);
  free(map);
  free(top_location);

  /* Return whether procedure returned a usable starting configuration */
  return(have_placement);
}
/***********************************************************************

compact_number  -  Compress number into shortest string for output

***********************************************************************/

int compact_number(float number,char *number_string,char *format)
{
  char temp_string[20];
  char *string_ptr;

  int ipos, len, nblanks, ppos;
  int done;

  /* Initialise variables */
  len = 0;
  for (ipos = 0; ipos < 20; ipos++)
    temp_string[ipos] = '\0';

  /* Write out the number using the given format string */
  sprintf(temp_string,format,number);
  strcpy(number_string,temp_string);

  /* Get the length of the string */
  len = strlen(number_string);
  ppos = len;

  /* If number contains a decimal point, remove any redundant zeros
     from the end */
  if ((string_ptr = strstr(number_string,".")) != NULL)
    {
      /* Get the position of the decimal point */
      ppos = string_ptr - number_string;

      /* Remove any trailing zeroes */
      done = FALSE;
      for (ipos = len - 1; ipos > ppos - 1 && done == FALSE; ipos--)
        {
          /* If this is a zero, then reduce length of string */
          if (number_string[ipos] == '0' || number_string[ipos] == '.')
            number_string[ipos] = '\0';

          /* Otherwise, set flag that we're done */
          else
            done = TRUE;
        }
    }

  /* Remove any redundant blank spaces from the front */
  done = FALSE;
  nblanks = 0;
  for (ipos = 0; ipos < ppos && done == FALSE; ipos++)
    {
      /* If this is a blank space, increment count */
      if (number_string[ipos] == ' ')
        nblanks++;

      /* Otherwise, set flag that we're done */
      else
        done = TRUE;
    }

  /* Transfer to temporary string */
  strcpy(temp_string,number_string + nblanks);

  /* Copy back */
  strcpy(number_string,temp_string);

  /* Get the string length */
  len = strlen(number_string);

  /* Return the string length */
  return(len);
}
/***********************************************************************

write_coords  -  Write out the coordinates of the chain circles to the
                 ppcont.dat file

***********************************************************************/

void write_coords(float *circle_coords,float *radius,float width,
                  float height,int nchains,float energy,
                  struct parameters params)
{
  char file_name[FILENAME_LEN], number_string[20];

  int ichain;

  FILE *file_ptr;

  /* Form name of output ppcont.dat file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(file_name,params.work_dir);
      strcat(file_name,DIR_SEP);
    }
  else
    file_name[0] = '\0';
  strcat(file_name,"ppcont.dat");

  /* Open the file */
  if ((file_ptr = fopen(file_name,"w")) !=  NULL)
    {
      /* Loop over the chains to write out the coords */
      for (ichain = 0; ichain < nchains; ichain++)
        {
          /* Write chain's x-coord */
          compact_number(circle_coords[2 * ichain],number_string,"%8.3f");
          fprintf(file_ptr,"%s",number_string);

          /* Write chain's y-coord */
          compact_number(circle_coords[2 * ichain + 1],number_string,"%8.3f");
          fprintf(file_ptr," %s",number_string);

          /* Write chain's radiu */
          compact_number(radius[ichain],number_string,"%8.3f");
          fprintf(file_ptr," %s\n",number_string);
        }

      /* Write out the width and height */
      compact_number(width,number_string,"%8.3f");
      fprintf(file_ptr,"%s",number_string);
      compact_number(height,number_string,"%8.3f");
      fprintf(file_ptr," %s\n",number_string);

      /* Write out the final energy */
      compact_number(energy,number_string,"%16.0f");
      fprintf(file_ptr,"Energy: %s\n",number_string);

      /* Close the file */
      fclose(file_ptr);
    }

  /* Show warning */
  else
    printf("*** Warning. Unable to open output file: %s\n",file_name);
}
/***********************************************************************

read_coords  -  Read in the coordinates of the chain circles from the
                ppcont.dat file

***********************************************************************/

int read_coords(float *circle_coords,float *radius,float *wd,float *ht,
                int nchains,char *pdbsum_dir,struct parameters params)
{
  char input_line[LINELEN + 1], file_name[FILENAME_LEN];

  int ichain, len, nlines;
  int have_coords;

  float energy, height, r, width, x, y;

  FILE *file_ptr;

  /* Initialise variables */
  have_coords = FALSE;
  height = 0.0;
  ichain = 0;
  width = 0.0;
  nlines = 0;
/* debug */
  printf("In read_coords ...\n");
  fflush(stdout);
/* debug */

  /* Form name of output ppcont.dat file */
  strcpy(file_name,pdbsum_dir);
  if (params.pdb_type == PDBSUM1)
    strcat(file_name,"protprot/");
  strcat(file_name,"ppcont.dat");
/* debug */
  printf("  Opening file %s\n",file_name);
  fflush(stdout);
/* debug */

  /* Open the file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Set flag */
      have_coords = TRUE;
/* debug */
      printf("   --- have file\n");
      fflush(stdout);
/* debug */

      /* Read in the data from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If this line contains chain coords, store */
          if (nlines < nchains)
            {
              /* Store the chain coords and radius */
              sscanf(input_line,"%f %f %f",&x,&y,&r);
              circle_coords[2 * ichain] = x;
              circle_coords[2 * ichain + 1] = y;
              radius[ichain] = r;

              /* Increment chain count */
              ichain++;
            }

          /* If this is the plot energy, store that */
          else if (!strncmp(input_line,"Energy:",6))
            sscanf(input_line+6,"%f",&energy);

          /* Otherwise, take to be the picture height and width */
          else
            sscanf(input_line,"%f %f",&width,&height);

          /* Increment count of lines read in */
          nlines++;
        }

      /* Close the file */
      fclose(file_ptr);
    }

  /* If no data read in from the file, unset flag */
  if (nlines < nchains || height == 0.0 || width == 0.0)
    have_coords = FALSE;

  /* Return picture height and width */
  *ht = height;
  *wd = width;

  /* Return whether we have data */
  return(have_coords);
}
/***********************************************************************

calc_force_vector  -  Calculate gradient of energy - ie force vector
                      - computed using the central difference
                      approximation

***********************************************************************/

float calc_force_vector(struct chain *first_chain_ptr,int nchains,
                        int nparams,float *radius,float *circle_coords,
                        float *force_vector,int *include)
{
  int iparam;

  float energy1, energy2, rms_force, save_coord;

  /* Initialise variables */
  rms_force = 0.0;

  /* Loop through all the coordinates, one by one */
  for (iparam = 0; iparam < nparams; iparam++)
    {
      /* Save the current coordinate */
      save_coord = circle_coords[iparam];

      /* Apply small shift to this coordinate */
      circle_coords[iparam] = save_coord + DELTA;

      /* Calculate total energy */
      energy2 = - calc_energy(first_chain_ptr,nchains,radius,circle_coords,
                              include);

      /* Shift this coordinate by DELTA in the opposite direction */
      circle_coords[iparam] = save_coord - DELTA;

      /* Calculate new energy */
      energy1 = - calc_energy(first_chain_ptr,nchains,radius,circle_coords,
                              include);

      /* Calculate force component for this parameter */
      force_vector[iparam] = - (energy2 - energy1) / (2 * DELTA);

      /* Accumulate terms for calculating rms force */
      rms_force = rms_force + force_vector[iparam] * force_vector[iparam];

      /* Restore coordinates starting value */
      circle_coords[iparam] = save_coord;
    }

  /* Calculate the rms force */
  rms_force = (float) sqrt(rms_force / (3 * nparams));

  /* Return the rms force */
  return(rms_force);
}
/***********************************************************************

shift_coords  -  Apply given shift to each of the coordinates

***********************************************************************/

void shift_coords(float *circle_coords,int nparams,float *force_vector,
                  float rms_force,float lambda)
{
  int iparam;

  float shift;

  /* Loop through all the coordinates, one by one */
  for (iparam = 0; iparam < nparams; iparam++)
    {
      /* Calculate the size of the shift */
      shift = force_vector[iparam] * lambda / rms_force;

      /* Apply small shift to this coordinate */
      circle_coords[iparam] = circle_coords[iparam] + shift;
    }
}
/***********************************************************************

energy_minimize  -  Use energy minimization to improve relative
                    dispositions of the chains

***********************************************************************/

float energy_minimize(struct chain *first_chain_ptr,int nchains,
                      float *radius,float *wd,float *ht,
                      float *circle_coords)
{
  int ichain, iseed, nparams, show_step, step, switch_on_step;
  int done, print_details;
  int include[NTERMS];

  float energy, lambda, last_energy, rms_force;
  float height, width;
  float minx, miny, maxx, maxy, midx, midy, x, y;

  float *force_vector;

  /* Initialise variables */
  done = FALSE;
  height = 0.0;
  lambda = (float) 1.0;
  energy = 0.0;
  nparams = 2 * nchains;
  print_details = TRUE;
  rms_force = 0.0;
  step = 0;
  show_step = 2 * MAXSTEPS;
  switch_on_step = 0;
  width = 0.0;

  /* Initialise energy terms to be included */
  include[CLASH_ENERGY] = FALSE;
  include[OVERLAP_ENERGY] = FALSE;
  include[SEPARATION_ENERGY] = TRUE;
  include[WEDGE_ENERGY] = FALSE;

  /* Initialise random number generator */
  iseed = -1;
  get_random_number(&iseed,RANDOM_START);

  /* Allocate memory for the force vector */
  force_vector = (float *) malloc(nparams * sizeof(float));

  /* Loop until convergence reached */
  while (done == FALSE)
    {
      /* Check whether this is the switch-on step */
      if (step == switch_on_step)
        {
          include[CLASH_ENERGY] = TRUE;
          include[OVERLAP_ENERGY] = TRUE;
          include[WEDGE_ENERGY] = TRUE;
        }

      /* Calculate total energy of system */
      energy = calc_energy(first_chain_ptr,nchains,radius,circle_coords,
                           include);

      /* If this is the first loop, then save starting energy */
      if (step == 0  || step == switch_on_step)
        last_energy = energy;

      /* Calculate gradient of energy - ie force vector */
      rms_force
        = calc_force_vector(first_chain_ptr,nchains,nparams,radius,
                            circle_coords,force_vector,include);

      /* Determine whether to stop here */
      if (rms_force < 0.3 || step >= MAXSTEPS)
        done = TRUE;

      /* Print status */
      if ((step % show_step) == 0 || done == TRUE)
        {
          printf("Step %5d. Energy = %10.3f, RMS force = %10.3f ",
                 step,energy,rms_force);
          printf("Lambda = %10.6f\n",lambda);
        }

      /* Otherwise, calculate new coordinates of whole system */
      if (done == FALSE)
        {
          /* Adjust step-size if necessary */
          if (energy > last_energy)
            lambda = (float) 1.2 * lambda;
          else if (energy < last_energy)
            lambda = (float) 0.75 * lambda;

          /* Apply all the coordinate shifts */
          shift_coords(circle_coords,nparams,force_vector,rms_force,
                       lambda);

          /* Store the current energy */
          last_energy = energy;

          /* Increment count of steps done */
          step++;
        }
    }

  /* Calculate the coordinates max and min coords of the final picture */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Get this circle's coords */
      x = circle_coords[2 * ichain];
      y = circle_coords[2 * ichain + 1];

      /* Adjust maximum and minimum coords, if necessary */
      if (ichain == 0)
        {
          minx = x - radius[ichain];
          maxx = x + radius[ichain];
          miny = y - radius[ichain];
          maxy = y + radius[ichain];
        }
      else
        {
          if (x - radius[ichain] < minx)
            minx = x - radius[ichain];
          if (x + radius[ichain] > maxx)
            maxx = x + radius[ichain];
          if (y - radius[ichain] < miny)
            miny = y - radius[ichain];
          if (y + radius[ichain] > maxy)
            maxy = y + radius[ichain];
        }
    }

  /* Calculate width and height of the current picture */
  width = maxx - minx;
  height = maxy - miny;

  /* Calculate mid-point of picture */
  midx = (maxx + minx) / 2;
  midy = (maxy + miny) / 2;

  /* Adjust the circle coords to place the origin in the centre of
     the picture */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Adjust x and y coords */
      circle_coords[2 * ichain] = circle_coords[2 * ichain] - midx;
      circle_coords[2 * ichain + 1]
        = circle_coords[2 * ichain + 1] - midy;
    }

  /* Return height and width of the picture */
  *ht = height;
  *wd = width;

  /* Return the final energy */
  return(energy);
}
/***********************************************************************

place_clusters  -  Arrange the clusters on the page

***********************************************************************/

void place_clusters(float *radius,float *circle_coords,int nchains,
                    int *chain_count,int nclusters,float *ht,float *wd)
{
  int first_chain, i, ichain, icluster;

  float height, width;
  float x, xmin, xmax, y, ymin, ymax, ybottom, yshift, ytop;

  /* Initialise variables */
  height = 0.0;
  first_chain = 0;
  ybottom = 0.0;
  ytop = 0.0;
  width = 0.0;

  /* Loop over all clusters */
  for (icluster = 0; icluster < nclusters; icluster++)
    {
      /* Initialise chain */
      ichain = first_chain;

      /* Initialise maximum and minimum values */
      xmin = circle_coords[2 * ichain] - radius[ichain];
      xmax = circle_coords[2 * ichain] + radius[ichain];
      ymin = circle_coords[2 * ichain + 1] - radius[ichain];
      ymax = circle_coords[2 * ichain + 1] + radius[ichain];
      ichain++;

      /* Loop over this cluster's other chains */
      for (i = 1; i < chain_count[icluster]; i++)
        {
          /* Adjust max and min, if required */
          x = circle_coords[2 * ichain] - radius[ichain];
          if (x < xmin)
            xmin = x;
          x = circle_coords[2 * ichain] + radius[ichain];
          if (x > xmax)
            xmax = x;

          /* Adjust max and min, if required */
          y = circle_coords[2 * ichain + 1] - radius[ichain];
          if (y < ymin)
            ymin = y;
          y = circle_coords[2 * ichain + 1] + radius[ichain];
          if (y > ymax)
            ymax = y;

          /* Increment chain counter */
          ichain++;
        }
      
      /* If this is the first cluster, then max and min y-values
         define the overall max and min values */
      if (icluster == 0)
        {
          /* Set top and bottom y-coords */
          ybottom = ymin;
          ytop = ymax;
        }

      /* Otherwise, shift whole cluster downwards */
      else
        {
          /* Determine y-shift to be applied */
          yshift = YGAP + ymax - ybottom;

          /* Loop over this cluster's chains to apply shift */
          for (i = 0; i < chain_count[icluster]; i++)
            {
              /* Get the chain identifier */
              ichain = i + first_chain;

              /* Apply shift */
              circle_coords[2 * ichain + 1]
                = circle_coords[2 * ichain + 1] - yshift;
            }

          /* Calculate the new bottom coords */
          ybottom = ymin - yshift;
        }

      /* Reset first chain */
      first_chain = first_chain + chain_count[icluster];
    }

  /* Calculate overall height and width */
  width = xmax - xmin;
  height = ytop - ybottom;
  yshift = (ytop + ybottom) / 2;

  /* Adjust all the y-coords to centre the clusters at the origin */
  for (ichain = 0; ichain < nchains; ichain++)
    circle_coords[2 * ichain + 1] = circle_coords[2 * ichain + 1] - yshift;
      
  /* Return height and width */
  *ht = height;
  *wd = width;
}
/***********************************************************************

arrange_clusters  -  Arrange the chains in each cluster to give the most
                     compact representation

***********************************************************************/

void arrange_clusters(struct chain *first_chain_ptr,int nchains,
                      struct cluster *first_cluster_ptr,int nclusters,
                      float **ccoords,float **rad,float *ht,float *wd,
                      char *pdbsum_dir,struct parameters params)
{
  int ichain, icluster;
  int have_coords, have_placement;

  int *chain_count;

  float energy, height, scale_factor, width;
  float frame_aspect, origin[2], picture_aspect;
  float global_radius, max_radius, tot_energy;

  float *circle_coords, *radius;

  struct cluster *cluster_ptr;

  /* Initialise variables */
  have_coords = FALSE;
  tot_energy = 0.0;
/* debug */
  printf("Chains %d\n",nchains);
  fflush(stdout);
/* debug */

  /* Grab memory for chain radii and circle coords */
  radius = (float *) malloc(nchains * sizeof(float));
  circle_coords = (float *) malloc(2 * nchains * sizeof(float));

  /* Initialise coords and radii */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      radius[ichain] = 0.0;
      circle_coords[2 * ichain] = 0.0;
      circle_coords[2 * ichain + 1] = 0.0;
    }

  /* Grab memory for chain counts */
  chain_count = (int *) malloc(nclusters * sizeof(int));

  /* Initialise counts */
  for (icluster = 0; icluster < nclusters; icluster++)
    chain_count[icluster] = 0;

  /* Initialise chain */
  ichain = 0;

  /* Compute optimal interaction wedges */
  compute_interaction_wedges(first_chain_ptr);

  /* If circle coords have already been computed, read in from ppcont.dat
     file */
  if (params.read == TRUE)
    have_coords
      = read_coords(circle_coords,radius,&width,&height,nchains,
                    pdbsum_dir,params);

  /* If have no circle coords, then need to compute and minimize them */
  if (have_coords == FALSE)
    {
      /* Get first of the clusters */
      cluster_ptr = first_cluster_ptr;

      /* Loop over all the clusters to plot */
      while (cluster_ptr != NULL)
        {
          /* Store the count of chains in this cluster */
          icluster = cluster_ptr->number - 1;
          chain_count[icluster] = cluster_ptr->nchains;

          /* Place the chains to give as clear a picture of the inter-chain
             interactions as possible */
          if (cluster_ptr->nchains > 2)
            have_placement
              = place_chains(cluster_ptr->first_chain_ptr,
                             cluster_ptr->nchains,radius + ichain,
                             &width,&height,
                             circle_coords + 2 * ichain);
          else
            have_placement = FALSE;

          /* If failed to place the chains on the grid, then place on
             a circle instead */
          if (have_placement == FALSE && cluster_ptr->nchains > 1)
            {
              /* Calculate overall radius of circle on which all the
                 chain-circles are to be placed */
              global_radius
                = calc_global_radius(cluster_ptr->first_chain_ptr,
                                     cluster_ptr->nchains,
                                     radius + ichain,&max_radius,
                                     circle_coords + 2 * ichain);

              /* Determine maximum extent of the diagram */
              width = 2 * (global_radius + max_radius);
              height = 2 * (global_radius + max_radius);
            }

          /* Energy-minimize current picture to improve appearance */
          energy
            = energy_minimize(cluster_ptr->first_chain_ptr,
                              cluster_ptr->nchains,radius + ichain,&width,
                              &height,circle_coords + 2 * ichain);
          tot_energy = tot_energy + energy;

          /* Increment chain */
          ichain = ichain + cluster_ptr->nchains;

          /* Get the next cluster */
          cluster_ptr = cluster_ptr->next_cluster_ptr;
        }

      /* If have more than one cluster, need to arrange them on a
         single page */
      if (nclusters > 1)
        {
          /* Place the clusters on the page */
          place_clusters(radius,circle_coords,nchains,chain_count,
                         nclusters,&height,&width);
        }

      /* If current chain coords to be saved, write to ppcont.dat file */
      if (params.save == TRUE)
        write_coords(circle_coords,radius,width,height,nchains,tot_energy,
                     params);
    }

  /* Return the chain radii and coords */
  *ccoords = circle_coords;
  *rad = radius;

  /* Return height and width of overall picture */
  *ht = height;
  *wd = width;
}
/***********************************************************************

plot_all_interactions  -  Plot all the interactions between the chains

***********************************************************************/

void plot_all_interactions(struct chain *first_chain_ptr,int nchains,
                           float scale_factor,float *origin,
                           float *circle_coords,struct pspage *pspage_ptr)
{
  char chain[2];

  static char *line_colour[] = { "sky_blue", "orange", "red", "yellow" };

  int i, ichain, jchain, nlines;
  int flip;

  float line_sep, line_width, theta;
  float dx, dy, x, xleft, xright, x1, x2, y, yleft, yright, y1, y2;
  float xfirst, xsecond, yfirst, ysecond;

  struct chain *chain_ptr;
  struct cluster *cluster_ptr;
  struct interface *interface_ptr;

  /* Get pointer to the first `chain */
  chain_ptr = first_chain_ptr;

  /* Define colours and line widths */
  line_width = scale_factor * LINE_WIDTH;
  if (line_width < 1.0)
    line_width = (float) 1.0;
  pscomm_(pspage_ptr->ps_file,"All interactions");
  pslwid_(pspage_ptr->ps_file,line_width);

  /* Determine the separation between interaction lines joining any two
     circles */
  line_sep = LINE_SEP * scale_factor;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get this chain */
      chain[0] = chain_ptr->chain;
      chain[1] = '\0';
      ichain = chain_ptr->chain_no;

      /* Get this chains's coords */
      xfirst = origin[0] + scale_factor * circle_coords[2 * ichain];
      yfirst = origin[1] + scale_factor * circle_coords[2 * ichain + 1];

      /* Get the first of this chain's interactions */
      interface_ptr = chain_ptr->first_interface_ptr;

      /* Loop over the interfaces */
      while (interface_ptr != NULL)
        {
          /* Get the other chain involved in the interface */
          jchain = interface_ptr->chain_ptr->chain_no;

          /* Get this chains's coords */
          xsecond = origin[0] + scale_factor * circle_coords[2 * jchain];
          ysecond = origin[1] + scale_factor * circle_coords[2 * jchain + 1];

          /* Determine whyich point is on the left and which on the right */
          if (xfirst < xsecond)
            {
              xleft = xfirst;
              yleft = yfirst;
              xright = xsecond;
              yright = ysecond;
            }
          else
            {
              xleft = xsecond;
              yleft = ysecond;
              xright = xfirst;
              yright = yfirst;
            }

          /* Calculate angle made by the line */
          x = xright - xleft;
          y = yright - yleft;
          theta = get_angle(x,y);

          /* Store this angle */
          if (xfirst < xsecond)
            interface_ptr->angle = theta;
          else
            interface_ptr->angle = theta - PI;

          /* Count the number of interaction lines we need to draw */
          nlines = 0;
          for (i = 0; i < NTYPES; i++)
            if (interface_ptr->count[i] > 0)
              nlines++;

          /* Initialise the flip factor (zero if an odd number of lines and
             1 if an even number) */
          flip = 0;
          if (nlines % 2 == 0)
            flip = 1;

          /* Loop over the lines showing interaction counts */
          for (i = 0; i < NTYPES; i++)
            {
              if (interface_ptr->count[i] > 0)
                {
                  /* Define this interaction's line colour */
                  pscolb_(pspage_ptr->ps_file,line_colour[i]);

                  /* Calculate coords of this line */
                  dx = flip * line_sep * (float) sin(theta) / 2;
                  dy = flip * line_sep * (float) cos(theta) / 2;
                  if (theta > PI)
                    dx = -dx;
                  else
                    dy = -dy;
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

          /* Get the next interface */
          interface_ptr = interface_ptr->next_interface_ptr;
        }

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

plot_all_circles  -  Plot all chains as circles in a ring

***********************************************************************/

void plot_all_circles(struct chain *first_chain_ptr,int nchains,
                      float scale_factor,float *origin,float *radius,
                      float *circle_coords,struct pspage *pspage_ptr)
{
  char chain[2], colour[COLNAMLEN + 1];

  int ichain;

  float line_width, r, x, y;

  struct chain *chain_ptr;

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Define colours and line widths */
  line_width = scale_factor * CIRCLE_BORDER_WIDTH;
  if (line_width < 1.0)
    line_width = (float) 1.0;
  pscomm_(pspage_ptr->ps_file,"All circles");
  pslwid_(pspage_ptr->ps_file,line_width);

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get this chain */
      chain[0] = chain_ptr->chain;
      chain[1] = '\0';
      ichain = chain_ptr->chain_no;

      /* Set line colour for circle border */
      pscolb_(pspage_ptr->ps_file,"black");

      /* Calculate coords of centre of circle */
      x = origin[0] + scale_factor * circle_coords[2 * ichain];
      y = origin[1] + scale_factor * circle_coords[2 * ichain + 1];

      /* Calculate radius of this circle */
      r = scale_factor * radius[ichain];

      /* Get this chain's colour */
      get_ppchain_colour(chain_ptr->chain,colour);

      /* Define circle colour */
      psccol_(pspage_ptr->ps_file,colour);

      /* Plot the circle */
      pscirc_(pspage_ptr->ps_file,x,y,r);

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;
    }
}
/***********************************************************************

plot_wedges  -  Plot the interaction wedges on each circle

***********************************************************************/

void plot_wedges(struct chain *first_chain_ptr,float scale_factor,
                 float *origin,float *radius,float *circle_coords,
                 struct pspage *pspage_ptr)
{
  char bgchain[2], chain[2], colour[COLNAMLEN + 1], bgcolour[COLNAMLEN + 1];

  int ichain;

  float dx, line_width, r, text_size, x, y;

  struct chain *chain_ptr;
  struct interface *interface_ptr;
  struct wedge *wedge_ptr;

  /* Initialise variables */
  pscomm_(pspage_ptr->ps_file,"Interface wedges");

  /* Define colours and line widths */
  line_width = scale_factor * WEDGE_LINE_WIDTH;
  if (line_width < 1.0)
    line_width = (float) 1.0;
  pslwid_(pspage_ptr->ps_file,line_width);

  /* Get pointer to the first chain */
  chain_ptr = first_chain_ptr;

  /* Loop over all the chains */
  while (chain_ptr != NULL)
    {
      /* Get this chain's colour */
      get_ppchain_colour(chain_ptr->chain,bgcolour);

      /* Get this chain identifier */
      bgchain[0] = chain_ptr->chain;
      bgchain[1] = '\0';
      ichain = chain_ptr->chain_no;

      /* Calculate coords of centre of circle */
      x = origin[0] + scale_factor * circle_coords[2 * ichain];
      y = origin[1] + scale_factor * circle_coords[2 * ichain + 1];

      /* Calculate radius of this circle */
      r = scale_factor * radius[ichain];

      /* Define wedge border colour */
      pscolb_(pspage_ptr->ps_file,"black");

      /* Get pointer to first of this circle's int5erface wedges */
      wedge_ptr = chain_ptr->first_wedge_ptr;

      /* Loop over the wedges to plot */
      while (wedge_ptr != NULL)
        {
          /* Get the interface that this wedge corresponds to */
          interface_ptr = wedge_ptr->interface_ptr;

          /* Get interaction angle */
          if (wedge_ptr->first_chain == TRUE)
            wedge_ptr->angle[CENTRE] = RADDEG * interface_ptr->angle;
          else
            wedge_ptr->angle[CENTRE] = RADDEG * (interface_ptr->angle - PI);

          /* Calculate start and end angles for this wedge */
          wedge_ptr->angle[START]
            = wedge_ptr->angle[START] + wedge_ptr->angle[CENTRE];
          wedge_ptr->angle[END]
            = wedge_ptr->angle[END] + wedge_ptr->angle[CENTRE];

          /* Get the chain with which current chain interacts */
          chain[0] = wedge_ptr->other_chain_ptr->chain;
          chain[1] = '\0';

          /* Get this chain's colour */
          get_ppchain_colour(chain[0],colour);

          /* Define circle colour */
          psccol_(pspage_ptr->ps_file,colour);

          /* Plot the wedge representing the contact area */
          pspcir_(pspage_ptr->ps_file,x,y,r,wedge_ptr->angle[START],
                  wedge_ptr->angle[END]);

          /* Plot inner line of wedge */
          pspcir_(pspage_ptr->ps_file,x,y,r * INNER_FRACTION + line_width,
                  wedge_ptr->angle[START],wedge_ptr->angle[END]);

          /* Get the next wedge */
          wedge_ptr = wedge_ptr->next_wedge_ptr;
        }

      /* Redefine circle colour */
      pscolb_(pspage_ptr->ps_file,bgcolour);
      psccol_(pspage_ptr->ps_file,bgcolour);

      /* Plot the inner circle to overwite wedge */
      pscirc_(pspage_ptr->ps_file,x,y,r * INNER_FRACTION);

      /* Plot the chain identifier */
      text_size = scale_factor * CHAINID1_HEIGHT;
      if (text_size > r)
        text_size = r;
      dx = text_size / (float) 20.0;
      pscolb_(pspage_ptr->ps_file,"black");
      psctxt_(pspage_ptr->ps_file,x + dx,y + dx,text_size,bgchain);

      /* Plot slightly shifted text */
      pscolb_(pspage_ptr->ps_file,"grey");
      psctxt_(pspage_ptr->ps_file,x,y,text_size,bgchain);

      /* Get the next chain */
      chain_ptr = chain_ptr->next_chain_ptr;

      /* Increment chain counter */
      ichain++;
    }
}
/***********************************************************************

plot_allchains_diagram  -  Plot the diagram of all chains in this cluster
                           and their mutual interactions

***********************************************************************/

void plot_allchains_diagram(struct cluster *first_cluster_ptr,
                            struct chain *first_chain_ptr,int nchains,
                            float *circle_coords,float *radius,
                            float height,float width,
                            struct psparams psparams_ptr,
                            struct pspage *pspage_ptr,
                            struct parameters params)
{
  int have_coords, have_placement;

  float default_height, energy, scale_factor;
  float frame_aspect, origin[2], picture_aspect;
  float global_radius, max_radius;

  /* Calculate default height */
  default_height = PICTURE_HEIGHT * DEFAULT_WIDTH / PICTURE_WIDTH;

  /* Calculate aspect ratio for picture and for frame within picture to
     be placed */
  picture_aspect = height / width;
  frame_aspect = PICTURE_HEIGHT / PICTURE_WIDTH;

  /* Determine optimal scaling factor for the picture */
  if (picture_aspect > frame_aspect)
    {
      /* Bounded by y: calculate the scaling factor */
      if (height < default_height)
        height = default_height;
      scale_factor = PICTURE_HEIGHT / height;

      /* Determine origin of central point */
      origin[0] = PICTURE_WIDTH / 2;
      origin[1] = scale_factor * height / 2;
    }
  else
    {
      /* Bounded by x: calculate the scaling factor */
      if (width < DEFAULT_WIDTH)
        width = DEFAULT_WIDTH;
      scale_factor = PICTURE_WIDTH / width;

      /* Determine origin of central point */
      origin[0] = scale_factor * width / 2;
      origin[1] = PICTURE_HEIGHT / 2;
    }

  /* Determine the picture origin */
  origin[0] = origin[0] + XMARGIN;
  origin[1] = origin[1] + YMARGIN;

  /* Form name of output PostScript file */
  if (params.work_dir[0] != '\0')
    {
      strcpy(pspage_ptr->ps_name,params.work_dir);
      strcat(pspage_ptr->ps_name,DIR_SEP);
    }
  strcat(pspage_ptr->ps_name,"ppcont.ps");

  /* Initialise PostScript page variables */
  init_pspage(psparams_ptr,pspage_ptr);

  /* Start the page holding the plot */
  start_new_page(psparams_ptr,pspage_ptr);

  /* Plot all the interactions between the circles */
  plot_all_interactions(first_chain_ptr,nchains,scale_factor,origin,
                        circle_coords,pspage_ptr);
/* debug */
  printf("scale_factor =  %8.3f\n",scale_factor);
  fflush(stdout);
  printf("radii =  %8.3f  and %8.3f\n",radius[0],radius[1]);
  fflush(stdout);
/* debug */

  /* Draw all the circles */
  plot_all_circles(first_chain_ptr,nchains,scale_factor,origin,radius,
                   circle_coords,pspage_ptr);

  /* Plot the interaction wedges on each circle */
  plot_wedges(first_chain_ptr,scale_factor,origin,radius,circle_coords,
              pspage_ptr);

  /* Close the current page */
  close_page(psparams_ptr,pspage_ptr);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  char pdbsum_dir[FILENAME_LEN];

  int icluster, nchains, nclusters;
  int dir_type;

  float height, width;

  float *circle_coords, *radius;

  struct chain *first_chain_ptr;
  struct cluster *cluster_ptr, *first_cluster_ptr;
  struct c_file_data file_name_ptr;
  struct interface interface_ptr;
  struct parameters params;
  struct psparams psparams_ptr;
  struct pspage pspage_ptr;

  /* Read in the command-line arguments */
  get_command_arguments(argv,argc - 1,&params);

  /* Get required paths and file names */
  c_get_paths(&file_name_ptr,params.run_64);
/* debug */
  printf("TYPE = %d\n",params.pdb_type);
  fflush(stdout);
/* debug */

  /* Form the name of the PDBsum directory */
  dir_type = find_pdbsum_dir(params.pdb_code,file_name_ptr.pdbsum_template,
                             params.pdb_type,pdbsum_dir);
  if (dir_type == NOT_FOUND)
    {
      printf("*** ERROR. Unable to locate PDBsum directory for %s\n",
             params.pdb_code);
      exit(-1);
    }

  /* Initialise PostScript parameters */
  initialise_psparams(&psparams_ptr,&pspage_ptr);

  /* If plotting a single interface between two chains, plot that */
  if (params.all_chains == FALSE)
    {
      /* Read in the data for this interface from the interfaces.dat file */
      read_database_file(pdbsum_dir,"interfaces.dat",&interface_ptr,
                         params);

      /* Plot the diagram */
      plot_diagram(interface_ptr,psparams_ptr,&pspage_ptr,params);
    }

  /* Otherwise, read in data for all chains and produce overall plot */
  else
    {
      /* Read in all the chains and their interface data */
      nchains = read_chains(pdbsum_dir,"interfaces.dat",&first_chain_ptr,
                            params);

      /* Plot the "all-chains" diagram */
      if (nchains > 1)
        {
          /* Determine number of interacting chain clusters there are */
          nclusters = count_clusters(first_chain_ptr,nchains,
                                     &first_cluster_ptr);

          /* Arrange the chains in each cluster */
          arrange_clusters(first_chain_ptr,nchains,first_cluster_ptr,
                           nclusters,&circle_coords,&radius,&height,
                           &width,pdbsum_dir,params);

          /* Plot each cluster diagram */
          plot_allchains_diagram(first_cluster_ptr,first_chain_ptr,
                                 nchains,circle_coords,radius,height,
                                 width,psparams_ptr,&pspage_ptr,
                                 params);
        }
    }

  return(0);
}
