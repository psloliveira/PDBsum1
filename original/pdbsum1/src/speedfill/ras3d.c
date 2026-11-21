/* ras3d.c  -  Commonly used Raster3d routines */

/**************   R A S T E R 3 D   R O U T I N E S   *************/

/* Routines:

   -> calculate_rascale
   -> generate_polygons
         -> binary_edge_search
         -> get_grid_value
               -> extract_grid_value
         -> make_line
         -> make_polygons
               -> generate_triangles
                     -> create_polygon_record
                           -> create_polist_record
                     -> check_vertex_order
                           -> get_vertex_coords
                     -> get_vertex_coords
                     -> update_normals
   -> order_polygon_vertices
         -> orien3
               -> vecprd
               -> dot3
         -> calculate_normal
   -> calculate_vertex_normals
   -> write_raster3d_polygons
         -> get_colour_cpk
               -> get_colour_rgb (in common.c)
         -> get_colour_restype
               -> get_colour_rgb (in common.c)
         -> get_colour_conserv
               -> get_colour_rgb (in common.c)
         -> get_chain_colour_numbers
         -> get_colour_rgb (in common.c)

*/


#include "common.h"
#include "ras3d.h"

/* Prototypes
   ---------- */

/* In common.c */
void get_colour_rgb(char *c_name,double rgb[3]);


/***********************************************************************

calculate_rascale  -  Calculate the scaling factor for Raster3D output

***********************************************************************/

double calculate_rascale(struct limits vertex_limits,double border)
{
  double maxwidth, rascale, xwidth, ywidth;

  /* Initialise variables */
  rascale =  0.0;
  maxwidth =  0.0;

  /* Calculate maximum width in x or y direction */
  xwidth = vertex_limits.valmax[0] - vertex_limits.valmin[0];
  ywidth = vertex_limits.valmax[1] - vertex_limits.valmin[1];
  if (xwidth > ywidth)
    maxwidth = xwidth;
  else
    maxwidth = ywidth;

  /* Calculate scaling factor */
  rascale =  ((1.0 - 2.0 * border) * 1.0 / maxwidth);

  /* Return the scaling factor */
  return(rascale);
}
/***********************************************************************
  
binary_edge_search  -  Find vertex corresponding to given edge number,
                       if there is one
  
***********************************************************************/

struct vertex *binary_edge_search(int edge_no,
                                  struct vertex **vertex_index_ptr,
                                  int nvertices)
{
  int ivertex, lower, upper;
  int done;

  struct vertex *vertex_ptr, *match_vertex_ptr;

  /* Initialise variables */
  done = FALSE;
  match_vertex_ptr = NULL;
  lower = 0;
  upper = nvertices - 1;

  /* Loop until search item located, or no more searching to be done */
  while (match_vertex_ptr == NULL && done == FALSE)
    {
      /* Get position of next vertex to be tried */
      ivertex = (lower + upper) / 2;

      /* Get this vertex pointer */
      vertex_ptr = vertex_index_ptr[ivertex];

      /* If edge number matches the search edge, then have found match */
      if (vertex_ptr->edge_no == edge_no)
        {
          /* Store the vertex pointer */
          match_vertex_ptr = vertex_ptr;
        }

      /* If this vertex's edge number is too high, then reduce the
         upper boundary */
      else if (vertex_ptr->edge_no > edge_no)
        {
          /* If upper boundary can't be reduced, then we are done */
          if (lower == upper)
            done = TRUE;

          /* Otherwise, set it to the current vertex position */
          else
            upper = ivertex;
        }

      /* Otherwise, increase the lower boundary */
      else
        {
          /* If lower boundary can't be increased, then we are done */
          if (lower == upper)
            done = TRUE;

          /* If we are already at the lowest position, then set
             lower boundary same as upper one */
          else if (lower == ivertex)
            lower = upper;

          /* Otherwise, set lower boundary to current vertex position */
          else
            lower = ivertex;
        }
    }

  /* Return the matching vertex pointer, if there is one */
  return(match_vertex_ptr);
}
/***********************************************************************
  
create_polist_record  -  Create and initialise a new polist record
  
***********************************************************************/

struct polist *create_polist_record(struct polist **fst_polist_ptr,
                                    struct polist **lst_polist_ptr,
                                    struct polygon *polygon_ptr)
{
  struct polist *polist_ptr, *first_polist_ptr, *last_polist_ptr;

  /* Initialise polist pointers */
  polist_ptr = NULL;
  first_polist_ptr = *fst_polist_ptr;
  last_polist_ptr = *lst_polist_ptr;

  /* Allocate memory for structure to hold polist info */
  polist_ptr = (struct polist *) malloc(sizeof(struct polist));
  if (polist_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct polist\n");
      exit (1);
    }

  /* If this is the very first polist, then store pointer */
  if (first_polist_ptr == NULL)
    first_polist_ptr = polist_ptr;

  /* Add link from previous polist to the current one */
  if (last_polist_ptr != NULL)
    last_polist_ptr->next_polist_ptr = polist_ptr;
  last_polist_ptr = polist_ptr;

  /* Store pointer to polygon */
  polist_ptr->polygon_ptr = polygon_ptr;

  /* Initialise the other elements of the structure */
  polist_ptr->next_polist_ptr = NULL;

  /* Return current pointers */
  *fst_polist_ptr = first_polist_ptr;
  *lst_polist_ptr = last_polist_ptr;

  /* Return new polist pointer */
  return(polist_ptr);
}
/***********************************************************************
  
create_polygon_record  -  Create and initialise a new polygon record
  
***********************************************************************/

struct polygon *create_polygon_record(struct polygon **fst_polygon_ptr,
                                      struct polygon **lst_polygon_ptr,
                                      int polygon_no,
                                      struct vertex *vertex1_ptr,
                                      struct vertex *vertex2_ptr,
                                      struct vertex *vertex3_ptr,
                                      int low_point)
{
  struct polygon *polygon_ptr, *first_polygon_ptr, *last_polygon_ptr;

  /* Initialise polygon pointers */
  polygon_ptr = NULL;
  first_polygon_ptr = *fst_polygon_ptr;
  last_polygon_ptr = *lst_polygon_ptr;

  /* Allocate memory for structure to hold polygon info */
  polygon_ptr = (struct polygon *) malloc(sizeof(struct polygon));
  if (polygon_ptr == NULL)
    {
      printf("*** Can't allocate memory for struct polygon\n");
      exit (1);
    }

  /* If this is the very first polygon, then store pointer */
  if (first_polygon_ptr == NULL)
    first_polygon_ptr = polygon_ptr;

  /* Add link from previous polygon to the current one */
  if (last_polygon_ptr != NULL)
    last_polygon_ptr->next_polygon_ptr = polygon_ptr;
  last_polygon_ptr = polygon_ptr;

  /* Store the 3 vertices */
  polygon_ptr->vertex_ptr[0] = vertex1_ptr;
  polygon_ptr->vertex_ptr[1] = vertex2_ptr;
  polygon_ptr->vertex_ptr[2] = vertex3_ptr;

  /* Add the current polygon to each of these 3 vertices' polygon
     lists */
  create_polist_record(&(vertex1_ptr->first_polist_ptr),
                       &(vertex1_ptr->last_polist_ptr),polygon_ptr);
  create_polist_record(&(vertex2_ptr->first_polist_ptr),
                       &(vertex2_ptr->last_polist_ptr),polygon_ptr);
  create_polist_record(&(vertex3_ptr->first_polist_ptr),
                       &(vertex3_ptr->last_polist_ptr),polygon_ptr);

  /* Initialise the other elements of the structure */
  polygon_ptr->polygon_no = polygon_no;
  polygon_ptr->polygon_colour = 0;
  polygon_ptr->outside_point = low_point;
  polygon_ptr->next_polygon_ptr = NULL;

  /* Return current pointers */
  *fst_polygon_ptr = first_polygon_ptr;
  *lst_polygon_ptr = last_polygon_ptr;

  /* Return new polygon pointer */
  return(polygon_ptr);
}
/***********************************************************************
  
make_line  -  Determine which edges are joined by lines linking
              contour cross-over points
  
***********************************************************************/

int make_line(int ncross_face,int nlines,int cell_line[12][2],
              struct vertex *line_vertex_ptr[12][2],
              struct vertex *side_vertex_ptr[4],int on_edge[4],
              int last_value,int contour_level)
{
  /* If only two cross-points found, then join them */
  if (ncross_face == 2)
    {
      /* Store the two vertices at either end of the line */
      line_vertex_ptr[nlines][0] = side_vertex_ptr[0];
      line_vertex_ptr[nlines][1] = side_vertex_ptr[1];

      /* Store which of the face's edges are involved */
      cell_line[nlines][0] = on_edge[0];
      cell_line[nlines][1] = on_edge[1];

      /* Increment count of lines created */
      nlines++;
    }

  /* If all four edges crossed, then order of joining depends on 
     whether last point encountered was above or below the contour
     level */
  else if (ncross_face == 4)
    {
      /* Check face's last value for order */
      if (last_value > contour_level)
        {
          /* Store the line's two vertices and edges */
          line_vertex_ptr[nlines][0] = side_vertex_ptr[0];
          line_vertex_ptr[nlines][1] = side_vertex_ptr[3];
          cell_line[nlines][0] = on_edge[0];
          cell_line[nlines][1] = on_edge[3];

          /* Increment count of lines created */
          nlines++;

          /* Repeat for other line */
          line_vertex_ptr[nlines][0] = side_vertex_ptr[1];
          line_vertex_ptr[nlines][1] = side_vertex_ptr[2];
          cell_line[nlines][0] = on_edge[1];
          cell_line[nlines][1] = on_edge[2];

          /* Increment count of lines created */
          nlines++;
        }
      else
        {
          /* Store the line's two vertices and edges */
          line_vertex_ptr[nlines][0] = side_vertex_ptr[0];
          line_vertex_ptr[nlines][1] = side_vertex_ptr[1];
          cell_line[nlines][0] = on_edge[0];
          cell_line[nlines][1] = on_edge[1];

          /* Increment count of lines created */
          nlines++;

          /* Repeat for other line */
          line_vertex_ptr[nlines][0] = side_vertex_ptr[2];
          line_vertex_ptr[nlines][1] = side_vertex_ptr[3];
          cell_line[nlines][0] = on_edge[2];
          cell_line[nlines][1] = on_edge[3];

          /* Increment count of lines created */
          nlines++;
        }
    }

  /* Return the current number of lines */
  return(nlines);
}
/***********************************************************************

get_vertex_coords  -  Get this polygon's 9 vertex coords

***********************************************************************/

void get_vertex_coords(struct polygon *polygon_ptr,double *v)
{
  int icoord, ipos, ivert;

  struct vertex *vertex_ptr;

  /* Initialise */
  ipos = 0;

  /* Loop over the polygon's three vertices */
  for (ivert = 0; ivert < 3; ivert++)
    {
      /* Get pointer to this vertex */
      vertex_ptr = polygon_ptr->vertex_ptr[ivert];

      /* Write out this vertex's coordinates */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Store the coordinate */
          v[ipos] = vertex_ptr->coords[icoord];
          ipos++;
        }
    }
}
/***********************************************************************

check_vertex_order  -  Check whether polygon vertices are in the correct
                       order

***********************************************************************/

int check_vertex_order(struct polygon *polygon_ptr,int npoints[3],
                       double valmin[3],double valmax[3],
                       double grid_sep,int low_point,double *v)
{
  int icoord, iorder;
  int calc, i, j, k;
  int swap;
  
  double vswap, x, y, z;
  
  struct vertex *vertex_ptr;
  
  /* Initialise variables */
  swap = FALSE;

  /* Convert low-point into array location */
  i = low_point / (npoints[1] * npoints[2]);
  calc = low_point - i * npoints[1] * npoints[2];
  j = calc / npoints[2];
  k = calc - j * npoints[2];

  /* Calculate the point's x-, y- and z-coordinates */
  x = valmin[0] + i * grid_sep;
  y = valmin[1] + j * grid_sep;
  z = valmin[2] + k * grid_sep;

  /* Get the vertex coords */
  get_vertex_coords(polygon_ptr,v);

  /* Check the vertex order */
  iorder = orien3(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],
                  x,y,z);

  /* If the order is anti-clockwise, then reverse order of vertices
     to make order clockwise */
  if (iorder == 1)
    {
      /* Set swap flag */
      swap = TRUE;

      /* Swap the first and third vertices */
      vertex_ptr = polygon_ptr->vertex_ptr[0];
      polygon_ptr->vertex_ptr[0] = polygon_ptr->vertex_ptr[2];
      polygon_ptr->vertex_ptr[2] = vertex_ptr;

      /* Swap the coordinates as well */
      for (icoord = 0; icoord < 3; icoord++)
        {
          vswap = v[icoord];
          v[icoord] = v[icoord + 6];
          v[icoord + 6] = vswap;
        }
    }

  /* Return flag indicating whether order has been swapped */
  return(swap);
}
/***********************************************************************

calculate_vertex_normals  -  Calculate vertex normals as an average of
                             the polygon normals or which the vertex is
                             a member

***********************************************************************/

void calculate_vertex_normals(struct vertex *first_vertex_ptr,
                              struct limits *vertex_limits)
{
  int icoord, npts;

  struct vertex *vertex_ptr;

  /* Initialise variables */
  for (icoord = 0; icoord < 3; icoord++)
    {
      vertex_limits->valmin[icoord] = 0.0;
      vertex_limits->valmax[icoord] = 0.0;
      vertex_limits->c_of_g[icoord] = 0.0;
    }
  npts = 0;
  vertex_limits->first = TRUE;

  /* Get first vertex pointer */
  vertex_ptr = first_vertex_ptr;

  /* Loop through all the vertices */
  while (vertex_ptr != NULL)
    {
      /* If vertex belongs to any of the polygons, then calculate an
         average normal */
      if (vertex_ptr->npolygons > 0)
        {
          /* Calculate averaged normal */
          vertex_ptr->normal[0]
            = vertex_ptr->normal[0] / vertex_ptr->npolygons;
          vertex_ptr->normal[1]
            = vertex_ptr->normal[1] / vertex_ptr->npolygons;
          vertex_ptr->normal[2]
            = vertex_ptr->normal[2] / vertex_ptr->npolygons;
        }

      /* If first vertex, then store its coords as maximum and minimum */
      if (vertex_limits->first == TRUE)
        {
          for (icoord = 0; icoord < 3; icoord++)
            {
              vertex_limits->valmin[icoord] = vertex_ptr->coords[icoord];
              vertex_limits->valmax[icoord] = vertex_ptr->coords[icoord];
            }

          /* Reset the flag */
          vertex_limits->first = FALSE;
        }

      /* Otherwise, see if limits need updating */
      else
        {
          /* Update the cordinate limits */
          for (icoord = 0; icoord < 3; icoord++)
            {
              if (vertex_ptr->coords[icoord]
                  < vertex_limits->valmin[icoord])
                vertex_limits->valmin[icoord]
                  = vertex_ptr->coords[icoord];
              if (vertex_ptr->coords[icoord]
                  > vertex_limits->valmax[icoord])
                vertex_limits->valmax[icoord]
                  = vertex_ptr->coords[icoord];
            }
        }

      /* Update centre of gravity accumulators */
      npts++;
      for (icoord = 0; icoord < 3; icoord++)
        {
          vertex_limits->c_of_g[icoord]
            = vertex_limits->c_of_g[icoord] + vertex_ptr->coords[icoord];
        }

      /* Get next vertex record */
      vertex_ptr = vertex_ptr->next_vertex_ptr;
    }

  /* Calculate overall centre of gravity */
  if (npts > 0)
    {
      for (icoord = 0; icoord < 3; icoord++)
        {
          vertex_limits->c_of_g[icoord]
            = vertex_limits->c_of_g[icoord] / npts;
        }
    }
}
/***********************************************************************
  
update_normals  -  For polygons with more than three vertices, break
                       up into a set of individual triangles
  
***********************************************************************/

void update_normals(struct polygon *polygon_ptr,double *v)
{
  int icoord, ivert;

  double normal[3];

  struct vertex *vertex_ptr;

  /* Calculate vector normal to polygon */
  calculate_normal(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],
                   normal);

  /* Add normals to each vertex's accumulator for deriving an
     average surface normal at that vector */
  for (ivert = 0; ivert < 3; ivert++)
    {
      /* Get pointer to this vertex */
      vertex_ptr = polygon_ptr->vertex_ptr[ivert];

      /* Update the normal */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Store the coordinate */
          vertex_ptr->normal[icoord]
            = vertex_ptr->normal[icoord] + normal[icoord];
        }

      /* Increment count of polygons to which this vertex belongs */
      vertex_ptr->npolygons++;
    }
}
/***********************************************************************
  
generate_triangles  -  For polygons with more than three vertices, break
                       up into a set of individual triangles
  
***********************************************************************/

int generate_triangles(struct polygon **fst_polygon_ptr,
                       struct polygon **lst_polygon_ptr,int npol,
                       struct vertex *polygon_vertex_ptr[12],int nvert,
                       int store_low_point[12],double corner_coords[3],
                       int corner_value[8],double grid_sep,
                       int contour_level,int npoints[3],double valmin[3],
                       double valmax[3])
{
  int icoord, iloop, ipoint, itriang, ivert, jvert, kvert;
  int bloop, navpoint, nloop, npolygons;
  int swap;

  double ave_score, min_score, score;
  double distance, mid_point[3], v[9], x, y, z;

  struct polygon *polygon_ptr;
  struct polygon *first_polygon_ptr, *last_polygon_ptr;
  struct vertex *vertex_ptr;

  static int cell_points[8][3]
      = { 0, 0, 0,
          0, 0, 1,
          0, 1, 0,
          0, 1, 1,
          1, 0, 0,
          1, 0, 1,
          1, 1, 0,
          1, 1, 1 };

  /* Initialise variables */
  npolygons = npol;

  /* Initialise pointers */
  first_polygon_ptr = *fst_polygon_ptr;
  last_polygon_ptr = *lst_polygon_ptr;

  /* If polygon has three vertices, then have a triangle, so can create
     a polygon record */
  if (nvert == 3)
    {
      /* Create a polygon record */
      polygon_ptr
        = create_polygon_record(&first_polygon_ptr,&last_polygon_ptr,
                                npolygons,polygon_vertex_ptr[0],
                                polygon_vertex_ptr[1],
                                polygon_vertex_ptr[2],
                                store_low_point[0]);

      /* Check the vertex order */
      swap = check_vertex_order(polygon_ptr,npoints,valmin,valmax,
                                grid_sep,store_low_point[0],v);

      /* Calculate normal to polygon and update vertex normal
         accumulators */
      update_normals(polygon_ptr,v);

      /* Increment count of polygons */
      npolygons++;
    }

  /* Otherwise, chop up the polygon into triangles so as to give a
     surface that most closely matches the surface being contoured */
  else if (nvert > 3)
    {
      /* Determine how many trial arrangements of triangles we need
         to make */
      nloop = nvert;

      /* If we have 4 vertices, then only need to make two trials
         as trials 3 and 4 are equivalent to trials 1 and 2
           
                            x            x
                          /   \        / | \
                        x ----- x    x   |   x
                         \    /        \ | /
                           x             x                    */
      if (nvert == 4)
        nloop = 2;

      /* Initialse minimum score */
      min_score =  0.0;

      /* Loop through all the trials */
      for (iloop = 0; iloop < nloop; iloop++)
        {
          /* Initialise score */
          score =  0.0;

          /* Loop through all possible triangles in this arrangement */
          for (itriang = 0; itriang < nvert - 2; itriang++)
            {
              /* Determine the vertices of the current triangle */
              ivert = iloop;
              jvert = ivert + itriang + 1;
              if (jvert > nvert - 1)
                jvert = jvert - nvert;
              kvert = jvert + 1;
              if (kvert > nvert - 1)
                kvert = kvert - nvert;
              
              /* Calculate mid-point of current triangle */
              for (icoord = 0; icoord < 3; icoord++)
                {
                  mid_point[icoord]
                    = (polygon_vertex_ptr[ivert]->coords[icoord]
                       + polygon_vertex_ptr[jvert]->coords[icoord]
                       + polygon_vertex_ptr[kvert]->coords[icoord])
                      / 3.0;
                }
              
              /* Initialise average score */
              ave_score =  0.0;
              navpoint = 0;
              
              /* Loop over all the cell vertices to calculate an
                 average score at the triangle's midpoint */
                  for (ipoint = 0; ipoint < 8; ipoint++)
                  {
                      /* Get this cell-point's coordinates */
                    x = corner_coords[0]
                      + grid_sep * cell_points[ipoint][0];
                    y = corner_coords[1]
                      + grid_sep * cell_points[ipoint][1];
                    z = corner_coords[2]
                      + grid_sep * cell_points[ipoint][2];
                    
                    /* Calculate the distance between the mid-point
                       and the cell-point */
                    x = x - mid_point[0];
                    y = y - mid_point[1];
                    z = z - mid_point[2];
                    distance =  sqrt(x * x + y * y + z * z);
                    
                    /* Apply distance as a factor to the value at
                       the given grid-point */
                    if (distance < 2.0 * grid_sep)
                      {
                        ave_score = ave_score
                          +  (2.0 * grid_sep - distance)
                            * corner_value[ipoint];
                        navpoint++;
                      }
                  }
              
              /* Calculate the averaged score at the triangle's
                 mid-point */
              if (navpoint > 0)
                {
                  /* Get average score */
                  ave_score = ave_score / navpoint;
                }
              else
                ave_score =  0.0;
              
              /* Calculate least-squares deviation from ideal value
                 (namely the contour level) */
              score = score
                + (ave_score - contour_level)
                  * (ave_score - contour_level);
            }
          
          /* If this is the best score so far, save it */
          if (iloop == 0 || score < min_score)
            {
              /* Store score and which loop produced it */
              min_score = score;
              bloop = iloop;
            }
        }
      
      /* Loop to generate the required number of triangles */
      for (itriang = 0; itriang < nvert - 2; itriang++)
        {
          /* Determine the vertices of the current triangle */
          ivert = bloop;
          jvert = ivert + itriang + 1;
          if (jvert > nvert - 1)
            jvert = jvert - nvert;
          kvert = jvert + 1;
          if (kvert > nvert - 1)
            kvert = kvert - nvert;

          /* Create a polygon record */
          polygon_ptr
            = create_polygon_record(&first_polygon_ptr,&last_polygon_ptr,
                                    npolygons,polygon_vertex_ptr[ivert],
                                    polygon_vertex_ptr[jvert],
                                    polygon_vertex_ptr[kvert],
                                    store_low_point[ivert]);

          /* If this is the first polygon, perform full check on vertex
             order */
          if (itriang == 0)
            swap
              = check_vertex_order(polygon_ptr,npoints,valmin,valmax,
                                   grid_sep,store_low_point[ivert],v);

          /* Otherwise, order in reverse order to previous triangle */
          else
            {
              /* If vertex order wrong, swap the 1st and 3rd vertices */
              if (swap == TRUE)
                {
                  /* Swap the first and third vertices */
                  vertex_ptr = polygon_ptr->vertex_ptr[0];
                  polygon_ptr->vertex_ptr[0] = polygon_ptr->vertex_ptr[2];
                  polygon_ptr->vertex_ptr[2] = vertex_ptr;
                }

              /* Get the coordinates of the vertices */
              get_vertex_coords(polygon_ptr,v);
            }

          /* Calculate normal to polygon and update vertex normal
             accumulators */
          update_normals(polygon_ptr,v);

          /* Increment count of polygons */
          npolygons++;
        }
    }

  /* Return updated polygon pointers */
  *fst_polygon_ptr = first_polygon_ptr;
  *lst_polygon_ptr = last_polygon_ptr;

  /* Return the current number of polygons */
  return(npolygons);
}
/***********************************************************************
  
make_polygons  -  Convert the lines between contour cross-over
                  points into polygons within current cell (There
                  can be up to 4 polygons crossing the cell).
  
***********************************************************************/

int make_polygons(int nlines,struct vertex *line_vertex_ptr[12][2],
                  int cell_line[12][2],int contour_level,
                  int low_point[12],double corner_coords[3],
                  int corner_value[8],double grid_sep,
                  int npoints[3],double valmin[3],double valmax[3],
                  struct polygon **fst_polygon_ptr,
                  struct polygon **lst_polygon_ptr,int *npol)
{
  int iline, nfound, non_null, old_edge, on_edge;
  int npolygons, nvert;
  int store_low_point[12];
  int done, match, ok;

  struct polygon *polygon_ptr, *first_polygon_ptr, *last_polygon_ptr;
  struct polygon *prev_polygon_ptr;
  struct vertex *polygon_vertex_ptr[12];
  struct vertex *end_vertex_ptr, *old_end_vertex_ptr;

  /* Initialise variables */
  done = FALSE;
  ok = TRUE;
  non_null = 1;
  npolygons = *npol;
  nvert = 0;

  /* Initialise pointers */
  end_vertex_ptr = NULL;
  polygon_ptr = NULL;
  first_polygon_ptr = *fst_polygon_ptr;
  last_polygon_ptr = *lst_polygon_ptr;
  prev_polygon_ptr = *lst_polygon_ptr;

  /* Loop until all lines have been processed and the 1 to 4
     polygons for this cell have been generated */
  while (non_null > 0)
    {
      /* Set flag giving whether anything found in next pass
         through all the lines */
      nfound = 0;
      non_null = 0;

      /* Loop through all the lines for current cell */
      for (iline = 0; iline < nlines; iline++)
        {
          /* If line OK to use, then continue */
          if (line_vertex_ptr[iline][0] != NULL)
            {
              /* Increment count of available lines remaining to be
                 processed */
              non_null++;

              /* Initialise flag */
              match = FALSE;

              /* If we're not looking for a continuation of the
                 current polygon, then start a new one */
              if (end_vertex_ptr == NULL)
                {
                  /* Store one end of the line as the current
                     search vertex for the new polygon */
                  end_vertex_ptr = line_vertex_ptr[iline][0];

                  /* Store which edge this vertex came from */
                  on_edge = cell_line[iline][0];

                  /* Initialise count of vertices for this polygon */
                  nvert = 0;
                }

              /* Initialise variables */
              old_end_vertex_ptr = end_vertex_ptr;
              old_edge = on_edge;

              /* Check if either end of the line matches one of the end
                 vertices of the current polygon */
              if (line_vertex_ptr[iline][0] == end_vertex_ptr)
                {
                  /* Store the other end of the line */
                  end_vertex_ptr = line_vertex_ptr[iline][1];
                  on_edge = cell_line[iline][1];

                  /* Set flag to indicate that match was found */
                  match = TRUE;
                }

              /* Check for match with other end of the current line */
              else if (line_vertex_ptr[iline][1] == end_vertex_ptr)
                {
                  /* Store the first end of the line */
                  end_vertex_ptr = line_vertex_ptr[iline][0];
                  on_edge = cell_line[iline][0];

                  /* Set flag to indicate that match was found */
                  match = TRUE;
                }

              /* If match found, create new line */
              if (match == TRUE)
                {
                  /* Store this polygon's next vertex */
                  polygon_vertex_ptr[nvert] = end_vertex_ptr;

                  /* Save the location of the best low-point
                     to use in checking polygon vertex-ordering */
                  store_low_point[nvert] = low_point[old_edge];

                  /* Increment count of polygon's vertices */
                  nvert++;

                  /* Blank out the line just done */
                  line_vertex_ptr[iline][0] = NULL;

                  /* Increment count of lines found */
                  nfound++;
                }
            }
        }

      /* If haven't found any more vertices for current polygon, reset
         so that next unused line starts a new polygon */
      if (nfound == 0)
        {
          /* Reset end vertex pointer */
          end_vertex_ptr = NULL;

          /* Check that the last polygon created has at least 3 vertices
             and, where it has more than 3, break it up into triangles */
          npolygons
            = generate_triangles(&first_polygon_ptr,&last_polygon_ptr,
                                 npolygons,polygon_vertex_ptr,nvert,
                                 store_low_point,corner_coords,
                                 corner_value,grid_sep,contour_level,
                                 npoints,valmin,valmax);
        }
    }
  
  /* Return pointer to first polygon */
  *fst_polygon_ptr = first_polygon_ptr;
  *lst_polygon_ptr = last_polygon_ptr;
  
  /* Return current number of polygons */
  *npol = npolygons;

  /* Return whether generation went well */
  return(ok);
}
/***********************************************************************

extract_grid_value  -  Extract the current grid value and mask-flag, if
                       any

***********************************************************************/

int extract_grid_value(int grid_value,int *masked)
{
  int actual_value;

  /* Decode value and flag */
  *masked = grid_value % 10;
  actual_value = grid_value / 10;

  /* Return the actual value */
  return(actual_value);
}
/***********************************************************************

get_grid_value  -  Get the grid value, adjusting if value happens to be
                   identical to the contouring level

***********************************************************************/

int get_grid_value(int grid_value,int *mskd,int contour_level)
{
  int value;
  int masked;

  /* Get mask flag */
  masked = *mskd;

  /* Get the grid value */
  value = extract_grid_value(grid_value,&masked);

  /* If value equals the contour level, apply tiny adjustment */
  if (value - contour_level < 1000)
    value = value + 1000;

  /* Return mask flag */
  *mskd = masked;

  /* Return adjusted value */
  return(value);
}
/***********************************************************************
  
generate_polygons  -  Generate all the polygons (triangles) to best
                      approximate the contour surface within each grid
                      cell
  
***********************************************************************/

int generate_polygons(int *grid,int npoints[3],double valmin[3],
                      double valmax[3],double grid_sep,int contour_level,
                      struct grid_limits glimits,
                      struct vertex **vertex_index_ptr,int nvertices,
                      struct polygon **fst_polygon_ptr,
                      struct atom *dummy_atom_ptr,int buried_vertices_only)
{
  int edge_number[12], face[12][2], ncross_face[6];
  int edge_no, edge_neighbour, i, iedge, iface, ipoint, iside, j, k;
  int cell_line[12][2], on_edge[4];
  int line_i, line_j, line_k;
  int check_edge[12], low_point[12];
  int icorner, istart, iend, iline, jstart, jend, kstart, kend;
  int ncross, nlines, nover, nunder, npolygons, nsides, nx, ny, nz;
  int corner_value[8], cross_add, grid_value, masked;
  int missing_atom_details, ok, wanted;

  double corner_coords[3];
  double x, y, z;
  
  struct polygon *first_polygon_ptr, *last_polygon_ptr;
  struct vertex *vertex_ptr, *edge_vertex_ptr[12], *side_vertex_ptr[4];
  struct vertex *line_vertex_ptr[12][2];

  static int grid_pos[12][2][3]
    = { 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        1, 0, 0, 1, 1, 0,
        1, 0, 0, 1, 0, 1,
        0, 1, 0, 1, 1, 0,
        0, 1, 0, 0, 1, 1,
        1, 1, 0, 1, 1, 1,
        0, 0, 1, 1, 0, 1,
        0, 0, 1, 0, 1, 1,
        1, 0, 1, 1, 1, 1,
        0, 1, 1, 1, 1, 1 };
  
  static int last_face_point[6][3]
    = { 1, 0, 0,
        1, 0, 0,
        0, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1 };
  
  static int edges_per_face[6][4]
    = { 0,  1,  5,  3,
        0,  2,  8,  4,
        1,  6,  9,  2,
        3,  4, 10,  7,
        5,  7, 11,  6,
        8, 10, 11,  9 };

  static int cell_points[8][3]
      = { 0, 0, 0,
          0, 0, 1,
          0, 1, 0,
          0, 1, 1,
          1, 0, 0,
          1, 0, 1,
          1, 1, 0,
          1, 1, 1 };

  static int corner_to_edge[8][3]
    = { 0,  1,  2,
        2,  8,  9,
        1,  5,  6,
        6,  9, 11,
        0,  3,  4,
        4,  8, 10,
        3,  5,  7,
        7, 10, 11 };
  
/*-------------------------------------------------------------------------
    
 Numbering of faces                     Numbering of edges
                                                                 10
    +                 +--------+           +                 +--------+
    |                /   5    /|           |             11 /      8 /|
    |               /        / |          7|               /   9    / |
    |    3         +--------+  |           |              +--------+  | 4
  4 |              |        | 1|           |    3         |        |  |
    +--------+     |   2    |  +           +--------+    6|       2|  +
   /               |        | /         5 /               |        | /
  /    0           |        |/           /                |    1   |/ 0
 +                 +--------+           +                 +--------+  
     Back             Front                 Back             Front
    
 Numbering of corners
                                                                 10
    +                 7--------5
    |                /        /|
    |               /        / |
    |              3--------1  |
    |              |        |  |
    6--------+     |        |  4
   /               |        | / 
  /                |        |/  
 +                 2--------0   
     Back             Front     
    
-------------------------------------------------------------------------*/
  
  /* Initialise polygon pointers */
  first_polygon_ptr = NULL;
  last_polygon_ptr = NULL;

  /* Initialise variables */
  iedge = 0;
  ipoint = 0;
  npolygons = 0;
  ok = TRUE;
  wanted = TRUE;

  /* Define range of non-zero array */
  istart = glimits.minpos[0];
  iend = glimits.maxpos[0] + 1;
  if (iend > npoints[0] - 1)
      iend = npoints[0] - 1;
  jstart = glimits.minpos[1];
  jend = glimits.maxpos[1] + 1;
  if (jend > npoints[1] - 1)
      jend = npoints[1] - 1;
  kstart = glimits.minpos[2];
  kend = glimits.maxpos[2] + 1;
  if (kend > npoints[2] - 1)
      kend = npoints[2] - 1;
  
  /* Determine size of reduced array region */
  nx = iend - istart;
  ny = jend - jstart;
  nz = kend - kstart;

  /* Initialise array determining the two faces to which each edge belongs */
  for (iedge = 0; iedge < 12; iedge++)
    {
      /* Initialise face numbers */
      face[iedge][0] = -1;
      face[iedge][1] = -1;

      /* Loop over the 6 faces */
      for (iface = 0; iface < 6; iface++)
        {
          /* Loop over this face's four sides */
          for (iside = 0; iside < 4; iside++)
            {
              /* If this face's current side corresponds to the current edge,
                 store the face as being associated with the current edge */
              if (edges_per_face[iface][iside] == iedge)
                {
                  /* If this is the first of the edge's two faces, store
                     as the first face */
                  if (face[iedge][0] == -1)
                    face[iedge][0] = iface;

                  /* Otherwise, store as the edge's second face */
                  else
                    face[iedge][1] = iface;
                }
            }
        }
    }
  
  /* Initialise the edge number */
  edge_no = 0;
  
  /* Loop through the x points */
  for (i = istart; i < iend; i++)
    {
      /* Get the x-coordinate */
      x = valmin[0] + i * grid_sep;

      /* Loop through the y points */
      for (j = jstart; j < jend; j++)
        {
          /* Get the y-coordinate */
          y = valmin[1] + j * grid_sep;

          /* Loop through the z points */
          for (k = kstart; k < kend && ok == TRUE; k++)
            {
              /* Get the z-coordinate */
              z = valmin[2] + k * grid_sep;

              /* Store the coordinates of this cell's origin */
              corner_coords[0] = x;
              corner_coords[1] = y;
              corner_coords[2] = z;

              /* Initialise counts of cross-overs in this cell */
              ncross = 0;
              for (iface = 0; iface < 6; iface++)
                ncross_face[iface] = 0;

              /* Initialise array showing which edges have a cross-over
                 and are worth checking out */
              for (iedge = 0; iedge < 12; iedge++)
                check_edge[iedge] = 0;

              /* Initialise counts of grid-points over and under the
                 contour level */
              nover = 0;
              nunder = 0;

              /* Initialise buried-only flag */
              wanted = TRUE;

              /* Loop over the 8 corners of this cell */
              for (icorner = 0; icorner < 8; icorner++)
                {
                  /* Get the array position of this corner */
                  line_i = i + cell_points[icorner][0];
                  line_j = j + cell_points[icorner][1];
                  line_k = k + cell_points[icorner][2];

                  /* Get the array location of this grid-point */
                  ipoint = line_i * npoints[1] * npoints[2]
                    + line_j * npoints[2] + line_k;

                  /* Get the grid-value at the current point */
                  corner_value[icorner]
                    = get_grid_value(grid[ipoint],&masked,contour_level);

                  /* Set cross-over indicator to reflect whether this
                     corner-point is above or below the contour level */
                  if (corner_value[icorner] > contour_level)
                    {
                      /* Corner's value is over the contour level */
                      cross_add = 1;
                      nover++;
                    }
                  else
                    {
                      /* Corner's value is below the contour level */
                      cross_add = -1;
                      nunder++;
                    }

                  /* Update the check-edge markers of the 3 edges this
                     corner belongs to */
                  for (edge_neighbour = 0; edge_neighbour < 3;
                       edge_neighbour++)
                    {
                      /* Get the appropriate edge */
                      iedge = corner_to_edge[icorner][edge_neighbour];

                      /* Update the check-edge marker */
                      check_edge[iedge] = check_edge[iedge] + cross_add;
                    }
                }

              /* If have at least one low-point or at least one
                 high-point in this cell, then need to check the crossed
                 edges for valid vertices */
              if (nover > 0 && nunder > 0)
                {
                  /* Calculate the first of the 12 edge-numbers making up
                     the cell commencing at the current grid-point */
                  edge_no = 3 * ((i - istart) * ny * nz
                                 + (j - jstart) * nz + k - kstart);

                  /* Calculate all 12 edge-numbers */
                  edge_number[0] = edge_no;
                  edge_number[1] = edge_no + 1;
                  edge_number[2] = edge_no + 2;
                  edge_number[3] = edge_no + 3 * ny * nz + 1;
                  edge_number[4] = edge_no + 3 * ny * nz + 2;
                  edge_number[5] = edge_no + 3 * nz;
                  edge_number[6] = edge_no + 3 * nz + 2;
                  edge_number[7] = edge_no + 3 * ny * nz + 3 * nz + 2;
                  edge_number[8] = edge_no + 3;
                  edge_number[9] = edge_no + 4;
                  edge_number[10] = edge_no + 3 * ny * nz + 4;
                  edge_number[11] = edge_no + 3 * nz + 3;

                  /* Loop over the edges in the current cell to check
                     which contain vertices */
                  for (iedge = 0; iedge < 12; iedge++)
                    {
                      /* Initialiase reference to nearest low-point
                         (ie nearest grid-point below the contour level)
                         to any cross-over on this edge */
                      low_point[iedge] = -1;

                      /* Initialise edge's vertex pointer */
                      vertex_ptr = NULL;

                      /* Only need to check this edge if its check-number
                         is zero (ie it had one corner over the contour
                         level and one below) */
                      if (check_edge[iedge] == 0)
                        {
                          /* Get the edge number */
                          edge_no = edge_number[iedge];

                          /* Perform a binary search on the vertices to see
                             if the current edge has a cross-over vertex
                             on it */
                          vertex_ptr
                            = binary_edge_search(edge_no,vertex_index_ptr,
                                                 nvertices);

                          /* If buried-only vertices are required,
                             check whether this is one */
                          if (vertex_ptr != NULL &&
                              buried_vertices_only == TRUE)
                            {
                              /* Check buried ratio */
                              if (vertex_ptr->ratio < 0.01)
                                {
                                  /* Set flags that this vertex, and
                                     indeed none of the ones in this
                                     cell are required */
                                  wanted = FALSE;
                                  vertex_ptr = NULL;
                                }
                            }

                          /* If edge contains vertex, store its location */
                          if (vertex_ptr != NULL)
                            {
                              /* Update total count of vertices and counts
                                 for each face */
                              ncross++;
                              iface = face[iedge][0];
                              ncross_face[iface]++;
                              iface = face[iedge][1];
                              ncross_face[iface]++;

                              /* Get the array position of one end of
                                 the line */
                              line_i = i + grid_pos[iedge][0][0];
                              line_j = j + grid_pos[iedge][0][1];
                              line_k = k + grid_pos[iedge][0][2];

                              /* Get the array location of this grid-point */
                              ipoint = line_i * npoints[1] * npoints[2]
                                + line_j * npoints[2] + line_k;

                              /* Get the grid-value at the current point */
                              grid_value
                                = get_grid_value(grid[ipoint],&masked,
                                                 contour_level);

                              /* If this is not a low-point, then want the
                                 array position of the other end of the line */
                              if (grid_value > contour_level)
                                {
                                  /* Get the array position of the other end
                                     of the line */
                                  line_i = i + grid_pos[iedge][1][0];
                                  line_j = j + grid_pos[iedge][1][1];
                                  line_k = k + grid_pos[iedge][1][2];

                                  /* Get the array location of this
                                     grid-point */
                                  ipoint = line_i * npoints[1] * npoints[2]
                                    + line_j * npoints[2] + line_k;
                                }

                              /* Store the array location of low-point
                                 just found */
                              low_point[iedge] = ipoint;
                            }
                        }

                      /* Store this edge's vertex pointer */
                      edge_vertex_ptr[iedge] = vertex_ptr;
                    }
                }

              /* If we have any cross-overs in this cell, then process */
              if (ncross > 0 && wanted == TRUE)
                {
                  /* Initialise count of lines crossing through this cell */
                  nlines = 0;

                  /* Process each of the 6 faces */
                  for (iface = 0; iface < 6; iface++)
                    {
                      /* Determine which edges make up the current face */
                      nsides = 0;
                      for (iside = 0; iside < 4; iside++)
                        {
                          /* Get the current edge number */
                          iedge = edges_per_face[iface][iside];

                          /* Get the pointer to the vertex on this
                             edge */
                          vertex_ptr = edge_vertex_ptr[iedge];

                          /* If vertex pointer not NULL, then have a
                             cross-over on this edge */
                          if (vertex_ptr != NULL)
                            {
                              /* Store details */
                              side_vertex_ptr[nsides] = vertex_ptr;
                              on_edge[nsides] = iedge;

                              /* Increment count of edges corssed */
                              nsides++;
                            }
                        }

                      /* Get the last point traced out by the edges
                         of the current face */
                      line_i = i + last_face_point[iface][0];
                      line_j = j + last_face_point[iface][1];
                      line_k = k + last_face_point[iface][2];

                      /* Get the array location of this grid-point */
                      ipoint = line_i * npoints[1] * npoints[2]
                        + line_j * npoints[2] + line_k;

                      /* Get the grid-value at the current point */
                      grid_value
                        = get_grid_value(grid[ipoint],&masked,contour_level);

                      /* If there are any cross-over points in the
                         current face, then join them by lines */
                      nlines
                        = make_line(ncross_face[iface],nlines,cell_line,
                                    line_vertex_ptr,side_vertex_ptr,
                                    on_edge,grid_value,contour_level);
                    }

                  /* If have any lines for the current cell convert into
                     polygons */
                  if (nlines > 0)
                    {
                      /* Initialise flag */
                      missing_atom_details = FALSE;

                      /* Check that all the vertices have associated
                         atom details */
                      for (iline = 0; iline < nlines &&
                           missing_atom_details == FALSE; iline++)
                        {
                          /* Get first of this line's vertices */
                          vertex_ptr = line_vertex_ptr[iline][0];

                          /* If this vertex does not have an associated
                             atom, then link to dummy atom */
                          if (vertex_ptr->atom_ptr == NULL)
                            vertex_ptr->atom_ptr = dummy_atom_ptr;
                        }

                      /* Create a set of polygons for the current cell */
                      if (missing_atom_details == FALSE)
                        {
                          ok
                            = make_polygons(nlines,line_vertex_ptr,
                                            cell_line,contour_level,low_point,
                                            corner_coords,corner_value,
                                            grid_sep,npoints,valmin,
                                            valmax,&first_polygon_ptr,
                                            &last_polygon_ptr,&npolygons);
                        }
                    }
                }
            }
        }
    }
  
  /* Return current pointers */
  *fst_polygon_ptr = first_polygon_ptr;
  
  /* Return number of polygons generated */
  return(npolygons);
}
/***********************************************************************

vecprd  -  Calculates the vector product (vx,vy,vz) of two vectors 
           (px,py,pz) and (qx,qy,qz)

***********************************************************************/

void vecprd(double px, double py, double pz, double qx, double qy, double qz,
            double *vx, double *vy, double *vz)
{
  *vx = py * qz - pz * qy;
  *vy = pz * qx - px * qz;
  *vz = px * qy - py * qx;
}
/***********************************************************************

dot3  -  Returns the scalar product of the two vectors (x1,y1,z1) 
         and (x2,y2,z2)

***********************************************************************/

double dot3(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double dot_prod;

  dot_prod = x1 * x2 + y1 * y2 + z1 * z2;
  return(dot_prod);
}
/***********************************************************************

orien3  -  Returns the orientation of the polygon with the consecutive
           vertices (x1,y1,z1), (x2,y2,z2), and (x3,y3,z3) as viewed
           from (ex,ey,ez).
                      -1  :  clockwise orientation
                      +1  :  anti-clockwise orientation (desired order)
                       0  :  degenerate (line or point)

 [Adapted from Angell I O & Griffith G (1987). High-resolution Computer
  Graphics using FORTRAN 77. Macmillan Education Ltd, Basingstoke, UK.]

***********************************************************************/

int orien3(double x1, double y1, double z1, double x2, double y2, double z2,
           double x3, double y3, double z3, double ex, double ey, double ez)
{
    double a, b, c, dx1, dy1, dz1, dx2, dy2, dz2, fe;
    int orientation;

    /* Initialise variables */
    dx1 = x2 - x1;
    dy1 = y2 - y1;
    dz1 = z2 - z1;
    dx2 = x3 - x2;
    dy2 = y3 - y2;
    dz2 = z3 - z2;

    /* Calculate triple scalar product */
    vecprd(dx1,dy1,dz1,dx2,dy2,dz2,&a,&b,&c);
    fe = dot3(a,b,c,ex-x1,ey-y1,ez-z1);
    if (fe < 0.0)
      orientation = -1;
    else
      orientation = 1;
    return(orientation);
}
/***********************************************************************

calculate_normal  -  Calculate the normal of the plane defined by the
                     three sets of coordinates supplied

***********************************************************************/

void calculate_normal(double x1,double y1,double z1,
                      double x2,double y2,double z2,
                      double x3,double y3,double z3,double normal[3])
{
  double ax, ay, az, bx, by, bz, nx, ny, nz;
  double diff, magnitude;

  /* Define two vectors in the plane of the triangle */
  ax = x2 - x1;
  ay = y2 - y1;
  az = z2 - z1;
  bx = x3 - x1;
  by = y3 - y1;
  bz = z3 - z1;

  /* Calculate the normal to this triangle */
  nx = ay * bz - az * by;
  ny = az * bx - ax * bz;
  nz = ax * by - ay * bx;

  /* Calculate magnitude of vector */
  magnitude =  sqrt(nx * nx + ny * ny + nz * nz);

  /* Calculate absolute value */
  diff =  fabs( magnitude);

  /* Normalise it */
  if (diff >  0.000001)
    {
      normal[0] = nx / magnitude;
      normal[1] = ny / magnitude;
      normal[2] = nz / magnitude;
    }
  else
    {
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 0.0;
    }
}
/***********************************************************************

order_polygon_vertices  -  Order each polygon's vertices so that they
                           are in an anti-clockwise sense as viewed
                           from outside

***********************************************************************/

void order_polygon_vertices(struct polygon *first_polygon_ptr,
                            int npoints[3],double valmin[3],
                            double valmax[3],double grid_sep)
{
  int icoord, iorder, ivert;
  int calc, i, ipoint, ipos, j, k;
  
  double normal[3], swap, v[9], x, y, z;
  
  struct polygon *polygon_ptr;
  struct vertex *vertex_ptr;
  
  /* Get first polygon pointer */
  polygon_ptr = first_polygon_ptr;
  
  /* Loop through all the vertices */
  while (polygon_ptr != NULL)
    {
      /* Get the reference to the outside point */
      ipoint = polygon_ptr->outside_point;

      /* Convert into array location */
      i = ipoint / (npoints[1] * npoints[2]);
      calc = ipoint - i * npoints[1] * npoints[2];
      j = calc / npoints[2];
      k = calc - j * npoints[2];

      /* Calculate the point's x-, y- and z-coordinates */
      x = valmin[0] + i * grid_sep;
      y = valmin[1] + j * grid_sep;
      z = valmin[2] + k * grid_sep;

      /* Initialise */
      ipos = 0;

      /* Loop over the polygon's three vertices */
      for (ivert = 0; ivert < 3; ivert++)
        {
          /* Get pointer to this vertex */
          vertex_ptr = polygon_ptr->vertex_ptr[ivert];

          /* Write out this vertex's coordinates */
          for (icoord = 0; icoord < 3; icoord++)
            {
              /* Store the coordinate */
              v[ipos] = vertex_ptr->coords[icoord];
              ipos++;
            }
        }

      /* Check the vertex order */
      iorder = orien3(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],
                      x,y,z);

      /* If the order is anti-clockwise, then reverse order of vertices
         to make order clockwise */
      if (iorder == 1)
        {
          /* Swap the first and third vertices */
          vertex_ptr = polygon_ptr->vertex_ptr[0];
          polygon_ptr->vertex_ptr[0] = polygon_ptr->vertex_ptr[2];
          polygon_ptr->vertex_ptr[2] = vertex_ptr;

          /* Swap the coordinates as well */
          for (icoord = 0; icoord < 3; icoord++)
            {
              swap = v[icoord];
              v[icoord] = v[icoord + 6];
              v[icoord + 6] = swap;
            }
        }

      /* Calculate normal to surface of polygon */
      calculate_normal(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],
                       normal);

      /* Add normals to each vertex's accumulator for deriving an
         average surface normal at that vector */
      for (ivert = 0; ivert < 3; ivert++)
        {
          /* Get pointer to this vertex */
          vertex_ptr = polygon_ptr->vertex_ptr[ivert];

          /* Update the normal */
          for (icoord = 0; icoord < 3; icoord++)
            {
              /* Store the coordinate */
              vertex_ptr->normal[icoord]
                = vertex_ptr->normal[icoord] + normal[icoord];
            }

          /* Increment count of polygons to which this vertex belongs */
          vertex_ptr->npolygons++;
        }

      /* Get next polygon record */
      polygon_ptr = polygon_ptr->next_polygon_ptr;
    }
}
/***********************************************************************

write_raster3d_headers  -  Write out the header records to the output
                           Raster3D file

***********************************************************************/

void write_raster3d_headers(FILE *outfile_ptr,char *prog_name,
                            char *background_colour,int sizex,int sizey)
{
  double rgb[3];

  /* Convert background colour to RGB values */
   get_colour_rgb(background_colour,rgb);

  /* Initialise variables */
  fprintf(outfile_ptr,"Generated by %s. R A Laskowski (2002).\n",
          prog_name);
  fprintf(outfile_ptr,"%d %d     tiles in x,y\n",sizex,sizey);
  fprintf(outfile_ptr,"24 24     pixels (x,y) per tile\n");
  fprintf(outfile_ptr,"3         3x3 virtual pixels -> 2x2 pixels\n");
  fprintf(outfile_ptr,"%5.3f %5.3f %5.3f   Background colour\n",
      rgb[0],rgb[1],rgb[2]);
  fprintf(outfile_ptr,"T         shadows cast\n");
  fprintf(outfile_ptr,"25        Phong power\n");
  fprintf(outfile_ptr,"0.15      secondary light contribution\n");
  fprintf(outfile_ptr,"0.05      ambient light contribution\n");
  fprintf(outfile_ptr,"0.25      specular reflection component\n");
  fprintf(outfile_ptr,"0.0       eye position\n");
  fprintf(outfile_ptr,"1.0 1.0 1.0  main light source position\n");
  fprintf(outfile_ptr,"1 0 0 0   input coordinate, radius ");
  fprintf(outfile_ptr,"transformation\n");
  fprintf(outfile_ptr,"0 1 0 0\n");
  fprintf(outfile_ptr,"0 0 1 0\n");
  fprintf(outfile_ptr,"0 0 0 1\n");
  fprintf(outfile_ptr,"3         mixed objects\n");
  fprintf(outfile_ptr,"*\n");
  fprintf(outfile_ptr,"*\n");
  fprintf(outfile_ptr,"*\n");
}
/***********************************************************************

get_colour_cpk  -  Get the rgb colour for the given element

***********************************************************************/

int get_colour_cpk(char *element,double rgb[3])
{
#define NELEMENTS        5

  char c_name[COLNAM_LEN];

  int colour, ielem;

  static char *element_list[] = {
    " X", " C", " O", " N", " S"
  };

  static char *colour_list[] = {
    "magenta", "white", "red", "blue", "yellow"
  };

  /* Initialise colour to first colour */
  colour = 0;

  /* Search elements list for given element */
  for (ielem = 0; ielem < NELEMENTS; ielem++)
    if (!strcmp(element,element_list[ielem]))
      colour = ielem;

  /* Get the corresponding RGB triplet for element's colour */
  strcpy(c_name,colour_list[colour]);
  get_colour_rgb(c_name,rgb);

  /* Return the corresponding colour number */
  return(colour);
}
/***********************************************************************

get_colour_restype  -  Get vertex colour from associated residue type

***********************************************************************/

int get_colour_restype(char *res_name,double rgb[3])
{
  char c_name[COLNAM_LEN];

  int colour, iamino, ires;

  static char *aa_name[] = {
    "XXX", "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
    "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
    "SER", "THR", "VAL", "TRP", "TYR" };

  static char *colour_name[] = { 
    "white", "skyblue", "red", "green", "grey", "purple", "brown", "yellow" };
  
  static char *colour_desc[] = { 
    "unknown", "positive", "negative", "neutral", "aliphatic", "aromatic",
    "Pro&Gly", "cysteine" };
  
  static int res_type[NAMINO + 1]
    = { 0, 4, 7, 2, 2, 5, 6, 1, 4, 1, 4, 4, 3, 6, 3, 1, 3, 3, 4, 5, 5 };
  
  /* Initialise colour to first colour */
  colour = 0;
  iamino = 0;

  /* Search residues list for given residue name */
  for (ires = 1; ires < NAMINO + 1; ires++)
    if (!strcmp(res_name,aa_name[ires]))
      iamino = ires;

  /* Get corresponding residue colour */
  colour = res_type[iamino];

  /* Get the corresponding RGB triplet for residue's colour */
  strcpy(c_name,colour_name[colour]);
  get_colour_rgb(c_name,rgb);

  /* Return the corresponding colour number */
  return(colour);
}
/***********************************************************************

get_colour_conserv  -  Get vertex colour from associated residue's
                       conservation score

***********************************************************************/

int get_colour_conserv(double conservation,double rgb[3],
                       char *conservation_colours)
{
  char c_name[COLNAM_LEN], number_string[20];

  int colour, iband, ipos;

  static char *colour_name[] = { 
    "white", "blue", "purple", "skyblue", "cyan", "grey", "green",
    "yellow", "orange", "pink", "red" };

  /* Initialise variables */
  iband = 0;
  colour = 0;

  /* Loop until find which band the conservation score falls into */
  for (ipos = 1; ipos < 11 && iband == 0; ipos++)
    {
      if (conservation <= RANGE[ipos])
        iband = ipos;
    }

  /* If have a set of user-defined colours, then retrieve */
  if (iband > 0 && conservation_colours[0] != '\0')
    {
      /* Get colour corresponding to given band */
      number_string[0] = conservation_colours[iband - 1];
      number_string[1] = '\0';

      /* If colour-number is an "A", replace by 10 */
      if (number_string[0] == 'A')
        colour = 10;
      
      /* Otherwise, retrieve the colour-number */
      else
        colour = atoi(number_string);
    }

  /* Otherwise, use colour corresponding to given band */
  else
    colour = iband;

  /* Get colour name */
  strcpy(c_name,colour_name[colour]);
  get_colour_rgb(c_name,rgb);

  /* Return the corresponding colour number */
  return(colour);
}
/***********************************************************************

get_chain_colour  -  Determine the colour for the given chain

***********************************************************************/

int get_chain_colour_numbers(char chain,char colour_name[COLNAM_LEN],
                             char colour_numbers[COLNO_LEN])
{
  int colour, ichain, nchains;

  static char chain_list[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890";
  static char *col_name[] = {
    "purple", "red", "brown", "pink", "blue",
    "green", "cyan", "yellow"
  };
  static char *col_numbers[] = {
    "[178,51,255]", "[255,0,0]", "[128,0,0]", "[255,128,255]", "[0,0,255]",
    "[0,255,0]", "[128,255,255]", "[255,255,0]"
  };

  /* Initialise variable */
  colour = 0;
  colour_name[0] = '\0';
  colour_numbers[0] = '\0';

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
          colour = colour % MAX_OBJ_COL;
          strcpy(colour_name,col_name[colour]);
          strcpy(colour_numbers,col_numbers[colour]);
        }
    }

  /* Return the colour number */
  return(colour);
}
/***********************************************************************

write_raster3d_polygons  -  Write out the polygons to the output
                            Raster3D file

***********************************************************************/

void write_raster3d_polygons(FILE *outfile_ptr,
                             struct polygon *first_polygon_ptr,
                             double c_of_g[3],double rascale,char chain,
                             int vertex_colouring,
                             char *conservation_colours,
                             char *col_name)
{
  char colour_name[COLNAM_LEN], colour_numbers[COLNO_LEN];
  char element[3];

  int icoord, ivert;
  int colour[3], pcolour;
  int vertices_equal;

  double coord, normal[3][3], rgb[3][3];

  struct polygon *polygon_ptr;
  struct residue *residue_ptr;
  struct vertex *vertex_ptr;

  /* Get first polygon pointer */
  polygon_ptr = first_polygon_ptr;

  /* Loop through all the vertices */
  while (polygon_ptr != NULL)
    {
      /* Write out object type */
      fprintf(outfile_ptr,"1\n");

      /* Loop over the polygon's three vertices */
      for (ivert = 0; ivert < 3; ivert++)
        {
          /* Get pointer to this vertex */
          vertex_ptr = polygon_ptr->vertex_ptr[ivert];

          /* Write out this vertex's coordinates */
          for (icoord = 0; icoord < 3; icoord++)
            {
              /* Get this coordinate value */
              coord = vertex_ptr->coords[icoord];

              /* Shift by centre of gravity and apply scaling factor */
              coord = rascale * (coord - c_of_g[icoord]);

              /* Write out */
              fprintf(outfile_ptr,"%8.4f ",coord);

              /* Store the vertex normal */
              normal[ivert][icoord] = vertex_ptr->normal[icoord];
            }

          /* Initialise colour */
          colour[ivert] = -1;

          /* Get the colour associated with this vertex */
          if (vertex_colouring == CPK)
            {
              /* Get the element associated with the given atom */
              strcpy(element,vertex_ptr->atom_ptr->element);

              /* Get element's RGB colour */
              colour[ivert] = get_colour_cpk(element,rgb[ivert]);
            }

          /* Colouring by residue type */
          else if (vertex_colouring == RESIDUE_TYPE)
            {
              /* Get residue pointer */
              residue_ptr = vertex_ptr->atom_ptr->residue_ptr;

              /* If not NULL, get residue's RGB colour */
              if (residue_ptr != NULL)
                colour[ivert]
                  = get_colour_restype(residue_ptr->res_name,rgb[ivert]);
            }

          /* Colouring by residue conservation */
          else if (vertex_colouring == CONSERVATION)
            {
              /* Get residue pointer */
              residue_ptr = vertex_ptr->atom_ptr->residue_ptr;

              /* If not NULL, get residue's RGB colour */
              if (residue_ptr != NULL)
                colour[ivert]
                  = get_colour_conserv(residue_ptr->conservation_score,
                                       rgb[ivert],
                                       conservation_colours);
            }

          /* Colouring by chain */
          else if (vertex_colouring == BY_CHAIN)
            {
              /* Get residue pointer */
              residue_ptr = vertex_ptr->atom_ptr->residue_ptr;

              /* If not NULL, get chain colour */
              if (residue_ptr != NULL)
                {
                  colour[ivert]
                    = get_chain_colour_numbers(chain,colour_name,
                                               colour_numbers);

                  /* Convert to RGB triplet */
                  get_colour_rgb(colour_name,rgb[ivert]);
                }
            }

          /* Colouring by cluster number */
          else if (vertex_colouring == STANDARD)
            {
              /* Convert supplied colour name to RGB triplet */
              get_colour_rgb(col_name,rgb[ivert]);
              colour[ivert] = 1;
            }

          /* If no colour assigned, set to default as white */
          if (colour[ivert] == -1)
            {
              colour[ivert] = 0;
              rgb[ivert][0] =  1.0;
              rgb[ivert][1] =  1.0;
              rgb[ivert][2] =  1.0;
            }
        }

      /* Close off the current line */
      fprintf(outfile_ptr,"\n");

      /* Determine polygon colour according to colours associated
         with each vertex */
      pcolour = 0;
      vertices_equal = FALSE;

      /* If all vertices are of the same colour, then only need to
         colour the polygon */
      if (colour[0] == colour[1] && colour[1] == colour[2])
        vertices_equal = TRUE;

      /* Otherwise, check which pair of vertices have the same
         colour */
      else if (colour[1] == colour[2])
        pcolour = 1;

      /* Write out the polygon colour */
      fprintf(outfile_ptr,"%7.4f %7.4f %7.4f\n",
              rgb[pcolour][0],rgb[pcolour][1],rgb[pcolour][2]);

      /* Write out the stored vertex normals */
      fprintf(outfile_ptr,"7\n");
      for (ivert = 0; ivert < 3; ivert++)
        {
          for (icoord = 0; icoord < 3; icoord++)
            {
              fprintf(outfile_ptr,"%8.3f ",normal[ivert][icoord]);
            }
        }

      /* Close off the current line */
      fprintf(outfile_ptr,"\n");

      /* If vertices are to be of differing colours, write these out */
      if (vertices_equal == FALSE)
        {
          /* Write out the object type */
          fprintf(outfile_ptr,"17\n");
      
          /* Loop over the polygon's three vertices */
          for (ivert = 0; ivert < 3; ivert++)
            {
              /* Write out this vertex's colour */
              fprintf(outfile_ptr,"%7.4f %7.4f %7.4f ",
                      rgb[ivert][0],rgb[ivert][1],rgb[ivert][2]);
            }

          /* Close off the current line */
          fprintf(outfile_ptr,"\n");
        }

      /* Get next polygon record */
      polygon_ptr = polygon_ptr->next_polygon_ptr;
    }
}
