/***********************************************************************

  ligps.c - PostScript functions called by LIGPLOT.

************************************************************************

Date:-         11 June 1999

Written by:-   Roman Laskowski

               Department of Crystallography
               Birkbeck College
               Malet Street, London WC1E 7HX, UK.


Version:-      Extracted from ligplot.c for v.4.2 of LIGPLOT.

Amendments:-   Amendments after v.4.2 are labelled by version number
               v.m.n--> and v.m.n<-- where m.n is the version
               number corresponding to the change

----------------------------------------------------------------------*/

/* Functions
   ---------

centre_of_mass
get_matrix
rotate_centres_of_mass
to_lower
to_upper
psarc_
psbbox_
psccol_
pscirc_
psclos_
pscolb_
pscomm_
pscrgb_
psctxt_
psdash_
pselip_
pshade_
psline_
pslwid_
psphcl_
pspher_
psrest_
psrot_
pssave_
pstext_
psubox_
psucir_
define_object_colours
define_text_colours
define_object_sizes
define_text_sizes
get_other_ring_bonds
psopen_
adjust_for_landscape
print_title
print_key_text
plot_spokes
plot_contact_group
accessibility_circle
print_symbol_key
extract_minmax
scale
centre_point
sort_accessibilities
initialise_atoms
shade_accessibilities
double_bond_split
get_atom_colour
get_subrings
plot_aromatics
plot_bond_lines
print_hbond_lengths
plot_hbond_lines
plot_ext_bond_lines
get_segment
position_names
print_current_atom
get_atom_text_colour
get_ring_atom_coords
print_atom_name
plot_atoms
convert_residue_name
print_name
plot_hydrophobics
plot_simple_ligand_residues
plot_simple_hgroups
extract_search_atoms
label_dist
perform_search
resname_grid_search
print_residue_names
plot_postscript_picture
unset_bond_plot_flags
plot_highlights
   -> extract_residue_name
         -> to_upper
         -> format_res_number
   -> create_highres_entry
   -> get_other_chains
         -> string_truncate
   -> get_variants
         -> string_truncate
         -> format_res_number
         -> find_highres_entry
         -> create_highres_entry
   -> calc_ellipse
         -> get_best_ellipse
               -> get_rot_matrix
               -> apply_atom_rotation
   -> find_highres_entry
   -> plot_ellipse
         -> pscomm_
         -> pssave_
         -> pslwid_
         -> pscolb_
         -> psdash_
         -> pselip_
         -> psrest_

*/


/***********************************************************************

string_truncate  -  Truncate the given string at the last non-blank
                    character, or at its maximum length

***********************************************************************/

int string_truncate(char *string,int max_length)
{
  char ch;

  int iend, ipos, len;

  /* Initialise variables */
  iend = -1;
  ch = string[0];
  for (ipos = 0; ipos < max_length && ch != '\0'; ipos++)
    {
      /* Get the current character and check for end of string */
      ch = string[ipos];
      if (ch == '\n' || ch == 13)
        ch = '\0';
      if (ch != '\0')
        {
          /* If not a space, then mark end of string */
          if (ch != ' ')
            iend = ipos;
        }
    }

  /* End the string after the last non-blank character */
  if (iend + 1 < max_length)
    {
      string[iend + 1] = '\0';
      len = iend + 1;
    }
  else
    {
      string[max_length - 1] = '\0';
      len = max_length - 1;
    }

  /* Return the string length */
  return(len);
}
/***********************************************************************

to_lower  -  Convert given character string to lower case

***********************************************************************/

void to_lower(char search_string[],int length)
{
  char ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'A';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'a';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

to_upper  -  Convert given character string to upper case

***********************************************************************/

void to_upper(char *search_string,int length)
{
  int ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'a';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'A';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

initialise_atoms  -  Initialise all the atom flags

***********************************************************************/

void initialise_atoms(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Initialise the flag and stack-pointer */
      atom_ptr->checked = FALSE;
      atom_ptr->next_stack_ptr = NULL;
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }               
}
/***********************************************************************

initialise_bonds  -  Initialise all the bond flags and stack pointers

***********************************************************************/

void initialise_bonds(void)
{
  struct bond *bond_ptr;

  /* Set pointer to first bond in linked-list */
  bond_ptr = first_bond_ptr;

  /* Loop through the stored bond-list to initialise flags and 
     stack-pointers */
  while (bond_ptr != NULL)
    {
      /* Initialise flag and stack-pointer */
      bond_ptr->checked = FALSE;
/* v.4.5.2--> */
      bond_ptr->ring_no = 0;
/* <--v.4.5.2 */
      bond_ptr->next_link_ptr = NULL;
      bond_ptr->next_stack_ptr = NULL;

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

centre_of_mass  -  Find the centre of mass of the supplied coordinates

***********************************************************************/

void centre_of_mass(float coord_store[MAXATOMS][3],int natoms,
                    float *centre_x,float *centre_y,float *centre_z)
{
  int iatom;
  float mass_x = 0.0, mass_y = 0.0, mass_z = 0.0, centre_mass_x,
  centre_mass_y, centre_mass_z;

  /* Loop through all the coordinates supplied */ 
  for (iatom = 0; iatom < natoms; iatom++)
    {
      mass_x = mass_x + coord_store[iatom][0];
      mass_y = mass_y + coord_store[iatom][1];
      mass_z = mass_z + coord_store[iatom][2];
    }

  /* Calculate centre of mass */
  centre_mass_x = mass_x / natoms;
  centre_mass_y = mass_y / natoms;
  centre_mass_z = mass_z / natoms;
  *centre_x = centre_mass_x;
  *centre_y = centre_mass_y;
  *centre_z = centre_mass_z;
}
/***********************************************************************

get_other_ring_bonds  -  Get all the other bonds in the current ring

***********************************************************************/

void get_other_ring_bonds(struct bond *test_bond_ptr)
{

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr;
  struct bond *next_stack_ptr, *next_link_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise pointers */
  next_link_ptr = NULL;
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Add the test-bond onto the start of the stack */
  bond_ptr = test_bond_ptr;
  next_stack_ptr = bond_ptr;
  next_free_stack_ptr = bond_ptr;
  bond_ptr->checked = TRUE;

  /* Initialise linked-list that will eventually contain all the
     bonds in the current ring */
  next_link_ptr = bond_ptr;

  /* Loop until stack of pointers to be processed has been exhausted
     or have determined that we have a ring */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
         pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Loop through all the bonds connected to the current one */

      /* Get the pointer to the first of the bonds coming off the 
         current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
        {
          /* Get the pointer to this bond */
          other_bond_ptr = bond_link_ptr->bond_ptr;

          /* Process only if this is a ring-bond or end-bond which
             hasn't already been processed */
          if (other_bond_ptr->checked == FALSE &&
              (other_bond_ptr->ring_bond == TRUE ||
               other_bond_ptr->end_bond == TRUE))
            {
              /* If this is a ring-bond, then add it to the stack */
              if (other_bond_ptr->ring_bond == TRUE)
                {
                  /* Add to end of stack of rings to be processed */
                  next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
                  next_free_stack_ptr = other_bond_ptr;

                  /* If the stack has run out, restart it at the bond
                     just entered */
                  if (next_stack_ptr == NULL)
                    next_stack_ptr = other_bond_ptr;
                }

              /* Add to the linked-list of bonds belonging to the
                 current ring (and end-bonds sprouting off it) */
              next_link_ptr->next_link_ptr = other_bond_ptr;
              next_link_ptr = other_bond_ptr;

              /* Mark this bond as checked */
              other_bond_ptr->checked = TRUE;
            }

          /* Go to the next bond-link in the linked list */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }
    }
}
/* v.4.5.2--> */
/***********************************************************************

create_path_record  -  Create a new path record

***********************************************************************/

struct path *create_path_record(struct path **fst_path_ptr,
                                struct path **lst_path_ptr,
                                struct bond *bond_ptr,int end)
{
  struct path *path_ptr;
  struct path *first_path_ptr, *last_path_ptr;

  /* Initialise path pointers */
  path_ptr = NULL;
  first_path_ptr = *fst_path_ptr;
  last_path_ptr = *lst_path_ptr;

  /* Create path entry */
  path_ptr = (struct path *) malloc(sizeof(struct path));
  if (path_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for struct path\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the bond pointer */
  path_ptr->bond_ptr = bond_ptr;
  path_ptr->end = end;

  /* Initialise pointer to next entry */
  path_ptr->next_path_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head
     of linked list */
  if (first_path_ptr == NULL)
    first_path_ptr = path_ptr;

  /* Otherwise, make previous entry point to the current one */
  else
    last_path_ptr->next_path_ptr = path_ptr;
                      
  /* Save the current residue's pointer */
  last_path_ptr = path_ptr;

  /* Return current pointers */
  *fst_path_ptr = first_path_ptr;
  *lst_path_ptr = last_path_ptr;
  return(path_ptr);
}
/***********************************************************************

get_subrings  -  Determine how many subrings there are in the given
                 ring system

***********************************************************************/

int get_subrings(struct bond *ring_bond_ptr)
{
  int ibond, ipath, jpath, npaths, nrings, pathlen[MAXPATHS];
  int end, other_end;
  int added_to, first;

  struct bond *bond_ptr, *link_bond_ptr, *other_bond_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *end_atom_ptr;
  struct path *first_path_ptr[MAXPATHS], *last_path_ptr[MAXPATHS];
  struct path *path_ptr;

  /* Initialise variables */
  nrings = 0;

  /* Get pointer to the first bond in the current ring */
  bond_ptr = ring_bond_ptr;

  /* Loop through all the bonds in the current ring system */
  while (bond_ptr != NULL)
    {
      /* If this is an aromatic ring-bond which hasn't yet been assigned
         to an aromatic ring, see if we can find a set of 6 or fewer
         aromatic bonds linking back to it */
      if (bond_ptr->ring_bond == TRUE &&
          bond_ptr->bond_order == AROMATIC && bond_ptr->ring_no == 0)
        {
          /* Store the end atom that we hope to reach */
          end_atom_ptr = bond_ptr->second_atom_ptr;

          /* Initialise paths from this bond */
          for (ipath = 0; ipath < MAXPATHS; ipath++)
            {
              first_path_ptr[ipath] = NULL;
              last_path_ptr[ipath] = NULL;
              pathlen[ipath] = 0;
            }
          npaths = 0;

          /* Start path at this bond */
          path_ptr = create_path_record(&first_path_ptr[npaths],
                                        &last_path_ptr[npaths],bond_ptr,
                                        FIRST);

          /* Increment path length */
          pathlen[npaths]++;

          /* Increment count of paths */
          npaths++;

          /* Initialise path to be checked */
          added_to = TRUE;
          ipath = 0;

          /* Process the paths until none left or have assigned current
             bond to a ring */
          while (npaths > 0 && bond_ptr->ring_no == 0 && added_to == TRUE)
            {
              /* Initialise flag */
              added_to = FALSE;

              /* If path length already too long to be 6-sided ring,
                 then delete path and free up its space */
              if (pathlen[ipath] > MAXPATH_LEN)
                {
                  /* Shift all subsequent pointers and path lengths up a
                     slot */
                  for (jpath = ipath; jpath < npaths - 1; jpath++)
                    {
                      /* Shift all details up a slot */
                      first_path_ptr[jpath] = first_path_ptr[jpath + 1];
                      last_path_ptr[jpath] = last_path_ptr[jpath + 1];
                      pathlen[jpath] = pathlen[jpath + 1];
                    }

                  /* Reinitialise last slot, now unused */
                  first_path_ptr[npaths - 1] = NULL;
                  last_path_ptr[npaths - 1] = NULL;
                  pathlen[npaths - 1] = 0;

                  /* Decrease count of paths left */
                  npaths--;
                }

              /* If we are off the end of the array, reset to the
                 beginning */
              if (ipath > npaths - 1)
                ipath = 0;

              /* Get this path's last bond */
              if (npaths > 0)
                {
                  other_bond_ptr = last_path_ptr[ipath]->bond_ptr;
                  end = last_path_ptr[ipath]->end;
                }
              else
                other_bond_ptr = NULL;

              /* If have a valid bond and path, then extend path */
              if (other_bond_ptr != NULL && pathlen[ipath] > 0)
                {
                  /* Get the first of this bond's links */
                  bond_link_ptr = other_bond_ptr->first_bond_link_ptr;

                  /* Initialise flag */
                  first = TRUE;

                  /* Work through all the bonds coming off one end */
                  while (bond_link_ptr != NULL)
                    {
                      /* Get the linked bond */
                      link_bond_ptr = bond_link_ptr->bond_ptr;

                      /* If this is our starting bond, check to see if
                         we have a closed ring, so can assign all bonds
                         in this path the same ring number */
                      if (link_bond_ptr == bond_ptr &&
                          (other_bond_ptr->first_atom_ptr == end_atom_ptr ||
                           other_bond_ptr->second_atom_ptr == end_atom_ptr))
                        {
                          /* Increment count of rings encountered */
                          nrings++;

                          /* If not the first offspring bond, then need
                             to stop short of previous bond that will have
                             been added to this path */
                          if (first == FALSE)
                            pathlen[ipath]--;

                          /* Get the first of this path's bonds */
                          path_ptr = first_path_ptr[ipath];
                          ibond = 0;

                          /* Loop through all the bonds in the current
                             path to assign a ring number */
                          while (path_ptr != NULL && ibond < pathlen[ipath])
                            {
                              /* Get this bond */
                              link_bond_ptr = path_ptr->bond_ptr;

                              /* If already assigned to a ring,
                                 add multiples of 10 of current
                                 ring number */
                              if (link_bond_ptr->ring_no != 0)
                                link_bond_ptr->ring_no
                                  = link_bond_ptr->ring_no + 10 * nrings;
                              else
                                link_bond_ptr->ring_no = nrings;

                              /* Get the next path pointer */
                              path_ptr = path_ptr->next_path_ptr;
                              ibond++;
                            }
                        }

                      /* Otherwise, if this is an aromatic ring bond,
                         add to current path */
                      else if (link_bond_ptr->ring_bond == TRUE &&
                               link_bond_ptr->bond_order == AROMATIC &&
                               bond_link_ptr->bond_end == end)
                        {
                          /* Determine which end is being attached */
                          if ((link_bond_ptr->first_atom_ptr ==
                               other_bond_ptr->first_atom_ptr) ||
                              (link_bond_ptr->first_atom_ptr ==
                               other_bond_ptr->second_atom_ptr))
                            other_end = SECOND;
                          else
                            other_end = FIRST;

                          /* Determine which of the current bond's
                             end this atom belongs to */

                          /* If this is the first such bond found,
                             add to current path */
                          if (first == TRUE)
                            {
                              jpath = ipath;
                              first = FALSE;
                            }

                          /* Otherwise, need to make a copy of the
                             current path */
                          else if (npaths < MAXPATHS)
                            {
                              /* Initialise length of new path */
                              pathlen[npaths] = 0;

                              /* Get the first path record */
                              path_ptr = first_path_ptr[ipath];

                              /* Loop through all bonds in the path
                                 excluding the last */
                              while (path_ptr != NULL &&
                                     path_ptr != last_path_ptr[ipath])
                                {
                                  /* Create a new path record to store
                                     the copy of the current bond */
                                  create_path_record(&first_path_ptr[npaths],
                                                     &last_path_ptr[npaths],
                                                     path_ptr->bond_ptr,
                                                     path_ptr->end);

                                  /* Increment path length */
                                  pathlen[npaths]++;

                                  /* Set flag */
                                  added_to= TRUE;

                                  /* Get the next path pointer */
                                  path_ptr = path_ptr->next_path_ptr;
                                }

                              /* Save identifier of new path */
                              jpath = npaths;

                              /* Increment count of paths */
                              npaths++;
                            }

                          /* If can't create a new path, show warning
                             message */
                          else
                            {
                              printf("*** Warning. Maximum no. of paths "
                                     "exceeded: %d\n",MAXPATHS);
                              jpath = -1;
                            }

                          /* Create a new path record for the new bond */
                          if (jpath > -1)
                            {
                              path_ptr
                                = create_path_record(&first_path_ptr[jpath],
                                                     &last_path_ptr[jpath],
                                                     link_bond_ptr,
                                                     other_end);

                              /* Increment path length */
                              pathlen[jpath]++;

                              /* Set flag */
                              added_to= TRUE;
                            }
                        }

                      /* Go to the next bond-link in the linked list */
                      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
                    }

                  /* Increment path being considered */
                  ipath++;
                  if (ipath >= npaths)
                    ipath = 0;
                }
            }
        }

      /* Get the next bond in the current ring */
      bond_ptr = bond_ptr->next_link_ptr;
    }

  /* Return number of subrings found */
  return(nrings);
}
/* <--v.4.5.2 */
/***********************************************************************

get_ring_atom_coords  -  Extract all the coordinates of the current
                         ring's atoms, so that can squash to make planar

***********************************************************************/

void get_ring_atom_coords(struct bond *ring_bond_ptr,
                          float coord_store[MAXATOMS][3],
/* v.4.4--> */
/*                          int *nring_atoms,int *nring_bonds) */
/* v.4.5.2--> */
                          // int *nring_atoms,int *nring_bonds,
                          int *nring_atoms,int *nring_bonds,int ring,
/* <--v.4.5.2 */
                          int aromatics_only)
/* <--v.4.4 */
{
/* v.4.5.2--> */
  // int loop, natoms, nbonds;
  int factor, i, iring, loop, natoms, nbonds, ring_no;
  int wanted;
/* <--v.4.5.2 */

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  natoms = 0;
  nbonds = 0;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Get pointer to the first bond in the current ring */
  bond_ptr = ring_bond_ptr;

  /* Loop through all the bonds in the current ring */
  while (bond_ptr != NULL)
    {
/* v.4.4--> */
      /* Only want coordinates if this is a ring-bond (as opposed
         to an associated end-bond) */
      if (bond_ptr->ring_bond == TRUE)
        {

          /* Check this bond has the correct ring number */
          wanted = FALSE;
          iring = bond_ptr->ring_no;
          factor = 100;
          for (i = 0; i < 3 && wanted == FALSE; i++)
            {
              /* Divide by current factor if ring-number high enough */
              if (iring > factor - 1)
                {
                  ring_no = iring / factor;
                  if (ring_no == ring)
                    wanted = TRUE;
                  iring = ring_no % factor;
                }

              /* Reduce factor by factor of 10 */
              factor = factor / 10;
            }

          /* If only interested in aromatics, make sure this bond is one */
          if (aromatics_only == FALSE ||
              (bond_ptr->bond_order == AROMATIC && wanted == TRUE))
            {
/* <--v.4.4 */
              /* Increment count of bonds in this ring */
              nbonds++;

              /* Loop over this bond's two atoms */
              for (loop = 0; loop < 2; loop++)
                {
                  if (loop == 0)
                    atom_ptr = bond_ptr->first_atom_ptr;
                  else
                    atom_ptr = bond_ptr->second_atom_ptr;

                  /* If this atom's coordinates have not already been stored,
                     then add to list for flattening */
                  if (atom_ptr->checked == FALSE)
                    {
                      /* Store the atom's coordinates */

                      /* Check that maximum storage not exceeded */
                      if (natoms > MAXSTACK)
                        {
                          printf("*** Maximum number of atoms");
                          printf(" exceeded in flatten_ring");
                          printf("  %d\n",MAXSTACK);
                          exit(1);
                        }
                      coord_store[natoms][0] = atom_ptr->x;
                      coord_store[natoms][1] = atom_ptr->y;
                      coord_store[natoms][2] = atom_ptr->z;
                      natoms++;

                      /* Mark atom as already included */
                      atom_ptr->checked = TRUE;
                    }
                }
/* v.4.4--> */
              /* If plotting aromatic ring bonds, mark this bond as
                 processed */
              if (aromatics_only == TRUE)
                bond_ptr->plotted = TRUE;
            }
/* <--v.4.4 */
        }

      /* Mark this bond as already flattened */
      bond_ptr->flattened = TRUE;

      /* Get the next bond in the current ring */
      bond_ptr = bond_ptr->next_link_ptr;
    }

  /* Return the number of atoms stored */
  *nring_atoms = natoms;
  *nring_bonds = nbonds;
}


/**************   P O S T S C R I P T   R O U T I N E S   *************/

/***********************************************************************

psarc_  -  Write arc out to PostScript file

***********************************************************************/

void psarc_(float x1,float y1,float radius,float angle_start,float angle_end)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f Arc\n",x1, y1, radius,
          angle_start, angle_end);
}
/***********************************************************************

psbbox_  -  Function to write out bounded box to Postscript file

**********************************************************************/

void psbbox_(float x1,float y1,float x2,float y2,float x3,float y3,
             float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pline4\n",
          x1,y1,x2,y2,x3,y3,x4,y4);
}
/***********************************************************************

psccol_  -  Write circle colour out to PostScript file

***********************************************************************/

void psccol_(char *colour)
{
  fprintf(ps_file,"/Circol { %s } def\n",colour);
}
/***********************************************************************

pscirc_  -  Write a circle out to PostScript file

***********************************************************************/

void pscirc_(float x,float y,char *radius)
{
  fprintf(ps_file,"%7.2f %7.2f %s Circle\n",x,y,radius);
}
/***********************************************************************

psclos_  -  Write final lines to PostScript file and close

***********************************************************************/

/* v.4.4.4--> */
/* void psclos_(void) */
void psclos_(int *bounds)
/* <--v.4.4.4 */
{

    /* Write out closing lines to PostScript file*/
    fprintf(ps_file,"\n\n");
    fprintf(ps_file,"LigplotSave restore\n");
    fprintf(ps_file,"showpage\n");
    fprintf(ps_file,"%%%%Trailer\n");
/* v.4.4.4--> */
/*    fprintf(ps_file,"%%%%BoundingBox: %d %d %d %d\n",(int) BBOXX1 - 1,
      (int) BBOXY1 - 1, (int) BBOXX2 + 1, (int) BBOXY2 + 1); */
    fprintf(ps_file,"%%%%BoundingBox: %d %d %d %d\n",bounds[0],
            bounds[1],bounds[2],bounds[3]);
/* <--v.4.4.4 */
    fprintf(ps_file,"%%%%EOF\n");

/* v.4.3--> */
    /* Close the PostScript file */
    fclose(ps_file);
/* <--v.4.3 */
}
/***********************************************************************

pscolb_  -  Write background colour level (for text, etc)

***********************************************************************/

void pscolb_(char *text_string)
{
    fprintf(ps_file,"%s\n",text_string);
}
/***********************************************************************

pscomm_  -  Write out a comment to the PostScript file

***********************************************************************/

void pscomm_(char *text_string)
{
  fprintf(ps_file,"\n");
  fprintf(ps_file,"%% %s\n",text_string);
  fprintf(ps_file,"\n");
}
/***********************************************************************

pscrgb_  -  Write RGB circle colour out to PostScript file

***********************************************************************/

void pscrgb_(float red, float green, float blue)
{
  fprintf(ps_file,"/Circol { %7.4f %7.4f %7.4f setrgbcolor } def\n",
          red, green, blue);
}
/***********************************************************************

psctxt_  -  Write out centred text to PostScript file

***********************************************************************/

void psctxt_(float x,float y,char *size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %s Center\n",text_string,size);
  fprintf(ps_file,"(%s) %s Print\n",text_string,size);
}
/***********************************************************************

psdash_  -  Switch dashed lines on/off

***********************************************************************/

void psdash_(int ONOFF)
{
  if (ONOFF == 0)
    fprintf(ps_file,"[]  0 setdash\n");
  else
    fprintf(ps_file,"[ %2d %2d ] 0 setdash\n",ONOFF,ONOFF);
}
/* v.5.3.8--> */
/***********************************************************************

pselip_  -  Write out a coloured ellipse to PostScript file

***********************************************************************/

void pselip_(float x,float y,float width,float ratio,float angle)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n",x,y);
  fprintf(ps_file,"%7.2f RotAngle\n",-angle);
  fprintf(ps_file,"0.0 0.0 %7.2f 1.0 %7.2f 0 0 Oellipse\n",width,ratio);
  fprintf(ps_file,"grestore\n");
}
/* <--v.5.3.8 */
/* v.4.4.1--> */
/***********************************************************************

psfont_  -  Write out text in appropriate font

***********************************************************************/

void psfont_(char *size,char *text_string,float text_height,int greek,
             int half)
{
  char size_string[PS_STRING_LENGTH];

  /* Copy text size */
  strcpy(size_string,size);

  /* If text needs to be half-height, then adjust size */
  if (half == TRUE)
    {
      /* Half the text height */
      text_height = text_height / 2;
      sprintf(size_string,"%7.2f",text_height);
    }

  /* Print text in appropriate font */
  if (greek == TRUE)
    fprintf(ps_file,"(%s) Gfont %s Fprint\n",text_string,size_string);
  else
    fprintf(ps_file,"(%s) Font %s Fprint\n",text_string,size_string);
}
/* <--v.4.4.1 */
/***********************************************************************

pshade_  -  Write colour/shade out to PostScript file

***********************************************************************/

void pshade_(char *colour)
{
    fprintf(ps_file,"G %s\n",colour);
}
/***********************************************************************

psline_  -  Write line out to PostScript file

***********************************************************************/

void psline_(float x1,float y1,float x2,float y2)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f L\n",x1, y1, x2, y2);
}
/***********************************************************************

pslwid_  -  Write line-width out to PostScript file

***********************************************************************/

void pslwid_(char *width)
{
    fprintf(ps_file,"%s W\n",width);
}
/* v.4.4.1--> */
/***********************************************************************

psmove_  -  Move to start of text

***********************************************************************/

void psmove_(float x,float y,char *size,char *text_string)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n",x,y);
  fprintf(ps_file,"(%s) %s Center\n",text_string,size);
}
/* <--v.4.4.1 */
/***********************************************************************

psphcl_  -  Write sphere colour out to PostScript file

***********************************************************************/

void psphcl_(char *colour)
{
    fprintf(ps_file,"/Sphcol { %s } def\n",colour);
}
/***********************************************************************

pspher_  -  Write a sphere out to PostScript file

***********************************************************************/

void pspher_(float x,float y,char *radius)
{
  fprintf(ps_file,"%7.2f %7.2f %s Sphere\n",x,y,radius);
}
/***********************************************************************

psrest_  -  Restore graphics state

***********************************************************************/

void psrest_(void)
{
    fprintf(ps_file,"grestore\n");
}
/***********************************************************************

psrot_  -  Rotate graphics page through 90 degrees

***********************************************************************/

void psrot_(float new_origin_x, float new_origin_y)
{
    fprintf(ps_file," %7.2f %7.2f moveto Rot90\n",new_origin_x,new_origin_y);
}
/***********************************************************************

pssave_  -  Save graphics state

***********************************************************************/

void pssave_(void)
{
    fprintf(ps_file,"gsave\n");
}
/***********************************************************************

pstext_  -  Write out text to PostScript file

***********************************************************************/

void pstext_(float x,float y,char *size,char *text)
{
  fprintf(ps_file,"%7.2f %7.2f moveto\n", x, y);
  fprintf(ps_file,"(%s) %s Print\n",text,size);
}
/***********************************************************************

psubox_  -  Write out unbounded box to PostScript file

***********************************************************************/

void psubox_(float x1,float y1,float x2,float y2,float x3,float y3,
        float x4,float y4)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f Pl4\n",
          x1, y1, x2, y2, x3, y3, x4, y4);
}
/***********************************************************************

psucir_  -  Write out a grey/coloured circle, with no border, to
            PostScript file

***********************************************************************/

void psucir_(float x,float y,float radius)
{
  fprintf(ps_file,"%7.2f %7.2f %7.2f Ucircle\n",x,y,radius);
}
/***********************************************************************

match_colour_names  -  Match up the colour names with the Colour
                       Names Table

***********************************************************************/

void match_colour_names(void)
{
  int icolour, iletter, jcolour, last_letter, match;

  /* Check that colour name can be used, adjusting where necessary */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
    {
      last_letter = 0;

      /* Loop through all the letters of the name */
      for (iletter = 0; iletter < COL_NAME_LEN; iletter++)
        {
          /* Replace blanks by underscores */
          if (Colour_Table_Name[icolour][iletter] == ' ' ||
              Colour_Table_Name[icolour][iletter] == '\'')
            Colour_Table_Name[icolour][iletter] = '_';
          else
            last_letter = iletter;
        }

      /* Blank out any trailing spaces */
      for (iletter = last_letter + 1; iletter < COL_NAME_LEN; iletter++)
        Colour_Table_Name[icolour][iletter] = ' ';

      /* Check for duplicate colours, discarding second instance */
      for (jcolour = 0; jcolour < icolour; jcolour++)
        {
          if (!strncmp(Colour_Table_Name[jcolour],
                       Colour_Table_Name[icolour],COL_NAME_LEN))
/* v.4.2--> */
/*            Colour_Table_Name[jcolour][0] = '\0'; */
            Colour_Table_Name[icolour][0] = '\0';
/* <--v.4.2 */
        }
    }

  /* Check for atom-colour option on ligand and non-ligand bonds */
  if (!strncmp(object_colour[1],"ATOM           ",15))
    Split_Colour_Ligand_Bonds = TRUE;
  if (!strncmp(object_colour[2],"ATOM           ",15))
    Split_Colour_Nonligand_Bonds = TRUE;
  
  /* Adjust the object colours */
  for (icolour = 0; Picture->In_Colour == TRUE && icolour < OBJECT_COLOURS;
       icolour++)
    {
      last_letter = 0;
      
      /* Loop through all the letters of the name */
      for (iletter = 0; iletter < COL_NAME_LEN; iletter++)
        {
          /* Replace blanks by underscores */
          if (object_colour[icolour][iletter] == ' ')
            object_colour[icolour][iletter] = '_';
          else
            last_letter = iletter;
        }

      /* Blank out any trailing spaces */
      for (iletter = last_letter + 1; iletter < COL_NAME_LEN; iletter++)
        object_colour[icolour][iletter] = ' ';
      
      /* If the colour of ligand or non-ligand bonds has been defined
         as ATOM colour, then don't need to match against the colour
         table */
      match = FALSE;
      if ((icolour == 1 && Split_Colour_Ligand_Bonds == TRUE) ||
          (icolour == 2 && Split_Colour_Nonligand_Bonds == TRUE))
        match = TRUE;

      /* Check that this colour has a matching entry in the table
         of colour names */
      for (jcolour = 0; jcolour < MAX_COLOURS && match == FALSE; jcolour++)
        {
          if (!strncmp(object_colour[icolour],
                       Colour_Table_Name[jcolour],COL_NAME_LEN))
            match = TRUE;
        }

      /* If failed to find a match, use whichever is the first colour */
      if (match == FALSE)
        {
          printf("*** Warning. Unidentified colour [%s]. Replaced by [%s]\n",
                 object_colour[icolour],Colour_Table_Name[0]);
          strcpy(object_colour[icolour],Colour_Table_Name[0]);
          Nwarnings++;
        }
    }

  /* Adjust the text colours */
  for (icolour = 0; Picture->In_Colour == TRUE && icolour < TEXT_COLOURS;
       icolour++)
    {
      last_letter = 0;
      
      /* Loop through all the letters of the name */
      for (iletter = 0; iletter < COL_NAME_LEN; iletter++)
        {
          /* Replace blanks by underscores */
          if (text_colour[icolour][iletter] == ' ')
            text_colour[icolour][iletter] = '_';
          else
            last_letter = iletter;
        }

      /* Blank out any trailing spaces */
      for (iletter = last_letter + 1; iletter < COL_NAME_LEN; iletter++)
        text_colour[icolour][iletter] = ' ';

      /* Check that this colour has a matching entry in the table
         of colour names */
      match = FALSE;
      for (jcolour = 0; jcolour < MAX_COLOURS; jcolour++)
        {
          if (!strncmp(text_colour[icolour],
                       Colour_Table_Name[jcolour],COL_NAME_LEN))
            match = TRUE;
        }

      /* If failed to find a match, use whichever is the first colour */
      if (match  == FALSE)
        {
          printf("*** Warning. Unidentified colour [%s]. Replaced by [%s]\n",
                 text_colour[icolour],Colour_Table_Name[0]);
          strcpy(text_colour[icolour],Colour_Table_Name[0]);
          Nwarnings++;
        }
    }
}
/***********************************************************************

define_object_colours  -  Write the object colours and definitions
                          to the PostScript file

***********************************************************************/

void define_object_colours(void)
{
    char colour[COL_NAME_LEN + 1];
    int icolour, jcolour;

    /* Write out the object colours */
    for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
    {
        /* Retrieve the appropriate description */
        switch (icolour)
        {
            case 0 : strcpy(colour,Colour->Background); break;
            case 1 : strcpy(colour,Colour->Ligand_Bonds); break;
            case 2 : strcpy(colour,Colour->Nonligand_Bonds); break;
            case 3 : strcpy(colour,Colour->Hydrogen_Bonds); break;
            case 4 : strcpy(colour,Colour->External_Bonds); break;
            case 5 : strcpy(colour,Colour->Hydrophobics); break;
            case 6 :
                strcpy(colour,Colour->Accessibility_Min);

                /* Save the RGB values for the minimum accessibility colour */
                for (jcolour = 0; jcolour < MAX_COLOURS; jcolour++)
                {
                    if (!strncmp(object_colour[icolour],
                        Colour_Table_Name[jcolour],COL_NAME_LEN))
                    {
                        Accessibility_Min[0] = Colour_Table[jcolour][0];
                        Accessibility_Min[1] = Colour_Table[jcolour][1];
                        Accessibility_Min[2] = Colour_Table[jcolour][2];
                    }
                }
                break;
            case 7 :
                strcpy(colour,Colour->Accessibility_Max);

                /* Save the RGB values for the maximum accessibility colour */
                for (jcolour = 0; jcolour < MAX_COLOURS; jcolour++)
                {
                    if (!strncmp(object_colour[icolour],
                        Colour_Table_Name[jcolour],COL_NAME_LEN))
                    {
                        Accessibility_Max[0] = Colour_Table[jcolour][0];
                        Accessibility_Max[1] = Colour_Table[jcolour][1];
                        Accessibility_Max[2] = Colour_Table[jcolour][2];
                    }
                }
                break;
            case 8 : strcpy(colour,Colour->Nitrogen); break;
            case 9 : strcpy(colour,Colour->Oxygen); break;
            case 10 : strcpy(colour,Colour->Carbon); break;
            case 11 : strcpy(colour,Colour->Sulphur); break;
            case 12 : strcpy(colour,Colour->Water); break;
            case 13 : strcpy(colour,Colour->Phosphorus); break;
            case 14 : strcpy(colour,Colour->Iron); break;
            case 15 : strcpy(colour,Colour->Other); break;
            case 16 : strcpy(colour,Colour->Atom_Edges); break;
            case 17 : strcpy(colour,Colour->Simple_Residues); break;
            default :
              printf("*** Program error. Invalid Colours parameter: %d\n",
                     icolour);
              exit(1);
        }

        /* Write out the colour and its definition to the PostScript file */
        if (!strncmp(object_colour[icolour],"ATOM           ",15))
          fprintf(ps_file,"%% %s defined as ATOM colour\n",colour);
        else
          fprintf(ps_file,"/%s { %s } def\n",colour,object_colour[icolour]);
    }
}
/***********************************************************************

define_text_colours  -  Write the text colours and definitions
                        to the PostScript file

***********************************************************************/

void define_text_colours(void)
{
  char colour[COL_NAME_LEN + 1];
  int icolour;

  /* Write out the text colours */
  for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
    {
      /* Retrieve the appropriate description */
      switch (icolour)
        {
        case 0 : strcpy(colour,Text_Colour->Title); break;
        case 1 : strcpy(colour,Text_Colour->Key_Text); break;
        case 2 : strcpy(colour,Text_Colour->Ligand_Residue_Names); break;
        case 3 : strcpy(colour,Text_Colour->Nonligand_Residue_Names); break;
        case 4 : strcpy(colour,Text_Colour->Water_Names); break;
        case 5 : strcpy(colour,Text_Colour->Hydrophobic_Names); break;
        case 6 : strcpy(colour,Text_Colour->Ligand_Atom_Names); break;
        case 7 : strcpy(colour,Text_Colour->Nonligand_Atom_Names); break;
        case 8 : strcpy(colour,Text_Colour->Hbond_Lengths); break;
        default :
          printf("*** Program error. Invalid Text Colours parameter: %d\n",
                 icolour);
          exit(1);
        }

      /* Write out the colour and its definition to the PostScript file */
      fprintf(ps_file,"/%s { %s } def\n",colour,text_colour[icolour]);
    }
}
/***********************************************************************

define_object_sizes  -  Write the object size definitions to the
                        PostScript file

***********************************************************************/

void define_object_sizes(float scale_factor,int print_out)
{
  char object[COL_NAME_LEN + 1];
  int isize;
  float simple_thick;

  /* Initialise maximum object sizes */
  Max_object_size = 0.0;

  /* Write out the object sizes */
  for (isize = 0; isize < OBJECT_SIZES; isize++)
    {
      /* Use scaling factor to change from Angstroms to PostScript
         sizes */
/* v.4.2--> */
/*      object_size[isize] = scale_factor * object_size[isize]; */
      object_size[isize] = scale_factor * param_object_size[isize];
/* <--v.4.2 */
      
      /* Retrieve the appropriate description */
      switch (isize)
        {
        case 0 : 
          strcpy(object,Size->Ligand_Atoms);
          Size_Val->Ligand_Atoms = object_size[isize];
          if (Include->Ligand_Atoms == TRUE &&
              object_size[isize] > Max_object_size)
            Max_object_size = object_size[isize];
          break;
        case 1 : 
          strcpy(object,Size->Nonligand_Atoms);
          Size_Val->Nonligand_Atoms = object_size[isize];
          if (Include->Nonligand_Atoms == TRUE &&
              object_size[isize] > Max_object_size)
            Max_object_size = object_size[isize];
          break;
        case 2 :
          strcpy(object,Size->Waters);
          Size_Val->Waters = object_size[isize];
          if (Include->Waters == TRUE &&
              object_size[isize] > Max_object_size)
            Max_object_size = object_size[isize];
          break;
        case 3 : 
          strcpy(object,Size->Hydrophobics);
          Size_Val->Hydrophobics = object_size[isize];
          if (Include->Hydrophobics == TRUE &&
              object_size[isize] > Max_object_size)
            Max_object_size = SPOKE_EXTENT * object_size[isize];
          break;
        case 4 :
          strcpy(object,Size->Simple_Residues);
          Size_Val->Simple_Residues = object_size[isize];
          if (Include->Simple_Ligand_Residues == TRUE &&
              object_size[isize] > Max_object_size)
            Max_object_size = object_size[isize];
          simple_thick = object_size[isize] / 10.0;
          break;
        case 5 : 
          strcpy(object,Size->Ligand_Bonds);
          Size_Val->Ligand_Bonds = object_size[isize];
          break;
        case 6 : 
          strcpy(object,Size->Nonligand_Bonds);
          Size_Val->Nonligand_Bonds = object_size[isize];
          break;
        case 7 : 
          strcpy(object,Size->Hydrogen_Bonds);
          Size_Val->Hydrogen_Bonds = object_size[isize];
          break;
        case 8 : 
          strcpy(object,Size->External_Bonds);
          Size_Val->External_Bonds = object_size[isize];
          break;
        default :
          printf("*** Program error. Invalid Size parameter: %d\n",
                 isize);
          exit(1);
        }

      /* Write out the size definition and its value to the PostScript
         file */
      if (print_out == TRUE)
        fprintf(ps_file,"/%s { %8.3f } def\n",object,object_size[isize]);

/* v.4.1--> */
      /* Save the description associated with the object size */
      strcpy(object_size_desc[isize],object);
/* <--v.4.1 */
    }
  /* Write out corresponding line-thicknesses */
  if (print_out == TRUE)
    {
      fprintf(ps_file,"/%s { %8.3f } def\n","Default_linewidth",0.2);
      fprintf(ps_file,"/%s { %8.3f } def\n","Arc_linewidth    ",0.75);
      fprintf(ps_file,"/%s { %8.3f } def\n","Simple_width     ",
              simple_thick);
    }
}
/***********************************************************************

define_text_sizes  -  Write the text size definitions to the
                      PostScript file

***********************************************************************/

void define_text_sizes(float scale_factor,int print_out)
{
  char text[COL_NAME_LEN + 1];
  int isize;

  /* Initialise maximum text height and width */
  Max_Textsize = 0.0;
  Max_Textwidth = 0.0;

  /* Write out the textt sizes */
  for (isize = 0; isize < TEXT_SIZES; isize++)
    {
      /* Use scaling factor to change from Angstroms to PostScript
         sizes */
/* v.4.2--> */
/*      text_size[isize] = scale_factor * text_size[isize]; */
      text_size[isize] = scale_factor * param_text_size[isize];
/* <--v.4.2 */
        
      /* Retrieve the appropriate description */
      switch (isize)
        {
        case 0 :
          strcpy(text,Text_Size->Ligand_Residue_Names);
          Text_Size_Val->Ligand_Residue_Names = text_size[isize];
          if (Include->Residue_Names == TRUE &&
              text_size[isize] > Max_Textsize)
            {
              Max_Textsize = text_size[isize];
              Max_Textwidth = 10 * Max_Textsize * CHAR_ASPECT;
            }
          break;
        case 1 :
          strcpy(text,Text_Size->Nonligand_Residue_Names);
          Text_Size_Val->Nonligand_Residue_Names = text_size[isize];
          if (Include->Residue_Names == TRUE &&
              text_size[isize] > Max_Textsize)
            {
              Max_Textsize = text_size[isize];
              Max_Textwidth = 10 * Max_Textsize * CHAR_ASPECT;
            }
          break;
        case 2 : 
          strcpy(text,Text_Size->Water_Names);
          Text_Size_Val->Water_Names = text_size[isize];
          if (Include->Residue_Names == TRUE &&
              text_size[isize] > Max_Textsize)
            {
              Max_Textsize = text_size[isize];
              Max_Textwidth = 10 * Max_Textsize * CHAR_ASPECT;
            }
          break;
        case 3 :
          strcpy(text,Text_Size->Hydrophobic_Names);
          Text_Size_Val->Hydrophobic_Names = text_size[isize];
          break;
        case 4 : 
          strcpy(text,Text_Size->Simple_Residue_Names); 
          Text_Size_Val->Simple_Residue_Names = text_size[isize];
          break;
        case 5 : 
          strcpy(text,Text_Size->Ligand_Atom_Names);
          Text_Size_Val->Ligand_Atom_Names = text_size[isize];
          if (Include->Atom_Names == TRUE &&
              text_size[isize] > Max_Textsize)
            {
              Max_Textsize = text_size[isize];
              Max_Textwidth = 4 * Max_Textsize * CHAR_ASPECT;
            }
          break;
        case 6 : 
          strcpy(text,Text_Size->Nonligand_Atom_Names);
          Text_Size_Val->Nonligand_Atom_Names = text_size[isize];
          if (Include->Residue_Names == TRUE &&
              text_size[isize] > Max_Textsize)
            {
              Max_Textsize = text_size[isize];
              Max_Textwidth = 4 * Max_Textsize * CHAR_ASPECT;
            }
          break;
        case 7 : 
          strcpy(text,Text_Size->Hbond_Lengths); 
          Text_Size_Val->Hbond_Lengths = text_size[isize];
          break;
        default :
          printf("*** Program error. Invalid Text Size parameter: %d\n",
                 isize);
          exit(1);
        }

      /* Write out the size definition and its value to the PostScript
         file */
      if (print_out == TRUE)
        fprintf(ps_file,"/%s { %8.3f } def\n",text,text_size[isize]);

/* v.4.1--> */
      /* Save the description associated with the text size */
      strcpy(text_size_desc[isize],text);
/* <--v.4.1 */
    }
}
/***********************************************************************

psopen_  -  Open PostScript file and write out header records

***********************************************************************/

void psopen_(void)
{
  char  colno[3], exclamation, percent;
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */

  float grey_shade, xx1, xx2, xy1, xy2;

  int   icolour, igrey;
/* v.5.2 --> */
  int i;
/* <-- v.5.2 */

  /* Initialise variables */
  percent = '%';
  exclamation = '!';

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

  /* Print the plot title */
  printf("\nWriting to PostScript output file: ligplot.ps\n\n");
  if (Print_Title[0] != '\0')
      printf("Title: %s\n\n",Print_Title);

  /* Open the PostScript file */
/* v.5.0--> */
  /* if ((ps_file = fopen("ligplot.ps","w")) == NULL) */
  filename[0] = '\0';
/* v.5.3.8--> */
  // if (wkdir[0] != '\0')
  if (wkdir[0] != '\0' && highlight_res[0] == '\0')
/* <--v.5.3.8 */
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.ps");
  if ((ps_file = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
      printf("\n*** Unable to open PostScript file: ligplot.ps\n");
      exit(1);
    }

  /* Check and adjust colour names in colour definitions */
/* v.4.2-->
  match_colour_names();
<--v.4.2 */

  /* Write out headings to PostScript file*/
  fprintf(ps_file,"%c%cPS-Adobe-3.0\n",percent,exclamation);
/* v.4.1--> */
/*  fprintf(ps_file,"%%%%Creator: Ligplot v.3.0\n"); */
  fprintf(ps_file,"%%%%Creator: Ligplot v.4.1\n");
/* <--v.4.1 */
  fprintf(ps_file,"%%%%DocumentNeededResources: font Times-Roman Symbol\n");
  fprintf(ps_file,"%%%%BoundingBox: (atend)\n");
  fprintf(ps_file,"%%%%Pages: 1\n");
  if (Print_Title[0] != '\0')
      fprintf(ps_file,"%%%%Title: %s\n",Print_Title);
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

  /* Loop to write out the different grey-shades */
  for (igrey = 0; igrey < 11; igrey++)
  {
      grey_shade = (float)igrey / 10.0;
      sprintf(colno,"%2d",igrey);
      colno[2] = '\0';

      if (colno[0] == ' ')
          colno[0] = '0';
      fprintf(ps_file,"/Grey%s { %6.2f setgray } def\n",colno,grey_shade);
  }

  /* Loop to write out colour definitions to PostScript file */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
  {
      if (Colour_Table_Name[icolour][0] != '\0')
      {
          fprintf(ps_file,"/%s { %8.4f %8.4f %8.4f",
              Colour_Table_Name[icolour],Colour_Table[icolour][0],
              Colour_Table[icolour][1],Colour_Table[icolour][2]);
          fprintf(ps_file," setrgbcolor } def\n");
      }
  }

/* v.4.1--> */
  /* If black-and-white option selected, overwrite all colours with their
     grey equivalents */
  if (Picture->In_Colour == FALSE)
    {
      /* Set object colours to their grey equivalents */
      for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
        {
          /* Get string equivalent of grey value */
          sprintf(colno,"%2d",object_grey[icolour]);
          if (colno[0] == ' ')
            colno[0] = '0';

          /* Form the string with the name of the grey colour */
          strcpy(object_colour[icolour],"Grey");
          strcat(object_colour[icolour],colno);
        }

      /* Repeat for text colours */
      for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
        {
          /* Get string equivalent of grey value */
          sprintf(colno,"%2d",text_grey[icolour]);
          if (colno[0] == ' ')
            colno[0] = '0';

          /* Form the string with the name of the grey colour */
          strcpy(text_colour[icolour],"Grey");
          strcat(text_colour[icolour],colno);
        }
    }
/* <--v.4.1 */

/* v.5.2 --> */
  /* For LigSearch annotation, write out residue-similarity colours */
  if (ligsearch_annot == TRUE)
    {
      /* Loop over the similarity levels */
      for (i = 0; i < NRES_SIMTYPES; i++)
        {
          /* Write out this colour definition */
          fprintf(ps_file,"/%s { %8.4f %8.4f %8.4f",
              SIMILARITY_DESC[i],SIMILARITY_COLOUR[i][0],
              SIMILARITY_COLOUR[i][1],SIMILARITY_COLOUR[i][2]);
          fprintf(ps_file," setrgbcolor } def\n");
        }

      /* Write out the spoke width */
      fprintf(ps_file,"/Spoke_width {    2.000 } def\n");
      fprintf(ps_file,"/Cut_width   {    2.000 } def\n");
    }
/* <-- v.5.2 */

  /* Write out the object- and text-colours */
  define_object_colours();
  define_text_colours();
  define_object_sizes(1.0,FALSE);
  define_text_sizes(1.0,FALSE);

  /* Print the remainder of the definitions */
  fprintf(ps_file,"/Title_size { %8.3f } def\n",TITLE_SIZE);
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
  fprintf(ps_file,"/Gprint { /Symbol findfont exch scalefont setfont show");
  fprintf(ps_file," } bind def\n");
/* v.4.4.1--> */
  fprintf(ps_file,"/Font { /Times-Roman findfont } bind def\n");
  fprintf(ps_file,"/Gfont { /Symbol findfont } bind def\n");
  fprintf(ps_file,"/Fprint { scalefont setfont show } bind def\n");
/* <--v.4.4.1 */
  fprintf(ps_file,"/Center {\n");
  fprintf(ps_file,"  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_file,"  exch stringwidth pop -2 div exch -3 div rmoveto\n");
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
  fprintf(ps_file,"/RotAngle { gsave currentpoint translate rotate } bind def\n");
/* <--v.5.3.8 */
  fprintf(ps_file,"%%%%EndProlog\n");
  fprintf(ps_file,"%%%%BeginSetup\n");
  fprintf(ps_file,"1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray");
  fprintf(ps_file," [ ] 0 setdash newpath\n");
  fprintf(ps_file,"%%%%EndSetup\n");
  fprintf(ps_file,"%%%%Page: 1 1 \n");
  fprintf(ps_file,"/LigplotSave save def\n");
/* v.5.1.7--> */
  fprintf(ps_file,"1.000 setgray\n");
/* <--v.5.1.7 */
  fprintf(ps_file,"%6.2f%6.2f moveto %7.2f%7.2f lineto%7.2f%7.2f lineto\n",
          xx1,xy1,xx2,xy1,xx2,xy2);
  fprintf(ps_file,"%7.2f%7.2f lineto closepath\n",xx1,xy2);
  fprintf(ps_file,"gsave 1.0000 setgray fill grestore\n");
  fprintf(ps_file,"stroke gsave\n");

  /*Draw in the background colour*/
  if (Picture->In_Colour == TRUE)
    {
      pscomm_("Background");
      pshade_(Colour->Background);
      psubox_(BBOXX1, BBOXY1, BBOXX2, BBOXY1, BBOXX2, BBOXY2, BBOXX1, BBOXY2);
    }

/* v.5.2 --> */
  /* For LigSearch annotation, write out the id */
  if (ligsearch_annot == TRUE)
    {
      /* Write out the identifier */
      strcpy(filename,"LigSearch id: ");
      strcat(filename,ligsearch_id);
      pscomm_(filename);
    }
/* <-- v.5.2 */
}

/****************   P L O T T I N G    R O U T I N E S   ***************/

/***********************************************************************

adjust_for_landscape  -  Perform boundary adjustments for printing
                         in Landscape, as opposed to Portrait,
                         orientation

***********************************************************************/

void adjust_for_landscape(void)
{

    /* Make necessary adjustments */
    pscomm_("Landscape orientation");
    Page_Min_x = BBOXY1;
    Page_Max_x = BBOXY2;
    Page_Min_y = BBOXX1;
    Page_Max_y = BBOXX2;
    Plot_Centre_x = (Page_Min_x + Page_Max_x) / 2.0;
    Plot_Centre_y = (Page_Min_y + Page_Max_y) / 2.0;
}
/***********************************************************************

print_title  -  Print the title for the plot

***********************************************************************/

void print_title(char string[TITLE_LEN + 1])
{
  float y;

  /* Define the margin taken up by the title area */
  y = Margin_y + TITLE_FRACTION * (BBOXY2 - BBOXY1) / 2.0;
  Margin_y = Margin_y + TITLE_FRACTION * (BBOXY2 - BBOXY1);

  /* Print the title */
  pscomm_("Plot heading");
  pscolb_(Text_Colour->Title);
  psctxt_(Plot_Centre_x,y,"Title_size",string);
}
/***********************************************************************

print_key_text  -  Print the text items for the key

***********************************************************************/

void print_key_text(float x_pos, float y_pos, float xmin, float ymin,
    float width, float height, char *text_size,char *text_string)
{
  float x, y;

  /* Print the text item */
  x = xmin + x_pos * width;
  y = ymin + y_pos * height;
  pstext_(x,y,text_size,text_string);
}
/***********************************************************************

plot_spokes  -  Plot the spokes to signify hydrophobic contacts

***********************************************************************/

void plot_spokes(float x, float y, float centre_x, float centre_y,
                 float scale_factor, float outer_radius, float radius,
                 float arc_angle, float direction_angle, float ang_step)
{
  int done, sign;
  float angle, angle_radians;
  float tot_step;
  float xc, x1, x2, yc, y1, y2;
  float ps_coord_x1, ps_coord_x2, ps_coord_y1, ps_coord_y2;

  /* Loop twice to plot the spokes in the clockwise and anti-clockwise
     directions */
  for (sign = -1; sign < 2; sign = sign + 2)
  {
      /* Initialise start for plotting in this direction */
      angle = ang_step * (int)(0.5 + direction_angle / ang_step);
      tot_step = 0.0;
      done = FALSE;
      while (done == FALSE)
      {
          /* Calculate angle in radians and start- and end-points of spoke */
          angle_radians = angle / RADDEG;
          xc = radius * cos(angle_radians) / scale_factor;
          yc = radius * sin(angle_radians) / scale_factor;
          x1 = x + xc;
          y1 = y + yc;
          xc = outer_radius * cos(angle_radians) / scale_factor;
          yc = outer_radius * sin(angle_radians) / scale_factor;
          x2 = x + xc;
          y2 = y + yc;

          /* Convert to PostScript coordinates */
          ps_coord_x1 = scale_factor * (x1 - centre_x) + Plot_Centre_x;
          ps_coord_x2 = scale_factor * (x2 - centre_x) + Plot_Centre_x;
          ps_coord_y1 = scale_factor * (y1 - centre_y) + Plot_Centre_y;
          ps_coord_y2 = scale_factor * (y2 - centre_y) + Plot_Centre_y;

          /* Plot this spoke */
          psline_(ps_coord_x1,ps_coord_y1,ps_coord_x2,ps_coord_y2);

          /* Increment the angle for the next spoke */
          angle = angle + sign * ang_step;
          if (angle < 0.0)
              angle = angle + 360.0;
          if (angle > 360.0)
              angle = angle - 360.0;

          /* Check whether have done all the spokes required in this
             direction */
          tot_step = tot_step + ang_step;
          if (tot_step > arc_angle / 2.0)
              done = TRUE;
      }
  }
}
/***********************************************************************

plot_contact_group  -  Plot a symbol for hydrophobic contacts - the
                       direction is defined as the start and end of
                       the virtual bond joining the contact group to
                       the atom it makes contact with

***********************************************************************/

void plot_contact_group(float x, float y, float dir_x, float dir_y,
                        float centre_x, float centre_y, float scale_factor,
                        float radius, float arc_angle,
/* v.4.0--> */
/*                        int spokes_only) */
                        int spokes_only,int interface)
/* <--v.4.0 */
{
  float angle_end, angle_start, ang_step, outer_radius;
  float direction_angle;
  float x1, y1;
  float ps_coord_x, ps_coord_y, ps_radius;

  /* Initialise values */
  ang_step = SPOKE_ANGLE;
  if (spokes_only == TRUE)
  {
      radius = radius * (1.75 / 1.4);
      ang_step = SPOKE_ANGLE * 3.0 / 2.0;
  }
  outer_radius = SPOKE_EXTENT * radius;

  /* Initialise PostScript output */
  pssave_();
/* v.4.0--> */
/*  pscolb_(Colour->Hydrophobics); */
  if (interface != 1)
    pscolb_(Colour->Hydrophobics);
  else
    pscolb_(Colour->Ligand_Bonds);
/* <--v.4.0 */

  /* Calculate the angle defining the direction of this group */
  x1 = dir_x - x;
  y1 = dir_y - y;
  if (fabs(x1) > 0.0001)
  {
      direction_angle = RADDEG * atan(y1 / x1);
      if (direction_angle >= 0.0)
      {
          if (x1 < 0.0 && y1 <= 0.0)
            direction_angle = direction_angle + 180.0;
      }
      else
      {
          direction_angle = direction_angle + 360.0;
          if (y1 > 0.0)
            direction_angle = direction_angle - 180.0;
      }
  }
  else if (fabs(y1) < 0.0001)
      direction_angle = 0.0;
  else if (y1 > 0.0)
      direction_angle = 90.0;
  else
      direction_angle = 270.0;

  /* Calculate start- and end-angles of the arc */
  angle_start = direction_angle - arc_angle / 2.0;
  if (angle_start < 0.0)
      angle_start = angle_start + 360.0;
  angle_end = direction_angle + arc_angle / 2.0;
  if (angle_end > 360.0)
      angle_end = angle_end - 360.0;

  /* Plot the arc */
  ps_coord_x = scale_factor * (x - centre_x) + Plot_Centre_x;
  ps_coord_y = scale_factor * (y - centre_y) + Plot_Centre_y;
  ps_radius = radius;
  if (spokes_only == FALSE)
  {
      pslwid_("Arc_linewidth");
      psarc_(ps_coord_x,ps_coord_y,ps_radius,angle_start,angle_end);
  }
  else
      pslwid_("Default_linewidth");

  /* Plot the spokes for the interacting residue */
  plot_spokes(x,y,centre_x,centre_y,scale_factor,outer_radius,radius,
              arc_angle,direction_angle,ang_step);

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

accessibility_circle  -  Plot the coloured/shaded circle representing
                         the current atom's accessibility

***********************************************************************/

void accessibility_circle(float ps_coord_x, float ps_coord_y,
    float min_circle_size, float max_circle_size, float scale_factor,
    float accessibility)
{
    int icol;
    float angle, circle_size, col[3], ratio, shade;

    /* Calculate the accessibility ratio */
    ratio = accessibility / Maximum_accessibility;
    if (ratio > 1.0)
        ratio = 1.0;
    if (ratio < 0.0)
        ratio = 0.0;
    if (xsite_file == TRUE)
    {
        if (ratio < 0.1)
            ratio = 0.0;
        else if (ratio > 0.25)
            ratio = 1.0;
        else
            ratio = (ratio - 0.1) / 0.15;
    }

    /* Compute the corresponding shade for the circle representing the
       accessibility */
    angle = PI_BY_2 * sqrt(ratio);
    shade = sin(angle);
    for (icol = 0; icol < 3; icol++)
    {
         col[icol] = Accessibility_Min[icol]
             + shade * (Accessibility_Max[icol] - Accessibility_Min[icol]);
         if (xsite_file == TRUE)
             col[icol] = sin(angle);
    }

    /* Draw in the circle of the appropriate size and shade */
    pscrgb_(col[0],col[1],col[2]);
    circle_size = min_circle_size + (max_circle_size - min_circle_size)
        * (1.0 - ratio);
    psucir_(ps_coord_x,ps_coord_y,circle_size * scale_factor);
}
/***********************************************************************

print_symbol_key  -  Print the key to the symbols used in the
                     figure
***********************************************************************/

void print_symbol_key(void)
{
  float margin;
  float yscale, width, text_size, x, x1, x2, xmid, xmin, y, y1, y2, ymin;
  float key_size, key_x, key_y;
  float ligand_atom_size, nonligand_atom_size;
  float atom_x1_1, atom_x1_2, atom_x1_3;
  float atom_x2_1, atom_x2_2, atom_x2_3;
  float atom_x3, atom_x4, atom_x5;
  float atom_y11, atom_y12, atom_y13, atom_y21, atom_y22, atom_y23,
        atom_y3;
  float bond1, bond2, bond3, ps_scale;
  float text_x1_1, text_x1_2, text_x1_3;
  float text_x2, text_x31, text_x32, text_x33;
  float text_y11, text_y12, text_y13, text_y21, text_y22, text_y23,
        text_y3, text_drop;

  /* Define positioning of elements */
  key_x = 0.03;
  key_y = 0.85;
  key_size = 0.17;
  text_size = 0.09;
  text_drop = text_size / 4.0;

  /* Define sizes, ensuring that not too large for key-area */
  ps_scale = 20.0;
  ligand_atom_size = ps_scale *  Size_Val->Ligand_Atoms;
  if (ligand_atom_size > 5.0)
    ligand_atom_size = 5.0;
  nonligand_atom_size = ps_scale *  Size_Val->Nonligand_Atoms;
  if (nonligand_atom_size > 5.0)
    nonligand_atom_size = 5.0;
  bond1 = ps_scale * Size_Val->Ligand_Bonds;
  bond2 = ps_scale * Size_Val->Nonligand_Bonds;
  bond3 = ps_scale * Size_Val->Hydrogen_Bonds;

  /* Bond positioning */
  atom_x1_1 = 0.06;
  atom_x2_1 = 0.12;
  atom_x1_2 = 0.307;
  atom_x2_2 = 0.367;
  atom_x1_3 = 0.59;
  atom_x2_3 = 0.65;
  text_x1_1 = 0.15;
  text_x1_2 = 0.397;
  text_x1_3 = 0.68;
  atom_y11 = 0.64;
  atom_y12 = 0.48;
  atom_y13 = 0.32;
  text_y11 = atom_y11 - text_drop;
  text_y12 = atom_y12 - text_drop;
  text_y13 = atom_y13 - text_drop;

  /* Hydrophobic contacts */
  atom_x3 = 0.46;
  text_x2 = 0.50;
  atom_y21 = 0.62;
  atom_y22 = 0.52;
  atom_y23 = 0.32;
  text_y21 = atom_y21 - text_drop;
  text_y22 = atom_y22 - text_drop;
  text_y23 = atom_y23 - text_drop;

  /* Accessibility shading */
  atom_x4 = 0.332;
  atom_x5 = 0.46;
  text_x31 = 0.052;
  text_x32 = 0.375;
  text_x33 = 0.50;
  atom_y3 = 0.09;
  text_y3 = atom_y3 - text_drop;

  /* Determine size of margin according to what is to be plotted */
  margin = 0.0;
  if (Include->Accessibilities == TRUE)
      margin = MARGIN_ACC * (BBOXY2 - BBOXY1);

  /* If no hydrophobic interactions, then show bonds on a single line */
  if (Include->Hydrophobics == FALSE)
  {
      margin = margin - MARGIN_NO_HPH * (BBOXY2 - BBOXY1);
      atom_y11 = atom_y13;
      atom_y12 = atom_y13;
      text_y11 = text_y13;
      text_y12 = text_y13;
      key_y = 0.53;
  }
  else
  {
      atom_x1_2 = atom_x1_1;
      atom_x1_3 = atom_x1_1;
      atom_x2_2 = atom_x2_1;
      atom_x2_3 = atom_x2_1;
      text_x1_2 = text_x1_1;
      text_x1_3 = text_x1_1;
  }

  /* Define border round the key area */
  x1 = (float)BBOXX1;
  x2 = (float)BBOXX2;
  width = x2 - x1;
  y1 = Margin_y;
  y2 = y1 + margin + KEY_HEIGHT * (BBOXY2 - BBOXY1);
  yscale = KEY_SCALE * (BBOXY2 - BBOXY1);
  Margin_y = y2;

  /* Determine minimum points for text */
  xmin = x1;
  if (Picture->Portrait == FALSE)
      xmin = x1 + (Page_Max_x - Page_Min_x - width) / 2.0;
  ymin = y1;
  if (Include->Accessibilities == TRUE)
      ymin = y1 + MARGIN_ACC * (BBOXY2 - BBOXY1);

  /* Write comment to PostScript file */
  pscomm_("KEY TO SYMBOLS");
  pssave_();

  /* Define the various sizes for text and objects */
  fprintf(ps_file,"/%s { %8.3f } def\n","Default_linewidth",0.2);
  fprintf(ps_file,"/%s { %8.3f } def\n","Arc_linewidth",0.75);
  fprintf(ps_file,"/%s { %8.3f } def\n","Key_size",key_size * yscale);
  fprintf(ps_file,"/%s { %8.3f } def\n","Text_size",text_size * yscale);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Ligand_Atoms,ligand_atom_size);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Nonligand_Atoms,
          nonligand_atom_size);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Ligand_Bonds,bond1);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Nonligand_Bonds,bond2);
  fprintf(ps_file,"/%s { %8.3f } def\n",Size->Hydrogen_Bonds,bond3);
  fprintf(ps_file,"/%s { %8.3f } def\n",Text_Size->Hbond_Lengths,
      text_size * yscale / 2.0);
  fprintf(ps_file,"/%s { %8.3f } def\n",Text_Size->Hydrophobic_Names,
      text_size * yscale * 3.0 / 4.0);

  /* Print key heading */
  pscolb_(Text_Colour->Key_Text);
  pslwid_("Default_linewidth");
  print_key_text(key_x,key_y,xmin,ymin,width,yscale,"Key_size","Key");

  /* Print ligand bond */
  pscomm_("Ligand bond");
  pssave_();
  x1 = xmin + atom_x1_1 * width;
  x2 = xmin + atom_x2_1 * width;
  xmid = (x1 + x2) / 2.0;
  y = ymin + atom_y11 * yscale;
  pslwid_(Size->Ligand_Bonds);

  /* If bond-split required, then print in two halves */
  if (Split_Colour_Ligand_Bonds == TRUE)
    {
      pscolb_(Colour->Nitrogen);
      psline_(x1,y,xmid,y);
      pscolb_(Colour->Carbon);
      psline_(x2,y,xmid,y);
    }
  else
    {
      pscolb_(Colour->Ligand_Bonds);
      psline_(x1,y,x2,y);
    }
  psrest_();

  /* Print the two atoms */
  pssave_();
  if (Include->Ligand_Atoms == TRUE)
    {
      pscolb_(Colour->Atom_Edges);
      psphcl_(Colour->Nitrogen);
      pspher_(x1,y,Size->Ligand_Atoms);
      psphcl_(Colour->Carbon);
      pspher_(x2,y,Size->Ligand_Atoms);
    }
  psrest_();

  /* Print the text */
  if (Include->Simple_Ligand_Residues == TRUE)
      print_key_text(text_x1_1,text_y11,xmin,ymin,width,yscale,
          "Text_size","Ligand m/c");
  else
/* v.4.0--> */
    {
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
        print_key_text(text_x1_1,text_y11,xmin,ymin,width,yscale,
                       "Text_size","Ligand bond");
/* v.4.0--> */
      else
        print_key_text(text_x1_1,text_y11,xmin,ymin,width,yscale,
                       "Text_size","Residues of first surface");
    }
/* <--v.4.0 */

  /* Print non-ligand bond */
  pscomm_("Nonligand bond");
  pssave_();
  x1 = xmin + atom_x1_2 * width;
  x2 = xmin + atom_x2_2 * width;
  xmid = (x1 + x2) / 2.0;
  y = ymin + atom_y12 * yscale;
  pslwid_(Size->Nonligand_Bonds);

  /* If bond-split required, then print in two halves */
  if (Split_Colour_Nonligand_Bonds == TRUE)
    {
      pscolb_(Colour->Nitrogen);
      psline_(x1,y,xmid,y);
      pscolb_(Colour->Carbon);
      psline_(x2,y,xmid,y);
    }
  else
    {
      pscolb_(Colour->Nonligand_Bonds);
      psline_(x1,y,x2,y);
    }
  psrest_();

  /* Print the two atoms */
  pssave_();
  if (Include->Nonligand_Atoms == TRUE)
    {
      pscolb_(Colour->Atom_Edges);
      psphcl_(Colour->Nitrogen);
      pspher_(x1,y,Size->Nonligand_Atoms);
      psphcl_(Colour->Carbon);
      pspher_(x2,y,Size->Nonligand_Atoms);
    }
  psrest_();

  /* Print the text */
  if (Include->Simple_Ligand_Residues == TRUE)
      print_key_text(text_x1_2,text_y12,xmin,ymin,width,yscale,
                     "Text_size","Ligand sidechain");
  else
/* v.4.0--> */
    {
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
        print_key_text(text_x1_2,text_y12,xmin,ymin,width,yscale,
                       "Text_size","Non-ligand bond");
/* v.4.0--> */
      else
        print_key_text(text_x1_2,text_y12,xmin,ymin,width,yscale,
                       "Text_size","Residues of second surface");
    }
/* <--v.4.0 */

  /* Print H-bond */
  if (Include->Hbonds == TRUE)
    {
      pscomm_("Hydrogen bond");
      pssave_();
      x1 = xmin + atom_x1_3 * width;
      x2 = xmin + atom_x2_3 * width;
      y = ymin + atom_y13 * yscale;
      pscolb_(Colour->Hydrogen_Bonds);
      pslwid_(Size->Hydrogen_Bonds);
      psdash_(2);
      psline_(x1,y,x1 + 0.3 * (x2 - x1),y);
      psline_(x2,y,x1 + 0.6 * (x2 - x1),y);
      psdash_(0);
      psrest_();

      /* Print the length */
      x = (x1 + x2) / 2.0;
      pssave_();
      pscolb_(Text_Colour->Hbond_Lengths);
      psctxt_(x,y,Text_Size->Hbond_Lengths,"3.0");
      psrest_();

      /* Print the two atoms */

      /* Print ligand nitrogen atom */
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
        {
          pscolb_(Colour->Atom_Edges);
          psphcl_(Colour->Nitrogen);
          pspher_(x1,y,Size->Ligand_Atoms);
        }

      /* If not plotting atoms, then show tip of the bond */
      else
        {
          if (Split_Colour_Ligand_Bonds == TRUE)
            psccol_(Colour->Nitrogen);
          else
            psccol_(Colour->Ligand_Bonds);
          psucir_(x1,y,bond1 / 2.0);
        }
      psrest_();

      /* Print non-ligand oxygen atom */
      pssave_();
      if (Include->Nonligand_Atoms == TRUE)
        {
          pscolb_(Colour->Atom_Edges);
          psphcl_(Colour->Oxygen);
          pspher_(x2,y,Size->Nonligand_Atoms);
        }

      /* If not plotting atoms, then show tip of the bond */
      else
        {
          if (Split_Colour_Nonligand_Bonds == TRUE)
            psccol_(Colour->Oxygen);
          else
            psccol_(Colour->Nonligand_Bonds);
          psucir_(x2,y,bond2 / 2.0);
        }
      psrest_();

      /* Print explanatory text */
      if (Include->Hbond_Lengths == TRUE)
        print_key_text(text_x1_3,text_y13,xmin,ymin,width,yscale,
                       "Text_size","Hydrogen bond and its length");
      else
        print_key_text(text_x1_3,text_y13,xmin,ymin,width,yscale,
                       "Text_size","Hydrogen bond");
    }

  /* If hydrophobic contacts have been shown, explain their symbols */
  if (Include->Hydrophobics == TRUE)
  {
      /* Print residues involved in hydrophobic contacts */
      pscomm_("Hydrophobic interaction");
      pssave_();
      x = xmin + atom_x3 * width;
      y = ymin + atom_y21 * yscale;
      pscolb_(Text_Colour->Hydrophobic_Names);
      psctxt_(x,y,Text_Size->Hydrophobic_Names,"His 53");
      psrest_();

      /* Print the carbon atom from the ligand */
      x1 = xmin + atom_x3 * width;
      y1 = ymin + atom_y23 * yscale;
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
        {
          pscolb_(Colour->Atom_Edges);
          psphcl_(Colour->Carbon);
          pspher_(x1,y1,Size->Ligand_Atoms);
        }

      /* If not plotting atoms, then show tip of the bond */
      else
        {
          if (Split_Colour_Ligand_Bonds == TRUE)
            psccol_(Colour->Carbon);
          else
            psccol_(Colour->Ligand_Bonds);
          psucir_(x1,y1,bond1 / 2.0);
        }
      psrest_();

      /* Plot the spokes radiating from atom */
      plot_contact_group(x1 - Plot_Centre_x,y1 - Plot_Centre_y,
                         x - Plot_Centre_x,y - Plot_Centre_y,
/* v.4.0--> */
/*                         0.0,0.0,1.0,ligand_atom_size,120.0,TRUE); */
                         0.0,0.0,1.0,ligand_atom_size,120.0,TRUE,0);
/* <--v.4.0 */

      /* Plot the spokes rediating from hydrophobic group */
      plot_contact_group(x - Plot_Centre_x,y - Plot_Centre_y,
                         x1 - Plot_Centre_x,y1 - Plot_Centre_y,
/* v.4.0--> */
/*                         0.0,0.0,1.0,2.2 * ligand_atom_size,120.0,FALSE); */
                         0.0,0.0,1.0,2.2 * ligand_atom_size,120.0,FALSE,0);
/* <--v.4.0 */

      /* Explanatory text: Line 1 */
/* v.4.0--> */
      if (Interface_Plot == FALSE)
/* <--v.4.0 */
        print_key_text(text_x2,text_y21,xmin,ymin,width,yscale,"Text_size",
                       "Non-ligand residues involved in hydrophobic");
/* v.4.0--> */
      else
        print_key_text(text_x2,text_y21,xmin,ymin,width,yscale,"Text_size",
                       "Residues involved in hydrophobic");
/* <--v.4.0 */

      /* Explanatory text: Line 2 */
      print_key_text(text_x2,text_y22,xmin,ymin,width,yscale,
          "Text_size","contact(s)");

      /* Explanatory text: Line 3 */
      print_key_text(text_x2,text_y23,xmin,ymin,width,yscale,
          "Text_size",
          "Corresponding atoms involved in hydrophobic contact(s)");

  }

  /* If accessibility-shading on the plot, explain the shading */
  if (Include->Accessibilities == TRUE)
  {
      /* Show shadings/colours for buried and accessible atoms */
      pscomm_("Accessibility shades/colours");
      pssave_();

      /* Explanatory text: Line 1 */
      print_key_text(text_x31,text_y3,xmin,ymin,width,yscale,
          "Text_size","Solvent accessibility shading:");

      /* Shading for completely buried atoms */
      x = xmin + atom_x4 * width;
      y = ymin + atom_y3 * yscale;
      accessibility_circle(x,y,2.1 * ligand_atom_size,
                           3.0 * ligand_atom_size,1.0,0.0);
      psccol_(Colour->Background);
      psucir_(x,y,2.1 * ligand_atom_size);

      /* Show the atom */
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
        {
          pscolb_(Colour->Atom_Edges);
          psphcl_(Colour->Nitrogen);
          pspher_(x,y,Size->Ligand_Atoms);
        }

      /* If not plotting atoms, then show tip of the bond */
      else
        {
          if (Split_Colour_Ligand_Bonds == TRUE)
            psccol_(Colour->Nitrogen);
          else
            psccol_(Colour->Ligand_Bonds);
          psucir_(x,y,bond1 / 2.0);
        }
      psrest_();

      /* Explanatory text: Line 2 - Buried shade/colour */
      print_key_text(text_x32,text_y3,xmin,ymin,width,yscale,
          "Text_size","Buried");

      /* Shading for completely accessible atoms */
      x = xmin + atom_x5 * width;
      y = ymin + atom_y3 * yscale;
      accessibility_circle(x,y,2.1 * ligand_atom_size,
                           3.0 * ligand_atom_size,1.0,
                           0.9 * Maximum_accessibility);
      psccol_(Colour->Background);
      psucir_(x,y,2.1 * ligand_atom_size);

      /* Show the atom */
      pssave_();
      if (Include->Ligand_Atoms == TRUE)
        {
          pscolb_(Colour->Atom_Edges);
          psphcl_(Colour->Nitrogen);
          pspher_(x,y,Size->Ligand_Atoms);
        }

      /* If not plotting atoms, then show tip of the bond */
      else
        {
          if (Split_Colour_Ligand_Bonds == TRUE)
            psccol_(Colour->Nitrogen);
          else
            psccol_(Colour->Ligand_Bonds);
          psucir_(x,y,bond1 / 2.0);
        }
      psrest_();
      psrest_();

      /* Explanatory text: Line 2 - Accessible shade/colour */
      print_key_text(text_x33,text_y3,xmin,ymin,width,yscale,
          "Text_size","Highly accessible");
  }

  /* Reset graphics defaults */
  psrest_();
}
/***********************************************************************

extract_minmax  -  Determine the max and min x and y coordinates of the
                   structure

***********************************************************************/

void extract_minmax(float *min_x,float *max_x,float *min_y,float *max_y)
{
  int first;

  float atom_size, minx, maxx, miny, maxy, x, y;

  struct coordinate *atom_ptr;

  /* Initialise variables */
  first = TRUE;
  minx = miny = 0.0;
  maxx = maxy = 1.0;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Get this atom's relative size */
      atom_size = atom_ptr->atom_size;

      /* Get atom's x- and y-coordinates */
      x = atom_ptr->x;
      y = atom_ptr->y;

      /* If this is the first atom, then use it to initialise the
         minimum and maximum values */
      if (first == TRUE)
        {
          minx = x - atom_size;
          maxx = x + atom_size;
          miny = y - atom_size;
          maxy = y + atom_size;
          first = FALSE;
        }

      /* Otherwise, update maximum and minimum values */
      else
        {
          if ((x - atom_size) < minx)
            minx = x - atom_size;
          if ((x + atom_size) > maxx)
            maxx = x + atom_size;
          if ((y - atom_size) < miny)
            miny = y - atom_size;
          if ((y + atom_size) > maxy)
            maxy = y + atom_size;
        }

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }

  /* Return the computed values */
  *min_x = minx;
  *max_x = maxx;
  *min_y = miny;
  *max_y = maxy;
}
/***********************************************************************

scale  -  Scale coordinates to the PostScript page

***********************************************************************/

void scale(float *scale_factor,float min_x,float max_x,float min_y,
    float max_y)
{
  float diff_x, diff_y, picture_height, picture_width;

  /* Calculate size of picture in real coordinates */
  diff_x = max_x - min_x;
  diff_y = max_y - min_y;

  /* If residue-names being printed, extend the size by the label-width
     and label-height */
  Max_object_width = Max_Textwidth;
  Max_object_height = Max_Textsize;
  if (Max_object_size > Max_object_width)
    Max_object_width = Max_object_size;
  if (Max_object_size > Max_object_height)
    Max_object_height = Max_object_size;

  /* Compute extra margin required */
  diff_x = diff_x + 2.0 * Max_object_width;
  diff_y = diff_y + 2.0 * Max_object_height;

  /* Determine extent of PostScript page that is available for the
     picture area */
  picture_height = Page_Max_y - Margin_y;
  picture_width = Page_Max_x - Page_Min_x;
  Plot_Centre_x = Page_Min_x + picture_width / 2.0;
  Plot_Centre_y = Margin_y + picture_height / 2.0;

  /* Determine whether picture bounded in the x- or y-direction */
  if ((diff_y / diff_x) < (picture_height / picture_width))
    {
      *scale_factor = BORDER_MARGIN * picture_width / diff_x;
    }
  else
    {
      *scale_factor = BORDER_MARGIN * picture_height / diff_y;
    }
}
/***********************************************************************

centre_point  -  Find the centre x- and y-coordinate on the
                 PostScript page

***********************************************************************/

void centre_point(float *centre_x,float *centre_y,float scale_factor,
/* v.4.4.4--> */
/*             float min_x,float max_x,float min_y,float max_y) */
                  float min_x,float max_x,float min_y,float max_y,
/* v.5.1.4--> */
/*                  int *bounds) */
                  int *bounds,int Clip_Diagram)
/* <--v.5.1.4 */
/* <--v.4.4.4 */
{
/* v.4.3--> */
#define SMALL_IMAGE_REDUCTION  0.33

/*  float ps_maxx, ps_minx, ps_maxy, ps_miny; */
  float ps_maxx, ps_minx, ps_maxy, ps_miny, swap;
  float f;
/* <--v.4.3 */

  *centre_x = (max_x + min_x) / 2.0;
  *centre_y = (max_y + min_y) / 2.0;

  /* Calculate PostScript limits of current picture */
  ps_minx = scale_factor * (min_x - *centre_x - Max_object_width)
    + Plot_Centre_x;
  ps_maxx = scale_factor * (max_x - *centre_x + Max_object_width)
    + Plot_Centre_x;
  ps_miny = scale_factor * (min_y - *centre_y - Max_object_height)
    + Plot_Centre_y;
  ps_maxy = scale_factor * (max_y - *centre_y + Max_object_height)
    + Plot_Centre_y;
  if (ps_minx < Page_Min_x || Include->Key == TRUE)
    ps_minx = Page_Min_x;
  if (ps_maxx > Page_Max_x || Include->Key == TRUE)
    ps_maxx = Page_Max_x;
  if (ps_miny < Page_Min_y || Include->Key == TRUE || Print_Title[0] != '\0')
    ps_miny = Page_Min_y;
  if (ps_maxy > Page_Max_y)
    ps_maxy = Page_Max_y;

  /* Flip the limits for landscape orientation */
  if (Picture->Portrait == FALSE)
    {
/* v.4.3--> */
/*      ps_minx = BBOXX1;
      ps_maxx = BBOXX2;
      ps_miny = BBOXY1;
      ps_maxy = BBOXY2; */
      swap = ps_minx;
      ps_minx = ps_miny;
      ps_miny = swap;
      swap = ps_maxx;
      ps_maxx = ps_maxy;
      ps_maxy = swap;
/* <--v.4.3 */
    }

  /* Write out "convert" command to PostScript file as a comment (used
     when using "convert" from within a script file to generate a .gif
     file from the ligplot.ps PostScript file) */
  fprintf(ps_file,"%%convert -clip %dx%d+%d-%d\n",(int)(ps_maxx - ps_minx),
          (int)(ps_maxy - ps_miny),(int) ps_minx,(int) ps_miny);

/* v.4.4.4--> */
/* v.5.1.4--> */
/*  bounds[0] = (int) ps_minx;
  bounds[1] = (int) ps_miny;
  bounds[2] = (int) ps_maxx;
  bounds[3] = (int) ps_maxy;
  bounds[0] = (int) BBOXX1 - 1;
  bounds[1] = (int) BBOXY1 - 1;
  bounds[2] = (int) BBOXX2 + 1;
  bounds[3] = (int) BBOXY2 + 1; */
  /* Store the picture boundaries */
  if (Clip_Diagram == TRUE)
    {
      bounds[0] = (int) ps_minx;
      bounds[1] = (int) ps_miny;
      bounds[2] = (int) ps_maxx;
      bounds[3] = (int) ps_maxy;
    }
  else
    {
      bounds[0] = (int) BBOXX1 - 1;
      bounds[1] = (int) BBOXY1 - 1;
      bounds[2] = (int) BBOXX2 + 1;
      bounds[3] = (int) BBOXY2 + 1;
    }
/* <--v.5.1.4 */
/* <--v.4.4.4 */

/* v.4.4.2--> */
  /* Add equivalent -crop line for newer versions of convert */
  fprintf(ps_file,"%%convert -crop %dx%d+%d+%d\n",(int)(ps_maxx - ps_minx),
          (int)(ps_maxy - ps_miny),(int) (ps_minx- BBOXX1 + 1),
          (int) (ps_miny - BBOXY1 + 1));
/* <--v.4.4.2 */

/* v.4.3--> */
  /* Write out the landscape small version (for use in PDBsum only) */
  if (Picture->Portrait == FALSE)
    {
      f = SMALL_IMAGE_REDUCTION;
      ps_minx = (1.0 - f) * (BBOXX2 + BBOXX1) + f * ps_minx;
      ps_maxx = (1.0 - f) * (BBOXX2 + BBOXX1) + f * ps_maxx;
      ps_miny = ps_miny * SMALL_IMAGE_REDUCTION;
      ps_maxy = ps_maxy * SMALL_IMAGE_REDUCTION;

      /* Write out convert command (disguised for undisguising by PDBsum
         script) */
      fprintf(ps_file,"%%honvert -clip %dx%d+%d-%d\n",
              (int)(ps_maxx - ps_minx),(int)(ps_maxy - ps_miny),
              (int) ps_minx,(int) ps_miny);
/* v.4.4.2--> */
      /* Repeat for -crop line for newer versions of convert */
      fprintf(ps_file,"%%honvert -crop %dx%d+%d+%d\n",
              (int)(ps_maxx - ps_minx),(int)(ps_maxy - ps_miny),
              (int) (ps_minx- BBOXX1 + 1),(int) (ps_miny - BBOXY1 + 1));
/* <--v.4.4.2 */

      fprintf(ps_file,"%%small %5.2f %5.2f scale\n",SMALL_IMAGE_REDUCTION,
              SMALL_IMAGE_REDUCTION);
    }
/* <--v.4.3 */
}
/***********************************************************************
 
sort_accessibilities  -  Sort routine using bubble sort method. Sorts
                         the atomic accessibilities in ascending or
                         descending order

***********************************************************************/

void sort_accessibilities(struct coordinate **first_stack_ptr,int order)
{
  int flipped;

  float accessibility, other_accessibility;

  struct coordinate *atom_ptr, *last_atom_ptr, *other_atom_ptr;

  /* Initialise variables */
  flipped = TRUE;

  /* Loop until all the values have been sorted */
  while (flipped == TRUE)
    {
      /* Initialise flag */
      flipped = FALSE;

      /* Get pointer to the first of the atoms */
      atom_ptr = *first_stack_ptr;
      last_atom_ptr = NULL;

      /* Loop through all the atoms */
      while (atom_ptr != NULL)
        {
          /* Get the atom's accessibility */
          accessibility = order * atom_ptr->accessibility;

          /* Get pointer to the next atom */
          other_atom_ptr = atom_ptr->next_stack_ptr;

          /* If have the next atom, then proceed */
          if (other_atom_ptr != NULL)
            {
              /* Get the other atom's accessibility */
              other_accessibility = order * other_atom_ptr->accessibility;

              /* If other atom's accessibility is smaller,
                 then swap the order of the atoms */
              if (other_accessibility < accessibility)
                {
                  /* Set flag to indicate that swap occurred */
                  flipped = TRUE;

                  /* If this is the first atom, get the first pointer
                     to point to the second atom */
                  if (last_atom_ptr == NULL)
                      *first_stack_ptr = other_atom_ptr;

                  /* Otherwise, set the last atom's pointer to the
                     other atom */
                  else
                    last_atom_ptr->next_stack_ptr = other_atom_ptr;

                  /* Swap the pointers */
                  atom_ptr->next_stack_ptr = other_atom_ptr->next_stack_ptr;
                  other_atom_ptr->next_stack_ptr = atom_ptr;

                  /* Set the other atom as the last atom */
                  last_atom_ptr = other_atom_ptr;
                }

              /* Otherwise, set the current atom as the last atom */
              else
                last_atom_ptr = atom_ptr;
            }
          else
            last_atom_ptr = atom_ptr;

          /* Get pointer to the next atom in the stack */
          atom_ptr = last_atom_ptr->next_stack_ptr;
        }
    }
}
/***********************************************************************

shade_accessibilities  -  Shade in the accessibility values for each
                          atom in the ligand

***********************************************************************/

void shade_accessibilities(float centre_x, float centre_y,
                           float scale_factor)
{
  int iatom, natoms, order;
/* v.4.0--> */
  int plot_accessibility;
/* <--v.4.0 */

  float ps_coord_x, ps_coord_y;
  float min_circle_size, max_circle_size, circle_size;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *bond_atom_ptr[2], *first_stack_ptr,
  *last_stack_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise plot of atom accessibilities */
  printf("   Shading accessibilities ...\n");
  min_circle_size = 2.8 * Size_Val->Ligand_Atoms / scale_factor;
  max_circle_size = 4.0 * Size_Val->Ligand_Atoms / scale_factor;
  pscomm_("Atomic accessibilities");

  /* Adjust for maximum accessibility read in */
  if (Maximum_accessibility == 0.0)
    Maximum_accessibility = 1.0;

  /* For an XSITE file, set the maximum accessibility to 90.0% */
  if (xsite_file == TRUE)
    Maximum_accessibility = 100.0;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise first and last stack pointers */
  first_stack_ptr = NULL;
  last_stack_ptr = NULL;
  natoms = 0;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
/* v.4.0--> */
      /* Determine whether accessibility to be plotted */
      plot_accessibility = TRUE;
      if (object_ptr->object_type == HYDROPHOBIC)
        plot_accessibility = FALSE;
      else if (object_ptr->object_type != LIGAND)
        {
          if (object_ptr->object_type == WATER)
            plot_accessibility = FALSE;
          else if (Include->Ligand_Accessibilities_Only == TRUE)
            plot_accessibility = FALSE;
        }
/* <--v.4.0 */
      /* If object is not a hydrophobic group, then plot all its atoms */
/* v.4.0--> */
/*      if (object_ptr->object_type != HYDROPHOBIC) */
      if (plot_accessibility == TRUE)
/* <--v.4.0 */
        {
          /* Get pointer to the first of this object's bonds */
          object_bond_ptr = object_ptr->first_object_bond_ptr;
      
          /* Loop through all object's bonds */
          while (object_bond_ptr != NULL)
            {
              /* Get the bond pointer */
              bond_ptr = object_bond_ptr->bond_ptr;

              /* Get pointers to the bond's two atoms */
              bond_atom_ptr[0] = bond_ptr->first_atom_ptr;
              bond_atom_ptr[1] = bond_ptr->second_atom_ptr;

              /* Loop to plot the bond's two atoms */
              for (iatom = 0; iatom < 2; iatom++)
                {
                  /* Check that this atom hasn't already been plotted */
                  if (bond_atom_ptr[iatom]->checked == FALSE)
                    {
                      /* Add atom to the stack */

                      /* If this is the first atom on the stack, set
                         first stack pointer */
                      if (first_stack_ptr == NULL)
                        first_stack_ptr = bond_atom_ptr[iatom];

                      /* Otherwise, get the last atom to point to
                         the current atom */
                      else
                        last_stack_ptr->next_stack_ptr = bond_atom_ptr[iatom];

                      /* Save the current atom's pointer */
                      last_stack_ptr = bond_atom_ptr[iatom];
                      natoms++;

                      /* Update flag to indicate atom has already been
                         plotted */
                      bond_atom_ptr[iatom]->checked = TRUE;
                    }
                }
      
              /* Get pointer to this object's next bond entry */
              object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
            }
        }
      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Process if there are any atoms to be done */
  if (natoms > 0)
    {
      /* Set the sort-order to be ascending order */
      order = 1;

      /* If this is an XSITE file, then plot the shades in the reverse
         order so that the darker ones come out on top */
      if (xsite_file == TRUE)
        order = -1;

      /* Sort all the stored accessibility values in ascending order */
      sort_accessibilities(&first_stack_ptr,order);

      /* Plot all the accessibility values as shaded circles */
      
      /* Get pointer to the first of the atoms to be plotted */
      atom_ptr = first_stack_ptr;

      /* Loop through all the atoms in order of increasing accessibility */
      while (atom_ptr != NULL)
        {
          /* Calculate PostScript coords of the atom's centre */
          ps_coord_x = scale_factor * (atom_ptr->x - centre_x)
            + Plot_Centre_x;
          ps_coord_y = scale_factor * (atom_ptr->y - centre_y)
            + Plot_Centre_y;

          /* Set the shade/colour of the circle according
             to this atom's accessibility */
          accessibility_circle(ps_coord_x,ps_coord_y,min_circle_size,
                               max_circle_size,scale_factor,
                               atom_ptr->accessibility);

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next_stack_ptr;
        }
    }

  /* Set the shade to the background colour */
  pscomm_("Blank out centres");
  psccol_(Colour->Background);

  /* Get pointer to the first of the atoms again */
  atom_ptr = first_stack_ptr;

  /* Loop through all the atoms a second time to blank out the
       central region */
  while (atom_ptr != NULL)
    {
      /* Calculate PostScript coords of the atom's centre */
      ps_coord_x = scale_factor * (atom_ptr->x - centre_x)
        + Plot_Centre_x;
      ps_coord_y = scale_factor * (atom_ptr->y - centre_y)
        + Plot_Centre_y;

      /* Draw in the circle */
      circle_size = min_circle_size;
      psucir_(ps_coord_x,ps_coord_y,circle_size * scale_factor);

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_stack_ptr;
    }
}
/* v.5.2 --> */
/***********************************************************************

ligsearch_annotations  -  Annotate the ligand atoms according to the
                          relevance of the H-bonded and non-bonded
                          contacts they make to their host protein

***********************************************************************/

void ligsearch_annotations(float centre_x, float centre_y,
                           float scale_factor)
{
  char comment[100], ps_radius[20];

  int best, i, iatom, natoms;

  float ang_step, arc_angle, dash_step, direction_angle;
  float outer_radius, radius, spoke_radius, zero_radius;
  float ps_coord_x, ps_coord_y;
  float min_circle_size, max_circle_size, circle_size;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *bond_atom_ptr[2], *first_stack_ptr,
  *last_stack_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise plot of atom accessibilities */
  printf("   Annotating ligand atoms ...\n");
  ang_step = SPOKE_ANGLE * 2.0;
  dash_step = SPOKE_ANGLE * 3.0;
  min_circle_size = Size_Val->Ligand_Atoms / scale_factor;
  max_circle_size = 1.5 * Size_Val->Ligand_Atoms / scale_factor;
  pscomm_("LigSearch annotations");
  arc_angle = 360;
  direction_angle = 0;
  radius = min_circle_size * scale_factor;
  spoke_radius = 1.4 * max_circle_size * scale_factor;
  outer_radius = 3 * (radius + spoke_radius) / 4;
  zero_radius = 0;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise first and last stack pointers */
  first_stack_ptr = NULL;
  last_stack_ptr = NULL;
  natoms = 0;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a ligand, then plot all annotations that we have */
      if (object_ptr->object_type == LIGAND)
        {
          /* Get pointer to the first of this object's bonds */
          object_bond_ptr = object_ptr->first_object_bond_ptr;
      
          /* Loop through all object's bonds */
          while (object_bond_ptr != NULL)
            {
              /* Get the bond pointer */
              bond_ptr = object_bond_ptr->bond_ptr;

              /* Get pointers to the bond's two atoms */
              bond_atom_ptr[0] = bond_ptr->first_atom_ptr;
              bond_atom_ptr[1] = bond_ptr->second_atom_ptr;

              /* Loop to plot the bond's two atoms */
              for (iatom = 0; iatom < 2; iatom++)
                {
                  /* Check that this atom hasn't already been plotted */
                  if (bond_atom_ptr[iatom]->checked == FALSE)
                    {
                      /* Add atom to the stack if it has any annotation */
                      if (bond_atom_ptr[iatom]->have_score == TRUE)
                        {
                          /* If this is the first atom on the stack, set
                             first stack pointer */
                          if (first_stack_ptr == NULL)
                            first_stack_ptr = bond_atom_ptr[iatom];

                          /* Otherwise, get the last atom to point to
                             the current atom */
                          else
                            last_stack_ptr->next_stack_ptr
                              = bond_atom_ptr[iatom];

                          /* Save the current atom's pointer */
                          last_stack_ptr = bond_atom_ptr[iatom];
                          natoms++;
                        }

                      /* Update flag to indicate atom has already been
                         plotted */
                      bond_atom_ptr[iatom]->checked = TRUE;
                    }
                }
      
              /* Get pointer to this object's next bond entry */
              object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
            }
        }
      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Process if there are any atoms to be done */
  if (natoms > 0)
    {
      /* Plot all the annotations */
      
      /* Get pointer to the first of the atoms to be plotted */
      atom_ptr = first_stack_ptr;

      /* Loop through all the atoms */
      while (atom_ptr != NULL)
        {
          /* Calculate PostScript coords of the atom's centre */
          ps_coord_x = scale_factor * (atom_ptr->x - centre_x)
            + Plot_Centre_x;
          ps_coord_y = scale_factor * (atom_ptr->y - centre_y)
            + Plot_Centre_y;

          /* Determine whether the atom's 'best' interaction is
             H-bonded or otherwise */
          if (atom_ptr->score[S_HBONDS] > -1 ||
              atom_ptr->score[S_SALT_BRIDGES] > -1 ||
              atom_ptr->score[S_DISULPHIDES] > -1)
            {
              /* Initialise the best colour */
              best = S_NTYPES;

              /* Loop over the types to see which has the best
                 similarity score */
              for (i = 0; i < S_NTYPES; i++)
                {
                  /* If this is the best so far, then store */
                  if (i != S_CONTACTS && atom_ptr->score[i] < best &&
                      atom_ptr->score[i] != -1)
                    best = atom_ptr->score[i];
                }

              /* If have a best, then set colour accordingly */
              if (best < RES_MISSING)
                {
                  /* Show atom being annotated */
                  sprintf(comment,"Spokes for atom %s",
                          atom_ptr->atom_name);
                  pscomm_(comment);

                  /* Set the colour according to the best similarity
                     score */
                  pscolb_(SIMILARITY_DESC[best]);

                  /* Set the spoke line width */
                  pslwid_("Spoke_width");

                  /* Plot the spokes around this atom */
                  plot_spokes(atom_ptr->x,atom_ptr->y,centre_x,centre_y,
                              scale_factor,spoke_radius,radius,
                              arc_angle,direction_angle,ang_step);
                }
            }

          /* Otherwise plot a dotted line around the atom */
          else
            {
              /* Initialise the best colour */
              best = atom_ptr->score[S_CONTACTS];

              /* If have a best, then set colour accordingly */
              if (best != -1 && best != RES_MISSING)
                {
                  /* Show atom being annotated */
                  sprintf(comment,"Dotted line for atom %s",
                          atom_ptr->atom_name);
                  pscomm_(comment);

                  /* Set the colour according to the best similarity
                     score */
                  //pscrgb_((float) SIMILARITY_COLOUR[best][0],
                  //        (float) SIMILARITY_COLOUR[best][1],
                  //        (float) SIMILARITY_COLOUR[best][2]);
                  psccol_(Colour->Background);
                  pscolb_(SIMILARITY_DESC[best]);
                  pslwid_("Cut_width");

                  /* Draw a circle around the atom */
                  sprintf(ps_radius,"%10.3f",
                          max_circle_size * scale_factor);
                  pscirc_(ps_coord_x,ps_coord_y,ps_radius);

                  /* Break it up with white spokes to make the circles
                     dashed ones */
                  
                  /* Set the colour according to the best similarity
                     score */
                  pscolb_("Grey10");

                  /* Set the spoke line width */
                  pslwid_("Cut_width");

                  /* Plot the spokes around this atom */
                  plot_spokes(atom_ptr->x,atom_ptr->y,centre_x,centre_y,
                              scale_factor,outer_radius,zero_radius,
                              arc_angle,direction_angle,dash_step);
                }
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next_stack_ptr;
        }
    }

  /* Set the shade to the background colour */
  pscomm_("Blank out centres");
  psccol_(Colour->Background);

  /* Get pointer to the first of the atoms again */
  atom_ptr = first_stack_ptr;

  /* Loop through all the atoms a second time to blank out the
       central region */
  while (atom_ptr != NULL)
    {
      /* Calculate PostScript coords of the atom's centre */
      ps_coord_x = scale_factor * (atom_ptr->x - centre_x)
        + Plot_Centre_x;
      ps_coord_y = scale_factor * (atom_ptr->y - centre_y)
        + Plot_Centre_y;

      /* Draw in the circle */
      circle_size = min_circle_size;
      psucir_(ps_coord_x,ps_coord_y,circle_size * scale_factor);

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next_stack_ptr;
    }
}
/* <-- v.5.2 */
/***********************************************************************

double_bond_split  -  Calculate the coordinate shifts required for
                      printing a double- or triple_bond

***********************************************************************/

void double_bond_split(float x1,float y1,float x2,float y2,
                       float line_width,float scale_factor,
/* v.4.4--> */
/*                       int bond_order,float *ps_dx,float *ps_dy) */
                       int bond_order,float *ps_dx,float *ps_dy,
                       int chemical_notation)
/* <--v.4.4 */
{
  float dx, dy, gap, len, x, y;

  /*Initialize */
  dx = 0.0;
  dy = 0.0;

  /* Check that bond is either a double- or triple-bond */
/* v.4.3--> */
/*  if (bond_order == 2 || bond_order == 3) */
  if (bond_order == DOUBLE || bond_order == TRIPLE)
/* <--v.4.3 */
    {
      /* Get the lengths of the bond in the x- and y-directions */
      x = x2 - x1;
      y = y2 - y1;

      /* Get the size of the gap between the two parallel lines defining
         the current double-bond */
/* v.4.4--> */
/* v.4.3--> */
/*      gap = 0.75 * line_width / scale_factor;
      if (bond_order == 't')
        gap = 1.5 * line_width / scale_factor;
      gap = 1.0 * line_width / scale_factor;
      if (bond_order == TRIPLE)
        gap = 1.0 * line_width / scale_factor; */
      if (chemical_notation == TRUE)
        {
/* v.4.5.4--> */
          gap = 0.5 * line_width / scale_factor;
/* v.5.0.5--> */
          if (line_width / scale_factor < 0.19)
            gap = 0.85 * line_width / scale_factor;
/* <--v.5.0.5 */
/* <--v.4.5.4 */
          if (bond_order == TRIPLE)
            gap = 1.0 * line_width / scale_factor;
        }
      else
        {
          gap = 0.75 * line_width / scale_factor;
          if (bond_order == 't')
            gap = 1.5 * line_width / scale_factor;
        }
/* <--v.4.4 */
/* <--v.4.3 */

      /* Calculate special cases where bond is close to horizontal or
         close to vertical */
      if (fabs(x) < 0.001)
        {
          dx = gap;
          dy = 0.0;
        }
      else if (fabs(y) < 0.001)
        {
          dx = 0.0;
          dy = gap;
        }

      /* Otherwise calculate the shift-sizes in the x- and y-directions */
      else
        {
          len = x * x + y * y;
          len = sqrt((double) len);
          dx = gap * y / len;
          dy = gap * x / len; 
        }
    }

  /* Compute the shift-lengths as scaled to the PostScript coordinates */
  *ps_dx = scale_factor * dx;
  *ps_dy = scale_factor * dy;
}
/***********************************************************************

get_atom_colour  -  Get the colour of the current atom

***********************************************************************/

void get_atom_colour(char atom_name[5],char res_name[4],
                     char colour[COL_NAME_LEN + 1])
{
  /* Determine the atom-type and its colour */

  /* Nitrogen */
  if (atom_name[1] == 'N' && atom_name[0] != 'Z' && atom_name[0] != 'M')
    strcpy(colour,Colour->Nitrogen);

  /* Carbon */
  else if (atom_name[1] == 'C')
    strcpy(colour,Colour->Carbon);

  /* Water */
  else if (!strncmp(res_name,"HOH",3))
    strcpy(colour,Colour->Water);

  /* Oxygen */
  else if (atom_name[1] == 'O' && atom_name[0] != 'C' && atom_name[0] != 'M')
    strcpy(colour,Colour->Oxygen);

  /* Sulphur */
  else if (atom_name[1] == 'S')
    strcpy(colour,Colour->Sulphur);

  /* Phosphorus */
  else if (atom_name[1] == 'P')
    strcpy(colour,Colour->Phosphorus);

  /* Iron */
  else if (!strncmp(atom_name,"FE",2))
    strcpy(colour,Colour->Iron);

  /* Default colour */
  else
    strcpy(colour,Colour->Other);
}
/***********************************************************************

plot_bond_lines  -  Print all the bond lines in each object in turn
  
***********************************************************************/

void plot_bond_lines(float centre_x,float centre_y,float scale_factor)
{
  char atom1_colour[COL_NAME_LEN + 1], atom2_colour[COL_NAME_LEN + 1];
/* v.4.5.4--> */
  char number_string[20];
/* <--v.4.5.4 */

  int iobject, width_change;
  int object_type;

  float ps_dx, ps_dy, line_width;
  float ps_coord_ax, ps_coord_ay, ps_coord_bx, ps_coord_by;
  float ps_coord_midx, ps_coord_midy;
/* v.4.5.4--> */
  float width;
/* <--v.4.5.4 */
  
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise print */
  pscomm_("Bond lines");
  pssave_();
  line_width = Size_Val->Ligand_Bonds;
  ps_dx = 0.0;
  ps_dy = 0.0;

  /* Initialise pointer to the first stored object */
  printf("   Plotting covalent bonds ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get this object's type */
      object_type = object_ptr->object_type;

      /* Process only if not a hydrophobic group, and not a simplified
         H-group */
      if (object_type != HYDROPHOBIC && object_type != SIMPLE_HGROUP)
        {
          /* Set default line colour and thickness according to the
             object type */
          if (object_type == LIGAND)
            {
              if (Split_Colour_Ligand_Bonds == FALSE)
                pscolb_(Colour->Ligand_Bonds);
              pslwid_(Size->Ligand_Bonds);
/* v.4.5.4--> */
              line_width = Size_Val->Ligand_Bonds;
/* <--v.4.5.4 */
            }

          /* Default colours and lines thicknesses for non-ligand bonds */
          else
            {
              if (Split_Colour_Nonligand_Bonds == FALSE)
                pscolb_(Colour->Nonligand_Bonds);
              pslwid_(Size->Nonligand_Bonds);
/* v.4.5.4--> */
              line_width = Size_Val->Nonligand_Bonds;
/* <--v.4.5.4 */
            }

          /* Get pointer to the first of this object's bonds */
          object_bond_ptr = object_ptr->first_object_bond_ptr;

          /* Loop through all this object's bonds, one by one */
          while (object_bond_ptr != NULL)
            {
              /* Get the current bond */
              bond_ptr = object_bond_ptr->bond_ptr;

/* v.4.5.4--> */
              /* Initialise width-change flag */
              width_change = FALSE;
/* <--v.4.5.4 */

/* v.4.4--> */
              /* Check that hasn't already been plotted */
              if (bond_ptr->plotted == FALSE)
                {
/* <--v.4.4 */
                  /* Get pointers to the two atoms at either end of the
                     bond */
                  atom1_ptr = bond_ptr->first_atom_ptr;
                  atom2_ptr = bond_ptr->second_atom_ptr;

                  /* Compute the PostScript coordinates of the two ends
                     of the bond */
                  ps_coord_ax = Plot_Centre_x
                    + scale_factor * (atom1_ptr->x - centre_x);
                  ps_coord_ay = Plot_Centre_y
                    + scale_factor * (atom1_ptr->y - centre_y);
                  ps_coord_bx = Plot_Centre_x
                    + scale_factor * (atom2_ptr->x - centre_x);
                  ps_coord_by = Plot_Centre_y
                    + scale_factor * (atom2_ptr->y - centre_y);

                  /* If double- and triple-bonds to be plotted, calculate
                     coordinate-shifts required */
                  if (Include->Double_Bonds == TRUE)
                    double_bond_split(atom1_ptr->x,atom1_ptr->y,atom2_ptr->x,
                                      atom2_ptr->y,line_width,scale_factor,
                                      bond_ptr->bond_order,
/* v.4.4--> */
/*                                      &ps_dx,&ps_dy); */
                                      &ps_dx,&ps_dy,
                                      Include->Chemical_Notation);
/* <--v.4.4 */

                  /* Plot the bond */

                  /* If bond of a single colour, then plot */
                  if ((object_type == LIGAND &&
                       Split_Colour_Ligand_Bonds == FALSE) ||
                      (object_type == HGROUP &&
                       Split_Colour_Nonligand_Bonds == FALSE))
                    {
                      /* If this is a double- or triple-bond, plot a
                         double-line */
                      width_change = FALSE;
                      if (Include->Double_Bonds == TRUE &&
/* v.4.4--> */
/*                          (bond_ptr->bond_order == 2 ||
                           bond_ptr->bond_order == 3)) */
                          (bond_ptr->bond_order == DOUBLE ||
                           bond_ptr->bond_order == TRIPLE))
/* <--v.4.4 */
                        {
                          /* If this is a ligand bond, then need to
                             make thinner */
/* v.4.5.4--> */
/*                          if (object_type == LIGAND)
                            {
                              pslwid_(Size->Nonligand_Bonds);
                              width_change = TRUE;
                            }
*/
                          /* Get the line's size */
                          if (object_type == LIGAND)
                            width = Size_Val->Ligand_Bonds;
                          else
                            width = Size_Val->Nonligand_Bonds;

                          /* Make line thinner */
                          width = 0.5 * width;
                          sprintf(number_string,"%5.2f",width);

                          /* Set line width to new value */
                          pslwid_(number_string);

                          /* Set flag indicating we have changed the
                             line width */
                          width_change = TRUE;
/* <--v.4.5.4 */

                          /* Draw the two lines representing the double bond */
                          psline_(ps_coord_ax - ps_dx,ps_coord_ay + ps_dy,
                                  ps_coord_bx - ps_dx,ps_coord_by + ps_dy);
                          psline_(ps_coord_ax + ps_dx,ps_coord_ay - ps_dy,
                                  ps_coord_bx + ps_dx,ps_coord_by - ps_dy);
                        }

                      /* For a single or triple-bond, print line between
                         atom centres */
                      if (Include->Double_Bonds == FALSE ||
                          (Include->Double_Bonds == TRUE &&
/* v.4.4--> */
/*                           (bond_ptr->bond_order == 1 ||
                            bond_ptr->bond_order == 3))) */
                           (bond_ptr->bond_order == SINGLE ||
                            bond_ptr->bond_order == TRIPLE ||
                            bond_ptr->bond_order == AROMATIC)))
/* <--v.4.4 */
                        psline_(ps_coord_ax,ps_coord_ay,ps_coord_bx,
                                ps_coord_by);

                      /* If ligand-line was made thin, revert to usual width */
                      if (width_change == TRUE)
                        pslwid_(Size->Ligand_Bonds);
                    }

                  /* If split-colour bonds required, then process */
                  else
                    {
                      /* Calculate mid-point coordinates */
                      ps_coord_midx = (ps_coord_ax + ps_coord_bx) / 2.0;
                      ps_coord_midy = (ps_coord_ay + ps_coord_by) / 2.0;

                      /* Get the colours of the two atoms at either end of
                         the bond */
                      get_atom_colour(atom1_ptr->atom_name,
                                      atom1_ptr->residue_ptr->res_name,
                                      atom1_colour);
                      get_atom_colour(atom2_ptr->atom_name,
                                      atom2_ptr->residue_ptr->res_name,
                                      atom2_colour);

                      /* Plot the two halves of the line */

                      /* If this is a double-bond, plot a double-line */
                      if (Include->Double_Bonds == TRUE &&
                          (bond_ptr->bond_order == 2 ||
                           bond_ptr->bond_order == 3))
                        {
                          /* If this is a ligand bond, then need to
                             make thinner */
/* v.4.5.4--> */
/*                          if (object_type == LIGAND)
                            {
                              pslwid_(Size->Nonligand_Bonds);
                              width_change = TRUE;
                            }
*/
                          /* Get the line's size */
                          if (object_type == LIGAND)
                            width = Size_Val->Ligand_Bonds;
                          else
                            width = Size_Val->Nonligand_Bonds;

                          /* Make line thinner */
                          width = 0.5 * width;
                          sprintf(number_string,"%5.2f",width);

                          /* Set line width to new value */
                          pslwid_(number_string);

                          /* Set flag indicating we have changed the
                             line width */
                          width_change = TRUE;
/* <--v.4.5.4 */

                          /* Draw the two lines representing the double bond */
                          pscolb_(atom1_colour);
                          psline_(ps_coord_ax - ps_dx,ps_coord_ay + ps_dy,
                                  ps_coord_midx - ps_dx,ps_coord_midy + ps_dy);
                          psline_(ps_coord_ax + ps_dx,ps_coord_ay - ps_dy,
                                  ps_coord_midx + ps_dx,ps_coord_midy - ps_dy);
                          pscolb_(atom2_colour);
                          psline_(ps_coord_midx - ps_dx,ps_coord_midy + ps_dy,
                                  ps_coord_bx - ps_dx,ps_coord_by + ps_dy);
                          psline_(ps_coord_midx + ps_dx,ps_coord_midy - ps_dy,
                                  ps_coord_bx + ps_dx,ps_coord_by - ps_dy);
                        }

                      /* For a single or triple-bond, print line between
                         atom centres */
                      if (Include->Double_Bonds == FALSE ||
                          (Include->Double_Bonds == TRUE &&
/* v.4.4--> */
/*                           (bond_ptr->bond_order == 1 ||
                            bond_ptr->bond_order == 3))) */
                           (bond_ptr->bond_order == SINGLE ||
                            bond_ptr->bond_order == TRIPLE ||
                            bond_ptr->bond_order == AROMATIC)))
/* <--v.4.4 */
                        {
                          pscolb_(atom1_colour);
                          psline_(ps_coord_ax,ps_coord_ay,ps_coord_midx,
                                  ps_coord_midy);
                          pscolb_(atom2_colour);
                          psline_(ps_coord_midx,ps_coord_midy,ps_coord_bx,
                                  ps_coord_by);
              
                          /* If ligand-line was made thin, revert to
                             usual width */
/* v.4.5.4--> */
/*                          if (width_change == TRUE)
                            pslwid_(Size->Ligand_Bonds); */
/* <--v.4.5.4 */
                        }
                    }

/* v.4.4--> */
                  /* Mark bond as plotted */
                  bond_ptr->plotted = TRUE;
                }
/* <--v.4.4 */

/* v.4.5.4--> */
              /* If line-width was changed, then reset */
              if (width_change == TRUE &&
                  object_type == LIGAND)
                pslwid_(Size->Ligand_Bonds);
              else if (width_change == TRUE)
                pslwid_(Size->Nonligand_Bonds);
/* <--v.4.5.4 */

              /* Get pointer to this object's next bond entry */
              object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
            }
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Restore the graphics state */
    psrest_();
}
/* v.4.4--> */
/***********************************************************************

unset_bond_plot_flags  -  Re-initialise all the bond plot flags

***********************************************************************/

void unset_bond_plot_flags(void)
{
  struct bond *bond_ptr;

  /* Set pointer to first bond in linked-list */
  bond_ptr = first_bond_ptr;

  /* Loop through the stored bond-list to initialise flags and 
     stack-pointers */
  while (bond_ptr != NULL)
    {
      /* Initialise plotted marker */
      bond_ptr->plotted = FALSE;

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

plot_aromatics  -  Plot circles within any rings that are defined as
                   aromatic
  
***********************************************************************/

void plot_aromatics(float centre_x,float centre_y,float scale_factor)
{
#define CIR_FACTOR 0.70

  char ps_radius[20];

/* v.4.5.2--> */
  // int iatom, natoms, nbonds;
  int iatom, iring, natoms, nbonds, nrings;
/* <--v.4.5.2 */

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
  float dist, radius;
  float ps_coord_x, ps_coord_y;

  struct bond *bond_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initiales all the plotted flags */
  unset_bond_plot_flags();

  /* Get pointer to the first object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;

      /* Loop through all object's bonds */
      while (object_bond_ptr != NULL)
        {
          /* Get the bond pointer */
          bond_ptr = object_bond_ptr->bond_ptr;

          /* Process bond if valid */
          if (bond_ptr != NULL)
            {
              /* If this is an aromatic ring-bond which hasn't already been
                 processed, then find the ring it belongs to */
              if (bond_ptr->ring_bond == TRUE &&
                  bond_ptr->bond_order == AROMATIC &&
                  bond_ptr->plotted == FALSE)
                {
                  /* Get a set of ring bonds which form a ring with the
                     current bond */
                  get_other_ring_bonds(bond_ptr);

/* v.4.5.2--> */
                  /* Determine how many subrings there are in the given
                     ring system */
                  nrings = get_subrings(bond_ptr);

                  /* Loop over all the aromatic rings identified */
                  for (iring = 0; iring < nrings; iring++)
                    {
/* <--v.4.5.2 */

                      /* Extract all the coordinates of the current
                         ring's atoms so can draw circle within them */
/* v.4.5.2--> */
                      // get_ring_atom_coords(bond_ptr,coord_store,&natoms,
                      //                     &nbonds,TRUE);
                      get_ring_atom_coords(bond_ptr,coord_store,&natoms,
                                           &nbonds,iring + 1,TRUE);
/* <--v.4.5.2 */

                      /* If have at least 3 atoms, continue */
                      if (natoms > 2)
                        {
                          /* Have the coords of all the atoms making up the
                             current ring, so calculate the centre of mass
                             of these atoms */
                          centre_of_mass(coord_store,natoms,&mass_x,&mass_y,
                                         &mass_z);

                          /* Determine average distance from centre-point
                             to atoms of the ring */
                          radius = 0.0;
                          for (iatom = 0; iatom < natoms; iatom++)
                            {
                              /* Get distance and save */
                              dist = (mass_x - coord_store[iatom][0])
                                * (mass_x - coord_store[iatom][0])
                                + (mass_y - coord_store[iatom][1])
                                * (mass_y - coord_store[iatom][1])
                                + (mass_z - coord_store[iatom][2])
                                * (mass_z - coord_store[iatom][2]);
                              dist = sqrt(dist);
                              radius = radius + dist;
                            }
                          
                          /* Calculate the average distance */
                          radius = radius / natoms;
                          
                          /* Adjust to give radius of inscribed circle */
                          radius = CIR_FACTOR * radius;

                          /* Set circle colour */
                          psccol_(Colour->Background);

                          /* Set circle colour */
                          if (object_ptr->object_type == LIGAND)
                            {
                              pscolb_(Colour->Ligand_Bonds);
                              pslwid_(Size->Ligand_Bonds);
                            }
                          else
                            {
                              pscolb_(Colour->Nonligand_Bonds);
                              pslwid_(Size->Nonligand_Bonds);
                            }

                          /* Calculate coords of centre and scaled radius */
                          ps_coord_x = Plot_Centre_x
                            + scale_factor * (mass_x - centre_x);
                          ps_coord_y = Plot_Centre_y
                            + scale_factor * (mass_y - centre_y);
                          sprintf(ps_radius,"%10.3f",scale_factor * radius);

                          /* Plot the circle */
                          pscirc_(ps_coord_x,ps_coord_y,ps_radius);
                        }
/* v.4.5.2--> */
                    }
/* <--v.4.5.2 */
                }
            }

          /* Get pointer to point to next bond for this object */
          object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* <--v.4.4 */
/***********************************************************************

print_hbond_lengths  -  Prints H-bond lengths on the H-bond lines

***********************************************************************/

/* v.4.1--> */
/* void print_hbond_lengths(float bond_length,float x1,float y1,
                         float x2,float y2,float scale_factor,
                         float centre_x,float centre_y) */
void print_hbond_lengths(struct bond *bond_ptr,float x1,float y1,
                         float x2,float y2,float scale_factor,
                         float centre_x,float centre_y)
/* <--v.4.1 */
{
  char hbond_dist[5];
  float sizex, sizey;
  float ps_coord_x, ps_coord_y, xc, yc;

  /* Get the bond length */
/* v.4.1--> */
/*  sprintf(hbond_dist,"%4.2f",bond_length); */
  sprintf(hbond_dist,"%4.2f",bond_ptr->bond_length);
/* <--v.4.1 */
  hbond_dist[4] = '\0';

/* v.4.1--> */
  /* Store the bond-length string */
  strcpy(bond_ptr->label,hbond_dist);
/* <--v.4.1 */

  /* Get the PostScript position of the label */
  xc = (x1 + x2) / 2;
  yc = (y1 + y2) / 2;

  /* Save the coordinates of the label */
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = xc;
      label_coord[n_labels][1] = yc;
      n_labels++;
    }
      
  /* Save the coordinates of the quarter- and three-quarter
     points along the H-bond, too */
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = x1 + 0.25 * (x2 - x1);
      label_coord[n_labels][1] = y1 + 0.25 * (y2 - y1);
      n_labels++;
    }
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = x1 + 0.75 * (x2 - x1);
      label_coord[n_labels][1] = y1 + 0.75 * (y2 - y1);
      n_labels++;
    }

  /* Convert coordinates to PostScript coordinates */
  ps_coord_x = scale_factor * (xc - centre_x) + Plot_Centre_x;
  ps_coord_y = scale_factor * (yc - centre_y) + Plot_Centre_y;
      
  /* Calculate size of rectangle required to blank out a region
     over the line such that can write the H-bond length clearly */
  sizex = 4 * Text_Size_Val->Hbond_Lengths * CHAR_ASPECT / 2.0;
  sizey = 0.6 * Text_Size_Val->Hbond_Lengths;

  /* Blank out the rectangle */
  pssave_();
  pscolb_(Colour->Background);
  pssave_();
  pslwid_("Zero_linewidth");
  pscolb_(Colour->Background);
  psubox_(ps_coord_x - sizex,ps_coord_y + sizey,
          ps_coord_x - sizex,ps_coord_y - sizey,
          ps_coord_x + sizex,ps_coord_y - sizey,
          ps_coord_x + sizex,ps_coord_y + sizey);

  /* Print the H-bond length */
  pscolb_(Colour->Hydrogen_Bonds);
  psctxt_(ps_coord_x,ps_coord_y,Text_Size->Hbond_Lengths,hbond_dist);

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_hbond_lines  -  Print all the H-bond lines
  
***********************************************************************/

void plot_hbond_lines(float centre_x,float centre_y,float scale_factor)
{
  int bond_type, last_hbond;

  float x1, x2, y1, y2;
  float ps_coord_ax, ps_coord_ay, ps_coord_bx, ps_coord_by;
  
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;

  /* Initialise print */
  pscomm_("Hydrogen bond lines");
  pssave_();
  pscolb_(Colour->Hydrogen_Bonds);
  pslwid_(Size->Hydrogen_Bonds);
  psdash_(3);
  last_hbond = TRUE;

  /* Initialise bond pointer to the first bond */
  printf("   Plotting H-bonds ...\n");
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get this bond's type */
      bond_type = bond_ptr->bond_type;

      /* Process only if this is a hydrogen bond */
      if (bond_type == HBOND ||
/* v.4.1.2--> */
/*          (Include->Internal_Hbonds == TRUE && bond_type == HBOND) || */
          (Include->Internal_Hbonds == TRUE && bond_type == INTERNAL) ||
/* <--v.4.1.2 */
          (Include->Hydrophobic_Bonds == TRUE && bond_type == CONTACT))
        {
/* v.4.1.2--> */
          /* If this is an internal H-bond, treat as ordinary H-bond */
          if (bond_type == INTERNAL)
            bond_type = HBOND;
/* <--v.4.1.2 */

          /* Get pointers to the two atoms at either end of the
             bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the atom coordinates */
          x1 = atom1_ptr->x;
          y1 = atom1_ptr->y;
          x2 = atom2_ptr->x;
          y2 = atom2_ptr->y;

          /* Compute the PostScript coordinates of the two ends
             of the bond */
          ps_coord_ax = Plot_Centre_x + scale_factor * (x1 - centre_x);
          ps_coord_ay = Plot_Centre_y + scale_factor * (y1 - centre_y);
          ps_coord_bx = Plot_Centre_x + scale_factor * (x2 - centre_x);
          ps_coord_by = Plot_Centre_y + scale_factor * (y2 - centre_y);

          /* Change colour, if necessary */
          if (bond_type == HBOND && last_hbond == FALSE)
            {
              pscolb_(Colour->Hydrogen_Bonds);
              last_hbond = TRUE;
            }
          else if (bond_type != HBOND && last_hbond == TRUE)
            {
              pscolb_(Colour->Hydrophobics);
              last_hbond = FALSE;
            }

          /* Plot the H-bond */
          psline_(ps_coord_ax,ps_coord_ay,ps_coord_bx,ps_coord_by);
              
          /* If Hbond-lengths are required, then print */
          if (bond_type == HBOND && Include->Hbond_Lengths == TRUE)
/* v.4.1--> */
/*            print_hbond_lengths(bond_ptr->bond_length,x1,y1,x2,y2,
                                scale_factor,centre_x,centre_y); */
            print_hbond_lengths(bond_ptr,x1,y1,x2,y2,
                                scale_factor,centre_x,centre_y);
/* <--v.4.1 */
        }

      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Restore the graphics state */
  psdash_(0);
  psrest_();
}
/***********************************************************************

plot_ext_bond_lines  -  Print any covalent bonds between objects
  
***********************************************************************/

void plot_ext_bond_lines(float centre_x,float centre_y,float scale_factor)
{
  int bond_type;

  float x1, x2, y1, y2;
  float ps_coord_ax, ps_coord_ay, ps_coord_bx, ps_coord_by;
  
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;

  /* Initialise print */
  pscomm_("External bond lines");
  pssave_();
  pscolb_(Colour->External_Bonds);
  pslwid_(Size->External_Bonds);
/* v.3.2--> */
  if (Include->External_Bonds_Solid == FALSE)
/* <--v.3.2 */
    psdash_(3);

  /* Initialise bond pointer to the first bond */
  printf("   Plotting additional bonds ...\n");
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get this bond's type */
      bond_type = bond_ptr->bond_type;

      /* Process only if this is a covalent, elastic bond */
      if (bond_type == COVALENT && bond_ptr->elastic == TRUE)
        {
          /* Get pointers to the two atoms at either end of the
             bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the atom coordinates */
          x1 = atom1_ptr->x;
          y1 = atom1_ptr->y;
          x2 = atom2_ptr->x;
          y2 = atom2_ptr->y;

          /* Compute the PostScript coordinates of the two ends
             of the bond */
          ps_coord_ax = Plot_Centre_x + scale_factor * (x1 - centre_x);
          ps_coord_ay = Plot_Centre_y + scale_factor * (y1 - centre_y);
          ps_coord_bx = Plot_Centre_x + scale_factor * (x2 - centre_x);
          ps_coord_by = Plot_Centre_y + scale_factor * (y2 - centre_y);

          /* Plot the bond */
          psline_(ps_coord_ax,ps_coord_ay,ps_coord_bx,ps_coord_by);

/* v.4.4--> */
          /* Mark bond as plotted */
          bond_ptr->plotted = TRUE;
/* <--v.4.4 */
        }

      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Restore the graphics state */
/* v.3.2--> */
  if (Include->External_Bonds_Solid == FALSE)
/* <--v.3.2 */
    psdash_(0);
  psrest_();
}
/***********************************************************************

get_segment  -  Determine which segment the given bond occupies

***********************************************************************/

void get_segment(float x1,float y1,float x2,float y2,int segment[8],
                 int first_time)
{
  int iseg;

  float angle, dx, dy;

  /* Initialise variables */
  if (first_time == TRUE)
    for (iseg = 0; iseg < 8; iseg++)
      segment[iseg] = FREE;

  /* Determine the direction of the current bond */
  dx =  x2 - x1;
  dy =  y2 - y1;
  if (fabs(dx) > 0.0001)
    angle = atan(fabs(dy) / fabs(dx));
  else
    angle = PI / 2.0;

  /* Convert angle into degrees */
  angle = angle * RADDEG;

  /* Adjust angle to get in right quadrant */
  if (dx < 0.0 && dy >= 0.0)
    angle = 180.0 - angle;
  else if (dx < 0.0 && dy < 0.0)
    angle = 180.0 + angle;
  else if (dx >= 0.0 && dy < 0.0)
    angle = 360.0 - angle;

  /* Determine which quadrant the bond lies in */
  if (angle > 10.0 && angle <= 80.0)
    segment[0] = NOT_FREE;
  else if (angle > 80.0 && angle <= 100.0)
    segment[1] = NOT_FREE;
  else if (angle > 100.0 && angle <= 170.0)
    segment[2] = NOT_FREE;
  else if (angle > 170.0 && angle <= 190.0)
    segment[3] = NOT_FREE;
  else if (angle > 190.0 && angle <= 260.0)
    segment[4] = NOT_FREE;
  else if (angle > 260.0 && angle <= 280.0)
    segment[5] = NOT_FREE;
  else if (angle > 280.0 && angle <= 350.0)
    segment[6] = NOT_FREE;
  else if ((angle > 350.0 && angle <= 360.0) ||
           (angle >= 0.0 && angle <= 10.0))
    segment[7] = NOT_FREE;
}
/***********************************************************************

position_names  -  Check the segment-regions around the current atom to
                   see which does not contain a bond and hence is
                   a suitable site for placement of the atom name

***********************************************************************/

void position_names(struct coordinate *atom_ptr,float *x_pos,float *y_pos)
{
  int first_time;
  int segment[8];

  float x1, x2, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *other_atom_ptr;

  /* Get the current atom's coordinates */
  x1 = atom_ptr->x;
  y1 = atom_ptr->y;
  first_time = TRUE;

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds, picking any involving the current atom */
  while (bond_ptr != NULL)
    {
      /* Process only if this is a covalent bond */
      if (bond_ptr->bond_type == COVALENT)
        {
          /* Initialise pointer to other atom */
          other_atom_ptr = NULL;

          /* Test whether this bond links to the current atom */
          if (bond_ptr->first_atom_ptr == atom_ptr)
            other_atom_ptr = bond_ptr->second_atom_ptr;
          if (bond_ptr->second_atom_ptr == atom_ptr)
            other_atom_ptr = bond_ptr->first_atom_ptr;
           
          /* If bond contains the atom at either end, then process */
          if (other_atom_ptr != NULL)
            {
              /* Get the other atom's coordinates */
              x2 = other_atom_ptr->x;
              y2 = other_atom_ptr->y;

              /* Determine which segment the current bond lies in */
              get_segment(x1,y1,x2,y2,segment,first_time);
              first_time = FALSE;
            }
        }

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Check which segments contain one or more bonds and which are
     free for placement of the atom name */
  if (segment[0] == FREE)
    {
      *x_pos =  1.0;
      *y_pos =  1.0;
    }
  else if (segment[2] == FREE)
    {
      *x_pos = -1.0;
      *y_pos =  1.0;
    }
  else if (segment[4] == FREE)
    {
      *x_pos = -1.0;
      *y_pos = -1.0;
    }
  else if (segment[6] == FREE)
    {
      *x_pos =  1.0;
      *y_pos = -1.0;
    }
  else if (segment[1] == FREE)
    {
      *x_pos =  0.0;
      *y_pos =  1.0;
    }
  else if (segment[3] == FREE)
    {
      *x_pos = -1.0;
      *y_pos =  0.0;
    }
  else if (segment[5] == FREE)
    {
      *x_pos =  0.0;
      *y_pos = -1.0;
    }
  else if (segment[7] == FREE)
    {
      *x_pos =  1.0;
      *y_pos =  0.0;
    }
}
/***********************************************************************

print_current_atom  -  Print out the current atom as a sphere in the
                       PostScript file

***********************************************************************/

void print_current_atom(struct coordinate *atom_ptr,int in_ligand,
                        float centre_x,float centre_y,float scale_factor)
{
  char atom_colour[COL_NAME_LEN + 1];

  int inligand, wanted;

  float ps_coord_ax, ps_coord_ay;

  /* Initialise flag */
  wanted = TRUE;
  inligand = atom_ptr->residue_ptr->inligand;

  /* Check whether this is a main-chain atom */
  if (atom_ptr->side_chain == FALSE)
    {
      /* If plotting the schematic plot, then don't want to plot
         the mainchain atoms (unless carbonyl oxygen) */
      if (Include->Simple_Ligand_Residues == TRUE &&
          inligand == TRUE &&
          (strncmp(atom_ptr->atom_name," O  ",4) &&
           strncmp(atom_ptr->atom_name," OXT",4)))
        wanted = FALSE;

      /* Don't plot atom if XSITE plot and occupancy is zero */
      if (xsite_file == TRUE && strncmp(xsite_probe,"  0.00",6) &&
          !strncmp(atom_ptr->occupancy,"  0.00",6))
        wanted = FALSE;
    }

  /* Check whether atom previously set as non-plottable */
  if (atom_ptr->plot_atom == FALSE)
    wanted = FALSE;

  /* If atom is to be printed, proceed */
  if (wanted == TRUE)
    {
      /* Calculate PostScript coordinates */
      ps_coord_ax = scale_factor * (atom_ptr->x - centre_x)
        + Plot_Centre_x;
      ps_coord_ay = scale_factor * (atom_ptr->y - centre_y)
        + Plot_Centre_y;

      /* Get this atom's colour and set sphere colour accordingly */
      pscolb_(Colour->Atom_Edges);
      get_atom_colour(atom_ptr->atom_name,atom_ptr->residue_ptr->res_name,
                      atom_colour);
      psphcl_(atom_colour);

      /* If this is an XSITE plot and occupancy does not
         correspond to the probe number read in at the start,
         then plot reduced sphere size */
      if (xsite_file == TRUE)
        pspher_(ps_coord_ax,ps_coord_ay,Size->Nonligand_Atoms);
            
      /* If water molecule, set its size */
      else if (!strncmp(atom_ptr->residue_ptr->res_name,"HOH",3))
        pspher_(ps_coord_ax,ps_coord_ay,Size->Waters);

      /* Otherwise, plot the standard sphere size according to 
         whether this is a ligand or non-ligand atom */
      else
        {
          if (in_ligand == TRUE)
            pspher_(ps_coord_ax,ps_coord_ay,Size->Ligand_Atoms);
          else
            pspher_(ps_coord_ax,ps_coord_ay,Size->Nonligand_Atoms);
        }
    }
}
/***********************************************************************

get_atom_text_colour  -  Get the colour of the text for the name of the
                         current atom

***********************************************************************/

void get_atom_text_colour(int in_ligand,char res_name[4],
                          char text_colour[COL_NAME_LEN + 1])
{
  /* Determine the atom-type and its text colour */

  /* Water */
  if (!strncmp(res_name,"HOH",3))
    strcpy(text_colour,Text_Colour->Water_Names);

  /* Ligand atom */
  else if (in_ligand == TRUE)
      strcpy(text_colour,Text_Colour->Ligand_Atom_Names);

  /* Ligand atom */
  else if (in_ligand == FALSE)
      strcpy(text_colour,Text_Colour->Nonligand_Atom_Names);
}
/***********************************************************************

print_atom_name  -  Print the current atom's name

***********************************************************************/

void print_atom_name(struct coordinate *atom_ptr,int in_ligand,
                     float centre_x,float centre_y,float scale_factor)
{
  char text_size[COL_NAME_LEN + 1], text_colour[COL_NAME_LEN + 1];
/* v.4.3--> */
  char text[5];
/* <--v.4.3 */
/* v.4.4.1--> */
  char text_string[5];
/* <--v.4.4.1 */

  int wanted;
/* v.4.3--> */
  int nchar;
/* <--v.4.3 */
/* v.4.4.1--> */
  int fpos, ipos, len;
/* <--v.4.4.1 */

  float atom_radius, x_pos, y_pos, text_height, text_x, text_y;
  float ps_coord_x, ps_coord_y;
/* v.4.3--> */
  float text_size_val, sizex, sizey;
/* <--v.4.3 */

/* v.4.3--> */
  /* Initialise variables */
  nchar = 2;
/* <--v.4.3 */

  /* Get the current atom's radius and text-size */
  atom_radius = atom_ptr->atom_size / scale_factor;
  if (in_ligand == TRUE)
    {
      strcpy(text_size,Text_Size->Ligand_Atom_Names);
      text_height = (Text_Size_Val->Ligand_Atom_Names / scale_factor) / 3.0;
/* v.4.3--> */
      text_size_val = Text_Size_Val->Ligand_Atom_Names;
/* <--v.4.3 */
    }
  else
    {
      strcpy(text_size,Text_Size->Nonligand_Atom_Names);
      text_height = (Text_Size_Val->Nonligand_Atom_Names / scale_factor) / 3.0;
/* v.4.3--> */
      text_size_val = Text_Size_Val->Nonligand_Atom_Names;
/* <--v.4.3 */
    }

  /* Get the atom's radius */
  atom_radius = atom_ptr->atom_size;

/* v.4.3--> */
  /* Transfer text to be printed */
  strcpy(text,atom_ptr->print_name);

  /* If printing element name only, then remove leading blank space, if
     there is one */
  if (Include->Chemical_Notation == TRUE && text[0] == ' ')
    {
      /* Transfer second letter of element name across */
      text[0] = atom_ptr->print_name[1];
      text[1] = '\0';

      /* Blank out carbons */
      if (text[0] == 'C')
        atom_ptr->print_name[0] = '\0';

      /* Set number of characters to 1 */
      nchar = 1;
    }
/* <--v.4.3 */

  /* Set default position for atom label diagonally above and to the
     right of the atom */
  x_pos = y_pos = 1.0;

  /* Determine whether this atom is required for printing */
  wanted = TRUE;
  if (atom_ptr->print_name[0] == '\0')
    wanted = FALSE;

  /* If atom-name is required, then proceed to print it */
  if (wanted == TRUE)
    {
/* v.4.3--> */
      /* If printing in chemical notation, show on top of atom centre */
      if (Include->Chemical_Notation == TRUE)
        {
          /* Make position the atom's centre */
          x_pos = 0.0;
          y_pos = 0.0;
        }

      /* Otherwise place relative to atom's bonds */
      else
        {
/* <--v.4.3 */
          /* Find the best position for the atom name */
          position_names(atom_ptr,&x_pos,&y_pos);
/* v.4.3--> */
        }
/* <--v.4.3 */
      text_x = (atom_radius + text_height) * x_pos;
      text_y = (atom_radius + text_height) * y_pos;

/* v.4.1--> */
      /* Store the atom name's printed location */
      atom_ptr->label_x = atom_ptr->x + text_x;
      atom_ptr->label_y = atom_ptr->y + text_y;
/* <--v.4.1 */

      /* Calculate the PostScript coordinates for name */
      ps_coord_x = scale_factor * (atom_ptr->x - centre_x + text_x)
        + Plot_Centre_x;
      ps_coord_y = scale_factor * (atom_ptr->y - centre_y + text_y)
        + Plot_Centre_y;

/* v.4.3--> */
      /* If printing in chemical notation, then blank out the region
         underneath the element name */
      if (Include->Chemical_Notation == TRUE)
        {
          /* Calculate size of rectangle required to blank out a region
             over the line such that can write the H-bond length clearly */
          sizex = nchar * text_size_val * CHAR_ASPECT / 1.5;
          sizey = 0.5 * text_size_val;

          /* Blank out the rectangle */
          pssave_();
          pslwid_("Zero_linewidth");
          pscolb_(Colour->Background);

          /* If two characters, blank out a rectangular region */
          if (nchar == 2)
            psubox_(ps_coord_x - sizex,ps_coord_y + sizey,
                    ps_coord_x - sizex,ps_coord_y - sizey,
                    ps_coord_x + sizex,ps_coord_y - sizey,
                    ps_coord_x + sizex,ps_coord_y + sizey);

          /* If only a single character, blank out a circular region */
          else
            {
              psccol_(Colour->Background);
              psucir_(ps_coord_x,ps_coord_y,sizey);
              psrest_();
            }

          /* Use the atom colour for the text colour */
          get_atom_colour(atom_ptr->atom_name,
                          atom_ptr->residue_ptr->res_name,
                          text_colour);
        }

      /* Otherwise, use defined text colour */
      else
        {
/* <--v.4.3 */

          /* Print the atom name */
          get_atom_text_colour(in_ligand,atom_ptr->residue_ptr->res_name,
                               text_colour);
/* v.4.3--> */
        }
/* <--v.4.3 */

      /* Set text colour */
      pscolb_(text_colour);
/* v.4.3--> */
/*      psctxt_(ps_coord_x,ps_coord_y,text_size,atom_ptr->print_name); */
/* v.4.4.1--> */

      /* Get the length of the text to be printed */
      len = strlen(text);

      /* Check whether conversion of atom names to Greek characters
         required */
      if (Greek_Names == TRUE)
        {
          /* If first character is blank, then start with second */
          fpos = 0;
          if (text[fpos] == ' ')
            fpos = 1;

          /* Move to appropriate, centred, position on the page */
          psmove_(ps_coord_x,ps_coord_y,text_size,text+fpos);

          /* Print first character in regular font */
          text_string[0] = text[fpos];
          text_string[1] = '\0';
          psfont_(text_size,text_string,3*text_height,FALSE,FALSE);

          /* Loop through remainder of characters printing as Greek
             if a letter, and as half-height if a number */
          for (ipos = fpos + 1; ipos < len; ipos++)
            {
              /* Get the letter */
              text_string[0] = text[ipos];

              /* Check for number */
              if (text_string[0] >= '0' && text_string[0] <= '9')
                {
                  /* Print at half-height */
                  psfont_(text_size,text_string,text_height,FALSE,TRUE);
                }

              /* Otherwise, print as the corresponding Greek character */
              else if (text_string[0] != ' ')
                {
                  /* Convert to lower case */
                  to_lower(text_string,1);

                  /* Print in Greek mode */
                  psfont_(text_size,text_string,text_height,TRUE,FALSE);
                }
            }
        }

      /* Otherwise, print as before */
      else
/* <--v.4.4.1 */
        psctxt_(ps_coord_x,ps_coord_y,text_size,text);
/* <--v.4.3 */
    }
}
/***********************************************************************

plot_atoms  -  Plot all the atoms in the picture

***********************************************************************/

void plot_atoms(float centre_x,float centre_y,float scale_factor)
{
  int iatom, iobject, iresid, natoms, nbonds, nresid;
  int deleted, in_ligand;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr[2];
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;
  struct residue *residue_ptr;
    
  /* Initialise */
  pscomm_("Atoms");
  pssave_();
  pslwid_("Default_linewidth");

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is not a hydrophobic group, and not a simplified
         H-group, then plot all its atoms */
/* v.5.0--> */
      // if (object_ptr->object_type != HYDROPHOBIC &&
      //     object_ptr->object_type != SIMPLE_HGROUP)
      if (object_ptr->object_type != HYDROPHOBIC &&
          object_ptr->object_type != SIMPLE_HGROUP &&
          object_ptr->object_type != DUMMY)
/* <--v.5.0 */
        {
          /* Get pointer to first residue to determine whether this
             object is a ligand or non-ligand one */
          residue_ptr = object_ptr->first_residue_ptr;
          in_ligand = residue_ptr->inligand;

          /* If object is a water, then plot its single atom */
          if (object_ptr->object_type == WATER)
            {
              /* Get its atom pointer */
              atom_ptr[0] = residue_ptr->first_atom_ptr;

              /* Plot the water atom */
/* v.4.0--> */
/*              if (Include->Nonligand_Atoms == TRUE) */
              if (Include->Water_Atoms == TRUE)
/* <--v.4.0 */
                print_current_atom(atom_ptr[0],FALSE,
                                   centre_x,centre_y,scale_factor);
            }

          /* Otherwise, loop through the object's bonds and plot all
             atoms on either end */
          else
            {
              /* Initialise count of bonds encountered */
              nbonds = 0;

              /* Get pointer to the first of this object's bonds */
              object_bond_ptr = object_ptr->first_object_bond_ptr;
      
              /* Loop through all object's bonds */
              while (object_bond_ptr != NULL)
                {
                  /* Get the bond pointer */
                  bond_ptr = object_bond_ptr->bond_ptr;

                  /* Get pointers to the bond's two atoms */
                  atom_ptr[0] = bond_ptr->first_atom_ptr;
                  atom_ptr[1] = bond_ptr->second_atom_ptr;

                  /* Check that neither atom has been deleted */
                  deleted = FALSE;
                  if (atom_ptr[0]->deleted == TRUE ||
                      atom_ptr[1]->deleted == TRUE)
                    deleted = TRUE;
                  else
                    nbonds++;

                  /* Loop to plot the bond's two atoms */
                  for (iatom = 0; iatom < 2 && deleted == FALSE; iatom++)
                    {
                      /* Check that this atom hasn't already been plotted */
                      if (atom_ptr[iatom]->checked == FALSE)
                        {
                          /* Plot the current atom */
                          if ((object_ptr->object_type == LIGAND &&
                               Include->Ligand_Atoms == TRUE) ||
/* v.4.0--> */
                              (object_ptr->object_type == WATER &&
                               Include->Water_Atoms == TRUE) ||
/* <--v.4.0 */
/* v.4.0.1--> */
/*                              Include->Nonligand_Atoms == TRUE) */
                              (object_ptr->object_type != LIGAND &&
                               object_ptr->object_type != WATER &&
                               Include->Nonligand_Atoms == TRUE))
/* <--v.4.0.1 */
                            print_current_atom(atom_ptr[iatom],in_ligand,
                                               centre_x,centre_y,
                                               scale_factor);

                          /* If atom name is required, print it by the atom */
/* v.4.3--> */
/*                          if (Include->Atom_Names == TRUE) */
                          if (Include->Atom_Names == TRUE ||
                              Include->Chemical_Notation == TRUE)
/* <--v.4.3 */
                            print_atom_name(atom_ptr[iatom],in_ligand,
                                            centre_x,centre_y,scale_factor);

                          /* Update flag to indicate atom has already been
                             plotted */
                          atom_ptr[iatom]->checked = TRUE;
                        }
                    }
      
                  /* Get pointer to this object's next bond entry */
                  object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
                }

              /* If the object contains no bonds, loop through all its
                 atoms and print them one by one */
              if (nbonds == 0)
                {
                  /* Get pointer to the first of this object's residues */
                  residue_ptr = object_ptr->first_residue_ptr;
                  nresid = object_ptr->nresidues;
                  iresid = 0;

                  /* Loop over all this object's residues */
                  while (iresid < nresid && residue_ptr != NULL)
                    {
                      /* Get pointer to the first of this residue's atoms */
                      atom_ptr[0] = residue_ptr->first_atom_ptr;
                      iatom = 0;
                      natoms = residue_ptr->natoms;

                      /* Loop over all the residue's atoms to plot */
                      while (iatom < natoms && atom_ptr[0] != NULL)
                        {
                          /* If atom not deleted, then plot */
                          if (atom_ptr[0]->deleted == FALSE)
                            {
                              /* Plot the current atom */
                              if ((object_ptr->object_type == LIGAND &&
                                   Include->Ligand_Atoms == TRUE) ||
/* v.4.0--> */
                                  (object_ptr->object_type == WATER &&
                                   Include->Water_Atoms == TRUE) ||
/* <--v.4.0 */
                                  Include->Nonligand_Atoms == TRUE)
                                print_current_atom(atom_ptr[0],in_ligand,
                                                   centre_x,centre_y,
                                                   scale_factor);

                              /* If atom name is required, print it by
                                 the atom */
/* v.4.3--> */
/*                              if (Include->Atom_Names == TRUE) */
                              if (Include->Atom_Names == TRUE ||
                                  Include->Chemical_Notation == TRUE)
/* <--v.4.3 */
                                print_atom_name(atom_ptr[0],in_ligand,
                                                centre_x,centre_y,
                                                scale_factor);
                            }

                          /* Get pointer to the next atom */
                          atom_ptr[0] = atom_ptr[0]->next;
                          iatom++;
                        }

                      /* Get pointer to the next residue in the list */
                      residue_ptr = residue_ptr->next_residue_ptr;
                      iresid++;
                    }
                }
            }
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

    /* Reset graphics state */
    psrest_();
}
/***********************************************************************

convert_residue_name  -  Convert residue name and number as it is to
                         appear on the plot

***********************************************************************/

void convert_residue_name(char res_name[4],char res_num[6],char chain,
                          char res_name_string[14],int *n_chars)
{
  int i, j, found, number;
  char converted_res_name[4];
  char upper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  char lower[] = "abcdefghijklmnopqrstuvwxyz";
  char digit[] = "0123456789";

  /* Initialise variables */
  *n_chars = 0;

  /* Check whether any digit is a number */
  number = FALSE;
/* v.4.0--> */
  for (i = 0; i < 3; i++)
    {
/* <--v.4.0 */
      for (j = 0; j < 10 && number == FALSE; j++)
        {
/* v.4.0--> */
/*          if (res_name[2] == digit[j]) */
          if (res_name[i] == digit[j])
/* <--v.4.0 */
            number = TRUE;
        }
/* v.4.0--> */
    }
/* <--v.4.0 */

  /* If residue name is a water, or last character is a number,
     or first character is a space, then keep in upper case */
/* v.4.0--> */
/*  if (!strncmp(res_name,"HOH",3) || !strncmp(res_name,"WAT",3) ||
      number == TRUE) */
  if (!strncmp(res_name,"HOH",3) || !strncmp(res_name,"WAT",3) ||
      res_name[0] == ' ' || number == TRUE)
/* <--v.4.0 */
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
  strncpy(res_name_string,converted_res_name,3);
  res_name_string[3] = ' ';
  *n_chars = 4;

  /* Add the residue number */
  if (res_num[2] == ' ')
    {
      strncpy(res_name_string+4,res_num+3,2);
      *n_chars = 6;
    }
  else if (res_num[1] == ' ')
    {
      strncpy(res_name_string+4,res_num+2,3);
      *n_chars = 7;
    }
  else if (res_num[0] == ' ')
    {
      strncpy(res_name_string+4,res_num+1,4);
      *n_chars = 8;
    }
  else
    {
      strncpy(res_name_string+4,res_num,5);
      *n_chars = 9;
    }
  if (res_num[4] == ' ')
    (*n_chars)--;

  /* Add the chain-id */
  if (chain == ' ')
    res_name_string[*n_chars]  = '\0';
  else
    {
      res_name_string[(*n_chars)] = '(';
      res_name_string[(*n_chars) + 1] = chain;
      res_name_string[(*n_chars) + 2] = ')';
      res_name_string[(*n_chars) + 3]  = '\0';
      *n_chars = (*n_chars) + 3;
    }
}
/***********************************************************************

print_name  -  Print residue name

***********************************************************************/

void print_name(char *res_name,char *res_num,char chain,
                float x,float y,float centre_x,float centre_y,
                float text_size,char *text_size_string,
                char *text_colour,float scale_factor,
                int blank_out,int border,int shift_x,int shift_y,
                int stack)
{
  char res_name_string[14], rname[4];

  int n_chars;

  float ps_coord_x, ps_coord_y, ps_y, sizex, sizey;

  /* Calculate PostScript coordinates */
  ps_coord_x = scale_factor * (x - centre_x) + Plot_Centre_x;
  ps_coord_y = scale_factor * (y - centre_y) + Plot_Centre_y;

  /* Prepare string containing residue name, number and chain-id */
  convert_residue_name(res_name,res_num,chain,res_name_string,&n_chars);

  /* Get size of box just bounding residue-name label */
  sizex = 10 * text_size * CHAR_ASPECT / 2.0;
  sizey = 0.6 * text_size;

  /* If shift is required, then apply it */
  ps_coord_x = ps_coord_x - shift_x * sizex;
  ps_coord_y = ps_coord_y - shift_y * sizey;

  /* Set the text colour */
  pscolb_(text_colour);
              
  /* Blank out region for residue name, if required */
  if (blank_out == TRUE)
    {
      pslwid_("Zero_linewidth");
      pshade_(Colour->Background);

      /* Plot a bounded or unbounded box, as required */
      if (border == TRUE)
        psbbox_(ps_coord_x - sizex,ps_coord_y + sizey,
                ps_coord_x - sizex,ps_coord_y - sizey,
                ps_coord_x + sizex,ps_coord_y - sizey,
                ps_coord_x + sizex,ps_coord_y + sizey);
      else
        psubox_(ps_coord_x - sizex,ps_coord_y + sizey,
                ps_coord_x - sizex,ps_coord_y - sizey,
                ps_coord_x + sizex,ps_coord_y - sizey,
                ps_coord_x + sizex,ps_coord_y + sizey);
    }

  /* Print the residue name, either with name stacked above number,
     or in one long string */
  if (stack == TRUE)
    {
      /* Extract just the residue name */
      strncpy(rname,res_name_string,3);
      rname[3] = '\0';

      /* Print the residue name above the number and chain-id */
      ps_y = ps_coord_y + text_size / 2.0;
      psctxt_(ps_coord_x,ps_y,text_size_string,rname);
      ps_y = ps_coord_y - text_size / 2.0;
      psctxt_(ps_coord_x,ps_y,text_size_string,res_name_string + 4);
    }
  else
    psctxt_(ps_coord_x,ps_coord_y,text_size_string,res_name_string);
}
/***********************************************************************

plot_hydrophobics  -  Plot all the symbols representing hydrophobic
                      interactions

***********************************************************************/

void plot_hydrophobics(float centre_x,float centre_y,float scale_factor)
{
/* v.4.0--> */
  char text_colour[PS_STRING_LENGTH];
/* <--v.4.0 */

  int iatom, iobject, jatom, spokes_only;
/* v.4.0--> */
  int interface;
/* <--v.4.0 */

  float radius, xc, x[2], yc, y[2];

  struct bond *bond_ptr;
  struct coordinate *atom_ptr[2];
  struct object *object_ptr;
  struct residue *residue_ptr;
    
  /* Initialise */
  pscomm_("Hydrophobic interactions");
  pssave_();
  pslwid_("Default_linewidth");
  printf("   Plotting hydrophobic contacts ...\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a hydrophobic group, then print residue name */
      if (object_ptr->object_type == HYDROPHOBIC)
        {
          /* Get the coordinates of the object */
          xc = object_ptr->minx;
          yc = object_ptr->miny;

          /* Get pointer to residue so that can pick up residue
             details */
          residue_ptr = object_ptr->first_residue_ptr;

/* v.4.0--> */
          /* Get the appropriate colour for the text of the residue name */
          if (Interface_Plot == FALSE)
            strcpy(text_colour,Text_Colour->Hydrophobic_Names);
          else if (object_ptr->interface == 1)
            strcpy(text_colour,Text_Colour->Ligand_Residue_Names);
          else
            strcpy(text_colour,Text_Colour->Nonligand_Residue_Names);
/* <--v.4.0 */

/* v.4.1--> */
          /* Save the label-position */
          residue_ptr->label_x = xc;
          residue_ptr->label_y = yc;
/* <--v.4.1 */

          /* Plot the residue name at this point */
          pscomm_("Hydrophobic contact");
          print_name(residue_ptr->res_name,residue_ptr->res_num,
                     residue_ptr->chain,xc,yc,centre_x,centre_y,
                     Text_Size_Val->Hydrophobic_Names,
                     Text_Size->Hydrophobic_Names,
/* v.4.0--> */
/*                     Text_Colour->Hydrophobic_Names,scale_factor, */
                     text_colour,scale_factor,
/* <--v.4.0 */
                     TRUE,FALSE,0,0,FALSE);
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Initialise pointer to the first stored bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds to pick out those representing
     hydrophobic contacts */
  while (bond_ptr != NULL)
    {
      /* If this represents a non-bonded hydrophobic contact, then process */
      if (bond_ptr->bond_type == CONTACT)
        {
          /* Get the two atoms involved in the contact */
          atom_ptr[0] = bond_ptr->first_atom_ptr;
          atom_ptr[1] = bond_ptr->second_atom_ptr;

          /* Get the atom coordinates */
          x[0] = atom_ptr[0]->x;
          y[0] = atom_ptr[0]->y;
          x[1] = atom_ptr[1]->x;
          y[1] = atom_ptr[1]->y;

          /* Loop over the two atoms to plot their interactions */
          for (iatom = 0; iatom < 2; iatom++)
            {
              /* Plot atom if not deleted */
              if (atom_ptr[iatom]->deleted == FALSE)
                {
                  /* Get index of the other atom */
                  jatom = 1 - iatom;

                  /* Determine whether this atom belongs to a hydrophobic
                     group or to a ligand/H-bond group */
                  residue_ptr = atom_ptr[iatom]->residue_ptr;
                  object_ptr = residue_ptr->object_ptr;

/* v.4.0--> */
                  /* Initialise interface side */
                  interface = 0;
/* <--v.4.0 */

                  /* Determine radius of arc from which spokes will emanate */
                  spokes_only = TRUE;
                  if (object_ptr->object_type == HYDROPHOBIC)
                    {
                      radius = Size_Val->Hydrophobics;
                      spokes_only = FALSE;

/* v.4.0--> */
                      /* For interface plot, determine which surface this
                         hydrophobic group belongs to */
                      if (Interface_Plot == TRUE)
                        interface = object_ptr->interface;
/* <--v.4.0 */
                    }
                  else
                    radius = scale_factor * atom_ptr[iatom]->atom_size;

/* Old ***          else if (residue_ptr->residue_type == SIMPLE_LIGAND &&
                           !strncmp(atom_ptr[iatom]->atom_name," CA ",4))
                    radius = Size_Val->Simple_Residues;
                  else if (residue_ptr->inligand == TRUE)
                    radius = Size_Val->Ligand_Atoms;
                  else
                    radius = Size_Val->Nonligand_Atoms; */

                  /* Plot the spokes */
                  plot_contact_group(x[iatom],y[iatom],x[jatom],y[jatom],
                                     centre_x,centre_y,scale_factor,
/* v.4.0--> */
/*                                     radius,120.0,spokes_only); */
                                     radius,120.0,spokes_only,interface);
/* <--v.4.0 */
                }
            }
        }

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_simple_ligand_residues  -  Plot all the residues in the simplified
                                representation

***********************************************************************/

void plot_simple_ligand_residues(float centre_x,float centre_y,
                                 float scale_factor)
{
  int iatom, iobject, iresid, natoms, nresid;

  float ps_coord_x, ps_coord_y, xc, yc;

  struct coordinate *atom_ptr, *ca_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;
  char radius[COL_NAME_LEN + 1];
    
  /* Initialise */
  pscomm_("Simplified ligand residues");
  pssave_();
  pslwid_("Default_linewidth");
  printf("   Plotting ligand residues ...\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a ligand, then process its residues */
      if (object_ptr->object_type == LIGAND)
        {
          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* If this residue is to be plotted in the simplified
                 representation, then do so */
              if (residue_ptr->residue_type == SIMPLE_LIGAND)
                {
                  /* Get pointer to the first of this residue's atoms */
                  atom_ptr = residue_ptr->first_atom_ptr;
                  iatom = 0;
                  natoms = residue_ptr->natoms;
                  ca_ptr = NULL;
                  
                  /* Loop over all the residue's atoms to locate
                     the CA */
                  while (iatom < natoms && atom_ptr != NULL &&
                         ca_ptr == NULL)
                    {
                      /* If this is the CA, then save its pointer */
                      if (!strncmp(atom_ptr->atom_name," CA ",4))
                        ca_ptr = atom_ptr;

                      /* Get pointer to the next atom */
                      atom_ptr = atom_ptr->next;
                      iatom++;
                    }

                  /* If have found the CA atom, then plot the residue
                     name at that position */
                  if (ca_ptr != NULL)
                    {
                      /* Get the coordinates of the CA */
                      xc = ca_ptr->x;
                      yc = ca_ptr->y;
                      ps_coord_x = scale_factor * (xc - centre_x)
                        + Plot_Centre_x;
                      ps_coord_y = scale_factor * (yc - centre_y)
                        + Plot_Centre_y;

                      /* Plot the circle representing the residue */
                      strcpy(radius,Size->Simple_Residues);
                      psccol_(Colour->Background);
                      pscolb_(Colour->Simple_Residues);
                      pslwid_("Simple_width");
                      pscirc_(ps_coord_x,ps_coord_y,radius);

                      /* Plot the residue name at this point */
                      pscomm_("Ligand residue");
                      print_name(residue_ptr->res_name,residue_ptr->res_num,
                                 residue_ptr->chain,xc,yc,centre_x,centre_y,
                                 Text_Size_Val->Simple_Residue_Names,
                                 Text_Size->Simple_Residue_Names,
                                 Text_Colour->Ligand_Residue_Names,
                                 scale_factor,FALSE,FALSE,0,0,TRUE);
                    }
                }

              /* Get pointer to the next residue in the list */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_simple_hgroups  -  Plot all the symbols representing the simplified
                        H-bonded residues

***********************************************************************/

void plot_simple_hgroups(float centre_x,float centre_y,float scale_factor)
{
  int iobject;
  int shift_x, shift_y;
  int first_time;
  int diff, iseg, jseg, max_score, segment[8], seg_score[8], sdist;

  float xc, yc;
  float x1, x2, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  static int score[] = { 10, 9, 7, 4, 0 };
    
  /* Initialise */
  pscomm_("Simplified H-groups");
  pssave_();
  pslwid_(Size->Nonligand_Bonds);
  printf("   Plotting simplified H-groups ...\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If object is a hydrophobic group, then plot it */
      if (object_ptr->object_type == SIMPLE_HGROUP)
        {
          /* Get the coordinates of the object */
          xc = object_ptr->minx;
          yc = object_ptr->miny;

          /* Get pointer to residue so that can pick up residue
             details */
          residue_ptr = object_ptr->first_residue_ptr;

          /* Determine best position of residue label so that it
             doesn't interfere with the bonds coming from it */
          shift_x = 1;
          shift_y = 0;
          first_time = TRUE;

          /* Get pointer to first bond */
          bond_ptr = first_bond_ptr;

          /* Loop through all the hydrogen bonds involving this H-group
             and determine the direction in which each bond lies */
          while (bond_ptr != NULL)
            {
              /* If this is a hydrogen bond, then get its two atoms */
              if (bond_ptr->bond_type == HBOND)
                {
                  /* Get the bond's two atoms */
                  atom1_ptr = bond_ptr->first_atom_ptr;
                  atom2_ptr = bond_ptr->second_atom_ptr;

                  /* If either atom involves the current object, then
                     get direction of the bond */
                  if (atom1_ptr->residue_ptr->object_ptr == object_ptr ||
                      atom2_ptr->residue_ptr->object_ptr == object_ptr)
                    {
                      /* Get the two atoms' coordinates */
                      if (atom1_ptr->residue_ptr->object_ptr == object_ptr)
                        {
                          x1 = atom1_ptr->x;
                          y1 = atom1_ptr->y;
                          x2 = atom2_ptr->x;
                          y2 = atom2_ptr->y;
                        }
                      else
                        {
                          x1 = atom2_ptr->x;
                          y1 = atom2_ptr->y;
                          x2 = atom1_ptr->x;
                          y2 = atom1_ptr->y;
                        }

                      /* Determine which segment the current bond lies in */
                      get_segment(x1,y1,x2,y2,segment,first_time);
                      first_time = FALSE;
                    }
                }

              /* Go to the next bond in the linked list */
              bond_ptr = bond_ptr->next_bond_ptr;
            }

          /* Initialise the segment scores */
          for (iseg = 0; iseg < 8; iseg++)
            seg_score[iseg] = 0;

          /* Loop through to see which segments are free and which are
             taken, and accumulate the scores accordingly */
          for (iseg = 0; iseg < 8; iseg++)
            {
              /* If this segment is taken, score the nearest segments
                 most highly */
              if (segment[iseg] == NOT_FREE)
                {
                  /* Loop through all segment positions, scoring
                     each one according to its proximity to the
                     current segment */
                  for (jseg = 0; jseg < 8; jseg++)
                    {
                      /* Get the distance between the segments and
                         score accordingly */
                      diff = jseg - iseg;
                      if (diff > 4)
                        sdist = 8 - diff;
                      else if (diff < -4)
                        sdist = 8 + diff;
                      else if (diff < 0)
                        sdist = - diff;
                      else
                        sdist = diff;

                      /* Compute the score */
                      seg_score[jseg] = seg_score[jseg] + score[sdist];
                    }
                }
            }

          /* Determine which is the maximum-scoring position */
          jseg = 0;
          max_score = 0;
          for (iseg = 0; iseg < 8; iseg++)
            {
              if (seg_score[iseg] > max_score)
                {
                  max_score = seg_score[iseg];
                  jseg = iseg;
                }
            }

          /* Convert the optimal segment position into x- and y-shifts */
          if (jseg == 0)
            {
              shift_x = 1;
              shift_y = 1;
            }
          else if (jseg == 1)
            {
              shift_x =  0;
              shift_y =  1;
            }
          else if (jseg == 2)
            {
              shift_x = -1;
              shift_y =  1;
            }
          else if (jseg == 3)
            {
              shift_x = -1;
              shift_y =  0;
            }
          else if (jseg == 4)
            {
              shift_x = -1;
              shift_y = -1;
            }
          else if (jseg == 5)
            {
              shift_x =  0;
              shift_y = -1;
            }
          else if (jseg == 6)
            {
              shift_x =  1;
              shift_y = -1;
            }
          else if (jseg == 7)
            {
              shift_x =  1;
              shift_y =  0;
            }
          
          /* Print the residue name, with a border around it */
          pscomm_("H-group");
          print_name(residue_ptr->res_name,residue_ptr->res_num,
                     residue_ptr->chain,xc,yc,centre_x,centre_y,
                     Text_Size_Val->Nonligand_Residue_Names,
                     Text_Size->Nonligand_Residue_Names,
                     Text_Colour->Nonligand_Residue_Names,scale_factor,
                     TRUE,TRUE,shift_x,shift_y,FALSE);
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

extract_search_atoms  -  Extract all the atoms that are within range
                         of the given range-limits

***********************************************************************/

void extract_search_atoms(struct residue *current_residue_ptr,
                          struct coordinate **first_stack_ptr,
                          float search_min_x,float search_max_x,
                          float search_min_y,float search_max_y,
                          float scale_factor)
{
  int current, iatom, natoms, wanted;
/* v.4.0--> */
  int plotted_accessibility;
/* <--v.4.0 */
/* v.5.3.6 --> */
  int object_type;
/* <-- v.5.3.6 */

  float acc_size, atom_size, x, y;

  struct coordinate *atom_ptr, *last_atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  acc_size = 4.0 * Size_Val->Ligand_Atoms / scale_factor;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Initialise stack pointers */
  last_atom_ptr = NULL;
  *first_stack_ptr = NULL;

  /* Initialise pointer to the first stored residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all residues */
  while (residue_ptr != NULL)
    {
      /* Determine whether this residue is in range of the search region */
      wanted = TRUE;
      if (residue_ptr->minx > search_max_x ||
          residue_ptr->maxx < search_min_x ||
          residue_ptr->miny > search_max_y ||
          residue_ptr->maxy < search_min_y)
        wanted = FALSE;

      /* If atoms from this residue likely to be wanted, then loop through
         them */
      if (wanted == TRUE)
        {
/* v.5.3.6 --> */
          /* Get the object type */
          object_type = residue_ptr->object_ptr->object_type;
/* <-- v.5.3.6 */

          /* If this is the same residue as the current one, then set
             flag */
          current = FALSE;
          if (residue_ptr == current_residue_ptr)
            current = TRUE;

          /* Get pointer to this residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all the residue's atoms, performing the
             transformation */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* Get this atom's coordinates and its size */
              x = atom_ptr->x;
              y = atom_ptr->y;
              atom_size = atom_ptr->atom_size;

/* v.4.0--> */
              /* Determine whether accessibility will have been plotted
                 for this atom */
              plotted_accessibility = FALSE;
              if (Include->Accessibilities == TRUE)
                {
                  if (residue_ptr->inligand == TRUE ||
                      Include->Ligand_Accessibilities_Only == FALSE)
                    plotted_accessibility = TRUE;
                }
/* <--v.4.0 */

              /* If have plotted accessibility shading, need to increase
                 atom's size (if not already done) */
/* v.4.0--> */
/*              if (Include->Accessibilities == TRUE) */
              if (plotted_accessibility == TRUE)
/* <--v.4.0 */
                {
                  if (atom_size < (acc_size))
                    {
                      atom_size = acc_size;
                      atom_ptr->atom_size = atom_size;
                    }
                }

              /* Determine whether this atom is within range */
              wanted = TRUE;
              if ((x - atom_size > search_max_x) ||
                  (x + atom_size < search_min_x) ||
                  (y - atom_size > search_max_y) ||
                  (y + atom_size < search_min_y))
                wanted = FALSE;

/* v.5.3.6 --> */
              /* If the residue is a hydrophobic or simple group, only
                 want a single representative atom for the energy
                 calculations */
              if ((object_type == HYDROPHOBIC ||
                   object_type == SIMPLE_HGROUP) &&
                  atom_ptr != residue_ptr->first_atom_ptr)
                wanted = FALSE;
/* <-- v.5.3.6 */

              /* If atom is within range, then add it to the stack */
              if (wanted == TRUE)
                {
                  /* If this is the first atom, get the first pointer
                     to point to the second atom */
                  if (*first_stack_ptr == NULL)
                    *first_stack_ptr = atom_ptr;

                  /* Otherwise, set the last atom's pointer to the
                     other atom */
                  else
                    last_atom_ptr->next_stack_ptr = atom_ptr;

                  /* Save pointer to the current atom */
                  last_atom_ptr = atom_ptr;

                  /* If this atom belongs to the current residue,
                     flag it as checked */
                  if (current == TRUE)
                    atom_ptr->checked = TRUE;
                }

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

label_dist  -  Calculate (squared) distance of point from residue-label
               position

***********************************************************************/

float label_dist(float x, float y, float label_coordx, float label_coordy,
                 float zone_x, float zone_y)
{
  int i, j;
  float corner_x, corner_y, distance_sqrd, dist_min, dist_x, dist_y;

  /* Check if the atom is overlapping any part of the label area */
  if (x < label_coordx + zone_x && x > label_coordx - zone_x && 
      y < label_coordy + zone_y && y > label_coordy - zone_y)
      distance_sqrd = 0.0;

  /* If the atom's x-coordinate falls between the
     corner x-coordinates, then use the straight-line
     y-distance */
  else if (x < label_coordx + zone_x && x > label_coordx - zone_x)
    {
      distance_sqrd = fabs(y - (label_coordy + zone_y));
      dist_y = fabs(y - (label_coordy - zone_y));
      if (dist_y < distance_sqrd)
        distance_sqrd = dist_y;
      distance_sqrd = distance_sqrd * distance_sqrd;
    }
  
  /* If the atom's y-coordinate falls between the
     corner y-coordinates, then use the straight-line
     x-distance */
  else if (y < label_coordy + zone_y && y > label_coordy - zone_y)
    {
      distance_sqrd = fabs(x - (label_coordx + zone_x));
      dist_x = fabs(x - (label_coordx - zone_x));
      if (dist_x < distance_sqrd)
        distance_sqrd = dist_x;
      distance_sqrd = distance_sqrd * distance_sqrd;
    }
  
  /* Otherwise, calculate the distance of the atom
     from the nearest corner */
  else
    {
      dist_min = 100000.0;

      /* Loop across the four corners */
      for (i = 0; i < 2; i++)
        {
          corner_x =  label_coordx + (i * 2 - 1) * zone_x;
          for (j = 0; j < 2; j++)
            {
              corner_y =  label_coordy + (j * 2 - 1) * zone_y;
              
              /* Calculate the distance of the atom from
                 this corner */
              dist_x = x - corner_x;
              dist_y = y - corner_y;
              distance_sqrd = dist_x * dist_x + dist_y * dist_y;
              
              /* Save this distance if it is the minimum
                 so far */
              if (distance_sqrd < dist_min)
                dist_min = distance_sqrd;
            }
        }
      distance_sqrd = dist_min;
    }
  
  /* Compute the effective energy from the closest approach */
  if (distance_sqrd < 0.00001)
    distance_sqrd = 0.00001;

  return(distance_sqrd);
}
/***********************************************************************

perform_search  -  Search for a suitable location for the placement of
                   the residue label, avoiding overlaps with atoms and
                   other labels as far as possible

***********************************************************************/

void perform_search(struct coordinate *first_stack_ptr,float grid_size,
                    float box_min_x,float box_max_x,
                    float box_min_y,float box_max_y,
                    float half_width,float half_height,
                    float *pos_x,float *pos_y)
{
  int ilabel;
  int npoint_x, npoint_y, point_x, point_y;
  int no_good;

  float x, y;
  float atom_size, length_x, length_y;
  float label_coordx, label_coordy;
  float energy, energy_diff, energy_store, label_energy;
  float lowest_energy, highest_energy, worst_energy, worst_other;
  float distance, distance_sqrd, closest_distance_sqrd;

  struct coordinate *atom_ptr;

  /* Initialise variables */
  lowest_energy = 100000.0;
  highest_energy = 0.0;
  energy_store = 0.0;

  /* Work out the dimensions of the box in which the residue name can
     be placed */
  length_x = (box_max_x - box_min_x);
  length_y = (box_max_y - box_min_y);
  *pos_x = box_min_x;
  *pos_y = box_min_y;

  /* Work out number of grid-points required for trial placement of
     residue name */
  npoint_x = 2 + length_x / grid_size;
  npoint_y = 2 + length_y / grid_size;

  /* Put the residue name at each grid-point and calculate an "energy"
     for it, depending on how close it is to its residue and whether
     it clashes with any other atoms */

  /* Loop over the grid-points in the y-direction */
  for (point_y = 0; point_y < npoint_y; point_y++)
    {
      label_coordy = box_min_y + point_y * grid_size;

      /* Loop over the grid-points in the x-direction */
      for (point_x = 0; point_x < npoint_x; point_x++)
        {
          label_coordx = box_min_x + point_x * grid_size; 
          worst_energy = 0.0;
          label_energy = 0.0;
          worst_other = 0.0;
/*          closest_atom_ptr = NULL; */
          closest_distance_sqrd = 1000000.0;
          no_good = FALSE;

          /* Get pointer to the first atom on the stack of stored atoms
             within range */
          atom_ptr = first_stack_ptr;

          /* Loop through all the stored atoms in the vicinity of
             the current residue */
          while (atom_ptr != NULL && no_good == FALSE)
            {
              /* Get the coordinates and size of this atom */
              x = atom_ptr->x;
              y = atom_ptr->y;
              atom_size = atom_ptr->atom_size;

              /* Calculate the nearest distance of atom from
                 any part of the residue label */
              distance_sqrd = label_dist(x,y,label_coordx,label_coordy,
                                         half_width,half_height);

              /* If atom overlaps the label, then discard this position */
/* v.5.3.6 --> */
//              if (distance_sqrd < (atom_size + 0.005) * (atom_size + 0.005))
              if (distance_sqrd <
                  (atom_size + LABEL_GAP) * (atom_size + LABEL_GAP))
/* <-- v.5.3.6 */
                {
                  no_good = TRUE;
                  energy_store = 10000000.0;
                  distance_sqrd = 0.0;
                }

                /* Otherwise, compute the effective energy from the closest
                   approach */
                else
                  {
                    distance = sqrt((double)distance_sqrd) - atom_size;
                    distance_sqrd = distance * distance;
                    energy_store = (1.0 / distance_sqrd);
                  }

              /* If this is the closest atom so far, save it */
              if (distance_sqrd < closest_distance_sqrd)
                {
/*                    closest_atom_ptr = atom_ptr; */
                  closest_distance_sqrd = distance_sqrd;
                }

              /* If the current atom belongs to the correct residue
                 save the energy only if it is the highest so far */
              if (atom_ptr->checked == TRUE)
                {
                  /* Save the worst energy encountered so far */
/* v.5.3.6 --> */
//                    energy = 10.0 * (IDEAL_ENERGY - energy_store);
                  energy = distance_sqrd;
                  if (no_good == TRUE)
                    energy = 10000000.0;
//                    if (energy < 0)
//                      energy = 10000000.0;
/* <-- v.5.3.6 */
                  worst_energy = worst_energy + energy;
                }

              /* Otherwise, save the worst energy from close approaches
                 to other residues */
              else
                worst_other = worst_other + energy_store;

              /* Get pointer to the next atom in the stack */
              atom_ptr = atom_ptr->next_stack_ptr;
            }

          /* Loop through all the labels written so far (ie H-bond
             distances and residue labels done to date) */
          for (ilabel = 0; ilabel < n_labels && no_good == FALSE; ilabel++)
            {
              /* Get the coordinates of this label */
              x = label_coord[ilabel][0];
              y = label_coord[ilabel][1];

              /* Calculate the nearest distance of this point from
                 any part of the residue label */
              distance = label_dist(x,y,label_coordx,label_coordy,
                                    half_width,half_height);
              if (distance < 0.0001)
                {
                  no_good = TRUE;
                  energy_store = 10000000.0;
                }

              /* Compute the effective energy from the closest
                 aproach */
              else
                energy_store = (1.0 / distance);

              /* Save the worst label-energy encountered so far */
              label_energy = label_energy + energy_store;
            }

          /* Check difference between the highest energy (ie coming from
             the atom that's closest to the label) and the "ideal"
             energy */
          energy_diff = worst_energy + worst_other + 10.0 * label_energy;
          if (energy_diff < 0.0)
            energy_diff = 10000000.0;

          /* If this is the very first loop, set current
             energy as the lowest so far encountered */
          if (point_x == 0 && point_y == 0)
            lowest_energy = energy_diff + 1.0;

          /* Now pick out the coordinates of the lowest energy */
          if (energy_diff < lowest_energy)
            {
              *pos_x = label_coordx;
              *pos_y = label_coordy;
              lowest_energy = energy_diff;
            }
          else if (energy_diff < 1000000.0 &&
                   energy_diff > highest_energy)
            highest_energy = energy_diff;
        }
    }

  /* Save the coords actually used */
  if (n_labels < MAXLABELS)
    {
      label_coord[n_labels][0] = *pos_x;
      label_coord[n_labels][1] = *pos_y;
      n_labels++;
    }
}
/***********************************************************************

resname_grid_search  -  Move the residue name so it doesn't overlap with
                        the atoms

***********************************************************************/

void resname_grid_search(struct residue *current_residue_ptr,
                         float text_size,int n_chars,
                         float *new_coordx,float *new_coordy,
                         float min_x,float min_y,
                         float max_x,float max_y,
                         float scale_factor)
{
  float grid_size;
  float box_min_x, box_max_x, box_min_y, box_max_y;
  float search_min_x, search_max_x, search_min_y, search_max_y;
  float x, y;
  float half_width, half_height, zone_x, zone_y;

  struct coordinate *first_stack_ptr;

  /* Initialise variables */
  grid_size = 0.3;
  first_stack_ptr = NULL;

  /* Compute width of residue name label, given its size and number
     of characters */
  half_width = n_chars * text_size * CHAR_ASPECT / 2.0;
  half_height = text_size / 2.0;
/*  upper_width = 15 * text_size * CHAR_ASPECT / 2.0; */

  /* Calculate the size of the exclusion zone for the residue label */
  zone_x = half_width + Max_object_size / scale_factor;
  zone_y = half_height + Max_object_size / scale_factor;

  /* Define the minimum and maximum coordinates of the box within which
     the residue name can be placed and of the search region for
     atoms */
  box_min_x = min_x - 2.0 * half_width;
  box_max_x = max_x + 2.0 * half_width;
  box_min_y = min_y - 4.0 * half_height;
  box_max_y = max_y + 4.0 * half_height;
/* v.5.3.6--> */
  if (Interface_Plot == TRUE)
    {
      box_min_y = min_y - 8.0 * half_height;
      box_max_y = max_y + 8.0 * half_height;
    }
/* <-- v.5.3.6 */
  search_min_x = box_min_x - zone_x;
  search_max_x = box_max_x + zone_x;
  search_min_y = box_min_y - zone_y;
  search_max_y = box_max_y + zone_y;

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Extract all the atoms that are within the search box */
  extract_search_atoms(current_residue_ptr,&first_stack_ptr,
                       search_min_x,search_max_x,search_min_y,
                       search_max_y,scale_factor);
  /* Perform the search for a suitable location for the residue name */
  perform_search(first_stack_ptr,grid_size,box_min_x,box_max_x,box_min_y,
                 box_max_y,half_width,half_height,&x,&y);

  /* Return the calculated coordinates for the residue label */
  *new_coordx = x;
  *new_coordy = y;
}
/***********************************************************************

print_object  -  Print current object (for debugging purposes)

***********************************************************************/

void print_object(struct object *object_ptr)
{
  int atom_no, iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Print object type */
  printf("Object type %d, number of residues = %d\n",
         object_ptr->object_type,object_ptr->nresidues);
  
  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;
  atom_no = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, printing each one */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Print atom details */
          printf("ATOM  %5d %s %s %c%s   %8.3f%8.3f%8.3f%s%s %7.3f\n",
                 atom_no + 1,
                 atom_ptr->atom_name,
                 atom_ptr->residue_ptr->res_name,
                 atom_ptr->residue_ptr->chain,
                 atom_ptr->residue_ptr->res_num,
                 atom_ptr->x,
                 atom_ptr->y,
                 atom_ptr->z,
                 atom_ptr->occupancy,
                 atom_ptr->bvalue,
                 atom_ptr->accessibility);

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
          atom_no++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

print_residue_names  -  Print all the residue names for the ligand and
                        nonligand sidechains

***********************************************************************/

void print_residue_names(float centre_x,float centre_y,float scale_factor)
{
  char chain, res_name[4], res_num[6], res_name_string[14];
  char text_size[COL_NAME_LEN + 1];

  int iobject, iresid, n_chars, nresid;
  int object_type;

  float min_x, min_y, max_x, max_y, new_coordx, new_coordy;
  float ps_coord_x, ps_coord_y, text_size_val;
  
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise print */
  pscomm_("Residue names");
  pssave_();

  /* Initialise pointer to the first stored object */
  printf("   Printing residue names ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get this object's type */
      object_type = object_ptr->object_type;

      /* Process only if not a hydrophobic group */
      if (object_type != HYDROPHOBIC && object_type != SIMPLE_HGROUP)
        {
          /* Set default text colour and size according to the
             object type */
          if (object_type == LIGAND)
            {
              pscolb_(Text_Colour->Ligand_Residue_Names);
              strcpy(text_size,Text_Size->Ligand_Residue_Names);
              text_size_val = Text_Size_Val->Ligand_Residue_Names
                / scale_factor;
            }

          /* Default text colour and size for waters */
          else if (object_type == WATER)
            {
/* v.4.0--> */
/*              pscolb_(Text_Size->Water_Names); */
              pscolb_(Text_Colour->Water_Names);
/* <--v.4.0 */
              strcpy(text_size,Text_Size->Water_Names);
              text_size_val = Text_Size_Val->Water_Names / scale_factor;
            }

          /* Default text colour and size for non-ligand bonds */
          else
            {
              pscolb_(Text_Colour->Nonligand_Residue_Names);
              strcpy(text_size,Text_Size->Nonligand_Residue_Names);
              text_size_val = Text_Size_Val->Nonligand_Residue_Names
                / scale_factor;
            }

          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop through all this object's residues, one by one */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Precess only if this is a standard residue (ie not a
                 schematic in the simplified ligand representation */
              if (residue_ptr->residue_type == STANDARD)
                {
                  /* Get the residue name, number and chain ID */
                  strncpy(res_name,residue_ptr->res_name,3);
                  res_name[3] = '\0';
                  strncpy(res_num,residue_ptr->res_num,5);
                  res_num[5] = '\0';
                  chain = residue_ptr->chain;

                  /* Prepare string containing residue name, number
                     and chain-id */
                  convert_residue_name(res_name,res_num,chain,
                                       res_name_string,&n_chars);

                  /* Get the minimum and maximum coordinates of the atoms
                     belonging to this residue */
                  min_x = residue_ptr->minx;
                  max_x = residue_ptr->maxx;
                  min_y = residue_ptr->miny;
                  max_y = residue_ptr->maxy;

                  /* If already have the residue label positions, then
                     use the stored values */
/* v.4.3--> */
                  if (Plot_from_drw == TRUE)
                    {
                      /* Get the label positions */
                      new_coordx = residue_ptr->label_x;
                      new_coordy = residue_ptr->label_y;
                    }

                  /* Otherwise, calculate the positions of the residue
                     names */
                  else
                    {
/* <--v.4.3 */
                      /* Perform the grid-search for the optimal positioning
                         of the residue name */
                      resname_grid_search(residue_ptr,text_size_val,n_chars,
                                          &new_coordx,&new_coordy,
                                          min_x,min_y,max_x,max_y,
                                          scale_factor);
/* v.4.3--> */
                    }
/* <--v.4.3 */

/* v.4.1--> */
                  /* Save the label-position */
                  residue_ptr->label_x = new_coordx;
                  residue_ptr->label_y = new_coordy;
/* <--v.4.1 */

                  /* Calculate PostScript coords for residue label and text
                     height */
                  ps_coord_x = scale_factor * (new_coordx - centre_x)
                    + Plot_Centre_x;
                  ps_coord_y = scale_factor * (new_coordy - centre_y)
                    + Plot_Centre_y;

                  /* Print out string containing residue name and number */
                  psctxt_(ps_coord_x,ps_coord_y,text_size,res_name_string);
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Reset graphics state */
  psrest_();
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

void extract_residue_name(char *res_name,char *res_num,char *chn,
                          int *vtype)
{
  char chain, number_string[20], tmp_string[20];

  char *string_ptr;

  int name_end, var_type;

  /* Initialise variables */
  res_name[0] = res_num[0] = '\0';
  chain = ' ';
  var_type = NO_VARIANT;
  
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

  /* Get the variant type being plotted */
  strcpy(number_string,highlight_res + name_end + 1);
  string_ptr = strstr(number_string,".");
  if (string_ptr != NULL)
    var_type = atoi(string_ptr + 1);
  if (var_type == -9)
    var_type = NO_VARIANT;

  /* Return the chain id and variant type */
  *chn = chain;
  *vtype = var_type;
}
/***********************************************************************

get_other_chains  -  Read in the chain equivalences from the chains.dat
                   file and return the other chain in the interface if
                   it is equivalent to the one with the highlighted
                   residue

***********************************************************************/

int get_other_chains(char chain,char *chains)
{
  char file_name[FILENAME_LEN + 1], input_line[LINELEN + 1];
  char chain_id, copy_of, other_chain, store_chain;

  int ichain, len, nchains, nlines;

  FILE *file_ptr;

  /* Initialise variables */
  nchains = 0;
  nlines = 0;
  other_chain = '\0';

  /* Store our chain */
  chains[nchains] = chain;
  nchains++;

  /* Form the name of the chains.dat file in the PDBsum directory */
  strcpy(file_name,pdbsum_dir);
  strcat(file_name,"chains.dat");

  /* Open the chains.dat input file */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Loop while reading in records from the input file */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Increment count of lines read in */
          nlines++;

          /* Truncate the string to strip off the newline */
          len = string_truncate(input_line,LINELEN);

          /* If this is a copy of our chain, store it */
          if (len > 2)
            {
              /* Get the two chain ids */
              chain_id = input_line[0];
              copy_of = input_line[2];
              store_chain = '\0';

              /* If either chain corresponds to our chain, store
                 the other one as one of the equivalents */
              if (chain_id == chain)
                store_chain = copy_of;
              else if (copy_of == chain)
                store_chain = chain_id;

              /* If have a chain to store, do so */
              if (store_chain != '\0' && nchains < MAX_CHAINS)
                {
                  /* Store this chain */
                  chains[nchains] = store_chain;

                  /* Increment count of stored chains */
                  nchains++;
                }
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }

  /* Show warning tnat chain.dat file not found */
  else
    printf("*** Warning. File not found: %s\n",file_name);

  /* Return number of chains */
  return(nchains);
}
/***********************************************************************

get_variants  -  Read in the DisaStr veriants on the current LIGPLOT
                 diagram

***********************************************************************/

int get_variants(struct highres **fst_highres_ptr,
                 struct highres **lst_highres_ptr,int nhigh,
                 char *chains,int nchains)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char chain, res_name[4], res_num[6];
  char number_string[20], tmp[LINELEN + 1];

  char *string_ptr;

  int ichain, len, line, type;
  int ref, wanted;

  struct highres *highres_ptr, *first_highres_ptr, *last_highres_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  line = 0;
  ref = FALSE;
  type = NO_VARIANT;
  wanted = FALSE;

  /* Initialise pointers */
  first_highres_ptr = *fst_highres_ptr;
  last_highres_ptr = *lst_highres_ptr;
  highres_ptr = NULL;

  /* Form name of input annotations.dat file */
  strcpy(file_name,wkdir);
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
          type = NO_VARIANT;
          wanted = FALSE;
        }

      /* If interaction type, save the type */
      else if (!strncmp(input_line,"INT_VAR_TYPE ",13))
        {
          /* Get the variant type */
          if (string_ptr != NULL)
            {
              /* Get the type */
              strcpy(number_string,string_ptr + 1);
              type = atoi(number_string);
              if (type == -9)
                type = NO_VARIANT;
            }
        }

      /* If this is the start of a new residue, then re-initialise */
      else if (!strncmp(input_line,"INT_PDB_CODE[",13))
        {
          /* Reinitialise flag */
          wanted = FALSE;

          /* Check if this is the right PDB code */
          if (string_ptr != NULL && !strcmp(string_ptr + 1,pdb_code))
            wanted = TRUE;
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

          /* If this is the LIGPLOT identifer, save the residue */
          else if (!strncmp(input_line,"INT_LIGPLOT[",12))
            {
              /* See if the current chain is one of the wanted ones */
              wanted = FALSE;
              for (ichain = 0; ichain < nchains; ichain++)
                {
                  /* If the right chain, set the flag */
                  if (chains[ichain] == chain)
                    wanted = TRUE;
                }
              
              /* Loop over all the equivalent chains, to create a
                 highres record for each */
              for (ichain = 0; ichain < nchains && wanted == TRUE;
                   ichain++)
                {
                  /* See if we already have this residue */
                  highres_ptr
                    = find_highres_entry(first_highres_ptr,res_name,
                                         res_num,chains[ichain]);

                  /* If not the reference residue, then create a record
                     for it */
                  if (highres_ptr == NULL)
                    {
                      /* Create record */
                      highres_ptr
                        = create_highres_entry(&first_highres_ptr,
                                               &last_highres_ptr,type,ref,
                                               res_name,res_num,
                                               chains[ichain]);

                      /* Increment count of variants stored */
                      nhigh++;
                    }

                  /* If this is the reference residue, then save the type */
                  else if (highres_ptr->ref == TRUE &&
                           highres_ptr->type == NO_VARIANT)
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

get_rot_matrix  -  Calculate the rotation matrix for a rotation about the
                   z-axis through the given angle

***********************************************************************/

void get_rot_matrix(float angle,float matrix[2][2])
{
  double theta;

  /* Convert angle into radians */
  theta = angle / RADDEG;

  /* Calculate rotation matrix */
  matrix[0][0] = cos(theta);
  matrix[0][1] = sin(theta);
  matrix[1][0] = -sin(theta);
  matrix[1][1] = cos(theta);
}
/***********************************************************************

apply_atom_rotation  -  Apply the given rotation to the atoms

***********************************************************************/

void apply_atom_rotation(float *coords[2],int iatom,float matrix[2][2])
{
  float x, y;

  /* Get the coordinates of this atom */
  x = coords[0][iatom];
  y = coords[1][iatom];

  /* Apply the rotation */
  coords[0][iatom] = (x * matrix[0][0] + y * matrix[0][1]);
  coords[1][iatom] = (x * matrix[1][0] + y * matrix[1][1]);
}
/***********************************************************************

get_best_ellipse  -  Calculate optimal orientation of bounding ellipse

***********************************************************************/

float get_best_ellipse(struct coordinate *first_atom_ptr,int natoms,
                       float *centre,int nstored,float *ellipse_radius)
{
  int iatom, istep, istored, nsteps;

  float angle, cmin[2], cmax[2], matrix[2][2], theta;
  float height, max_ratio, ratio, width;

  float *coords[2];

  struct coordinate *atom_ptr;

  /* Initialise */
  iatom = istored = 0;
  max_ratio = theta = 0;

  /* Create two arrays for the coordinates */
  coords[0] = (float *) malloc(nstored * sizeof(float));
  coords[1] = (float *) malloc(nstored * sizeof(float));

  /* Get residue's atom */
  atom_ptr = first_atom_ptr;

  /* Loop over all the residue's atoms to plot */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* If atom not deleted, then process */
      if (atom_ptr->deleted == FALSE)
        {
          /* Stored the coords, adjusted to the origin */
          coords[0][istored] = atom_ptr->x - centre[0];
          coords[1][istored] = atom_ptr->y - centre[1];

          /* Increment count of stored atoms */
          istored++;
        }
                          
      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
      iatom++;
    }

  /* Calculate number of angle steps to try */
  nsteps = (int) (180.0 / ANGLE_STEP);
  angle = ANGLE_STEP;

  /* Calculate the rotation matrix for the given angle step */
  get_rot_matrix(angle,matrix);

  /* Loop over all the rotation steps */
  for (istep = 0; istep < nsteps; istep++)
    {
      /* Loop over ths stored atom coords to get the max and min values */
      for (iatom = 0; iatom < nstored; iatom++)
        {
          /* Apply rotation to this atom's coords */
          apply_atom_rotation(coords,iatom,matrix);

          /* If first atom, store the maximum and minimum */
          if (iatom == 0)
            {
              cmin[0] = cmax[0] = coords[0][iatom];
              cmin[1] = cmax[1] = coords[1][iatom];
            }
          else
            {
              if (coords[0][iatom] < cmin[0])
                cmin[0] = coords[0][iatom];
              if (coords[0][iatom] > cmax[0])
                cmax[0] = coords[0][iatom];
              if (coords[1][iatom] < cmin[1])
                cmin[1] = coords[1][iatom];
              if (coords[1][iatom] > cmax[1])
                cmax[1] = coords[1][iatom];
            }
        }

      /* Calculate effective ellipse diameters */
      width = cmax[0] - cmin[0] + 2 * ELLIPSE_MARGIN;
      height = cmax[1] - cmin[1] + 2 * ELLIPSE_MARGIN;

      /* Calculate x:y ratio */
      ratio = width / height;

      /* If this is the highest ratio so far, store the angle and
         diameters */
      if (ratio > max_ratio)
        {
          /* Store this angle */
          if (angle > 0)
            theta = -angle;
          else
            theta = 0;

          /* Store the ellipse diameters */
          ellipse_radius[0] = width / 2;
          ellipse_radius[1] = height / 2;

          /* Store the maximum ratio */
          max_ratio = ratio;
        }

      /* Increment the angle */
      angle = angle + ANGLE_STEP;
    }

  /* Return the optimal rotation angle */
  return(theta);
}
/***********************************************************************

calc_ellipse  -  Calculate the ellipse parameters for highlighting the
                 given residue

***********************************************************************/

float calc_ellipse(struct residue *residue_ptr,int type,float *centre,
                   float *ellipse_radius,float scale_factor)
{
  int iatom, natoms, nstored;

  float max_radius, radius, theta;
  
  struct coordinate *atom_ptr;

  /* Initialise */
  centre[0] = centre[1] = 0;
  max_radius = 0;
  nstored = 0;
  theta = 0;

  /* Get residue's atom */
  atom_ptr = residue_ptr->first_atom_ptr;
  iatom = 0;
  natoms = residue_ptr->natoms;

  /* Loop over all the residue's atoms to plot */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* If atom not deleted, then plot */
      if (atom_ptr->deleted == FALSE)
        {
          /* Update centre of mass accumulator */
          centre[0] = centre[0] + atom_ptr->x;
          centre[1] = centre[1] + atom_ptr->y;
          nstored++;

          /* Get this atom's rtadius */
          if (type == HYDROPHOBIC)
            radius = Size_Val->Hydrophobics / scale_factor;
          else
            radius = atom_ptr->atom_size;

          /* Update maximum radius */
          if (radius > max_radius)
            max_radius = radius;
        }
                          
      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
      iatom++;
    }

  /* Calculate the centre of mass */
  if (nstored > 0)
    {
      centre[0] = centre[0] / nstored;
      centre[1] = centre[1] / nstored;
    }

  /* If we only have one atom, or the residue is a hydrophobic gourp,
     then define ellipse as a circle */
  if (nstored < 2 || type == HYDROPHOBIC)
    ellipse_radius[0] = ellipse_radius[1] = max_radius + ELLIPSE_MARGIN;

  /* Otherwise, compute the best ellipse to surround the atoms */
  else
    theta
      = get_best_ellipse(residue_ptr->first_atom_ptr,natoms,centre,
                         nstored,ellipse_radius);
    
  /* Return the ellipse angle */
  return(theta);
}
/***********************************************************************

plot_ellipse  -  Plot the ellipse around the highlighted residue

***********************************************************************/

void plot_ellipse(float centre_x,float centre_y,float scale_factor,
                  float *centre,float *ellipse_radius,float angle,
                  struct highres *highres_ptr)
{
  char colour[60], comment[60], number_string[30];

  int type;

  float line_width, ratio, width, x, y;

  /* Initialise variables */
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
  pscomm_(comment);
  pssave_();

  /* Define ellipse's line width */
  line_width = scale_factor * ELLIPSE_LINEWIDTH;
  sprintf(number_string,"%8.2f",line_width);
  if (line_width < NEAR_ZERO)
    strcpy(number_string,"Default_linewidth");
  pslwid_(number_string);

  /* Set the line colour */
  sprintf(colour,"Ellipse_col_%s",ELLIPSE_COL_DEF[type + 2]);
  pscolb_(colour);

  /* If not the reference residue, make the ellipse dashed */
  if (highres_ptr->ref == FALSE)
    psdash_(5);

  /* Calculate the height to width ratio */
  ratio = ellipse_radius[1] / ellipse_radius[0];

  /* Convert the ellipse-coordinates to PostScript coords */
  x = scale_factor * (centre[0] - centre_x) + Plot_Centre_x;
  y = scale_factor * (centre[1] - centre_y) + Plot_Centre_y;

  /* Convert the ellipse width to PostScript equivalent */
  width = scale_factor * ellipse_radius[0];

  /* Draw the ellipse */
  pselip_(x,y,width,ratio,angle);
  
  /* If not the reference residue, unset dashed lines */
  if (highres_ptr->ref == FALSE)
    psdash_(0);

  /* Reset graphics state */
  psrest_();
}
/***********************************************************************

plot_highlights  -  Draw an ellipse around the residue to be highlighted

***********************************************************************/

void plot_highlights(float centre_x,float centre_y,float scale_factor)
{
  char chain, chains[MAX_CHAINS], res_name[4], res_num[6];

  int ichain, nchains, nhigh, type, var_type;
  int ref;

  float centre[2], ellipse_radius[2], theta;
  
  struct object *object_ptr;
  struct residue *residue_ptr;
  struct highres *highres_ptr, *first_highres_ptr, *last_highres_ptr;
    
  /* Initialise variables */
  ellipse_radius[0] = ellipse_radius[1] = theta = 0;
  nhigh = 0;
  ref = TRUE;
  type = NO_VARIANT;

  /* Initialise pointers */
  highres_ptr = first_highres_ptr = last_highres_ptr = NULL;
 
  /* Unpack the name of the highlighted residue */
  extract_residue_name(res_name,res_num,&chain,&var_type);

  /* Pick up the chains from the chains.dat file, if present */
  nchains = get_other_chains(chain,chains);

  /* Loop over all the equivalent chains, to create a highres record
     for each */
  for (ichain = 0; ichain < nchains; ichain++)
    {
      /* Create a highlight record for this residue */
      highres_ptr
        = create_highres_entry(&first_highres_ptr,&last_highres_ptr,type,
                               ref,res_name,res_num,chains[ichain]);
      nhigh++;

      /* Save the variant type */
      highres_ptr->type = var_type;
    }

  /* Get all the other residues to be highlighted */
  nhigh = get_variants(&first_highres_ptr,&last_highres_ptr,nhigh,
                       chains,nchains);

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If not a dummy object, see if it is the residue to be highlighted */
      if (object_ptr->object_type != DUMMY)
        {
          /* Get pointer to first residue to determine whether this
             object is a ligand or non-ligand one */
          residue_ptr = object_ptr->first_residue_ptr;

          /* See if this residue is to be highlighted */
          highres_ptr
            = find_highres_entry(first_highres_ptr,residue_ptr->res_name,
                                 residue_ptr->res_num,
                                 residue_ptr->chain);

          /* If it is a variant residue, work out parameters of ellipse
             to be plotted */
          if (highres_ptr != NULL && 
              (highres_ptr->type != NO_VARIANT || highres_ptr->ref == TRUE))
            {
              /* Calculate the parameters of the ellipse for this
                 residue */
              theta = calc_ellipse(residue_ptr,
                                   object_ptr->object_type,centre,
                                   ellipse_radius,scale_factor);

              /* Plot the ellipse */
              plot_ellipse(centre_x,centre_y,scale_factor,centre,
                           ellipse_radius,theta,highres_ptr);
            }
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* <--v.5.3.8 */
/***********************************************************************

plot_postscript_picture  -  Plot the PostScript picture

***********************************************************************/

void plot_postscript_picture(float centre_x, float centre_y,
                             float scale_factor)
{
  /* If accessibility shading required, do it */
  if (Include->Accessibilities == TRUE)
    shade_accessibilities(centre_x,centre_y,scale_factor);

/* v.5.2 --> */
  /* If showing LigSearch annotations, plot them for all ligand atoms */
  if (ligsearch_annot == TRUE)
    ligsearch_annotations(centre_x,centre_y,scale_factor);
/* <-- v.5.2 */

  /* Draw all the H-bond lines */
  plot_hbond_lines(centre_x,centre_y,scale_factor);
    
  /* Draw any additional covalent bonds */
  if (Include->External_Bonds == TRUE)
    plot_ext_bond_lines(centre_x,centre_y,scale_factor);
    
  /* Draw all the covalent bond lines */
  plot_bond_lines(centre_x,centre_y,scale_factor);
    
/* v.4.4--> */
  /* Plot circles within any rings identified as aromatic */
  if (Include->Double_Bonds == TRUE)
    plot_aromatics(centre_x,centre_y,scale_factor);
/* <--v.4.4 */

  /* Plot the atoms */
  plot_atoms(centre_x,centre_y,scale_factor);

  /* Plot the hydrophobic symbols */
  plot_hydrophobics(centre_x,centre_y,scale_factor);

  /* If this is a schematic ligand plot, plot the individual components */
  if (Include->Simple_Ligand_Residues == TRUE)
    {
      /* Print the residue names for the ligand */
      plot_simple_ligand_residues(centre_x,centre_y,scale_factor);
    }

  /* If the non-ligand residues are schematic, plot them */
  if (Include->Simple_Nonligand_Residues == TRUE)
    {
      /* Print the residue names for the H-groups */
      plot_simple_hgroups(centre_x,centre_y,scale_factor);
    }

  /* Print the residue names, if required */
/* v.4.3--> */
  if (Include->Residue_Names == TRUE)
/* <--v.4.3 */
    print_residue_names(centre_x,centre_y,scale_factor);

/* v.5.3.8--> */
  /* If we have a highlighted residue, draw ellipse around it */

  if (highlight_res[0] != '\0')
    plot_highlights(centre_x,centre_y,scale_factor);
/* <--v.5.3.8 */
}

