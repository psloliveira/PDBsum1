/***********************************************************************

  ligxmasx.c - Dummy version of ligxmas.c for compilation without the
               full XMAS routines

************************************************************************

Date:-         10 June 1999

Written by:-   Roman Laskowski

               Department of Crystallography
               Birkbeck College
               Malet Street, London WC1E 7HX, UK.


Version:-      Original version was part of v.4.2 of LIGPLOT.

Amendments:-   Amendments after v.4.2 are labelled by version number
               v.m.n--> and v.m.n<-- where m.n is the version
               number corresponding to the change

----------------------------------------------------------------------*/

/* Dummy Functions
   ---------------

have_xmas_link
read_in_xmas_file
get_nlig_xmas_molecs
get_xmaslig_start_end
form_plot_title
get_xmas2pdb_line
get_xmas_atom
get_xmas_conect
get_xmas2hb_line
get_xmas_hbdata
get_next_xmas_ligand
free_up_memory

*/

typedef int XMAS;

/***********************************************************************

have_xmas_link  -  Return flag indicating that only dummy functions
                   are present

***********************************************************************/

int have_xmas_link(void)
{
  return(FALSE);
}
/***********************************************************************

read_in_xmas_file  -  Read in the required atom coordinates from the
                      XMAS file

***********************************************************************/

XMAS *read_in_xmas_file(char xmas_name[FILENAME_LEN])
{
  XMAS *xmas;

  xmas = NULL;

  return(xmas);
}
/***********************************************************************

get_xmas_atom  -  Read in the next atom record from the cached XMAS file

***********************************************************************/

int get_xmas_atom(XMAS *xmas,int *atom_num,char atom_type[5],char res_name[4],
		  char res_num[6],char *chain,float *x,float *y,float *z,
		  float *occup,float *bvalue,int *id)
{
  int end_of_file = FALSE;

  return(end_of_file);
}
/***********************************************************************

get_xmas_conect  -  Read in the next CONECT record from the cached
                    XMAS file

***********************************************************************/

int get_xmas_conect(XMAS *xmas,int *atom_num1,int *atom_num2)
{
  int end_of_file = FALSE;
  return(end_of_file);
}
/***********************************************************************

get_nlig_xmas_molecs  -  Get the number of ligand molecules in the XMAS
                         file, and the number for the first one

***********************************************************************/

int get_nlig_xmas_molecs(XMAS *xmas,int *first_lig_molec)
{
  return(0);
}
/***********************************************************************

get_xmaslig_start_end  -  Given the molecule-id, get the residue range
                          this represents

***********************************************************************/

int get_xmaslig_start_end(XMAS *xmas,char res_name1[4],char res_num1[6],
			  char res_name2[4],char res_num2[6],
			  char *chain_identi)
{
  return(0);
}
/***********************************************************************

get_xmas2pdb_line  -  Return the next ATOM/HETATM record from the XMAS
                      cache in PDB format

***********************************************************************/

char *get_xmas2pdb_line(XMAS *xmas,char line[LINELEN + 1],int first,
/* v.4.2.2--> */
/*			int rec_types[NTYPES],int *current_type) */
			int rec_types[NTYPES],int *current_type,
			int *mol_id)
/* <--v.4.2.2 */
{
  return(NULL);
}
/***********************************************************************

get_xmas_hbdata  -  Read in the next H-bond or non-bonded contact
                    from the cached XMAS file

***********************************************************************/

int get_xmas_hbdata(XMAS *xmas,int file_type,int *interaction_type,
		    char atom_type1[5],
		    char res_name1[4], char res_num1[6],char *chain1,
		    char atom_type2[5],
		    char res_name2[4], char res_num2[6],char *chain2)
{
  int end_of_file = FALSE;

  return(end_of_file);
}
/***********************************************************************

form_plot_title  -  Form the plot title from PDB code and ligand residue
                    range

***********************************************************************/

int form_plot_title(XMAS *xmas,char pdb_code[5],char res_name1[4],
		    char res_num1[6], char res_name2[4],char res_num2[6],
		    char chain_identi)
{
  return(0);
}
/***********************************************************************

get_xmas2hb_line  -  Return the next H-bond or non-bonded contact record
                     from the XMAS cache in HBPLUS format

***********************************************************************/

char *get_xmas2hb_line(XMAS *xmas,char line[LINELEN + 1],int first,
		       int file_type)
{
  return(NULL);
}
/***********************************************************************

get_next_xmas_ligand  -  Get the number of next ligand molecule in the
                         XMAS file

***********************************************************************/

int get_next_xmas_ligand(XMAS *xmas,int current_lig_molec)
{
  return(0);
}
/* v.4.3--> */
/***********************************************************************

free_up_memory  -  Free up the memory used by the structures in the
                   last run of LIGPLOT, in preparation for plotting
                   the next ligand

***********************************************************************/

void free_up_memory(void)
{
  struct hhb_info *hhb_info_ptr, *next_hhb_info_ptr;
  struct object *object_ptr, *next_object_ptr;
  struct residue *residue_ptr, *next_residue_ptr;
  struct coordinate *atom_ptr, *next_atom_ptr;
  struct bond *bond_ptr, *next_bond_ptr;

  /* HBPLUS info data */

  /* Set pointer to the first of the links */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Search through all the stored items */
  while (hhb_info_ptr != NULL)
    {
      /* Get pointer next one */
      next_hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;

      /* Free the memory */
      free(hhb_info_ptr);

      /* Go to next one */
      hhb_info_ptr = next_hhb_info_ptr;
    }

  /* Object data */

  /* Set pointer to the first of the links */
  object_ptr = first_object_ptr;

  /* Search through all the stored items */
  while (object_ptr != NULL)
    {
      /* Get pointer next one */
      next_object_ptr = object_ptr->next_object_ptr;

      /* Free the memory */
      free(object_ptr);

      /* Go to next one */
      object_ptr = next_object_ptr;
    }

  /* Residue data */

  /* Set pointer to the first of the links */
  residue_ptr = first_residue_ptr;

  /* Search through all the stored items */
  while (residue_ptr != NULL)
    {
      /* Get pointer next one */
      next_residue_ptr = residue_ptr->next_residue_ptr;

      /* Free the memory */
      free(residue_ptr);

      /* Go to next one */
      residue_ptr = next_residue_ptr;
    }

  /* Atom data */

  /* Set pointer to the first of the links */
  atom_ptr = first_atom_ptr;

  /* Search through all the stored items */
  while (atom_ptr != NULL)
    {
      /* Get pointer next one */
      next_atom_ptr = atom_ptr->next;

      /* Free the memory */
      free(atom_ptr);

      /* Go to next one */
      atom_ptr = next_atom_ptr;
    }

  /* Bond data */

  /* Set pointer to the first of the links */
  bond_ptr = first_bond_ptr;

  /* Search through all the stored items */
  while (bond_ptr != NULL)
    {
      /* Get pointer next one */
      next_bond_ptr = bond_ptr->next_bond_ptr;

      /* Free the memory */
      free(bond_ptr);

      /* Go to next one */
      bond_ptr = next_bond_ptr;
    }
}
/* <--v.4.3 */
