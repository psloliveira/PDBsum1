/***********************************************************************

  ligplot.c - Program to produce a schematic 2D plot of a given protein
              ligand, showing the protein residues to which it is
              hydrogen-bonded and those with which it makes hydrophobic
              contact

************************************************************************

Version:-      v.5.4.6  -  13 Sep 2024

Written by:-   Andrew Wallace & Roman Laskowski

               Department of Biochemistry and Molecular Biology
               University College London
               Gower Street, London WC1E 6BT, UK.

Additional
routines:-     qikfit and matfit - least-squares fitting routines by
(v.3.0)        Andrew Martin, Dept. Biochem. & Mol. Biol., University
               College London.

<               Eigen value calculation routines by Ian Tickle & David
               Moss, Dept. Crystallography, Birkbeck College.

Reference:-    Wallace A C, Laskowski R A & Thornton J M .
               LIGPLOT: A program to generate schematic diagrams
               of protein-ligand interactions. Prot. Eng., 8,
               127-134.

Further info:- http://www.biochem.ucl.ac.uk/bsm/ligplot/ligplot.html

Requirements:- The list of hydrogen bonds and non-bonded contacts
               required by the program can either be user-generated
               in the format described in the Operating Manual, or
               produced by the HBPLUS program, written by Ian McDonald.
               To obtain a copy of the latter program, see

               http://www.biochem.ucl.ac.uk/~mcdonald/hbplus/home.html

               The ligand atoms can also be shaded (coloured) according
               to their solvent accessibility. Accessibilities can be
               calculated using the NACCESS program, written by
               Simon Hubbard. To obtain a copy, see

               http://sjh.bi.umist.ac.uk/naccess.html

Amendments:-   Amendments after v.4.0 are labelled by version number
               v.m.n--> and v.m.n<-- where m.n is the version
               number corresponding to the change

v.1.0          Original release version.                 September 1994
                                     (Andrew Wallace & Roman Laskowski)

v.2.0          Revised input parameters, including many new parameters,
               particularly those defining the colours and sizes of the
               different items and labels plotted. Portrait/landscape
               option added. Name of input parameter file changed from
               ligplot.prm to ligplot.par to prevent confusion with
               former version.
               Changes to PostScript output routines to make the output
               PostScript file simpler to understand and edit.
               Several bugs in v.1.0 found and fixed.
               Double- and triple-bond option added (these have to be
               manually defined, though).
               Read-in and use of CONECT records in input PDB file.
               Output of ligplot.frm file giving original coordinates,
               in PDB format, of all the residues shown on the LIGPLOT
               diagram.
                                                         Jan 1996 (RAL)

v.2.0.1        Addition of routine to check for hydrogen atoms so that
               they can be discarded.                  2 Feb 1996 (RAL)

               Replacement of "connect" array by the bond structure
               to improve program's storage efficiency. New routines:
               update_bond, get_bond_type.
               Extra option (only used when parameter Include_Hydrogens
               set to TRUE and program recompiled) to include hydrogen
               atoms on the plot (Note: corresponding H-bonds have to be
               defined to H-atoms in the .hhb file).
                                                      26 Feb 1996 (RAL)

v.2.1          Bug-fix to incorrect conversion of atom-names from ligplot-
               hbond file format to .hhb hbond file format.
               Addition of read_conec_records and verify_conec_records
               routines to pick up all the non-ligand residues that
               are covalently bonded to the ligand (as per the CONECT
               records), and to store these bonds as pseudo H-honds.
                                                   11-14 Nov 1996 (RAL)
               Routines to store all ligand residue-names and numbers,
               as encountered in get_ligand_start_end, to compare against
               when reading in the H-bonds.
                                                      21 Nov 1996 (RAL)
v.3.0          Major revisions, including change of ligplot.par back
               to ligplot.prm.
                                             29 Jan - 31 Mar 1997 (RAL)
v.3.0.1        Bug-fix to routine determining whether line read in from
               HBPLUS file is valid. Lines with blanks at the start were
               ignored, but HBPLUS fills line-starts with blanks if it
               cannot identify a PDB 4-letter code in the name of the
               input PDB file.
                                                       8 Apr 1997 (RAL)
v.3.1          Extension of combine_objects routine to ligand objects,
               specifically to handle phosphotyrosines in the ligand. In
               the process, also fixed a bug in the join_residues function
               which goes wrong if block of residues or atoms to be moved
               are already adjacent in the list.
                                                      11 Apr 1997 (RAL)
               Amendments to hbadd.c program.
                                                      13 Apr 1997 (RAL)
               Improvements to some of the most computationally
               intensive routines to speed up whole program.
                                                      14 Apr 1997 (RAL)
               Check for change in residue name as well as number in
               read_pdb_file function, otherwise H-bond info not picked
               up for HET group attachements to residues
               (eg phosphotyrosines)
                                                      15 Apr 1997 (RAL)
               Improved program's ability to cope with duplicate
               residue numbering in PDB files.
                                                      16 Apr 1997 (RAL)
v.3.1.1        Bug-fix for certain cases where cyclic peptides identified
               (eg covalently-linked sugars in 1gya). New function to
               check for such cases prior to the split_object function
               so that the different units are split into separate
               objects for ease of flattening.
               Check for ENDMDL added to function get_ligand_start_end
               so that warning messages not repeated unnecessarily.
                                                      20 Apr 1997 (RAL)
v.3.1.2        Amendment primarily for PDBsum. Where ligand residues
               are duplicated and last "apparent" ligand is consists of
               a single atom only, then take the "true" ligand to be the
               penultimate one, rather than the last.
                                                      22 Apr 1997 (RAL)
               Speeded up processing of CONECT records for NMR
               structures having large numbers of them.
               Additional flag to specify that ligand is a metal.
                                                      23 Apr 1997 (RAL)
v.3.2          Bug-fix in stretch_mainchain function causing crash when
               no carbonyl oxygen present.
                                                      28 Apr 1997 (RAL)
               Addition of graph-matching routines to hbadd.c.
                                                       1 May 1997 (RAL)
               Bug-fix. When two residues were combined into a hydrophobic
               group, and then one of the residues' atoms all deleted as
               unwanted, gave problems in the plot with the residue name
               displaced from the "eye-lash" representing it.
               "Hydrophobic contacts" extended to include non-bonded
               contact between any pair of atoms where at least one of
               the pair is a carbon or sulphur.
                                                       2 May 1997 (RAL)
               Addition of two new parameters to the ligplot.prm file,
               one dealing with the definition of non-bonded contacts and
               the other with a final rotation of the picture.
               Removal of one of the warning messages.
                                                       8 May 1997 (RAL)
               Amendment to distance-energy for non-bonded contacts.
               Amendment to read-in of CONECT records to check for early
               line-end (as in output from the Babel program).
                                                      12 May 1997 (RAL)
               Amendment to initially place all hydrophobic interactions
               a long way from the ligand to give the H-bonded groups a
               better chance of finding better positions.
               New option to show covalent bonds to external groups as
               solid lines rather than as dotted lines.
                                                      13 May 1997 (RAL)
               Bug-fix. Residues belonging to deleted objects now marked
               as deleted to prevent them being written out to the
               ligplot.frm file.
               Bug-fix. Problems in transformation routines when ligand
               consisting of a single atom only.
                                                      15 May 1997 (RAL)
               New output file, ligplot.res, listing all the residues
               included in the plot (for use by the raslig program in
               the PDBsum LIGPLOT pages).
               Improvements to function that interprets the command-line
               parameters.
                                                      16 May 1997 (RAL)
               Bug-fix to lig_search routines which were failing to find
               and display the first "cofactor".
                                                      19 May 1997 (RAL)
               Bug-fix to remove any bonds between hydrophobic groups 
               and any non-ligand residues.
                                                      25 May 1997 (RAL)
               Amendment to cope with real problem cases where squashing
               flat of the ligand is required (eg as in 2cmm). Function
               align_object added to orientate the ligand as best as
               possible before it is brutally forced flat.
               Added check for Hg and Ho in hydrogen-check.
                                                      28 May 1997 (RAL)
v.3.2.1        Addition of marker to identify the ligand residues in the
               output ligplot.res file.
                                                      19 Jun 1997 (RAL)
v.4.0          Addition of routines to read in interactions details across
               a dimer or domain interface, as generated by Alice
               Wuensche's program, dimer.c
                                             28 Jun - 10 Aug 1997 (RAL)
               Addition of output ligplot.rcm, which holds the residue
               centres of mass for use in generating similar LIGPLOT
               layouts for similar input structures. Routines to read in
               the ligplot.rcm file and use its data for the minimization
               procedure.
                                                      10 Aug 1997 (RAL)
               Amendments to hbadd (rewrite of clique-detection routine).
                                                   16-18 Sep 1997 (RAL)
               Addition of check that maximum accessibility non-zero,
               otherwise led to divide by zero. Spotted by Dave Renouf.
                                                      21 Sep 1997 (RAL)
               Amendment of .rcm file routines to have centres of mass
               of residues, rather than object, and to allow restraint
               of individual atoms also.
                                                      16 Oct 1997 (RAL)
               Increase number of atoms stored from 1200 to 5000 as was
               not enough for some structures (eg NAD 336 of 1gd1).
               Fix to exclude CONECT records giving covalent bonds
               between non-ligand residues (as in 1gd1).
               Minor fixes to ligplot.res file (only used by PDBsum).
                                                      10 Nov 1997 (RAL)
               Additional command-line input of residue names, as well as
               residue numbers, for defining the range of the ligand's
               residues. Particularly useful where PDB file contains
               multiple residues having exactly the same residue number.
               Amendment to prevent protein residues being incorporated
               into the ligand if they have the same residue number and
               chain ID as a residue in the ligand.
                                                      11 Nov 1997 (RAL)
               Addition of list of covalent bonds between ligand and
               protein residues to .res file.
                                                      12 Nov 1997 (RAL)
               Addition of option to parameters to allow water atoms to
               be shown as spheres, or not.
                                                      16 Nov 1997 (RAL)
               Storage of maximum accessibility value in the ligplot.pdb
               file so that shading preserved when plot regenerated.
                                                      18 Nov 1997 (RAL)
v.4.0.1        Bug-fix. Ligand atoms were being printed even if option
               not to print them had been selected.
                                                       3 Dec 1997 (RAL)
v.4.0.2        Bug-fix to deal with cases where ligand residues are
               numbered backwards (as in 9lpr).
               Bug-fix to increase maximum number of command-line
               parameters that can be stored (MAXTOKENS) as giving
               problems with certain structures (eg MG 435 -> HAD 438
               in 1gin).
                                                       8 Dec 1997 (RAL)
v.4.1          Addition of output to ligplot.drw file for input into
               the Java LIGPLOT editor, LigEd.
                                                      28 Jul 1998 (RAL)
               Discard covalent links from CONECT records between ligand
               and non-ligand residues if the plot external bonds option
               is "N".
                                                       3 Sep 1998 (RAL)
v.4.1.1        More metals added to list in metal-check.
                                                      22 Sep 1998 (RAL)
               Fix to ligand-search routine to cope with cases where
               residue name contains leading blanks.
                                                       1 Oct 1998 (RAL)
               Addition of summary interactions file ligplot.sum (for
               internal use only).
                                                      19 Oct 1998 (RAL)
v.4.1.2        Addition of debug flag to help locate errors.
                                                      10 Dec 1998 (RAL)
               Bug-fix. Internal H-bonds not being shown, even when
               selected.
                                                      16 Dec 1998 (RAL)
               Bug-fix. Internal H-bonds not being printed as normal
               H-bonds.
                                                      28 Mar 1999 (RAL)
               Bug-fix (Malcolm Gillies). If first atom record in the
               file is a water, program goes into infinite loop.
                                                      14 May 1999 (RAL)
               Bug-fix (Malcolm Gillies). Amendment to read-in of
               CONECT records in read_pdb_file to check for early
               line-end.
v.4.2          Addition of routines to read in Andrew Martin's XMAS
               file formats containing the coordinates, H-bonds and
               van der Waals contacts.
                                                    4-11 Jun 1999 (RAL)
               Split off of PostScript routines to separate file,
               ligps.c.
                                                      11 Jun 1999 (RAL)
v.4.2.1        Bug-fix (Andrew Martin). Amendment put into v.4.2 in 
               calling read_hbplus_file contained error in parameters
               when called for reading non-bonded contacts.
                                                      11 Jun 1999 (RAL)
               Bug-fix (Martyn Pearce). One of the loops in the 
               principal_components function going off end of 
               transformation array. Possibly cause of LIGPLOT hanging
               in print_residue_names function when compiled with -O2
               under linux.
                                                      22 Sep 1999 (RAL)
v.4.2.2        Fixes to XMAS-format reads to deal with ligands
               consisting of several "residues" which are interspersed
               in the original PDB file with non-ligand residues.
                                                       3 Mar 2000 (RAL)
               Fix to errors in array definitions in parameters for
               function calculate_matrix.
                                                      13 Jun 2000 (RAL)

v.4.3          Amendments to read in bond order info from hbadd.bonds
               file created by the modified version of hbadd.c using
               the new CIF-format Het Group Dictionary.
                                                      17 Nov 2000 (RAL)
               Additional "C" option in the ligplot.prm file for
               generating plots in "chemical notation".

               New check_for_metal routine.
                                                       1 Dec 2000 (RAL)

               Modification for dealing with new .asa format generated
               by NACCESS.
                                                      14 Dec 2000 (RAL)
               Minor fixes to deal with crashes under linux.
                                                      20 Dec 2000 (RAL)
               Modifications to -clip line for Landscape orientation.
               Addition of write_ligand_fit_coords routine for PDBsum
               generation run.
                                                       8 Feb 2001 (RAL)
v.4.4          Modification to get_ligand_start_end algorithm such that
               if exact range is supplied, then all residues within the
               range are ligand residues. (This to overcome aberrant
               residue-numberings such as that of the residues of chain
               R in 1a2c).
                                                      16 Feb 2001 (RAL)
               Check for molecules interpenetrating ligand so as to
               discard them.
                                                       6 Mar 2001 (RAL)
               Check for residues names that are shorter than 3
               characters long. Halt on too many ligand residues. If
               ligand residue is a metal and is "covalently" bonded to
               other ligand residues, make these bonds elastic.
                                                      12 Mar 2001 (RAL)
               Modification to "chemical notation plot" to keep atoms
               during energy minimization, but drop them right
               afterwards.
                                                      13 Mar 2001 (RAL)
               Addition of circles within aromatic rings.
                                                      20 Mar 2001 (RAL)
               Reinstatement of old double-bond print for old cases
               where chemical notation not used.
                                                       9 May 2001 (RAL)
v.4.4.1        Bug-fix to dimer.c causing it to fail under latest SGI
               Origin compiler. Spotted by Atsushi Nakagawa.
                                                      10 Sep 2001 (RAL)
               Bug-fix on the .frm and .rcm files being closed a second
               time.
               Addition of -hyd run-time option to include hydrogens,
               and -gk option to convert atom names to Greek letters.
                                                      20 Feb 2002 (RAL)
v.4.4.2        Addition of -crop line to ligplot.ps output for newer
               versions of convert.
                                                      10 Dec 2002 (RAL)
v.4.4.3        ligplot.sum file being closed twice and causing core
               dump on the second close, so first close removed.
                                                      17 Apr 2003 (RAL)
v.4.4.4        Bounding box altered to be equivalent to clipping region.
                                                      12 Jan 2004 (RAL)
               Bug-fix on Print_as_is option, crashing when trying to
               close unopened files,
                                                      10 Feb 2004 (RAL)
               New user-definable scaling option (only for use in
               generating plots for PDBsum and the Enzyme Structures
               Database.
                                                      25 Mar 2004 (RAL)
               Call to get_element to deal with atom names read in from
               the ligplot.drw file.
                                                      23 Jul 2004 (RAL)
v.4.5          Addition of dimhtlm to dimplot and dimonly scripts to
               give new schematic diagrams of interactions across
               chain or domain interfaces.
                                                       6 Feb 2006 (RAL)
               Aromatic bond type written out to ligplot.drw file and
               read in again (for use with PDBsum).
                                                       6 Feb 2006 (RAL)
v.4.5.1        Bug-fix. Schematic residue-plot was showing atoms rather
               than boxes when coords read in from ligplot.pdb file.
                                                       2 Jun 2006 (RAL)
v.4.5.2        Bug-fix. Where have aromatic rings as part of a joined
               set of rings (eg 1y5d and ) the circle representing
               the aromatic rings was being incorrectly drawn.
                                                      14 Sep 2006 (RAL)
v.4.5.3        Commenting out of free_up_memory function to prevent
               glibc error in the newer versions of linux.
                                                      17 Apr 2007 (RAL)
v.4.5.4        Additional flag to ensure that bond types are written out
               to the ligplot.drw file.
                                                       5 Nov 2008 (RAL)

v.5.0          Addition of -wkdir and -prm command-line arguments
               giving working directory and parameter file. The former
               defines where all files are to be written and read, while
               the latter gives the full path of the parameter file,
               ligplot.prm. Change to get_filename routine so that
               directories not stripped from file name.
                                                       9 Mar 2009 (RAL)
               Write out atom details (as in PDB file) for each atom
               written out to the .drw file.
                                                      31 Jul 2009 (RAL)
               Addition of -wat command-line argument to specify that
               waters are to be included on the plot.
                                                       7 Aug 2009 (RAL)
               Allow "missing H-bonds" to be added to the PDB file using
               keywords HHB and NNB. Format the same as for ligplot.hhb
               records.

               HHB...RRR.C.NNNNI.AAAA.....RRR.C.NNNNI.AAAA....D.DD
               HHB   HIS A  231   N       ASP A  226   OD2    2.77
                                                       8 Aug 2009 (RAL)

               Moved routine for reading in .rcm file prior to the
               flatten_all_objects routine so that initial flattening
               matches the atom position constraints.
                                                      15 Sep 2009 (RAL)

               If all object's atoms have anchors, set its starting
               coords to those anchor positions. 
                                                      21 Sep 2009 (RAL)

               Addition of "blocking" atom positions to .rcm file
               (identified by residue name "...") so that corresponding
               region of plot is kept clear (eg this is the ligand
               location for an earlier ligand in an ensemble from LigPlus)
                                                      10 Mar 2010 (RAL)

               Finally removed lig_search routine and the -h command-line
               parameter so as not to get that silly compiler warning
               about gets being dangerous.
                                                      11 Mar 2010 (RAL)

               Addition of -ligfit option to indicate program to minimize
               ligand to its atom-position restraints. Addition of extra
               minimzation routines specifically for this minimization,
                                                      13 Apr 2010 (RAL)

               Change to prevent right-justificiation of residue name
               if entered with "-n" prefix. Strip out quotes from
               command-line tokens.
                                                       4 May 2010 (RAL)

v.5.0.1        Addition of -conect option to include all CONECT records
               irrespective of their length (used for program schemat in
               PDBsum).
                                                      24 May 2010 (RAL)

               Bug-fix involving use of HHB/NNB records in the PDB file
               for holding additional H-bonds or non-bonded contacts.
                                                       9 Jun 2010 (RAL)

               Allow removal of H-bonds and non-bonded contacts via
               new -HB and -NB records in the PDB file

               -HB...RRR.C.NNNNI.AAAA.....RRR.C.NNNNI.AAAA....D.DD
               -HB   HIS A  231   N       ASP A  226   OD2    2.77
                                                      16 Jun 2010 (RAL)

v.5.0.2        One-line change to fitting algorithm for -ligfit option..
                                                      26 Jul 2010 (RAL)

v.5.0.3        Bug-fix. Addition for v.5.0 making LIGPLOT fail for PDBsum
               because program expecting .hhb and .nnb files in same
               directory as the PDB file, rather than default directory.
                                                      31 Aug 2010 (RAL)

               Bug-fix. Parameters from ligplot.prm file being used with
               default arguments before being read in!
                                                       1 Sep 2010 (RAL)

v.5.0.4        Removal of pare_down_hydrophobics routine so that all
               atoms are written to the .drw file (for use by ligfit)
                                                       5 Oct 2010 (RAL)

               Change so that bonds between atoms in any hydrophobic
               group are not written out to the ligplot.drw file.
                                                       6 Oct 2010 (RAL)

v.5.0.5        Minor adjustment to gaps for double bonds.
                                                       1 Dec 2010 (RAL)

v.5.1          Change to hbadd program so that it can read in a Het
               Group Dictionary in SDF format.
                                                      17 Jan 2011 (RAL)
               Addition of -wcut option to reduce the number of waters
               on the plot by excluding waters with fewer than 2 H-bonds
               to protein/ligand.
                                                      10 Feb 2011 (RAL)
               Extra weighting for anchored ring atoms with the -ligfit
               option so that rings are fitted better for molecules
               with floppy phosphate groups such as ATP.
                                                      11 Mar 2011 (RAL)
v.5.1.1        Bug-fix to atom_linkage function which was storing all
               bonds, rather than just covalent ones. This had an impact
               in the calc_internal_energy routine which needs to skip
               covalently bonded atoms. Hence some atom clashes were not
               scored - leading to bad internal overlaps.
                                                       3 May 2011 (RAL)
v.5.1.2        Bug-fix to length of Print_Title field, which was causing
               program to crash when called by DIMPLOT on macs with a
               very long string in the Title: line of dimplot.pdb.
                                                      19 May 2011 (RAL)
v.5.1.3        Bug-fix to dimer.c to deal with HHB and NNB record in
               input PDB file.
                                                       5 Aug 2011 (RAL)
v.5.1.4        New flag, -clip, which write out the picture limits into
               the PostScript file's BoundingBox parameters.
                                                       3 Oct 2011 (RAL)
v.5.1.5        Additional check for format of .asa file in case it was
               generated by running naccess with the "-f" flag.
                                                      17 Oct 2011 (RAL)
v.5.1.6        Modifications to strcpy(field,field + n) commands which
               fail under some linux operating systems.
                                                      24 Jan 2012 (RAL)
v.5.1.7        Change to PostScript header records in ligps.c to allow
               the -trim option in convert to trim the images correctly.
                                                      17 Feb 2012 (RAL)
v.5.2          Alteration for LigSearch to allow annotations in
               ligplot.drw to show up as shadings on the LIGPLOT diagram.
                                                       5 Mar 2012 (RAL)
v.5.2.1        Bug-fix - LigFit option switched off all ligand flips,
               so it wasn not fitting atoms at all well to the anchor
               positions.Reduced the maximum number of rotatable bonds
               from 20 to 14 for the brute force flip_ligand_bonds
               calculation, otherwise it was taking far too long.
                                                       9 May 2012 (RAL)
               Alteration to flip_ligand_bonds routine to deal with
               large ligands by performing brute force minimization on
               subsets of the bonds to increase chances of finding a
               good fit.
                                                      17 May 2012 (RAL)
v.5.2.2        Bug fix to read of domain file in dimer.c.
                                                      31 May 2012 (RAL)
v.5.3          Change to DIMPLOT operation when residues restraints are
               read in from the dimplot.rcm file.
                                                      26 Jun 2012 (RAL)
               Bug fix. Losing residues from DIMPLOT.
                                                       2 Jul 2012 (RAL)
               Addition of -ctype command-line parameter so that calls
               from LigPlot+ can define the contact type.
                                                       3 Jul 2012 (RAL)
v.5.3.1        Bug-fix to hbadd.c.
                                                       8 Oct 2012 (RAL)
v.5.3.2        Alteration to get_drw_atom to enable it to read in .drw
               files created by LigPlot+.
                                                       5 Mar 2013 (RAL)
v.5.3.3        Bug-fix to hbadd.c for cases where empty Het Group
               dictionary empty causing program to crash.
                                                        RAL 13 May 2013
               Bug-fix to find_hhb_info routine which wasn't checking
               both atom names.
                                                      26 Jun 2013 (RAL)
v.5.3.4        Minor change such that when a .rcm file is used to drive
               a plot, no final rotation is then applied.
                                                       8 Sep 2014 (RAL)
               New routine to close up residues at end in a DIMPLOT
               diagram.
                                                      10 Nov 2014 (RAL)
               Amendments to dimhtml to tighten up the layout of the
               initial dimplot diagram.
                                                      12 Nov 2014 (RAL)
               Addition of -dimfit flag so the tightening routine is
               omitted when fitting to a previous DIMPLOT.
                                                      14 Nov 2014 (RAL)

v.5.3.5        Bug-fix for object types in interface plots
                                                       2 Feb 2015 (RAL)

v.5.3.6        Allow for larger hydrophobic groups in dimhtml
                                                      28 Aug 2015 (RAL)
               Reverse residue order in dimhtml for plotting in portrait
               mode.
                                                       3 Sep 2015 (RAL)
               Make search-size for placing residue names larger for
               interface plots as no space for label found when an H-bond
               group is sandwiched between two non-bonded groups.
                                                      20 Oct 2015 (RAL)
v.5.3.7        Bug-fix in hbadd to identify flourine as a metal.
                                                       4 Jan 2016 (RAL)
               Bug-fix in get_hbplus_info which was skipping residues
               with negative residue numbers.
                                                      22 Jan 2016 (RAL)

v.5.3.8        Addition of plot of red ellipse around selected residue.
                                                      23 Jan 2017 (RAL)
               Added extra variant types.
                                                      14 Feb 2019 (RAL)
v.5.4          Addition of salt-bridges and disulphides to DIMPLOT.
                                                       7 Feb 2020 (RAL)
v.5.4.1        Bug-fix for H-bonds and non-bonded contacts added as HHB
               and NNB records in PDB file not being plotted due to
               changes made in v.5.4.
                                                      11 May 2020 (RAL)
v.5.4.2        New command-line parameter to override abort on bad
               energy of plot, so that squished plots can go through to
               LigPlot+ for manual fixing.
                                                       8 Jul 2020 (RAL)
               Second bug-fix for H-bonds and non-bonded contacts added
               as -HB and -NB records in PDB file not being deleted due to
               changes made in v.5.4.
                                                      15 Jul 2020 (RAL)
v.5.4.3        Change to hbadd.c for new format of the Het Group Dictionary
                                                       8 Nov 2020 (RAL)
v.5.4.4        Change to dimhtml.c to reduce processing time for very
               large interfaces.
               Change to ligplot.c for the same reason.
                                                      25 Jan 2022 (RAL)
               Output various error messages to ligplot.drw file so that
               LigPlus can report the program's failure.
                                                      28 Jan 2022 (RAL)
v.5.4.5        Additions for running in PDBsum1.
                                                      23 Mar 2022 (RAL)
v.5.4.6        Bug fix for special residue interactions not being
               plotted due to changes made in v.5.4.
                                                      13 Sep 2024 (RAL)

----------------------------------------------------------------------*/

/* Compilation
   -----------

Standard compile command:-

   cc -o ligplot ligplot.c -lm

Optimization flags:-  cc -O2 -funroll-loops -o ligplot ligplot.c -lm

To compile with the XMAS-format routines (UCL only):-

   cc -DPLUSXMAS -o ligplot ligplot.c $LIBXMAS $LIBACRM -lm

where $LIBXMAS and $LIBACRM point to Andrew Martin's libxmas.a and
libacrm.a XMAS libraries, respectively.

*/

/* Datafiles
   ---------

Actual name       File pointer    Description
-----------       ------------    -----------
<filename>.pdb    fil_pdb         Input PDB file holding coords of
                                  protein and ligand
<filename>.hhb    fil_hbplus      Input files produced by HBPLUS (or
<filename>.nnb                    otherwise). The .hhb file gives all
                                  the hydrogen bonds in the structure,
                                  while the .nnb file gives all the
                                  non-bonded contacts. These files can
                                  either be in LIGPLOT format or in
                                  HBPLUS format.
<filename>.xmas   fil_xmas        Input file holding XMAS version of
                                  PDB file, H-bonds and non-bonded contacts
                                  as calculated by Andrew Martin's XMAS
                                  programs.
ligplot.bonds     ligplot_bonds   Output file listing of bonds and bond
                  (fil_bonds)     types, for the molecules shown in the
                                  plot.                                 [*]
ligplot.frm       ligplot_frm     Output file in PDB format of the
                                  molecules shown in the plot, prior to
                                  flattening.
ligplot.hhb       ligplot_hhb     Output file containing just the
                                  hydrogen bonds from the original .hhb
                                  file actually used in the plot.       [*]
ligplot.nnb       ligplot_nnb     Output file containing just the
                                  non-bonded contacts from the original
                                  .nnb file actually used in the plot.  [*]
ligplot.prm       ligplot_par     Input parameter file governing
                                  appearance of plot produced.
ligplot.pdb       ligplot_pdb     Output file in PDB format of the
                                  flattened molecules shown in the
                                  plot.                                 [*]
ligplot.ps        ps_file         Output PostScript file of the plot.
ligplot.rcm       ligplot_rcm     Output file listing all the residue
                                  names and numbers in the plot, plus their
                                  centre of mass position on the plot. For
                                  use when generating a plot that is to
                                  have a similar residue layout to a plot
                                  for a related protein/complex. File is
                                  used if LIGPLOT detects <filename>.rcm
                                  file, where <filename> is the name of
                                  the PDB file being plotted.
ligplot.res       ligplot_res     Output file listing all the residue
                                  names and numbers in the plot. For use
                                  by the raslig program in PDBsum LIGPLOT
                                  pages. (File only generated with -r
                                  or -l command-line options).
ligplot.lig       file_ptr        Output file giving coords of ligand
                                  after principal components transformation
                                  (File only generated with -l
                                  command-line option).
ligplot.drw       ligplot_drw     Output file listing all the drawn objects
                                  for import into the Java LIGPLOT Editor, 
                                  LigEd.
ligplot.sum       ligplot_sum     Output file summarising the interactions 
                                  between residue-pairs.

[*] The asterisked output files are taken as input when the program
    is run with ligplot.pdb as the input PDB file (in which case the
    Print_as_is flag is set to TRUE).

*/

/* Function calling-tree
   ---------------------

main
   -> initialise_parameters
         -> have_xmas_link
   -> read_ligplus_params
         -> create_parameter_record
         -> retrieve_parameters
               -> find parameter
   -> read_in_parameters
         -> update_current_parameter
   -> match_colour_names
   -> get_random_number
   -> get_command_arguments
         -> check_for_dimplot
         -> lig_search  (no longer used)
               -> chk_cofactor
               -> first_res
               -> print_cofactor
               -> lig_size
               -> print_chains
               -> ligand_input
         -> get_residue_number
         -> get_filename
               -> getnam
   -> read_in_xmas_file
   -> get_nlig_xmas_molecs
         -> ResetXmasData
         -> ReadXmasData
   -> initialise_variables
   -> define_object_sizes
   -> define_text_sizes
   -> get_xmaslig_start_end
   -> form_plot_title
         -> ResetXmasData
         -> ReadXmasData
         -> to_lower
   -> get_ligand_start_end
         -> get_xmas2pdb_line
               -> ResetXmasData
               -> get_xmas_atom
                     -> ReadXmasData
                     -> remove_dots
         -> check_if_in_ligand
   -> write_drw_headers
         -> to_lower
         -> write_colour_item
               -> get_object_colour
                     -> get_equiv_colour
               -> get_text_colour
                     -> get_equiv_text_colour
   -> read_drw_file
         -> get_line_type
         -> get_drw_residue
               -> create_residue_record
         -> get_drw_atom
               -> create_atom_record
               -> get_element
         -> get_drw_bond
               -> get_tokens
               -> create_bond_record
         -> get_drw_title
         -> create_object_record
   -> read_conec_records
         -> get_xmas2pdb_line
               -> ResetXmasData
               -> get_xmas_atom
                     -> ReadXmasData
                     -> remove_dots
   -> verify_conec_records
         -> get_xmas2pdb_line
               -> ResetXmasData
               -> get_xmas_atom
                     -> ReadXmasData
                     -> remove_dots
   -> read_hbplus_file
         -> get_hbplus_info
               -> check_atom_pair
         -> check_for_metal
         -> check_if_in_ligand
   -> special_resinc
         -> get_hbplus_info
   -> read_pdb_file
         -> get_xmas2pdb_line
               -> ResetXmasData
               -> get_xmas_atom
                     -> ReadXmasData
                     -> remove_dots
         -> get_hbplus_info
         -> find_hhb_info
         -> delete_hhb_info
         -> create_hhb_info_record
         -> check_if_in_ligand
         -> get_element
               -> is_hydrogen_atom
               -> check_for_metal
                     -> is_hydrogen_atom
         -> create_atom_record
         -> create_residue_record
         -> update_bond
               -> create_bond_record
   -> define_objects
         -> create_object_record
   -> get_bonds
         -> create_bond_record
   -> update_all_boundaries
         -> update_residue_boundaries
   -> open_pdb_output
   -> open_bonds_output
   -> covalent_connectivity
         -> update_bond
               -> create_bond_record
   -> interpenetration_check
         -> count_joins
         -> delete_unwanted_objects
               -> delete_object
   -> combine_objects
         -> check_for_join
         -> check_standard_aa
         -> join_residues
               -> get_last
               -> get_pointer_to_residue
               -> get_pointer_to_atom
   -> delete_false_ligands
         -> check_for_metal
   -> add_hbplus_bonds
         -> check_if_in_ligand
         -> update_bond
               -> create_bond_record
   -> get_bond_orders
   -> get_standard_bond_orders
   -> check_metal_bonds
         -> check_for_metal
   -> bond_linkage
   -> check_connections
         -> get_test_bond
         -> remove_unconnected_atoms
               -> initialise_atoms
               -> mark_downstream_atoms
                    -> initialise_bonds
         -> delete_atoms_and_bonds
               -> delete_unwanted_atoms
               -> delete_unwanted_bonds
               -> delete_unwanted_bond_links
               -> delete_unwanted_bond_connections
   -> check_for_cyclics
   -> split_objects
         -> get_last_object
         -> check_for_attachment
   -> assign_bonds_to_objects
         -> delete_object
   -> delete_unwanted_bond_connections
   -> list_bonds
   -> mark_unwanted_waters
   -> mark_reachable_atoms
         -> initialise_atoms
   -> delete_unreachable_residues
         -> delete_object
   -> show_deleted_residues_message
   -> get_atom_sizes
   -> atom_linkage
   -> identify_rotatable_bonds
         -> ring_check
               -> initialise_bonds
   -> write_frm_file
   -> write_lig_prot_bonds
   -> pare_down_hydrophobics
         -> initialise_atoms
   -> delete_atoms_and_bonds
         -> delete_unwanted atoms
         -> delete_unwanted_bonds
         -> delete_unwanted_bond_links
         -> delete_unwanted_bond_connections
   -> calc_fit_coordinates
         -> centre_of_mass
         -> adjust_stored_coords
         -> principal_components
               -> calc_eigen_values
               -> eigen_sort
               -> orien3
                     -> vecprd
                     -> dot3
         -> apply_rotation
   -> flatten_all_objects
         -> create_simple_hgroup
               -> initialise_atoms
         -> delete_atoms_and_bonds
               -> delete_unwanted_atoms
               -> delete_unwanted_bonds
               -> delete_unwanted_bond_links
               -> delete_unwanted_bond_connections
         -> remove_mainchains
               -> initialise_atoms
         -> update_residue_boundaries
         -> place_group_at_origin
         -> flatten_all_rings
               -> get_other_ring_bonds
                     -> initialise_bonds
               -> flatten_ring
                     -> get_ring_atom_coords
                           -> initialise_atoms
                     -> centre_of_mass
                     -> move_object
                     -> adjust_stored_coords
                     -> principal_components
                           -> calc_eigen_values
                           -> eigen_sort
                           -> orien3
                                 -> vecprd
                                 -> dot3
                     -> rotate_object
                     -> squash_z
         -> correct_bond_lengths
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
                     -> apply_transformation
               -> transform_object
                     -> apply_transformation
               -> stretch_bond
                     -> initialise_bonds
                     -> initialise_atoms
         -> flatten_ring_offshoots
               -> offshoot_splay
                     -> plane_transformation
                           -> get_angle
                           -> calculate_matrix
                           -> apply_transformation
                     -> transform_object
                           -> apply_transformation
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> apply_transformation
                     -> transform_downstream_atoms
                           -> initialise_bonds
                           -> initialise_atoms
                           -> apply_transformation
         -> unroll_structure
               -> initialise_atoms
               -> flatten_object
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> transform_object
                     -> flatten_bonds
                           -> plane_transformation
                                 -> get_angle
                                 -> calculate_matrix
                                 -> apply_transformation
                           -> transform_downstream_atoms
                                 -> initialise_bonds
                                 -> initialise_atoms
                                 -> apply_transformation
                           -> line_transformation
                                 -> get_angle
                                 -> calculate_matrix
                           -> apply_transformation
                     -> move_object
                     -> flip_object
               -> move_object
         -> align_object
               -> centre_of_mass
               -> adjust_object_to_origin
               -> adjust_stored_coords
               -> principal_components
                     -> calc_eigen_values
                     -> eigen_sort
                     -> orien3
                           -> vecprd
                           -> dot3
               -> rotate_object
         -> force_flat
         -> create_simple_ligand
               -> stretch_mainchain
                     -> line_transformation
                     -> transform_object
                     -> transform_downlist_atoms
                           -> transform_downstream_atoms
                                 -> initialise_bonds
                                 -> initialise_atoms
                                 -> apply_transformation
                     -> move_object
                     -> apply_transformation
               -> mark_sidechain_atoms
         -> adjust_for_schematic
               -> initialise_atoms
         -> delete_atoms_and_bonds
               -> delete_unwanted atoms
               -> delete_unwanted_bonds
               -> delete_unwanted_bond_links
               -> delete_unwanted_bond_connections
         -> fix_object
               -> update_residue_boundaries
         -> untangle_object
               -> get_random_number
               -> save_restore_coords
               -> calc_internal_energy
                     -> update_residue_boundaries
                     -> fit_object
                           -> matfit
                                 -> qikfit
                           -> get_object_distance_match
                           -> move_object
                           -> flip_object
                           -> rotate_object
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
               -> transform_object
               -> flip_about_bond
                     -> initialise_bonds
               -> flip_back
   -> score_objects
   -> sort_objects
   -> calc_residue_means
   -> read_rcm_file
         -> get_last_object
         -> get_last_residue
         -> get_last_atom
         -> get_filename
   -> lay_interface_objects_out
         -> get_object_means
         -> centre_of_mass
         -> adjust_stored_coords
         -> principal_components
         -> rotate_centres_of_mass
         -> place_interacting_residues
   -> lay_objects_out
         -> fit_object
               -> matfit
                     -> qikfit
               -> get_object_distance_match
               -> move_object
               -> flip_object
               -> rotate_object
         -> update_residue_boundaries
   -> flip_ligand_bonds
         -> save_restore_coords
         -> get_ligand_conformation
               -> bond_flip (see below)
         -> calculate_energy (see below)
   -> energy_minimize
         -> update_all_boundaries (as above)
         -> get_random_number
         -> get_total_energy
               -> calculate_energy
                     -> get_atom_clash_energy
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> get_bond_clash_energy
                           -> apply_transformation
                           -> get_distance_energy
                     -> get_bond_overlap_energy
                           -> apply_transformation
                     -> get_relative_anchor_energy
                           -> compare_anchor_dists
                     -> get_anchor_energy
                     -> get_rel_pstn_energy
                     -> get_boundary_energy
               -> get_closeup_energy
                     -> update_all_boundaries
         -> save_restore_all_objects
               -> save_restore_coords
         -> save_restore_coords
         -> calculate_energy
               -> get_atom_clash_energy
               -> line_transformation
                     -> get_angle
                     -> calculate_matrix
               -> get_bond_clash_energy
                     -> apply_transformation
                     -> get_distance_energy
               -> get_bond_overlap_energy
                     -> apply_transformation
               -> get_anchor_energy
               -> get_rel_pstn_energy
               -> get_boundary_energy
         -> get_closeup_energy
               -> update_all_boundaries
         -> update_residue_boundaries
         -> make_random_move
               -> get_random_number
               -> move_object
         -> make_random_rotation
               -> get_random_number
               -> move_object
               -> rotate_object
         -> make_random_flip
               -> get_random_number
               -> bond_flip
                     -> move_object
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> transform_object
                     -> flip_about_bond
                     -> calc_internal_energy
                           -> update_residue_boundaries
                           -> fit_object
                                 -> matfit
                                       -> qikfit
                                 -> move_object
                                 -> flip_object
                                 -> rotate_object
                      -> flip_attached_groups
                            -> move_object
                            -> transform_object
                            -> flip_object
         -> make_random_atom_flip
               -> get_random_number
               -> move_object
               -> flip_object
         -> get_total_energy
               -> calculate_energy
                     -> get_atom_clash_energy
                     -> line_transformation
                           -> get_angle
                           -> calculate_matrix
                     -> get_bond_clash_energy
                           -> apply_transformation
                           -> get_distance_energy
                     -> get_bond_overlap_energy
                           -> apply_transformation
                     -> get_anchor_energy
                     -> get_rel_pstn_energy
                     -> get_boundary_energy
               -> get_closeup_energy
                     -> update_all_boundaries
   -> ligand_orientate
         -> centre_of_mass
         -> adjust_to_origin
         -> adjust_stored_coords
         -> principal_components
               -> calc_eigen_values
               -> eigen_sort
               -> orien3
                     -> vecprd
                     -> dot3
         -> rotate_whole_structure
         -> swap_x_and_y
         -> check_ligand_direction
               -> flip_about_x
               -> flip_about_y
         -> update_residue_boundaries
   -> tighten_dimplot
         -> update_all_boundaries (as above)
         -> sort_iface_objects
         -> pair_up_objects
         -> pack_paired_objects
         -> update_residue_boundaries
   -> calc_rotation_matrix
   -> rotate_whole_structure
   -> update_all_boundaries (as above)
   -> write_pdb_file
   -> write_bonds
   -> write_ligand_fit_coords
   -> write_rcm_file
   -> write_sum_file
   -> psopen_
         -> define_object_colours
         -> define_text_colours
         -> define_object_sizes
         -> define_text_sizes
         -> pscomm_
         -> pshade_
         -> psubox_
   -> write_drw_file
   -> adjust_for_landscape
   -> psrot_
   -> print_title
         -> pscomm_
         -> pscolb_
         -> psctxt_
   -> print_symbol_key
         -> pscomm_
         -> pssave_
         -> pscolb_
         -> pslwid_
         -> print_key_text
               -> pstext_
         -> psline_
         -> psrest_
         -> psphcl_
         -> pspher_
         -> psdash_
         -> psctxt_
         -> plot_contact_group
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> psarc_
               -> plot_spokes
                     -> psline_
               -> psrest_
         -> accessibility_circle
               -> pscrgb_
               -> psucir_
         -> psccol_
         -> psucir_
   -> extract_minmax
   -> scale
   -> centre_point
   -> define_object_sizes
   -> define_text_sizes
   -> plot_postscript_picture
         -> shade_accessibilities
               -> pscomm_
               -> initialise_atoms
               -> sort_accessibilities
               -> accessibility_cirle
                     -> pscrgb_
                     -> psucir_
               -> psccol_
               -> psucir_
         -> ligsearch_annotations
               -> pscomm_
               -> initialise_atoms
         -> plot_hbond_lines
               -> pscomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> psdash_
               -> psline_
               -> print_hbond_lengths
                     -> pssave_
                     -> pscolb_
                     -> pslwid_
                     -> psubox_
                     -> psctxt_
                     -> psrest_
               -> psrest_
         -> plot_ext_bond_lines
               -> pscomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> psdash_
               -> psline_
               -> psrest_
         -> plot_bond_lines
               -> pscomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> double_bond_split
               -> get_atom_colour
               -> psrest_
         -> plot_aromatics
               -> unset_bond_plot_flags
               -> get_other_ring_bonds
               -> get_subrings
               -> get_ring_atom_coords
                     -> initialise_atoms
               -> centre_of_mass
               -> psomm_
               -> pssave_
               -> pscolb_
               -> pslwid_
               -> pscirc_
               -> psrest_
         -> plot_atoms
               -> pscomm_
               -> pssave_
               -> initialise_atoms
               -> print_current_atom
                     -> pscolb_
                     -> get_atom_colour
                     -> psphcl_
                     -> pspher_
               -> print_atom_name
                     -> position_names
                     -> pssave_
                     -> pscolb_
                     -> pslwid_
                     -> psubox_
                     -> psccol_
                     -> psucir_
                     -> psrest_
                     -> get_atom_colour
                     -> get_atom_text_colour
                     -> psctxt_
               -> psctxt_
               -> psrest_
         -> plot_hydrophobics
               -> pscomm_
               -> pssave_
               -> pslwid_
               -> print_name
                     -> convert_residue_name
                     -> pssave_
                     -> pslwid_
                     -> pscolb_
                     -> pshade_
                     -> psbbox_
                     -> psubox_
                     -> psctxt_
               -> plot_contact_group
                     -> pssave_
                     -> pscolb_
                     -> pslwid_
                     -> psarc_
                     -> plot_spokes
                           -> psline_
                     -> psrest_
               -> psrest_
         -> plot_simple_ligand_residues
               -> pscomm_
               -> pssave_
               -> pslwid_
               -> psccol_
               -> pscolb_
               -> pscirc_
               -> print_name
                     -> convert_residue_name
                     -> pssave_
                     -> pslwid_
                     -> pscolb_
                     -> pshade_
                     -> psbbox_
                     -> psubox_
                     -> psctxt_
               -> psrest_
         -> plot_simple_hgroups
               -> pscomm_
               -> pssave_
               -> pslwid_
               -> print_name
                     -> convert_residue_name
                     -> pssave_
                     -> pslwid_
                     -> pscolb_
                     -> pshade_
                     -> psbbox_
                     -> psubox_
                     -> psctxt_
               -> psrest_
         -> print_residue_names
               -> pssave_
               -> pscolb_
               -> convert_residue_name
                     -> find_ca_position
                     -> coords_resname
               -> resname_grid_search
                     -> initialise_atoms
                     -> extract_search_atoms
                           -> initialise_atoms
                     -> perform_search
                           -> label_dist
               -> psrest_
         -> plot_highlights (in ligps.c)
   -> psclos_
   -> close_output_files
   -> free_up_memory
   -> get_next_xmas_ligand
         -> ResetXmasData
         -> ReadXmasData

*/

#include "ligplot.h"
#include "ligps.c"

/* v.4.2--> */
#ifdef PLUSXMAS
#include "ligxmas.c"
#else
#include "ligxmasx.c"
#endif
/* <--v.4.2 */


/***********************************************************************

initialise_parameters  -  Initialise the plot parameters

***********************************************************************/

void initialise_parameters(void)
{
  int icol, icolour, iline, isize;

  /* Initialise counts */
  n_labels = 0;
  verbose = FALSE;

/* v.5.0--> */
  /* Initialise working directory and name of parameter file */
  wkdir[0] = '\0';
  strcpy(ligplot_prm,"ligplot.prm");
  nensemble = 0;
/* v.5.3.8--> */
  highlight_res[0] = '\0';
  pdb_code[0] = '\0';
  pdbsum_dir[0] = '\0';
/* <--v.5.3.8 */
  ligfit = FALSE;
/* <--v.5.0 */
/* v.5.3.4--> */
  dimfit = FALSE;
/* <--v.5.3.4 */
/* v.5.2 --> */
  ligsearch_id[0] = '\0';
  ligsearch_annot = FALSE;
/* v.5.3 --> */
  contact_type = NOT_DEFINED;
/* <--v.5.3 */
/* <-- v.5.2 */

/* v.5.0.1--> */
  conect = FALSE;
/* <--v.5.0.1 */

  /* Initialise global variables */
  Page_Min_x = BBOXX1;
  Page_Max_x = BBOXX2;
  Page_Min_y = BBOXY1;
  Page_Max_y = BBOXY2;
  Margin_y = (float)Page_Min_y;
  Split_Colour_Ligand_Bonds = FALSE;
  Split_Colour_Nonligand_Bonds = FALSE;
/* v.4.0--> */
  Have_Anchors = FALSE;
/* <--v.4.0 */

/* v.4.2--> */
  /* Check whether XMAS file-reading routines are linked in */
  Have_Xmas = have_xmas_link();

  /* Initialise other flags and parameters */
  Ligand_Molecule_id = -1;
  Output_drw_File_Only = FALSE;
  Read_Xmas = FALSE;
  Plot_All_Xmas_Ligands = FALSE;
/* <--v.4.2 */
  
  /* Initialise pointers */
  first_bond_ptr = NULL;
  first_hhb_info_ptr = NULL;
  last_hhb_info_ptr = NULL;
/* v.5.0--> */
  first_ensemble_ptr = NULL;
/* <--v.5.0 */
  
  /* Initialise picture options */
  Picture = (struct Picture_Options *)
    malloc(sizeof(struct Picture_Options));
  Picture->In_Colour = FALSE;
  Picture->Portrait = TRUE;
/* v.5.3.4--> */
  Picture->Allow_Rotation = TRUE;
/* <--v.5.3.4 */

  /* Title of plot */
  root_name[0] = '\0';
  Print_Title[0] = '\0';
/* v.5.4.5 --> */
  Drw_Title[0] = '\0';
/* <-- v.5.4.5 */
  
  /* Special flags and values */
  xsite_file = FALSE;
  Max_atom_radius = 0.0;
  Maximum_accessibility = (float) 100.0;
  Nobjects = 0;
  Plot_Centre_x = (float) ((Page_Min_x + Page_Max_x) / 2.0);
  Plot_Centre_y = (float) ((Page_Min_y + Page_Max_y) / 2.0);
  
  /* Parameters determining what's to be included in the plot */
  Include = (struct Include_Parameters *)
    malloc(sizeof(struct Include_Parameters));
  Include->Waters = FALSE;
  Include->Mainchain_Atoms = TRUE;
  Include->Hydrophobics = TRUE;
  Include->Hbonds = TRUE;
  Include->External_Bonds = TRUE;
  Include->Hydrophobic_Bonds = FALSE;
  Include->Internal_Hbonds = FALSE;
  Include->Simple_Ligand_Residues = FALSE;
  Include->Simple_Nonligand_Residues = FALSE;
  Include->Accessibilities = FALSE;
  Include->Ligand_Atoms = TRUE;
  Include->Nonligand_Atoms = TRUE;
/* v.4.0--> */
  Include->Water_Atoms = TRUE;
  Include->Ligand_Accessibilities_Only = TRUE;
/* <--v.4.0 */
/* v.4.3--> */
  Include->Chemical_Notation = FALSE;
/* <--v.4.3 */
  Include->Atom_Names = TRUE;
  Include->Residue_Names = TRUE;
  Include->Hbond_Lengths = TRUE;
  Include->Linked_Residues = FALSE;
  Include->Key = TRUE;
  Include->Double_Bonds = FALSE;
  Include->Filename_for_Title = TRUE;
  Include->External_Bonds_Solid = FALSE;
  Include->Contact_Type = 0;
  Include_Hydrogens = FALSE;
/* v.4.4.1--> */
  Greek_Names = FALSE;
/* <--v.4.4.1 */
  
  /* Residues H-bonded to non-ligand residues */
  for (iline = 0; iline < MAXSPECIAL_RES; iline++)
    strcpy(special_res[iline],"       ");
  
  /* Sizes of objects on plot */
  for (isize = 0; isize < OBJECT_SIZES; isize++)
/* v.4.2--> */
/*    object_size[isize] = 0.0; */
    param_object_size[isize] = 0.0;
/* <--v.4.2 */
  Size = (struct Sizes *) malloc(sizeof(struct Sizes));
  Size_Val = (struct Size_Vals *) malloc(sizeof(struct Size_Vals));
  strcpy(Size->Ligand_Atoms,           "Ligatom_radius ");
  strcpy(Size->Nonligand_Atoms,        "Nligatom_radius");
  strcpy(Size->Waters,                 "Water_radius   ");
  strcpy(Size->Hydrophobics,           "Hphobic_radius ");
  strcpy(Size->Simple_Residues,        "Simple_radius  ");
  strcpy(Size->Ligand_Bonds,           "Ligbond_width  ");
  strcpy(Size->Nonligand_Bonds,        "Nligbond_width ");
  strcpy(Size->Hydrogen_Bonds,         "Hbond_width    ");
  strcpy(Size->External_Bonds,         "Ext_bond_width ");

  /* Text sizes on plot */
  for (isize = 0; isize < TEXT_SIZES; isize++)
/* v.4.2--> */
/*    text_size[isize] = 0.0; */
    param_text_size[isize] = 0.0;
/* <--v.4.2 */
  Text_Size = (struct Text_Sizes *) malloc(sizeof(struct Text_Sizes));
  Text_Size_Val 
    = (struct Text_Size_Vals *) malloc(sizeof(struct Text_Size_Vals));
  strcpy(Text_Size->Ligand_Residue_Names,   "Ligresname_size");
  strcpy(Text_Size->Nonligand_Residue_Names,"Nligresnam_size");
  strcpy(Text_Size->Water_Names,            "Watername_size ");
  strcpy(Text_Size->Hydrophobic_Names,      "Hydrophnam_size");
  strcpy(Text_Size->Simple_Residue_Names,   "Simpletext_size");
  strcpy(Text_Size->Ligand_Atom_Names,      "Ligatmname_size");
  strcpy(Text_Size->Nonligand_Atom_Names,   "Nligatmnam_size");
  strcpy(Text_Size->Hbond_Lengths,          "HBtext_size    ");

/* v.4.1--> */
  /* Initialise the grey values */
  for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
    object_grey[icolour] = 0;
  object_grey[0] = 10;
  object_grey[6] = 10;
  object_grey[8] = 8;
  object_grey[9] = 4;
  object_grey[11] = 9;
  object_grey[12] = 2;
  object_grey[13] = 4;
  object_grey[14] = 4;
  object_grey[15] = 4;

  /* Initialise the inverse colours */
  for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
    object_inverse[icolour] = 0;
  object_inverse[1] = 5;
  object_inverse[2] = 6;
  object_inverse[3] = 9;
  object_inverse[4] = 8;
  object_inverse[5] = 7;
  object_inverse[8] = 11;
  object_inverse[9] = 2;
  object_inverse[10] = 12;
  object_inverse[11] = 5;
  object_inverse[12] = 8;
  object_inverse[13] = 10;
  object_inverse[14] = 8;
  object_inverse[15] = 10;
  object_inverse[16] = 0;
  object_inverse[17] = 7;
/* <--v.4.1 */

  /* Object colours */
  for (icolour = 0; icolour < OBJECT_COLOURS; icolour++)
    strcpy(object_colour[icolour],"Grey00");
/* v.4.1-->
  strcpy(object_colour[0],"Grey10");
  strcpy(object_colour[6],"Grey10");
  strcpy(object_colour[8],"Grey08");
  strcpy(object_colour[9],"Grey04");
  strcpy(object_colour[11],"Grey09");
  strcpy(object_colour[12],"Grey02");
  strcpy(object_colour[13],"Grey04");
  strcpy(object_colour[14],"Grey04");
  strcpy(object_colour[15],"Grey04");
<--v.4.1 */
  Colour = (struct Colours *) malloc(sizeof(struct Colours));
  strcpy(Colour->Background,       "Background_col    ");
  strcpy(Colour->Ligand_Bonds,     "Ligbond_colour    ");
  strcpy(Colour->Nonligand_Bonds,  "Nligbond_colour   ");
  strcpy(Colour->Hydrogen_Bonds,   "Hbond_colour      ");
  strcpy(Colour->External_Bonds,   "Ext_bond_colour   ");
  strcpy(Colour->Hydrophobics,     "Hydrophobic_col   ");
  strcpy(Colour->Accessibility_Min,"Accessib_min      ");
  strcpy(Colour->Accessibility_Max,"Accessib_max      ");
  strcpy(Colour->Nitrogen,         "Nitrogen_colour   ");
  strcpy(Colour->Oxygen,           "Oxygen_colour     ");
  strcpy(Colour->Carbon,           "Carbon_colour     ");
  strcpy(Colour->Sulphur,          "Sulphur_colour    ");
  strcpy(Colour->Water,            "Water_colour      ");
  strcpy(Colour->Phosphorus,       "Phosphorus_colour ");
  strcpy(Colour->Iron,             "Iron_colour       ");
  strcpy(Colour->Other,            "Other_colour      ");
  strcpy(Colour->Atom_Edges,       "Atom_edge_col     ");
  strcpy(Colour->Simple_Residues,  "Simpleres_col     ");
  
  /* Accessibility maximum and minimum */
  for (icol = 0; icol < 3; icol++)
    {
      Accessibility_Max[icol] = (float) (1.0);
      Accessibility_Min[icol] = (float) (0.2);
    }

/* v.4.1--> */
  /* Initialise the grey values */
  for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
    text_grey[icolour] = 0;

  /* Initialise the inverse colors */
  for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
    text_inverse[icolour] = 1;
  text_inverse[0] = 5;
  text_inverse[1] = 5;
  text_inverse[2] = 5;
  text_inverse[3] = 6;
  text_inverse[4] = 8;
  text_inverse[5] = 1;
  text_inverse[6] = 1;
  text_inverse[7] = 9;
  text_inverse[8] = 9;
/* <--v.4.1 */
  
  /* Text colours */
  for (icolour = 0; icolour < TEXT_COLOURS; icolour++)
    strcpy(text_colour[icolour],"Grey00");
  Text_Colour = (struct Text_Colours *) malloc(sizeof(struct Text_Colours));
  strcpy(Text_Colour->Title,                  "Title_colour      ");
  strcpy(Text_Colour->Key_Text,               "Key_Text_colour   ");
  strcpy(Text_Colour->Ligand_Residue_Names,   "Ligresname_col    ");
  strcpy(Text_Colour->Nonligand_Residue_Names,"Nligresname_col   ");
  strcpy(Text_Colour->Water_Names,            "Watername_col     ");
  strcpy(Text_Colour->Hydrophobic_Names,      "Hydrophname_col   ");
  strcpy(Text_Colour->Ligand_Atom_Names,      "Ligatmname_col    ");
  strcpy(Text_Colour->Nonligand_Atom_Names,   "Nligatmname_col   ");
  strcpy(Text_Colour->Hbond_Lengths,          "HBlength_col      ");
  
  /* Colour definitions */
  for (icolour = 0; icolour < MAX_COLOURS; icolour++)
    {
      strcpy(Colour_Table_Name[icolour],"BLACK          ");
      Colour_Table[icolour][0] = 0.0000;
      Colour_Table[icolour][1] = 0.0000;
      Colour_Table[icolour][2] = 0.0000;
    }
  
  /* Minimization parameters */
  Atom_Atom_Clash = (float) (10.0);
  Bond_Atom_Clash = (float) (0.20);
  Overlap_Score = (float) (50.0);
  HB_Weight = (float) (0.5);
  Nonbond_Weight = (float) (0.2);
  Internal_Energy_Weight = (float) (10.0);
  Max_HBdist = (float) (15.0);
  Bond_Stretch = (float) (1.0);
  Energy_Drop_Maximum = 0.0;
  Max_Loops = 1000;
  Random_Start = FALSE;
/* v.4.0--> */
  Anchor_Weight = (float) (2.0);
  Boundary_Energy_Weight = (float) (1.0);
  Min_Boundary_Dist = (float) (2.0);
  Rel_Pstn_Energy_Weight = (float) (-1.0);
/* <--v.4.0 */
/* v.5.3.4--> */
  Closeup_Energy_Weight = (float) (1.0);
/* <--v.5.3.4 */
}
/***********************************************************************

update_current_parameter  -  Update the current plot parameter with
                             the value just read in

***********************************************************************/

void update_current_parameter(int parameter_group, int parameter,
                              int true_false, float size,
                              char colour[COL_NAME_LEN + 1],
                              float red, float green, float blue,
                              char residue_pair[8],
                              char special_res[MAXSPECIAL_RES][15])
{
  /* Process according to the parameter-group */
  switch (parameter_group)
    {
      /* Colour option */
    case 0 :
      
      /* Update the appropriate parameter */
      switch (parameter)
        {
        case 0 : Picture->In_Colour = true_false; break;
        case 1 : Picture->Portrait = true_false; break;
        case 2 : Picture->Rotation_Angle = size; break;
        default :
          printf("*** Program error. Invalid Picture parameter: %d\n",
                 parameter);
          exit(1);
        }
      break;

      /* Include options */
    case 1 :

      /* Update the appropriate parameter */
      switch (parameter)
        {
        case  0 : Include->Hydrophobics = true_false; break;
/* v.5.0--> */
        /* case  1 : Include->Waters = true_false; break; */
        case  1 :
          Include->Waters = true_false;
          if (Include_Waters == TRUE)
            Include->Waters = TRUE;
          break;
/* <--v.5.0 */
        case  2 : Include->Mainchain_Atoms = true_false; break;
        case  3 : Include->Linked_Residues = true_false; break;
        case  4 : Include->Hbonds = true_false; break;
        case  5 : Include->Internal_Hbonds = true_false; break;
        case  6 : Include->External_Bonds = true_false; break;
        case  7 : Include->Hydrophobic_Bonds = true_false; break;
        case  8 : Include->Simple_Ligand_Residues = true_false; break;
        case  9 : Include->Simple_Nonligand_Residues = true_false; break;
        case 10 : Include->Accessibilities = true_false; break;
        case 11 : Include->Ligand_Atoms = true_false; break;
        case 12 : Include->Nonligand_Atoms = true_false; break;
        case 13 : Include->Double_Bonds = true_false; break;
        case 14 : Include->Key = true_false; break;
        case 15 : Include->Residue_Names = true_false; break;
        case 16 : Include->Atom_Names = true_false; break;
        case 17 : Include->Hbond_Lengths = true_false; break;
        case 18 : Include->Filename_for_Title = true_false; break;
        case 19 : Include->External_Bonds_Solid = true_false; break;
        case 20 : Include->Contact_Type = (int) size; break;
/* v.4.0--> */
        case 21 : Include->Water_Atoms = true_false; break;
        case 22 : Include->Ligand_Accessibilities_Only = true_false; break;
/* <--v.4.0 */
/* v.4.3--> */
        case 23 : Include->Chemical_Notation = true_false; break;
/* <--v.4.3 */
        default :
          printf("*** Program error. Invalid Include parameter: %d\n",
                 parameter);
          exit(1);
        }
      break;

      /* Linked residues */
    case 2 :
      
      /* Update the appropriate parameter */
      strncpy(special_res[parameter],residue_pair,8);
      special_res[parameter][8] = '\0';
      break;

      /* Sizes */
    case 3 :
      
      /* Update the appropriate parameter */
/* v.4.2--> */
/*      object_size[parameter] = size; */
      param_object_size[parameter] = size;
/* <--v.4.2 */
      break;

      /* Text sizes */
    case 4 :
      
      /* Update the appropriate parameter */
/* v.4.2--> */
/*      text_size[parameter] = size; */
      param_text_size[parameter] = size;
/* <--v.4.2 */
      break;

      /* Object colours */
    case 5 :

      /* Update the appropriate object colour */
/* v.4.1-->
      if (Picture->In_Colour == TRUE)
<--v.4.1 */
      strcpy(object_colour[parameter],colour);

      /* If black-and-white, update only if this is an ATOM option
         on ligand or non-ligand bonds */
/* v.4.1-->
      else if ((parameter == 1 || parameter == 2) &&
               !strncmp(colour,"ATOM           ",15))
        strcpy(object_colour[parameter],colour);
<--v.4.1 */
      break;

      /* Text colours */
    case 6 :
      
      /* Update the appropriate text colour */
/* v.4.1-->
      if (Picture->In_Colour == TRUE)
<--v.4.1 */
      strcpy(text_colour[parameter],colour);
      break;

      /* Colour definitions */
    case 7 :

      /* Update the appropriate parameter */
      Colour_Table[parameter][0] = red;
      Colour_Table[parameter][1] = green;
      Colour_Table[parameter][2] = blue;
      strcpy(Colour_Table_Name[parameter],colour);
      break;

      /* Minimization parameters */
    case 8 :

      /* Update the appropriate parameter */
      switch (parameter)
        {
        case  0 : Atom_Atom_Clash = size; break;
        case  1 : Bond_Atom_Clash = size; break;
        case  2 : Overlap_Score = size; break;
        case  3 : HB_Weight = size; break;
        case  4 : Nonbond_Weight = size; break;
        case  5 : Internal_Energy_Weight = size; break;
        case  6 : Max_HBdist = size; break;
        case  7 : Bond_Stretch = size; break;
        case  8 : Max_Loops = (int) size; break;
        case  9 : Energy_Drop_Maximum = size; break;
        case 10 : Random_Start = true_false; break;
/* v.4.0--> */
        case 11 : Anchor_Weight = size; break;
        case 12 : Boundary_Energy_Weight = size; break;
        case 13 : Min_Boundary_Dist = size; break;
        case 14 : Rel_Pstn_Energy_Weight = size; break;
/* <--v.4.0 */
        default :
          printf("*** Program error. \n");
          printf("*** Invalid Minimization parameter: %d\n",
                 parameter);
          exit(1);
        }
      break;

      /* Error message */
    default :
      printf("*** Program error. Invalid parameter group: %d\n",
             parameter_group);
      exit(1);
    }
}
/* v.5.0--> */
/***********************************************************************

store_string  -  Allocate memory for given string and store pointer to
                 it

***********************************************************************/

void store_string(char **string_ptr,char *string)
{
  int len;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for the string */
  *string_ptr = (char *) malloc(sizeof(char)*(len + 1));

  /* Store the string */
  strcpy(*string_ptr,string);
}
/***********************************************************************

create_parameter_record  -  Create a new parameter record

***********************************************************************/

struct parameter *create_parameter_record(struct parameter **fst_parameter_ptr,
                                          struct parameter **lst_parameter_ptr,
                                          char *input_line,char *string_ptr)
{
  struct parameter *parameter_ptr;
  struct parameter *first_parameter_ptr, *last_parameter_ptr;

  /* Initialise parameter pointers */
  parameter_ptr = NULL;
  last_parameter_ptr = *lst_parameter_ptr;
  first_parameter_ptr = *fst_parameter_ptr;

  /* Create parameter entry */
  parameter_ptr = (struct parameter*)malloc(sizeof(struct parameter));
  if (parameter_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct parameter\n");
      printf("***        Program terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  store_string(&(parameter_ptr->value),string_ptr+1);
  string_ptr[0] = '\0';
  store_string(&(parameter_ptr->key),input_line);

  /* Initialise pointer to next parameter record */
  parameter_ptr->next_parameter_ptr = NULL;

  /* If this is the first entry stored, then update pointer to head of
     linked list */
  if (first_parameter_ptr == NULL)
    first_parameter_ptr = parameter_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_parameter_ptr->next_parameter_ptr = parameter_ptr;
         
  /* Save the current parameter's pointer */
  last_parameter_ptr = parameter_ptr;

  /* Return current pointers */
  *fst_parameter_ptr = first_parameter_ptr;
  *lst_parameter_ptr = last_parameter_ptr;
  return(parameter_ptr);
}
/***********************************************************************

find_parameter  -  Find the given parameter entry

***********************************************************************/

struct parameter *find_parameter(struct parameter *first_parameter_ptr,
                                 char *key,int *true_false,float *size,
                                 char *colour_name,float *red,
                                 float *green,float *blue)
{
  struct parameter *parameter_ptr, *found_parameter_ptr;

  /* Initialise variables */
  found_parameter_ptr = NULL;

  /* Get pointer to the first parameter */
  parameter_ptr = first_parameter_ptr;

  /* Loop over all the parameters to find the one of interest */
  while (parameter_ptr != NULL && found_parameter_ptr == NULL)
    {
      /* If this is the right parameter, store pointer */
      if (!strcmp(key,parameter_ptr->key))
        found_parameter_ptr = parameter_ptr;

      /* Get pointer to the next parameter */
      parameter_ptr = parameter_ptr->next_parameter_ptr;
    }

  /* Return the parametering parameter pointer */
  return(found_parameter_ptr);
}
/***********************************************************************

retrieve_parameters  -  Store the LigPlus parameters just read in

***********************************************************************/

void retrieve_parameters(struct parameter *first_parameter_ptr)
{
  char colour_name[COL_NAME_LEN + 1], residue_pair[8];
  char special_res[MAXSPECIAL_RES][15];

  int ikey, parameter, parameter_group;
  int true_false;

  float blue, green, red, size;

  struct parameter *parameter_ptr;

  /* Loop over the LigPlus parameters */
  for (ikey = 0; ikey < NPARAM_KEYS; ikey++)
    {
      /* Use the key to locate this parameter */
      /*
      parameter_ptr
        = find_parameter(first_parameter_ptr,PARAM_KEY[ikey].key,
                         &true_false,&size,&colour_name,&red,&green,
                         &blue);
      */
      /* If found, then process */
      if (parameter_ptr != NULL)
        {

          /* Update the current parameter */
          update_current_parameter(parameter_group,parameter,true_false,
                                   size,colour_name,red,green,blue,
                                   residue_pair,special_res);
        }
    }

}
/***********************************************************************

read_ligplus_params  -  Read in the plot parameters output by LigPlot+

***********************************************************************/

void read_ligplus_params(void)
{
  char line[LINELEN + 1];

  char *string_ptr;

  int nparams;

  struct parameter *parameter_ptr;
  struct parameter *first_parameter_ptr, *last_parameter_ptr;

  /* Initialise variables */
  nparams = 0;

  /* Initialise pointers */
  parameter_ptr = first_parameter_ptr = last_parameter_ptr = NULL;

  /* Open the parameter file */
  if ((ligplot_par = fopen(ligplot_prm,"r")) == NULL)
    {
      printf("\n*** Unable to open parameter file: %s\n\n",ligplot_prm);
      exit(1);
    }

  /* Loop while reading records from the parameter file */
  while (fgets(line,LINELEN,ligplot_par)!= NULL)
    {
      /* See if line contains an equals sign */
      string_ptr = strchr(line,'=');

      /* If it does, then store the keyword/value pair */
      if (string_ptr != NULL)
        {
          /* Create a parameter record */
          parameter_ptr
            = create_parameter_record(&first_parameter_ptr,
                                      &last_parameter_ptr,line,string_ptr);

          /* Increment count of parameters read in */
          nparams++;
        }
    }

  /* Close the parameter file */
  fclose(ligplot_par);

  /* Retrieve the parameters */
  retrieve_parameters(first_parameter_ptr);

}
/* <--v.5.0 */
/* v.4.1--> */
/***********************************************************************

get_equiv_colour  -  Get the equivalent colour position for the given
                     colour string

***********************************************************************/

int get_equiv_colour(int icol)
{
  int icolour, jcol;
  int match;

  /* Initialise variables */
  icolour = 0;

  /* Locate this colour in the table of colour names */
  match = FALSE;
  for (jcol = 0; jcol < MAX_COLOURS && match == FALSE; jcol++)
    {
      if (!strncmp(object_colour[icol],Colour_Table_Name[jcol],
                   COL_NAME_LEN))
        {
          match = TRUE;
          icolour = jcol;
        }
    }

  /* If no match found, check whether colour defined as "ATOM" */
  if (match == FALSE)
    icolour = ncols;

  /* Return the colour number */
  return(icolour);
}
/***********************************************************************

get_equiv_text_colour  -  Get the equivalent colour position for the
                          given colour string

***********************************************************************/

int get_equiv_text_colour(int icol)
{
  int icolour, jcol;
  int match;

  /* Initialise variables */
  icolour = 0;

  /* Locate this colour in the table of colour names */
  match = FALSE;
  for (jcol = 0; jcol < MAX_COLOURS && match == FALSE; jcol++)
    {
      if (!strncmp(text_colour[icol],
                   Colour_Table_Name[jcol],COL_NAME_LEN))
        {
          match = TRUE;
          icolour = jcol;
        }
    }

  /* Return the colour number */
  return(icolour);
}
/***********************************************************************

get_object_colour  -  Get the equivalent colour position for the given
                      object's colour string

***********************************************************************/

int get_object_colour(char colour[COL_NAME_LEN + 1],int *igrey,int *inverse)
{
  int icol, icolour;

  /* Initialise variables */
  icol = 0;

  /* Check for each string in turn */
  if (!strcmp(colour,Colour->Background))
    icol = 0;
  else if (!strcmp(colour,Colour->Ligand_Bonds))
    icol = 1;
  else if (!strcmp(colour,Colour->Nonligand_Bonds))
    icol = 2;
  else if (!strcmp(colour,Colour->Hydrogen_Bonds))
    icol = 3;
  else if (!strcmp(colour,Colour->External_Bonds))
    icol = 4;
  else if (!strcmp(colour,Colour->Hydrophobics))
    icol = 5;
  else if (!strcmp(colour,Colour->Accessibility_Min))
    icol = 6;
  else if (!strcmp(colour,Colour->Accessibility_Max))
    icol = 7;
  else if (!strcmp(colour,Colour->Nitrogen))
    icol = 8;
  else if (!strcmp(colour,Colour->Oxygen))
    icol = 9;
  else if (!strcmp(colour,Colour->Carbon))
    icol = 10;
  else if (!strcmp(colour,Colour->Sulphur))
    icol = 11;
  else if (!strcmp(colour,Colour->Water))
    icol = 12;
  else if (!strcmp(colour,Colour->Phosphorus))
    icol = 13;
  else if (!strcmp(colour,Colour->Iron))
    icol = 14;
  else if (!strcmp(colour,Colour->Other))
    icol = 15;
  else if (!strcmp(colour,Colour->Atom_Edges))
    icol = 16;
  else if (!strcmp(colour,Colour->Simple_Residues))
    icol = 17;

  /* Locate this colour in the table of colour names */
  icolour = get_equiv_colour(icol);

  /* Get its grey equivalent */
  *igrey = object_grey[icol];

  /* Get its inverse colour equivalent */
  *inverse = object_inverse[icol];

  /* Return the colour number */
  return(icolour);
}
/***********************************************************************

get_text_colour  -  Get the equivalent colour position for the given
                    text item's colour string

***********************************************************************/

int get_text_colour(char colour[COL_NAME_LEN + 1],int *igrey,int *inverse)
{
  int icol, icolour;

  /* Initialise variables */
  icol = 0;

  /* Check for each string in turn */
  if (!strcmp(colour,Text_Colour->Title))
    icol = 0;
  else if (!strcmp(colour,Text_Colour->Key_Text))
    icol = 1;
  else if (!strcmp(colour,Text_Colour->Ligand_Residue_Names))
    icol = 2;
  else if (!strcmp(colour,Text_Colour->Nonligand_Residue_Names))
    icol = 3;
  else if (!strcmp(colour,Text_Colour->Water_Names))
    icol = 4;
  else if (!strcmp(colour,Text_Colour->Hydrophobic_Names))
    icol = 5;
  else if (!strcmp(colour,Text_Colour->Ligand_Atom_Names))
    icol = 6;
  else if (!strcmp(colour,Text_Colour->Nonligand_Atom_Names))
    icol = 7;
  else if (!strcmp(colour,Text_Colour->Hbond_Lengths))
    icol = 8;

  /* Locate this colour in the table of colour names */
  icolour = get_equiv_text_colour(icol);

  /* Get its grey equivalent */
  *igrey = text_grey[icol];

  /* Get its inverse colour equivalent */
  *inverse = text_inverse[icol];

  /* Return the colour number */
  return(icolour);
}
/***********************************************************************

write_colour_item  -  Write out all the colour-number associated with
                      the colour of this object to the ligplot.drw file

***********************************************************************/

void write_colour_item(char item_name[COL_NAME_LEN + 1],int object_type)
{
  int icol, igrey, inverse;

  /* Get the associated colour number, and write out */
  if (object_type == 0)
    icol = get_object_colour(item_name,&igrey,&inverse);
  else
    icol = get_text_colour(item_name,&igrey,&inverse);
  fprintf(ligplot_drw,"%2d %2d %2d %s\n",icol,inverse,igrey,item_name);
}
/***********************************************************************

write_drw_headers  -  Write out the plot parameters to the ligplot.drw
                      file

***********************************************************************/

/* v.4.2--> */
/* void write_drw_headers(void) */
void write_drw_headers(char pdb_code[5])
/* <--v.4.2 */
{
  int isize;
  int i, icol, ipos, len, slen;
  int done;
/* v.4.2--> */
  int last_letter;
/* <--v.4.2 */

  char colour_name[COL_NAME_LEN + 1];
  char number_string[10];
/* v.4.2--> */
  char filename[FILENAME_LEN];
/* <--v.4.2 */

  printf("\nWriting out the drw file...\n");

/* v.5.0--> */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
/* <--v.5.0 */

/* v.4.2--> */
  /* If plotting all ligands in XMAS file, then add molcule identifier
     to name of file */
  if (Plot_All_Xmas_Ligands == TRUE)
    {
      sprintf(number_string,"%4d",Ligand_Molecule_id);

      /* Zero-fill leading spaces */
      for (ipos = 0; ipos < 4; ipos++)
        if (number_string[ipos] == ' ')
          number_string[ipos] = '0';

      /* Form the filename of the ligplot.drw file */
/* v.5.0--> */
      /* strcpy(filename,"lig"); */
      strcat(filename,"lig");
/* <--v.5.0 */
      strcat(filename,pdb_code);
      strcat(filename,number_string);
      strcat(filename,".drw");
    }
  else
/* v.5.0--> */
    /* strcpy(filename,"ligplot.drw"); */
    strcat(filename,"ligplot.drw");
/* <--v.5.0 */
/* <--v.4.2 */

  /* Open drawn objects file, ligplot.drw */
/* v.4.2--> */
/*  if ((ligplot_drw = fopen("ligplot.drw","w")) ==  NULL) */
  if ((ligplot_drw = fopen(filename,"w")) ==  NULL)
/* <--v.4.2 */
    {
      printf("\n*** Unable to open drawn objects file,");
/* v.4.2--> */
/*      printf(" ligplot.drw\n"); */
      printf(" %s\n",filename);
/* <--v.4.2 */
      exit(1);
    }

  /* Write first record signifying start of plot parameters */
  fprintf(ligplot_drw,"#P\n");

  /* Write out the plot parameters */

  /* Colour/black-and-white option */
  if (Picture->In_Colour == TRUE)
    fprintf(ligplot_drw,"C\n");
  else
    fprintf(ligplot_drw,"B\n");

  /* Portrait/Landscape option */
  if (Picture->Portrait == TRUE)
    fprintf(ligplot_drw,"P\n");
  else
    fprintf(ligplot_drw,"L\n");

  /* Ligand atoms on/off */
/* v.4.4--> */
/*  if (Include->Ligand_Atoms == TRUE) */
  if (Include->Ligand_Atoms == TRUE && Include->Chemical_Notation == FALSE)
/* <--v.4.4 */
    fprintf(ligplot_drw,"Y\n");
  else
    fprintf(ligplot_drw,"N\n");

  /* Non-ligand atoms on/off */
/* v.4.4--> */
/*  if (Include->Nonligand_Atoms == TRUE) */
  if (Include->Nonligand_Atoms == TRUE && Include->Chemical_Notation == FALSE)
/* <--v.4.4 */
    fprintf(ligplot_drw,"Y\n");
  else
    fprintf(ligplot_drw,"N\n");

  /* Water atoms on/off */
/* v.4.4--> */
/*  if (Include->Water_Atoms == TRUE) */
  if (Include->Water_Atoms == TRUE && Include->Chemical_Notation == FALSE)
/* <--v.4.4 */
    fprintf(ligplot_drw,"Y\n");
  else
    fprintf(ligplot_drw,"N\n");

  /* Write out #M record as start of colour-map definitions */
  fprintf(ligplot_drw,"#M\n");
  ncols = 0;

  /* Loop to write out colour definitions */
  for (icol = 0; icol < MAX_COLOURS; icol++)
  {
    if (Colour_Table_Name[icol][0] != '\0')
      {
        /* Write out the colour number */
        fprintf(ligplot_drw,"%2d ",icol);

        /* Process each of the RGB numbers to compress */
        for (i = 0; i < 3; i++)
          {
            /* Get this RGB value */
            sprintf(number_string,"%6.4f",Colour_Table[icol][i]);
            number_string[6] = '\0';

            /* Remove any trailing zeroes */
            done = FALSE;
            for (ipos = 5; ipos > 2 && done == FALSE; ipos--)
              {
                /* If this is a zero, then reduce length of string */
                if (number_string[ipos] == '0')
                  number_string[ipos] = '\0';
                else
                  done = TRUE;
              }

            /* Write out the value */
            fprintf(ligplot_drw,"%s ",number_string);
          }

        /* Convert all but first character of colour name to lower case */
        strcpy(colour_name,Colour_Table_Name[icol]);

        /* Get the length of the colour name */
        len = 0;
        slen = strlen(colour_name);
        for (i = 0; i < slen; i++)
          if (colour_name[i] != ' ' && colour_name[i] != '\0')
            len = i;
        to_lower(colour_name+1,len);
/* v.4.2--> */
        last_letter = 0;
/* <--v.4.2 */

        /* Replace any spaces by underscores */
        for (i = 0; i < len; i++)
          {
            if (colour_name[i] == ' ')
              colour_name[i] = '_';
/* v.4.2--> */
            else if (colour_name[i] != '_')
              last_letter = i;
/* <--v.4.2 */
          }

/* v.4.2--> */
        /* Terminate the colour name */
        colour_name[last_letter + 1] = '\0';
/* <--v.4.2 */

        /* Write out converted colour name */
        fprintf(ligplot_drw,"%s\n",colour_name);
        ncols++;
      }
  }

  /* Write out the ATOM colour for bonds */
  fprintf(ligplot_drw,"%2d 0.0 0.0 0.0 ATOM_colour\n",ncols);

  /* Write out the various object colours */
  fprintf(ligplot_drw,"#C\n");

  /* Get each object's position in the colour table */
  write_colour_item(Colour->Background,0);
  write_colour_item(Colour->Ligand_Bonds,0);
  write_colour_item(Colour->Nonligand_Bonds,0);
  write_colour_item(Colour->Hydrogen_Bonds,0);
  write_colour_item(Colour->External_Bonds,0);
  write_colour_item(Colour->Hydrophobics,0);
  write_colour_item(Colour->Atom_Edges,0);
  write_colour_item(Colour->Nitrogen,0);
  write_colour_item(Colour->Oxygen,0);
  write_colour_item(Colour->Carbon,0);
  write_colour_item(Colour->Sulphur,0);
  write_colour_item(Colour->Water,0);
  write_colour_item(Colour->Phosphorus,0);
  write_colour_item(Colour->Iron,0);
  write_colour_item(Colour->Other,0);
  write_colour_item(Colour->Simple_Residues,0);

  /* Get each text item's position in the colour table */
  write_colour_item(Text_Colour->Title,1);
/*  write_colour_item(Text_Colour->Key_Text,1); */
  write_colour_item(Text_Colour->Ligand_Residue_Names,1);
  write_colour_item(Text_Colour->Nonligand_Residue_Names,1);
  write_colour_item(Text_Colour->Water_Names,1);
  write_colour_item(Text_Colour->Hydrophobic_Names,1);
  write_colour_item(Text_Colour->Ligand_Atom_Names,1);
  write_colour_item(Text_Colour->Nonligand_Atom_Names,1);
  write_colour_item(Text_Colour->Hbond_Lengths,1);

  /* Write out the various object sizes */
  fprintf(ligplot_drw,"#S\n");
  for (isize = 0; isize < OBJECT_SIZES; isize++)
    fprintf(ligplot_drw,"%8.3f %s\n",object_size[isize],
            object_size_desc[isize]);

  /* Insert Title and key-text sizes */
  fprintf(ligplot_drw,"   1.000 Title_size\n");
/*  fprintf(ligplot_drw,"   0.090 Key_Text_size\n"); */

  /* Write out the various text sizes */
  for (isize = 0; isize < TEXT_SIZES; isize++)
    fprintf(ligplot_drw,"%8.3f %s\n",text_size[isize],
            text_size_desc[isize]);

  /* If this is a DIMPLOT, write out the DIMPLOT flag */
  if (Interface_Plot == TRUE)
    fprintf(ligplot_drw,"#D\n");
}
/* <--v.4.1 */
/* v.4.3--> */
/***********************************************************************

get_line_type  -  Determine the line-type of the .drw lines that follow

***********************************************************************/

int get_line_type(char identifier)
{
  int icode, line_type;

  /* Initialise variables */
  line_type = -1;

  /* Determine which hash code this is */
  for (icode = 0; icode < (int) strlen(hash_codes); icode++)
    {
      if (identifier == hash_codes[icode])
        line_type = icode;
    }

  /* Return the line-type */
  return(line_type);
}
/***********************************************************************

create_object_record  -  Create a new object record

***********************************************************************/

struct object *create_object_record(struct object **lst_object_ptr)
{
  struct object *object_ptr;
  struct object *last_object_ptr;


  /* Initialise object pointers */
  object_ptr = NULL;
  last_object_ptr = *lst_object_ptr;

  /* Create object entry */
  object_ptr = (struct object*)malloc(sizeof(struct object));
  if (object_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct object\n");
      printf("***        Program ligplot terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */

  /* Initialise the object details */
  object_ptr->nresidues = 0;
  object_ptr->nbonds = 0;
  object_ptr->nrot_bonds = 0;
  object_ptr->nmain_chain = 0;
  object_ptr->n_ca = 0;
  object_ptr->weight = 0;
  object_ptr->minx = 0.0;
  object_ptr->miny = 0.0;
  object_ptr->maxx = 0.0;
  object_ptr->maxy = 0.0;
  object_ptr->interface = 0;
  object_ptr->have_anchors = FALSE;
/* v.5.0--> */
  object_ptr->fixed = FALSE;
/* <--v.5.0 */
/* v.5.3.4--> */
  object_ptr->coords_fixed = FALSE;
/* <--v.5.3.4 */
  object_ptr->max_atom_size = 0.0;
  object_ptr->place_x = 0.0;
  object_ptr->place_y = 0.0;
  object_ptr->internal_energy = 0.0;
  object_ptr->total_energy = 0.0;
  object_ptr->first_residue_ptr = NULL;
  object_ptr->first_object_bond_ptr = NULL;
  object_ptr->next_object_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_object_ptr == NULL)
    first_object_ptr = object_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_object_ptr->next_object_ptr = object_ptr;
                      
  /* Save the current residue's pointer */
  last_object_ptr = object_ptr;

  /* Return current pointers */
  *lst_object_ptr = last_object_ptr;
  return(object_ptr);
}
/***********************************************************************

create_residue_record  -  Create a new residue record

***********************************************************************/

struct residue *create_residue_record(struct residue **lst_residue_ptr,
                                      char res_name[4],char res_num[6],
                                      char chain)
{
  struct residue *residue_ptr;
  struct residue *last_residue_ptr;


  /* Initialise residue pointers */
  residue_ptr = NULL;
  last_residue_ptr = *lst_residue_ptr;

  /* Create residue entry */
  residue_ptr = (struct residue*)malloc(sizeof(struct residue));
  if (residue_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct residue\n");
      printf("***        Program ligplot terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(residue_ptr->res_name,res_name);
  strcpy(residue_ptr->res_num,res_num);
  residue_ptr->chain = chain;

  /* Initialise the remaining residue details */
  residue_ptr->residue_type = STANDARD;
  residue_ptr->natoms = 0;
  residue_ptr->inligand = FALSE;
  residue_ptr->minx = 0.0;
  residue_ptr->maxx = 0.0;
  residue_ptr->miny = 0.0;
  residue_ptr->maxy = 0.0;
  residue_ptr->max_atom_size = 0.0;
  residue_ptr->deleted = FALSE;
  residue_ptr->original_mean_x = 0.0;
  residue_ptr->original_mean_y = 0.0;
  residue_ptr->original_mean_z = 0.0;
  residue_ptr->flattened_mean_x = 0.0;
  residue_ptr->flattened_mean_y = 0.0;
  residue_ptr->have_anchor = FALSE;
  residue_ptr->anchor_pstn_x = 0.0;
  residue_ptr->anchor_pstn_y = 0.0;
  residue_ptr->nanchored_atoms = 0;
  residue_ptr->label_x = 0.0;
  residue_ptr->label_y = 0.0;
  residue_ptr->object_ptr = NULL;
  residue_ptr->first_atom_ptr = NULL;
  residue_ptr->next_residue_ptr = NULL;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_residue_ptr == NULL)
    first_residue_ptr = residue_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_residue_ptr->next_residue_ptr = residue_ptr;
                      
  /* Save the current residue's pointer */
  last_residue_ptr = residue_ptr;

  /* Return current pointers */
  *lst_residue_ptr = last_residue_ptr;
  return(residue_ptr);
}
/***********************************************************************

get_drw_residue  -  Extract residue details from the given .drw record

***********************************************************************/

void get_drw_residue(char input_line[LINELEN + 1],
                     struct residue **lst_residue_ptr,
                     struct object *object_ptr)
{
  char chain, number_string[20], res_name[4], res_num[6], type[2];

  int len, object_type;

  float label_x, label_y;

  struct residue *residue_ptr, *last_residue_ptr;

  /* Initialise variables */
  chain = ' ';
  label_x = 0.0;
  label_y = 0.0;

  /* Initialise pointers */
  residue_ptr = NULL;
  last_residue_ptr = *lst_residue_ptr;

  /* Get the line length */
  len = strlen(input_line);

  /* Proceed only if line is long enough to contain the required data */
  if (len > 9)
    {
      /* Get the residue type */
      type[0] = input_line[0];
      type[1] = '\0';
      object_type = atoi(type); 

      /* Extract the residue details */
      strncpy(res_name,input_line+1,3);
      res_name[3] = '\0';
      strncpy(res_num,input_line+4,5);
      res_num[5] = '\0';

      /* Extract chain identifier */
      chain = input_line[9];
      if (chain == '-')
        chain = ' ';

      /* If have residue name details, then store */
      if (len > 15)
        {
          /* Extract x-coordinate for residue name */
          strncpy(number_string,input_line+11,8);
          number_string[8] = '\0';
          label_x = (float) atof(number_string);

          /* Repeat for y-coordinate */
          strncpy(number_string,input_line+20,8);
          number_string[8] = '\0';
          label_y = (float) atof(number_string);
        }

      /* Create a new residue record */
      residue_ptr = create_residue_record(&last_residue_ptr,res_name,
                                          res_num,chain);

      /* Determine whether residue is a ligand residue */
/* v.5.3.8--> */
      // if (object_ptr->object_type == LIGAND)
      if (object_type == LIGAND)
/* <--v.5.3.8 */
        residue_ptr->inligand = TRUE;

      /* Store the label coordinates */
      residue_ptr->label_x = label_x;
      residue_ptr->label_y = label_y;

      /* Store pointer to residue's object */
      residue_ptr->object_ptr = object_ptr;

      /* Update the details of the object to which this residue
         belongs */
      if (object_ptr->first_residue_ptr == NULL)
        object_ptr->first_residue_ptr = residue_ptr;

      /* Increment count of this object's residues */
      object_ptr->nresidues++;

      /* Store the object type */
      object_ptr->object_type = object_type;
    }

  /* Return current pointer */
  *lst_residue_ptr = last_residue_ptr;
}
/***********************************************************************

create_atom_record  -  Create a new atom record

***********************************************************************/

struct coordinate *create_atom_record(struct coordinate **lst_atom_ptr,
                                      char atom_name[5],
                                      char atom_number[6],float x,
                                      float y,float z)
{
/* v.5.2 --> */
  int i;
/* <-- v.5.2 */

  struct coordinate *atom_ptr;
  struct coordinate *last_atom_ptr;


  /* Initialise atom pointers */
  atom_ptr = NULL;
  last_atom_ptr = *lst_atom_ptr;

  /* Create atom entry */
  atom_ptr = (struct coordinate*)malloc(sizeof(struct coordinate));
  if (atom_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct coordinate\n");
      printf("***        Program ligplot terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  strcpy(atom_ptr->atom_name,atom_name);
  strcpy(atom_ptr->print_name,atom_name);
  strcpy(atom_ptr->atom_number,atom_number);
/* v.4.3--> */
  strcpy(atom_ptr->element,"  ");
/* <--v.4.3 */
  atom_ptr->x = x;
  atom_ptr->y = y;
  atom_ptr->z = z;
  atom_ptr->original_x = x;
  atom_ptr->original_y = y;
  atom_ptr->original_z = z;
  atom_ptr->fit_x = x;
  atom_ptr->fit_y = y;
  atom_ptr->fit_z = z;
  atom_ptr->save_x = x;
  atom_ptr->save_y = y;
  atom_ptr->accessibility = 0.0;
  atom_ptr->atom_size = 0.0;
  atom_ptr->deleted = FALSE;

  /* Mark all residue side chain atoms*/
  if (!strncmp(atom_ptr->atom_name," C ",3)||
      !strncmp(atom_ptr->atom_name," O ",3)||
      !strncmp(atom_ptr->atom_name," CA ",4)||
      !strncmp(atom_ptr->atom_name," N ",3)||
      !strncmp(atom_ptr->atom_name," P ",3)||
      !strncmp(atom_ptr->atom_name," ON ",4)||
      !strncmp(atom_ptr->atom_name," OXT",4))
    atom_ptr->side_chain = FALSE;
  else
    atom_ptr->side_chain = TRUE;

  /* Initialise the remaining atom details */
  strcpy(atom_ptr->occupancy,"  1.00");
  strcpy(atom_ptr->bvalue," 10.00");
  atom_ptr->checked = FALSE;
  atom_ptr->natom_links = 0;
  atom_ptr->plot_atom = TRUE;
  atom_ptr->have_anchor = FALSE;
  atom_ptr->anchor_pstn_x = 0.0;
  atom_ptr->anchor_pstn_y = 0.0;
/* v.5.0--> */
  atom_ptr->anchor_weight = 1.0;
/* v.5.2.1 --> */
  atom_ptr->anchor_weight = 0.0;
/* <-- v.5.2.1 */
/* <--v.5.0 */
/* v.5.1--> */
  atom_ptr->ring_atom = FALSE;
/* <--v.5.1 */
  atom_ptr->atom_count = 0;
  atom_ptr->label_x = 0.0;
  atom_ptr->label_y = 0.0;
  atom_ptr->next_stack_ptr = NULL;
  atom_ptr->first_atom_link_ptr = NULL;
  atom_ptr->next = NULL;
/* v.4.3--> */
  atom_ptr->residue_ptr = NULL;
/* <--v.4.3 */
/* v.5.2 --> */
  for (i = 0; i < S_NTYPES; i++)
    atom_ptr->score[i] = -1;
  atom_ptr->have_score = FALSE;
/* <-- v.5.2 */

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_atom_ptr == NULL)
    first_atom_ptr = atom_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_atom_ptr->next = atom_ptr;
                      
  /* Save the current atom's pointer */
  last_atom_ptr = atom_ptr;

  /* Return current pointers */
  *lst_atom_ptr = last_atom_ptr;
  return(atom_ptr);
}
/* v.4.3--> */
/***********************************************************************

check_for_metal  -  Check whether the given atom is a metal ion

***********************************************************************/

int check_for_metal(char atom_name[5],char element[3])
{
  int imetal;
  int metal;

  static int nmetals = 100;

  /* Note "AC  " changed to "ac  " below as carbon atoms in certain
     molecules labelled AC1, AC1, etc (eg in FAD). Similarly NP, NO
     and PO */
  static char *metal_name[] = {
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

  /* Initialise variables */
  metal = FALSE;

  /* Loop through all the metals */
  for (imetal = 0; imetal < nmetals && metal == FALSE; imetal++)
    {
      /* If name matches in first two characters, and third character
         is either a space or a number, then have a metal */
      if (!strncmp(atom_name,metal_name[imetal],2) &&
          (atom_name[2] == ' ' ||
           (atom_name[2] >= '0' && atom_name[2] <= '9')))
        {
          /* Set flag */
          metal = TRUE;

          /* Return the metal name matched */
          strncpy(element,metal_name[imetal],2);
          element[2] = '\0';
        }
    }

  /* Return the answer */
  return(metal);
}
/***********************************************************************

is_hydrogen_atom  -  Check whether this is a hydrogen or a pseudo atom
                     from any of the many possible representations of
                     such in the PDB

***********************************************************************/

int is_hydrogen_atom(char atname[5])
{
  char atom_name[5], element[3];

  int is_hydrogen, is_metal;

  /* Initialise */
  is_hydrogen = FALSE;
  strncpy(atom_name,atname,4);
  atom_name[4] = '\0';

  /* Check the atom name for any indicators to suggest this is either
     a hydrogen or a pseudo atom */
  if (atom_name[0] == 'H' || atom_name[1] == 'H' || atom_name[1] == 'Q' ||
      (atom_name[1] == 'D' && (atom_name[0] == ' ' || atom_name[0] == '1' ||
                               atom_name[0] == '2' || atom_name[0] == '3')))
    is_hydrogen = TRUE;

  /* Check for possible strange hydrogens */
  if (atom_name[0] >= '0' && atom_name[0] <= '9' && atom_name[2] == 'H')
    is_hydrogen = TRUE;

  /* If this looks likely to be a hydrogen, check that it isn't in
     fact a metal */
  if (is_hydrogen == TRUE)
    {
      /* Check whether this is a metal */
      is_metal = check_for_metal(atom_name,element);

      /* If a metal, then isn't a hydrogen */
      if (is_metal == TRUE)
        is_hydrogen = FALSE;
    }

  return is_hydrogen;
}
/***********************************************************************

get_element  -  Determine element from atom's name

***********************************************************************/

void get_element(char atom_name[5],char element[3])
{
  char metal_name[3];

  int ielem;
  int found, is_hydrogen, is_metal;

  static char *valid_elem[] = { " C", " N", " O", " P", " S" };

  static int nelem = 5;

  /* Initialise variables */
  strcpy(element,"  ");

  /* Check whether this is a hydrogen atom */
  is_hydrogen = is_hydrogen_atom(atom_name);

  /* If atom is a hydrogen, then have its element type */
  if (is_hydrogen == TRUE)
    strcpy(element," H");

  /* Otherwise, identify non-hydrogen element */
  else
    {
      /* First, whether this might be a metal */
      is_metal = check_for_metal(atom_name,metal_name);

      /* If it is, then store the element */
      if (is_metal == TRUE)
        strcpy(element,metal_name);

      /* Otherwise, take character in 2nd position to be representative
         of the element */
      else
        {
          /* Form the proposed element name */
          element[0] = ' ';
          element[1] = atom_name[1];
          element[2] = '\0';

          /* Check that this is a valid element name */
          found = FALSE;
          for (ielem = 0; ielem < nelem && found == FALSE; ielem++)
            {
              /* If have a match, then can stop looking */
              if (!strcmp(element,valid_elem[ielem]))
                found = TRUE;
            }

          /* If not found, then show warning message */
          if (found == FALSE)
            {
              /* Show warning */
              printf("*** Warning. Can't identify element for atom");
              printf(" [%s]\n",atom_name);
            }
        }
    }
}
/* <--v.4.3 */
/***********************************************************************

get_tokens  -  Extract all the token_strings (ie separate char-strings)
               from the line just read in

***********************************************************************/

int get_tokens(char input_line[LINELEN + 1],
               char token_string[MAXTOKENS][TOKENLEN + 1])
{
  int end_line, ipos, ipstn, token, ntokens;


  /* Initialise variables for extraction of information from
     the line just read in */
  end_line = FALSE;
  ipstn = 0;
  ipos = 0;
  for (token = 0; token < MAXTOKENS; token++)
    token_string[token][0] = '\0';
  token = -1;

  /* Search the line for each token_string and store */
  for (ipos = 0; ipos < LINELEN && end_line == FALSE; ipos++)
    {
      /* Check for end of line or end of string */
      if (input_line[ipos] == '\n' || input_line[ipos] == '\0')
        {
          input_line[ipos] = ' ';
          end_line = TRUE;
        }

      /* If have a space, then check if have just ended a token_string */
      if (input_line[ipos] == ' ' || input_line[ipos] == '\t')
        {
          /* If currently saving a token_string, then have reached its end */
          if (ipstn > 0)
            token_string[token][ipstn] = '\0';
          ipstn = 0;
        }

      /* Otherwise, if this is a character, then store in the current
         token_string */
      else
        {
          /* If this is the start of a new token_string, increment
             token_string-count */
          if (ipstn == 0)
            token = token + 1;

          /* Process providing we haven't exceeded number of token_strings
             allowed */
          if (token < MAXTOKENS)
            {
              /* Store the current character in the token_string */
              if (ipstn < TOKENLEN)
                token_string[token][ipstn] = input_line[ipos];

              /* Increment position in the current token_string */
              ipstn++;
            }
        }
    }

  /* Save the number of token_strings read in fom this line */
  ntokens = token + 1;
  if (ntokens > MAXTOKENS)
    ntokens = MAXTOKENS;
  return(ntokens);
}
/***********************************************************************

get_drw_atom  -  Extract atom details from the given .drw record

***********************************************************************/

void get_drw_atom(char input_line[LINELEN + 1],
                  struct coordinate **lst_atom_ptr,int natoms,
/* v.5.3.2 --> */
//                  struct residue *residue_ptr)
                  struct residue *residue_ptr,int ligplus)
/* <--v.5.3.2 */
{
  char atom_number[6], atom_name[5], element[3], number_string[20];
/* v.5.2 --> */
  char token[MAXTOKENS][TOKENLEN + 1];

  char *string_ptr;

  int itoken, ntokens, score;
/* <-- v.5.2 */
  int len;
/* v.5.3.2 --> */
  int ok;
/* <--v.5.3.2 */

  float atom_x, atom_y, label_x, label_y;

  struct coordinate *atom_ptr, *last_atom_ptr;

  /* Initialise variables */
  atom_x = 0.0;
  atom_y = 0.0;
  label_x = 0.0;
  label_y = 0.0;
/* v.5.3.2 --> */
  ntokens = 0;
  ok = FALSE;
/* <--v.5.3.2 */

  /* Initialise pointers */
  atom_ptr = NULL;
  last_atom_ptr = *lst_atom_ptr;

  /* Get the line length */
  len = strlen(input_line);

/* v.5.3.2 --> */
  /* If file created by LigPlus, need to read in the space-separated
     tokens from the line */
  if (ligplus == TRUE)
    {
      /* Extract the tokens */
      ntokens = get_tokens(input_line,token);

      /* If at least 3 tokens, then store */
      if (ntokens > 2)
        {
          /* Extract the atom name */
          strncpy(atom_name,input_line,4);
          atom_name[4] = '\0';

          /* Extract x- and y-coordinates for atom */
          atom_x = (float) atof(token[1]);
          atom_y = (float) atof(token[2]);

          /* If have more than 4 tokens, take them to be atom label
             and label coordinates */
          if (ntokens > 4 && strcmp(token[3],"#"))
            {
              /* Get label's x- and y-coordinates */
              label_x = (float) atof(token[3]);
              label_x = (float) atof(token[4]);
            }

          /* Set flag that atom details OK */
          ok = TRUE;
        }
    }
/* <--v.5.3.2 */

  /* Proceed only if line is long enough to contain the required data */
/* v.5.3.2 --> */
  //if (len > 21)
  else if (len > 21)
/* <--v.5.3.2 */
    {
      /* Extract the atom name */
      strncpy(atom_name,input_line,4);
      atom_name[4] = '\0';

      /* Extract x-coordinate for atom */
      strncpy(number_string,input_line+5,8);
      number_string[8] = '\0';
      atom_x = (float) atof(number_string);

      /* Repeat for y-coordinate */
      strncpy(number_string,input_line+14,8);
      number_string[8] = '\0';
      atom_y = (float) atof(number_string);

      /* If line contains coords for atom label, read these in */
      if (len > 38)
        {
          /* Extract x-coordinate for atom label */
          strncpy(number_string,input_line+23,8);
          number_string[8] = '\0';
          label_x = (float) atof(number_string);

          /* Repeat for y-coordinate */
          strncpy(number_string,input_line+32,8);
          number_string[8] = '\0';
          label_y = (float) atof(number_string);
        }

/* v.5.3.2 --> */
      /* Set flag that atom details OK */
      ok = TRUE;
    }

  /* If atom details OK, then save */
  if (ok == TRUE)
    {
/* <--v.5.3.2 */
      /* Convert atom number to a string */
      sprintf(atom_number,"%5d",natoms);

      /* Create a new atom record */
      create_atom_record(&last_atom_ptr,atom_name,atom_number,atom_x,
                         atom_y,0.0);
      atom_ptr = last_atom_ptr;

/* v.4.4.4--> */
      /* Get the element corresponding to this atom name */
      get_element(atom_name,element);

      /* Store element */
      strcpy(atom_ptr->element,element);

      /* If printing in chemical notation, replace atom name
         by element */
      if (Include->Chemical_Notation == TRUE)
        {
          /* If carbon, then no name required */
          if (!strcmp(element," C"))
            atom_ptr->print_name[0] = '\0';

          /* Otherwise, use element name */
          else
            strcpy(atom_ptr->print_name,element);
        }
/* <--v.4.4.4 */

      /* Update atom data */
      atom_ptr->label_x = label_x;
      atom_ptr->label_y = label_y;
      atom_ptr->atom_count = natoms;

      /* Store pointer to atom's residue */
      atom_ptr->residue_ptr = residue_ptr;

      /* Update residue details */
      if (residue_ptr->first_atom_ptr == NULL)
        residue_ptr->first_atom_ptr = atom_ptr;
      residue_ptr->natoms++;

/* v.5.2 --> */
      /* Check for any LigSearch annotations at the end of the line */
      string_ptr = strchr(input_line,'%');

      /* If have LigSearch annotations, store the details */
      if (string_ptr != NULL)
        {
          /* Extract the tokens at the end of the line */
          ntokens = get_tokens(string_ptr + 1,token);

          /* If have the right number of tokens, then store their
             values */
          if (ntokens >= S_NTYPES)
            {
              /* Loop over the tokens */
              for (itoken = 0; itoken < S_NTYPES; itoken++)
                {
                  /* Get the score */
                  score = atoi(token[itoken]);

                  /* Store the score */
                  atom_ptr->score[itoken] = score;
                  atom_ptr->have_score = TRUE;
                }
            }
        }
/* <-- v.5.2 */
    }

  /* Return current pointer */
  *lst_atom_ptr = last_atom_ptr;
}
/***********************************************************************

create_bond_record  -  Create a new bond record

***********************************************************************/

struct bond *create_bond_record(struct bond *first_bond_ptr,
                                int bond_type,int bond_source,
                                float bond_length,
                                struct coordinate *first_atom_ptr,
                                struct coordinate *second_atom_ptr)
{
  struct bond *bond_ptr;


  /* Initialise bond pointers */
  bond_ptr = NULL;

  /* Create bond entry */
  bond_ptr = (struct bond*)malloc(sizeof(struct bond));
  if (bond_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct bond\n");
      printf("***        Program ligplot terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  bond_ptr->first_atom_ptr = first_atom_ptr;
  bond_ptr->second_atom_ptr = second_atom_ptr;
  bond_ptr->bond_type = bond_type;

  /* Initialise the remaining bond details */
  if (bond_type == COVALENT)
    bond_ptr->elastic = FALSE;
  else
    bond_ptr->elastic = TRUE;
  bond_ptr->rotatable_bond = FALSE;
  bond_ptr->ring_bond = FALSE;
  bond_ptr->end_bond = FALSE;
  bond_ptr->bond_order = SINGLE;
  bond_ptr->bond_source = bond_source;
  bond_ptr->bond_length = bond_length;
  bond_ptr->checked = FALSE;
  bond_ptr->flattened = FALSE;
  bond_ptr->nfirst_atom_links = 0;
  bond_ptr->nsecond_atom_links = 0;
  bond_ptr->nbonds_from_first_atom = 0;
  bond_ptr->nbonds_from_second_atom = 0;
  bond_ptr->first_bond_link_ptr = NULL;
  bond_ptr->next_stack_ptr = NULL;
  bond_ptr->next_link_ptr = NULL;
  bond_ptr->next_bond_ptr = first_bond_ptr;
  bond_ptr->label[0] = '\0';
/* v.4.4--> */
  bond_ptr->plotted = FALSE;
/* <--v.4.4 */
/* v.4.5.2--> */
  bond_ptr->ring_no = 0;
/* <--v.4.5.2 */

  /* Return current pointers */
  return(bond_ptr);
}
/***********************************************************************

get_drw_bond  -  Extract bond details from the given .drw record

***********************************************************************/

void get_drw_bond(char input_line[LINELEN + 1],
                  struct bond **lst_bond_ptr)
{
  char type[2], token[MAXTOKENS][TOKENLEN + 1];

  int atom1, atom2, bond_order, bond_source, bond_type, len, ntokens;

  float bond_length;

  struct bond *bond_ptr, *last_bond_ptr;
  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;

  /* Initialise variables */
  bond_length = 0.0;
  bond_order = SINGLE;
  bond_source = DRW_FILE;

  /* Initialise pointers */
  atom1_ptr = NULL;
  atom2_ptr = NULL;
  bond_ptr = NULL;
  last_bond_ptr = *lst_bond_ptr;

  /* Get the line length */
  len = strlen(input_line);

  /* Proceed only if line is long enough to contain the required data */
  if (len > 4)
    {
      /* Extract the bond type */
      type[0] = input_line[0];
      type[1] = '\0';
      bond_type = atoi(type);
      bond_order = SINGLE;

      /* Extract all the tokens on the input line */
      ntokens = get_tokens(input_line,token);

      /* Extract bond's two atoms */
      atom1 = atoi(token[1]);
      atom2 = atoi(token[2]);

      /* Get the pointers for these two atoms */
      atom_ptr = first_atom_ptr;

      /* Loop through all the atoms to find the matching ones */
      while (atom_ptr != NULL && (atom1_ptr == NULL || atom2_ptr == NULL))
        {
          /* If have a match for the first atom, store pointer */
          if (atom_ptr->atom_count == atom1)
            atom1_ptr = atom_ptr;

          /* Or for second atom */
          else if (atom_ptr->atom_count == atom2)
            atom2_ptr = atom_ptr;

          /* Get the next atom record */
          atom_ptr = atom_ptr->next;
        }

      /* If line contains coords bond label, read it in */
      if (ntokens > 3)
        {
          /* Check whether label simply gives the bond order */
          if (bond_type != HBOND && bond_type != INTERNAL)
            {
              /* Store the bond order */
              if (Include->Chemical_Notation == TRUE)
                {
                  if (token[3][0] == 'd')
                    bond_order = DOUBLE;
                  else if (token[3][0] == 't')
                    bond_order = TRIPLE;
/* v.4.5--> */
                  else if (token[3][0] == 'a')
                    bond_order = AROMATIC;
/* <--v.4.5 */
                }
            }

          /* Otherwise, label gives the H-bond length */
          else
            {
              /* Extract bond length */
              bond_length = (float) atof(token[3]);
            }
        }

      /* If have a bond between two atoms, store it */
      if (atom1_ptr != NULL && atom2_ptr != NULL)
        {
          /* Create a new bond record */
          bond_ptr = create_bond_record(first_bond_ptr,bond_type,
                                        bond_source,bond_length,
                                        atom1_ptr,atom2_ptr);

          /* Store the bond order */
          bond_ptr->bond_order = bond_order;

          /* Store this pointer as the last-encountered bond */
          first_bond_ptr = bond_ptr;
        }
    }

  /* Return current pointer */
  *lst_bond_ptr = last_bond_ptr;
}
/***********************************************************************

get_drw_title  -  Extract title from the given .drw record

***********************************************************************/

void get_drw_title(char input_line[LINELEN + 1])
{
  char token[MAXTOKENS][TOKENLEN + 1];

  int itoken, len, ntokens;

/*  float title_x, title_y; */

  /* Initialise variables */
/*  title_x = 0.0;
  title_y = 0.0; */

  /* Get the line length */
  len = strlen(input_line);

  /* Proceed only if line is long enough to contain the required data */
  if (len > 4)
    {
      /* Extract all the tokens on the input line */
      ntokens = get_tokens(input_line,token);

      /* If line contains coords title label, read it in */
      if (ntokens > 2)
        {
          /* Take last two strings to be the title coords */
/*          title_x = (float) atof(token[ntokens - 2]);
          title_y = (float) atof(token[ntokens - 1]); */

          /* Concatenate remaining tokens to make up the plot title */
          Print_Title[0] = '\0';

          /* Loop through the tokens and add on */
          for (itoken = 0; itoken < ntokens - 2; itoken++)
            {
              strcat(Print_Title,token[itoken]);
              strcat(Print_Title," ");
            }

          /* Remove trailing blank */
          len = strlen(Print_Title);
          if (len > 0)
            Print_Title[len - 1] = '\0';
        }
    }
}
/***********************************************************************

read_drw_file  -  Read in the ligplot.drw file to extract all the objects,
                  residues, atom and bonds to be plotted

***********************************************************************/

int read_drw_file(char drw_name[FILENAME_LEN])
{
  int natoms, nbonds, nobjects, nresid;

  char input_line[LINELEN + 1];

  int len, line, line_type;
/* v.5.3.2 --> */
  int ligplus;
/* <--v.5.3.2 */

  struct coordinate *last_atom_ptr;
  struct bond *last_bond_ptr;
  struct object *object_ptr, *last_object_ptr;
  struct residue *last_residue_ptr;

  FILE *file_ptr;

  /* Initialise variables */
/* v.5.3.2 --> */
  ligplus = FALSE;
/* <--v.5.3.2 */
  line = 0;
  line_type = -1;
  natoms = 0;
  nbonds = 0;
  nobjects = 0;
  nresid = 0;

  /* Initialise pointers */
  first_bond_ptr = NULL;
  last_atom_ptr = NULL;
  last_object_ptr = NULL;
  last_residue_ptr = NULL;
  object_ptr = NULL;

  /* Open .drw file to see if it is present */
  if ((file_ptr = fopen(drw_name,"r")) ==  NULL)
    {
      printf("\n*** Unable to open .drw file [%s]\n",drw_name);
      exit(1);
    }

  /* Read in the data from the file */
  printf("\nReading in the .drw file [%s] ...\n",drw_name);

  /* Loop through the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Truncate the string to strip off the newline */
      len = string_truncate(input_line,LINELEN);

      /* Process input line only if not empty */
      if (len > 0)
        {
          /* If this is a keyword line, set flag according to the
             keyword given by #x, where x is one of the hash
             codes: M, C, S, R, A, B */
          if (input_line[0] == '#')
            {
              /* Get the keyword identifier if there is one */
              if (len > 1)
                {
                  /* Identify the identifier */
                  line_type = get_line_type(input_line[1]);

/* v.5.3.2 --> */
                  /* Check whether line shows that file was created
                     by LigPlot+ */
                  if (!strncmp(input_line,"#--- Created by: LigPlus",24))
                    ligplus = TRUE;
/* <--v.5.3.2 */
                }
            }

          /* Otherwise, process the relevant line */
          else
            {
              /* Residue name, number and type */
              if (line_type == RES_DATA)
                {
                  /* Extract residue details and increment count */
                  get_drw_residue(input_line,&last_residue_ptr,object_ptr);
                  nresid++;
                }

              /* Atom details */
              else if (line_type == ATOM_DATA)
                {
                  /* Extract atom details and increment count */
                  get_drw_atom(input_line,&last_atom_ptr,natoms,
/* v.5.3.2 --> */
//                               last_residue_ptr);
                               last_residue_ptr,ligplus);
/* <--v.5.3.2 */
                  natoms++;
                }

              /* Bond details */
              else if (line_type == BOND_DATA)
                {
                  /* Extract bond details and increment count */
                  get_drw_bond(input_line,&last_bond_ptr);
                  nbonds++;
                }

              /* If this is the plot title, extract it and its coords */
              else if (line_type == PLOT_TITLE)
                get_drw_title(input_line);
            }

          /* If line-type is a new object, the create object record */
          if (line_type == OBJECT_DATA)
            {
              /* Create object record and increment count */
              object_ptr = create_object_record(&last_object_ptr);
              nobjects++;
            }
        }

      /* Increment line-count */
      line++;
    }

  /* Close the .drw file */
  fclose(file_ptr);

  /* Return count of objects read in */
  return(nobjects);
}
/* <--v.4.3 */
/***********************************************************************

getnam  -  Peel off the directory path and extension from the full name
           of the .pdb file

***********************************************************************/

int getnam(char pdbfil[FILENAME_LEN],int *istart,int *iend)
{
  char pchar;
  int  finish, gotdot, ipos, istate;


  /* Initialise variables */
  finish = FALSE;
  *iend = -1;
  *istart = 0;
  istate = 1;
  ipos = strlen(pdbfil) - 1;
  gotdot = FALSE;

  /* Check through the filename from right to left */
  while (finish == FALSE && ipos > -1)
    {
      /* Pick off next character */
      pchar = pdbfil[ipos];

      /* State 1: Searching for first non-blank character */
      if (istate == 1)
        {
          if (pchar == '/' || pchar == '\\' || pchar == ']')
            {
              printf("*** ERROR in supplied name of file: [%s\n",
                     pdbfil);
              return(-1);
            }

          if (pchar != ' ' && pchar != '.')
            {
              *iend = ipos;
              istate = 2;
            }
        }

      /* State 2: Searching for end of extension, or end of
         directory path */
      else if (istate == 2)
        {
          /* If character is a dot, and is the first dot, then
             note position */
          if (pchar == '.' && gotdot == FALSE)
            {
              *iend = ipos - 1;
              gotdot = TRUE;
            }

          /* If character signifies the end of a directory path,
             note pstn */
          else if (pchar == '/' || pchar == '\\' || pchar == ']')
            {
              *istart = ipos + 1;
              finish = TRUE;
            }
        }

      /* Step back a character */
      ipos = ipos - 1;
    }

  /* Check whether file name is sensible */
  if (*istart > *iend)
    {
      printf("*** ERROR in supplied name of file: %s\n", pdbfil);
      printf("    No name found\n");
      return(-1);
    }

  return(0);
}
/***********************************************************************

get_filename  -  Strip off directory path of PDB filename
                 and append required extention to root

***********************************************************************/

void get_filename(char *pdb_name, char *file_name,char *extension)
{
    int iend, ierror, inew, ipos, istart;

    /* Initialise returned filename */
    file_name[0] = '\0';

    /* Check that PDB filename is valid */
    if (pdb_name[0] != '\0' && pdb_name[0] != '\n')
    {
        /* Peel off the directory path and extension, returning just
           the root file name */
        ierror = getnam(pdb_name,&istart,&iend);
        if (ierror != -1)
        {
            /* Transfer just the root name into new file name */
            for (ipos = istart, inew = 0; ipos < iend + 1;
                 ipos++, inew++)
                file_name[inew] = pdb_name[ipos];
            file_name[inew] = '\0';

            /* Store the root name */
            strcpy(root_name,file_name);

/* v.5.0--> */
            /* Get the file's full pathname, up to the extension */
/* v.5.0.3--> */
/*            strncpy(file_name,pdb_name,iend + 1);
              file_name[iend + 1] = '\0'; */
/* v.5.3.8--> */
            // if (wkdir[0] != '\0')
            if (wkdir[0] != '\0' && highlight_res[0] == '\0')
/* <--v.5.3.8 */
              {
                strcpy(file_name,wkdir);
                strcat(file_name,root_name);
              }
/* <--v.5.0.3 */
/* <--v.5.0 */

            /* Append the extension */
            strcat(file_name,extension);
        }
    }
}
/***********************************************************************

get_residue_number  -  Interpret the residue-number entered by the
                       user

***********************************************************************/

/* v.4.0--> */
/* void get_residue_number(char tmpstr[LINELEN],char res_num[6]) */
int get_residue_number(char tmpstr[LINELEN],char res_num[6],char res_name[4])
/* <--v.4.0 */
{
  char insertion_code, res_number[6], new_number[6];
/* v.4.0--> */
  char temp_string[LINELEN];
/* <--v.4.0 */
  int cpos, got_number, idigit, ipos, ndigits;
/* v.4.0--> */
  int at_end, residue_name, residue_number;
  int len, nchars;
/* <--v.4.0 */

  /* Initialise position of insertion code */
  cpos = 0;
  insertion_code = ' ';
  got_number = FALSE;
  ndigits = 0;
  for (idigit = 0; idigit < 5; idigit++)
    new_number[idigit] = ' ';
/* v.4.0--> */
  at_end = FALSE;
  nchars = 0;
  residue_name = FALSE;
  residue_number = FALSE;

  /* Check whether this has a -n prefix indicating a residue name */
  if (!strncmp(tmpstr,"-n",2))
    {
      residue_name = TRUE;

      /* Excise the -n */
      strcpy(temp_string,tmpstr+2);
      strcpy(tmpstr,temp_string);
    }
/* <--v.4.0 */

  /* Loop through the characters in the entered residue-number
     searching for end of number or insertion-code */
/* v.4.0--> */
/*  for (ipos = 0; ipos < 5 && cpos == 0; ipos++) */
  for (ipos = 0; ipos < 5 && at_end == FALSE; ipos++)
/* <--v.4.0 */
    {
      /* If this is the end of the string, store end-position */
      if (tmpstr[ipos] == '\0')
/* v.4.0--> */
        {
/* <--v.4.0 */
          cpos = ipos;
/* v.4.0--> */
          at_end = TRUE;
        }
/* <--v.4.0 */

      /* If this is a digit, then have part of a residue number */
/* v.4.0--> */
/*      else if (tmpstr[ipos] >= '0' && tmpstr[ipos] <= '9') */
      else if ((tmpstr[ipos] >= '0' && tmpstr[ipos] <= '9') ||
               tmpstr[ipos] == '-')
/* <--v.4.0 */
        {
          res_number[ndigits] = tmpstr[ipos];
          ndigits++;
          got_number = TRUE;

/* v.4.0--> */
          /* If already have a character, then this must be a residue name */
          if (nchars > 0)
            {
              got_number = FALSE;
              residue_name = TRUE;
            }
/* <--v.4.0 */
        }

      /* If this is a character it might be a prefix to the
         residue number, or the insertion code */
      else if ((tmpstr[ipos] >= 'A' && tmpstr[ipos] <= 'Z') ||
               tmpstr[ipos] == '\'')
        {
          /* If we already have a number, this must be an
             insertion code */
          if (got_number == TRUE)
            {
              cpos = ipos;
              insertion_code = tmpstr[ipos];
            }

          /* Otherwise it must be a residue name */
          else
            {
/* v.4.0--> */
/*            res_number[ndigits] = tmpstr[ipos];
              ndigits++; */
              residue_name = TRUE;
              got_number = FALSE;
/* <--v.4.0 */
            }
/* v.4.0--> */
          /* Increment count of characters */
          nchars++;
/* <--v.4.0 */
        }

      /* If this is a space, have reached the end of the number */
      else if (tmpstr[ipos] == ' ')
        {
          if (got_number == TRUE)
            cpos = ipos;
        }      /* Otherwise, have an invalid character in the residue
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

  /* If entry doesn't resemble a residue number, then see whether it
     might be a residue name instead */
/* v.4.0--> */
/*  if (cpos < 0 || got_number == FALSE) */
  if (residue_name == TRUE || cpos < 0 || got_number == FALSE ||
      nchars > 1)
/* <--v.4.0 */
    {
/* v.4.0--> */
/*        printf("*** Invalid residue-number: [%s]\n",tmpstr);
        exit(1);
      } */
      /* Get length of residue name */
      len = strlen(tmpstr);

      /* If too long, then have an error */
      if (len > 3)
        {
          printf("*** Invalid residue-number or name:");
          printf(" [%s]\n",tmpstr);
          exit(1);
        }

      /* Store the name */
      strcpy(res_name,tmpstr);

/* v.4.4--> */
      /* Check for blanks at end */
/* v.5.0--> */
      /* if (len == 3 && !strncmp(tmpstr+1,"  ",2)) */
      if (len == 3 && !strncmp(tmpstr+1,"  ",2) && residue_name == FALSE)
/* <--v.5.0 */
        {
          len = 1;
          tmpstr[1] = '\0';
        }
/* v.5.0--> */
      /* if (len == 3 && tmpstr[2] == ' ') */
      if (len == 3 && tmpstr[2] == ' ' && residue_name == FALSE)
/* <--v.5.0 */
        {
          len = 2;
          tmpstr[2] = '\0';
        }

      /* Check for short residue names */
      if (len == 1)
        {
          strcpy(res_name,"  ");
          strcat(res_name,tmpstr);
        }
      else if (len == 2)
        {
          strcpy(res_name," ");
          strcat(res_name,tmpstr);
        }
/* <--v.4.4 */

      residue_number = FALSE;
    }

  /* Otherwise, reformat the number */
  else
    {
/* <--v.4.0 */
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
/* v.4.0--> */

      /* Set flag to say that residue number is OK */
      residue_number = TRUE;
    }

  /* Return flag indicating whether residue number identified */
  return(residue_number);
/* <--v.4.0 */

}
/***********************************************************************

chk_cofactor  -  Check for a cofactor

***********************************************************************/

int chk_cofactor(int a,char data[MAXCYCLE][4],char pdb_line[81])
{
    int 
        i,
        cofactor;
    cofactor = -1;
    for (i = 0;i < a;i++) 
        if (!strncmp(data[i],pdb_line+17,3))
            cofactor = i;
    return cofactor;
}
/***********************************************************************

print_cofactor  -  Print out the cofactor information

***********************************************************************/

void print_cofactor(char pdb_line[81],char name[4],char number[6],
                    int cofactor,char infoc[MAXCYCLE][41])
{
    infoc[cofactor][40] = '\0';
    printf("\t%s %s %c\t%s\n",name,number,
           pdb_line[21],infoc[cofactor]);
}
/***********************************************************************

first_res  -  Save the first residues details

***********************************************************************/

void first_res(char pdb_line[81],char num_start[6],char *chain_start,
               char name_start[3])
{
    strncpy(num_start,pdb_line+22,5);
    num_start[5] = '\0';
    *chain_start = pdb_line[21];
    strncpy(name_start,pdb_line+17,3);
    name_start[3] = '\0';
}
/***********************************************************************

print_chains  -  Print out the information

***********************************************************************/

void print_chains(char num_start[6],
             char chain_start,  
             char name_start[3],
             char res_num[6],
             char chain_id,
             char res_name[4])
{
    num_start[5] = '\0';
    name_start[3] = '\0';
    res_num[5] = '\0';
    res_name[3] = '\0';

    if ((atoi(res_num) - atoi(num_start)) > 10)
    {
        printf("\t%s %s %c to",
               name_start,num_start,chain_start);
        printf(" %s %s %c \n",res_name,
               res_num,chain_id);
    }
}
/***********************************************************************

ligand_input  -  For the ligand input string 

***********************************************************************/

void ligand_input(char character,
              char ligand2[6],
              char ligand3[6],
              char ligand4[6],
              char ligand5[6])
{
    static int
        num = 2,
        i = 0, 
        copied = OFF;

    if (character != ' ')
    {
        switch (num)
        {
        case(2):
        {
            ligand2[i] = character;
            i++;
            ligand2[i] = '\0';
            copied = ON;
        }
            break;
        case(3):
        {
            ligand3[i] = character;
            i++;
            ligand3[i] = '\0';
        }
            break;
        case(4):
        {
            ligand4[i] = character;
            i++;
            ligand4[i] = '\0';
        }
            break;
        case(5):
        {
            ligand5[i] = character;
            i++;
            ligand5[i] = '\0';
        }
            break;
        }
    }
    else if (character == ' ' &&
            (num != 2 || copied == ON))
    {
        i = 0;
        num++;
    }
    if (!strncmp(ligand4,"-h",2))
    {
        strncpy(ligand5,"-h",2);
        ligand5[2] = '\0';
        ligand4[0] = '\0';
    }
    if (!strncmp(ligand4,"-w",2))
    {
        strncpy(ligand5,"-w",2);
        ligand5[2] = '\0';
        ligand4[0] = '\0';
    }
}

/***********************************************************************

lig_size  -  See how big a ligand chain is

***********************************************************************/

int lig_size(FILE *fil_pdb,char num[5],char chain)
{
    char
        num_temp[6],
        pdb_line[81];
    int
        count;
    static int
        once = 1;
    if (once == 1)
    {
        once++;
        count = 0;
        strncpy(num_temp,"#####",5);
        while (fgets (pdb_line,80,fil_pdb)!= NULL)    
        {
            if (strncmp(num_temp,pdb_line+22,5))
            {
                strncpy(num_temp,pdb_line+22,5);
                if (pdb_line[21] == chain && 
                   strncmp(num,pdb_line+22,5))
                    count++;
                else
                    count = 0;
            }
            if (count == 12)
                return OFF;
        }
        return ON;
    }
    return OFF;
}
/* v.4.0--> */
/***********************************************************************

check_for_dimplot  -  Check the ligplot.pdb file to test whether it
                      contains a flattened DIMPLOT from a previous run

***********************************************************************/

int check_for_dimplot(void)
{
  char line[LINELEN + 1];
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */

  int dimplot, keep_reading;

  FILE *fil_pdb;

  /* Open the ligplot.pdb file */
/* v.5.0--> */
  /* if ((fil_pdb = fopen("ligplot.pdb","r")) == NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.pdb");
  if ((fil_pdb = fopen(filename,"r")) == NULL)
/* <--v.5.0 */
    {
      printf("\n*** Unable to open the ligplot.pdb file\n");
      exit(1);
    }

  /* Initialise variables */
  dimplot = FALSE;
  keep_reading = TRUE;

  while (fgets(line,LINELEN,fil_pdb) != NULL && keep_reading == TRUE)
    {
      /* If this is a DIMPLOT record indicating that this is a flattened
         DIMPLOT, then can stop reading */
      if (!strncmp(line,"DIMPLOT",7))
        {
          /* Set the flags */
          dimplot = TRUE;
          keep_reading = FALSE;
        }

      /* Otherwise, if have hit the first ATOM or HETATM record, stop
         reading */
      else if (!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6))
        keep_reading = FALSE;
    }

  /* Close the ligplot.pdb file */
  fclose(fil_pdb);

  /* Return whether this is a DIMPLOT */
  return(dimplot);
}
/* <--v.4.0 */
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *strptrs[],int num,
                           char pdb_name[FILENAME_LEN],
/* v.4.0--> */
                           char res_name1[4],
/* <--v.4.0 */
                           char res_num1[6],
/* v.4.0--> */
                           char res_name2[4],
/* <--v.4.0 */
                           char res_num2[6],
                           char *chain_identi,
                           char hhb_name[FILENAME_LEN],
                           char nnb_name[FILENAME_LEN],
                           char bonds_name[FILENAME_LEN],
/* v.4.1.2--> */
/*                         char title[TITLE_LEN]) */
/* v.5.4.2--> */
//                           char title[TITLE_LEN],int *debug)
                           char title[TITLE_LEN],int *debug,
                           int *no_abort)
/* <--v.5.4.2 */
/* <--v.4.1.2 */
{
/* v.4.0-->
  int
    lig_name;
<--v.4.0 */
  int irange, itoken;
/* v.4.0--> */
  int res_num;
/* <--v.4.0 */
  char token[MAXTOKENS][FILENAME_LEN];
  char
/* v.4.0--> */
/*    ligand2[5],
    ligand3[5],
    ligand4[5],
    ligand5[5], */
/* <--v.4.0 */
    tmpstr1[LINELEN];
/* v.4.1.1--> */
  char res_name[4], temp[4];
/* v.4.2--> */
  char number_string[20];
/* <--v.4.2 */
/* v.5.4.5 --> */
  char ligmet[20];
/* <-- v.5.4.5 */
/* v.5.1.6 --> */
  char tmp[FILENAME_LEN];
/* <-- v.5.1.6 */
  int fpos, len, loop, ipos;
/* <--v.4.1.1 */
/* v.5.3 --> */
  int number;
/* <--v.5.3 */

  /* Initialise */
  Print_as_is = FALSE;
  title[0] = '\0';
  pdb_name[0] = '\0';
  res_num1[0] = '\0';
  res_num2[0] = '\0';
/* v.4.0--> */
  res_name1[0] = '\0';
  res_name2[0] = '\0';
/* <--v.4.0 */
  *chain_identi = ' ';
  for (itoken = 0; itoken < MAXTOKENS; itoken++)
    token[itoken][0] = '\0';
/* v.5.0--> */
  Include_Waters = FALSE;
  run_type = STANDARD_RUN;
  first_ensemble_ptr = NULL;
/* <--v.5.0 */
/* v.5.1--> */
  Filter_Waters = FALSE;
/* <--v.5.1 */
/* v.5.1.4--> */
  Clip_Diagram = FALSE;
/* <--v.5.1.4 */
/* v.5.4.2--> */
  *no_abort = FALSE;
/* <--v.5.4.2 */
/* v.5.4.5 --> */
  PDBsum1 = FALSE;
/* <-- v.5.4.5 */

  /* If no parameters present, not even filename, then stop */
  if (num == 0)
    {
      printf("*** ERROR. No command-line arguments entered!\n");
      exit(-1);
    }

  /* If only a single parameter entered, then take this to be the filename
     and prompt user for remaining residue information */
  else if (num == 1)
    {
/* v.4.0--> */
/*      lig_name = ON; */
/* <--v.4.0 */

      strcpy(pdb_name,strptrs[1]);
      strcpy(token[0],strptrs[1]);

/* v.4.0--> */
      /* Check whether the file is the output file from dimer, giving
         the residue interactions across a dimer or domain interface */
/* v.5.0--> */
/*      if (!strcmp(pdb_name,"dimplot.pdb")) */
      if (strstr(pdb_name,"dimplot.pdb") != NULL)
/* <--v.5.0 */
        {
          /* Have dimplot file, so set appropriate flag */
          Interface_Plot = TRUE;
/* v.5.3.4--> */
          From_Dimer = TRUE;
/* <--v.5.3.4 */
        }

      /* If this is a ligplot.pdb file, then check whether it contains
         a flattened DIMPLOT from a previous run */
      if (!strcmp(pdb_name,"ligplot.pdb"))
        Interface_Plot = check_for_dimplot();

      /* If this is not a DIMPLOT, ask user for the ligand details */
      if (Interface_Plot == FALSE)
        {
          // lig_search(pdb_name,token,&num);
          printf("*** ERROR. Not enough command-line arguments entered!\n");
          exit(-1);
        }
/* <--v.5.0 */
    }

  /* If two or more parameters, transfer the command-line parameters
     into the array to be processed */
  else
    {
      for (itoken = 0; itoken < num; itoken++)
/* v.5.0--> */
        {
/* <--v.5.0 */
          strcpy(token[itoken],strptrs[itoken + 1]);
/* v.5.0--> */

          /* Get the length of the token */
          len = strlen(token[itoken]);

          /* Check for bounding quotes */
          if (len > 1 && token[itoken][0] == '"' &&
              token[itoken][len - 1] == '"')
            {
              /* Strip out the quotes */
/* v.5.1.6 --> */
              // strcpy(token[itoken],token[itoken]+1);
              strcpy(tmp,token[itoken]+1);
              strcpy(token[itoken],tmp);
/* <-- v.5.1.6 */
              token[itoken][len - 2] = '\0';
            }
        }
/* <--v.5.0 */
    }

  /* Process the entered tokens */

  /* Initialise count of non-flag tokens */
  irange = 0;

  /* Loop over all the command-line arguments */
  for (itoken = 0; itoken < num; itoken++)
    {
      /* Check for the -h option */
      if (!strcmp(token[itoken],"-h"))
        {
/* v.5.0--> */
          printf("*** ERROR. The -h option is no longer supported\n");
          exit(1);

          /* Title for plot to be entered by user - accept it */
          // printf("\nEnter title of your plot\n");
          // gets(title);

          /* Blank out this token */
          // token[itoken][0] = '\0';
/* <--v.5.0 */
        }

/* v.5.3.8--> */
      /* Check for the -res option identifying residue to be highlighted */
      else if (!strcmp(token[itoken],"-res"))
        {
          /* Retrieve the residue name from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(highlight_res,token[itoken]);
            }
        }

      /* Check for the -pdb option giving PDB code */
      else if (!strcmp(token[itoken],"-pdb"))
        {
          /* Retrieve the UniProt accession from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(pdb_code,token[itoken]);
            }
        }

      /* Check for the -pdbsum option giving the PDB entry's PDBsum directory */
      else if (!strcmp(token[itoken],"-pdbsum"))
        {
          /* Retrieve the directory name from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(pdbsum_dir,token[itoken]);

              /* Check that final character is a file separator */
              len = strlen(pdbsum_dir);
              if (pdbsum_dir[len - 1] != '/' && pdbsum_dir[len - 1] != '\\')
                strcat(pdbsum_dir,"/");
            }

          /* If no directory supplied, show error message */
          else
            {
              printf("*** ERROR. No working directory given for -pdbsum");
              printf(" option\n");
              exit(1);
            }
        }
/* <--v.5.3.8 */

/* v.5.0--> */
      /* Check for the -wkdir option */
      else if (!strcmp(token[itoken],"-wkdir"))
        {
          /* Retrieve the directory name from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(wkdir,token[itoken]);

              /* Check that final character is a file separator */
              len = strlen(wkdir);
              if (wkdir[len - 1] != '/' && wkdir[len - 1] != '\\')
                strcat(wkdir,"/");
            }

          /* If no directory supplied, show error message */
          else
            {
              printf("*** ERROR. No working directory given for -wkdir");
              printf(" option\n");
              exit(1);
            }
        }

      /* Check for the -prm option */
      else if (!strcmp(token[itoken],"-prm"))
        {
          /* Retrieve the file name from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(ligplot_prm,token[itoken]);
            }

          /* If no file supplied, show error message */
          else
            {
              printf("*** ERROR. No parameter file name given for -prm");
              printf(" option\n");
              exit(1);
            }
        }

      /* Check for the -ligfit option indicating fitting ligand only (from
         ligfit program in LigPlus) */
      else if (!strcmp(token[itoken],"-ligfit"))
        ligfit = TRUE;

/* v.5.3.4--> */
      /* Check for the -dimfit option indicating fitting DIMPLOT to a
         previous DIMPLOT diagram */
      else if (!strcmp(token[itoken],"-dimfit"))
        dimfit = TRUE;
/* <--v.5.3.4 */

      /* Check for the -ligsearch option indicating atom annotations
         to be shown in the ligand */
      else if (!strcmp(token[itoken],"-ligsearch"))
        {
          /* Retrieve the LigSearch id from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(ligsearch_id,token[itoken]);

              /* Set flag */
              ligsearch_annot = TRUE;
            }
        }

/* v.5.3 --> */
      /* Check for the -ctype option defining the non-bonded contacts to
         be included */
      else if (!strcmp(token[itoken],"-ctype"))
        {
          /* Get the contact type from the next token, if there is one */
          if (itoken < num - 1)
            {
              itoken++;
              number = atoi(token[itoken]);
              if (number > -1 && number < 3)
                contact_type = number;
            }          
        }
/* <--v.5.3 */

      /* Check for the -ps+ option indicating that PostScript output
         required from LigPlot+ input files */
      else if (!strcmp(token[itoken],"-ps+"))
        {
          /* Retrieve the name of the .drw file from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(ligplot_drw_name,token[itoken]);

              /* Set run-type and plot flags */
              run_type = LIGPLUS_PLOT;
              Plot_from_drw = TRUE;
            }

          /* If no file supplied, show error message */
          else
            {
              printf("*** ERROR. No .drw file name given for -ps+");
              printf(" option\n");
              exit(1);
            }
        }
/* <--v.5.0 */

/* v.5.0.1--> */
      /* Check for the -conect option indicating that all CONECT records
         required, irrespective of their length (for schemat program in
         PDBsum) */
      else if (!strcmp(token[itoken],"-conect"))
        conect = TRUE;

/* <--v.5.0.1 */

/* v.4.2--> */
      /* Check for the -xmas option signifying that input file is in
         XMAS format */
      else if (!strcmp(token[itoken],"-xmas"))
        {
          /* If not linked to the XMAS functions, then show warning
             message  */
          if (Have_Xmas == FALSE)
            {
              printf("*** ERROR. -xmas option invalid as program not linked");
                  printf(" to XMAS routines\n");
              exit(1);
            }

          /* Set flag that input file is in XMAS format */
          Read_Xmas = TRUE;
        }

      /* Check for the -xall option signifying that all ligands in the
         XMAS format are to be plotted one after another */
      else if (!strcmp(token[itoken],"-xall"))
        {
          /* If not linked to the XMAS functions, then show warning
             message  */
          if (Have_Xmas == FALSE)
            {
              printf("*** ERROR. -xall option invalid as program not linked");
                  printf(" to XMAS routines\n");
              exit(1);
            }

          /* Set flag that all XMAS ligands to be processed */
          Plot_All_Xmas_Ligands = TRUE;
          Ligand_Molecule_id = 0;

          /* Set flag that input file is in XMAS format */
          Read_Xmas = TRUE;
        }

      /* Check for the -xmol option signifying the molecule number in
         the XMAS file corresponding to the ligand to be plotted */
      else if (!strcmp(token[itoken],"-xmol"))
        {
          /* If not linked to the XMAS functions, then show warning
             message  */
          if (Have_Xmas == FALSE)
            {
              printf("*** ERROR. -xmol option invalid as program not linked");
                  printf(" to XMAS routines\n");
              exit(1);
            }

          /* Retrieve the molecule number from the next token */
          if (itoken < num - 1)
            {
              /* Increment token-number and get next token */
              itoken++;
              strcpy(number_string,token[itoken]);

              /* Get the number corresponding to this string */
              Ligand_Molecule_id = atoi(number_string);

              /* Set flag that input file must be in XMAS format */
              Read_Xmas = TRUE;
            }

          /* If no molecule-id supplied, show error message */
          else
            {
              printf("*** ERROR. No molecule number given for -xmol");
              printf(" option\n");
              exit(1);
            }
        }
/* <--v.4.2 */

/* v.4.4.1--> */
      /* Check for the -hyd option signifying hydrogens to be included */
      else if (!strcmp(token[itoken],"-hyd"))
        Include_Hydrogens = TRUE;

      /* Check for the -gk option signifying atom names to be converted
         to Greek characters, where appropriate */
      else if (!strcmp(token[itoken],"-gk"))
        Greek_Names = TRUE;
/* <--v.4.4.1 */

/* v.5.0--> */
      /* Check for the -wat option signifying waters to be included */
      else if (!strcmp(token[itoken],"-wat"))
          Include_Waters = TRUE;
/* <--v.5.0 */

/* v.5.1--> */
      /* Check for the -wcut option signifying waters to be filtered to
         retain only those that make at least 2 H-bonds to protein/ligand */
      else if (!strcmp(token[itoken],"-wcut"))
        {
          Include_Waters = TRUE;
          Filter_Waters = TRUE;
        }
/* <--v.5.1 */

      /* Check for the -d debug option */
      else if (!strcmp(token[itoken],"-d"))
        *debug = TRUE;

/* v.5.4.2--> */
      /* Check for the -no_abort option */
      else if (!strcmp(token[itoken],"-no_abort"))
        *no_abort = TRUE;
/* <--v.5.4.2 */

      /* Check for the -m option signifying that ligand is a metal */
      else if (!strcmp(token[itoken],"-m"))
        Metal_as_Ligand = TRUE;

      /* Check for the -w option signifying that ligand is a water */
      else if (!strcmp(token[itoken],"-w"))
        Water_as_Ligand = TRUE;

      /* Check for the -r option signifying that ligplot.res file
         to be generated */
      else if (!strcmp(token[itoken],"-r"))
        Write_Res_File = TRUE;

/* v.4.3--> */
      /* Check for the -l option signifying that ligplot.lig and
         ligplot.res files to be generated */
      else if (!strcmp(token[itoken],"-l"))
        {
          Write_Lig_File = TRUE;
          Write_Res_File = TRUE;
        }
/* <--v.4.3 */

/* v.4.5.4--> */
      /* Check for the -t option signifying that bond types to be
         written out to the ligplot.drw file */
      else if (!strcmp(token[itoken],"-t"))
        Write_Bond_Types = TRUE;
/* <--v.4.5.4 */

/* v.4.4.4--> */
      /* Check for the -scale option signifying that molecules on plot
         to be drawn to fixed scale (used for Enzyme Structures
         Database and PDBsum only) */
      else if (!strcmp(token[itoken],"-scale"))
        Scale_Plot = TRUE;
/* <--v.4.4.4 */

/* v.5.1.4--> */
      /* Check for the -clip option signifying that the PostScript file's
         bounding box is to clip the picture */
      else if (!strcmp(token[itoken],"-clip"))
        Clip_Diagram = TRUE;
/* <--v.5.1.4 */

      /* Check for the -x option signifying that ligand to be plotted
         exactly as it is, without any flattening */
      else if (!strcmp(token[itoken],"-x"))
        Print_as_is = TRUE;

/* v.5.3.4--> */
      /* Check for the -norot option signifying that no final rotation
         to be applied to the plot */
      else if (!strcmp(token[itoken],"-norot"))
        Picture->Allow_Rotation = FALSE;
/* <--v.5.3.4 */

/* v.5.4.5 --> */
      /* Check for the -pdbsum1 option signifying program part of
         PDBsum1 */
      else if (!strcmp(token[itoken],"-pdbsum1"))
        PDBsum1 = TRUE;
/* <-- v.5.4.5 */

      /* Otherwise, take the token to be the filename or part of
         the residue-range, according to order of entry */
      else if (token[itoken][0] != '\0')
        {
          /* If this is the first non-flag parameter, then it must
             be the file name */
          if (irange == 0)
            {
              strcpy(pdb_name,token[itoken]);

/* v.4.0--> */
              /* Check whether the file is dimplot.pdb */
/* v.5.0--> */
/*              if (!strcmp(pdb_name,"dimplot.pdb")) */
              if (strstr(pdb_name,"dimplot.pdb") != NULL)
/* <--v.5.0 */
/* v.5.3.4--> */
                {
/* <--v.5.3.4 */
                  Interface_Plot = TRUE;
/* v.5.3.4--> */
                  From_Dimer = TRUE;
                }
/* <--v.5.3.4 */

              /* If this is a ligplot.pdb file, then check whether it
                 contains a flattened DIMPLOT from a previous run */
              if (!strcmp(pdb_name,"ligplot.pdb"))
                Interface_Plot = check_for_dimplot();
/* <--v.4.0 */
            }

/* v.4.0--> */
          /* If not an interface plot, check for residue range and
             chain-id information */
          if (Interface_Plot == FALSE)
            {
/* <--v.4.0 */
          /* If this is the first of the residue-range parameters,
             then take it to be the first residue in the range */
/* v.4.0--> */
/*        else if (irange == 1) */
              if (irange == 1)
/* <--v.4.0 */
                {
                  /* Interpret the entered residue number */
                  strcpy(tmpstr1,token[itoken]);
/* v.4.0--> */
/*                get_residue_number(tmpstr1,res_num1); */
                  res_num = get_residue_number(tmpstr1,res_num1,res_name1);

                  /* If not a residue number, then decrement range as
                     still need to pick up the residue number */
                  if (res_num == FALSE)
                    irange--;
/* <--v.4.0 */
                }

              /* If this is the second of the residue-range parameters,
                 then take it to be the second residue in the range */
              else if (irange == 2)
                {
                  /* Interpret the entered residue number */
                  strcpy(tmpstr1,token[itoken]);
/* v.4.0--> */
/*                get_residue_number(tmpstr1,res_num2); */
                  res_num = get_residue_number(tmpstr1,res_num2,res_name2);

                  /* If not a residue number, then decrement range as
                     still need to pick up the residue number */
                  if (res_num == FALSE)
                    irange--;
/* <--v.4.0 */
                }

              /* The third residue-range parameter will be the chain-id */
              else if (irange == 3)
                *chain_identi = token[itoken][0];
/* v.4.0--> */
            }
/* <--v.4.0 */

          /* Increment count of residue-range parameters */
          irange++;
        }
    }

  /* If only one residue entered for the range, then make the end-residue
     number the same */
  if (res_num2[0] == '\0')
/* v.4.0--> */
    {
/* <--v.4.0 */
      strcpy(res_num2,res_num1);
/* v.4.0--> */
      strcpy(res_name2,res_name1);
    }
/* <--v.4.0 */

  /* If no PDB filename entered, then abort */
  if (pdb_name[0] == '\0')
    {
      printf("*** ERROR. No filename entered!\n");
      exit(-1);
    }

  /* Determine what the names of the corresponding .hhb, .nnb and
     .bonds files are from the entered PDB filename */
  printf("Input file names:-\n");
  printf("   PDB file              %s\n",pdb_name);
  get_filename(pdb_name,hhb_name,".hhb");
/* v.5.4.5 --> */
  if (PDBsum1 == TRUE)
    {
      strcpy(tmp,hhb_name);
      sprintf(hhb_name,"../hbplus/%s",tmp);
    }
/* <-- v.5.4.5 */
  printf("   H-bonds file          %s\n",hhb_name);
  get_filename(pdb_name,nnb_name,".nnb");
/* v.5.4.5 --> */
  if (PDBsum1 == TRUE)
    {
      strcpy(tmp,nnb_name);
      sprintf(nnb_name,"../hbplus/%s",tmp);
    }
/* <-- v.5.4.5 */
  printf("   Non-bonded contacts   %s\n",nnb_name);
  get_filename(pdb_name,bonds_name,".bonds");
  printf("\n");

  /* The .bonds option only currently works where the .bonds file
     has been output as ligplot.bonds */
  if (strncmp(bonds_name,"ligplot.bonds",13))
    bonds_name[0] = '\0';

/* v.4.1.1--> */
  /* Loop over both residue names to right-justify, if necessary */
  strcpy(res_name,res_name1);
  for (loop = 0; loop < 2; loop++)
    {
      /* Get the length of the residue name */
      len = strlen(res_name);

      /* If not blank and not full */
      if (len > 0 && len != 3)
        {
          /* Fill temporary array with blanks */
          strcpy(temp,"   ");
          temp[3] = '\0';

          /* Get the start position for the transfer of characters */
          fpos = 3 - len;

          /* Copy across the required characters */
          for (ipos = 0; ipos < len; ipos++)
            temp[fpos + ipos] = res_name[ipos];

          /* Copy back into appropriate residue name */
          if (loop == 0)
            strcpy(res_name1,temp);
          else
            strcpy(res_name2,temp);
        }

      /* Get the second residue name for the second time around */
      strcpy(res_name,res_name2);
    }
/* <--v.4.1.1 */

  /* Print confirmation of entered details */
/* v.4.0--> */
/*      printf("Residue range:  %s to %s   ",res_num1,res_num2); */
/* v.4.3--> */
/*  if (Interface_Plot == FALSE) */
  if (Interface_Plot == FALSE && Plot_from_drw == FALSE)
/* <--v.4.3 */
    {
/* v.4.2--> */
      /* If plotting all ligands in an XMAS-format file, print info */
      if (Ligand_Molecule_id == 0)
        {
          printf("Plotting all ligands in XMAS-format file\n");
        }

      /* If have a molecule-id for an XMAS-format file, show just the id */
      else if (Ligand_Molecule_id > 0)
        {
          printf("Ligand defined as molecule %d in XMAS-format file\n",
                 Ligand_Molecule_id);
        }

      /* Otherwise, show usual residue range */
      else
        {
/* <--v.4.2 */
          printf("Residue range:  ");
          if (res_name1[0] != '\0')
            printf("[%s] ",res_name1);
          printf("[%s] to ",res_num1);
          if (res_name2[0] != '\0')
            printf("[%s] ",res_name2);
          printf("[%s]   ",res_num2);
/* <--v.4.0 */
          if (*chain_identi != '#')
            printf("Chain [%c]\n",*chain_identi);
          else
            printf("\n");
/* v.4.2--> */
        }
/* <--v.4.2 */

/* v.4.0--> */
    }
/* <--v.4.0 */

/* v.5.4.5 --> */
  /* For PDBsum1, form the plot title */
  if (PDBsum1 == TRUE)
    {
      /* Get the plot type */
      if (Metal_as_Ligand == TRUE)
        strcpy(ligmet,"Metal ion:");
      else
        strcpy(ligmet,"Ligand:");

      if (res_name2[0] != '\0')
        sprintf(Drw_Title,"%s %s %s(%c) - %s %s(%c)",ligmet,
                res_name1,res_num1,*chain_identi,
                res_name2,res_num2,*chain_identi);
      else
        sprintf(Drw_Title,"%s %s %s(%c)",ligmet,
                res_name1,res_num1,*chain_identi);
    }
/* <-- v.5.4.5 */

  /* If the input file is called 'ligplot.pdb' then want
     to leave the structure as it is */
  if (!strncmp(pdb_name,"ligplot.pdb",11))
      Print_as_is = TRUE;

/* v.4.3--> */
  /* Similarly, if the file is called 'ligplot.drw' then plot as is */
  if (!strcmp(pdb_name,"ligplot.drw"))
    {
      Print_as_is = TRUE;
      Plot_from_drw = TRUE;
    }
/* <--v.4.3 */

/* v.4.2--> */
  /* If printing from a ligplot.pdb file, then it won't be in XMAS 
     format */
  if (Print_as_is == TRUE)
    {
      Read_Xmas = FALSE;
      Ligand_Molecule_id = -1;
      Plot_All_Xmas_Ligands = FALSE;
    }
/* <--v.4.2 */

/* v.5.0.3--> */
//   /* If no title entered, and option to use filename as the default
//      title, then store the filename as the title */
////   if (Include->Filename_for_Title == TRUE && Print_Title[0] == '\0' &&
///* v.4.0--> */
///*       Print_as_is == FALSE) */
//       Print_as_is == FALSE && Interface_Plot == FALSE)
//* <--v.4.0 */
//        strcpy(Print_Title,root_name);
//
//  /* If the ligand is a water molecule, then make sure that water
//     molecules are allowed */
//  if (Water_as_Ligand == TRUE)
//    Include->Waters = TRUE;
/* <--v.5.0.3 */
}
/***********************************************************************

read_in_parameters  -  Read in the plot parameters from the
                       ligplot.par parameter file

***********************************************************************/

void read_in_parameters(char special_res[MAXSPECIAL_RES][15])
{
  char line[LINELEN + 1];
  char colour_name[COL_NAME_LEN + 1], residue_pair[8];
  char show_string[14];
/* v.5.4.6--> */
  char ligplot_version[] = "LIGPLOT v.5.4.6";
/* <--v.5.4.6 */
  char old_version[16];
  char parameter_type, parameter_types[] = {'P', 'Y', 'L', 'S', 'S',
                                              'C', 'C', 'D', 'S'};
  int iline, last_group, parameter, parameter_group, parameter_lines;
  int number_of_parameters[] = { PICTURE_OPTIONS, INCLUDE_PARAMETERS,
                                   MAXSPECIAL_RES, OBJECT_SIZES,
                                   TEXT_SIZES, OBJECT_COLOURS,
                                   TEXT_COLOURS, MAX_COLOURS,
                                   MINIMIZATION_PARAMETERS };
  int true_false;
  float blue, green, red, size;

/* v.4.0--> */
  /* Print version number */
  printf("\n");
  printf("Running %s ...\n",ligplot_version);
  printf("\n");
/* <--v.4.0 */

  /* Open the parameter file */
/* v.5.0--> */
  /* if ((ligplot_par = fopen("ligplot.prm","r")) == NULL) */
  if ((ligplot_par = fopen(ligplot_prm,"r")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open parameter file: ligplot.prm\n\n"); */
      printf("\n*** Unable to open parameter file: %s\n\n",ligplot_prm);
/* <--v.5.0 */
      exit(1);
    }

  /* Initialise flags and counts */
  strcpy(colour_name,"WHITE");
  iline = 0;
  last_group = -1;
  parameter = 0;
  parameter_group = -1;
  parameter_lines = -1;
  size = 0.0;
  true_false = FALSE;

  /* Loop while reading records from the parameter file */
  while (fgets(line,LINELEN,ligplot_par)!= NULL)
    {
      /* If this is the first line, then check that the parameter
         file corresponds to the current program version number */
      if (iline == 0)
        {
          if (strncmp(line,ligplot_version,13))
            {
              /* Check for a version which may work OK */
/* v.4.0--> */
/*            if (!strncmp(line,ligplot_version,11)) */
              if (!strncmp(line,"LIGPLOT v.3",11))
/* <--v.4.0 */
                {
                  printf("*** Warning. Your ligplot.prm file is for ");
                  strncpy(old_version,line,15);
                  old_version[15] = '\0';
                  printf("%s\n",old_version);
                  printf("***          You will get better results if you");
                  printf(" delete the file and re-run the\n");
/* v.4.0--> */
/*                printf("***          program which is %s\n",
                         ligplot_version); */
                  printf("***          program to generate a new ");
                  printf("parameter file for %s\n",ligplot_version);
                  printf("\n");
/* <--v.4.0 */
                  Nwarnings++;
                }
/* v.4.1--> */
              /* Versions 4.0 and 4.1 parameter files will work OK */
/* v.4.1.1--> */
/*            else if (!strncmp(line,"LIGPLOT v.4.0",13)) */
              else if (!strncmp(line,"LIGPLOT v.4.0",13) ||
/* v.4.2--> */
/*                     !strncmp(line,"LIGPLOT v.4.1",14)) */
                       !strncmp(line,"LIGPLOT v.4.1",13))
/* <--v.4.2 */
/* <--v.4.1.1 */
                {
                  printf("\n");
                }
/* <--v.4.1 */

              /* If this is a ligplot file for any other version, then
                 unlikely to work OK */
              else if (!strncmp(line,"LIGPLOT",7))
                {
                  printf("*** ERROR. Incorrect version of parameter file");
                  printf(", ligplot.prm\n");
                  strncpy(show_string,line,13);
                  show_string[13] = '\0';
                  printf("           File is for %s, whereas program ",
                         show_string);
                  printf("requires %s\n",ligplot_version);
                  exit(1);
                }

              /* Otherwise, if not a LIGPLOT parameter file at all, then
                 forget it */
              else
                {
                  printf("*** ERROR. Unrecognised first line in parameter");
                  printf(" file, ligplot.prm\n");
                  printf("%s\n",line);
                  exit(1);
                }
            }
        }

      /* Check for one of the keyword lines */
      if (!strncmp(line,"PRINT OPTIONS",13))
        parameter_group = 0;
      else if (!strncmp(line,"PLOT PARAMETERS",15))
        parameter_group = 1;
      else if (!strncmp(line,"LINKED RESIDUES",15))
        parameter_group = 2;
      else if (!strncmp(line,"SIZES",5))
        parameter_group = 3;
      else if (!strncmp(line,"TEXT SIZES",10))
        parameter_group = 4;
      else if (!strncmp(line,"COLOURS",7))
        parameter_group = 5;
      else if (!strncmp(line,"TEXT COLOURS",12))
        parameter_group = 6;
      else if (!strncmp(line,"COLOUR DEFINITIONS",18))
        parameter_group = 7;
      else if (!strncmp(line,"MINIMIZATION PARAMETERS",23))
        parameter_group = 8;

      /* If this is a new parameter-group, then initialise
         line-counts */
      if (parameter_group != last_group)
        {
          parameter = -2;
          last_group = parameter_group;
          parameter_lines = number_of_parameters[parameter_group];
          parameter_type = parameter_types[parameter_group];
        }

      /* If this line is one of the current group of parameters,
         then extract the values and store */
      if (parameter > -1 && parameter < parameter_lines)
        {
          /* Check for a Yes/No response */
          if (line[0] == 'Y' || line[0] == 'y')
            true_false = TRUE;
          else
            true_false = FALSE;

          /* If this is a yes/no response, then get its value */
          if (parameter_type == 'Y')
            {
              if (line[0] == 'Y' || line[0] == 'y')
                true_false = TRUE;
              else
                true_false = FALSE;
              if (line[0] != 'Y' && line[0] != 'y' &&
                  line[0] != 'N' && line[0] != 'n')
                sscanf(line,"%f ",&size);
            }

          /* If this is one of the print options, interpret it */
          else if (parameter_type == 'P')
            {
              if (parameter == 0 && (line[0] == 'Y' || line[0] == 'y'))
                true_false = TRUE;
              else if (parameter == 1 && (line[0] == 'P' || line[0] == 'p'))
                true_false = TRUE;
              else
                true_false = FALSE;
              if (parameter == 2)
                sscanf(line,"%f ",&size);
            }

          /* If this is a colour, then store its name */
          else if (parameter_type == 'C')
            {
/* v.4.1--> */
/*            strncpy(colour_name,line,COL_NAME_LEN);
              colour_name[COL_NAME_LEN] = '\0'; */
              strncpy(colour_name,line,PRM_COL_NAME_LEN);
              colour_name[PRM_COL_NAME_LEN] = '\0';
/* <--v.4.1 */
            }

          /* If this is a size, then retrieve the value */
          else if (parameter_type == 'S')
            {
              /* Read in the value of the size */
              sscanf(line,"%f ",&size);
            }

          /* If this is a link-residue pair, retrieve the residue names */
          else if (parameter_type == 'L')
            {
              strncpy(residue_pair,line,7);
              residue_pair[7] = '\0';
            }

          /* If this is a colour definition line, read in the RGB
             values and the colour name */
          else if (parameter_type == 'D')
            {
              /* Read in the RGB values */
              sscanf(line,"%f %f %f ",&red,&green,&blue);
              
              /* Read in the colour name */
/* v.4.1--> */
/*            strncpy(colour_name,line+22,COL_NAME_LEN);
              colour_name[COL_NAME_LEN] = '\0'; */
              strncpy(colour_name,line+22,PRM_COL_NAME_LEN);
              colour_name[PRM_COL_NAME_LEN] = '\0';
/* <--v.4.1 */
            }

          /* Update the current parameter */
          update_current_parameter(parameter_group,parameter,true_false,
                                   size,colour_name,red,green,blue,
                                   residue_pair,special_res);
        }
      /* Increment the parameter line-count */
      parameter++;
      iline++;
    }

  /* Correct any incompatible parameters */
/* v.4.0--> */
/*  if (Include->Simple_Ligand_Residues == TRUE ||
      Include->Simple_Nonligand_Residues == TRUE)
    {
      Include->Accessibilities = FALSE;
      Include->External_Bonds = FALSE;
    } */
  if (Include->Simple_Ligand_Residues == TRUE)
    {
      Include->Accessibilities = FALSE;
      Include->External_Bonds = FALSE;
    }
  if (Include->Simple_Nonligand_Residues == TRUE)
    {
      Include->External_Bonds = FALSE;
      Include->Ligand_Accessibilities_Only = TRUE;
    }
/* <--v.4.0 */

  /* For interface plots, need to ensure that certain options aren't set */
  if (Interface_Plot == TRUE)
    {
      Include->Simple_Ligand_Residues = FALSE;
      Include->Simple_Nonligand_Residues = FALSE;

      /* Initialise weight for boundary energy calcs */
      if (Rel_Pstn_Energy_Weight < 0.0)
        Rel_Pstn_Energy_Weight = (float) (0.3);
    }

/* v.4.3--> */
  /* For chemical notation, switch off "sphere" representations and
     printing of atom names */
  if (Include->Chemical_Notation == TRUE)
    {
      /* Switch off "sphere" representations */
/* v.4.4--> */
/*      Include->Ligand_Atoms = FALSE;
      Include->Nonligand_Atoms = FALSE;
      Include->Water_Atoms = FALSE; */
/* <--v.4.4 */

      /* Switch on double-bond representation */
      Include->Double_Bonds = TRUE;

      /* Switch off atom names */
      Include->Atom_Names = FALSE;
    }
/* <--v.4.3 */
/* v.5.0--> */
  if (ligfit == TRUE)
    {
      Anchor_Weight = 1;
/* v.5.2.1 --> */
//      Max_Loops = 2000;
      Max_Loops = 3000;
/* <-- v.5.2.1 */
    }
/* <--v.5.0 */

/* v.5.3 --> */
  /* If contact type defined on command line, override the value in the
     parameter file */
  if (contact_type !=  NOT_DEFINED)
    Include->Contact_Type = contact_type;
/* <--v.5.3 */

/* v.4.2--> */
  /* Close the parameter file */
  fclose(ligplot_par);
/* <--v.4.2 */
}
/* v.4.2--> */
/***********************************************************************

initialise_variables  -  Initialise the variables for the current plot

***********************************************************************/

void initialise_variables(void)
{
  /* Re-initialise page variables */
  Page_Min_x = BBOXX1;
  Page_Max_x = BBOXX2;
  Page_Min_y = BBOXY1;
  Page_Max_y = BBOXY2;
  Margin_y = (float)Page_Min_y;

  /* Re-initialise pointers */
  first_bond_ptr = NULL;
  first_hhb_info_ptr = NULL;
  last_hhb_info_ptr = NULL;

  /* Special flags and values */
  Max_atom_radius = 0.0;
  Maximum_accessibility = (float) (100.0);
  Nobjects = 0;
  Plot_Centre_x = (float) ((Page_Min_x + Page_Max_x) / 2.0);
  Plot_Centre_y = (float) ((Page_Min_y + Page_Max_y) / 2.0);
}
/* <--v.4.2 */
/***********************************************************************

check_if_in_ligand  -  Determine whether the given residue is a ligand
                       residue

***********************************************************************/

int check_if_in_ligand(char chain, char res_name[4], char res_num[6],
/* v.4.2.2--> */
/*                     int partial_match,char duplicate_res_name[4]) */
                       int partial_match,char duplicate_res_name[4],
                       int mol_id)
/* <--v.4.2.2 */
{
  int ires, inligand;

  /* Initialise answer */
  inligand = FALSE;
  ires = 0;
  duplicate_res_name[0] = '\0';

/* v.4.2.2--> */
  /* If reading from a XMAS file, need only to verify the molecule-id */
  if (partial_match == FALSE && Ligand_Molecule_id > -1 && mol_id > -1)
    {
      /* Check the molecule identifier */
      if (mol_id == Ligand_Molecule_id)
        inligand = TRUE;
    }

  /* Otherwise, use ligand range */
  else
    {
/* <--v.4.2.2 */

      /* Loop through the stored names of the residues making up the ligand
         to see if given residue is amongst them */
      while (ires < ligand_residues && inligand == FALSE)
        {
          /* First check if have match on residue number and chain-id */
          if (chain == lig_chain_store[ires] &&
              !strncmp(res_name,lig_res_name_store[ires],3) &&
              !strncmp(res_num,lig_res_num_store[ires],5))
            inligand  = TRUE;

          /* If only a partial match required, check for match on residue
             number and chain-id only */
          if (partial_match == TRUE)
            if (chain == lig_chain_store[ires] &&
                !strncmp(res_num,lig_res_num_store[ires],5))
              {
                inligand  = TRUE;
                strncpy(duplicate_res_name,lig_res_name_store[ires],3);
                duplicate_res_name[3] = '\0';
              }

          /* Increment ligand count */
          ires++;
        }
/* v.4.2.2--> */
    }
/* <--v.4.2.2 */

  /* Return answer */
  return(inligand);
}
/***********************************************************************

get_ligand_start_end  -  Read through PDB file to get the start
                         and end residue positions of the ligand

***********************************************************************/

/* v.4.0--> */
/* void get_ligand_start_end(char pdb_name[FILENAME_LEN],char res_num1[6],
                          char res_num2[6],
                          char chain_id) */
/* v.4.2--> */
/* void get_ligand_start_end(char pdb_name[FILENAME_LEN],char res_name1[4],
                          char res_num1[6],char res_name2[4],char res_num2[6],
                          char chain_id) */
void get_ligand_start_end(char pdb_name[FILENAME_LEN],XMAS *xmas,
                          char res_name1[4],char res_num1[6],
                          char res_name2[4],char res_num2[6],
                          char chain_id)
/* <--v.4.2 */
/* <--v.4.0 */
{
  char line[LINELEN + 1];
  char duplicate_resnam[4], resnam[4];
  char chain, last_chain, last_resnum[6], resnum[6];
  char last_resnam[4];
  char icode, icode_from, icode_to, temp_num[5];
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */

  int resno, resno_from, resno_to;
  int inligand, residue_count;
  int duplicate;
  int fstpos, ipos, ires;
  int keep_reading;
/* v.4.0--> */
  int first_residue, got_end, got_ligand, last_residue;
/* <--v.4.0 */
/* v.4.0.2--> */
/* v.4.4--> */
/*  int backwards; */
  int backwards, exact;
/* <--v.4.4 */
/* <--v.4.0.2 */
/* v.4.2--> */
  int rec_types[NTYPES], type;
/* <--v.4.2 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */

  /* Open the pdb file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
    {
/* <--v.4.2 */
      if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
        {
          printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
          exit(1);
        }
/* v.4.2--> */
    }
/* <--v.4.2 */

  /* Initialise variables */
  keep_reading = TRUE;
/* v.4.0--> */
  got_end = FALSE;
  got_ligand = FALSE;
/* <--v.4.0 */
  inligand = FALSE;
  ligand_residues = 0;
  residue_count = 0;

  /* Read in the data from the PDB file */
  printf("\nSearching for ligand residues in PDB file ...\n");

  /* Initialise variables */
  last_resnam[0] = '\0';
  last_resnum[0] = '\0';
  last_chain = '\0';

  /* Split the residue numbers giving the range of residues making up the
     ligand into their residue number and insertion code components */

  /* First residue of range */
  strncpy(temp_num,res_num1,4);
  temp_num[4] = '\0';
  resno_from = atoi(temp_num);
  icode_from = res_num1[4];
    
  /* Second residue of range */
  strncpy(temp_num,res_num2,4);
  temp_num[4] = '\0';
  resno_to = atoi(temp_num);
  icode_to = res_num2[4];
    
/* v.4.0.2--> */
  /* Check whether residues are numbered backwards (!) */
  backwards = FALSE;
  if (resno_from > resno_to ||
      (resno_from ==resno_to && icode_from > icode_to))
    backwards = TRUE;
/* <--v.4.0.2 */

/* v.4.4--> */
  /* Check whether ligand range has been exactly specified (ie with residue
     names as well as numbers) */
  exact = FALSE;
  if (res_name1[0] != '\0' && res_name2[0] != '\0')
    exact = TRUE;
/* <--v.4.4 */

/* v.4.2--> */
  /* Get the first line from the PDB or XMAS file */
  if (Read_Xmas == FALSE)
    eofstatus = fgets(line,LINELEN,fil_pdb);
  else
    {
      /* Form the record types to be retrieved */
      rec_types[REC_ATOM] = TRUE;
      rec_types[REC_CONECT] = FALSE;
      type = REC_ATOM;

      /* Get first record from the XMAS cache */
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type); */
      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
    }
/* <--v.4.2 */

/* v.4.2--> */
/*  while (fgets(line,LINELEN,fil_pdb) != NULL && keep_reading == TRUE) */
  /* Search through the PDB file to find residue-range defining the
     ligand */
  while (eofstatus != NULL && keep_reading == TRUE)
/* <--v.4.2 */
    {
      /* If this is an ATOM or HETATM record, process it */
      if (!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6))
        {
          /* Get the residue name, number and chain-ID */
          strncpy(resnam,line+17,3);
          resnam[3] = '\0';
          strncpy(resnum,line+22,5);
          resnum[5] = '\0';
          chain = line[21];

          /* Check whether this is a new residue */
          if (chain != last_chain || strncmp(resnum,last_resnum,5) ||
              strncmp(resnam,last_resnam,3))
            {
              /* Split the residue number in an integer and insertion code */
              strncpy(temp_num,resnum,4);
              temp_num[4] = '\0';
              resno = atoi(temp_num);
              icode = resnum[4];

              /* Determine whether this residue falls in the range defining
                 the ligand */
              inligand = FALSE;
/* v.4.0--> */
              first_residue = FALSE;
              last_residue = FALSE;
/* <--v.4.0 */

/* v.4.4--> */
              /* If residue range has been exactly specified, check whether
                 this is the start residue, the end residue, or somewhere
                 in between */
              if (exact == TRUE)
                {
                  /* Check for a match with the first residue */
                  if (!strcmp(resnam,res_name1) && !strcmp(resnum,res_num1) &&
                      chain == chain_id)
                    {
                      /* Set flag that this residue is in the ligand and
                         that we've hit the ligand-start */
                      inligand = TRUE;
                      first_residue = TRUE;
                      got_ligand = TRUE;
                    }

                  /* Check for a match with the final residue */
                  if (!strcmp(resnam,res_name2) && !strcmp(resnum,res_num2) &&
                      chain == chain_id)
                    {
                      /* Set flag that ligand-end reached */
                      if (got_ligand == TRUE)
                        {
                          inligand = TRUE;
                          last_residue = TRUE;
                          got_end = TRUE;
                        }
                    }                 

                  /* Otherwise, if have had a ligand-start, without a
                     ligand end, then take this residue to be in the
                     ligand */
                  else if (got_ligand == TRUE && got_end == FALSE)
                    {
                      /* Check that not a water */
                      if (!strncmp(resnam,"HOH",3) &&
                          Water_as_Ligand == FALSE)
                        inligand = FALSE;
                      else
                        inligand = TRUE;
                    }
                }

              /* Otherwise, perform all the old checks */
              else
                {
/* <--v.4.4 */
              
/* v.4.0.2--> */
/*            if (chain == chain_id && resno >= resno_from &&
                  resno <= resno_to) */
                  if (chain == chain_id &&
                      ((backwards == FALSE && resno >= resno_from &&
                        resno <= resno_to) ||
                       (backwards == TRUE && resno <= resno_from &&
                        resno >= resno_to)))
/* <--v.4.0.2 */
                    {
                      /* Set flag to say residue is in ligand */
                      inligand = TRUE;

                      /* If residue number is start of range, check the
                         insertion code is within the range, too */
/* v.4.0.2--> */
/*                if (resno == resno_from && icode < icode_from) */
                      if (resno == resno_from &&
                          ((backwards == FALSE && icode < icode_from) ||
                           (backwards == TRUE && icode > icode_from)))
/* <--v.4.0.2 */
                        inligand = FALSE;

/* v.4.0--> */
                      /* Check whether this is the start of the residue
                         range */
                      if (resno == resno_from && icode == icode_from)
                        {
                          first_residue = TRUE;
                          got_ligand = TRUE;
                        }

                      /* Check whether this is the end of the residue range */
                      if (resno == resno_to && icode == icode_to)
                        last_residue = TRUE;
/* <--v.4.0 */

                      /* If residue number is end of range, check the
                         insertion code is within the range, too */
/* v.4.0.2--> */
/*                if (resno == resno_to && icode > icode_to) */
                      if (resno == resno_to &&
                          ((backwards == FALSE && icode > icode_to) ||
                           (backwards == TRUE && icode < icode_to)))
/* <--v.4.0.2 */
                        inligand = FALSE;
                    }
/* v.4.0--> */
                  else
                    got_ligand = FALSE;
/* <--v.4.0 */

                  /* If this is a water molecule, then cannot be in the
                     ligand unless -w option selected in command-line
                     parameters */
                  if (!strncmp(resnam,"HOH",3) && Water_as_Ligand == FALSE)
                    inligand = FALSE;

/* v.4.0--> */
                  /* If this is the start-residue and have a name for the
                     start-residue, then check whether this name corresponds */
                  if (first_residue == TRUE && res_name1[0] != '\0')
                    {
                      /* If the residue name doesn't match, then it is not the
                         first ligand residue */
                      if (strncmp(resnam,res_name1,3))
                        got_ligand = FALSE;
                    }

                  /* If haven't encountered a proper start for ligand, or have
                     hit the recognised end, then even if residue is in the
                     right range, then don't want it */
                  if (got_ligand == FALSE || got_end == TRUE)
                    inligand = FALSE;

                  /* Check whether this is the predefined residue end */
                  if (last_residue == TRUE && res_name2[0] != '\0')
                    {
                      /* If the residue name matches, then have the last
                         ligand residue */
                      if (!strncmp(resnam,res_name2,3))
                        got_end = TRUE;
                    }
/* <--v.4.0 */

/* v.4.4--> */
                }
/* <--v.4.4 */

              /* If residue is in the ligand, then process */
              if (inligand == TRUE)
                {
                  /* Check whether a residue with the same residue
                     number has already been stored */
                  duplicate = check_if_in_ligand(chain,resnam,resnum,
/* v.4.2.2--> */
/*                                               TRUE,duplicate_resnam); */
                                                 TRUE,duplicate_resnam,
                                                 mol_id);
/* <--v.4.2.2 */

                  /* If have a duplicate, then print warning message */
                  if (duplicate == TRUE)
                    {
                      printf("\n");
                      printf("*** Warning. Duplicate res. no. in");
                      printf(" PDB file for ligand res: %s %c (%s & %s)\n",
                             resnum,chain,resnam,duplicate_resnam);
                      printf("\n");
                      inligand = TRUE;
                      Nwarnings++;
                    }

                  /* Check that maximum number of ligand residues
                     not exceeded */
                  if (ligand_residues > MAXLIGRES - 1)
                    {
/* v.4.4--> */
/*                      printf("*** Warning. Maximum number of ligand ");
                      printf("residues (%d) exceeded!\n",MAXLIGRES);
                      Nwarnings++; */
                      printf("*** ERROR. Maximum number of ligand ");
                      printf("residues (%d) exceeded!\n",MAXLIGRES);
/* v.5.4.4 --> */
                      /* Write message out to the ligplot.drw file */
                      fprintf(ligplot_drw,
                              "ERROR: Maximum number of ligand "
                              "residues (%d) exceeded\n",MAXLIGRES);

                      /* Close the .drw file */
                      fclose(ligplot_drw);
/* <-- v.5.4.4 */
                      exit(-1);
/* <--v.4.4 */
                    }

                  /* Store the residue details */
                  else
                    {
                      strncpy(lig_res_name_store[ligand_residues],resnam,3);
                      lig_res_name_store[ligand_residues][3] = '\0';
                      strncpy(lig_res_num_store[ligand_residues],resnum,5);
                      lig_res_num_store[ligand_residues][5] = '\0';
                      lig_chain_store[ligand_residues] = chain;
                      ligand_residues++;
                    }
                }

              /* Save the residue details */
              strncpy(last_resnam,resnam,3);
              last_resnam[3] = '\0';
              strncpy(last_resnum,resnum,5);
              last_resnum[5] = '\0';
              last_chain = chain;

              /* Increment the residue-count */
              residue_count++;
            }
        }

      /* If this is an ENDMDL record, stop reading here */
      else if (!strncmp(line,"ENDMDL",6))
        keep_reading = FALSE;

/* v.4.2--> */
      /* Get the next atom record */
      if (keep_reading == TRUE)
        {
          /* Get the next line from the PDB or XMAS file */
          if (Read_Xmas == FALSE)
            eofstatus = fgets(line,LINELEN,fil_pdb);
          else
/* v.4.2.2--> */
/*          eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type); */
            eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type,
                                          &mol_id);
        }
/* <--v.4.2.2 */
/* <--v.4.2 */
    }

  /* Print counts of data read in */
  printf("   Number of residues read in        = %7d\n",residue_count);
  printf("   Number of ligand residues         = %7d\n",ligand_residues);

  if (ligand_residues == 0)
    {
      printf("\n*** No ligand residues found. Nothing to plot\n");

/* v.5.4.4 --> */
      /* Write message out to the ligplot.drw file */
      fprintf(ligplot_drw,
              "ERROR: No ligand residues found\n");

      /* Close the .drw file */
      fclose(ligplot_drw);
/* <-- v.5.4.4 */
      exit(1);
    }
  else
    {
      /* Print the ligand residues identified */
      printf("\n");
      printf("   Ligand residues: ");
      ipos = 20;
      for (ires = 0; ires < ligand_residues; ires++)
        {
          fstpos = 0;
          if (lig_res_num_store[ires][0] == ' ')
            fstpos = 1;
          if (lig_res_num_store[ires][1] == ' ')
            fstpos = 2;
          if (lig_res_num_store[ires][2] == ' ')
            fstpos = 3;
          printf("%s %s",lig_res_name_store[ires],
                 lig_res_num_store[ires]+fstpos);
          ipos = ipos + 9 - fstpos;
          if (lig_chain_store[ires] != ' ')
            {
              printf("(%c)",lig_chain_store[ires]);
              ipos = ipos + 3;
            }
          if (ires < ligand_residues - 1)
            {
              printf(" - ");
              ipos = ipos + 3;
            }
          if (ipos > 60)
            {
              printf("\n                    ");
              ipos = 20;
            }
        }
      printf("\n\n");
    }
  
  /* Close the PDB file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
/* <--v.4.2 */
    fclose(fil_pdb);
}
/* v.5.0--> */
/***********************************************************************

create_hhb_info_record  -  Create a new hhb_info record

***********************************************************************/

struct hhb_info *create_hhb_info_record(int atom_number1,char *atom_name1,
                                        char *res_name1,char *res_num1,
                                        char chain1,
                                        int atom_number2,char *atom_name2,
                                        char *res_name2,char *res_num2,
                                        char chain2,int source,
                                        float bond_length)
{
  struct hhb_info *hhb_info_ptr;

  /* Initialise hhb_info pointers */
  hhb_info_ptr = NULL;

  /* Create hhb_info entry */
  hhb_info_ptr = (struct hhb_info*)malloc(sizeof(struct hhb_info));
  if (hhb_info_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for");
      printf(" struct hhb_info\n");
      printf("***        Program ligplot terminated with error.\n");
      exit (1);
    }

  /* Store the details for this entry */
  hhb_info_ptr->source = source;
  hhb_info_ptr->atom_number1 = atom_number1;
  strcpy(hhb_info_ptr->atom_name1,atom_name1);
  strcpy(hhb_info_ptr->res_name1,res_name1);
  strcpy(hhb_info_ptr->res_num1,res_num1);
  hhb_info_ptr->chain1 = chain1;
  hhb_info_ptr->atom_number2 = atom_number2;
  strcpy(hhb_info_ptr->atom_name2,atom_name2);
  strcpy(hhb_info_ptr->res_name2,res_name2);
  strcpy(hhb_info_ptr->res_num2,res_num2);
  hhb_info_ptr->chain2 = chain2;
  hhb_info_ptr->bond_length = bond_length;

  /* Initialise the remaining hhb_info details */
/* v.5.4--> */
  hhb_info_ptr->bond_type = INVALID;
/* <--v.5.4 */
  hhb_info_ptr->x1 = 0.0;
  hhb_info_ptr->y1 = 0.0;
  hhb_info_ptr->z1 = 0.0;
  hhb_info_ptr->x2 = 0.0;
  hhb_info_ptr->y2 = 0.0;
  hhb_info_ptr->z2 = 0.0;

  /* If this is the first entry stored, then update
     pointer to head of linked list */
  if (first_hhb_info_ptr == NULL)
    first_hhb_info_ptr = hhb_info_ptr;

  /* Otherwise, make previous entry point to the current
     one */
  else
    last_hhb_info_ptr->next_hhb_info_ptr = hhb_info_ptr;
                      
  /* Save the current hhb_info's pointer */
  last_hhb_info_ptr = hhb_info_ptr;

  /* Return current pointer */
  return(hhb_info_ptr);
}
/* <--v.5.0 */
/***********************************************************************

read_conec_records  -  Read through CONECT records to pick up any
                       non-ligand residues that are covalently
                       bonded to the ligand, and to store these.

***********************************************************************/

/* v.4.2--> */
/* void read_conec_records(char pdb_name[FILENAME_LEN],int *nlinks) */
void read_conec_records(char pdb_name[FILENAME_LEN],int *nlinks,XMAS *xmas)
/* <--v.4.2 */
{
/* v.4.2--> */
/*  char line[101]; */
  char line[LINELEN + 1];
/* <--v.4.2 */
  char atmnum[6];
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */

  int done, have_link;
  int iatom, icon, ipos, jatom, ncon;
  int store_atom[8];
  int atom_number, last_atom_number;
  int keep_reading;
/* v.4.2--> */
  int rec_types[NTYPES], type;
/* <--v.4.2 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */

  struct hhb_info *hhb_info_ptr;

  /* Open the pdb file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
    {
/* <--v.4.2 */
      if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
        {
          printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
          exit(1);
        }
/* v.4.2--> */
    }
/* <--v.4.2 */
    
  /* Initialise variables */
  *nlinks = 0;
  last_atom_number = 0;
  keep_reading = TRUE;

  /* Read in the data from the PDB file */
  printf("\nChecking CONECT records in the PDB file ...\n");

/* v.4.2--> */
  /* Get the first line from the PDB or XMAS file */
  if (Read_Xmas == FALSE)
    eofstatus = fgets(line,LINELEN,fil_pdb);
  else
    {
      /* Form the record types to be retrieved */
      rec_types[REC_ATOM] = FALSE;
      rec_types[REC_CONECT] = TRUE;
      type = REC_CONECT;
      last_atom_number = 99999;

      /* Get first record from the XMAS cache */
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type); */
      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
    }
/* <--v.4.2 */

  /* Search through the PDB file to scan the CONECT records */
/* v.4.2--> */
/*  while (fgets(line,100,fil_pdb) != NULL) */
  while (eofstatus != NULL)
/* <--v.4.2 */
    {
      /* If this is an ATOM or HETATM record, store the atom number */
      if ((!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)) &&
               keep_reading == TRUE)
        {
          /* Store atom-number of last atom read in (for use in whittling
             out which CONECT records are to be processed) */
          strncpy(atmnum,line+6,5);
          atmnum[5] = '\0';
          last_atom_number = atoi(atmnum);
        }

      /* If this is an ENDMDL record, then not interested in any
         of the following ATOM and HETATM records */
      else if (!strncmp("ENDMDL",line,6))
        keep_reading = FALSE;

      /* If this is a CONECT record, process it */
      else if (!strncmp(line,"CONECT",6))
        {
          /* Get the atom-number of the first number in the CONECT record */
          strncpy(atmnum,line+6,5);
          atmnum[5] = '\0';
          atom_number = atoi(atmnum);

          /* If the atom number is within the range of atoms read in,
             then process */
          if (atom_number <= last_atom_number)
            {
              /* Extract all the atom numbers in the CONECT record */
              ipos = 6;
              ncon = 0;
              done = FALSE;
              while (done == FALSE && ncon < 5)
                {
                  /* Extract this atom number and store */
                  strncpy(atmnum,line+ipos,5);
                  if (line[ipos] == '\0' || line[ipos] == '\n')
                    strncpy(atmnum,"     ",5);
                  atmnum[5] = '\0';

                  /* If have blank, then end of connections reached */
                  if (!strncmp(atmnum,"     ",5))
                    done = TRUE;

                  /* Otherwise, store the number */
                  else
                    {
                      store_atom[ncon] = atoi(atmnum);

                      /* Increment variables */
                      ncon++;
                      ipos = ipos + 5;
                    }
                }

              /* Consider each connection in turn */
              iatom = store_atom[0];
              for (icon = 1; icon < ncon; icon++)
                {
                  jatom = store_atom[icon];

                  /* Check whether we already have this link */
                  have_link = FALSE;

                  /* Set pointer to the first of the links */
                  hhb_info_ptr = first_hhb_info_ptr;

                  /* Search through all the stored links */
                  while (hhb_info_ptr != NULL)
                    {
                      /* Check if we have this link already */
                      if ((hhb_info_ptr->atom_number1 == iatom &&
                           hhb_info_ptr->atom_number2 == jatom) ||
                          (hhb_info_ptr->atom_number1 == jatom &&
                           hhb_info_ptr->atom_number2 == iatom))
                        have_link = TRUE;

                      /* Get the pointer to the next link */
                      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
                    }
                      
                  /* Add if this is a new link */
                  if (have_link == FALSE)
                    {
                      /* Allocate memory for structure to hold
                         current interaction */
                      hhb_info_ptr = (struct hhb_info *)
                        malloc(sizeof(struct hhb_info));

                      /* If this is the first interaction stored, then
                         update pointer to head of linked list */
                      if (first_hhb_info_ptr == NULL)
                        first_hhb_info_ptr = hhb_info_ptr;
              
                      /* Otherwise, make previous interaction point to
                         the current one */
                      else
                        last_hhb_info_ptr->next_hhb_info_ptr
                          = hhb_info_ptr;
              
                      /* Save the current interaction's pointer */
                      last_hhb_info_ptr = hhb_info_ptr;

                      /* If either atom is in the ligand store this
                         connection as a proper CONECT */

                      /* Otherwise, store it as an EXTRA, possibly
                         for later use */
                      hhb_info_ptr->source = CONECT;
/* v.5.4--> */
                      hhb_info_ptr->bond_type = S_COVALENT;
/* <--v.5.4 */

                      /* Store the data */
                      hhb_info_ptr->atom_number1 = iatom;
                      hhb_info_ptr->atom_number2 = jatom;

                      /* Initialise the other fields in the
                         structure */

                      /* First atom's details */
                      strncpy(hhb_info_ptr->atom_name1,"    ",4);
                      hhb_info_ptr->atom_name1[4] = '\0';
                      strncpy(hhb_info_ptr->res_name1,"   ",3);
                      hhb_info_ptr->res_name1[3] = '\0';
                      strncpy(hhb_info_ptr->res_num1,"     ",5);
                      hhb_info_ptr->res_num1[5] = '\0';
                      hhb_info_ptr->x1 = 0.0;
                      hhb_info_ptr->y1 = 0.0;
                      hhb_info_ptr->z1 = 0.0;

                      /* Second atom's details */
                      strncpy(hhb_info_ptr->atom_name2,"    ",4);
                      hhb_info_ptr->atom_name2[4] = '\0';
                      strncpy(hhb_info_ptr->res_name2,"   ",3);
                      hhb_info_ptr->res_name2[3] = '\0';
                      strncpy(hhb_info_ptr->res_num2,"     ",5);
                      hhb_info_ptr->res_num2[5] = '\0';
                      hhb_info_ptr->chain1 = ' ';
                      hhb_info_ptr->chain2 = ' ';
                      hhb_info_ptr->x2 = 0.0;
                      hhb_info_ptr->y2 = 0.0;
                      hhb_info_ptr->z2 = 0.0;

                      /* Additional data fields */
                      hhb_info_ptr->bond_length = 0.0;
                      hhb_info_ptr->next_hhb_info_ptr = NULL;

                      /* Increment number of connections stored */
                      (*nlinks)++;
                    }
                }
            }
        }

/* v.4.2--> */
      /* Get the next line from the PDB or XMAS file */
      if (Read_Xmas == FALSE)
        eofstatus = fgets(line,LINELEN,fil_pdb);
      else
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type); */
        eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
/* <--v.4.2 */
    }

  /* Write out the number of pseudo H-bonds */
  printf("   Total number of CONECT records read in = %7d\n",*nlinks);

  /* Close the PDB file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
/* <--v.4.2 */
    fclose(fil_pdb);
}
/***********************************************************************

verify_conec_records  -  Verify the stored CONECT links by calculating
                         the appropriate distances

***********************************************************************/

/* v.4.2--> */
/* void verify_conec_records(char pdb_name[FILENAME_LEN]) */
void verify_conec_records(char pdb_name[FILENAME_LEN],XMAS *xmas)
/* <--v.4.2 */
{
/* v.4.2--> */
/*  char line[101]; */
  char line[LINELEN + 1];
/* <--v.4.2 */
  char atmnum[6];
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */

  int iatom;
  int delete_link, wanted;
/* v.4.2--> */
  int rec_types[NTYPES], type;
/* <--v.4.2 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */

  float dist, x, x1, x2, y, y1, y2, z, z1, z2;

  struct hhb_info *hhb_info_ptr, *previous_hhb_info_ptr;

  /* Initialise variables */
  previous_hhb_info_ptr = NULL;

  /* Open the pdb file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
    {
/* <--v.4.2 */
      if ((fil_pdb = fopen(pdb_name,"r")) == NULL)
        {
          printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
          exit(1);
        }
/* v.4.2--> */
    }
/* <--v.4.2 */
    
  /* Read in the data from the PDB file */
  printf("\nVerifying covalent bonds to ligand ...\n");

/* v.4.2--> */
  /* Get the first line from the PDB or XMAS file */
  if (Read_Xmas == FALSE)
    eofstatus = fgets(line,LINELEN,fil_pdb);
  else
    {
      /* Form the record types to be retrieved */
      rec_types[REC_ATOM] = TRUE;
      rec_types[REC_CONECT] = FALSE;
      type = REC_ATOM;

      /* Get first record from the XMAS cache */
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type); */
      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
    }
/* <--v.4.2 */

  /* Search through the PDB file to pick up the atoms involved in the
     CONECT records */
/* v.4.2--> */
/*  while (fgets(line,100,fil_pdb) != NULL) */
  while (eofstatus != NULL)
/* <--v.4.2 */
    {
      /* If this is an ATOM or HETATM record, process it */
      if ((!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)))
        {
          /* Extract the atom number */
          strncpy(atmnum,line+6,5);
          atmnum[5] = '\0';
          iatom = atoi(atmnum);

          /* Set pointer to the first of the links */
          hhb_info_ptr = first_hhb_info_ptr;

          /* Search through all the stored links to see if this atom
             is amongst them */
          while (hhb_info_ptr != NULL)
            {
              /* Check if the atom belongs to the current link */
              if (iatom == hhb_info_ptr->atom_number1 ||
                  iatom == hhb_info_ptr->atom_number2)
                {
                  /* Retrieve the atomic coordinates */
                  sscanf(line+30," %f %f %f ",&x,&y,&z);

                  /* Determine whether this is the first or second
                     atom of the pair */
                  if (iatom == hhb_info_ptr->atom_number1)
                    {
                      /* Transfer the appropriate details from the PDB
                         record to the H-bond info record */
                      strncpy(hhb_info_ptr->atom_name1,line+12,4);
                      hhb_info_ptr->atom_name1[4] = '\0';
                      strncpy(hhb_info_ptr->res_name1,line+17,3);
                      hhb_info_ptr->res_name1[3] = '\0';
                      strncpy(hhb_info_ptr->res_num1,line+22,5);
                      hhb_info_ptr->res_num1[5] = '\0';
                      hhb_info_ptr->chain1 = line[21];

                      /* Store this atom's coordinates */
                      hhb_info_ptr->x1 = x;
                      hhb_info_ptr->y1 = y;
                      hhb_info_ptr->z1 = z;
                    }

                  /* If this is the second atom of the link, store the
                     details */
                  else
                    {
                      /* Transfer the appropriate details from the PDB
                         record to the H-bond info record */
                      strncpy(hhb_info_ptr->atom_name2,line+12,4);
                      hhb_info_ptr->atom_name2[4] = '\0';
                      strncpy(hhb_info_ptr->res_name2,line+17,3);
                      hhb_info_ptr->res_name2[3] = '\0';
                      strncpy(hhb_info_ptr->res_num2,line+22,5);
                      hhb_info_ptr->res_num2[5] = '\0';
                      hhb_info_ptr->chain2 = line[21];

                      /* Store this atom's coordinates */
                      hhb_info_ptr->x2 = x;
                      hhb_info_ptr->y2 = y;
                      hhb_info_ptr->z2 = z;
                    }
                }

              /* Get the pointer to the next link */
              hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
            }
        }

/* v.4.2--> */
      /* Get the next line from the PDB or XMAS file */
      if (Read_Xmas == FALSE)
        eofstatus = fgets(line,LINELEN,fil_pdb);
      else
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type); */
        eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type,
                                      &mol_id);
/* <--v.4.2.2 */
/* <--v.4.2 */
    }

  /* Loop through the stored links to calculate distances and discard
     any links with missing atoms or which are clearly wrong */

  /* Set pointer to the first of the links */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Loop through all the stored links */
  while (hhb_info_ptr != NULL)
    {
      /* Initialise variables */
      delete_link = FALSE;

      /* If have located both atoms, then use stored coordinates to
         compute the distance between them */
      if (strncmp(hhb_info_ptr->atom_name1,"    ",4) &&
          strncmp(hhb_info_ptr->atom_name2,"    ",4))
        {
          /* Retrieve both atoms' coordinates */
          x1 = hhb_info_ptr->x1;
          y1 = hhb_info_ptr->y1;
          z1 = hhb_info_ptr->z1;
          x2 = hhb_info_ptr->x2;
          y2 = hhb_info_ptr->y2;
          z2 = hhb_info_ptr->z2;

          /* Calculate the distance between them */
          dist = (x1 - x2) * (x1 - x2)
            + (y1 - y2) * (y1 - y2)
              + (z1 - z2) * (z1 - z2);
          dist = (float) sqrt((double) dist);

          /* Check the bond-length */
          wanted = TRUE;
/* v.5.0.1--> */
/*          if (dist > 3.0) */
          if (dist > 3.0 && conect == FALSE)
/* <--v.5.0.1 */
            {
              printf("*** Warning. Long bond from CONECT ");
              printf("[%5d ->%5d] Length %7.2f",
                     hhb_info_ptr->atom_number1,
                     hhb_info_ptr->atom_number2,dist);
              if (dist > 4.0)
                {
                  printf(" - discarded\n");
                  wanted = FALSE;
                }
              else
                printf("\n");
              Nwarnings++;
            }

          /* Store the distance if not too long */
          if (wanted == TRUE)
            {
              hhb_info_ptr->bond_length = dist;
            }
          else
            delete_link = TRUE;
        }

      /* Otherwise, if one or both atoms not located, show warning message */
      else
        {
          /* Print warning message */
          printf("* Warning. Error in CONECT record joining atoms ");
          printf("%d and %d. \n",hhb_info_ptr->atom_number1,
                 hhb_info_ptr->atom_number2);
          printf("*               ");
          if (!strncmp(hhb_info_ptr->atom_name1,"    ",4))
            {
              printf("Atom %d not found. ",hhb_info_ptr->atom_number1);
            }
          if (!strncmp(hhb_info_ptr->atom_name2,"    ",4))
            {
              printf("Atom %d not found. ",hhb_info_ptr->atom_number2);
            }
          printf("\n");
          Nwarnings++;

          /* Mark for deletion */
          delete_link = TRUE;
        }

      /* If link to be deleted, do so by making pointers by-pass it */
      if (delete_link == TRUE)
        {
          /* If this is the first link, then set the next link as the
             first */
          if (previous_hhb_info_ptr == NULL)
            first_hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;

          /* Otherwise, make pointers bypass the current link */
          else
            previous_hhb_info_ptr->next_hhb_info_ptr
              = hhb_info_ptr->next_hhb_info_ptr;

          /* If this was the last in the list, then set the previous
             pointer as the new terminus */
          if (hhb_info_ptr->next_hhb_info_ptr == NULL)
            last_hhb_info_ptr = previous_hhb_info_ptr;
        }

      /* Save the current pointer */
      previous_hhb_info_ptr = hhb_info_ptr;

      /* Get the pointer to the next link */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }

  /* Close the PDB file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
/* <--v.4.2 */
    fclose(fil_pdb);
}
/* v.5.4--> */
/***********************************************************************

check_atom_pair  -  Check whether the two given atoms might form a salt
                    bridge

***********************************************************************/

int check_atom_pair(char *atom_name1,char *atom_name2,
                    char *res_name1,char *res_name2)
{
  char atom_name[5], res_name[4];

  int iatom, type[2];
  int salt_bridge, wanted;

  double dist_sqrd, min_dist_sqrd, x1, x2, y1, y2, z1, z2;

  /* Initialise variables */
  salt_bridge = FALSE;
  type[0] = type[1] = NOT_WANTED;

  /* Assume atom pair are wanted */
  wanted = TRUE;

  /* Loop over the two atoms to see if they are of the right type */
  for (iatom = 0; iatom < 2 && wanted == TRUE; iatom++)
    {
      /* Get this atom's name */
      if (iatom == 0)
        strcpy(atom_name,atom_name1);
      else
        strcpy(atom_name,atom_name2);

      /* Check whether this is the right kind of atom */
      if (!strcmp(atom_name," OD1") ||
          !strcmp(atom_name," OD2") ||
          !strcmp(atom_name," OE1") ||
          !strcmp(atom_name," OE2") ||
          !strcmp(atom_name," NE ") ||
          !strcmp(atom_name," NH1") ||
          !strcmp(atom_name," NH2") ||
          !strcmp(atom_name," NZ ") ||
          !strcmp(atom_name," ND1") ||
          !strcmp(atom_name," NE2"))
        wanted = TRUE;
      else
        wanted = FALSE;
    }

  /* If both atoms are wanted, check their residue types */
  if (wanted == TRUE)
    {
      /* Loop over the two atoms to get their residue types */
      for (iatom = 0; iatom < 2; iatom++)
        {
          /* Get this residue's name */
          if (iatom == 0)
            strcpy(res_name,res_name1);
          else
            strcpy(res_name,res_name2);

          /* Get residue type */
          if (!strcmp(res_name,"ASP") ||
              !strcmp(res_name,"GLU"))
            type[iatom] = NEGATIVE;
          else if (!strcmp(res_name,"ARG") ||
                   !strcmp(res_name,"LYS") ||
                   !strcmp(res_name,"HIS"))
            type[iatom] = POSITIVE;
          else
            type[iatom] = NOT_WANTED;
        }

      /* If the residues are oppositely charged, then we have a salt
         bridge */
      if ((type[0] == POSITIVE && type[1] == NEGATIVE) ||
          (type[1] == POSITIVE && type[0] == NEGATIVE))
        salt_bridge = TRUE;
    }
        
  /* Return whether atoms are of the correct types for salt-bridge */
  return(salt_bridge);
}
/* <--v.5.4 */
/***********************************************************************

get_hbplus_info  -  Extract the atom details from the interaction just
                    read in from the .hhb or .nnb file

***********************************************************************/

int get_hbplus_info(char *line,char *atom_name1,char *atom_name2,
                    char *chain1,char *chain2,char *res_name1,
                    char *res_name2,char *res_num1,char *res_num2,
                    float *bond_length,int file_type,int ligplot_format)
{
  char bond_lenstr[5];

/* v.5.4--> */
//  int ipos, nonzero, valid;
  int ipos, len, nonzero, bond_type;
  int salt_bridge;
/* <--v.5.4 */
/* v.5.3.7--> */
  int minus;
/* <--v.5.3.7 */

  /* Initialise variables */
/* v.5.4--> */
//  valid = TRUE;
  len = strlen(line);
  if (file_type == NNB_FILE)
    bond_type = S_CONTACTS;
  else
    bond_type = S_HBONDS;
/* <--v.5.4 */

  /* Check for obviously unwanted records */
/* v.4.0--> */
/*  if (!strncmp(line,"ligplot",7) || !strncmp(line + 10,"   ",10)) */
/* v.5.4--> */
//  if ((int) strlen(line) < 10 ||
  if ((int) strlen(line) < 25 ||
/* <--v.5.4 */
      (ligplot_format == TRUE && !strncmp(line + 21,"   ",3)) ||
      (ligplot_format == FALSE && !strncmp(line + 10,"   ",3)))
/* <--v.4.0 */
/* v.5.4--> */
//    valid = FALSE;
    bond_type = INVALID;
/* <--v.5.4 */

  /* Otherwise, extract the data according to the file format */
  else
    {
      /* LIGPLOT-format of the .hhb and .nnb files */
      if (ligplot_format == TRUE)
        {
          /* Atom name */
          strncpy(atom_name1,line+12,4);
          atom_name1[4] = '\0';
          strncpy(atom_name2,line+33,4);
          atom_name2[4] = '\0';
                
          /* Chain-id */
          *chain1 = line[4];
          *chain2 = line[25];
          if ((*chain1) == '-')
            *chain1 = ' ';
          if ((*chain2) == '-')
            *chain2 = ' ';
          
          /* Residue details */
          strncpy(res_num1,line+6,5);
          res_num1[5] = '\0';
          strncpy(res_num2,line+27,5);
          res_num2[5] = '\0';
          strncpy(res_name1,line,3);
          res_name1[3] = '\0';
          strncpy(res_name2,line+21,3);
          res_name2[3] = '\0';

          /* Bond length */
          strncpy(bond_lenstr,line+41,4);
          bond_lenstr[4] = '\0';
          *bond_length = (float) atof(bond_lenstr);

/* v.5.4--> */
          /* Check for salt-bridge or disulphide */
          if (len > 49)
            {
              /* Salt bridge */
              if (!strncmp(line + 47,"+/-",3))
                bond_type = S_SALT_BRIDGES;
              else if (!strncmp(line + 47,"S-S",3))
                bond_type = S_DISULPHIDES;
            }
/* <--v.5.4 */
        }

      /* HBPLUS-format hydrogen bonds file */
      else if (file_type == HHB_FILE || file_type == NNB_FILE)
        {
          /* Atom name */
          strncpy(atom_name1,line+19,4);
          atom_name1[4] = '\0';
          strncpy(atom_name2,line+39,4);
          atom_name2[4] = '\0';
                
          /* Chain-id */
          *chain1 = line[10];
          *chain2 = line[30];
          if ((*chain1) == '-')
            *chain1 = ' ';
          if ((*chain2) == '-')
            *chain2 = ' ';
          
          /* Residue details */
          strncpy(res_num1,line+11,5);
          res_num1[5] = '\0';
          strncpy(res_num2,line+31,5);
          res_num2[5] = '\0';
          strncpy(res_name1,line+16,3);
          res_name1[3] = '\0';
          strncpy(res_name2,line+36,3);
          res_name2[3] = '\0';

          /* Bond length */
          strncpy(bond_lenstr,line+45,4);
          bond_lenstr[4] = '\0';
          *bond_length = (float) atof(bond_lenstr);

/* v.5.4--> */
          /* If this is a non-bonded contact, Check whether it might
             correspond to a salt bridge */
          if (file_type == NNB_FILE && *bond_length <= SALT_BRIDGE_DIST)
            {
              /* Check this atom pair for a salt bridge */
              salt_bridge
                = check_atom_pair(atom_name1,atom_name2,
                                  res_name1,res_name2);
              if (salt_bridge == TRUE)
                bond_type = S_SALT_BRIDGES;
            }
/* <--v.5.4 */
        }

      /* For the residue numbers, may need to replace any leading
         zeros by spaces */

      /* First residue */
      nonzero = FALSE;
/* v.5.3.7--> */
      minus = FALSE;
//      for (ipos = 0; ipos < 3 && nonzero == FALSE; ipos ++)
      for (ipos = 0; ipos < 4 && nonzero == FALSE; ipos ++)
/* <--v.5.3.7 */
        {
/* v.5.3.7--> */
//          if (res_num1[ipos] == '0')
          if (res_num1[ipos] == '0' && ipos < 3)
/* <--v.5.3.7 */
            res_num1[ipos] = ' ';
/* v.5.3.7--> */
          else if (res_num1[ipos] == '-')
            {
              minus = TRUE;
              res_num1[ipos] = ' ';
            }
/* <--v.5.3.7 */
          else
/* v.5.3.7--> */
            {
/* <--v.5.3.7 */
              nonzero = TRUE;
/* v.5.3.7--> */
              if (minus == TRUE && ipos > 0)
                res_num1[ipos - 1] = '-';
            }
/* <--v.5.3.7 */
        }

      /* Replace hyphen in insertion code with a space */
      if (res_num1[4] == '-')
        res_num1[4] = ' ';
                
      /* Repeat for second residue */
      nonzero = FALSE;
/* v.5.3.7--> */
      minus = FALSE;
//      for (ipos = 0; ipos < 3 && nonzero == FALSE; ipos ++)
      for (ipos = 0; ipos < 4 && nonzero == FALSE; ipos ++)
/* <--v.5.3.7 */
        {
/* v.5.3.7--> */
//          if (res_num2[ipos] == '0')
          if (res_num2[ipos] == '0' && ipos < 3)
/* <--v.5.3.7 */
            res_num2[ipos] = ' ';
/* v.5.3.7--> */
          else if (res_num2[ipos] == '-')
            {
              minus = TRUE;
              res_num2[ipos] = ' ';
            }
/* <--v.5.3.7 */
          else
/* v.5.3.7--> */
            {
/* <--v.5.3.7 */
              nonzero = TRUE;
/* v.5.3.7--> */
              if (minus == TRUE && ipos > 0)
                res_num2[ipos - 1] = '-';
            }
/* <--v.5.3.7 */
        }

      /* Replace hyphen in insertion code with a space */
      if (res_num2[4] == '-')
        res_num2[4] = ' ';
    }

  /* Return interaction type */
/* v.5.4--> */
//  return(valid);
  return(bond_type);
/* <--v.5.4 */
}
/***********************************************************************

read_hbplus_file  -  Read in the required H-bond and nonbonded contacts
                     data from the appropriate HBPLUS output file

***********************************************************************/

void read_hbplus_file(char file_name[FILENAME_LEN],int file_type,
/* v.4.1.2--> */
/*                    int *nbonds) */
/* v.4.2--> */
/*                    int *nbonds,int debug) */
                      int *nbonds,int debug,XMAS *xmas)
/* <--v.4.2 */
/* <--v.4.1.2 */
{
  char file_desc[20], line[LINELEN + 1];
  char atom_name1[5], atom_name2[5], chain1, chain2, res_name1[4],
  res_name2[4], res_num1[6], res_num2[6];
  char dummy[4];
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */
/* v.4.3--> */
  char element[3];
/* <--v.4.3 */

/* v.5.4--> */
//  int have_file, inrange1, inrange2, ligplot_format, valid, wanted;
  int have_file, inrange1, inrange2, ligplot_format, wanted;
/* <--v.5.4 */
/* v.4.1.2--> */
  int nline;
/* <--v.4.1.2 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */
/* v.5.4--> */
  int bond_type;
/* <--v.5.4 */

  float bond_length;

  struct hhb_info *hhb_info_ptr;

  /* Initialise variables */
  have_file = TRUE;
/* v.5.4--> */
//  *nbonds = 0;
/* <--v.5.4 */
  ligplot_format = FALSE;
/* v.4.1.2--> */
  nline = 0;
/* <--v.4.1.2 */
/* v.4.2.2--> */
  mol_id = -1;
/* <--v.4.2.2 */
  if (file_type == HHB_FILE)
    strcpy(file_desc,"hydrogen bonds");
  else
    strcpy(file_desc,"non-bonded contacts");

  /* Open the HBPLUS input file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
    {
/* <--v.4.2 */
      if ((fil_hbplus = fopen(file_name,"r")) == NULL)
        {
/* v.4.3--> */
/*          printf("\n*** Unable to open your %s file [%s]\n",file_desc,
                 file_name);
          printf("*** No %s will be shown on the plot\n",file_desc); */
          printf("\n* Unable to open your %s file [%s]\n",file_desc,
                 file_name);
          printf("* No %s will be shown on the plot\n",file_desc);
/* <--v.4.3 */
          have_file = FALSE;
          Nwarnings++;
        }
/* v.4.2--> */
    }
/* <--v.4.2 */

  /* If input file exists, read in the interactions */
/* v.4.2--> */
/*  if (have_file == TRUE) */
  if (have_file == TRUE || Read_Xmas == TRUE)
/* <--v.4.2 */
    {
/* v.4.2--> */
/*      printf("\nReading in %s from file %s ...\n",file_desc,file_name); */

      /* Get the first line from the HBPLUS or XMAS file */
      if (Read_Xmas == FALSE)
        {
          printf("\nReading in %s from file %s ...\n",file_desc,file_name);
          eofstatus = fgets(line,LINELEN,fil_hbplus);
        }
      else
        {
          printf("\nGetting XMAS %s from file ...\n",file_desc);
          eofstatus = get_xmas2hb_line(xmas,line,TRUE,file_type);
        }
/* <--v.4.2 */

      /* Loop while reading through the file */
/* v.4.2--> */
/*      while (fgets(line,LINELEN,fil_hbplus) != NULL) */
      while (eofstatus != NULL)
/* <--v.4.2 */
        {
/* v.4.1.2--> */
          /* Increment line-count */
          nline++;
/* <--v.4.1.2 */

          /* Check if this line identifies the file as being in
             LIGPLOT format */
          if (!strncmp(line,"ligplot",7))
            ligplot_format = TRUE;

          /* Extract the relevant information from the line just
             read in */
/* v.5.4--> */
//          valid = get_hbplus_info(line,atom_name1,atom_name2,&chain1,
          bond_type = get_hbplus_info(line,atom_name1,atom_name2,&chain1,
/* <--v.5.4 */
                                      &chain2,res_name1,res_name2,
                                      res_num1,res_num2,&bond_length,
                                      file_type,ligplot_format);

          /* If this is a valid interaction check whether it is
             wanted */
/* v.5.4--> */
//          if (valid == TRUE)
          if (bond_type != INVALID)
/* <--v.5.4 */
            {
              /* Initialise flag which determines whether this
                 interaction is required */
              wanted = TRUE;

              /* For non-bonded contacts, check whether this is a desired
                 contact */
/* v.5.4--> */
//              if (file_type == NNB_FILE)
              if (file_type == NNB_FILE && bond_type != S_SALT_BRIDGES)
/* <--v.5.4 */
                {
                  /* If hydrophobics only (as in the old version) check
                     that both atoms are either a carbon or a sulphur */
                  if (Include->Contact_Type == HYDROPHOBIC_ONLY)
                    {
                      if ((atom_name1[1] != 'C' && atom_name1[1] != 'S') ||
                          (atom_name2[1] != 'C' && atom_name2[1] != 'S'))
                        wanted = FALSE;
                    }

                  /* If at least one of the atoms has to be a hydrophobic
                     then check for this */
                  else if (Include->Contact_Type == HYDROPHOBIC_ANY)
                    {
                      if (atom_name1[1] == 'C' || atom_name1[1] == 'S' ||
                          atom_name2[1] == 'C' || atom_name2[1] == 'S')
                        wanted = TRUE;
                      else
                        wanted = FALSE;
                    }

                  /* Otherwise, accept any contact */
                  else
                    wanted = TRUE;

                  /* If either atom is a water, then don't want a
                     non-bonded contact */
                  if (!strncmp(res_name1,"HOH",3) ||
                      !strncmp(res_name2,"HOH",3))
                    wanted = FALSE;

                  /* If either is a metal, then again don't want this
                     contact */
                  if (wanted == TRUE)
                    {
/* v.4.3--> */
/*                    if (check_for_metal(atom_name1) == TRUE)
                        wanted = FALSE;
                      else if (check_for_metal(atom_name2) == TRUE)
                        wanted = FALSE; */
                      if (check_for_metal(atom_name1,element) == TRUE)
                        wanted = FALSE;
                      else if (check_for_metal(atom_name2,element) == TRUE)
                        wanted = FALSE;
/* <--v.4.3 */
                    }
                }

              /* If waters not required and one of the atoms is a
                 water, then exclude */
              if (Include->Waters == FALSE &&
                  (!strncmp("HOH",res_name1,3) ||
                   !strncmp("HOH",res_name2,3)))
                wanted = FALSE;

              /* Check if either residue is in the ligand */
              else if (wanted == TRUE)
                {
/* v.4.0--> */
                  /* If plotting interface interactions, then check whether
                     this interaction is wanted */
                  if (Interface_Plot == TRUE)
                    {
                      /* Set interaction as wanted */
                      wanted = TRUE;
                    }

                  /* Otherwise check whether one of the residues belongs
                     to a ligand */
                  else
                    {
/* <--v.4.0 */
                      /* Check if first residue belongs to the ligand */
                      inrange1 = check_if_in_ligand(chain1,res_name1,
/* v.4.2.2--> */
/*                                                  res_num1,FALSE,dummy); */
                                                    res_num1,FALSE,dummy,
                                                    mol_id);
/* <--v.4.2.2 */
                      /* Repeat for second residue */
                      inrange2 = check_if_in_ligand(chain2,res_name2,
/* v.4.2.2--> */
/*                                                  res_num2,FALSE,dummy); */
                                                    res_num2,FALSE,dummy,
                                                    mol_id);
/* <--v.4.2.2 */
                    
                      /* If either residue falls within the range defining the
                         ligand, then want to save this interaction */
                      if (inrange1 == FALSE && inrange2 == FALSE)
                        wanted = FALSE;
                    }
                }

              /* If this interaction is still wanted, then store it */
              if (wanted == TRUE)
                {
                  /* Allocate memory for structure to hold current 
                     interaction */
                  hhb_info_ptr = (struct hhb_info *)
                    malloc(sizeof(struct hhb_info));

                  /* If this is the first interaction stored, then
                      update pointer to head of linked list */
                  if (first_hhb_info_ptr == NULL)
                    first_hhb_info_ptr = hhb_info_ptr;
              
                  /* Otherwise, make previous interaction point to the
                     current one */
                  else
                    last_hhb_info_ptr->next_hhb_info_ptr = hhb_info_ptr;
              
                  /* Save the current interaction's pointer */
                  last_hhb_info_ptr = hhb_info_ptr;

                  /* Store the data */
                  hhb_info_ptr->source = file_type;
/* v.5.4--> */
                  hhb_info_ptr->bond_type = bond_type;
/* <--v.5.4 */
                  hhb_info_ptr->atom_number1 = 0;
                  hhb_info_ptr->atom_number2 = 0;

                  /* Initialise the other fields in the structure */

                  /* First atom's details */
                  strncpy(hhb_info_ptr->atom_name1,atom_name1,4);
                  hhb_info_ptr->atom_name1[4] = '\0';
                  strncpy(hhb_info_ptr->res_name1,res_name1,3);
                  hhb_info_ptr->res_name1[3] = '\0';
                  strncpy(hhb_info_ptr->res_num1,res_num1,5);
                  hhb_info_ptr->res_num1[5] = '\0';
                  hhb_info_ptr->chain1 = chain1;
                  hhb_info_ptr->x1 = 0.0;
                  hhb_info_ptr->y1 = 0.0;
                  hhb_info_ptr->z1 = 0.0;

                  /* Second atom's details */
                  strncpy(hhb_info_ptr->atom_name2,atom_name2,4);
                  hhb_info_ptr->atom_name2[4] = '\0';
                  strncpy(hhb_info_ptr->res_name2,res_name2,3);
                  hhb_info_ptr->res_name2[3] = '\0';
                  strncpy(hhb_info_ptr->res_num2,res_num2,5);
                  hhb_info_ptr->res_num2[5] = '\0';
                  hhb_info_ptr->chain2 = chain2;
                  hhb_info_ptr->x2 = 0.0;
                  hhb_info_ptr->y2 = 0.0;
                  hhb_info_ptr->z2 = 0.0;

                  /* Distance between atoms */
                  hhb_info_ptr->bond_length = bond_length;

                  /* Additional data fields */
                  hhb_info_ptr->next_hhb_info_ptr = NULL;

                  /* Increment count of bonds of this type */
/* v.5.4--> */
//                  (*nbonds)++;
                  nbonds[bond_type]++;
/* <--v.5.4 */

/* v.4.1.2--> */
                  /* If in debug mode, print out accepted line */
                  if (debug == TRUE)
                    {
                      printf("  ------ Stored ");
/* v.5.4--> */
//                      if (file_type == HHB_FILE)
//                        printf("H-bond ");
//                      else
//                        printf("contact ");
                      if (bond_type == S_HBONDS)
                        printf("H-bond ");
                      else if (bond_type == S_CONTACTS)
                        printf("contact ");
                      else if (bond_type == S_SALT_BRIDGES)
                        printf("salt br ");
                      else if (bond_type == S_DISULPHIDES)
                        printf("disulph ");
/* <--v.5.4 */
                      printf("[%s %s %s %c] -> [%s %s %s %c]\n",
                             hhb_info_ptr->atom_name1,
                             hhb_info_ptr->res_name1,
                             hhb_info_ptr->res_num1,
                             hhb_info_ptr->chain1,
                             hhb_info_ptr->atom_name2,
                             hhb_info_ptr->res_name2,
                             hhb_info_ptr->res_num2,
                             hhb_info_ptr->chain2);
                    }
/* <--v.4.1.2 */
                }
            }

/* v.4.2--> */
          /* Get the next line from the HBPLUS or XMAS file */
          if (Read_Xmas == FALSE)
            eofstatus = fgets(line,LINELEN,fil_hbplus);
          else
            eofstatus = get_xmas2hb_line(xmas,line,FALSE,file_type);
/* <--v.4.2 */
        }

      /* Close the input file */
/* v.4.2--> */
      if (Read_Xmas == FALSE)
/* <--v.4.2 */
        fclose(fil_hbplus);
    }
}
/***********************************************************************

special_resinc  -  Include specially defined residues in the parameter
                   file - ie  residues that are connected to other
                   residues connected to the ligand

***********************************************************************/

/* v.4.2--> */
/* void special_resinc(char file_name[FILENAME_LEN],int file_type,
                    char special_res[MAXSPECIAL_RES][15]) */
void special_resinc(char file_name[FILENAME_LEN],int file_type,
                    char special_res[MAXSPECIAL_RES][15],XMAS *xmas)
/* <--v.4.2 */
{
  char file_desc[20], line[LINELEN + 1];
  char atom_name1[5], atom_name2[5], chain1, chain2, res_name1[4],
  res_name2[4], res_num1[6], res_num2[6];
/* v.5.1--> */
  char metal_name[3];
/* <--v.5.1 */
  char check_chain, check_res_name[4], check_res_num[6], other_res_name[4];
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */

/* v.5.4--> */
//  int have_file, ligplot_format, valid;
  int bond_type;
  int have_file, ligplot_format;
/* <--v.5.4 */
  int match1, match2, res_to_ligand, nspecial_res, sp_res;
/* v.5.1--> */
  int is_metal1, is_metal2;
/* <--v.5.1 */

  float bond_length;

  struct hhb_info *hhb_info_ptr;

  /* Initialise variables */
  have_file = TRUE;
  ligplot_format = FALSE;
  if (file_type == HHB_FILE)
    strcpy(file_desc,"hydrogen bonds");
  else
    strcpy(file_desc,"non-bonded contacts");

  /* Open the HBPLUS input file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
    {
/* <--v.4.2 */
      if ((fil_hbplus = fopen(file_name,"r")) == NULL)
        {
/* v.4.3--> */
/*        printf("\n*** Unable to open your %s file [%s]\n",file_desc,
                 file_name);
          printf("*** No %s will be shown on the plot\n",file_desc); */
        printf("\n* Unable to open your %s file [%s]\n",file_desc,
                 file_name);
          printf("* No %s will be shown on the plot\n",file_desc);
/* <--v.4.3 */
          have_file = FALSE;
          Nwarnings++;
        }
/* v.4.2--> */
    }
/* <--v.4.2 */

  /* If input file exists, read in the interactions */
/* v.4.2--> */
/*  if (have_file == TRUE) */
  if (have_file == TRUE || Read_Xmas == TRUE)
/* <--v.4.2 */
    {
/* v.4.2--> */
/*      printf("\nReading in %s from file %s ...\n",file_desc,file_name); */
/* <--v.4.2 */

      /* Determine how many special-residue records there are */
      nspecial_res = 0;
      for (sp_res = 0; sp_res < MAXSPECIAL_RES; sp_res++)
        {
          /* If entry non-blank, then save its position */
          if (strncmp(special_res[sp_res],"       ",7))
            nspecial_res = sp_res + 1;
        }

/* v.4.2--> */
      /* Get the first line from the HBPLUS or XMAS file */
      if (Read_Xmas == FALSE)
        {
          printf("\nReading in %s from file %s ...\n",file_desc,file_name);
          eofstatus = fgets(line,LINELEN,fil_hbplus);
        }
      else
        {
          printf("\nGetting XMAS %s from file ...\n",file_desc);
          eofstatus = get_xmas2hb_line(xmas,line,TRUE,file_type);
        }
/* <--v.4.2 */

      /* Loop while reading through the file */
/* v.4.2--> */
/*      while (fgets(line,LINELEN,fil_hbplus) != NULL) */
      while (eofstatus != NULL)
/* <--v.4.2 */
        {
          /* Check if this line identifies the file as being in
             LIGPLOT format */
          if (!strncmp(line,"ligplot",7))
            ligplot_format = TRUE;

          /* Extract the relevant information from the line just
             read in */
/* v.5.4--> */
//          valid = get_hbplus_info(line,atom_name1,atom_name2,&chain1,
          bond_type = get_hbplus_info(line,atom_name1,atom_name2,&chain1,
/* <--v.5.4 */
                                      &chain2,res_name1,res_name2,
                                      res_num1,res_num2,&bond_length,
                                      file_type,ligplot_format);

          /* If this is a valid interaction check whether it is
             wanted */
/* v.5.4--> */
//          if (valid == TRUE)
          if (bond_type != INVALID)
/* <--v.5.4 */
            {
/* v.5.1--> */
              /* Check whether either residue is a metal ion */
              is_metal1 = is_metal2 = FALSE;
              if (file_type == HHB_FILE)
                {
                  /* Check each atom to see if it a metal */
                  is_metal1 = check_for_metal(atom_name1,metal_name);
                  is_metal2 = check_for_metal(atom_name2,metal_name);
                }
/* <--v.5.1 */
              /* Check if this interaction matches any of the selection
                 criteria */
              for (sp_res = 0; sp_res < nspecial_res; sp_res++)
                {
                  /* For the case where two residues that are non-ligand 
                     and hydrogen bonded to the ligand and the user 
                     wishes to H-bond these together */
                  match1 = FALSE;
/* v.5.1--> */
/*                  if ((!strncmp(special_res[sp_res],res_name1,3) ||
                       !strncmp(special_res[sp_res],"***",3)) &&
                      (!strncmp(special_res[sp_res]+4,res_name2,3) ||
                      !strncmp(special_res[sp_res]+4,"***",3))) */
                  if ((!strncmp(special_res[sp_res],res_name1,3) ||
                       !strncmp(special_res[sp_res],"***",3) ||
                       (is_metal1 == TRUE &&
                        !strncmp(special_res[sp_res],"met",3))) &&
                      (!strncmp(special_res[sp_res]+4,res_name2,3) ||
                       !strncmp(special_res[sp_res]+4,"***",3) ||
                       (is_metal2 == TRUE &&
                        !strncmp(special_res[sp_res]+4,"met",3))))
/* <--v.5.1 */
                    match1 = TRUE;
                    
                  match2 = FALSE;
/* v.5.1--> */
/*                  if ((!strncmp(special_res[sp_res],res_name2,3) ||
                       !strncmp(special_res[sp_res],"***",3)) &&
                      (!strncmp(special_res[sp_res]+4,res_name1,3) ||
                      !strncmp(special_res[sp_res]+4,"***",3))) */
                  if ((!strncmp(special_res[sp_res],res_name2,3) ||
                       !strncmp(special_res[sp_res],"***",3) ||
                       (is_metal2 == TRUE &&
                        !strncmp(special_res[sp_res],"met",3))) &&
                      (!strncmp(special_res[sp_res]+4,res_name1,3) ||
                       !strncmp(special_res[sp_res]+4,"***",3) ||
                       (is_metal1 == TRUE &&
                        !strncmp(special_res[sp_res]+4,"met",3))))
/* <--v.5.1 */
                    match2 = TRUE;
                  
                  /* If have a possible match, need to check if either
                     of the residues is H-bonded to the ligand */
                  if (match1 == TRUE || match2 == TRUE)
                    {
                      /* Need to loop through the stored H-bond links
                         to check if either residue is H-bonded to the
                         ligand */
                      res_to_ligand = FALSE;

                      /* Store which of the two residues needs to be
                         H-bonded to the ligand */
                      if (match1 == TRUE)
                        {
                          strncpy(check_res_name,res_name1,3);
                          check_res_name[3] = '\0';
                          strncpy(check_res_num,res_num1,5);
                          check_res_num[5] = '\0';
                          check_chain = chain1;
                          strncpy(other_res_name,res_name2,3);
                          other_res_name[3] = '\0';
                        }
                      else
                        {
                          strncpy(check_res_name,res_name2,3);
                          check_res_name[3] = '\0';
                          strncpy(check_res_num,res_num2,5);
                          check_res_num[5] = '\0';
                          check_chain = chain2;
                          strncpy(other_res_name,res_name1,3);
                          other_res_name[3] = '\0';
                        }

                      /* Set pointer to the first of the H-bond links
                         stored so far */
                      hhb_info_ptr = first_hhb_info_ptr;
                      
                      /* Search through all the stored links */
                      while (hhb_info_ptr != NULL && res_to_ligand == FALSE)
                        {
                          /* Check if first residue is H-bonded to
                             the ligand */
                          if (hhb_info_ptr->source == HHB_FILE &&
                              ((!strncmp(hhb_info_ptr->res_name1,
                                         check_res_name,3) &&
                                !strncmp(hhb_info_ptr->res_num1,
                                         check_res_num,5) &&
                                hhb_info_ptr->chain1 == check_chain) ||
                               (!strncmp(hhb_info_ptr->res_name2,
                                         check_res_name,3) &&
                                !strncmp(hhb_info_ptr->res_num2,
                                         check_res_num,5) &&
                                hhb_info_ptr->chain2 == check_chain)))
                            res_to_ligand = TRUE;
                          
                          /* Get the pointer to the next link */
                          hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
                        }

                      /* If the outer residue is a water, then only
                         include if waters to be shown on plot */
                      if (!strncmp(other_res_name,"HOH",3) &&
                          Include->Waters == FALSE)
                        res_to_ligand = FALSE;

                      /* If have found an interaction between the
                         residue and the ligand, want to add the current
                         HBPLUS interaction to the list of hydrogen bonds */
                      if (res_to_ligand == TRUE)
                        {
                          /* Allocate memory for structure to hold current 
                             interaction */
                          hhb_info_ptr = (struct hhb_info *)
                            malloc(sizeof(struct hhb_info));

                          /* If this is the first interaction stored, then
                             update pointer to head of linked list */
                          if (first_hhb_info_ptr == NULL)
                            first_hhb_info_ptr = hhb_info_ptr;

                          /* Otherwise, make previous interaction point
                             to the current one */
                          else
                            last_hhb_info_ptr->next_hhb_info_ptr
                              = hhb_info_ptr;
              
                          /* Save the current interaction's pointer */
                          last_hhb_info_ptr = hhb_info_ptr;

                          /* Store the data */
                          hhb_info_ptr->source = file_type;
                          hhb_info_ptr->atom_number1 = 0;
                          hhb_info_ptr->atom_number2 = 0;
/* v.5.4.6--> */
                          hhb_info_ptr->bond_type = S_HBONDS;
/* <--v.5.4.6 */

                          /* Initialise the other fields in the structure */

                          /* First atom's details */
                          strncpy(hhb_info_ptr->atom_name1,atom_name1,4);
                          hhb_info_ptr->atom_name1[4] = '\0';
                          strncpy(hhb_info_ptr->res_name1,res_name1,3);
                          hhb_info_ptr->res_name1[3] = '\0';
                          strncpy(hhb_info_ptr->res_num1,res_num1,5);
                          hhb_info_ptr->res_num1[5] = '\0';
                          hhb_info_ptr->chain1 = chain1;
                          hhb_info_ptr->x1 = 0.0;
                          hhb_info_ptr->y1 = 0.0;
                          hhb_info_ptr->z1 = 0.0;

                          /* Second atom's details */
                          strncpy(hhb_info_ptr->atom_name2,atom_name2,4);
                          hhb_info_ptr->atom_name2[4] = '\0';
                          strncpy(hhb_info_ptr->res_name2,res_name2,3);
                          hhb_info_ptr->res_name2[3] = '\0';
                          strncpy(hhb_info_ptr->res_num2,res_num2,5);
                          hhb_info_ptr->res_num2[5] = '\0';
                          hhb_info_ptr->chain2 = chain2;
                          hhb_info_ptr->x2 = 0.0;
                          hhb_info_ptr->y2 = 0.0;
                          hhb_info_ptr->z2 = 0.0;

                          /* Distance between atoms */
                          hhb_info_ptr->bond_length = bond_length;

                          /* Additional data fields */
                          hhb_info_ptr->next_hhb_info_ptr = NULL;
                        }
                    }
                }
            }

/* v.4.2--> */
          /* Get the next line from the HBPLUS or XMAS file */
          if (Read_Xmas == FALSE)
            eofstatus = fgets(line,LINELEN,fil_hbplus);
          else
            eofstatus = get_xmas2hb_line(xmas,line,FALSE,file_type);
/* <--v.4.2 */
        }

      /* Close the input file */
/* v.4.2--> */
      if (Read_Xmas == FALSE)
/* <--v.4.2 */
      fclose(fil_hbplus);
    }
}
/***********************************************************************

update_bond  -  Add/remove bond to/from the bond-list, or update its
                bond-type

***********************************************************************/

void update_bond(struct coordinate *first_atom_ptr, 
                 struct coordinate *second_atom_ptr,int bond_type,
                 int bond_source,float bond_length)
{
  int got_it_already;

  struct bond *bond_ptr;

  /* Initialise variables */
/*  bond_order = SINGLE; */

  /* Check that bond not already stored */

  /* Loop through the stored bond-list to see if already have bond */
  bond_ptr = first_bond_ptr;
  got_it_already = FALSE;
  while (bond_ptr != NULL && got_it_already == FALSE)
    {
      /* If this bond contains both atoms, then have it already */
      if ((bond_ptr->first_atom_ptr == first_atom_ptr &&
           bond_ptr->second_atom_ptr == second_atom_ptr) ||
          (bond_ptr->first_atom_ptr == second_atom_ptr &&
           bond_ptr->second_atom_ptr == first_atom_ptr))
        {
          got_it_already = TRUE;

          /* If bond to be deleted, mark it so */
          if (bond_type == DELETED)
            bond_ptr->bond_type = DELETED;

/* v.5.0.1--> */
          /* If new bond is an H-bond, and old version is a non-bonded
             contact, then upgrade to H-bond */
/* v.5.4--> */
//          if (bond_type == HBOND && bond_ptr->bond_type == CONTACT)
//            bond_ptr->bond_type = HBOND;
          if ((bond_type == HBOND ||
               bond_type == SALT_BRIDGE ||
               bond_type == DISULPHIDE) && bond_ptr->bond_type == CONTACT)
            bond_ptr->bond_type = bond_type;
/* <--v.5.4 */
/* <--v.5.0.1 */
        }

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* If bond is a new bond to be added, then do so */
  if (bond_type != DELETED && got_it_already == FALSE)
    {
      /* Create a new bond record */
/* v.4.3--> */
      bond_ptr = create_bond_record(first_bond_ptr,bond_type,
                                    bond_source,bond_length,
                                    first_atom_ptr,second_atom_ptr);
/* <--v.4.3 */

      /* Store this pointer as the last-encountered bond */
      first_bond_ptr = bond_ptr;
    }
}
/* v.5.0--> */
/***********************************************************************

find_hhb_info  -  Check whether we already have the give H-bond or
                  non-bonded contact

***********************************************************************/

/* v.5.0.1--> */
/* int find_hhb_info(char *atom_name1,char *atom_name2,char chain1,
                  char chain2,char *res_name1,char *res_name2,
                  char *res_num1,char *res_num2,int bond_type) */
struct hhb_info *find_hhb_info(char *atom_name1,char *atom_name2,
                               char chain1,char chain2,char *res_name1,
                               char *res_name2,char *res_num1,
                               char *res_num2,int bond_type)
/* <--v.5.0.1 */
{
  int found = FALSE;

/* v.5.0.1--> */
/*  struct hhb_info *hhb_info_ptr; */
  struct hhb_info *hhb_info_ptr, *found_hhb_info_ptr;
/* <--v.5.0.1 */

/* v.5.0.1--> */
  /* Initialise variables */
  found_hhb_info_ptr = NULL;
/* <--v.5.0.1 */

  /* Set pointer to the first of the hhb_info records */
  hhb_info_ptr = first_hhb_info_ptr;
                     
  /* Loop through all the bonds to see if we have this one */
  while (hhb_info_ptr != NULL && found == FALSE)
    {
      /* Check if this is our bond */
/* v.5.4--> */
//      if (hhb_info_ptr->source == bond_type &&
      if (hhb_info_ptr->bond_type == bond_type &&
/* <--v.5.4 */
          ((!strcmp(hhb_info_ptr->res_num1,res_num1) &&
            !strcmp(hhb_info_ptr->res_num2,res_num2) &&
            !strcmp(hhb_info_ptr->atom_name1,atom_name1) &&
/* v.5.3.3 --> */
            !strcmp(hhb_info_ptr->atom_name2,atom_name2) &&
/* <--v.5.3.3 */
            !strcmp(hhb_info_ptr->res_name1,res_name1) &&
            !strcmp(hhb_info_ptr->res_name2,res_name2) &&
            hhb_info_ptr->chain1 == chain1 &&
            hhb_info_ptr->chain2 == chain2) ||
           (!strcmp(hhb_info_ptr->res_num1,res_num2) &&
            !strcmp(hhb_info_ptr->res_num2,res_num1) &&
            !strcmp(hhb_info_ptr->atom_name1,atom_name2) &&
/* v.5.3.3 --> */
            !strcmp(hhb_info_ptr->atom_name2,atom_name1) &&
/* <--v.5.3.3 */
            !strcmp(hhb_info_ptr->res_name1,res_name2) &&
            !strcmp(hhb_info_ptr->res_name2,res_name1) &&
            hhb_info_ptr->chain1 == chain2 &&
            hhb_info_ptr->chain2 == chain1)))
        {
          /* Set flag that bond found */
          found = TRUE;

/* v.5.0.1--> */
          /* Save the bond */
          found_hhb_info_ptr = hhb_info_ptr;
/* <--v.5.0.1 */
        }
                          
      /* Get the pointer to the next record */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }

  /* Return whether bond was found */
/* v.5.0.1--> */
/*  return(found); */
  return(found_hhb_info_ptr);
/* <--v.5.0.1 */
}
/* <--v.5.0 */
/* v.5.0.1--> */
/***********************************************************************

delete_hhb_info  -  Delete the given hhb_info record

***********************************************************************/

void delete_hhb_info(struct hhb_info *delete_hhb_info_ptr,
                     struct hhb_info **fst_hhb_info_ptr)
{
  int found;

  struct hhb_info *hhb_info_ptr, *first_hhb_info_ptr;
  struct hhb_info *next_hhb_info_ptr, *previous_hhb_info_ptr;

  /* Initialise pointer to the first record */
  hhb_info_ptr = first_hhb_info_ptr = *fst_hhb_info_ptr;
  previous_hhb_info_ptr = NULL;
  found = FALSE;

  /* Loop through all hhb_infos to find the one to be deleted */
  while (hhb_info_ptr != NULL && found == FALSE)
    {
      /* Save pointer to the next record */
      next_hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;

      /* If this is the one to be deleted, then get rid of it */
      if (hhb_info_ptr == delete_hhb_info_ptr)
        {
          /* Get previous record to point round this one */
          if (previous_hhb_info_ptr != NULL)
            previous_hhb_info_ptr->next_hhb_info_ptr
              = next_hhb_info_ptr;

          /* If no previous, then adjust first pointer */
          else
            first_hhb_info_ptr = next_hhb_info_ptr;

          /* Free up the memory taken by this record */
          free(hhb_info_ptr);

          /* Set flag that we're done */
          found = TRUE;
        }

      /* Save the current hhb_info pointer */
      previous_hhb_info_ptr = hhb_info_ptr;

      /* Get pointer to point to next hhb_info in the list and increment
         hhb_info count */
      hhb_info_ptr = next_hhb_info_ptr;
    }

  /* Return the first pointer */
  *fst_hhb_info_ptr = first_hhb_info_ptr;
}
/* <--v.5.0.1 */
/***********************************************************************

read_pdb_file  -  Read in the required atom coordinates from the
                  PDB file

***********************************************************************/

/* v.4.2--> */
/* int read_pdb_file(char pdb_name[FILENAME_LEN]) */
int read_pdb_file(char pdb_name[FILENAME_LEN],XMAS *xmas)
/* <--v.4.2 */
{
  char line[LINELEN + 1];
/* v.5.4.4 --> */
  char interface_name[LINELEN + 1];
/* <-- v.5.4.4 */
/* v.4.0--> */
/*  char accstring[8], res_name[4], res_num[6]; */
  char accstring[10], res_name[4], res_num[6];
/* <--v.4.0 */
  char atmnum[6];
  char dummy[4], last_res_name[4], occup[4];
  char last_res_num[6];
  char atom_partner[6];
  char atom_name[5], chain, last_chain;
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */
/* v.4.3--> */
  char atom_number[6], element[3];
/* <--v.4.3 */
  int asa_file, atom_counter, residue_counter, counter;
  int cpos, ipos, inligand, new_residue, found, keep_reading,
  want_atom, residue_in_ligand, residue_in_group, natoms_in_ligand,
  wanted;
  int iamino, ncopies, residues_stored;
  int hydrogen;
  int atom_no, last_atom_no, nconnect;
/* v.4.0--> */
  int in_lig1, in_lig2, ligplot_pdb_file, surface1;
  int len;
/* <--v.4.0 */
/* v.4.2--> */
  int rec_types[NTYPES], type;
/* <--v.4.2 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */
/* v.4.3--> */
  int asa_header;
/* <--v.4.3 */
/* v.5.0.1--> */
  int delete;
/* <--v.5.0.1 */
/* v.5.0--> */
  char atom_name1[5], atom_name2[5], chain1, chain2;
  char res_name1[4], res_name2[4], res_num1[6], res_num2[6];

/* v.5.4.4 --> */
  char *string_ptr;
/* <-- v.5.4.4 */

  int atom_number1, atom_number2;
  int file_type;
/* v.5.4--> */
  int bond_type;
//  int valid;
/* <--v.5.4 */

  float bond_length;
/* <--v.5.0 */

  float accessibility, length, length2, x1, x2, y1, y2, z1, z2;
/* v.4.3--> */
  float x, y, z;
/* <--v.4.3 */

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr, *last_atom_ptr,
  *other_atom_ptr;
  struct hhb_info *hhb_info_ptr;
  struct residue *residue_ptr, *last_residue_ptr, *other_residue_ptr;

  static char *amino[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL"
  };
    
  FILE *fil_pdb;

  /* Open PDB file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
    {
/* <--v.4.2 */
      if ((fil_pdb = fopen(pdb_name,"r")) ==  NULL)
        {
          printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
          exit(1);
        }
/* v.4.2--> */
    }
/* <--v.4.2 */

  /* Initialise variables */
/* v.4.3--> */
  asa_header = FALSE;
/* <--v.4.3 */
  atom_counter = 0;
  residue_counter = 0;
  residues_stored = 0;
  natoms_in_ligand = 0;
  Maximum_accessibility = 0.0;
/* v.5.4.4 --> */
  strcpy(interface_name,"No interactions across protein-protein interface");
/* <-- v.5.4.4 */

  /* Initialise pointers */
  first_atom_ptr = NULL;
  last_atom_ptr = NULL;
  first_residue_ptr = NULL;
  last_residue_ptr = NULL;
  residue_ptr = NULL;

  /* Read in the data from the PDB file */
  printf("\nReading in atom coordinates from PDB file ...\n");

  /* Initialise variables */
  asa_file = FALSE;
/* v.4.0--> */
  ligplot_pdb_file = FALSE;
/* <--v.4.0 */
  xsite_file = FALSE;
  xsite_probe[0] = '\0';
  counter = 0;
  keep_reading = TRUE;
  last_atom_no = 0;
  last_res_num[0] = '\0';
  last_res_name[0] = '\0';
  last_chain = '\0';
  inligand = FALSE;
  nconnect = 0;
  new_residue = TRUE;
  residue_in_ligand = FALSE;
  residue_in_group = FALSE;
/* v.4.0--> */
  surface1 = TRUE;
/* <--v.4.0 */
  want_atom = FALSE;

/* v.4.0--> */
  /* Check whether the file is a .asa file output by NACCESS */
  len = strlen(pdb_name);
  if (!strncmp(pdb_name+(len - 4),".asa",4))
    asa_file = TRUE;
  if (!strncmp(pdb_name,"ligplot.pdb",11))
    ligplot_pdb_file = TRUE;
/* <--v.4.0 */

/* v.4.2--> */
  /* Get the first line from the PDB or XMAS file */
  if (Read_Xmas == FALSE)
    eofstatus = fgets(line,LINELEN,fil_pdb);
  else
    {
      /* Form the record types to be retrieved */
      rec_types[REC_ATOM] = TRUE;
      rec_types[REC_CONECT] = TRUE;
      type = REC_ATOM;

      /* Get first record from the XMAS cache */
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type); */
      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
    }
/* <--v.4.2 */

  /* Loop through the PDB file, picking up all the atoms that
     are required */
/* v.4.2--> */
/*  while (fgets(line,LINELEN,fil_pdb) != NULL) */
  while (eofstatus != NULL)
/* <--v.4.2 */
    {
/* v.5.0.4--> */
      /* Truncate the string to strip off the newline */
      len = string_truncate(line,LINELEN);
/* <--v.5.0.4 */
/* v.4.0--> */
      /* Check for key words in a dimplot file */
      if (Interface_Plot == TRUE)
        {
          /* If this is the first surface, then treat as ligand */
/* v.5.4.4 --> */
//          if (!strncmp(line,"Title: ",7) &&
//              Include->Filename_for_Title == TRUE)
          if (!strncmp(line,"Title: ",7))
            {
              /* Save name of interface */
              string_ptr = strstr(line,"ligplus:");
              if (string_ptr != NULL)
                {
                  strcpy(interface_name,"No interactions between ");
                  strcat(interface_name,string_ptr + 9);
                }

              /* If using title for plot, save it */
              if (Include->Filename_for_Title == TRUE)
/* <-- v.5.4.4 */
/* v.5.0.4--> */
/*            strcpy(Print_Title,line+7); */
                {
                  strcpy(Print_Title,line+7);
                  if (strstr(Print_Title,"ligplus") != NULL)
                    strcpy(Print_Title,"ligplus");
                }
/* <--v.5.0.4 */
/* v.5.4.4 --> */
            }
/* <-- v.5.4.4 */
          else if (!strncmp(line,"SURFACE 1",9))
            surface1 = TRUE;
          else if (!strncmp(line,"SURFACE 2",9))
            surface1 = FALSE;
        }
/* <--v.4.0 */

      /* If this is a HEADER record, check whether the PDB file is
         a .asa file, output by ACCESS, containing all the atom
         accessibilities */
      if (!strncmp(line,"HEADER",6))
        {
          if (!strncmp(line+13,"ATOMIC ASAs",11))
/* v.4.3--> */
            {
/* <--v.4.3 */
              asa_file = TRUE;
/* v.4.3--> */
              asa_header = TRUE;
            }
/* <--v.4.3 */
        }
        
      /* Check for the end of an NMR structure model */
      else if (!strncmp("ENDMDL",line,6))
        keep_reading = FALSE;

      /* If this is the TITLE record in a ligplot.pdb file, then
         store it */
      else if (!strncmp("TITLE ",line,6) && Print_as_is == TRUE)
        {
          if (Include->Filename_for_Title == TRUE
              && Print_Title[0] == '\0')
            {
              /* Get the title and chop off at last non-blank character */
              ipos = 0;
              for (cpos = 0; cpos < TITLE_LEN && line[cpos + 7] != '\n';
                   cpos++)
                {
                  if (line[cpos + 7] != ' ')
                    {
                      ipos = cpos + 8;
                    }
                }
              line[ipos] = '\0';
              strcpy(Print_Title,line+7);
            }
        }
      
      /* Check whether this is an XSITE file */
      else if (!strncmp(line,"XSITE",5))
        {
          xsite_file = TRUE;
          strncpy(xsite_probe,line+5,6);
          xsite_probe[6] = '\0';
        }

/* v.5.0--> */
      /* If this is an HHB or NNB record, read in the HHB info */
/* v.5.0.1--> */
/*      else if (!strncmp(line,"HHB   ",6) || !strncmp(line,"NNB   ",6)) */
      else if (!strncmp(line,"HHB   ",6) || !strncmp(line,"NNB   ",6) ||
               !strncmp(line,"-HB   ",6) || !strncmp(line,"-NB   ",6))
/* <--v.5.0.1 */
        {
          /* Set the bond type */
/* v.5.0.1--> */
/*          if (!strncmp(line,"HHB",3)) */
          if (!strncmp(line,"HHB",3) || !strncmp(line,"-HB",3))
/* <--v.5.0.1 */
            file_type = HHB_FILE;
          else
            file_type = NNB_FILE;

/* v.5.0.1--> */
          delete = FALSE;

          /* If have a bond to be deleted, then set flag */
          if (!strncmp(line,"-HB   ",6) || !strncmp(line,"-NB   ",6))
            delete = TRUE;
/* <--v.5.0.1 */

          /* Extract the H-bond or non-bonded contact information from
             this line */
/* v.5.4--> */
//          valid = get_hbplus_info(line + 6,atom_name1,atom_name2,&chain1,
          bond_type = get_hbplus_info(line + 6,atom_name1,atom_name2,&chain1,
/* <--v.5.4 */
                                      &chain2,res_name1,res_name2,
                                      res_num1,res_num2,&bond_length,
                                      file_type,TRUE);

          /* If this is a valid bond, check whether we already have it */
/* v.5.4--> */
//          if (valid == TRUE)
          if (bond_type != INVALID)
/* <--v.5.4 */
            {
              /* Check for this bond */
/* v.5.0.1--> */
/*              found = find_hhb_info(atom_name1,atom_name2,chain1,chain2, */
              hhb_info_ptr
                = find_hhb_info(atom_name1,atom_name2,chain1,chain2,
/* <--v.5.0.1 */
                                res_name1,res_name2,res_num1,res_num2,
/* v.5.4.2--> */
//                                file_type);
                                bond_type);
/* <--v.5.4.2 */

/* v.5.0.1--> */
              /* If have found a bond that needs to be deleted, then
                 delete */
/*              if (found == TRUE) */
              if (hhb_info_ptr != NULL && delete == TRUE)
                {
                  /* Delete this record */
                  delete_hhb_info(hhb_info_ptr,&first_hhb_info_ptr);
                }

              /* If not found, then create a new bond */
/*              if (found == FALSE) */
              else if (hhb_info_ptr == NULL && delete == FALSE)
/* <--v.5.0.1 */
                {
                  /* Initialise variables */
                  atom_number1 = atom_number2 = 0;

                  /* Create a new record */
                  hhb_info_ptr
                    = create_hhb_info_record(atom_number1,atom_name1,
                                             res_name1,res_num1,chain1,
                                             atom_number2,atom_name2,
                                             res_name2,res_num2,chain2,
                                             file_type,bond_length);

/* v.5.4.1--> */
                  /* Define the bond type */
                  if (file_type == HHB_FILE)
                    hhb_info_ptr->bond_type = S_HBONDS;
                  else
                    hhb_info_ptr->bond_type = S_CONTACTS;
/* <--v.5.4.1 */
                }
            }
        }
/* <--v.5.0 */

      /* If this is an ATOM or HETATM record, process it */
      else if ((!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)) &&
               keep_reading == TRUE)
        {
          /* Get the residue name, number, chain-ID and atom type */
          strncpy(res_name,line+17,3);
          res_name[3] = '\0';
          strncpy(res_num,line+22,5);
          res_num[5] = '\0';
          chain = line[21];
          strncpy(atom_name,line+12,4);
          atom_name[4] = '\0';

          /* Store atom-number of last atom read in (for use in whittling
             out which CONECT records are to be stored */
          strncpy(atmnum,line+6,5);
          atmnum[5] = '\0';
          last_atom_no = atoi(atmnum);

          /* Initialise variables for this atom */
          ncopies = 0;

          /* Get the accessibility to update the maximum value */
          
/* v.4.0--> */
          /* If this is a ligplot.pdb file, then accessibility is beyond
             the B-value field */
          if (ligplot_pdb_file == TRUE)
            {
              strncpy(accstring,line+66,8);
              accstring[8] = '\0';
            }
/* <--v.4.0 */

          /* If this is an ordinary PDB file, store the B-value
             as the atom's accessibility */
/* v.4.0--> */
/*        if (asa_file == FALSE) */
          else if (asa_file == FALSE)
/* <--v.4.0 */
            {
              strncpy(accstring,line+60,6);
              accstring[6] = '\0';
            }
/* v.4.3--> */
          /* If no HEADER record in .asa file, then assume to be new
             NACCESS format */
          else if (asa_header == FALSE)
            {
/* v.5.1.5--> */
              /* Check if this is the short format */
              if (strlen(line) < 69)
                {
/* <--v.5.1.5 */
                  strncpy(accstring,line+55,7);
/* v.5.1.5--> */
                }

              /* Otherwise, take to be the "full" format, generated by
                 running naccess with the -f flag */
              else
                {
                  strncpy(accstring,line+65,7);
                }
/* <--v.5.1.5 */
              accstring[7] = '\0';
            }
/* <--v.4.3 */
            
          /* Otherwise, if this is a .asa file output by Simon Hubbard's
             NACCESS program, then pick up the accessibility from the
             appropriate location */
          else
            {
              strncpy(accstring,line+65,7);
              accstring[7] = '\0';
            }
          accessibility = (float) atof(accstring);
          if (accessibility > Maximum_accessibility)
            Maximum_accessibility = accessibility;

          /* Check whether atom belongs to a new residue */
          if (chain != last_chain || strncmp(res_num,last_res_num,5) ||
              strncmp(res_name,last_res_name,3))
            {
              /* If residue is a new one, need to check whether it
                 either belongs to the ligand or to one of the
                 residues H-bonded to the ligand */
              new_residue = TRUE;
              residue_in_ligand = FALSE;
              residue_in_group = FALSE;

              /* Check whether the residue falls within the range
                 defining the ligand */
              inligand = TRUE;

/* v.4.0--> */
              /* For interface plot, check whether this residue belongs
                 to the first surface */
              if (Interface_Plot == TRUE)
                {
                  inligand = surface1;

                  /* Check whether residue from the other surface
                     or water */
                  if (surface1 == FALSE)
                    residue_in_group = TRUE;
                }

              /* Otherwise, check whether residue falls within the
                 sequential range defining the ligand */
              else
/* <--v.4.0 */

                inligand = check_if_in_ligand(chain,res_name,res_num,
/* v.4.2.2--> */
/*                                            FALSE,dummy); */
                                              FALSE,dummy,mol_id);
/* <--v.4.2.2 */

              /* If the residue is a ligand residue, then want to
                 store all its atom coords */
              if (inligand == TRUE)
                residue_in_ligand = TRUE;
                  
              /* If residue is not one of the ligand residues, check
                 whether it is one of the residues H-bonded to the
                 ligand, or involved in non-bonded contacts */
/* v.4.0--> */
/*            else */
              else if (Interface_Plot == FALSE)
/* <--v.4.0 */
                {
                  inligand = FALSE;

                  /* If printing structure as it stands (ie reading
                     in from the ligplot.pdb file), then want
                     to read in all the atom records present */ 
                  if (Print_as_is == TRUE)
                    residue_in_group = TRUE;

                  /* Otherwise, check the HBPLUS information */
                  else
                    {
                      /* Set pointer to the first of the links */
                      hhb_info_ptr = first_hhb_info_ptr;

                      /* Search through all the stored links */
                      while (hhb_info_ptr != NULL)
                        {
                          /* Initialise flag */
                          wanted = FALSE;

                          /* Check if this residue is in the current
                             interaction */
                          if (chain == hhb_info_ptr->chain1 &&
                              !strncmp(res_name,hhb_info_ptr->res_name1,3) &&
                              !strncmp(res_num,hhb_info_ptr->res_num1,5))
                            wanted = TRUE;
                          else if (chain == hhb_info_ptr->chain2 &&
                                   !strncmp(res_name,
                                            hhb_info_ptr->res_name2,3) &&
                                   !strncmp(res_num,
                                            hhb_info_ptr->res_num2,5))
                            wanted = TRUE;

/* v.4.0--> */
                          /* If this is a covalent bond, then check
                             that is is between the ligand and a non-ligand
                             residue */
                          if (wanted == TRUE &&
                              hhb_info_ptr->source == CONECT)
                            {
                              /* Check if first residue is in the ligand */
                              in_lig1
                                = check_if_in_ligand(hhb_info_ptr->chain1,
                                                     hhb_info_ptr->res_name1,
                                                     hhb_info_ptr->res_num1,
/* v.4.2.2--> */
/*                                                   FALSE,dummy); */
                                                     FALSE,dummy,mol_id);
/* <--v.4.2.2 */

                              /* Check if first residue is in the ligand */
                              in_lig2
                                = check_if_in_ligand(hhb_info_ptr->chain2,
                                                     hhb_info_ptr->res_name2,
                                                     hhb_info_ptr->res_num2,
/* v.4.2.2--> */
/*                                                   FALSE,dummy); */
                                                     FALSE,dummy,mol_id);
/* <--v.4.2.2 */

                              /* If neither residue is in the ligand,
                                 discount this bond */
                              if (in_lig1 == FALSE && in_lig2 == FALSE)
                                wanted = FALSE;
                            }
/* <--v.4.0 */

                          /* If residue is included, then need to store
                             all its coordinates */
                          if (wanted == TRUE)
                            residue_in_group = TRUE;

                          /* Get the pointer to the next link */
                          hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
                        }
                    }
                }

              /* Save the residue details */
              last_chain = chain;
              strncpy(last_res_num,res_num,5);
              last_res_num[5] = '\0';
              strncpy(last_res_name,res_name,3);
              last_res_name[3] = '\0';
              
              /* Increment residue-count */
              residue_counter++;
            }

          /* Determine whether current atom is wanted for a ligand
             or H-group */
          if (residue_in_ligand == TRUE || residue_in_group == TRUE)
            want_atom = TRUE;
          else
            want_atom = FALSE;

          if (residue_in_ligand == TRUE)
            inligand = TRUE;

/* v.4.3--> */
          /* Get the atom type and number */
          strncpy(atom_name,line+12,4);
          atom_name[4] = '\0';

          /* Get the length of the input line */
          len = strlen(line);

          /* Get the element identifier, if there is one */
          if (len > 77)
            {
              strncpy(element,line+76,2);
              element[2] = '\0';
            }
          else
            strcpy(element,"  ");

          /* Check whether this is a valid element */
          if ((element[0] >= '0' && element[0] <= '9') ||
              (element[1] >= '0' && element[1] <= '9'))
            strcpy(element,"  ");

          /* If don't have the element, have to extract it from
             the atom name */
          if (!strcmp(element,"  "))
            get_element(atom_name,element);
/* <--v.4.3 */

          /* NMR structures in the pdb have hydrogens and pseudo-atoms, 
             so we need to remove these */
/* v.4.3--> */
/*        hydrogen = is_hydrogen_atom(line+12); */
          if (!strcmp(element," H"))
            hydrogen = TRUE;
          else
            hydrogen = FALSE;
/* <--v.4.3 */
          if (hydrogen == TRUE && Include_Hydrogens != TRUE)
            want_atom = FALSE;

          /* If this is an alternate-occupancy atom, check whether
             we already have it */
          strncpy(occup,line+56,3);
          occup[3] = '\0';
          if (want_atom == TRUE && (line[16] != ' ' || strncmp(occup,"1.0",3)))
            {
              /* Initialise pointer to start of linked list of atom coords */
              other_atom_ptr = first_atom_ptr;

              /* Loop through all the atoms to check whether we
                 already have this one */
              while (other_atom_ptr != NULL && want_atom == TRUE)
                {
                  /* Get pointer to the atom's residue */
                  other_residue_ptr = other_atom_ptr->residue_ptr;

                  /* If we already have this atom, then don't want
                     the second copy */
                  if (!strncmp(other_residue_ptr->res_name,res_name,3) &&
                      !strncmp(other_residue_ptr->res_num,res_num,5) &&
                      other_residue_ptr->chain == chain &&
                      !strncmp(other_atom_ptr->atom_name,atom_name,4))
                    want_atom = FALSE;

                  /* Get pointer to next atom in linked-list */
                  other_atom_ptr = other_atom_ptr->next;
                }              
            }

          /* If the current atom is required, save its details */
          if (want_atom == TRUE)
            {
              /* Get the atom's coordinates */
              sscanf(line+30," %f %f %f ",&x,&y,&z);

              /* Store the atom number */
              strncpy(atom_number,line+6,5);
              atom_number[5] = '\0';

              /* Create new atom record and store this atom's
                 details */
/* v.4.3--> */
              atom_ptr = create_atom_record(&last_atom_ptr,atom_name,
                                            atom_number,x,y,z);
/* <--v.4.3 */

              /* Get occupancy, B-balue and accessibility, where
                 present */
/* v.4.3--> */
              if (len > 60)
              {
/* <--v.4.3 */
                strncpy(atom_ptr->occupancy,line+54,6);
                atom_ptr->occupancy[6] = '\0';
/* v.4.3--> */
              }
              if (len > 66)
              {
/* <--v.4.3 */
              strncpy(atom_ptr->bvalue,line+60,6);
              atom_ptr->bvalue[6] = '\0';
/* v.4.3--> */
              }
/* <--v.4.3 */
              atom_ptr->accessibility = accessibility;

              /* If atom in ligand, increment count of ligand atoms */
              if (inligand == TRUE)
                natoms_in_ligand++;
                
              /* Check that this is a standard amino acid */
              found = FALSE;
              for (iamino = 0; iamino < 20 && found == FALSE; iamino++)
                if (!strncmp(res_name,amino[iamino],3))
                  found = TRUE;

              /* If this is not a standard amino acid, then
                 mark atom as a main-chain atom */
              if (Include->Simple_Ligand_Residues == TRUE && found == FALSE)
                atom_ptr->side_chain = FALSE;

/* v.4.3--> */
              /* Store element */
              strcpy(atom_ptr->element,element);

              /* If printing in chemical notation, replace atom name
                 by element */
              if (Include->Chemical_Notation == TRUE)
                {
                  /* If carbon, then no name required */
                  if (!strcmp(element," C"))
                    atom_ptr->print_name[0] = '\0';

                  /* Otherwise, use element name */
                  else
                    strcpy(atom_ptr->print_name,element);
                }
/* <--v.4.3 */

              /* If this is an XSITE file, and the print-field is
                 asterisked, read in the character-string to be
                 printed in place of the atom name on the final plot */
              if (xsite_file == TRUE && line[67] == '*')
                {
                  strncpy(atom_ptr->print_name,line+68,4);
                  atom_ptr->print_name[4] = '\0';
                }

/* v.4.3--> */
              /* If printing in chemical notation, then show only the
                 element name for non-carbon atoms */
              if (Include->Chemical_Notation == TRUE)
                {
                  /* If this is a metal atom, get element name */
                }
/* <--v.4.3 */

              /* Increment counts */
              counter++;
              if (atom_counter >= MAXATOMS)
                {
                  printf("*** Maximum number of atoms");
                  printf(" (%d) exceeded. Program aborted.\n",
                         MAXATOMS);
/* v.5.4.4 --> */
                  /* Write message out to the ligplot.drw file */
                  fprintf(ligplot_drw,
                          "ERROR: Exceeded maximum number of atoms "
                          "on plot (%d)\n",MAXATOMS);

                  /* Close the .drw file */
                  fclose(ligplot_drw);
/* <-- v.5.4.4 */
                  exit(1);
                }
              ncopies++;
              atom_counter++;

              /* If this is a new residue, create a residue record,
                 store residue data, and increment the residue-count */
              if (new_residue == TRUE)
                {
                /* Create new residue record and store this residue's
                     details */
/* v.4.3--> */
                  residue_ptr
                    = create_residue_record(&last_residue_ptr,res_name,
                                            res_num,chain);
/* <--v.4.3 */

                  /* Store the residue details */
                  if (inligand == TRUE &&
                      Include->Simple_Ligand_Residues == TRUE)
                    residue_ptr->residue_type = SIMPLE_LIGAND;
                  if (inligand == FALSE &&
                      Include->Simple_Nonligand_Residues == TRUE)
                    residue_ptr->residue_type = SIMPLE_HGROUP;
                  residue_ptr->inligand = inligand;
                  residue_ptr->first_atom_ptr = atom_ptr;

                  /* Increment count of residues */
                  residues_stored++;
                  new_residue = FALSE;
                }

              /* Increment count of atoms in this residue */
              residue_ptr->natoms++;

              /* Get atom record to point to its residue record */
              atom_ptr->residue_ptr = residue_ptr;
            }
        }
      
      /* If this is a CONECT record, store the covalent bonds it reports */
      else if (!strncmp(line,"CONECT",6))
        {
          /* Get the atom-number of the first number in the CONECT record */
          strncpy(atmnum,line+6,5);
          atmnum[5] = '\0';
          atom_no = atoi(atmnum);

          /* If the atom number is within the range of atoms read in,
             then process */
          if (atom_no <= last_atom_no)
            {
              /* Search through all the stored atoms for a match with the
                 first atom-number in the CONECT record */
              found = FALSE;

              /* Initialise pointer to start of linked list of atom coords */
              atom1_ptr = first_atom_ptr;

              /* Loop through all the atoms in the linked list */
              while (atom1_ptr != NULL && found == FALSE)
                {
                  /* Check if this atom's number matches */
                  if (!strncmp(atom1_ptr->atom_number,line+6,5))
                    {
                      /* If match found, store the atom's coordinates */
                      x1 = atom1_ptr->x;
                      y1 = atom1_ptr->y;
                      z1 = atom1_ptr->z;
                      found = TRUE;
                    }

                  /* If have found a match for the first atom, locate all
                     the other atoms it is reported as being bonded to */
                  if (found == TRUE)
                    {
                      /* Get the next atom-number */
                      ipos = 0;
                      strncpy(atom_partner,line+(5 * ipos + 11),5);
/* v.4.1.2--> */
                      if (line[5 * ipos + 11] == '\0' ||
                          line[5 * ipos + 11] == '\n')
                        strncpy(atom_partner,"     ",5);
/* <--v.4.1.2 */
                      atom_partner[5] = '\0';
                  
                      /* Loop through all possible bonding partners */
                      while (strncmp(atom_partner,"     ",5) && ipos < 4)
                        {
                          /* Search for this atom-partner in the atom
                             records stored */
                          found = FALSE;

                          /* Initialise pointer to start of linked list
                             of atom coords */
                          atom2_ptr = first_atom_ptr;

                          /* Loop through all the atoms in the linked list */
                          while (atom2_ptr != NULL && found == FALSE)
                            {
                              /* Initialise flag indicating whether bond
                                 is required */
                              wanted = FALSE;

                              /* Check if this atom's number matches */
                              if (!strncmp(atom2_ptr->atom_number,
                                           atom_partner,5))
                                {
                                  found = TRUE;

                                  /* Want this covalent bond only if it
                                     belongs entirely to the ligand (any
                                     ligand-nonligand bonds will already have
                                     been picked up in a previous pass) */
                                  if (atom1_ptr->residue_ptr->inligand
                                      == TRUE &&
                                      atom2_ptr->residue_ptr->inligand == TRUE)
                                    wanted = TRUE;
                                }

                              /* If desired match found, store the atom's 
                                 coordinates */
                              if (wanted == TRUE)
                                {
                                  x2 = atom2_ptr->x;
                                  y2 = atom2_ptr->y;
                                  z2 = atom2_ptr->z;

                                  /* Check the bond-length */
                                  length2 = (x2 - x1) * (x2 - x1)
                                    + (y2 - y1) * (y2 - y1)
                                      + (z2 - z1) * (z2 - z1);
/* v.5.0.1--> */
/*                                  if (length2 > 9.0) */
                                  if (length2 > 9.0 && conect == FALSE)
/* <--v.5.0.1 */
                                    {
                                      printf("*** Warning. Long bond from ");
                                      printf("CONECT [%5s -> %5s] %7.2f",
                                             atom1_ptr->atom_number,
                                             atom_partner,
                                             sqrt((double) length2));
                                      if (length2 > 16.0)
                                        {
                                          printf(" -- discarded\n");
                                          found = FALSE;
                                        }
                                      else
                                        printf("\n");
                                      Nwarnings++;
                                    }
                              
                                  /* Store this covalent bond */
                                  if (found == TRUE)
                                    {
                                      length = (float) sqrt((double) length2);
                                      update_bond(atom1_ptr,atom2_ptr,
                                                  COVALENT,CONECT,length);
                                      nconnect++;
                                    }
                                }

                              /* Get pointer to next atom in linked-list */
                              atom2_ptr = atom2_ptr->next;
                            }

                          /* Get the next atom-number */
                          ipos++;
                          strncpy(atom_partner,line+(5 * ipos + 11),5);
                          atom_partner[5] = '\0';
                        }
                    }

                  /* Get pointer to next atom in linked-list */
                  atom1_ptr = atom1_ptr->next;
                }
            }
        }

/* v.4.0--> */
      /* If this is an ENDLIN record in a ligplot.pdb file, pick up the
         maximum accessibility value */
      else if (ligplot_pdb_file == TRUE && !strncmp(line,"ENDLIN",6))
        {
          /* Get the accessibility value */
          strncpy(accstring,line+66,8);
          accstring[8] = '\0';
          
          /* Save as the maximum value, if indeed higher than current
             maximum */
          accessibility = (float) atof(accstring);
          if (accessibility > Maximum_accessibility)
            Maximum_accessibility = accessibility;
        }
/* <--v.4.0 */

/* v.4.2--> */
      /* Get the next line from the PDB or XMAS file */
      if (Read_Xmas == FALSE)
        eofstatus = fgets(line,LINELEN,fil_pdb);
      else
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type); */
        eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type,
                                      &mol_id);
/* <--v.4.2.2 */
/* <--v.4.2 */
    }
    
  /* Print counts of data read in */
  printf("   Number of residues stored         = %7d\n",residues_stored);
  printf("   Number of atoms stored            = %7d\n",atom_counter);
  printf("   Number of ligand atoms            = %7d\n",natoms_in_ligand);
  printf("   Number of CONECT bonds stored     = %7d\n",nconnect);
  if (atom_counter == 0)
    {
      printf("\n*** No atoms found. Nothing to plot\n");
/* v.5.4.4 --> */
/* <--v.4.0 */
      /* Write message out to the ligplot.drw file */
      if (Interface_Plot == TRUE)
        fprintf(ligplot_drw,"ERROR: %s\n",interface_name);
      else
        fprintf(ligplot_drw,"ERROR: No atoms found\n");

      /* Close the .drw file */
      fclose(ligplot_drw);
/* <-- v.5.4.4 */
      exit(1);
    }
  if (natoms_in_ligand == 0)
    {
      printf("\n*** No ligand atoms found. Nothing to plot\n");
/* v.5.4.4 --> */
      /* Write message out to the ligplot.drw file */
      fprintf(ligplot_drw,"ERROR: No ligand atoms found\n");

      /* Close the .drw file */
      fclose(ligplot_drw);
/* <-- v.5.4.4 */
      exit(1);
    }

  /* Close the PDB file */
/* v.4.2--> */
  if (Read_Xmas == FALSE)
/* <--v.4.2 */
    fclose(fil_pdb);

/* v.4.0--> */
  /* Check that maximum accessibility is non-zero */
  if (Maximum_accessibility == 0.0)
    Maximum_accessibility = 100.0;
/* <--v.4.0 */

  /* Return number of atoms read in */
  return(atom_counter);
}
/***********************************************************************

define_objects  -  Define the objects that will appear on the final
                   LIGPLOT diagram

***********************************************************************/

void define_objects(int *nobjects,int *nligand_objects)
{
  int iatom, natoms;
  int got_ligand_object;
/* v.4.0--> */
  int new_object;
/* <--v.4.0 */

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;
  struct object *object_ptr, *last_object_ptr;

  /* Initialise variables */
  got_ligand_object = FALSE;
  first_object_ptr = NULL;
  last_object_ptr = NULL;
  *nobjects = 0;
  *nligand_objects = 0;

  /* Initialise pointer to the first stored residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all residues */
  while (residue_ptr != NULL)
    {
/* v.4.0--> */
      /* Determine whether a new object is to be created */
      new_object = FALSE;
      
      /* If not a ligand residue, or the first of the ligand residues,
         need to create a new object record */
      if (residue_ptr->inligand == FALSE || got_ligand_object == FALSE)
        new_object = TRUE;

      /* If plotting a DIMPLOT from a flattened ligplot.pdb file, then
         each residue needs to be a separate object */
      if (Interface_Plot == TRUE && Print_as_is == TRUE)
        new_object = TRUE;
/* <--v.4.0 */

/* v.4.0--> */
      /* If this is a new object, create a new object record */
/*      if (residue_ptr->inligand == FALSE || got_ligand_object == FALSE) */
      if (new_object == TRUE)
/* <--v.4.0 */
        {
          /* Store this object's details */

/* v.4.3--> */
          /* Create an entry for this object */
          object_ptr = create_object_record(&last_object_ptr);
/* <--v.4.3 */

          /* Store pointer to this, the first, residue */
          object_ptr->first_residue_ptr = residue_ptr;

          /* Define whether object is the ligand or a peripheral
             residue */
          if (residue_ptr->inligand == TRUE)
            {
              object_ptr->object_type = LIGAND;
              got_ligand_object = TRUE;
              (*nligand_objects)++;
            }
          else if (!strncmp(residue_ptr->res_name,"HOH",3) ||
                   !strncmp(residue_ptr->res_name,"WAT",3))
            object_ptr->object_type = WATER;
/* v.5.0--> */
          else if (!strncmp(residue_ptr->res_name,"...",3))
            object_ptr->object_type = DUMMY;
/* <--v.5.0 */
          else
            object_ptr->object_type = HGROUP;

          /* If this is a non-ligand residue, reset the got-ligand
             flag */
          if (residue_ptr->inligand == FALSE)
            got_ligand_object = FALSE;

          /* Increment count of objects */
          (*nobjects)++;
        }

      /* Increment count of residues in this object */
      object_ptr->nresidues++;

      /* Store pointer from residue to object containing it */
      residue_ptr->object_ptr = object_ptr;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
              
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, counting number of
         mainchain atoms */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* If this is a water, then mark as non-mainchain */
          if (object_ptr->object_type == WATER)
            atom_ptr->side_chain = TRUE;

          /* If not a sidechain atom, increment count of this object's
             mainchain atoms */
          if (atom_ptr->side_chain == FALSE)
            object_ptr->nmain_chain++;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Print number of objects created */
  printf("   Number of objects                 = %7d\n",*nobjects);
}
/***********************************************************************

update_residue_boundaries  -  Calculate the maximum and minimum
                              coordinates of the given object's residues

***********************************************************************/

void update_residue_boundaries(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  int first_residue;
/* v.4.0--> */
  int ncoords;
  float flattened_mean_x, flattened_mean_y;
/* <--v.4.0 */

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  first_residue = TRUE;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
/* v.4.0--> */
      /* Initialise mean residue coordinates */
      flattened_mean_x = 0.0;
      flattened_mean_y = 0.0;
      ncoords = 0;
/* <--v.4.0 */

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Initialise residue limits to the coords of this first atom,
         if there is one(!) */
      if (natoms > 0)
        {
          residue_ptr->minx = atom_ptr->x;
          residue_ptr->maxx = atom_ptr->x;
          residue_ptr->miny = atom_ptr->y;
          residue_ptr->maxy = atom_ptr->y;
        }

      /* Loop over all the residue's atoms, checking the x- and
         y-coordinates */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* If this is a hydrophobic group, or simple residue, need
             to fix all the atom-coords to be equal */
          if (object_ptr->object_type == HYDROPHOBIC ||
              object_ptr->object_type == SIMPLE_HGROUP)
            {
              atom_ptr->x = residue_ptr->minx;
              atom_ptr->y = residue_ptr->miny;
/* v.4.0--> */
              residue_ptr->flattened_mean_x = residue_ptr->minx;
              residue_ptr->flattened_mean_y = residue_ptr->miny;
/* <--v.4.0 */
            }

          /* Otherwise check this atom's coordinates against the
             maximum and minimum values */
          else
            {
              if (atom_ptr->x < residue_ptr->minx)
                residue_ptr->minx = atom_ptr->x;
              if (atom_ptr->x > residue_ptr->maxx)
                residue_ptr->maxx = atom_ptr->x;
              if (atom_ptr->y < residue_ptr->miny)
                residue_ptr->miny = atom_ptr->y;
              if (atom_ptr->y > residue_ptr->maxy)
                residue_ptr->maxy = atom_ptr->y;

/* v.4.0--> */
              /* Update mean coords of the current x-y coords of the
                 atom in the residue's flattened state */
              flattened_mean_x = flattened_mean_x + atom_ptr->x;
              flattened_mean_y = flattened_mean_y + atom_ptr->y;
              ncoords++;
/* <--v.4.0 */
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* If this is the object's first residue, set its boundaries
         equal to that of the residue */
      if (natoms > 0)
        {
          if (first_residue == TRUE)
            {
              object_ptr->minx = residue_ptr->minx;
              object_ptr->maxx = residue_ptr->maxx;
              object_ptr->miny = residue_ptr->miny;
              object_ptr->maxy = residue_ptr->maxy;
              first_residue = FALSE;
            }

          /* Otherwise, adjust if necessary */
          else
            {
              if (residue_ptr->minx < object_ptr->minx)
                object_ptr->minx = residue_ptr->minx;
              if (residue_ptr->maxx > object_ptr->maxx)
                object_ptr->maxx = residue_ptr->maxx;
              if (residue_ptr->miny < object_ptr->miny)
                object_ptr->miny = residue_ptr->miny;
              if (residue_ptr->maxy > object_ptr->maxy)
                object_ptr->maxy = residue_ptr->maxy;
            }
        }

/* v.4.0--> */
      /* Get this residue's mean coordinates */
      if (ncoords > 0)
        {
          /* Calculate the flattened-out residue's mean position, and
             store the coords */
          flattened_mean_x = flattened_mean_x / ncoords;
          flattened_mean_y = flattened_mean_y / ncoords;

          /* Store the object's original and flattened mean coords */
          residue_ptr->flattened_mean_x = flattened_mean_x;
          residue_ptr->flattened_mean_y = flattened_mean_y;
        }
/* <--v.4.0 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

update_all_boundaries  -  Update all object and residue boundaries

***********************************************************************/

void update_all_boundaries(void)
{
  struct object *object_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Update this object's boundaries */
      update_residue_boundaries(object_ptr);

      /* Get pointer to the next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* v.5.0--> */
/***********************************************************************

compare_anchor_dists  -  Compare the current anchor position and
                         corresponding coords against all other anchor
                         positions

***********************************************************************/

float compare_anchor_dists(float x1,float y1,float anchor1_x,
                           float anchor1_y,struct residue *residue_ptr,
                           int iresid,int nresid,int start_atom,
                           float anchor_energy)
{
  int iatom, natoms;

  float anchor2_x, anchor2_y, x2, y2;
  float anchor_dist, anchor_dist2, atom_dist, atom_dist2, dist2;

  struct coordinate *atom_ptr;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* If any of the residue's atoms have anchor positions, loop through
         these */
      if (residue_ptr->nanchored_atoms > 0)
        {
          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all this residue's atoms */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* If this is an anchored atom, process it */
              if (atom_ptr->have_anchor == TRUE &&
                  (iatom > start_atom || start_atom == -1))
                {
                  /* Get the atom's anchor coordinates */
                  anchor2_x = atom_ptr->anchor_pstn_x;
                  anchor2_y = atom_ptr->anchor_pstn_y;

                  /* Store the atom's coordinates */
                  x2 = atom_ptr->x;
                  y2 = atom_ptr->y;

                  /* Calculate distance between the two anchor positions */
                  anchor_dist2
                    = (anchor1_x - anchor2_x) * (anchor1_x - anchor2_x)
                    + (anchor1_y - anchor2_y) * (anchor1_y - anchor2_y);
                  anchor_dist = sqrt(anchor_dist2);

                  /* Calculate distance between the corresponding
                     atom positions */
                  atom_dist2
                    = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
                  atom_dist = sqrt(atom_dist2);

                  /* Compare the two distances */
                  dist2 = (anchor_dist - atom_dist)
                    * (anchor_dist - atom_dist);

                  /* Calculate the equivalent energy */
                  anchor_energy = anchor_energy + Anchor_Weight * dist2;
                }

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;

      /* Re-initialise start atom */
      start_atom = -1;
    }

  /* Return this object's energy */
  return(anchor_energy);
}
/***********************************************************************

get_relative_anchor_energy  -  Calculate the energy of the atoms with
                               anchor positions, by comparing their
                               relative distances from one another with
                               the relative distances of the
                               corresponding anchor positions

***********************************************************************/

float get_relative_anchor_energy(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;

  float anchor_x, anchor_y, x, y;
  float anchor_energy;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  anchor_energy = 0.0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* If any of the residue's atoms have anchor positions, loop through
         these */
      if (residue_ptr->nanchored_atoms > 0)
        {
          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all this residue's atoms */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* If this is an anchored atom, process it */
              if (atom_ptr->have_anchor == TRUE)
                {
                  /* Get the atom's anchor coordinates */
                  anchor_x = atom_ptr->anchor_pstn_x;
                  anchor_y = atom_ptr->anchor_pstn_y;

                  /* Store the atom's coordinates */
                  x = atom_ptr->x;
                  y = atom_ptr->y;

                  /* Compare against all other anchor positions */
                  anchor_energy
                    = compare_anchor_dists(x,y,anchor_x,anchor_y,
                                           residue_ptr,iresid,nresid,
                                           iatom,anchor_energy);
                }

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return this object's energy */
  return(anchor_energy);
}
/* <--v.5.0 */
/* v.4.0--> */
/***********************************************************************

get_anchor_energy  -  Calculate the given object's distance from its
                      anchor-position and covert into an anchor energy

***********************************************************************/

float get_anchor_energy(struct object *object_ptr)
{
  int anchors, iatom, iresid, nanchored, natoms, nresid;

  float anchor_x, anchor_y, x, y;
  float dist2;
  float anchor_energy;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  anchor_energy = 0.0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Process only if this residue has an anchor position */
      if (residue_ptr->have_anchor == TRUE)
        {
          /* Get the residue's anchor coordinates */
          anchor_x = residue_ptr->anchor_pstn_x;
          anchor_y = residue_ptr->anchor_pstn_y;

          /* Store the residue's coordinates */
          x = residue_ptr->flattened_mean_x;
          y = residue_ptr->flattened_mean_y;

          /* Calculate distance between the residue's current mean
             position and its anchoring position */
          dist2 = (x - anchor_x) * (x - anchor_x)
            + (y - anchor_y) * (y - anchor_y);

          /* Calculate the equivalent energy */
          anchor_energy = anchor_energy + Anchor_Weight * dist2;
        }

      /* If any of the residue's atoms have anchor positions, loop through
         these */
      if (residue_ptr->nanchored_atoms > 0)
        {
          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;
          nanchored = 0;
          anchors = residue_ptr->nanchored_atoms;

          /* Loop over all this residue's atoms */
          while (iatom < natoms && atom_ptr != NULL && nanchored < anchors)
            {
              /* If this is an anchored atom, process it */
              if (atom_ptr->have_anchor == TRUE)
                {
                  /* Get the atom's anchor coordinates */
                  anchor_x = atom_ptr->anchor_pstn_x;
                  anchor_y = atom_ptr->anchor_pstn_y;

                  /* Store the atom's coordinates */
                  x = atom_ptr->x;
                  y = atom_ptr->y;

                  /* Calculate distance between the atom's current mean
                     position and its anchoring position */
                  dist2 = (x - anchor_x) * (x - anchor_x)
                    + (y - anchor_y) * (y - anchor_y);

                  /* Calculate the equivalent energy */
                  anchor_energy
/* v.5.0--> */
/*                    = anchor_energy + Anchor_Weight * dist2; */
                    = anchor_energy
                    + atom_ptr->anchor_weight * Anchor_Weight * dist2;
/* <--v.5.0 */

                  /* Increment count of anchors encountered */
                  nanchored++;
                }

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return this object's energy */
  return(anchor_energy);
}
/***********************************************************************

get_rel_pstn_energy  -  Calculate the relative distances between the
                        given object's residues and all other objects'
                        residues to get an energy from a comparison of the
                        distances in the current flattened layout with
                        the distances in the original coordinates

***********************************************************************/

float get_rel_pstn_energy(struct object *object_ptr)
{
  int iresid, jresid, nother_resid, nresid;

  float other_x, other_y, other_z, x, y, z;
  float other_flattened_x, other_flattened_y, flattened_x, flattened_y;
  float dist2, original_dist2;
  float energy, rel_pstn_energy, weight;

  struct object *other_object_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Initialise variables */
  rel_pstn_energy = 0.0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Store the residue's coordinates */
      x = residue_ptr->original_mean_x;
      y = residue_ptr->original_mean_y;
      z = residue_ptr->original_mean_z;
      flattened_x = residue_ptr->flattened_mean_x;
      flattened_y = residue_ptr->flattened_mean_y;

      /* Initialise pointer to the first stored object */
      other_object_ptr = first_object_ptr;

      /* Loop through all objects */
      while (other_object_ptr != NULL)
        {
          /* Process if not the same as the object of interest */
          if (other_object_ptr != object_ptr)
            {
              /* Set the pointer to the first of this object's residues */
              other_residue_ptr = other_object_ptr->first_residue_ptr;
              nother_resid = other_object_ptr->nresidues;
              jresid = 0;

              /* Loop over all this object's residues */
              while (jresid < nother_resid && other_residue_ptr != NULL)
                {
                  /* Get this residue's original and current mean coords */
                  other_x = other_residue_ptr->original_mean_x;
                  other_y = other_residue_ptr->original_mean_y;
                  other_z = other_residue_ptr->original_mean_z;
                  other_flattened_x = other_residue_ptr->flattened_mean_x;
                  other_flattened_y = other_residue_ptr->flattened_mean_y;

                  /* Calculate distance between these two residues in the
                     original structure */
                  original_dist2 = (x - other_x) * (x - other_x)
                    + (y - other_y) * (y - other_y)
                      + (z - other_z) * (z - other_z);

                  /* Calculate the weighting factor (such that nearby residues
                     in the original structure maintain their proximity better
                     than distant residues) */
                  if (original_dist2 > 0.01)
                    weight = Rel_Pstn_Energy_Weight / original_dist2;
                  else
                    weight = 0.0;

                  /* Repeat for distance between them for the current
                     flattened objects */
                  dist2 = (flattened_x - other_flattened_x)
                    * (flattened_x - other_flattened_x)
                      + (flattened_y - other_flattened_y)
                        * (flattened_y - other_flattened_y);

                  /* Calculate the difference, and hence the effective
                     energy */
                  energy = (float) (weight * fabs(original_dist2 - dist2));
                  rel_pstn_energy = rel_pstn_energy + energy;

                  /* Get pointer to the next residue */
                  other_residue_ptr = other_residue_ptr->next_residue_ptr;
                  jresid++;
                }
            }

          /* Get pointer to the next object */
          other_object_ptr = other_object_ptr->next_object_ptr;
        }

      /* Get pointer to the next residue of the current object */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return this object's energy */
  return(rel_pstn_energy);
}
/***********************************************************************

get_boundary_energy  -  Calculate the given object's boundary energy
                        according to how far its atoms have strayed into
                        the wrong side of the plot (only relevant for
                        interface plots when objects are placed on either
                        side of a boundary, with any waters lying on,
                        or close to, the boundary)

***********************************************************************/

float get_boundary_energy(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid, object_type;
  int interface;
/* v.5.3.4--> */
  int first;
/* <--v.5.3.4 */

/* v.5.3.4--> */
//  float x;
  float x, xdist, xmin, xmax, y, ydist, ymin, ymax;
/* <--v.5.3.4 */
  float dist, object_size;
  float boundary_energy, energy;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  boundary_energy = 0.0;
/* v.5.3.4--> */
  first = TRUE;
  xmin = xmax = ymin = ymax = 0;
/* <--v.5.3.4 */

  /* Get this object type */
  object_type = object_ptr->object_type;
  interface = object_ptr->interface;

  /* Get the object size for calculating distance from boundary */
  object_size = 0.0;
  if (object_type == HYDROPHOBIC)
    object_size = Size_Val->Hydrophobics;
  else if (object_type == WATER)
    object_size = Size_Val->Waters;
  else if (object_type == LIGAND)
    object_size = Size_Val->Ligand_Atoms;
  else if (object_type == HGROUP)
    object_size = Size_Val->Nonligand_Atoms;
  else if (object_type == SIMPLE_HGROUP)
    object_size = Size_Val->Simple_Residues;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, checking whether
         they are involved in interactions */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Get this atom's x- and y-coordinates */
          x = atom_ptr->x;
/* v.5.3.4--> */
          y = atom_ptr->y;

          /* If the first, save as the maximum and minum coords */
          if (first == TRUE)
            {
              xmin = xmax = x;
              ymin = ymax = y;
              first = FALSE;
            }

          /* Otherwise, update maximum and minimum values */
          else
            {
              if (x < xmin)
                xmin = x;
              if (x > xmax)
                xmax = x;
              if (y < ymin)
                ymin = y;
              if (y > ymax)
                ymax = y;
            }
/* <--v.5.3.4 */

          /* Initialise "bad" distance from boundary */
          dist = 0.0;

          /* If object is a water molecule, want it fairly close to the
             boundary */
          if (object_type == WATER)
             dist = (float) (fabs(x));

          /* If residue belongs to first surface, want all its
             atoms on the left-hand side of the boundary */
          else if (interface == 1 && x > - (Min_Boundary_Dist + object_size))
            dist = x - (Min_Boundary_Dist + object_size);

          /* If residue belongs to second surface, want all its
             atoms on the right-hand side of the boundary */
          else if (interface == 2 && x < (Min_Boundary_Dist + object_size))
            dist = (Min_Boundary_Dist + object_size) - x;

          /* If distance is non-zero, calculate a corresponding
             boundary-violation energy for this atom */
          if (dist > 0.0)
            {
              energy = Boundary_Energy_Weight * dist * dist;

              /* Add to total energy for this object */
              boundary_energy = boundary_energy + energy;
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

/* v.5.3.4--> */
  /* Calculate whether object is more vertical than horizontal */
  xdist = xmax - xmin;
  ydist = ymax - ymin;
//  if (xdist > 0)
//    factor = ydist / xdist;
//  else
//    factor = MAX_FACTOR;
/* <--v.5.3.4 */

  /* Return this object's energy */
  return(boundary_energy);
}
/* <--v.4.0 */
/***********************************************************************

open_pdb_output  -  Open output PDB files, ligplot.pdb, ligplot.frm and
                    ligplot.res

***********************************************************************/

/* v.4.0--> */
/* void open_pdb_output(char res_num1[6],char res_num2[6],char chain) */
void open_pdb_output(char res_name1[4],char res_num1[6],char res_name2[4],
                     char res_num2[6],char chain)
/* <--v.4.0 */
{
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */
  int fstpos, ipos, ires;

  /* Open output PDB file containing the final flattened coordinates
     of the residues on the plot */
/* v.5.0--> */
  /* if ((ligplot_pdb = fopen("ligplot.pdb","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.pdb");
  if ((ligplot_pdb = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open output PDB file, ligplot.pdb\n"); */
      printf("\n*** Unable to open output PDB file, %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }

  /* Open output PDB file for SURFNET link */
/* v.5.0--> */
  /* if ((ligplot_frm = fopen("ligplot.frm","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.frm");
  if ((ligplot_frm = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open output PDB file, ligplot.frm\n"); */
      printf("\n*** Unable to open output PDB file, %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }

/* v.4.0--> */
  /* Open residue centres of mass file, ligplot.rcm */
/* v.5.0--> */
  /* if ((ligplot_rcm = fopen("ligplot.rcm","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.rcm");
  if ((ligplot_rcm = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
      printf("\n*** Unable to open output residue centres of mass file,");
/* v.5.0--> */
      /* printf(" ligplot.rcm\n"); */
      printf(" %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }
/* <--v.4.0 */

  /* If writing out the ligplot.res file, then open it */
  if (Write_Res_File == TRUE)
    {
/* v.5.0--> */
      /* if ((ligplot_res = fopen("ligplot.res","w")) ==  NULL) */
      filename[0] = '\0';
      if (wkdir[0] != '\0')
        strcpy(filename,wkdir);
      strcat(filename,"ligplot.res");
      if ((ligplot_res = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
        {
/* v.5.0--> */
          /* printf("\n*** Unable to open output file, ligplot.res\n"); */
          printf("\n*** Unable to open output file, %s\n",filename);
/* <--v.5.0 */
          exit(1);
        }

      /* Write the header records to the ligplot.res file */
      else
        {
          /* Print the residue range */
/* v.4.0--> */
/*        fprintf(ligplot_res,"#Residue range: %s to %s",res_num1,res_num2); */
          fprintf(ligplot_res,"#Residue range: ");
          if (res_name1[0] != '\0')
            fprintf(ligplot_res,"%s ",res_name1);
          fprintf(ligplot_res,"%s to ",res_num1);
          if (res_name2[0] != '\0')
            fprintf(ligplot_res,"%s ",res_name2);
          fprintf(ligplot_res,"%s",res_num2);
/* <--v.4.0 */
          if (chain != ' ')
            fprintf(ligplot_res,", chain %c",chain);
          fprintf(ligplot_res,"\n");

          /* Print the ligand residues identified */
          fprintf(ligplot_res,"#Ligand residues: ");
          ipos = 20;
          for (ires = 0; ires < ligand_residues; ires++)
            {
              /* If not been blanked out, then write out */
              if (lig_chain_store[ires] != '\0')
                {
                  fstpos = 0;
                  if (lig_res_num_store[ires][0] == ' ')
                    fstpos = 1;
                  if (lig_res_num_store[ires][1] == ' ')
                    fstpos = 2;
                  if (lig_res_num_store[ires][2] == ' ')
                    fstpos = 3;
                  fprintf(ligplot_res,"%s %s",lig_res_name_store[ires],
                          lig_res_num_store[ires]+fstpos);
                  ipos = ipos + 9 - fstpos;
                  if (lig_chain_store[ires] != ' ')
                    {
                      fprintf(ligplot_res,"(%c)",lig_chain_store[ires]);
                      ipos = ipos + 3;
                    }
                  if (ires < ligand_residues - 1)
                    {
                      fprintf(ligplot_res," - ");
                      ipos = ipos + 3;
                    }
                  if (ipos > 60)
                    {
/* v.4.0--> */
/*                    fprintf(ligplot_res,"\n                    "); */
                      fprintf(ligplot_res,"\n");
                      fprintf(ligplot_res,"#                 ");
/* <--v.4.0 */
                      ipos = 20;
                    }
                }
            }
          fprintf(ligplot_res,"\n");
        }
    }
}
/***********************************************************************
  
get_bonds  -  Read in the bond connectivities from the .bonds file

***********************************************************************/

void get_bonds(char bonds_name[FILENAME_LEN])
{
  char line[LINELEN + 1];
  char atom_name[2][5], chain[2], res_name[2][5], res_num[2][6];
  char bond_type_char, ebond_char, rbond_char;
  char bond_no[5], bond_len_str[8];

  int bond_order, bond_source, bond_type;
  int ipos, lpos;
  int end_bond, elastic, rotatable_bond, wanted;

  float bond_length;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *bond_atom_ptr[2];
  struct residue *residue_ptr;
    
  /* Initialise variables */
  printf("\nReading in bonds from file  %s ...\n",bonds_name);

  /* Open the ligplot.bonds file */
  if ((fil_bonds = fopen(bonds_name,"r")) == NULL)
    {
      printf("\n*** Unable to open bonds file, %s\n",bonds_name);
      exit(1);
    }

  /* Read through the .bonds file */
  while (fgets(line,100,fil_bonds) != NULL)
    {
      /* Initialise flag */
      wanted = FALSE;

      /* Check this is a bond record */
      if (!strncmp(line,"Bond",4))
        {
          /* Get the bond number */
          strncpy(bond_no,line+4,4);
          bond_no[4] = '\0';

          /* Assume this bond is going to be OK, unless it turns out
             otherwise */
          wanted = TRUE;

          /* Extract the bond details from the record */
          lpos = 22;

          /* Store both atoms */
          for (ipos = 0; ipos < 2; ipos++)
            {
              /* Atom name */
              strncpy(atom_name[ipos],line+lpos,4);
              atom_name[ipos][4] = '\0';
              lpos = lpos + 5;

              /* Residue name */
              strncpy(res_name[ipos],line+lpos,3);
              res_name[ipos][3] = '\0';
              lpos = lpos + 4;

              /* Residue number */
              strncpy(res_num[ipos],line+lpos,5);
              res_num[ipos][5] = '\0';
              lpos = lpos + 6;

              /* Chain */
              chain[ipos] = line[lpos];
              lpos = lpos + 12;
            }

          /* Get bond type */
          bond_type_char = line[66];
          if (bond_type_char == 'H')
            bond_type = HBOND;
          else if (bond_type_char == 'c')
            bond_type = COVALENT;
          else if (bond_type_char == 'n')
            bond_type = CONTACT;
          else if (bond_type_char == 'i')
            bond_type = INTERNAL;
/* v.5.4--> */
          else if (bond_type_char == 'S')
            bond_type = SALT_BRIDGE;
          else if (bond_type_char == 'D')
            bond_type = DISULPHIDE;
/* <--v.5.4 */
          else
            {
              printf("*** Warning. Invalid bond-type in ligplot.bonds");
              printf(" file  [%c]\n",bond_type_char);
              wanted = FALSE;
              Nwarnings++;
            }

          /* Determine whether a rotatable bond */
          end_bond = FALSE;
          rotatable_bond = FALSE;
          rbond_char = line[67];
          if (rbond_char == 'r')
            rotatable_bond = TRUE;
          else if (rbond_char == 'e')
            end_bond = TRUE;

          /* Determine whether an elastic bond */
          ebond_char = line[68];
          elastic = FALSE;
          if (ebond_char == 'E')
            elastic = TRUE;

          /* Get the bond length */
          strncpy(bond_len_str,line+69,7);
          bond_len_str[7] = '\0';
          bond_length = (float) atof(bond_len_str);

          /* Get the bond order */
          bond_order = SINGLE;
          if (line[77] == 'd')
            bond_order = DOUBLE;
          else if (line[77] == 't')
            bond_order = TRIPLE;

          /* Get the bond source */
          bond_source = CALCULATED;
          if (!strncmp(line+79,"HBPLUS",6))
            bond_source = HBPLUS;
          else if (!strncmp(line+79,"CONECT",6))
            bond_source = CONECT;

          /* Initialise pointers to the two bond atoms */
          bond_atom_ptr[0] = NULL;
          bond_atom_ptr[1] = NULL;

          /* Initialise pointer to start of linked list of atom coords */
          atom_ptr = first_atom_ptr;
              
          /* Loop through all the atoms to find the two making up the
             current bond */
          while (atom_ptr != NULL)
            {
              /* Get this atom's residue */
              residue_ptr = atom_ptr->residue_ptr;

              /* Check if this matches either of the bond atoms */
              for (ipos = 0; ipos < 2; ipos++)
                {
                  /* If have match, then store this atom's pointer */
                  if (!strncmp(res_name[ipos],residue_ptr->res_name,3) &&
                      !strncmp(res_num[ipos],residue_ptr->res_num,5) &&
                      chain[ipos] == residue_ptr->chain &&
                      !strncmp(atom_name[ipos],atom_ptr->atom_name,4))
                    bond_atom_ptr[ipos] = atom_ptr;
                }

              /* Get pointer to next atom in linked-list */
              atom_ptr = atom_ptr->next;
            }          

          /* Check whether matches for both the bond's atoms were found */
          if (bond_atom_ptr[0] == NULL || bond_atom_ptr[1] == NULL)
            {
              /* Print warning message */
              if (bond_atom_ptr[0] == NULL)
                {
                  printf("*** WARNING. Failed to find atom %s %s %s %c for ",
                         atom_name[0],res_name[0],res_num[0],chain[0]);
                  printf("bond number %s\n",bond_no);
                  Nwarnings++;
                }
              if (bond_atom_ptr[1] == NULL)
                {
                  printf("*** WARNING. Failed to find atom %s %s %s %c for ",
                         atom_name[1],res_name[1],res_num[1],chain[1]);
                  printf("bond number %s\n",bond_no);
                  Nwarnings++;
                }
              wanted = FALSE;
            }
        }

      /* If this bond is wanted, then create a new bond record and store
         the details */
      if (wanted == TRUE)
        {
          /* Create a new bond record */
/* v.4.3--> */
          bond_ptr
            = create_bond_record(first_bond_ptr,bond_type,bond_source,
                                 bond_length,bond_atom_ptr[0],
                                 bond_atom_ptr[1]);
/* <--v.4.3 */

          /* Store the bond details */
          bond_ptr->elastic = elastic;
          bond_ptr->rotatable_bond = rotatable_bond;
          bond_ptr->ring_bond = FALSE;
          bond_ptr->end_bond = end_bond;
/* v.4.4.4--> */
          bond_ptr->bond_order = bond_order;
/* <--v.4.4.4 */

          /* Store this pointer as the last-encountered bond */
          first_bond_ptr = bond_ptr;
        }
    }

/* v.4.2--> */
  /* Close the bonds file */
  fclose(fil_bonds);
/* <--v.4.2 */

  /*    printf("   Number of bonds in structure      = %7d\n",nbonds);
    if (bond_stack == 0)
      {
        printf("*** No bonds found in structure. Nothing to plot\n");
        exit(1);
      }
    printf("   No. of hydrogen bonds             = %7d\n",nhbonds);
    printf("   No. of linked bonds on stack      = %7d\n",bond_stack - 1);
    printf("   No. of ligand bonds on stack      = %7d\n",nligbonds);
    return bond_stack;  */
    
}
/***********************************************************************
  
covalent_connectivity  -  Find the covalent connectivity using a simple
                          distance cut-off

***********************************************************************/

void covalent_connectivity(void)
{
  int object_type1, object_type2, wanted1, wanted2;

  float x1, x2, y1, y2, z1, z2;
  float distance_real, distance_x, distance_y, distance_z, dist_sqrd;
  float length;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  printf("\nCalculating bonds and bond-connectivity ...\n");

  /* Initialise pointer to start of linked list of atom coords */
  atom1_ptr = first_atom_ptr;
  length = 0.0;
    
  /* Loop through all the atoms in the linked list */
  while (atom1_ptr != NULL)
    {
      /* Get the atom's coordinates */
      x1 = atom1_ptr->x;
      y1 = atom1_ptr->y;
      z1 = atom1_ptr->z;

      /* Get the residue and object that the atom belongs to */
      residue1_ptr = atom1_ptr->residue_ptr;
      object_type1 = residue1_ptr->object_ptr->object_type;

      /* Determine whether it is worth proceeding with this atom */
      wanted1 = TRUE;
      if (object_type1 == WATER || object_type1 == HYDROPHOBIC)
        wanted1 = FALSE;

      /* Initialise pointer to start of linked list of atom coords */
      atom2_ptr = first_atom_ptr;
    
      /* Loop through all the atoms in the linked list */
      while (atom2_ptr != NULL && wanted1 == TRUE)
        {
          /* Get the residue and object that the atom belongs to */
          residue2_ptr = atom2_ptr->residue_ptr;
          object_type2 = residue2_ptr->object_ptr->object_type;

          /* Determine whether it is worth proceeding with this atom */
          wanted2 = TRUE;
          if (object_type2 == WATER || object_type2 == HYDROPHOBIC)
            wanted2 = FALSE;

          /* If atoms both belong to H-bond groups, but are not in
             the same residue, then don't want any bonds between them */
          if (object_type1 == HGROUP && object_type2 == HGROUP &&
              residue1_ptr != residue2_ptr)
            wanted2 = FALSE;

          /* If one atom in the ligand, and the other isn't, then only
             want this bond if external groups parameter is set */
          if ((object_type1 == HGROUP && object_type2 == LIGAND) ||
              (object_type1 == LIGAND && object_type2 == HGROUP))
            {
              if (Include->External_Bonds == FALSE)
                wanted2 = FALSE;
            }

          /* If both atoms are the same atom, then don't want distance
             calculation */
          if (atom1_ptr == atom2_ptr)
            wanted2 = FALSE;

/* v.4.0--> */
          /* If this is an interface plot and have different residues
             belonging to the same surface, then not interested in any
             covalent bonds between them */
          if (Interface_Plot == TRUE && residue1_ptr != residue2_ptr)
            wanted2 = FALSE;
/* <--v.4.0 */

          /* Check for a covalent bond only if both atoms are wanted */
          if (wanted2 == TRUE)
            {
              /* Get the second atom's coordinates */
              x2 = atom2_ptr->x;
              y2 = atom2_ptr->y;
              z2 = atom2_ptr->z;

              /* Calculate distance between atoms 1 & 2 */
              distance_x = (x1 - x2) * (x1 - x2);
              distance_y = (y1 - y2) * (y1 - y2);
              distance_z = (z1 - z2) * (z1 - z2);
              distance_real = distance_x + distance_y + distance_z;

              /* Calculate cut-off distance for covalent connectivity,
                 based on atom types involved */
              dist_sqrd = (float) (3.6);
              if (atom1_ptr->atom_name[1] == 'P' ||
                  atom2_ptr->atom_name[1] == 'P')
                dist_sqrd = (float) (3.8);
              if (!strncmp(atom1_ptr->atom_name,"FE",2) ||
                  !strncmp(atom2_ptr->atom_name,"FE",2))
                dist_sqrd = (float) (5.0);
              if (atom1_ptr->atom_name[1] == 'S' ||
                  atom2_ptr->atom_name[1] == 'S')
                dist_sqrd = (float) (5.0);
              
              /* If distance satisfies the distance criteria for covalent
                 bonding, then atoms are bonded */
              if (distance_real < dist_sqrd)
                {
                  length = (float) sqrt((double) distance_real);

                  /* If this is a disulphide bond, store it as such */
/* v.5.4--> */
//                  if (!strncmp(atom1_ptr->atom_name," S  ",4) &&
//                      !strncmp(atom2_ptr->atom_name," S  ",4))
//                    update_bond(atom1_ptr,atom2_ptr,HBOND,CALCULATED,
                  if (!strncmp(atom1_ptr->atom_name," SG ",4) &&
                      !strncmp(atom2_ptr->atom_name," SG ",4) &&
                      distance_real <= DISULPH_DIST_SQRD)
/* <--v.5.4 */
                    update_bond(atom1_ptr,atom2_ptr,DISULPHIDE,CALCULATED,
                                length);

                  /* Otherwise, store as a covalent bond */
                  else
                    update_bond(atom1_ptr,atom2_ptr,COVALENT,CALCULATED,
                                length);
                }
            }
          atom2_ptr = atom2_ptr->next;
        }
      atom1_ptr = atom1_ptr->next;
    }
}
/* v.4.4--> */
/***********************************************************************

count_joins  -  Count the number of covalent links between the two
                residues, checking for any very close ones

***********************************************************************/

int count_joins(struct residue *residue1_ptr,
                struct residue *residue2_ptr,int *nclose)
{
  int njoin;
  int wanted;

  struct bond *bond_ptr;
  struct residue *res1_ptr, *res2_ptr;

  /* Initialise variables */
  njoin = 0;
  *nclose = 0;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Check that covalent bond */
      if (bond_ptr->bond_type == COVALENT)
        {
          /* Get this bond's two residue details */
          res1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
          res2_ptr = bond_ptr->second_atom_ptr->residue_ptr;

          /* Check if this pair of residues matches our residues of
             interest */
          wanted = FALSE;
          if ((res1_ptr == residue1_ptr && res2_ptr == residue2_ptr) ||
              (res2_ptr == residue1_ptr && res1_ptr == residue2_ptr))
            wanted = TRUE;

          /* If bond spans the two residues in question, then have a
            link between them */
          if (wanted == TRUE)
            {
              /* Increment count of joins */
              njoin++;

              /* Check length of bond to ensure it is not too short */
              if (bond_ptr->bond_length < TOO_CLOSE)
                (*nclose)++;
            }
        }

      /* Get the pointer to the next link */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return number of covalent connections between the two residues */
  return(njoin);
}
/***********************************************************************

delete_object  -  Delete the given object, marking all its atoms for
                  deletion

***********************************************************************/

void delete_object(struct object *object_ptr,
                   struct object *previous_object_ptr,int show_deletions)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Show message that object to be deleted */
  if (show_deletions == TRUE)
    {
      printf("+++ Residue %s %s %c not connected to any other",
             object_ptr->first_residue_ptr->res_name,
             object_ptr->first_residue_ptr->res_num,
             object_ptr->first_residue_ptr->chain);
      printf("  ... deleting it!\n");
      Nwarnings++;
    }

  /* Mark all the object's atoms for deletion */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop through all the object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
                  
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, deleting any unwanted
         ones */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Mark the atom for deletion */
          atom_ptr->deleted = TRUE;

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Mark the residue itself as deleted */
      residue_ptr->deleted = TRUE;

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Delete the current object by making the object pointers
     bypass it */
  if (previous_object_ptr == NULL)
    first_object_ptr = object_ptr->next_object_ptr;
              
  /* Otherwise, make pointers bypass the current link */
  else
    previous_object_ptr->next_object_ptr = object_ptr->next_object_ptr;
}
/***********************************************************************

delete_unwanted_objects  -  Delete any unwanted objects

***********************************************************************/

void delete_unwanted_objects(int *ndeleted,int show_deletions)
{
  int iobject;

  struct object *object_ptr;
  struct object *next_object_ptr, *previous_object_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  previous_object_ptr = NULL;
  iobject = 0;

  /* Loop through all objects to check which are to be deleted */
  while (object_ptr != NULL)
    {
      /* Save pointer to the next object */
      next_object_ptr = object_ptr->next_object_ptr;

      /* If object marked for deletion, then get rid of it */
      if (object_ptr->object_type == DELETED)
        {
          /* Delete the current object */
          delete_object(object_ptr,previous_object_ptr,show_deletions);
          (*ndeleted)++;

          /* Set current pointer to the previous object */
          object_ptr = previous_object_ptr;
        }

      /* Save the current object pointer */
      previous_object_ptr = object_ptr;

      /* Get pointer to point to next object in the list and increment
         object count */
      object_ptr = next_object_ptr;
      iobject++;
    }
}
/***********************************************************************
  
interpenetration_check  -  Check whether any residues interpenetrate

***********************************************************************/

int interpenetration_check(int ndel)
{
  int iresid1, iresid2, nresid1, nresid2;
  int nclose, ndeleted, njoin;
  int deleted1, deleted2, have_deletions, show_deletions;

  struct object *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  have_deletions = FALSE;
  ndeleted = ndel;
  show_deletions = TRUE;

  /* Initialise pointer to the first stored object */
  printf("\nChecking for residue interpenetration ...\n");
  object1_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object1_ptr != NULL)
    {
     /* Get object's first residue and number of residues */
      residue1_ptr = object1_ptr->first_residue_ptr;
      nresid1 = object1_ptr->nresidues;
      iresid1 = 0;
      deleted1 = FALSE;

      /* Loop over all this object's residues */
      while (iresid1 < nresid1 && residue1_ptr != NULL && deleted1 == FALSE)
        {
          /* Loop through all other objects beyond this one */
          object2_ptr = object1_ptr->next_object_ptr;

          /* Loop through all other objects */
          while (object2_ptr != NULL && deleted1 == FALSE)
            {
              /* Get object's first residue and number of residues */
              residue2_ptr = object2_ptr->first_residue_ptr;
              nresid2 = object2_ptr->nresidues;
              iresid2 = 0;
              deleted2 = FALSE;

              /* Loop over all this object's residues */
              while (iresid2 < nresid2 && residue2_ptr != NULL &&
                     deleted2 == FALSE)
                {
                  /* Check through the list of bonds to count the
                     number of covalent links between the current
                     pair of residues */
                  njoin = count_joins(residue1_ptr,residue2_ptr,
                                      &nclose);

                  /* If there are just too many bonds between these
                     two residues, then they probably interpenetrate */
                  if (nclose > 1 && njoin > 2)
                    {
                      /* Print warning */
                      printf("*** Warning. Residues [%s %s %c] and [%s %s %c]",
                             residue1_ptr->res_name,
                             residue1_ptr->res_num,
                             residue1_ptr->chain,
                             residue2_ptr->res_name,
                             residue2_ptr->res_num,
                             residue2_ptr->chain);
                      printf(" seem to interpenetrate\n");
                      printf("***          Number of joins = %d, ",njoin);
                      printf("  No. of close contacts = %d\n",nclose);
                      Nwarnings++;

                      /* Determine which object to delete */
                      deleted1 = FALSE;
                      deleted2 = FALSE;

                      /* Whichever object is the ligand, delete the other */
                      if (object1_ptr->object_type == LIGAND)
                        deleted2 = TRUE;
                      else if (object2_ptr->object_type == LIGAND)
                        deleted1 = TRUE;

                      /* Otherwise, delete the more minor object */
                      else if (object1_ptr->object_type == WATER)
                        deleted1 = TRUE;
                      else if (object2_ptr->object_type == WATER)
                        deleted2 = TRUE;
                      else if (object1_ptr->object_type == HYDROPHOBIC)
                        deleted1 = TRUE;
                      else if (object2_ptr->object_type == HYDROPHOBIC)
                        deleted2 = TRUE;
                      else
                        deleted2 = TRUE;

                      /* Mark the object to be deleted */
                      if (deleted1 == TRUE)
                        object1_ptr->object_type = DELETED;
                      else
                        object2_ptr->object_type = DELETED;

                      /* Set flag to indicate that we need to delete
                         at least one object */
                      have_deletions = TRUE;
                    }

                  /* Get pointer to the next residue */
                  residue2_ptr = residue2_ptr->next_residue_ptr;
                  iresid2++;
                }

              /* Get pointer to next object */
              object2_ptr = object2_ptr->next_object_ptr;
            }

          /* Get pointer to the next residue */
          residue1_ptr = residue1_ptr->next_residue_ptr;
          iresid1++;
        }

      /* Get pointer to next object */
      object1_ptr = object1_ptr->next_object_ptr;
    }

  /* If any objects have been marked for deletion, get rid of them now */
  if (have_deletions == TRUE)
    delete_unwanted_objects(&ndeleted,show_deletions);

  /* Return the number of object deleted */
  return(ndeleted);
}
/* <--v.4.4 */
/***********************************************************************

get_last  -  Get the last residue and atom of the given object

***********************************************************************/

void get_last(struct object *object_ptr,struct residue **last_residue_ptr,
              struct coordinate **last_atom_ptr)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise pointers */
  *last_atom_ptr = NULL;
  *last_residue_ptr = NULL;

  /* Loop through the first object's residues to get to the last
     residue and last atom */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop through all this object's residues, one by one */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Store this residue pointer */
      *last_residue_ptr = residue_ptr;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, counting number of
         mainchain atoms */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Store this atom pointer */
          *last_atom_ptr = atom_ptr;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

get_pointer_to_residue  -  Get the residue that points to the given
                             residue

***********************************************************************/

void get_pointer_to_residue(struct residue *current_residue_ptr,
                            struct residue **to_residue_ptr)
{
  struct residue *residue_ptr;

  /* Initialise pointer */
  *to_residue_ptr = NULL;

  /* Get pointer to very first residue */
  residue_ptr = first_residue_ptr;
  
  /* Loop through all residues */
  while (residue_ptr != NULL && *to_residue_ptr == NULL)
    {
      /* If this residue points to the given residue then store
         its pointer */
      if (residue_ptr->next_residue_ptr == current_residue_ptr)
        *to_residue_ptr = residue_ptr;

      /* Get pointer to next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

get_pointer_to_atom  -  Get the atom that points to the given
                             atom

***********************************************************************/

void get_pointer_to_atom(struct coordinate *current_atom_ptr,
                         struct coordinate **to_atom_ptr)
{
  struct coordinate *atom_ptr;

  /* Initialise pointer */
  *to_atom_ptr = NULL;

  /* Get pointer to very first atom */
  atom_ptr = first_atom_ptr;
  
  /* Loop through all atoms */
  while (atom_ptr != NULL && *to_atom_ptr == NULL)
    {
      /* If this atom points to the given atom then store
         its pointer */
      if (atom_ptr->next == current_atom_ptr)
        *to_atom_ptr = atom_ptr;

      /* Get pointer to next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

join_residues  -  Combine the given objects' residues and atoms into a
                  single object

***********************************************************************/

void join_residues(struct object *object1_ptr,struct object *object2_ptr,
                   int inligand)
{
  int iresid, nresid;

  struct coordinate *swap_atom_ptr;
  struct coordinate *last_atom_object1_ptr, *last_atom_object2_ptr;
  struct coordinate *to_first_atom_object1_ptr;
  struct coordinate *to_first_atom_object2_ptr;
  struct object *object_ptr, *to_object2_ptr;
  struct residue *residue_ptr, *swap_residue_ptr;
  struct residue *last_residue_object1_ptr, *last_residue_object2_ptr;
  struct residue *to_first_residue_object1_ptr;
  struct residue *to_first_residue_object2_ptr;

  /* Print that about to connect the residues */
  printf("\n");
  printf("*** Combining residues %s %s %c and %s %s %c",
         object1_ptr->first_residue_ptr->res_name,
         object1_ptr->first_residue_ptr->res_num,
         object1_ptr->first_residue_ptr->chain,
         object2_ptr->first_residue_ptr->res_name,
         object2_ptr->first_residue_ptr->res_num,
         object2_ptr->first_residue_ptr->chain);
  printf(" into a single object\n");
  Nwarnings++;

  /* Get the first object's last residue and last atom */
  get_last(object1_ptr,&last_residue_object1_ptr,&last_atom_object1_ptr);

  /* Repeat for second object */
  get_last(object2_ptr,&last_residue_object2_ptr,&last_atom_object2_ptr);

  /* Get the object that points to the second object */
  to_object2_ptr = NULL;
  object_ptr = first_object_ptr;
  
  /* Loop through all objects */
  while (object_ptr != NULL && to_object2_ptr == NULL)
    {
      /* If this object points to the second object, store its pointer */
      if (object_ptr->next_object_ptr == object2_ptr)
        to_object2_ptr = object_ptr;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Get the residue that points to the first object's first residue */
  get_pointer_to_residue(object1_ptr->first_residue_ptr,
                         &to_first_residue_object1_ptr);

  /* Get the residue that points to the second object's first residue */
  get_pointer_to_residue(object2_ptr->first_residue_ptr,
                         &to_first_residue_object2_ptr);

  /* Get the atom that points to the first object's first atom */
  get_pointer_to_atom(object1_ptr->first_residue_ptr->first_atom_ptr,
                      &to_first_atom_object1_ptr);

  /* Get the atom that points to the second object's first atom */
  get_pointer_to_atom(object2_ptr->first_residue_ptr->first_atom_ptr,
                      &to_first_atom_object2_ptr);

  /* If adding a residue to a ligand, make sure the first object's
     residues are all identified as ligand residues */
  if (inligand == TRUE)
    {
      residue_ptr = object1_ptr->first_residue_ptr;
      nresid = object1_ptr->nresidues;
      iresid = 0;
 
      /* Loop through all this object's residues, one by one */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Set ligand marker */
          residue_ptr->inligand = TRUE;

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }
    }

  /* Loop through all the residues of the second object and make them
     point to the first object */
  residue_ptr = object2_ptr->first_residue_ptr;
  nresid = object2_ptr->nresidues;
  iresid = 0;
 
  /* Loop through all this object's residues, one by one */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Amend the object pointer */
      residue_ptr->object_ptr = object1_ptr;

      /* If in ligand, then mark as such */
      if (inligand == TRUE)
        residue_ptr->inligand = TRUE;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Cut out the second object altogether */
  if (to_object2_ptr != NULL)
    to_object2_ptr->next_object_ptr = object2_ptr->next_object_ptr;
  else if (first_object_ptr == object2_ptr)
    first_object_ptr = object2_ptr->next_object_ptr;

  /* If residues are already adjacent, then don't need to do anything */
  if (last_residue_object1_ptr->next_residue_ptr
      != object2_ptr->first_residue_ptr)
    {
      /* Move the second object's residues to their new location */
      swap_residue_ptr = last_residue_object1_ptr->next_residue_ptr;
      last_residue_object1_ptr->next_residue_ptr
        = object2_ptr->first_residue_ptr;
      if (to_first_residue_object2_ptr != NULL)
        to_first_residue_object2_ptr->next_residue_ptr
          = last_residue_object2_ptr->next_residue_ptr;
      else if (first_residue_ptr == object2_ptr->first_residue_ptr)
        first_residue_ptr = last_residue_object2_ptr->next_residue_ptr;
      last_residue_object2_ptr->next_residue_ptr = swap_residue_ptr;
    }

  /* If atoms are already adjacent, then don't need to do anything */
  if (last_atom_object1_ptr->next
      != object2_ptr->first_residue_ptr->first_atom_ptr)
    {
      /* Move the second object's atoms to their new location */
      swap_atom_ptr = last_atom_object1_ptr->next;
      last_atom_object1_ptr->next
        = object2_ptr->first_residue_ptr->first_atom_ptr;
      if (to_first_atom_object2_ptr != NULL)
        to_first_atom_object2_ptr->next
          = last_atom_object2_ptr->next;
      else if (first_atom_ptr
               == object2_ptr->first_residue_ptr->first_atom_ptr)
        first_atom_ptr = last_atom_object2_ptr->next;
      last_atom_object2_ptr->next = swap_atom_ptr;
    }

  /* Update the first object's residue-count */
  object1_ptr->nresidues = object1_ptr->nresidues + object2_ptr->nresidues;
  object1_ptr->nmain_chain = object1_ptr->nmain_chain
    + object2_ptr->nmain_chain;
  object1_ptr->n_ca = object1_ptr->n_ca + object2_ptr->n_ca;
}
/***********************************************************************

check_for_join  -  Check whether the given two residues are connected by
                   a covalent bond

***********************************************************************/

int check_for_join(struct residue *residue1_ptr,
                   struct residue *residue2_ptr)
{
  char chain1, chain2, res_name1[4], res_name2[4], res_num1[6], res_num2[6];

  int connected, match11, match12, match21, match22;

  struct hhb_info *hhb_info_ptr;

  /* Initialise flag */
  connected = FALSE;

  /* Get pointer to first stored bond item in the list */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Loop through all the bonds in the linked list */
  while (hhb_info_ptr != NULL)
    {
      /* If this is a covalent bond, check whether
         it links the two residues */
      if (hhb_info_ptr->source == CONECT)
        {
          /* Get this bond's two residue details */
          strncpy(res_name1,hhb_info_ptr->res_name1,3);
          res_name1[3] = '\0';
          strncpy(res_num1,hhb_info_ptr->res_num1,5);
          res_num1[5] = '\0';
          chain1 = hhb_info_ptr->chain1;
          strncpy(res_name2,hhb_info_ptr->res_name2,3);
          res_name2[3] = '\0';
          strncpy(res_num2,hhb_info_ptr->res_num2,5);
          res_num2[5] = '\0';
          chain2 = hhb_info_ptr->chain2;

          /* Check if the two residues are represented */

          /* First residue matches res1 of bond */
          match11 = FALSE;
          if (!strncmp(residue1_ptr->res_name,
                       res_name1,3) &&
              !strncmp(residue1_ptr->res_num,
                       res_num1,5) &&
              residue1_ptr->chain == chain1)
            match11 = TRUE;

          /* First residue matches res2 of bond */
          match12 = FALSE;
          if (!strncmp(residue1_ptr->res_name,
                       res_name2,3) &&
              !strncmp(residue1_ptr->res_num,
                       res_num2,5) &&
              residue1_ptr->chain == chain2)
            match12 = TRUE;

          /* Second residue matches res1 of bond */
          match21 = FALSE;
          if (!strncmp(residue2_ptr->res_name,
                       res_name1,3) &&
              !strncmp(residue2_ptr->res_num,
                       res_num1,5) &&
              residue2_ptr->chain == chain1)
            match21 = TRUE;

          /* Second residue matches res2 of bond */
          match22 = FALSE;
          if (!strncmp(residue2_ptr->res_name,
                       res_name2,3) &&
              !strncmp(residue2_ptr->res_num,
                       res_num2,5) &&
              residue2_ptr->chain == chain2)
            match22 = TRUE;

          /* If it spans the two residues in question,
             then have a link between them */
          if ((match11 && match22) ||
              (match12 && match21))
            {
              /* Set flag on */
              connected = TRUE;

              /* Upgrade this bond so that included
                 later on */
              hhb_info_ptr->source = CONECT;
/* v.5.4--> */
              hhb_info_ptr->bond_type = S_COVALENT;
/* <--v.5.4 */
            }
        }

      /* Get the pointer to the next link */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }

  /* Return flag indicating whether connection between the two residues
     found or not */
  return(connected);
}
/* v.4.0--> */
/***********************************************************************

check_standard_aa  -  Check whether given residue is a standard amino
                      acid

***********************************************************************/

int check_standard_aa(char res_name[4])
{
  int standard;
  int iamino;

  static char *amino[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL"
  };
    
  /* Initialise return-value */
  standard = FALSE;

  /* Check whether this is a standard amino acid */
  for (iamino = 0; iamino < 20 && standard == FALSE; iamino++)
    if (!strncmp(res_name,amino[iamino],3))
      standard = TRUE;

  /* Return whether given residue is a standard amino acid */
  return(standard);
}
/* <--v.4.0 */
/***********************************************************************

combine_objects  -  Check for any H-group objects that need to be
                    combined (ie one is an attachment for another)

***********************************************************************/

void combine_objects(void)
{
  int connected, inligand, joined;
  int iresid1, iresid2, nresid1, nresid2;
/* v.4.0--> */
  int standard;
/* <--v.4.0 */

  struct object *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to the first stored object */
  object1_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object1_ptr != NULL)
    {
/* v.4.1.2--> */
      /* Initialise flag */
      joined = FALSE;
/* <--v.4.1.2 */

      /* Only process if this is not a water */
      if (object1_ptr->object_type != WATER)
        {
          /* Get object's first residue and number of residues */
          residue1_ptr = object1_ptr->first_residue_ptr;
          nresid1 = object1_ptr->nresidues;
          iresid1 = 0;
          joined = FALSE;

          /* Loop over all this object's residues */
          while (iresid1 < nresid1 && residue1_ptr != NULL &&
                 joined == FALSE)
            {
              /* Loop through all other non-ligand objects except this one */
              object2_ptr = object1_ptr->next_object_ptr;

              /* Loop through all other objects */
              while (object2_ptr != NULL && joined == FALSE)
                {
                  /* Only consider if this is not a water */
                  if (object2_ptr->object_type != WATER)
                    {
                      /* Get object's first residue and number of residues */
                      residue2_ptr = object2_ptr->first_residue_ptr;
                      nresid2 = object2_ptr->nresidues;
                      iresid2 = 0;

                      /* Loop over all this object's residues */
                      while (iresid2 < nresid2 && residue2_ptr != NULL &&
                             joined == FALSE)
                        {
                          /* Check whether residue number and chain match,
                             or both are ligand residues */
/* v.4.2.2--> */
/*                        if (!strncmp(residue1_ptr->res_num,
                                       residue2_ptr->res_num,5) &&
                              residue1_ptr->chain == residue2_ptr->chain) */
                          if ((!strncmp(residue1_ptr->res_num,
                                       residue2_ptr->res_num,5) &&
                              residue1_ptr->chain == residue2_ptr->chain) ||
                              (object1_ptr->object_type == LIGAND &&
                               object2_ptr->object_type == LIGAND))
/* <--v.4.2.2 */
                            {
                              /* Check through the list of bonds to confirm
                                 that there is indeed a link between these
                                 objects */
                              connected = check_for_join(residue1_ptr,
                                                         residue2_ptr);

/* v.4.0--> */
                              /* If one object is a ligand, and the other
                                 isn't, then check that the residue from
                                 the non-ligand object is not a standard
                                 amino acid */
                              if (object1_ptr->object_type == LIGAND ||
                                  object2_ptr->object_type == LIGAND)
                                {
                                  standard = FALSE;

                                  /* For whichever object is the non-ligand
                                     one, check whether the residue is a
                                     standard amino acid */
                                  if (object1_ptr->object_type != LIGAND)
                                    standard = 
                                      check_standard_aa(residue1_ptr->res_name);
                                  else if (object2_ptr->object_type != LIGAND)
                                    standard = 
                                      check_standard_aa(residue2_ptr->res_name);

                                  /* If non-ligand residue is a standard
                                     amino acid, then don't want to join
                                     the two objects */
                                  if (standard == TRUE)
                                    connected = FALSE;
                                }
/* <--v.4.0 */

                              /* If have a link between the two residues,
                                 then may be an attachment */
                              if (connected == TRUE)
                                {
                                  /* If either object is a ligand, then
                                     make sure both are marked as such */
                                  if (object1_ptr->object_type == LIGAND ||
                                      object2_ptr->object_type == LIGAND)
/* v.4.0--> */
                                    {
                                      object1_ptr->object_type = LIGAND;
                                      object2_ptr->object_type = LIGAND;
/* <--v.4.0 */
                                      inligand = TRUE;
/* v.4.0--> */
                                    }
/* <--v.4.0 */
                                  else
                                    inligand = FALSE;

                                  /* Join the residues into a single object */
/* v.4.0--> */
/*                                join_residues(object1_ptr,object2_ptr,
                                                inligand); */
                                  join_residues(object1_ptr,object2_ptr,
                                                inligand);
/* <--v.4.0 */
                                  joined = TRUE;
                                }
                            }

                          /* Get pointer to the next residue */
                          residue2_ptr = residue2_ptr->next_residue_ptr;
                          iresid2++;
                        }
                    }

                  /* Get pointer to next object */
                  object2_ptr = object2_ptr->next_object_ptr;
                }

              /* Get pointer to the next residue */
              residue1_ptr = residue1_ptr->next_residue_ptr;
              iresid1++;
            }
        }

      /* Get pointer to next object */
      if (joined == FALSE)
        object1_ptr = object1_ptr->next_object_ptr;
    }
}
/***********************************************************************

delete_false_ligands  -  Delete false ligands (eg protein residues with
                         the same residue numbers as the ligand residues)

***********************************************************************/

void delete_false_ligands(void)
{
/* v.4.3--> */
  char element[3];
/* <--v.4.3 */

  int iresid, nligand, nligand_objects, nobject, nresid;
  int delete_last, delete_ligand, is_metal, metal, water;
  int natoms;
  int ires;

  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  nligand_objects = 0;
  nobject = 0;
  delete_last = FALSE;
  natoms = 0;
  nresid = 0;
  metal = FALSE;
  water = FALSE;

  /* Loop through all objects to count the number of ligand objects
     and make sure that still have more than one */
  while (object_ptr != NULL)
    {
      /* If this is a ligand, then increment count of ligand objects */
      if (object_ptr->object_type == LIGAND)
        {
          /* Get the object's first residue pointer */
          residue_ptr = object_ptr->first_residue_ptr;

          /* If object is a water, and are plotting waters, then this
             is OK as the ligand */
          water = FALSE;
          if (!strncmp(residue_ptr->res_name,"HOH",3) ||
              !strncmp(residue_ptr->res_name,"WAT",3))
            if (Water_as_Ligand == TRUE)
              water = TRUE;

          /* If object is a metal, and are plotting metals, then this
             is OK as the ligand */
          metal = FALSE;
          if (object_ptr->nresidues == 1 && residue_ptr->natoms == 1)
            {
              /* Check if this might be a metal */
/* v.4.3--> */
/*            is_metal
                = check_for_metal(residue_ptr->first_atom_ptr->atom_name); */
              is_metal
                = check_for_metal(residue_ptr->first_atom_ptr->atom_name,
                                  element);
/* <--v.4.3 */
              if (Metal_as_Ligand == TRUE && is_metal == TRUE)
                metal = TRUE;
            }

          /* Increment number of apparent ligands and get number of
             residues in this object */
          nligand_objects++;
          nresid = object_ptr->nresidues;

          /* Get the number of atoms in the object's first residue */
          natoms = residue_ptr->natoms;
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* If only have one ligand left, then don't need to delete any objects */
  if (nligand_objects < 2)
    return;

  /* If last apparent ligand consisted of a single residue containing a
     single atom (and this wasn't a water), then take the penultimate
     ligand as the "true" one, rather than the last one (ie probably just
     a metal ion with the same residue number as a HET group, and probably
     the same number as some residue in the protein!) */
  if (nresid == 1 && natoms == 1 && water == FALSE && metal == FALSE)
    delete_last = TRUE;

  /* Initialise variables */
  nligand = 0;
  printf("*** Warning. Number of ligand objects identified: %d\n",
         nligand_objects);
  printf("***          Retaining last one as ligand only\n");
  Nwarnings++;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If this is a ligand, then increment count of ligand objects */
      if (object_ptr->object_type == LIGAND)
        nligand++;

      /* If a ligand object, then delete or keep it */
      if (object_ptr->object_type == LIGAND)
        {
          printf("    Ligand object:  %s %s %c",
                 object_ptr->first_residue_ptr->res_name,
                 object_ptr->first_residue_ptr->res_num,
                 object_ptr->first_residue_ptr->chain);

          /* If this is not the last ligand object, then delete */
          nobject++;

          /* Determine whether to delete this particular object */
          delete_ligand = FALSE;
          if (delete_last == FALSE && nobject < nligand_objects)
            delete_ligand = TRUE;
          else if (delete_last == TRUE && nobject != nligand_objects - 1)
            delete_ligand = TRUE;

          if (delete_ligand == TRUE)
            {
              /* Reclassify this object as a H-bond group */
              object_ptr->object_type = HGROUP;
              printf("  *** Marked as non-ligand\n");

              /* Get pointer to object's first residue */
              residue_ptr = object_ptr->first_residue_ptr;
              nresid = object_ptr->nresidues;
              iresid = 0;
              
              /* Loop over all this object's residues to classify them
                 as non-ligand residues */
              while (iresid < nresid && residue_ptr != NULL)
                {
                  /* Define residue as non-ligand */
                  residue_ptr->inligand = FALSE;

              /* Loop through the list of ligand residues to remove this
                 entry from it */
              for (ires = 0; ires < ligand_residues; ires++)
                {
                  /* If this is the one just deleted, then remove it
                     from the list */
                  if (lig_chain_store[ires] == residue_ptr->chain &&
                      !strncmp(lig_res_name_store[ires],
                               residue_ptr->res_name,3) &&
                      !strncmp(lig_res_num_store[ires],
                               residue_ptr->res_num,3))
                    {
                      lig_chain_store[ires] = '\0';
                      lig_res_name_store[ires][0] = '\0';
                      lig_res_num_store[ires][0] = '\0';
                    }
                }
                  /* Get pointer to the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }
            }

          /* Keep the last ligand object */
          else
            printf("      Retained as ligand\n");
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

add_hbplus_bonds  -  Match the atoms listed as H-bonded in the .hhb file
                     (or as involved in non-bonded contacts in the .nnb
                     file) with those read in from the PDB file

***********************************************************************/

void add_hbplus_bonds(void)
{
  int bond_source, bond_type, match, wanted;
/* v.4.0--> */
  int in_lig1, in_lig2;

  char dummy[4];
/* <--v.4.0 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */

  float length, x1, x2, y1, y2, z1, z2;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct hhb_info *hhb_info_ptr;
  struct residue *residue_ptr;

/* v.4.2.2--> */
  /* Initialise */
  mol_id = -1;
/* <--v.4.2.2 */

  /* Set pointer to the first of the HBPLUS interactions */
  hhb_info_ptr = first_hhb_info_ptr;

  /* Loop through each of the atom-pairs corresponding to the
     hydrogen-bonds read in from the .hhb file and the non-bonded
     contacts read in from the .nnb file */
  while (hhb_info_ptr != NULL)
    {
      /* Determine bond type */
      wanted = TRUE;
      if (hhb_info_ptr->source == CONECT)
        bond_type = COVALENT;
/* v.5.4--> */
//      else if (hhb_info_ptr->source == HHB_FILE)
//        bond_type = HBOND;
//      else if (hhb_info_ptr->source == NNB_FILE)
//        bond_type = CONTACT;
/* <--v.5.4.6 */
      else if (hhb_info_ptr->source == HHB_FILE)
        bond_type = HBOND;
/* v.5.4.6--> */
      else if (hhb_info_ptr->bond_type == S_HBONDS)
        bond_type = HBOND;
      else if (hhb_info_ptr->bond_type == S_CONTACTS)
        bond_type = CONTACT;
      else if (hhb_info_ptr->bond_type == S_SALT_BRIDGES)
        bond_type = SALT_BRIDGE;
      else if (hhb_info_ptr->bond_type == S_DISULPHIDES)
        bond_type = DISULPHIDE;
/* <--v.5.4 */
      else
        wanted = FALSE;

      /* Initialise pointer to start of linked list of atom coords */
      atom_ptr = first_atom_ptr;
      atom1_ptr = NULL;
      atom2_ptr = NULL;
      match = FALSE;

      /* Loop through all the atoms in the linked list */
      while (atom_ptr != NULL && match == FALSE && wanted == TRUE)
        {
          /* Get this atom's residue */
          residue_ptr = atom_ptr->residue_ptr;

          /* Check whether atom corresponds to the first of the H-bond
             atoms */
          if (residue_ptr->chain == hhb_info_ptr->chain1 &&
              !strncmp(residue_ptr->res_name,hhb_info_ptr->res_name1,3) &&
              !strncmp(residue_ptr->res_num,hhb_info_ptr->res_num1,5) &&
              !strncmp(atom_ptr->atom_name,hhb_info_ptr->atom_name1,4))
            atom1_ptr = atom_ptr;

          /* Check whether atom corresponds to the second of the H-bond
             atoms */
          else if (residue_ptr->chain == hhb_info_ptr->chain2 &&
              !strncmp(residue_ptr->res_name,hhb_info_ptr->res_name2,3) &&
              !strncmp(residue_ptr->res_num,hhb_info_ptr->res_num2,5) &&
              !strncmp(atom_ptr->atom_name,hhb_info_ptr->atom_name2,4))
            atom2_ptr = atom_ptr;

/* v.4.0--> */
          /* If this is a CONECT record, then check that is between
             ligand and non-ligand residue */
          if (atom1_ptr != NULL && atom2_ptr != NULL &&
              hhb_info_ptr->source == CONECT)
            {
              /* Check if first residue is in the ligand */
              in_lig1 = check_if_in_ligand(hhb_info_ptr->chain1,
                                           hhb_info_ptr->res_name1,
                                           hhb_info_ptr->res_num1,
/* v.4.2.2--> */
/*                                         FALSE,dummy); */
                                           FALSE,dummy,mol_id);
/* <--v.4.2.2 */

              /* Check if first residue is in the ligand */
              in_lig2 = check_if_in_ligand(hhb_info_ptr->chain2,
                                           hhb_info_ptr->res_name2,
                                           hhb_info_ptr->res_num2,
/* v.4.2.2--> */
/*                                         FALSE,dummy); */
                                           FALSE,dummy,mol_id);
/* <--v.4.2.2 */

              /* If neither residue is in the ligand, discount
                 this bond */
              if (in_lig1 == FALSE && in_lig2 == FALSE &&
                  atom1_ptr->residue_ptr->object_ptr !=
                  atom2_ptr->residue_ptr->object_ptr)
                {
                  atom1_ptr = NULL;
                  atom2_ptr = NULL;
                }
/* v.4.1--> */
              /* If the external bonds option is not set then don't
                 want this if between ligand and nonligand */
              if (Include->External_Bonds == FALSE)
                {
                  /* Check whether between ligand and non-ligand */
                  if ((in_lig1 == TRUE && in_lig2 == FALSE) ||
                      (in_lig1 == FALSE && in_lig2 == TRUE))
                    {
                      atom1_ptr = NULL;
                      atom2_ptr = NULL;
                    }
                }
/* <--v.4.1 */
            }
/* <--v.4.0 */

          /* If have found the two atoms corresponding to both ends
             of the bond, then add bond to list */
          if (atom1_ptr != NULL && atom2_ptr != NULL)
            {
              /* Set flag to denote that atoms matched */
              match = TRUE;

              /* Calculate H-bond distance */
              x1 = atom1_ptr->x;
              y1 = atom1_ptr->y;
              z1 = atom1_ptr->z;
              x2 = atom2_ptr->x;
              y2 = atom2_ptr->y;
              z2 = atom2_ptr->z;
              length = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
                + (z1 - z2) * (z1 - z2);
              length = (float) sqrt((double) length);

              /* Determine bond's source */
              if (hhb_info_ptr->source == CONECT)
                bond_source = CONECT;
              else
                bond_source = HBPLUS;

              /* Add bond to list */
              update_bond(atom1_ptr,atom2_ptr,bond_type,bond_source,
                          length);
            }

           /* Get pointer to next atom in linked-list */
           atom_ptr = atom_ptr->next;
         }

      /* Get the pointer to the next HBPLUS interaction */
      hhb_info_ptr = hhb_info_ptr->next_hhb_info_ptr;
    }
}
/***********************************************************************

open_bonds_output  -  Open output bonds files, ligplot.hhb, ligplot.nnb
                      and ligplot.bonds

***********************************************************************/

void open_bonds_output(void)
{
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */

  /* Open the ligplot.bonds file for output */
/* v.5.0--> */
  /* if ((ligplot_bonds_out = fopen("ligplot.bonds","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.bonds");
  if ((ligplot_bonds_out = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open output bonds file, ligplot.bonds\n"); */
      printf("\n*** Unable to open output bonds file, %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }

  /* Open the ligplot.hhb file for output */
/* v.5.0--> */
  /* if ((ligplot_hhb_out = fopen("ligplot.hhb","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.hhb");
  if ((ligplot_hhb_out = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open output H-bonds file, ligplot.hhb\n"); */
      printf("\n*** Unable to open output H-bonds file, %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }

  /* Write header records to ligplot.hhb */
  fprintf(ligplot_hhb_out,"ligplot.hhb output:\n");    
  fprintf(ligplot_hhb_out,"\n");    
  fprintf(ligplot_hhb_out,"   Donor                  Acceptor     ");
  fprintf(ligplot_hhb_out,"Distance\n");

  /* Open the ligplot.nnb file for output */
/* v.5.0--> */
  /* if ((ligplot_nnb_out = fopen("ligplot.nnb","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.nnb");
  if ((ligplot_nnb_out = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open output contacts file, ligplot.nnb\n"); */
      printf("\n*** Unable to open output contacts file, %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }

  /* Write header records to ligplot.nnb */
  fprintf(ligplot_nnb_out,"ligplot.nnb output:\n");    
  fprintf(ligplot_nnb_out,"\n");    
  fprintf(ligplot_nnb_out,"   Atom 1               Atom 2        ");
  fprintf(ligplot_nnb_out,"Distance\n");

/* v.4.1.1--> */
  /* Open the ligplot.sum file for output */
/* v.5.0--> */
  /* if ((ligplot_sum = fopen("ligplot.sum","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.sum");
  if ((ligplot_sum = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open output summary file, ligplot.sum\n"); */
      printf("\n*** Unable to open output summary file, %s\n",filename);
/* <--v.5.0 */
      exit(1);
    }

  /* Write header records to ligplot.sum */
  fprintf(ligplot_sum,"%s\n",Print_Title);
  fprintf(ligplot_sum,"\n");
  fprintf(ligplot_sum,"                                             ");
  fprintf(ligplot_sum,"Hydrogen bonds  Non-bonded\n");
  fprintf(ligplot_sum,"               Interacting residues          ");
  fprintf(ligplot_sum,"M-M   S-S   M-S  contacts\n");
  fprintf(ligplot_sum,"          ------------------------------     ");
  fprintf(ligplot_sum,"---   ---   ---    -----  \n");
/* <--v.4.1.1 */
}
/* v.4.3--> */
/***********************************************************************
      
get_bond_orders  -  Read in the hbadd.bonds file, if present, to pick up
                    the bond order information written out by hbadd.c
                    from the components.cif file

***********************************************************************/

void get_bond_orders(void)
{
  char atom_name1[5], atom_name2[5], bond_desc[5];
  char res_name1[4], res_num1[6];
  char res_name2[4], res_num2[6];
  char input_line[LINELEN + 1];
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */

  int bond_order;
  int found, have_file;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  have_file = TRUE;

  /* Open hbadd.bonds file to see if it is present */
/* v.5.0--> */
  /* if ((file_ptr = fopen("hbadd.bonds","r")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"hbadd.bonds");
/* v.5.4.5 --> */
  if (PDBsum1 == TRUE)
    strcpy(filename,"../hbplus/hbadd.bonds");
/* <-- v.5.4.5 */
  if ((file_ptr = fopen(filename,"r")) == NULL)
/* <--v.5.0 */
    have_file = FALSE;

  /* If file is present, read through all the records */
  if (have_file == TRUE)
    {
      /* Read in the data from the PDB file */
      printf("\nReading in the bond order data from hbadd.bonds ...\n");

      /* Loop through the file, picking up the bond order data */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* Only interested in this line if the bond order is anything
             other than a single bond */
          if (strncmp(input_line+52,"SING",4))
            {
              /* Extract the data for atom 1 */
              strncpy(atom_name1,input_line+6,4);
              atom_name1[4] = '\0';
              strncpy(res_name1,input_line+11,3);
              res_name1[3] = '\0';
              strncpy(res_num1,input_line+15,5);
              res_num1[5] = '\0';
          
              /* Repeat for atom 2 */
              strncpy(atom_name2,input_line+33,4);
              atom_name2[4] = '\0';
              strncpy(res_name2,input_line+38,3);
              res_name2[3] = '\0';
              strncpy(res_num2,input_line+42,5);
              res_num2[5] = '\0';

              /* Extract the bond order */
              strncpy(bond_desc,input_line+52,4);
              bond_desc[4] = '\0';
              if (!strncmp(bond_desc,"SING",4))
                bond_order = SINGLE;
              else if (!strncmp(bond_desc,"DOUB",4))
                bond_order = DOUBLE;
              else if (!strncmp(bond_desc,"TRIP",4))
                bond_order = TRIPLE;
              else if (!strncmp(bond_desc,"AROM",4))
                bond_order = AROMATIC;

              /* Initialise pointer to start of linked list of bonds */
              bond_ptr = first_bond_ptr;
              found = FALSE;
    
              /* Loop through all the bonds in the linked list,
                 excluding any elastic bonds */
              while (bond_ptr != NULL && found == FALSE)
                {
                  /* Get the two atoms at either end of this bond */
                  atom1_ptr = bond_ptr->first_atom_ptr;
                  atom2_ptr = bond_ptr->second_atom_ptr;

                  /* Get the residues to which these two atoms belong */
                  residue1_ptr = atom1_ptr->residue_ptr;
                  residue2_ptr = atom2_ptr->residue_ptr;

                  /* Check for a match */
                  if ((!strcmp(residue1_ptr->res_name,res_name1) &&
                       !strcmp(atom1_ptr->atom_name,atom_name1) &&
                       !strcmp(residue2_ptr->res_name,res_name2) &&
                       !strcmp(atom2_ptr->atom_name,atom_name2)) ||
                      (!strcmp(residue2_ptr->res_name,res_name2) &&
                       !strcmp(atom2_ptr->atom_name,atom_name2) &&
                       !strcmp(residue1_ptr->res_name,res_name1) &&
                       !strcmp(atom1_ptr->atom_name,atom_name1)))
                    found = TRUE;

                  /* If have a match, then store this bond's bond order */
                  if (found == TRUE)
                    {
                      bond_ptr->bond_order = bond_order;
                    }

                  /* Get pointer to next bond in linked-list */
                  bond_ptr = bond_ptr->next_bond_ptr;
                }

            }
        }

      /* Close the file */
      fclose(file_ptr);
    }
}
/***********************************************************************
      
get_order  -  Get bond order for bond in standard residue type

***********************************************************************/

int get_order(char res_name[4],char atom_name1[5],char atom_name2[5])
{
#define NENTRY 41

  int bond_order, ientry;
  int found;

  static char dbond[NENTRY][3][5] = 
    {
      {"***", " C  ", " O  "},
      {"ARG", " CZ ", " NH1"},
      {"ASN", " CG ", " OD1"},
      {"ASP", " CG ", " OD1"},
      {"GLN", " CD ", " OE1"},
      {"GLU", " CD ", " OE1"},
      {"HIS", " CG ", " CD2"},
      {"HIS", " ND1", " CE1"},
      {"PHE", " CG ", " CD1"},
      {"PHE", " CE1", " CZ "},
      {"PHE", " CD2", " CE2"},
      {"TRP", " CG ", " CD1"},
      {"TRP", " CD2", " CE2"},
      {"TRP", " CZ2", " CH2"},
      {"TRP", " CE3", " CZ3"},
      {"TYR", " CG ", " CD1"},
      {"TYR", " CE1", " CZ "},
      {"TYR", " CD2", " CE1"},
      {"  A", " C6 ", " N1 "},
      {"  A", " C2 ", " N3 "},
      {"  A", " C4 ", " C5 "},
      {"  A", " C8 ", " N7 "}, {"  A", " P  ", " O1P"},
      {"  C", " C2 ", " O2 "},
      {"  C", " C4 ", " N3 "},
      {"  C", " C5 ", " C6 "}, {"  C", " P  ", " O1P"},
      {"  G", " C6 ", " O6 "},
      {"  G", " C2 ", " N3 "},
      {"  G", " C4 ", " C5 "},
      {"  G", " C8 ", " N7 "}, {"  G", " P  ", " O1P"},
      {"  T", " C2 ", " O2 "},
      {"  T", " C4 ", " O4 "},
      {"  T", " C5 ", " C6 "}, {"  T", " P  ", " O1P"},
      {"  U", " C2 ", " O2 "},
      {"  U", " C4 ", " O4 "},
      {"  U", " C5 ", " O6 "}, {"  U", " P  ", " O1P"}
    };

  /* Initialise variables */
  bond_order = SINGLE;
  found = FALSE;

  /* Loop through the bond entries until find the right one */
  for (ientry = 0; ientry < NENTRY && found == FALSE; ientry++)
    {
      /* Check whether have a match */
      if (!strcmp(res_name,dbond[ientry][0]) ||
          !strcmp(dbond[ientry][0],"***"))
        {
          /* Check atom-name match */
          if ((!strcmp(atom_name1,dbond[ientry][1]) &&
               !strcmp(atom_name2,dbond[ientry][2])) ||
              (!strcmp(atom_name2,dbond[ientry][1]) &&
               !strcmp(atom_name1,dbond[ientry][2])))
            {
              /* Set bond order to double */
              bond_order = DOUBLE;

              /* Set flag that found */
              found = TRUE;
            }
        }
    }

  /* Return the bond order */
  return(bond_order);
}
/***********************************************************************
      
get_standard_bond_orders  -  Apply bond orders to standard residue types

***********************************************************************/

void get_standard_bond_orders(void)
{
  char atom_name1[5], atom_name2[5];
  char res_name[4];

  int bond_order;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
    
  /* Loop through all the bonds in the linked list, excluding any
     elastic bonds and any already done */
  while (bond_ptr != NULL)
    {
      /* Check that not an elastic bond and not already a double */
      if (bond_ptr->bond_type == COVALENT && bond_ptr->bond_order == SINGLE)
        {
          /* Get the two atoms at either end of this bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the residues to which these two atoms belong */
          residue1_ptr = atom1_ptr->residue_ptr;
          residue2_ptr = atom2_ptr->residue_ptr;

          /* If both atoms belong to the same residue, then process */
          if (residue1_ptr == residue2_ptr)
            {
              /* Get the residue name */
              strcpy(res_name,residue1_ptr->res_name);

              /* Get the atom names */
              strcpy(atom_name1,atom1_ptr->atom_name);
              strcpy(atom_name2,atom2_ptr->atom_name);

              /* Get bond order for this bond */
              bond_order = get_order(res_name,atom_name1,atom_name2);

              /* Store the bond order */
              bond_ptr->bond_order = bond_order;
            }
        }

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/* <--v.4.3 */
/***********************************************************************
      
bond_linkage  -  Determine which other bonds are linked to each bond's
                 two ends

***********************************************************************/

void bond_linkage(void)
{
  int link_to_first, link_to_second;

  struct bond *bond_ptr, *other_bond_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
    
  /* Loop through all the bonds in the linked list, excluding any elastic
     bonds */
  while (bond_ptr != NULL)
    {
      /* Loop through all other bonds to find the ones that have an atom
         in common with the current bond */
      other_bond_ptr = first_bond_ptr;
      
      /* Loop through all the bonds in the linked list, excluding current
         one */
      while (bond_ptr->elastic == FALSE && other_bond_ptr != NULL)
        {
          /* Check that this isn't the current bond */
          if (other_bond_ptr->elastic == FALSE && other_bond_ptr != bond_ptr)
            {
              /* Initialise flags identifying links to current bond */
              link_to_first = FALSE;
              link_to_second = FALSE;

              /* Check if other bond is linked to current bond's first
                 atom */
              if (bond_ptr->first_atom_ptr ==
                  other_bond_ptr->first_atom_ptr ||
                  bond_ptr->first_atom_ptr ==
                  other_bond_ptr->second_atom_ptr)
                link_to_first = TRUE;

              /* Check if other bond is linked to current bond's second
                 atom */
              if (bond_ptr->second_atom_ptr ==
                  other_bond_ptr->first_atom_ptr ||
                  bond_ptr->second_atom_ptr ==
                  other_bond_ptr->second_atom_ptr)
                link_to_second = TRUE;

              /* If the two bonds are linked, then create a link record */
              if (link_to_first == TRUE || link_to_second == TRUE)
                {
                  /* Allocate memory for structure to hold this bond link */
                  bond_link_ptr
                    = (struct bond_link*)malloc(sizeof(struct bond_link));
                  if (bond_link_ptr == NULL)
                    {
                      printf("*** ERROR. Unable to allocate memory for");
                      printf(" struct bond_link\n");
                      printf("***        Program ligplot terminated with");
                      printf(" error.\n");
                      exit (1);
                    }

                  /* Store the link and appropriate data */
                  bond_link_ptr->bond_ptr = other_bond_ptr;
                  if (link_to_first == TRUE)
                    {
                      bond_link_ptr->bond_end = FIRST;
                      bond_ptr->nfirst_atom_links++;
                    }
                  else
                    {
                      bond_link_ptr->bond_end = SECOND;
                      bond_ptr->nsecond_atom_links++;
                    }
                  bond_link_ptr->next_bond_link_ptr = NULL;
                  
                  /* If this is not the first link for current bond, make
                     new link point to previous link */
                  if (bond_ptr->first_bond_link_ptr != NULL)
                    bond_link_ptr->next_bond_link_ptr
                      = bond_ptr->first_bond_link_ptr;

                  /* Point from current bond to current link */
                  bond_ptr->first_bond_link_ptr = bond_link_ptr;
                }
            }

          /* Get pointer to next bond in linked-list */
          other_bond_ptr = other_bond_ptr->next_bond_ptr;
        }

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/* v.4.4--> */
/***********************************************************************
      
check_metal_bonds  -  Correct for "covalent bonds" to metal ions which
                      are likely to make plot difficult to unravel. Make
                      these bonds elastic

***********************************************************************/

void check_metal_bonds(void)
{
  char atom_name[5], element[3];
  int is_metal;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
    
  /* Loop through all covalent bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Process if this is a covelent bond */
      if (bond_ptr->bond_type == COVALENT)
        {
          /* Get the two residues involved */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;
          residue1_ptr = atom1_ptr->residue_ptr;
          residue2_ptr = atom2_ptr->residue_ptr;

          /* If residues aren't the same one, check further */
          if (residue1_ptr != residue2_ptr)
            {
              /* Initialise flag */
              residue_ptr = NULL;

              /* Check if first or second residue has only a single atom */
              if (residue1_ptr->natoms == 1)
                residue_ptr = residue1_ptr;
              else if (residue2_ptr->natoms == 1)
                residue_ptr = residue2_ptr;

              /* If either is a singleton, check whether it is also a
                 metal */
              if (residue_ptr != NULL)
                {
                  /* Get the atom */
                  atom_ptr = residue_ptr->first_atom_ptr;
                  strcpy(atom_name,atom_ptr->atom_name);
                  strcpy(element,atom_ptr->element);

                  /* Check whether this is a metal */
                  is_metal = check_for_metal(atom_name,element);

                  /* If it's a metal, then change the current bond
                     type to be elastic */
                  if (is_metal == TRUE ||
                      !strcmp(residue_ptr->res_name,"HOH"))
                    {
                      /* Show message that bond made elastic */
                      printf("*** Note. Bond to metal/water made elastic:");
                      printf(" %s %s %s %c",
                             atom1_ptr->atom_name,residue1_ptr->res_name,
                             residue1_ptr->res_num,residue1_ptr->chain);
                      printf(" -> %s %s %s %c\n",
                             atom2_ptr->atom_name,residue2_ptr->res_name,
                             residue2_ptr->res_num,residue2_ptr->chain);
                      bond_ptr->elastic = TRUE;
                      Nwarnings++;
                    }
                }
            }
        }

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/* <--v.4.4 */
/***********************************************************************

get_test_bond  -  Pick a representative bond for the given residue from
                  which to check connectivities to all atoms in the
                  residue

***********************************************************************/

struct bond *get_test_bond(struct residue *residue_ptr)
{
  int desirability, stored_desirability;

  struct bond *bond_ptr, *test_bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  test_bond_ptr = NULL;
  stored_desirability = 0;

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get the residues to which the bond's two atoms belong
         to */
      residue1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
      residue2_ptr = bond_ptr->second_atom_ptr->residue_ptr;

      /* Process only if this bond is a covalent bond entirely
         within the current residue */
      if (bond_ptr->bond_type == COVALENT &&
          residue1_ptr == residue_ptr &&
          residue2_ptr == residue_ptr)
        {
          /* Get the bond's two atoms */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Calculate this bond's desirability as the test
             bond */

          /* Best bond is N-CA bond */
          if ((!strncmp(atom1_ptr->atom_name," N  ",4) &&
               !strncmp(atom2_ptr->atom_name," CA ",4)) ||
              (!strncmp(atom2_ptr->atom_name," N  ",4) &&
               !strncmp(atom1_ptr->atom_name," CA ",4)))
            desirability = 14;

          /* Next-best is any mainchain bond */
          else if (atom1_ptr->side_chain == FALSE &&
                   atom2_ptr->side_chain == FALSE)
            {
              /* If rotatable, then better than not */
              if (bond_ptr->rotatable_bond == TRUE)
                desirability = 12;
              else
                desirability = 10;
            }

          /* Next any bond involving a mainchain atom */
          else if (atom1_ptr->side_chain == FALSE ||
                   atom2_ptr->side_chain == FALSE)
            {
              /* If rotatable, then better than not */
              if (bond_ptr->rotatable_bond == TRUE)
                desirability = 8;
              else
                desirability = 6;
            }

          /* Finally, any bond will do, preferably a
             rotatable one */
          else if (bond_ptr->rotatable_bond == TRUE)
            desirability = 4;
          else
            desirability = 2;

          /* If this bond is more desirable than the one we
             already have, then store it */
          if (desirability > stored_desirability)
            {
              /* Store the current bond */
              test_bond_ptr = bond_ptr;
              stored_desirability = desirability;
            }
        }

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return the test bond */
  return(test_bond_ptr);
}
/* v.4.1.2--> */
/***********************************************************************

list_bonds  -  List all the bonds in the structure

***********************************************************************/

void list_bonds(void)
{
  int ibond;
  char bond_type_char, ebond_char, rbond_char, single_double, source[7];
  struct bond *bond_ptr;
  struct coordinate *atom_ptr1, *atom_ptr2;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
  ibond = 0;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Get character for bond type */
      bond_type_char = '.';
      if (bond_ptr->bond_type == DELETED)
        bond_type_char = 'X';
      else if (bond_ptr->bond_type == HBOND)
        bond_type_char = 'H';
      else if (bond_ptr->bond_type == COVALENT)
        bond_type_char = 'c';
      else if (bond_ptr->bond_type == CONTACT)
        bond_type_char = 'n';
      else if (bond_ptr->bond_type == INTERNAL)
        bond_type_char = 'i';
/* v.5.4--> */
      else if (bond_ptr->bond_type == SALT_BRIDGE)
        bond_type_char = 'S';
      else if (bond_ptr->bond_type == DISULPHIDE)
        bond_type_char = 'D';
/* <--v.5.4 */

      /* Determine whether a rotatable bond */
      rbond_char = ' ';
      if (bond_ptr->rotatable_bond == TRUE)
        rbond_char = 'r';
      else if (bond_ptr->end_bond == TRUE)
        rbond_char = 'e';

      /* Determine whether an elastic bond */
      ebond_char = ' ';
      if (bond_ptr->elastic == TRUE)
        ebond_char = 'E';

      /* Get the bond's two atoms */
      atom_ptr1 = bond_ptr->first_atom_ptr;
      atom_ptr2 = bond_ptr->second_atom_ptr;

      /* Determine bond-order, single/double/treble */
      single_double = '-';
      if (bond_ptr->bond_order == DOUBLE)
        single_double = 'd';
      else if (bond_ptr->bond_order == TRIPLE)
        single_double = 't';

      /* Determine where bond info came from */
      if (bond_ptr->bond_source == CONECT)
        strcpy(source,"CONECT");
      else if (bond_ptr->bond_source == CALCULATED)
        strcpy(source,"Calc  ");
      else if (bond_ptr->bond_source == HBPLUS)
        strcpy(source,"HBPLUS");
      else
        strcpy(source,"      ");
          
      /* Write the bond information out */
      printf("Bond %3d. Atoms%s [%s %s %s %c]",ibond,
             atom_ptr1->atom_number,atom_ptr1->atom_name,
             atom_ptr1->residue_ptr->res_name,
             atom_ptr1->residue_ptr->res_num,
             atom_ptr1->residue_ptr->chain);
      printf(" ->%s [%s %s %s %c]",
             atom_ptr2->atom_number,atom_ptr2->atom_name,
             atom_ptr2->residue_ptr->res_name,
             atom_ptr2->residue_ptr->res_num,
             atom_ptr2->residue_ptr->chain);
      printf("%c%c%c%7.4f %c %s\n",bond_type_char,
             rbond_char,ebond_char,bond_ptr->bond_length,
             single_double,source);

      /* Increment bond counter */
      ibond++;

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/* <--v.4.1.2 */
/***********************************************************************

mark_downstream_atoms  -  Mark all the atoms downstream of the given bond

***********************************************************************/

void mark_downstream_atoms(struct bond *reference_bond_ptr,int end)
{
  int nbonds;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Mark both this bond's atoms as already checked */
  reference_bond_ptr->first_atom_ptr->checked = TRUE;
  reference_bond_ptr->second_atom_ptr->checked = TRUE;

  /* Set pointer to the first of the bonds sprouting off the given
     bond */
  bond_link_ptr = reference_bond_ptr->first_bond_link_ptr;

  /* Set the reference-bond as already checked */
  reference_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the reference-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
         end */
      if (bond_link_ptr->bond_end == end)
        {
          /* Get the bond this link is pointing to and add to stack */
          bond_ptr = bond_link_ptr->bond_ptr;

          /* Add to stack provided it is not an elastic bond */
          if (bond_ptr->elastic == FALSE)
            {
              /* If this is the first to be added to the stack, set both
                 stack-pointers to point to it */
              if (next_free_stack_ptr == NULL)
                {
                  next_free_stack_ptr = bond_ptr;
                  next_stack_ptr = bond_ptr;
                }

              /* Otherwise, make last bond on stack point to this one
                 and set this one as nopw the last on the stack */
              else
                {
                  next_free_stack_ptr->next_stack_ptr = bond_ptr;
                  next_free_stack_ptr = bond_ptr;
                }
            }

          /* Mark bond as added to the stack */
          bond_ptr->checked = TRUE;
        }
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Mark the bond's two atoms as checked */
      bond_ptr->first_atom_ptr->checked = TRUE;
      bond_ptr->second_atom_ptr->checked = TRUE;

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
         current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
        {
          /* Get the pointer to this bond */
          other_bond_ptr = bond_link_ptr->bond_ptr;

          /* If this bond hasn't already been processed, then move
             its atoms (if they need moving) and add the bond to the
             stack */
          if (other_bond_ptr->checked == FALSE &&
              other_bond_ptr->elastic == FALSE)
            {
              next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
              next_free_stack_ptr = other_bond_ptr;

              /* If the stack has run out, restart it at the bond
                 just entered */
              if (next_stack_ptr == NULL)
                next_stack_ptr = other_bond_ptr;

              /* Mark bond as added to the stack */
              other_bond_ptr->checked = TRUE;
            }

          /* Go to the next bond-link in the linked list */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }
    }
}
/***********************************************************************

remove_unconnected atoms  -  Find, and mark for deletion, all atoms in
                             the current residue that cannot be reached
                             from the given test-bond

***********************************************************************/

void remove_unconnected_atoms(struct bond *bond_ptr,
/* v.4.1.2--> */
/*                            struct residue *residue_ptr,int *marked) */
                              struct residue *residue_ptr,int *marked,
                              int debug)
/* <--v.4.1.2 */
{
  int iatom, natoms;
/* v.4.1.2--> */
  int nmarked;
/* <--v.4.1.2 */

  struct coordinate *atom_ptr;

  /* Initialise all atoms */
  initialise_atoms();
/* v.4.1.2--> */
  nmarked = 0;
/* <--v.4.1.2 */

  /* Mark all atoms downstream of each end of the given bonds */
  mark_downstream_atoms(bond_ptr,FIRST);
  mark_downstream_atoms(bond_ptr,SECOND);

  /* Get pointer to the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;

  /* Initialise atom count */
  iatom = 0;
  natoms = residue_ptr->natoms;

  /* Loop over all the residue's atoms checking which atoms in the
     current residue haven't been marked */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* Check which atoms in the current residue haven't been marked */
      if (atom_ptr->residue_ptr == residue_ptr &&
          atom_ptr->checked == FALSE)
        {
          /* Have an unconnected atom, so mark it for deletion */
          atom_ptr->deleted = TRUE;
          *marked = TRUE;
/* v.4.1.2--> */
          nmarked++;
/* <--v.4.1.2 */
        }

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }

/* v.4.1.2--> */
  /* If in debug mode, show residue being deleted */
  if (debug == TRUE && nmarked > 0)
    {
      printf("----> In remove_unconnected_atoms. Marked ");
       printf("%d unreachable atoms in residue ",nmarked);
      printf("[%s %s %c]\n",
             residue_ptr->res_name,
             residue_ptr->res_num,
             residue_ptr->chain);
      printf("       for deletion (from bond ");
      printf("[%s %s %s %c] -> [%s %s %s %c]\n",
             bond_ptr->first_atom_ptr->atom_name,
             bond_ptr->first_atom_ptr->residue_ptr->res_name,
             bond_ptr->first_atom_ptr->residue_ptr->res_num,
             bond_ptr->first_atom_ptr->residue_ptr->chain,
             bond_ptr->second_atom_ptr->atom_name,
             bond_ptr->second_atom_ptr->residue_ptr->res_name,
             bond_ptr->second_atom_ptr->residue_ptr->res_num,
             bond_ptr->second_atom_ptr->residue_ptr->chain);
    }
/* <--v.4.1.2 */

}
/***********************************************************************

delete_unwanted_atoms  -  Delete any atoms marked for deletion

***********************************************************************/

/* v.4.1.2--> */
/* void delete_unwanted_atoms(void) */
void delete_unwanted_atoms(int debug)
/* <--v.4.1.2 */
{
  int iatom, natoms;

  struct coordinate *atom_ptr, *last_atom_ptr;
  struct residue *residue_ptr;

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;
  last_atom_ptr = NULL;

  /* Loop over all the residues */
  while (residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
              
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, deleting any unwanted
         ones */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* If atom is marked for deletion, then make pointer of previous
             atom by-pass it */
          if (atom_ptr->deleted == TRUE)
            {
/* v.4.1.2--> */
              /* If in debug mode, show atom being deleted */
              if (debug == TRUE)
                printf("       * Deleted atom  [%4s %3s %5s %c]\n",
                       atom_ptr->atom_name,
                       atom_ptr->residue_ptr->res_name,
                       atom_ptr->residue_ptr->res_num,
                       atom_ptr->residue_ptr->chain);
/* <--v.4.1.2 */

              /* If this is the first atom, then set pointer to
                 next one to bypass it */
              if (atom_ptr == first_atom_ptr)
                first_atom_ptr = atom_ptr->next;

              /* Otherwise, set last pointer around current atom */
              else
                last_atom_ptr->next = atom_ptr->next;

              /* If this is the first atom of the current residue,
                 then alter the residue's pointer */
              if (atom_ptr == residue_ptr->first_atom_ptr)
                residue_ptr->first_atom_ptr = atom_ptr->next;

              /* Reduce count of this residue's atoms */
              residue_ptr->natoms--;
            }

          /* Otherwise, store pointer to this as the last undeleted
             atom encountered */
          else
            last_atom_ptr = atom_ptr;
      
          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }
}
/***********************************************************************

delete_unwanted_bonds  -  Delete any bonds involving deleted atoms

***********************************************************************/

void delete_unwanted_bonds(void)
{
  struct bond *bond_ptr, *last_bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;

  /* Initialise */
  last_bond_ptr = NULL;

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop over all the bonds */
  while (bond_ptr != NULL)
    {
      /* Get the pointers to this bond's two atoms */
      atom1_ptr = bond_ptr->first_atom_ptr;
      atom2_ptr = bond_ptr->second_atom_ptr;

      /* If either atom has been deleted, then want to delete this bond */
      if (atom1_ptr->deleted == TRUE || atom2_ptr->deleted == TRUE ||
           bond_ptr->bond_type == DELETED)
        {
          /* Mark the bond as deleted */
          bond_ptr->bond_type = DELETED;

          /* Remove the bond from the linked list */
         
          /* If this is the first bond, then set pointer to
             next one to bypass it */
          if (bond_ptr == first_bond_ptr)
            first_bond_ptr = bond_ptr->next_bond_ptr;

          /* Otherwise, set last pointer around current bond */
          else
            last_bond_ptr->next_bond_ptr = bond_ptr->next_bond_ptr;
        }

      /* Otherwise, store pointer to this as the last undeleted
         bond encountered */
      else
        last_bond_ptr = bond_ptr;

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

delete_unwanted_bond_links  -  Delete any unwanted bond links

***********************************************************************/

void delete_unwanted_bond_links(void)
{
  struct bond *bond_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr, *last_object_bond_ptr;

  /* Initialise */
  last_object_bond_ptr = NULL;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;
      
      /* Loop through all this object's bonds, one by one */
      while (object_bond_ptr != NULL)
        {
          /* Get the current bond */
          bond_ptr = object_bond_ptr->bond_ptr;

          /* If this bond has been deleted, then need to delete
             this link also */
          if (bond_ptr->bond_type == DELETED)
            {
              /* If this is the first link, then set pointer to
                 next one to bypass it */
              if (object_bond_ptr == object_ptr->first_object_bond_ptr)
                object_ptr->first_object_bond_ptr
                  = object_bond_ptr->next_object_bond_ptr;

              /* Otherwise, set last pointer around current bond */
              else
                last_object_bond_ptr->next_object_bond_ptr
                  = object_bond_ptr->next_object_bond_ptr;
            }

          /* Otherwise, store pointer to this as the last undeleted
             link encountered */
          else
            last_object_bond_ptr = object_bond_ptr;

          /* Get pointer to the next link for this object */
          object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
        }

      /* Get pointer to the next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

delete_unwanted_bond_connections  -  Delete any unwanted bond-to-bond
                                     connections used in chasing downstream
                                     atoms

***********************************************************************/

void delete_unwanted_bond_connections(void)
{
  struct bond *bond_ptr, *other_bond_ptr;
  struct bond_link *bond_link_ptr, *last_bond_link_ptr;

  /* Initialise */
  last_bond_link_ptr = NULL;

  /* Initialise pointer to the first stored bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all bonds */
  while (bond_ptr != NULL)
    {
      /* Get pointer to the first of the bonds springing off this one */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Loop through all this bond's adjacent bonds, one by one */
      while (bond_link_ptr != NULL)
        {
          /* Get the connected bond */
          other_bond_ptr = bond_link_ptr->bond_ptr;

          /* If this bond has been deleted, or has been made into
             an elastic bond, then need to delete this link */
          if (other_bond_ptr->bond_type == DELETED ||
              other_bond_ptr->elastic == TRUE)
            {

              /* If this is the first link, then set pointer to
                 next one to bypass it */
              if (bond_link_ptr == bond_ptr->first_bond_link_ptr)
                bond_ptr->first_bond_link_ptr
                  = bond_link_ptr->next_bond_link_ptr;

              /* Otherwise, set last pointer around current bond */
              else
                last_bond_link_ptr->next_bond_link_ptr
                  = bond_link_ptr->next_bond_link_ptr;
            }

          /* Otherwise, store pointer to this as the last undeleted
             link encountered */
          else
            last_bond_link_ptr = bond_link_ptr;

          /* Get pointer to the next link for this bond */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

delete_atoms_and_bonds  -  Delete any atoms marked for deletion and any
                           bonds and bond-links that involve them

***********************************************************************/

/* v.4.1.2--> */
/* void delete_atoms_and_bonds(void) */
void delete_atoms_and_bonds(int debug,char *call_from)
/* <--v.4.1.2 */
{
/* v.4.1.2--> */
  static int call_time = 0;

  /* Increment number of times this routine has been called and, if in
   debug mode, display */
  call_time++;
  if (debug == TRUE)
    printf("--> In delete_atoms_and_bonds (%2d) ... Called from %s.\n",
           call_time,call_from);
/* <--v.4.1.2 */

  /* Delete any unwanted atoms, previously marked for deletion */
/* v.4.1.2--> */
/*  delete_unwanted_atoms(); */
  delete_unwanted_atoms(debug);
/* <--v.4.1.2 */

  /* Delete any unwanted bonds */
  delete_unwanted_bonds();

  /* Delete any unwanted links */
  delete_unwanted_bond_links();

  /* Delete any unwanted bond-to-bond connections */
  delete_unwanted_bond_connections();
}
/***********************************************************************

fix_object  -  If all object's atoms have anchors, check if we can use
               these anchor coords as the object's starting position

***********************************************************************/

int fix_object(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
/* v.5.3.4--> */
//  int fixed, have_anchors, ok;
  int fixed, have_atom_anchors, have_res_anchor, ok;
/* <--v.5.3.4 */

  float dist_sqrd, x1, x2, y1, y2;
/* v.5.3.4--> */
  float shift_x, shift_y;
/* <--v.5.3.4 */

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct object_bond *object_bond_ptr;
  struct residue *residue_ptr;

  /* Initialise */
  fixed = FALSE;
  have_atom_anchors = FALSE;
  ok = FALSE;
/* v.5.3.4--> */
  shift_x = shift_y = 0;
/* <--v.5.3.4 */

  /* If object has anchors, check whether all its atoms have anchors */
  if (object_ptr->have_anchors == TRUE)
    {
      /* Assume all atoms have anchors */
      have_atom_anchors = TRUE;

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
/* v.5.3.4--> */
//      while (iresid < nresid && residue_ptr != NULL &&
//             have_atom_anchors == TRUE)
      while (iresid < nresid && residue_ptr != NULL)
/* <--v.5.3.4 */
        {
/* v.5.3.4--> */
          /* See if this residue has a CofM anchor */
          have_res_anchor = FALSE;
          if (residue_ptr->have_anchor == TRUE && nresid == 1)
            {
              /* Update this residue's boundaries */
              update_residue_boundaries(object_ptr);

              /* Calculate the shift to place mean position at the
                 anchor point */
              shift_x = residue_ptr->anchor_pstn_x
                - residue_ptr->flattened_mean_x;
              shift_y = residue_ptr->anchor_pstn_y
                - residue_ptr->flattened_mean_y;

              /* Set flags */
              have_res_anchor = TRUE;
              fixed = TRUE;
              object_ptr->fixed = TRUE;
            }
/* <--v.5.3.4 */

          /* Get pointer to this residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
              
          /* Initialise atom count */
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all the residue's atoms to see if any do not have
             anchor positions defined */
/* v.5.3.4--> */
//          while (iatom < natoms && atom_ptr != NULL &&
//                 have_atom_anchors == TRUE)
          while (iatom < natoms && atom_ptr != NULL)
/* <--v.5.3.4 */
            {
              /* If atom doesn't have an anchor, unset flag */
              if (atom_ptr->have_anchor == FALSE)
                have_atom_anchors = FALSE;

/* v.5.3.4--> */
              /* If we have a residue anchor, then shift this atom's
                 position */
              if (have_res_anchor == TRUE)
                {
                  atom_ptr->x = atom_ptr->x + shift_x;
                  atom_ptr->y = atom_ptr->y + shift_y;
                }
/* <--v.5.3.4 */

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }

          /* Get pointer to the next residue in the list */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }
    }
  
  /* If we have all atom anchors, check whether all resultant bond lengths
     seem OK */
  if (have_atom_anchors == TRUE)
    {
      /* Initialise flag */
      ok = TRUE;

      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;
      
      /* Loop through all object's bonds */
      while (object_bond_ptr != NULL && ok == TRUE)
        {
          /* Get the bond pointer */
          bond_ptr = object_bond_ptr->bond_ptr;

          /* Get pointers to the bond's two atoms */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Store the two atoms' coordinates */
          x1 = atom1_ptr->anchor_pstn_x;
          y1 = atom1_ptr->anchor_pstn_y;
          x2 = atom2_ptr->anchor_pstn_x;
          y2 = atom2_ptr->anchor_pstn_y;

          /* Calculate the bond length squared */
          dist_sqrd = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);

          /* If bond looks wrong, then don't use anchor positions */
          if (dist_sqrd < MIN_ANCHOR_LEN_SQRD ||
              dist_sqrd > MAX_ANCHOR_LEN_SQRD)
            {
              /* Unset flag */
              ok = FALSE;
            }

          /* Get pointer to this object's next bond entry */
          object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
        }
    }

  /* If OK to use anchor coords, apply them to the object's atom positions */
  if (ok == TRUE)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues to assign anchor positions
         to atom coords */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Get pointer to this residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
              
          /* Initialise atom count */
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all the residue's atoms to fix coords */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* Fix atom's position to its anchor coords */
              atom_ptr->x = atom_ptr->anchor_pstn_x;
              atom_ptr->y = atom_ptr->anchor_pstn_y;

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }

          /* Get pointer to the next residue in the list */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

/* v.5.0--> */
      /* Set flags that coords fixed */
      fixed = TRUE;
      object_ptr->fixed = TRUE;
/* <--v.5.0 */
/* v.5.3.4--> */
      object_ptr->coords_fixed = TRUE;
/* <--v.5.3.4 */
    }

  /* Return whether anchors used for fixing object */
  return(fixed);
}
/***********************************************************************

check_connections  -  Check that all atoms within each residue are
                      connected

***********************************************************************/

/* v.4.1.2--> */
/* void check_connections(void) */
void check_connections(int debug)
/* <--v.4.1.2 */
{
  int iresid, nresid;
  int iobject;
  int marked;

  struct bond *test_bond_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  marked = FALSE;

  /* Initialise pointer to the first stored object */
  printf("\nChecking for unconnected atoms ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects to find a representative bond */
  while (object_ptr != NULL)
    {
      /* Skip hydrophobic groups, waters and simple H-bonded groups */
      if (object_ptr->object_type != HYDROPHOBIC &&
          object_ptr->object_type != SIMPLE_HGROUP &&
          object_ptr->object_type != WATER)
        {
          /* Set the pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get a test bond for this residue */
              test_bond_ptr = get_test_bond(residue_ptr);

              /* If have a valid test bond then check that all atoms in
                 the current residue are reachable from it */
              if (test_bond_ptr != NULL)
                {
                  /* Mark any unwanted atoms for deletion */
                  remove_unconnected_atoms(test_bond_ptr,residue_ptr,
/* v.4.1.2--> */
/*                                         &marked); */
                                           &marked,debug);
/* <--v.4.1.2 */
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

  /* If any atoms have been marked for deletion, then remove them and
     their associated bonds, etc */
  if (marked == TRUE)
    {
      /* Delete any unwanted atoms */
/* v.4.1.2--> */
/*      delete_atoms_and_bonds(); */
      delete_atoms_and_bonds(debug,"check_connections");
/* <--v.4.1.2 */
    }
}
/***********************************************************************

check_for_attachment  -  Check whether the given residue is an attachment
                         to a previous residue in the given object

***********************************************************************/

int check_for_attachment(struct object *object_ptr,
                         struct residue *current_residue_ptr)
{
  int iresid, nresid;
  int connected;

  struct bond *bond_ptr;
  struct residue *residue_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise flag */
  connected = FALSE;

  /* Get number of residues making up this object */
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Get pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;

  /* Loop through all the residues until hit the current one */
  while (iresid < nresid && residue_ptr != NULL &&
         residue_ptr != current_residue_ptr)
    {
      /* Get pointer to first bond in the list */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds in the linked list */
      while (bond_ptr != NULL && connected == FALSE)
        {
          /* If this is a covalent bond, check whether
             it links the two residues */
          if (bond_ptr->bond_type == COVALENT &&
              bond_ptr->elastic == FALSE)
            {
              /* Get this bond's two residues */
              residue1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
              residue2_ptr = bond_ptr->second_atom_ptr->residue_ptr;

              /* If it spans the two residues in question,
                 then have a link between them */
              if ((residue1_ptr == residue_ptr &&
                   residue2_ptr == current_residue_ptr) ||
                  (residue2_ptr == residue_ptr &&
                   residue1_ptr == current_residue_ptr))
                connected = TRUE;
            }

          /* Get pointer to next atom in linked-list */
          bond_ptr = bond_ptr->next_bond_ptr;
        }

      /* Get pointer to the next residue in the list */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return flag indicating whether the given residue is an attachment */
  return(connected);
}
/***********************************************************************

check_for_cyclics  -  Check whether any objects are cyclic (ie residue
                      m is covalent bound to residue n where m and n are
                      not adjacent)

***********************************************************************/

void check_for_cyclics(void)
{
  int iobject;
  int cyclic;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first stored bond */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds, identifying any belonging
         to the current object */
      while (bond_ptr != NULL)
        {
          /* Get the two atoms on either end of this bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the objects containing these atoms via their respective
             residues */
          object1_ptr = atom1_ptr->residue_ptr->object_ptr;
          object2_ptr = atom2_ptr->residue_ptr->object_ptr;
          
          /* If both atoms belong to the current object, then check
             which residues the atoms belong to */
          if (object1_ptr == object_ptr && object2_ptr == object_ptr)
            {
              /* Get which residues the two atoms belong to */
              residue1_ptr = atom1_ptr->residue_ptr;
              residue2_ptr = atom2_ptr->residue_ptr;

              /* If the atoms belong to two different residues, then
                 check for cyclic peptides */
              if (bond_ptr->bond_type == COVALENT &&
                  bond_ptr->elastic == FALSE &&
                  residue1_ptr != residue2_ptr)
                {
                  cyclic = FALSE;

                  /* If residues are not adjacent, then looks like a
                     cyclic peptide */
                  if (residue1_ptr->next_residue_ptr != residue2_ptr &&
                      residue2_ptr->next_residue_ptr != residue1_ptr)
                    {
                      cyclic = TRUE;

                      /* If they are the same residue (ie one is an
                         attachment to the other) then not cyclic */
                      if (!strncmp(residue1_ptr->res_num,
                                  residue2_ptr->res_num,5) &&
                          residue1_ptr->chain == residue2_ptr->chain)
                        cyclic = FALSE;
                    }

                  /* If have what appears to be a cyclic peptide, then
                     need to make the bond elastic */
                  if (cyclic == TRUE)
                    {
                      /* Print warning message and make the bond elastic */
                      printf("*** Warning. Seems to be a cyclic peptide\n");
                      printf("***          Bond made elastic: %s %s %s %c",
                             atom1_ptr->atom_name,residue1_ptr->res_name,
                             residue1_ptr->res_num,residue1_ptr->chain);
                      printf(" ->  %s %s %s %c\n",
                             atom2_ptr->atom_name,residue2_ptr->res_name,
                             residue2_ptr->res_num,residue2_ptr->chain);
                      bond_ptr->elastic = TRUE;
                      Nwarnings++;
                    }
                }
            }

          /* Get pointer to point to next bond in the list */
          bond_ptr = bond_ptr->next_bond_ptr;
        }

      /* Get pointer to point to next object in the list and increment
         object count */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/* v.5.0--> */
/***********************************************************************

get_last_object  -  Loop through all the objects to get the last one

***********************************************************************/

struct object *get_last_object(struct object *first_object_ptr)
{
  struct object *object_ptr, *last_object_ptr;

  /* Initialise variables */
  last_object_ptr = NULL;

  /* Get pointer to first object */
  object_ptr = first_object_ptr;

  /* Loop through all objects to find last one */
  while (object_ptr != NULL)
    {
      /* Save this object's pointer as the last encountered */
      last_object_ptr = object_ptr;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Return the last object */
  return(last_object_ptr);
}
/***********************************************************************

get_last_residue  -  Loop through all the residues to get the last one

***********************************************************************/

struct residue *get_last_residue(struct residue *first_residue_ptr)
{
  struct residue *residue_ptr, *last_residue_ptr;

  /* Initialise variables */
  last_residue_ptr = NULL;

  /* Get pointer to first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all residues to find last one */
  while (residue_ptr != NULL)
    {
      /* Save this residue's pointer as the last encountered */
      last_residue_ptr = residue_ptr;

      /* Get pointer to next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Return the last residue */
  return(last_residue_ptr);
}
/***********************************************************************

get_last_atom  -  Loop through all the atoms to get the last one

***********************************************************************/

struct coordinate *get_last_atom(struct coordinate *first_atom_ptr)
{
  struct coordinate *atom_ptr, *last_atom_ptr;

  /* Initialise variables */
  last_atom_ptr = NULL;

  /* Get pointer to first atom */
  atom_ptr = first_atom_ptr;

  /* Loop through all atoms to find last one */
  while (atom_ptr != NULL)
    {
      /* Save this atom's pointer as the last encountered */
      last_atom_ptr = atom_ptr;

      /* Get pointer to next atom */
      atom_ptr = atom_ptr->next;
    }

  /* Return the last atom */
  return(last_atom_ptr);
}
/* <--v.5.0 */
/***********************************************************************

split_objects  -  Check for breaks between adjacent residues in objects
                  having more than one residue, and split into two, if
                  necessary

***********************************************************************/

void split_objects(int *nobjects,int *nligands)
{
  int iobject, iresid, nresid;
  int connected;

  struct bond *bond_ptr;
  struct object *last_object_ptr, *object_ptr, *new_object_ptr;
  struct residue *residue_ptr, *next_residue_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  *nligands = 0;

  /* Loop through all the objects to get the last one in case need
     to create any new ones in the routine */
/* v.5.0--> */
  last_object_ptr = get_last_object(first_object_ptr);
/* <--v.5.0 */

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;
  
  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get number of residues making up this object */
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* If this is a ligand, increment count */
      if (object_ptr->object_type == LIGAND)
        (*nligands)++;

      /* Get pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;

      /* Loop through all pairs of adjacent residues */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Update residue's pointer to its parent object */
          residue_ptr->object_ptr = object_ptr;

          /* Get pointer to the next residue */
          next_residue_ptr = residue_ptr->next_residue_ptr;

          /* If next residue exists, then check for covalent bonds
             common to both residues */
          if (next_residue_ptr != NULL && iresid < nresid - 1 &&
              object_ptr->object_type == LIGAND)
            {
              /* Update residue's pointer to its parent object */
              next_residue_ptr->object_ptr = object_ptr;

              /* Initialise flag */
              connected = FALSE;

              /* Get pointer to first bond in the list */
              bond_ptr = first_bond_ptr;

              /* Loop through all the bonds in the linked list */
              while (bond_ptr != NULL && connected == FALSE)
                {
                  /* If this is a covalent bond, check whether
                     it links the two residues */
                  if (bond_ptr->bond_type == COVALENT &&
                      bond_ptr->elastic == FALSE)
                    {
                      /* Get this bond's two residues */
                      residue1_ptr
                        = bond_ptr->first_atom_ptr->residue_ptr;
                      residue2_ptr
                        = bond_ptr->second_atom_ptr->residue_ptr;

                      /* If it spans the two residues in question,
                         then have a link between them */
                      if ((residue1_ptr == residue_ptr &&
                           residue2_ptr == next_residue_ptr) ||
                          (residue2_ptr == residue_ptr &&
                           residue1_ptr == next_residue_ptr))
                        connected = TRUE;
                    }

                  /* Get pointer to next atom in linked-list */
                  bond_ptr = bond_ptr->next_bond_ptr;
                }

              /* If have no link between the two residues, then may
                 be an attachment */
              if (connected == FALSE)
                {
                  /* Check whether the residue is an attachment
                     to a previous residue in this object */
                  connected = check_for_attachment(object_ptr,
                                                   next_residue_ptr);
                }

              /* If either residue is a heme, then split anyway */
              if (!strncmp(residue_ptr->res_name,"HEM",3) ||
                  !strncmp(next_residue_ptr->res_name,"HEM",3))
                connected = FALSE;

/* v.4.0--> */
              /* If plotting interactions across an interface, then want
                 all residues to be separate objects */
              if (Interface_Plot == TRUE)
                connected = FALSE;
/* <--v.4.0 */

              /* If still no link between the two residues, then may have
                 to split object into two */
              if (connected == FALSE)
                {
/* v.4.0--> */
                  if (Interface_Plot == FALSE)
                    {
/* <--v.4.0 */
                      printf("\n");
                      printf("*** Splitting ligand into two between");
                      printf("  %s %s %c and %s %s %c\n",
                             residue_ptr->res_name,
                             residue_ptr->res_num,
                             residue_ptr->chain,
                             next_residue_ptr->res_name,
                             next_residue_ptr->res_num,
                             next_residue_ptr->chain);
                      printf("\n");
                      Nwarnings++;
/* v.4.0--> */
                    }
/* <--v.4.0 */

                  /* Allocate memory for structure to hold new object's
                     details */
                  new_object_ptr = (struct object*)
                    malloc(sizeof(struct object));
                  if (new_object_ptr == NULL)
                    {
                      printf("*** ERROR. Unable to allocate memory for");
                      printf(" struct object\n");
                      printf("***        Program ligplot terminated with");
                      printf(" error.\n");
                      exit (1);
                    }

                  /* Fill in initial details for the new object */
                  new_object_ptr->object_type = object_ptr->object_type;
                  new_object_ptr->nresidues = nresid - iresid - 1;
                  new_object_ptr->first_residue_ptr = next_residue_ptr;
                  new_object_ptr->first_object_bond_ptr = NULL;
                  new_object_ptr->next_object_ptr = NULL;

                  /* Initialise remaining details (which will be
                     filled in later) */
                  new_object_ptr->nbonds = 0;
                  new_object_ptr->nrot_bonds = 0;
                  new_object_ptr->nmain_chain = 0;
                  new_object_ptr->n_ca = 0;
                  new_object_ptr->weight = 0;
                  new_object_ptr->minx = 0.0;
                  new_object_ptr->miny = 0.0;
                  new_object_ptr->maxx = 0.0;
                  new_object_ptr->maxy = 0.0;
/* v.4.0--> */
                  new_object_ptr->interface = 0;
/* <--v.4.0 */
                  new_object_ptr->max_atom_size = 0.0;
                  new_object_ptr->place_x = 0.0;
                  new_object_ptr->place_y = 0.0;
                  new_object_ptr->internal_energy = 0.0;
                  new_object_ptr->total_energy = 0.0;

                  /* Update number of residues now belonging to
                     the object that has been split */
                  object_ptr->nresidues = iresid + 1;

                  /* Add link from last object to the new one just
                     created */
                  last_object_ptr->next_object_ptr = new_object_ptr;
                  last_object_ptr = new_object_ptr;

                  /* Increment count of objects */
                  (*nobjects)++;

                  /* Finish with the current object */
                  iresid = nresid + 1;
                }
            }

          /* Get pointer to the next residue in the list */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

assign_bonds_to_objects  -  Determine which bonds relate to each object

***********************************************************************/

void assign_bonds_to_objects(int *ndeleted,int show_deletions)
{
  int iobject, iresid, nbonds, nresid;
  int ncovalent, nhbonds, ncontacts;
  int cyclic;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr, *next_object_ptr,
  *previous_object_ptr;
  struct object_bond *object_bond_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  previous_object_ptr = NULL;
  iobject = 0;
/* v.4.4--> */
/*  *ndeleted = 0; */
/* <--v.4.4 */

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Initialise counts of bonds between this object and all others */
      nbonds = 0;
      ncovalent = 0;
      nhbonds = 0;
      ncontacts = 0;

      /* Save pointer to the next object */
      next_object_ptr = object_ptr->next_object_ptr;

      /* Get pointer to the first stored bond */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds, identifying any belonging
         to the current object */
      while (bond_ptr != NULL)
        {
          /* Get the two atoms on either end of this bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the objects containing these atoms via their respective
             residues */
          object1_ptr = atom1_ptr->residue_ptr->object_ptr;
          object2_ptr = atom2_ptr->residue_ptr->object_ptr;

          /* If both atoms belong to the current object, then assign
             current bond, if covalent, to current object */
          if (object1_ptr == object_ptr && object2_ptr == object_ptr)
            {
              /* Get which residues the two atoms belong to */
              residue1_ptr = atom1_ptr->residue_ptr;
              residue2_ptr = atom2_ptr->residue_ptr;

              /* If the atoms belong to two different residues, then
                 check for cyclic peptides */
              if (bond_ptr->bond_type == COVALENT &&
                  residue1_ptr != residue2_ptr)
                {
                  cyclic = FALSE;

                  /* If residues are not adjacent, then looks like a
                     cyclic peptide */
                  if (residue1_ptr->next_residue_ptr != residue2_ptr &&
                      residue2_ptr->next_residue_ptr != residue1_ptr)
                    {
                      cyclic = TRUE;

                      /* If they are the same residue (ie one is an
                         attachment to the other) then not cyclic */
                      if (!strncmp(residue1_ptr->res_num,
                                  residue2_ptr->res_num,5) &&
                          residue1_ptr->chain == residue2_ptr->chain)
                        cyclic = FALSE;
                    }

                  /* If have what appears to be a cyclic peptide, then
                     need to make the bond elastic */
                  if (cyclic == TRUE)
                    {
                      /* Print warning message and make the bond elastic */
                      printf("*** Warning. Seems to be a cyclic peptide\n");
                      printf("***          Bond made elastic: %s %s %s %c",
                             atom1_ptr->atom_name,residue1_ptr->res_name,
                             residue1_ptr->res_num,residue1_ptr->chain);
                      printf(" ->  %s %s %s %c\n",
                             atom2_ptr->atom_name,residue2_ptr->res_name,
                             residue2_ptr->res_num,residue2_ptr->chain);
                      bond_ptr->elastic = TRUE;
                      Nwarnings++;
                    }
                }

              /* If bond is a covalent bond then perform assignment */
              if (bond_ptr->bond_type == COVALENT &&
                  bond_ptr->elastic == FALSE)
                {
                  /* Increment count of covalent bonds belonging to
                     the current object */
                  object_ptr->nbonds++;

                  /* Allocate memory for structure to hold this
                     object-bond link */
                  object_bond_ptr
                    = (struct object_bond*)malloc(sizeof(struct object_bond));
                  if (object_bond_ptr == NULL)
                    {
                      printf("*** ERROR. Unable to allocate memory for");
                      printf(" struct object_bond\n");
                      printf("***        Program ligplot terminated with");
                      printf(" error.\n");
                      exit (1);
                    }

                  /* Store the link and appropriate data */
                  object_bond_ptr->bond_ptr = bond_ptr;
                  object_bond_ptr->next_object_bond_ptr = NULL;
                  
                  /* If this is not the first link for current bond, make
                     new link point to previous link */
                  if (object_ptr->first_object_bond_ptr != NULL)
                    object_bond_ptr->next_object_bond_ptr
                      = object_ptr->first_object_bond_ptr;
              
                  /* Point from current bond to current link */
                  object_ptr->first_object_bond_ptr = object_bond_ptr;
                }

              /* If this is a hydrogen bond within the object, mark it
                 as internal */
/* v.5.4--> */
//              else if (bond_ptr->bond_type == HBOND)
              else if (bond_ptr->bond_type == HBOND ||
                       bond_ptr->bond_type == SALT_BRIDGE ||
                       bond_ptr->bond_type == DISULPHIDE)
/* <--v.5.4 */
                bond_ptr->bond_type = INTERNAL;

              /* Otherwise, if it is an internal non-bonded contact
                 interaction, then delete it altogether */
              else if (bond_ptr->bond_type != COVALENT)
                bond_ptr->bond_type = DELETED;
            }

          /* Otherwise, if only one of the ends of the bond belongs to
             this object, need to make the bond an elastic one */
          else if (object1_ptr == object_ptr || object2_ptr == object_ptr)
            {
              bond_ptr->elastic = TRUE;

              /* Increment counts of different bond types between this
                 object and any others */
              nbonds++;
/* v.5.4--> */
//              if (bond_ptr->bond_type == HBOND)
              if (bond_ptr->bond_type == HBOND ||
                  bond_ptr->bond_type == SALT_BRIDGE ||
                  bond_ptr->bond_type == DISULPHIDE)
/* <--v.5.4 */
                {
                  nhbonds++;

                  /* If both objects are ligand objects, then make this an
                     internal hydrogen bond */
                  if (object1_ptr->object_type == LIGAND &&
                      object2_ptr->object_type == LIGAND)
                    bond_ptr->bond_type = INTERNAL;
                }
              else if (bond_ptr->bond_type == CONTACT)
                ncontacts++;
              else if (bond_ptr->bond_type == COVALENT)
                {
                  /* If this covalent bond is between the current object
                     and the ligand, then increment count of covalent
                     bonds */
                  if ((object1_ptr == object_ptr &&
                       object2_ptr->object_type == LIGAND) ||
                      (object2_ptr == object_ptr &&
                       object1_ptr->object_type == LIGAND))
                    ncovalent++;
                }
            }

          /* Get pointer to point to next bond in the list */
          bond_ptr = bond_ptr->next_bond_ptr;
        }

      /* If this is not the ligand object, then determine what type it
         is */
/* v.4.0--> */
/*      if (object_ptr->object_type != LIGAND) */
      if (object_ptr->object_type != LIGAND || Interface_Plot == TRUE)
/* <--v.4.0 */
        {
          /* If there are no bonds whatsoever between this object and
             any others, then might as well delete it */
          if (nbonds == 0)
            {
              /* Delete the current object */
/* v.4.0--> */
              if (Print_as_is == FALSE)
                {
/* <--v.4.0 */
                  delete_object(object_ptr,previous_object_ptr,show_deletions);
                  (*ndeleted)++;
/* v.4.0--> */
                }
/* <--v.4.0 */

              /* Set current pointer to the previous object */
              object_ptr = previous_object_ptr;
            }

          /* Otherwise, if no hydrogen bonds or covalent bonds, then
             mark this object as being a hydrophobic group */
          else if (ncontacts > 0 && nhbonds == 0 && ncovalent == 0)
/* v.4.0--> */
            {
              /* Make object hydrophobic, or set interface flag */
              if (Interface_Plot == FALSE)
/* <--v.4.0 */
                object_ptr->object_type = HYDROPHOBIC;
/* v.4.0--> */
              else
                object_ptr->interface = 9;
            }
/* <--v.4.0 */
        }

      /* Save the current object pointer */
      previous_object_ptr = object_ptr;

      /* Get pointer to point to next object in the list and increment
         object count */
      object_ptr = next_object_ptr;
      iobject++;
    }

/* v.4.5.1--> */
  /* If printing as is, and using simplified representation, identify
     residues for simple representation */
  if (Print_as_is == TRUE)
    {
      /* Reinitialise to first object */
      object_ptr = first_object_ptr;

      /* Loop through all objects */
      while (object_ptr != NULL)
        {
          /* If ligand residue, and using simple representation of
             ligands, rassign object type */
          if (object_ptr->object_type == LIGAND &&
              Include->Simple_Ligand_Residues == TRUE)
            {
              /* Set the pointer to the first of this object's residues */
              residue_ptr = object_ptr->first_residue_ptr;
              nresid = object_ptr->nresidues;
              iresid = 0;

              /* Loop over all this object's residues */
              while (iresid < nresid && residue_ptr != NULL)
                {
                  /* Mark the residue as being of the simplified type */
                  residue_ptr->residue_type = SIMPLE_LIGAND;

                  /* Get pointer to the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }
            }

          /* Ditto for non-ligand residues */
          if (object_ptr->object_type == HGROUP &&
              Include->Simple_Nonligand_Residues == TRUE)
            object_ptr->object_type = SIMPLE_HGROUP;

          /* Get pointer to point to next object in the list */
          object_ptr = object_ptr->next_object_ptr;
        }
    }
/* <--v.4.5.1 */

/* v.4.0--> */
  /* For interface plots, loop through the object again setting the
     appropriate ones as hydrophobic objects */
  if (Interface_Plot == TRUE)
    {
      /* Reinitialise to first object */
      object_ptr = first_object_ptr;

      /* Loop through all objects */
      while (object_ptr != NULL)
        {
          /* If this should be a hydrophobic object, then make it so */
          if (object_ptr->interface == 9)
            {
              /* Mark the appropriate interface */
              if (object_ptr->object_type == LIGAND)
                object_ptr->interface = 1;
              else
                object_ptr->interface = 2;
              object_ptr->object_type = HYDROPHOBIC;
            }

          /* Otherwise set interface type according to the object type */
          else
            {
              /* Mark the appropriate interface */
              if (object_ptr->object_type == LIGAND)
                object_ptr->interface = 1;
              else if (object_ptr->object_type != WATER)
                object_ptr->interface = 2;
/* v.5.3.5--> */
              if (object_ptr->object_type != WATER)
                object_ptr->object_type = HGROUP;
/* <--v.5.3.5 */
            }

          /* Get pointer to point to next object in the list */
          object_ptr = object_ptr->next_object_ptr;
        }
    }
/* <--v.4.0 */
}
/***********************************************************************

show_deleted_residues_message  -  Show explanatory note about residues
                                  deleted

***********************************************************************/

void show_deleted_residues_message(int ndeleted)
{
  /* Print the message */
  printf("\n");
  if (ndeleted == 1)
    {
      printf("NOTE. The above residue has been deleted as it is not ");
      printf("connected to \n");
      printf("any otherin the diagram. It may have originally been ");
      printf("read in because it \n");
      printf("appears in the CONECT records of the PDB file.\n");
    }
  else
    {
      printf("NOTE. The above residues have been deleted as they ");
      printf("are not connected to \n");
      printf("any others in the diagram. They may have originally ");
      printf("been read in because\n");
      printf("they appear in the CONECT records of the PDB file.\n");
    }
}
/* v.5.1--> */
/***********************************************************************

mark_unwanted_waters  -  Mark any waters with fewer than 2 H-bonds to
                         protein/ligand

***********************************************************************/

void mark_unwanted_waters(void)
{
  int iatom, iresid, natoms, nbonds, nresid;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr, *other_atom_ptr;
  struct bond *bond_ptr;
  struct object *object_ptr, *other_object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects to mark the ligand atoms as checked */
  while (object_ptr != NULL)
    {
      /* If the object is a water, check all H-bonds that involve it */
      if (object_ptr->object_type == WATER)
        {
          /* Initialise count of significant bonds */
          nbonds = 0;

          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              
              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Get pointer to the first stored bond */
                  bond_ptr = first_bond_ptr;

                  /* Loop through all the bonds, to find any H-bonds
                     involving this water */
                  while (bond_ptr != NULL)
                    {
                      /* Get the two atoms on either end of this bond */
                      atom1_ptr = bond_ptr->first_atom_ptr;
                      atom2_ptr = bond_ptr->second_atom_ptr;

                      /* If either one is our water, get the other atom */
                      other_atom_ptr = NULL;
                      if (atom_ptr == atom1_ptr)
                        other_atom_ptr = atom2_ptr;
                      else if (atom_ptr == atom2_ptr)
                        other_atom_ptr = atom1_ptr;

                      /* Determine what the object at the other end of
                         the bond is */
/* v.5.4--> */
//                      if (bond_ptr->bond_type == HBOND &&
                      if ((bond_ptr->bond_type == HBOND ||
                           bond_ptr->bond_type == SALT_BRIDGE ||
                           bond_ptr->bond_type == DISULPHIDE) &&
/* <--v.5.4 */
                          other_atom_ptr != NULL)
                        {
                          /* Get the other atom's object */
                          other_object_ptr
                            = other_atom_ptr->residue_ptr->object_ptr;

                          /* If object is protein or ligand, increment
                             count of significant atoms */
                           if (other_object_ptr->object_type == LIGAND ||
                               other_object_ptr->object_type == HGROUP)
                             nbonds++;
                        }

                      /* Get pointer to point to next bond in the list */
                      bond_ptr = bond_ptr->next_bond_ptr;
                    }

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* If this water doesn't have at least 2 bonds to protein/ligand
             then mark it for deletion */
          if (nbonds < 2 && object_ptr->first_residue_ptr)
            object_ptr->first_residue_ptr->deleted = TRUE;
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* <--v.5.1 */
/***********************************************************************

mark_reachable_atoms  -  Mark all atoms that are reachable from the
                         ligand residue(s) via covalent bonds, H-bonds
                         or hydrophobic contacts

***********************************************************************/

void mark_reachable_atoms(void)
{
  int iatom, iresid, natoms, nresid;
  int checked;
  int wanted;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue_ptr;

  /* Inialise all atoms */
  initialise_atoms();

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects to mark the ligand atoms as checked */
  while (object_ptr != NULL)
    {
      /* If the object is a ligand, mark all its atoms */
      if (object_ptr->object_type == LIGAND)
        {
          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              
              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Mark this atom as wanted */
                  atom_ptr->checked = TRUE;

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Initialise flag */
  checked = TRUE;

  /* Loop until no new atoms checked as being connected */
  while (checked == TRUE)
    {
      /* reset flag */
      checked = FALSE;

      /* Get pointer to the first stored bond */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds, marking any unchecked atoms 
         connected to checked atoms */
      while (bond_ptr != NULL)
        {
          /* Get the two atoms on either end of this bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the two objects these atoms belong to */
          object1_ptr = atom1_ptr->residue_ptr->object_ptr;
          object2_ptr = atom2_ptr->residue_ptr->object_ptr;

          /* If the link is to a HYDROPHOBIC group, then ignore this
             bond */
          wanted = TRUE;
          if ((object1_ptr->object_type == HYDROPHOBIC &&
               object2_ptr->object_type != LIGAND) ||
              (object2_ptr->object_type == HYDROPHOBIC &&
               object1_ptr->object_type != LIGAND))
            wanted = FALSE;

          /* Process only if this is a valid bond */
          if (wanted == TRUE && (bond_ptr->bond_type == COVALENT ||
              bond_ptr->bond_type == HBOND || 
/* v.5.4--> */
              bond_ptr->bond_type == SALT_BRIDGE || 
              bond_ptr->bond_type == DISULPHIDE || 
/* <--v.5.4 */
              bond_ptr->bond_type == CONTACT))
            {
              /* If one is checked and other not, then mark the one that
                 isn't as reachable */
              if (atom1_ptr->checked == TRUE &&
                  atom2_ptr->checked == FALSE)
                {
                  atom2_ptr->checked = TRUE;
                  checked = TRUE;
                }
              else if (atom2_ptr->checked == TRUE &&
                       atom1_ptr->checked == FALSE)
                {
                  atom1_ptr->checked = TRUE;
                  checked = TRUE;
                }
            }

          /* Get pointer to point to next bond in the list */
          bond_ptr = bond_ptr->next_bond_ptr;
        }
    }
}
/***********************************************************************

delete_unreachable_residues  -  Delete any unreachable residues

***********************************************************************/

void delete_unreachable_residues(int *ndeleted,int show_deletions)
{
  int iatom, iobject, iresid, natoms, nresid;
  int got_atom;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct object *next_object_ptr, *previous_object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  previous_object_ptr = NULL;
  iobject = 0;

  /* Loop through all objects to check which are reachable */
  while (object_ptr != NULL)
    {
      /* Save pointer to the next object */
      next_object_ptr = object_ptr->next_object_ptr;

      /* Only want to look at non-ligand objects */
      if (object_ptr->object_type != LIGAND)
        {
          /* Set flag indicating whether any reachable atoms located */
          got_atom = FALSE;

          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL &&
                 got_atom == FALSE)
            {
              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              
              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL &&
                     got_atom == FALSE)
                {
                  /* If atom is reachable, then set flag */
                  if (atom_ptr->checked == TRUE)
                    got_atom = TRUE;

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

/* v.5.1--> */
          /* If this is a water that has been marked for deletion,
             the flag object for deletion */
          if (object_ptr->object_type == WATER &&
              object_ptr->first_residue_ptr != NULL &&
              object_ptr->first_residue_ptr->deleted == TRUE)
            got_atom = FALSE;
/* <--v.5.1 */

          /* If object is unreachable, then need to delete it and all its
             residues */
          if (got_atom == FALSE)
            {
              /* Delete the current object */
              delete_object(object_ptr,previous_object_ptr,show_deletions);
              (*ndeleted)++;

              /* Set current pointer to the previous object */
              object_ptr = previous_object_ptr;
            }
        }

      /* Save the current object pointer */
      previous_object_ptr = object_ptr;

      /* Get pointer to point to next object in the list and increment
         object count */
      object_ptr = next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

get_atom_sizes  -  Get all the atom-sizes, according to the type of
                   object the atoms belong to

***********************************************************************/

void get_atom_sizes(void)
{
  int object_type;

  float atom_size;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;
  struct object *object_ptr;

  /* Initialise all the maximum object- and residue- atom sizes */
  
  /* Get pointer to the first object */
  object_ptr = first_object_ptr;

  /* Loop through all the stored objects to initialise */
  while (object_ptr != NULL)
    {
      /* Initialise */
      object_ptr->max_atom_size = 0.0;

      /* Get pointer to the next object in the linked list */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Get pointer to the first residue */
  residue_ptr = first_residue_ptr;

  /* Loop through all the stored residues to initialise */
  while (residue_ptr != NULL)
    {
      /* Initialise */
      residue_ptr->max_atom_size = 0.0;

      /* Get pointer to the next residue in the linked list */
      residue_ptr = residue_ptr->next_residue_ptr;
    }

  /* Get pointer to the first atom */
  atom_ptr = first_atom_ptr;

  /* Loop over all the atoms */
  while (atom_ptr != NULL)
    {
      /* Get the object and residue to which the atom belongs */
      residue_ptr = atom_ptr->residue_ptr;
      object_ptr = residue_ptr->object_ptr;

      /* Get the object type */
      object_type = object_ptr->object_type;

      /* Initialise size */
      atom_size = 0.0;

      /* Water atom */
      if (object_type == WATER)
/* v.4.0--> */
/*      if (Include->Nonligand_Atoms == TRUE) */
        if (Include->Water_Atoms == TRUE)
/* <--v.4.0 */
          atom_size = Size_Val->Waters;
        else
          atom_size = Size_Val->Nonligand_Bonds;

      /* Hydrophobic contact residue */
      else if (object_type == HYDROPHOBIC)
        atom_size = SPOKE_EXTENT * Size_Val->Hydrophobics;

      /* Simplified ligand representation */
      else if (residue_ptr->residue_type == SIMPLE_LIGAND &&
               !strncmp(atom_ptr->atom_name," CA ",4))
        atom_size = Size_Val->Simple_Residues;

      /* Ligand atom */
      else if (object_type == LIGAND)
        if (Include->Ligand_Atoms == TRUE)
          atom_size = Size_Val->Ligand_Atoms;
        else
          atom_size = Size_Val->Ligand_Bonds;

      /* Non-ligand atom */
      else
        if (Include->Nonligand_Atoms == TRUE)
          atom_size = Size_Val->Nonligand_Atoms;
        else
          atom_size = Size_Val->Nonligand_Bonds;

      /* Store the atom's size */
      atom_ptr->atom_size = atom_size;

      /* Update maximum atom size for this atom's residue and object */
      if (atom_size > object_ptr->max_atom_size)
        object_ptr->max_atom_size = atom_size;
      if (atom_size > residue_ptr->max_atom_size)
        residue_ptr->max_atom_size = atom_size;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************
      
atom_linkage  -  Determine which atoms are covalently bonded to each other

***********************************************************************/

void atom_linkage(void)
{
  int loop;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *other_atom_ptr;
  struct atom_link *atom_link_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
    
  /* Loop through all the bonds in the linked list, excluding any elastic
     bonds */
  while (bond_ptr != NULL)
    {
/* v.5.1.1--> */
      if (bond_ptr->bond_type == COVALENT)
        {
/* <--v.5.1.1 */
          /* Loop over this bond's two atoms */
          for (loop = 0; loop < 2; loop++)
            {
              if (loop == 0)
                {
                  atom_ptr = bond_ptr->first_atom_ptr;
                  other_atom_ptr = bond_ptr->second_atom_ptr;
                }
              else
                {
                  atom_ptr = bond_ptr->second_atom_ptr;
                  other_atom_ptr = bond_ptr->first_atom_ptr;
                }

              /* Add the other atom to the current atom's list of covalent
                 partners */

              /* Allocate memory for structure to hold this atom link */
              atom_link_ptr
                = (struct atom_link*)malloc(sizeof(struct atom_link));
              if (atom_link_ptr == NULL)
                {
                  printf("*** ERROR. Unable to allocate memory for");
                  printf(" struct atom_link\n");
                  printf("***        Program ligplot terminated with");
                  printf(" error.\n");
                  exit (1);
                }

              /* Store the link, pointing at the other atom */
              atom_link_ptr->atom_ptr = other_atom_ptr;
              atom_link_ptr->next_atom_link_ptr = NULL;

              /* Increment count of this atom's partners */
              atom_ptr->natom_links++;
                  
              /* If this is not the first link for current atom, make
                 new link point to previous link */
              if (atom_ptr->first_atom_link_ptr != NULL)
                atom_link_ptr->next_atom_link_ptr
                  = atom_ptr->first_atom_link_ptr;

              /* Point from current bond to current link */
              atom_ptr->first_atom_link_ptr = atom_link_ptr;
            }
/* v.5.1.1--> */
        }
/* <--v.5.1.1 */

      /* Get pointer to next bond in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

ring_check  -  Checks whether a given bond is in a ring - a path traced
               from one end of the bond eventually returns to the bond
               via the other end

***********************************************************************/

int ring_check(struct bond *test_bond_ptr)
{
  int got_ring;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;

  /* Initialise variables */
  got_ring = FALSE;
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Get the pointer to the first of the bonds coming off the 
     test_bond */
  bond_link_ptr = test_bond_ptr->first_bond_link_ptr;

  /* Initialise stack pointers for ring-search by putting the bonds
     coming off one end of the test-bond onto the stack of bonds to
     be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the first
         atom of the test-bond */
      if (bond_link_ptr->bond_end == FIRST)
        {
          /* Get the bond this link is pointing to and add to stack */
          bond_ptr = bond_link_ptr->bond_ptr;

          /* If this is the first to be added to the stack, set both
             stack-pointers to point to it */
          if (next_free_stack_ptr == NULL)
            {
              next_free_stack_ptr = bond_ptr;
              next_stack_ptr = bond_ptr;
            }

          /* Otherwise, make last bond on stack point to this one
             and set this one as now the last on the stack */
          else
            {
              next_free_stack_ptr->next_stack_ptr = bond_ptr;
              next_free_stack_ptr = bond_ptr;
            }

          /* Mark this bond as checked */
          bond_ptr->checked = TRUE;
        }
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted
     or have determined that we have a ring */
  while (next_stack_ptr != NULL && got_ring == FALSE)
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
          
          /* If this is the test-bond, check whether it is connected
             to the current bond via the atom opposite to the one we
             started with, implying we have gone round in a ring */
          if (other_bond_ptr == test_bond_ptr)
            {
              /* Check if the current bond has the atom we're after */
              if (bond_ptr->first_atom_ptr
                  == test_bond_ptr->second_atom_ptr ||
                  bond_ptr->second_atom_ptr
                  == test_bond_ptr->second_atom_ptr)
                {
                  /* Have gone round, so must have a ring */
                  got_ring = TRUE;
                }
            }

          /* Otherwise, if this bond hasn't already been processed
             or added to the stack, add it now */
          else if (other_bond_ptr->checked == FALSE)
            {
              next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
              next_free_stack_ptr = other_bond_ptr;

              /* If the stack has run out, restart it at the bond
                 just entered */
              if (next_stack_ptr == NULL)
                next_stack_ptr = other_bond_ptr;

              /* Mark this bond as checked */
              other_bond_ptr->checked = TRUE;
            }

          /* Go to the next bond-link in the linked list */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }
    }
    
  /* Return the flag showing whether a ring has been found or not */
  return(got_ring);
}
/***********************************************************************

identify_rotatable_bonds  -  Determine which bonds are rotatable in
                             each object

***********************************************************************/

void identify_rotatable_bonds(void)
{
  int in_ring, iobject;

  struct bond *bond_ptr;
  struct object *object_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise pointer to the first stored object */
  printf("\nIdentifying rotatable bonds ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;

      /* Loop through all this object's bonds, one by one */
      while (object_bond_ptr != NULL)
        {
          /* Get the current bond */
          bond_ptr = object_bond_ptr->bond_ptr;

          /* To be a rotatable bond, has to have other bonds linked to
             it at either end */
          if (bond_ptr->nfirst_atom_links > 0 &&
              bond_ptr->nsecond_atom_links > 0)
            {
              /* Check whether the current bond is in some sort of ring
                 - ie can start at one end and trace through the bonds to
                 the other */
              in_ring = ring_check(bond_ptr);

              /* If bond is a ring-bond, then set marker */
              if (in_ring == TRUE)
/* v.5.1--> */
                {
/* <--v.5.1 */
                  bond_ptr->ring_bond = TRUE;
/* v.5.1--> */
                  /* Mark this two bond's atoms as ring atoms */
                  bond_ptr->first_atom_ptr->ring_atom = TRUE;
                  bond_ptr->second_atom_ptr->ring_atom = TRUE;
                }
/* <--v.5.1 */

              /* Otherwise, must be a rotatable bond */
              else
                {
                  bond_ptr->rotatable_bond = TRUE;

                  /* Increment count of rotatable bonds in the current
                     object */
                  object_ptr->nrot_bonds++;
                }
            }

            /* If the current bond has other bonds linked only to one
               of its ends, then mark it as an "end" bond */
          else if (bond_ptr->nfirst_atom_links > 0 ||
                   bond_ptr->nsecond_atom_links > 0)
            bond_ptr->end_bond = TRUE;
          
          /* Get pointer to this object's next bond entry */
          object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
        }

      /* If object has bonds, but none are rotatable, make the first
         one a rotatable one */
      if (object_ptr->nbonds > 0 && object_ptr->nrot_bonds == 0)
        {
          bond_ptr = object_ptr->first_object_bond_ptr->bond_ptr;
          bond_ptr->rotatable_bond = TRUE;
          object_ptr->nrot_bonds++;
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

write_frm_file  -  Read through the PDB file, writing out all the
                   residues (complete atom records) that will be used
                   in the final LIGPLOT picture

***********************************************************************/

/* v.4.2--> */
/* void write_frm_file(char pdb_name[FILENAME_LEN]) */
void write_frm_file(char pdb_name[FILENAME_LEN],XMAS *xmas)
/* <--v.4.2 */
{
  char line[LINELEN + 1];
  char res_number[6], residue[4], last_residue[4], last_resnum[6];
  char atom_name[5], chain, last_chain;
/* v.4.2--> */
  char *eofstatus;
/* <--v.4.2 */

  int keep_reading;
  int cpos, ipos, want_residue;
  int atom_counter;
  int got_residue;
/* v.4.2--> */
  int rec_types[NTYPES], type;
/* <--v.4.2 */
/* v.4.2.2--> */
  int mol_id;
/* <--v.4.2.2 */

  struct residue *residue_ptr;
    
/* v.4.2--> */
/*  FILE *fa; */
  FILE *fil_pdb;
/* <--v.4.2 */

  /* Open PDB file */
/* v.4.2--> */
/*  if ((fa = fopen(pdb_name,"r")) ==  NULL) */
  if (Read_Xmas == FALSE)
    {
      if ((fil_pdb = fopen(pdb_name,"r")) ==  NULL)
/* <--v.4.2 */
        {
          printf("\n*** Unable to open your PDB file [%s]\n",pdb_name);
          exit(1);
        }
/* v.4.2--> */
    }
/* <--v.4.2 */
    
  /* Initialise variables */
  atom_counter = 0;
  keep_reading = TRUE;
  last_residue[0] = '\0';
  last_resnum[0] = '\0';
  last_chain = '\0';
  want_residue = FALSE;

  /* Read in the data from the PDB file */
  printf("\nWriting out the ligplot.frm file ...\n");

/* v.4.2--> */
  /* Get the first line from the PDB or XMAS file */
  if (Read_Xmas == FALSE)
    eofstatus = fgets(line,LINELEN,fil_pdb);
  else
    {
      /* Form the record types to be retrieved */
      rec_types[REC_ATOM] = TRUE;
      rec_types[REC_CONECT] = FALSE;
      type = REC_ATOM;

      /* Get first record from the XMAS cache */
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type); */
      eofstatus = get_xmas2pdb_line(xmas,line,TRUE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
    }
/* <--v.4.2 */

  /* Loop through the PDB file, picking up all the atoms that
     are required */
/* v.4.2--> */
/*  while (fgets(line,LINELEN,fa)!= NULL) */
  while (eofstatus != NULL)
/* <--v.4.2 */
    {
      /* Check for the end of an NMR structure model */
      if (!strncmp("ENDMDL",line,6))
        keep_reading = FALSE;

      /* If this is an ATOM or HETATM record, process it */
      else if (keep_reading == TRUE &&
               (!strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6)))
        {
          /* Get the residue name, number, chain-ID and atom type */
          strncpy(residue,line+17,3);
          residue[3] = '\0';
          strncpy(res_number,line+22,5);
          res_number[5] = '\0';
          chain = line[21];
          strncpy(atom_name,line+12,4);
          atom_name[4] = '\0';

          /* Check whether atom belongs to a new residue */
          if (chain != last_chain || strncmp(res_number,last_resnum,5) ||
              strncmp(residue,last_residue,3))
            {
              /* Residue is a new one, need to check whether it
                 is one of those used for the LIGPLOT and hence is to
                 be written out in full */
              got_residue = FALSE;
              want_residue = FALSE;

              /* Initialise search pointer to first residue in the
                 linked list */
              residue_ptr = first_residue_ptr;

              /* Loop through all the stored residues to check for a
                 match */
              while (residue_ptr != NULL && got_residue == FALSE)
                {
                  /* Check whether this is the correct residue */
                  if (!strncmp(res_number,residue_ptr->res_num,5) &&
                      !strncmp(residue,residue_ptr->res_name,3) &&
                      chain == residue_ptr->chain)
                    got_residue = TRUE;

                  /* If have the residue, check that it hasn't been
                     deleted */
                  if (got_residue == TRUE && residue_ptr->deleted == FALSE)
                    want_residue = TRUE;

                  /* If residue is wanted and the ligplot.res file is
                     required, then write out the residue details */
                  if (want_residue == TRUE && Write_Res_File == TRUE)
                    {
                      fprintf(ligplot_res,"%s%s%c",
                              residue_ptr->res_name,
                              residue_ptr->res_num,
                              residue_ptr->chain);

                      /* If this is a ligand residue, then write out
                         marker */
                      if (residue_ptr->inligand == TRUE)
                        fprintf(ligplot_res," *\n");
                      else
                        fprintf(ligplot_res,"  \n");
                    }

                  /* Get pointer to the next residue in the linked list */
                  residue_ptr = residue_ptr->next_residue_ptr;
                }

              /* Save the current residue number and chain */
              strncpy(last_residue,residue,3);
              last_residue[3] = '\0';
              strncpy(last_resnum,res_number,5);
              last_resnum[5] = '\0';
              last_chain = chain;
            }

          /* If this residue is wanted, then write the atom-record
             out to the .frm file */
          if (want_residue == TRUE)
            {
              /* Truncate the record at last non-blank character */
              ipos = 0;
              for (cpos = 0; cpos < 101 && line[cpos] != '\n'; cpos++)
                {
                  if (line[cpos] != ' ')
                    {
                      ipos = cpos + 1;
                    }
                }
              line[ipos] = '\0';
              
              /* Write the record out */
              fprintf(ligplot_frm,"%s\n",line);
              atom_counter++;
            }
        }

/* v.4.2--> */
      /* Get the next line from the PDB or XMAS file */
      if (Read_Xmas == FALSE)
        eofstatus = fgets(line,LINELEN,fil_pdb);
      else
/* v.4.2.2--> */
/*      eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type); */
        eofstatus = get_xmas2pdb_line(xmas,line,FALSE,rec_types,&type,&mol_id);
/* <--v.4.2.2 */
/* <--v.4.2 */
    }

  /* Print counts of records written out */
  printf("   Number of atoms written out       = %7d\n",atom_counter);

/* v.4.2--> */
  /* Close the PDB file */
  if (Read_Xmas == FALSE)
    fclose(fil_pdb);

  /* Close the output file */
  fclose(ligplot_frm);
/* <--v.4.2 */
}
/* v.4.0--> */
/***********************************************************************

write_lig_prot_bonds  -  Write out any covalent bonds between the ligand
                         and the protein to the .res file

***********************************************************************/

void write_lig_prot_bonds(void)
{
  int nbonds;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object1_ptr, *object2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;    

  /* Initialise variables */
  nbonds = 0;

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is a covalent, elastic bond */
      if (bond_ptr->bond_type == COVALENT && bond_ptr->elastic == TRUE)
        {
          /* Get pointers to the two atoms at either end of the
             bond, and their respective residues and object */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;
          residue1_ptr = atom1_ptr->residue_ptr;
          residue2_ptr = atom2_ptr->residue_ptr;
          object1_ptr = residue1_ptr->object_ptr;
          object2_ptr = residue2_ptr->object_ptr;

          /* Write out the bond to the .res file only if it is across
             two different objects, of which one is a ligand */
          if (object1_ptr != object2_ptr &&
              (object1_ptr->object_type == LIGAND ||
               object2_ptr->object_type == LIGAND))
            {
              /* Write out the bond */
              fprintf(ligplot_res,"#Covalent links: ");
              if (object1_ptr->object_type == LIGAND)
                fprintf(ligplot_res,"%s %s %s %c -> %s %s %s %c\n",
                        atom1_ptr->atom_name,
                        residue1_ptr->res_name,
                        residue1_ptr->res_num,
                        residue1_ptr->chain,
                        atom2_ptr->atom_name,
                        residue2_ptr->res_name,
                        residue2_ptr->res_num,
                        residue2_ptr->chain);
              else
                fprintf(ligplot_res,"%s %s %s %c -> %s %s %s %c\n",
                        atom2_ptr->atom_name,
                        residue2_ptr->res_name,
                        residue2_ptr->res_num,
                        residue2_ptr->chain,
                        atom1_ptr->atom_name,
                        residue1_ptr->res_name,
                        residue1_ptr->res_num,
                        residue1_ptr->chain);

              /* Increment count of bonds written out */
              nbonds++;
            }
        }

      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Print counts of records written out */
  printf("   Number of covalent ligand-protein bonds = %7d\n",nbonds);
}
/* <--v.4.0 */
/***********************************************************************

pare_down_hydrophobics  -  Mark any unwanted atoms in hydrophobic groups
                           for deletion

***********************************************************************/

void pare_down_hydrophobics(void)
{
  int iatom, iobject, iresid, natoms, nresid;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct object *object1_ptr, *object2_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  printf("\nMarking unwanted atoms ...\n");
  initialise_atoms();

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is a hydrophobic contact */
      if (bond_ptr->bond_type == CONTACT)
        {
          /* Mark both atoms as being involved in the bond */
          bond_ptr->first_atom_ptr->checked = TRUE;
          bond_ptr->second_atom_ptr->checked = TRUE;
        }

      /* Get the two objects these atoms belong to */
      object1_ptr = bond_ptr->first_atom_ptr->residue_ptr->object_ptr;
      object2_ptr = bond_ptr->second_atom_ptr->residue_ptr->object_ptr;

      /* If this is a bond between a hydrophobic group and anything
         other than a ligand residue, then delete (except for case of
         interface plots) */
      if (object1_ptr != object2_ptr)
        {
/* v.4.0--> */
          if (Interface_Plot == FALSE)
            {
/* <--v.4.0 */
              if ((object1_ptr->object_type == HYDROPHOBIC &&
                   object2_ptr->first_residue_ptr->inligand == FALSE) ||
                  (object2_ptr->object_type == HYDROPHOBIC &&
                   object1_ptr->first_residue_ptr->inligand == FALSE))
                {
                  /* Mark this bond for deletion */
                  bond_ptr->bond_type = DELETED;
                }
/* v.4.0--> */
            }
/* <--v.4.0 */
        }
      /* Get pointer for next object */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Process only if this is a hydrophobic group */
      if (object_ptr->object_type == HYDROPHOBIC)
        {
          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              
              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* If atom not involved in any hydrophobic contacts
                     then can mark it for deletion */
                  if (atom_ptr->checked == FALSE)
                    atom_ptr->deleted = TRUE;

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
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
}
/***********************************************************************

adjust_stored_coords  -  Adjust the stored coordinates to origin

***********************************************************************/

void adjust_stored_coords(int natoms,float coord_store[MAXATOMS][3],
                          float centre_x,float centre_y,float centre_z)
{
  int iatom;

  for (iatom = 0; iatom < natoms; iatom++)
    {
      coord_store[iatom][0] = coord_store[iatom][0] - centre_x;
      coord_store[iatom][1] = coord_store[iatom][1] - centre_y;
      coord_store[iatom][2] = coord_store[iatom][2] - centre_z;
    }
}
/***********************************************************************
  
calc_eigen_values  -  Subroutine to compute eigenvalues & eigenvectors
                      of a real symmetric matrix, from IBM SSP manual
                      (see p165).

                      Ian Tickle April 1992
 
                      (modified by David Moss February 1993)

Eigenvalues & vectors of real symmetric matrix stored in triangular form.
 
Arguments:
 
mv = 0 to compute eigenvalues only.
mv = 1 to compute eigenvalues & eigenvectors.
 
n - supplied dimension of n*n matrix.
 
a - supplied array of size n*(n+1)/2 containing n*n matrix in the order:
 
           1      2      3    ...
     1    a[0]
     2    a[1]   a[2]
     3    a[3]   a[4]   a[5]   ...
     ...

     NOTE a is used as working space and is overwritten on return.
     Eigenvalues are written into diagonal elements of a
     i.e.  a[0]  a[2]  a[5]  for a 3*3 matrix.
 
r - Resultant matrix of eigenvectors stored columnwise in the same
     order as eigenvalues.

***********************************************************************/

void calc_eigen_values(int mv, int n, double* a, double* r)
{
  int i, il, ilq, ilr, im, imq, imr, ind, iq, j, l, ll, lm, lq, m, mm, mq;
  double anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;

  if (mv)
    { for (i=1; i<n*n; i++) r[i]=0.;
      for (i=0; i<n*n; i+=n+1) r[i]=1.;
    }

  /* Initial and final norms (anorm & anrmx). */
  anorm=0.;
  iq=0;
  for (i=0; i<n; i++) for (j=0; j<=i; j++)
    { if (j!=i) anorm+=a[iq]*a[iq];
      ++iq;
    }

  if (anorm>0.)
    { anorm=sqrt(2.0*anorm);
      anrmx=1.0e-6*anorm/n;

      /* Compute threshold and initialise flag. */
      thr=anorm;
      do
        { thr/=n;
      do
        { ind=0;
          l=0;
          do
            { lq=l*(l+1)/2;
              ll=l+lq;
              m=l+1;
              ilq=n*l;
              do

                /* Compute sin & cos. */
                { mq=m*(m+1)/2;
                  lm=l+mq;
                  if (fabs(a[lm])>=thr)
                    { ind=1;
                      mm=m+mq;
                      x=0.5*(a[ll]-a[mm]);
                      y=-a[lm]/sqrt(a[lm]*a[lm]+x*x);
                      if (x<0.0) y=-y;
                      sinx=y/sqrt(2.0*(1.0+(sqrt(1.0-y*y))));
                      sinx2=sinx*sinx;
                      cosx=sqrt(1.0-sinx2);
                      cosx2=cosx*cosx;
                      sincs=sinx*cosx;

                      /* Rotate l & m columns. */
                      imq=n*m;
                      for (i=0; i<n; i++)
                        { iq=i*(i+1)/2;
                          if (i!=l && i!=m)
                            { if (i<m) im=i+mq;
                            else im=m+iq;
                              if (i<l) il=i+lq;
                              else il=l+iq;
                              x=a[il]*cosx-a[im]*sinx;
                              a[im]=a[il]*sinx+a[im]*cosx;
                              a[il]=x;
                            }
                          if (mv)
                            { ilr=ilq+i;
                              imr=imq+i;
                              x=r[ilr]*cosx-r[imr]*sinx;
                              r[imr]=r[ilr]*sinx+r[imr]*cosx;
                              r[ilr]=x;
                            }
                        }

                      x=2.*a[lm]*sincs;
                      y=a[ll]*cosx2+a[mm]*sinx2-x;
                      x=a[ll]*sinx2+a[mm]*cosx2+x;
                      a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2);
                      a[ll]=y;
                      a[mm]=x;
                    }

                  /* Tests for completion.
                     Test for m = last column. */
                } while (++m!=n);

              /* Test for l =penultimate column. */
            } while (++l!=n-1);
        } while (ind);
          
          /* Compare threshold with final norm. */
        } while (thr>anrmx);
    }
}

/***********************************************************************

eigen_sort  -  Sort eigenvalues & eigenvectors (if mv=1) in order of
               descending eigenvalue.
               Arguments are as for eigen, which must have been called.

***********************************************************************/

void eigen_sort(int mv, int n, double* a, double* r)
{
  int i, im, j, k, km, l, m;
  double am;

  k=0;
  for (i=0; i<n-1; i++)
    { im=i;
      km=k;
      am=a[k];
      l=0;
      for (j=0; j<n; j++)
        { if (j>i && a[l]>am)
            { im=j;
              km=l;
              am=a[l];
            }
          l+=j+2;
        }
      if (im!=i)
        { a[km]=a[k];
          a[k]=am;
          if (mv)
            { l=n*i;
              m=n*im;
              for (j=0; j<n; j++)
                { am=r[l];
                  r[l++]=r[m];
                  r[m++]=am;
                }
            }
        }
      k+=i+2;
    }
}
/***********************************************************************

vecprd  -  Calculates the vector product (vx,vy,vz) of two vectors 
           (px,py,pz) and (qx,qy,qz)

***********************************************************************/

void vecprd(float px, float py, float pz, float qx, float qy, float qz,
            float *vx, float *vy, float *vz)
{
  *vx = py * qz - pz * qy;
  *vy = pz * qx - px * qz;
  *vz = px * qy - py * qx;
}
/***********************************************************************

dot3  -  Returns the scalar product of the two vectors (x1,y1,z1) 
         and (x2,y2,z2)

***********************************************************************/

float dot3(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float dot_prod;

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

int orien3(float x1, float y1, float z1, float x2, float y2, float z2,
           float x3, float y3, float z3, float ex, float ey, float ez)
{
    float a, b, c, dx1, dy1, dz1, dx2, dy2, dz2, fe;
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

principal_components  -  Calculate the transformation to align the
                          supplied atom coords with their principal axes
                          along the x-, y-, and z-directions

***********************************************************************/

void principal_components(int natoms,float transformation[3][3],
                          float coord_store[MAXATOMS][3])
{
  int iatom, iorder;
  int nzeroy, nzeroz;
  int array_size, i, ipos, j, nsize;
  double amatrix[6], matrix[3][3], rmatrix[9];
  double x, y, z;
/*  float eigen_value[3]; */
  double xxsum, xysum, xzsum, yysum, yzsum, zzsum;

  /* Initialise counts */
  nzeroy = 0;
  nzeroz = 0;
  xxsum = 0.0;
  xysum = 0.0;
  xzsum = 0.0;
  yysum = 0.0;
  yzsum = 0.0;
  zzsum = 0.0;

  /* Initialise transformation matrix */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (i == j)
        transformation[i][j] = (float) (1.0);
      else
        transformation[i][j] = 0.0;

  /* If only have a single atom, then return */
  if (natoms < 2)
    return;

  /* Loop through all data points */
  for (iatom = 0; iatom < natoms; iatom++)
    {
      /* Retrieve the coordinates for this point */
      x = coord_store[iatom][0];
      y = coord_store[iatom][1];
      z = coord_store[iatom][2];

      /* Count instances where the y- and z-coords are zero */
      if (fabs(y) < 0.000001)
        nzeroy++;
      if (fabs(z) < 0.000001)
        nzeroz++;

      /* Compute the terms of the scalar products matrix */
      xxsum = xxsum+x*x;
      xysum = xysum+x*y;
      xzsum = xzsum+x*z;
      yysum = yysum+y*y;
      yzsum = yzsum+y*z;
      zzsum = zzsum+z*z;
    }
    
  /* Form the scalar products matrix */
  matrix[0][0] = xxsum;
  matrix[0][1] = xysum;
  matrix[0][2] = xzsum;
  matrix[1][0] = xysum;
  matrix[1][1] = yysum;
  matrix[1][2] = yzsum;
  matrix[2][0] = xzsum;
  matrix[2][1] = yzsum;
  matrix[2][2] = zzsum;

  /* If there are only two points, and both their z-coords are zero, or
     both their y- and z-coords are zero, adjust the order of the scalar
     products matrix whose eigenvalues are required */
  nsize = 3;
  if ((natoms == 2) && (nzeroz == natoms))
    {
      if (nzeroy == natoms)
        nsize = 1;
      else
        nsize = 2;
    }

  /* Calculate the eigen vectors and eigen values providing order of the
     matrix is not 1 */
  if (nsize > 1)
    {
      /* Convert matrix into form required by the eigen routine */
      ipos = 0;
      for (i = 0; i < nsize; i++)
        for (j = 0; j <= i; j++)
          {
            amatrix[ipos] = matrix[i][j];
            ipos++;
          }

      /* Calculate the eigen values */
      array_size = nsize * (nsize + 1) / 2;
      calc_eigen_values(1,nsize,amatrix,rmatrix);

      /* Sort the eigen vectors according to the eigen values */
      eigen_sort(1,nsize,amatrix,rmatrix);

      /* Retrieve the eigen values */
      i = j = 0;
      for (ipos = 0; ipos < array_size; ipos++)
        {
          /* If this is a diagonal term, then is one of the eigen values */
/*          if (i == j)
            eigen_value[i] = amatrix[ipos]; */

          /* Increment array indices */
          j++;
          if (j > i)
            {
              i++;
              j = 0;
            }
        }

      /* Retrieve the eigen vectors */
      i = j = 0;
      for (ipos = 0; ipos < nsize * nsize; ipos++)
        {
          /* Store the transformation matrix */
          transformation[j][i] = (float) rmatrix[ipos];

          /* Increment array indices */
          j++;
          if (j > nsize - 1)
            {
              i++;
              j = 0;
            }
        }

      /* If the matrix was of order 2, fix up row 3 and column 3 of the
         output transformation matrix */
      if (nsize == 2)
        {
          for (i = 0; i < 3; i++)
            {
              for (j = 0; j < 3; j++)
                {
                  if ((i == 2) && (j == 2))
                    transformation[i][j] = (float) (1.0);
                  else if ((i == 2) || (j == 2))
                    transformation[i][j] = 0.0;
                }
            }
        }

      /* If the matrix was of order 2, fix up row 3 and column 3 of the
         output transformation matrix */
      if (nsize == 2)
        {
          for (i = 0; i < 3; i++)
            {
              for (j = 0; j < 3; j++)
                {
                  if ((i == 3) && (j == 3))
                    transformation[i][j] = (float) (1.0);
                  else if ((i == 3) || (j == 3))
                    transformation[i][j] = 0.0;
                }
            }
        }
    }

    /* If the matrix was of order 1, make the transormation matrix the
       identity matrix */
    else
      {
        for (i = 0; i < 3; i++)
          {
/* v.4.2.1--> */
/*          for (j = 1; j < 4; j++) */
            for (j = 0; j < 3; j++)
/* <--v.4.2.1 */
              {
                if (i == j)
                  transformation[i][j] = (float) (1.0);
                else
                  transformation[i][j] = 0.0;
              }
          }
      }
    
  /* Check whether the sorted eigenvectors make a right-handed set */
  iorder = orien3(transformation[0][0],transformation[0][1],
                  transformation[0][2],0.0,0.0,0.0,
                  transformation[1][0],transformation[1][1],
                  transformation[1][2],transformation[2][0],
                  transformation[2][1],transformation[2][2]);
  
  /* If the sorted vectors don't make a right-handed set, then reverse the
     signs on the third vector */
  if (iorder == 1)
    {
      transformation[0][2] = - transformation[0][2];
      transformation[1][2] = - transformation[1][2];
      transformation[2][2] = - transformation[2][2];
    }
}
/***********************************************************************

apply_rotation  -  Apply the given rotation to the supplied coordinates

***********************************************************************/

void apply_rotation(float x,float y,float z,float *xnew,
                    float *ynew,float *znew,float mat[3][3])
{
  /* Apply the tranformation */
  *xnew = mat[0][0] * x + mat[1][0] * y + mat[2][0] * z;
  *ynew = mat[0][1] * x + mat[1][1] * y + mat[2][1] * z;
  *znew = mat[0][2] * x + mat[1][2] * y + mat[2][2] * z;
}
/***********************************************************************

calc_fit_coordinates  -  Align principal axes of original coords with
                         the x-, y- and z-axes. Repeat for each object

***********************************************************************/

void calc_fit_coordinates(int nligands)
{
  int iobject, iatom, iresid, natoms, nresid, nstored;
  int wanted;

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
  float new_x, new_y, new_z, x, y, z;
  float transformation[3][3];

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  printf("\nCalculating fits ...\n");
  nstored = 0;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms and store their original coords */
  while (atom_ptr != NULL)
    {
      wanted = FALSE;

      /* If there is no ligand, then use all atoms */
      if (nligands == 0)
        wanted = TRUE;

      /* Otherwise, only want this atom if it belongs to the ligand */
      else
        {
          object_ptr = atom_ptr->residue_ptr->object_ptr;
          if (object_ptr->object_type == LIGAND)
            wanted = TRUE;
        }

      /* Store this atom's coords */
      if (wanted == TRUE && nstored < MAXATOMS)
        {
          coord_store[nstored][0] = atom_ptr->original_x;
          coord_store[nstored][1] = atom_ptr->original_y;
          coord_store[nstored][2] = atom_ptr->original_z;
          nstored++;
        }

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }          

  /* If have some atoms to fit onto, then do the fit */
  if (nstored > 0)
    {
      /* Calculate the centre of mass of these atoms */
      centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

      /* Adjust coordinates to place centre of mass at the origin */
      adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

      /* Calculate transformation to orient the principal axes of the
         original structure along the x-, y- and z-axes */
      principal_components(nstored,transformation,coord_store);

      /* Get pointer to first atom again */
      atom_ptr = first_atom_ptr;

      /* Loop through all the atoms to compute their fitted coords */
      while (atom_ptr != NULL)
        {
          /* Get the atom's original coords */
          x = atom_ptr->original_x - mass_x;
          y = atom_ptr->original_y - mass_y;
          z = atom_ptr->original_z - mass_z;

          /* Apply rotation matrix */
          apply_rotation(x,y,z,&new_x,&new_y,&new_z,transformation);
          atom_ptr->fit_x = new_x;
          atom_ptr->fit_y = new_y;
          atom_ptr->fit_z = new_z;

          /* Get pointer to next atom in linked-list */
          atom_ptr = atom_ptr->next;
        }
    }

  /* Get the translation that puts each object's fitted coords at the
     origin */

  /* Get pointer to first object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all object's residues to calculate its centre of mass */
  while (object_ptr != NULL)
    {
      /* Initialise count of stored coordinates */
      nstored = 0;

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all the residue's atoms */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* Store this atom's fitted coords */
              if (nstored < MAXATOMS)
                {
                  coord_store[nstored][0] = atom_ptr->fit_x;
                  coord_store[nstored][1] = atom_ptr->fit_y;
                  coord_store[nstored][2] = 0.0;
                  nstored++;
                }

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Calculate the centre of mass of these atoms */
      centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

      /* Store this centre of mass */
      object_ptr->place_x = mass_x;
      object_ptr->place_y = mass_y;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

move_object  -  Move the given object by the supplied translation

***********************************************************************/

void move_object(struct object *object_ptr,float move_x,
                 float move_y,float move_z)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         transformation */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Apply the move */
          atom_ptr->x = atom_ptr->x - move_x;
          atom_ptr->y = atom_ptr->y - move_y;
          atom_ptr->z = atom_ptr->z - move_z;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Adjust the residue boundaries */
      residue_ptr->minx = residue_ptr->minx - move_x;
      residue_ptr->maxx = residue_ptr->maxx - move_x;
      residue_ptr->miny = residue_ptr->miny - move_y;
      residue_ptr->maxy = residue_ptr->maxy - move_y;

/* v.4.0--> */
      /* Adjust the residue's mean position */
      residue_ptr->flattened_mean_x = residue_ptr->flattened_mean_x - move_x;
      residue_ptr->flattened_mean_y = residue_ptr->flattened_mean_y - move_y;
/* <--v.4.0 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Adjust the object boundaries */
  object_ptr->minx = object_ptr->minx - move_x;
  object_ptr->maxx = object_ptr->maxx - move_x;
  object_ptr->miny = object_ptr->miny - move_y;
  object_ptr->maxy = object_ptr->maxy - move_y;
}
/***********************************************************************

place_group_at_origin  -  Place all the atoms of the current hydrophobic
                          group at the origin

***********************************************************************/

void place_group_at_origin(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Initialise object boundaries */
  object_ptr->minx = 0.0;
  object_ptr->maxx = 0.0;
  object_ptr->miny = 0.0;
  object_ptr->maxy = 0.0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, zeroing the coordinates */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Zero the coordinates */
          atom_ptr->x = 0.0;
          atom_ptr->y = 0.0;
          atom_ptr->z = 0.0;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Initialise residue boundaries */
      residue_ptr->minx = 0.0;
      residue_ptr->maxx = 0.0;
      residue_ptr->miny = 0.0;
      residue_ptr->maxy = 0.0;

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

rotate_object  -  Apply the given rotation to the given object

***********************************************************************/

void rotate_object(struct object *object_ptr,float v[3][3])
{
  int iatom, icoord, iresid, natoms, nresid;
  float x, y, z;
  float new_coords[3];

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         transformation */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Extract this atom's coordinates */
          x = atom_ptr->x;
          y = atom_ptr->y;
          z = atom_ptr->z;
        
          /* Apply the transformation */
          for (icoord = 0; icoord < 3; icoord++)
            {
              /* Calculate transformed coordinates */
              new_coords[icoord]
                =  ((x * v[0][icoord])
                    + (y * v[1][icoord])
                    + (z * v[2][icoord]));
            }

          /* Save the new coordinates */
          atom_ptr->x = new_coords[0];
          atom_ptr->y = new_coords[1];
          atom_ptr->z = new_coords[2];

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

squash_z  -  Set the z-coord of the ring-atoms (and attached end-bond
             atoms) to zero to completely flatten current ring

***********************************************************************/

void squash_z(struct bond *ring_bond_ptr)
{
  struct bond *bond_ptr;

  /* Get pointer to the first bond in the current ring */
  bond_ptr = ring_bond_ptr;

  /* Loop through all the bonds in the current ring */
  while (bond_ptr != NULL)
    {
      /* Set z-coordinate of both this bond's atoms to zero */
      bond_ptr->first_atom_ptr->z = 0.0;
      bond_ptr->second_atom_ptr->z = 0.0;

      /* Get the next bond in the current ring */
      bond_ptr = bond_ptr->next_link_ptr;
    }
}
/***********************************************************************

flatten_ring  -  Flatten the current ring, given the starting bond

***********************************************************************/

void flatten_ring(struct object *object_ptr,struct bond *ring_bond_ptr)
{
  int natoms, nbonds;

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
  float transformation[3][3];

  struct residue *residue_ptr;

  /* Extract all the coordinates of the current ring's atoms, so that
     can squash to make whole ring perfectly planar */
/* v.4.4--> */
/*  get_ring_atom_coords(ring_bond_ptr,coord_store,&natoms,&nbonds); */
/* v.4.5.2--> */
  // get_ring_atom_coords(ring_bond_ptr,coord_store,&natoms,&nbonds,FALSE);
  get_ring_atom_coords(ring_bond_ptr,coord_store,&natoms,&nbonds,-1,FALSE);
/* <--v.4.5.2 */
/* <--v.4.4 */

  /* Have the coords of all the atoms making up the current ring,
     so calculate the centre of mass of these atoms */
  centre_of_mass(coord_store,natoms,&mass_x,&mass_y,&mass_z);

  /* Transform the object to place the ring's centre of mass at the
     origin */
  move_object(object_ptr,mass_x,mass_y,mass_z);

  /* Apply the same transformation to the stored ring-atom coordinates
     prior to calculating the best-fit plane */
  adjust_stored_coords(natoms,coord_store,mass_x,mass_y,mass_z);

  /* Calculate transformation to orient the plane of the ring
     approximately in the x-y plane */
  principal_components(natoms,transformation,coord_store);

  /* Apply this transformation to the entire object */
  rotate_object(object_ptr,transformation);

  /* Zero the z-coordinates of the atoms in the ring and any end bonds
     springing off it */
  squash_z(ring_bond_ptr);

  /* Display message to show which ring has been flattened */
  residue_ptr = ring_bond_ptr->first_atom_ptr->residue_ptr;
  printf("   Flattened %2d-bond ring of:  %s %s %c\n",
         nbonds,residue_ptr->res_name,residue_ptr->res_num,
         residue_ptr->chain);
}
/***********************************************************************

flatten_all_rings  -  Flatten any rings in the current object

***********************************************************************/

int flatten_all_rings(struct object *object_ptr)
{
  int n_flattened;

  struct bond *bond_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise variables */
  n_flattened = 0;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a ring-bond which isn't part of a ring already
         flattened, then process */
      if (bond_ptr->ring_bond == TRUE && bond_ptr->flattened == FALSE)
        {
          /* Get all the other bonds involved in this ring */
          get_other_ring_bonds(bond_ptr);

          /* Flatten the ring */
          flatten_ring(object_ptr,bond_ptr);
          n_flattened++;
        }

      /* Get pointer to point to next bond for this object */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }

  /* Return the number of rings flattened */
  return(n_flattened);
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
            ang =       PI + ang;
        }
    }
  return ang;
}
/***********************************************************************

calculate_matrix  -  Calculate the 4 x 3 transformation matrix from the
                     supplied translation and angles to rotate about
                     the z-x-z axes

***********************************************************************/

/* v.4.2.2--> */
/* void calculate_matrix(float transl[3], float cos_theta[2], float sin_theta[2]
,
            float mat[4][3]) */
void calculate_matrix(float transl[3],float cos_theta[3],
                      float sin_theta[3],float mat[4][3])
/* <--v.4.2.2 */
{
  float xcys, xsyc;

  /* Calculate each of the terms in the transformation matrix */
  xcys = (transl[0] * cos_theta[0] + transl[1] * sin_theta[0]);
  xsyc = (transl[0] * sin_theta[0] - transl[1] * cos_theta[0]);
  mat[0][0] =   cos_theta[0] * cos_theta[2]
              - sin_theta[0] * cos_theta[1] * sin_theta[2];
  mat[1][0] =   sin_theta[0] * cos_theta[2]
              + cos_theta[0] * cos_theta[1] * sin_theta[2];
  mat[2][0] =   sin_theta[1] * sin_theta[2];
  mat[3][0] = - xcys * cos_theta[2]
              + xsyc * cos_theta[1] * sin_theta[2]
              - transl[2] * sin_theta[1] * sin_theta[2];
  mat[0][1] = -(cos_theta[0] * sin_theta[2]
              + sin_theta[0] * cos_theta[1] * cos_theta[2]);
  mat[1][1] = -(sin_theta[0] * sin_theta[2]
              - cos_theta[0] * cos_theta[1] * cos_theta[2]);
  mat[2][1] =   sin_theta[1] * cos_theta[2];
  mat[3][1] =   xcys * sin_theta[2]
              + xsyc * cos_theta[1] * cos_theta[2]
              - transl[2] * sin_theta[1] * cos_theta[2];
  mat[0][2] =   sin_theta[0] * sin_theta[1];
  mat[1][2] = - cos_theta[0] * sin_theta[1];
  mat[2][2] =   cos_theta[1];
  mat[3][2] = - xsyc * sin_theta[1] - transl[2] * cos_theta[1];
}              
/***********************************************************************

line_transformation  -  Calculate the full transformation that will
                        put the given 2 atoms along the x-axis, with the
                        first atom at the origin and the other atom
                        along the +ve (or -ve) x-axis

***********************************************************************/

void line_transformation(float coord_store[3][3], float matrix[4][3],
                         float splay_angle)
{
  int   icoord, itheta;
  float  calc, cos_theta[3],
        nx, ny, nz, sin_theta[3], theta[3],
        transl[3], xyz[3][3];

  /* Initialise variables */
  for (itheta = 0; itheta < 3; itheta++)
    {
      theta[itheta] = 0.0;
      sin_theta[itheta] = 0.0;
      cos_theta[itheta] = (float) (1.0);
    }

  /* Translation to origin is just the reverse of the coordinates of the
     first atom */
  transl[0] = coord_store[0][0];
  transl[1] = coord_store[0][1];
  transl[2] = coord_store[0][2];

  /* Translate the 2 atoms to the origin */
  for (icoord = 0; icoord < 3; icoord++)
    {
      xyz[0][icoord] = coord_store[0][icoord] - transl[icoord];
      xyz[1][icoord] = coord_store[1][icoord] - transl[icoord];
    }

  /* Calculate theta1 - the angle to rotate about the z-axis */
  nx = xyz[1][0];
  ny = xyz[1][1];
  nz = xyz[1][2];
  theta[0] = get_angle(ny,-nx);

  /* Calculate theta2 - the angle to rotate about the x-axis to put the
     line along the +ve y-axis */
  calc = nx * nx + ny * ny;
  theta[1] = - get_angle(nz,(float) sqrt(calc));
  theta[1] = (float) (theta[1] + PI / 2.0);

  /* Third rotation is thro' 90 degrees about z-axis to put from +ve
     y-axis to +ve x-axis */
  theta[2] = (float) (PI / 2.0);

  /* If want the line pointing at some splay-angle round from the +ve
     x-axis, then add that to the final rotation */
  theta[2] = theta[2] - splay_angle;

  /* Precalculate sine and cosine values for theta1 and theta2 */
  cos_theta[0] = (float) cos(theta[0]);
  cos_theta[1] = (float) cos(theta[1]);
  sin_theta[0] = (float) sin(theta[0]);
  sin_theta[1] = (float) sin(theta[1]);
  if (fabs(cos_theta[0]) < 0.000001)
    cos_theta[0] = 0.0;
  if (fabs(cos_theta[1]) < 0.000001)
    cos_theta[1] = 0.0;
  if (fabs(sin_theta[0]) < 0.000001)
    sin_theta[0] = 0.0;
  if (fabs(sin_theta[1]) < 0.000001)
    sin_theta[1] = 0.0;

  /* Precalculate sine and cosine values for theta3 */
  cos_theta[2] = (float) cos(theta[2]);
  sin_theta[2] = (float) sin(theta[2]);
  if (fabs(cos_theta[2]) < 0.000001)
    cos_theta[2] = 0.0;
  if (fabs(sin_theta[2]) < 0.000001)
    sin_theta[2] = 0.0;

  /* Calculate the transformation matrix */
  calculate_matrix(transl,cos_theta,sin_theta,matrix);
}
/***********************************************************************

apply_transformation  -  Apply the tranformation to the supplied
                         coordinates

***********************************************************************/

void apply_transformation(float x,float y,float z,float *xnew,
                          float *ynew,float *znew,float mat[4][3])
{
  /* Apply the tranformation */
  *xnew = mat[0][0] * x + mat[1][0] * y + mat[2][0] * z + mat[3][0];
  *ynew = mat[0][1] * x + mat[1][1] * y + mat[2][1] * z + mat[3][1];
  *znew = mat[0][2] * x + mat[1][2] * y + mat[2][2] * z + mat[3][2];
}
/***********************************************************************

apply_xy_transformation  -  Apply the tranformation to the supplied
                            x and y coordinates

***********************************************************************/

void apply_xy_transformation(float x,float y,float *xnew,float *ynew,
                          float mat[4][3])
{
  /* Apply the tranformation */
  *xnew = mat[0][0] * x + mat[1][0] * y + mat[3][0];
  *ynew = mat[0][1] * x + mat[1][1] * y + mat[3][1];
}
/***********************************************************************

transform_object  -  Apply the given transformation to the given object

***********************************************************************/

void transform_object(struct object *object_ptr,float matrix[4][3])
{
  int iatom, iresid, natoms, nresid;
  float x, y, z;
  float new_x, new_y, new_z;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         transformation */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Extract this atom's coordinates */
          x = atom_ptr->x;
          y = atom_ptr->y;
          z = atom_ptr->z;
        
          /* Apply the transformation */
          apply_transformation(x,y,z,&new_x,&new_y,&new_z,matrix);

          /* Save the new coordinates */
          atom_ptr->x = new_x;
          atom_ptr->y = new_y;
          atom_ptr->z = new_z;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

stretch_bond  -  Stretch the current bond, moving all atoms off one of
                 its ends by the required stretch distance

***********************************************************************/

void stretch_bond(struct bond *stretch_bond_ptr,float move_dist)
{
  int loop, nbonds, end;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Initialise all the atom flags */
  initialise_atoms();

  /* If first atom has fewer bonds downstream of it, then want to
     process the bonds coming off it */
  if (stretch_bond_ptr->nbonds_from_first_atom
      < stretch_bond_ptr->nbonds_from_second_atom)
    {
      end = FIRST;
      move_dist = -1 * move_dist;
    }

  /* Otherwise, use the other end of the bond */
  else
    end = SECOND;

  /* Set pointer to the first of the bonds sprouting off the stretch
     bond */
  bond_link_ptr = stretch_bond_ptr->first_bond_link_ptr;

  /* Set the stretch-bond as already checked */
  stretch_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the stretch-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
         atom of the test-bond */
      if (bond_link_ptr->bond_end == end)
        {
          /* Get the bond this link is pointing to and add to stack */
          bond_ptr = bond_link_ptr->bond_ptr;

          /* Add to stack provided it is not an elastic bond */
          if (bond_ptr->elastic == FALSE)
            {
              /* If this is the first to be added to the stack, set both
                 stack-pointers to point to it */
              if (next_free_stack_ptr == NULL)
                {
                  next_free_stack_ptr = bond_ptr;
                  next_stack_ptr = bond_ptr;
                }

              /* Otherwise, make last bond on stack point to this one
                 and set this one as nopw the last on the stack */
              else
                {
                  next_free_stack_ptr->next_stack_ptr = bond_ptr;
                  next_free_stack_ptr = bond_ptr;
                }
            }

          /* Mark bond as added to the stack */
          bond_ptr->checked = TRUE;
        }
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Loop over the bond's two atoms */
      for (loop = 0; loop < 2; loop++)
        {
          /* Get pointer to appropriate atom */
          if (loop == 0)
            atom_ptr = bond_ptr->first_atom_ptr;
          else
            atom_ptr = bond_ptr->second_atom_ptr;

          /* If first atom has not been moved, then move it */
          if (atom_ptr->checked == FALSE)
            {
              atom_ptr->x = atom_ptr->x + move_dist;
              atom_ptr->checked = TRUE;
            }
        }

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
         current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
        {
          /* Get the pointer to this bond */
          other_bond_ptr = bond_link_ptr->bond_ptr;

          /* If this bond hasn't already been processed, then move
             its atoms (if they need moving) and add the bond to the
             stack */
          if (other_bond_ptr->checked == FALSE &&
              other_bond_ptr->elastic == FALSE)
            {
              next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
              next_free_stack_ptr = other_bond_ptr;

              /* If the stack has run out, restart it at the bond
                 just entered */
              if (next_stack_ptr == NULL)
                next_stack_ptr = other_bond_ptr;

              /* Mark bond as added to the stack */
              other_bond_ptr->checked = TRUE;
            }

          /* Go to the next bond-link in the linked list */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }
    }
    
  /* Update the counter of bonds downstream on stretched bond */
  if (end == FIRST)
    stretch_bond_ptr->nbonds_from_first_atom = nbonds;
  else
    stretch_bond_ptr->nbonds_from_second_atom = nbonds;
}
/***********************************************************************

correct_bond_lengths  -  Correct bond-lengths in case any have become
                         deformed during the flattening of the ring groups

***********************************************************************/

void correct_bond_lengths(struct object *object_ptr)
{
  float bond_length, coord_store[3][3], diff, matrix[4][3], move_dist;
  float x1, x2, y1, y2, z1, z2;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object_bond *object_bond_ptr;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a rotatable bond, or an end bond, then process */
      if (bond_ptr->rotatable_bond == TRUE || bond_ptr->end_bond == TRUE)
        {
          /* Get the pointers to the two atoms at either end of the bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get the atom coordinates */
          x1 = atom1_ptr->x;
          y1 = atom1_ptr->y;
          z1 = atom1_ptr->z;
          x2 = atom2_ptr->x;
          y2 = atom2_ptr->y;
          z2 = atom2_ptr->z;

          /* Calculate length of bond */
          bond_length = (float) sqrt((x2 - x1) * (x2 - x1)
                             + (y2 - y1) * (y2 - y1)
                             + (z2 - z1) * (z2 - z1));

          /* If this bond is significantly different from its original
             value, need to correct it */
          diff = (float) (fabs(bond_length - bond_ptr->bond_length));
          if (diff > 0.01)
            {
              /* Get the distance by which the bond needs to be stretched
                 or shrunk */
              move_dist = bond_ptr->bond_length - bond_length;

              /* Get the coordinates of the two atoms defining this bond */
              coord_store[0][0] = x1;
              coord_store[0][1] = y1;
              coord_store[0][2] = z1;
              coord_store[1][0] = x2;
              coord_store[1][1] = y2;
              coord_store[1][2] = z2;
              
              /* Calculate the transformation that will put this bond
                 along the x-axis, with atom 1 at the origin and atom 2
                 pointing in the +ve x-direction */
              line_transformation(coord_store,matrix,0.0);

              /* Apply the transformation to all atoms in the current
                 object */
              transform_object(object_ptr,matrix);

              /* If this is an end-bond, then only need to move a single
                 atom */
              if (bond_ptr->end_bond == TRUE)
                {
                  /* If the bond's first atom is at the free end */
                  if (bond_ptr->nfirst_atom_links == 0)
                    {
                      /* Move first atom (which should currently be
                         at the origin) back in the -ve x-direction */
                      atom1_ptr->x = atom1_ptr->x - move_dist;
                    }

                  /* Otherwise, move the other atom */
                  else
                    {
                      /* Move second atom (which should currently be
                         on the +ve x-axis) in the +ve x-direction */
                      atom2_ptr->x = atom2_ptr->x + move_dist;
                    }
                }
              /* Otherwise, if this is a rotatable bond, need to
                 move all atoms reachable from one of its ends
                 (the one with fewest connections) */
              else
                stretch_bond(bond_ptr,move_dist);
            }
        }
      
      /* Get pointer to point to next bond for this object */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }
}
/***********************************************************************

plane_transformation - Calculate the full transformation that will
                       put the given 3 atoms in the x-y plane, with
                       the middle atom at the origin and the other
                       two atoms placed either side of the positive
                       x-axis s.t. the angle 1-2-3 is bisected by
                       that axis

***********************************************************************/

void plane_transformation(float coord[3][3], float matrix[4][3],
                          int in_neg_direction)
{
  int   icoord, itheta;
  float bisect[2], calc, cos_theta[3], mag1, mag2, normal[3],
        newx, newy, newz, nx, ny, nz, sin_theta[3], theta[3],
        transl[3], xyz[3][3];

  /* Initialise variables */
  for (itheta = 0; itheta < 3; itheta++)
    {
      theta[itheta] = 0.0;
      sin_theta[itheta] = 0.0;
      cos_theta[itheta] = (float) (1.0);
    }

  /* Translation to origin is just the reverse of the coordinates of the
     middle atom */
  transl[0] = coord[1][0];
  transl[1] = coord[1][1];
  transl[2] = coord[1][2];

  /* Translate the 3 atoms to the origin */
  for (icoord = 0; icoord < 3; icoord++)
    {
      xyz[0][icoord] = coord[0][icoord] - transl[icoord];
      xyz[1][icoord] = coord[1][icoord] - transl[icoord];
      xyz[2][icoord] = coord[2][icoord] - transl[icoord];
    }

  /* Calculate normal to plane holding the 3 atoms */
  normal[0] = xyz[2][1] * xyz[0][2] - xyz[2][2] * xyz[0][1];
  normal[1] = xyz[2][2] * xyz[0][0] - xyz[2][0] * xyz[0][2];
  normal[2] = xyz[2][0] * xyz[0][1] - xyz[2][1] * xyz[0][0];
  nx = normal[0];
  ny = normal[1];
  nz = normal[2];

  /* Calculate theta1 - the angle to rotate about the z-axis */
  theta[0] = get_angle(ny,-nx);

  /* Calculate theta2 - the angle to rotate about the x-axis */
  calc = nx * nx + ny * ny;
  theta[1] = - get_angle(nz,(float) sqrt(calc));

  /* Precalculate sine and cosine values for theta1 and theta2 */
  cos_theta[0] = (float) cos(theta[0]);
  cos_theta[1] = (float) cos(theta[1]);
  sin_theta[0] = (float) sin(theta[0]);
  sin_theta[1] = (float) sin(theta[1]);
  if (fabs(cos_theta[0]) < 0.000001)
    cos_theta[0] = 0.0;
  if (fabs(cos_theta[1]) < 0.000001)
    cos_theta[1] = 0.0;
  if (fabs(sin_theta[0]) < 0.000001)
    sin_theta[0] = 0.0;
  if (fabs(sin_theta[1]) < 0.000001)
    sin_theta[1] = 0.0;

  /* Calculate the transformation matrix */
  calculate_matrix(transl,cos_theta,sin_theta,matrix);

  /* Transform coordinates of atoms 1 and 3 */
  apply_transformation(coord[0][0],coord[0][1],coord[0][2],
                       &newx,&newy,&newz,matrix);
  xyz[0][0] = newx;
  xyz[0][1] = newy;
  xyz[0][2] = newz;
  apply_transformation(coord[2][0],coord[2][1],coord[2][2],
                       &newx,&newy,&newz,matrix);
  xyz[2][0] = newx;
  xyz[2][1] = newy;
  xyz[2][2] = newz;

  /* Calculate theta3 - the get_angle to rotate about the z-axis s.t.
     the get_angle between the fragment's atoms 1-2-3 is bisected by
     the +ve x-axis */

  /* Find vector in direction of the bisector */
  mag1 = xyz[0][0] * xyz[0][0] + xyz[0][1] * xyz[0][1];
  mag2 = xyz[2][0] * xyz[2][0] + xyz[2][1] * xyz[2][1];
  if (mag1 > 0.000001 && mag2 > 0.000001)
    {
      bisect[0] = (float) (xyz[0][0] / sqrt(mag1) + xyz[2][0] / sqrt(mag2));
      bisect[1] = (float) (xyz[0][1] / sqrt(mag1) + xyz[2][1] / sqrt(mag2));
    }
  else
    printf("+WARNING: Bond-length zero!!!\n");

  /* Calculate get_angle the bisector makes with the x-axis */
  theta[2] = get_angle(bisect[0],bisect[1]);

  /* If want the bisector pointing along the -ve, rather than the +ve,
     x-axis, then turn through a further 180 degrees */
  if (in_neg_direction == TRUE)
    theta[2] = theta[2] - PI;

  /* Precalculate sine and cosine values for theta3 */
  cos_theta[2] = (float) cos(theta[2]);
  sin_theta[2] = (float) sin(theta[2]);
  if (fabs(cos_theta[2]) < 0.000001)
    cos_theta[2] = 0.0;
  if (fabs(sin_theta[2]) < 0.000001)
    sin_theta[2] = 0.0;

  /* Calculate the transformation matrix */
  calculate_matrix(transl,cos_theta,sin_theta,matrix);
}
/***********************************************************************

transform_downstream_atoms  -  Apply given transformation to all the
                               atoms downstream of the given bond

***********************************************************************/

void transform_downstream_atoms(struct bond *reference_bond_ptr,int end,
                                float matrix[4][3])
{
  int loop, nbonds;

  float newx, newy, newz;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Initialise all the atom flags */
  initialise_atoms();

  /* Mark both this bond's atoms as already checked */
  reference_bond_ptr->first_atom_ptr->checked = TRUE;
  reference_bond_ptr->second_atom_ptr->checked = TRUE;

  /* Set pointer to the first of the bonds sprouting off the reference
     bond */
  bond_link_ptr = reference_bond_ptr->first_bond_link_ptr;

  /* Set the reference-bond as already checked */
  reference_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the reference-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
         end */
      if (bond_link_ptr->bond_end == end)
        {
          /* Get the bond this link is pointing to and add to stack */
          bond_ptr = bond_link_ptr->bond_ptr;

          /* Add to stack provided it is not an elastic bond */
          if (bond_ptr->elastic == FALSE)
            {
              /* If this is the first to be added to the stack, set both
                 stack-pointers to point to it */
              if (next_free_stack_ptr == NULL)
                {
                  next_free_stack_ptr = bond_ptr;
                  next_stack_ptr = bond_ptr;
                }

              /* Otherwise, make last bond on stack point to this one
                 and set this one as nopw the last on the stack */
              else
                {
                  next_free_stack_ptr->next_stack_ptr = bond_ptr;
                  next_free_stack_ptr = bond_ptr;
                }
            }

          /* Mark bond as added to the stack */
          bond_ptr->checked = TRUE;
        }
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Loop over the bond's two atoms */
      for (loop = 0; loop < 2; loop++)
        {
          /* Get pointer to appropriate atom */
          if (loop == 0)
            atom_ptr = bond_ptr->first_atom_ptr;
          else
            atom_ptr = bond_ptr->second_atom_ptr;

          /* If first atom has not been transformed, then apply the
             matrix transformation */
          if (atom_ptr->checked == FALSE)
            {
              /* Apply the transformation to the current atom */
              apply_transformation(atom_ptr->x,atom_ptr->y,atom_ptr->z,
                                   &newx,&newy,&newz,matrix);
              atom_ptr->x = newx;
              atom_ptr->y = newy;
              atom_ptr->z = newz;
              atom_ptr->checked = TRUE;
            }
        }

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
         current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
        {
          /* Get the pointer to this bond */
          other_bond_ptr = bond_link_ptr->bond_ptr;

          /* If this bond hasn't already been processed, then move
             its atoms (if they need moving) and add the bond to the
             stack */
          if (other_bond_ptr->checked == FALSE &&
              other_bond_ptr->elastic == FALSE)
            {
              next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
              next_free_stack_ptr = other_bond_ptr;

              /* If the stack has run out, restart it at the bond
                 just entered */
              if (next_stack_ptr == NULL)
                next_stack_ptr = other_bond_ptr;

              /* Mark bond as added to the stack */
              other_bond_ptr->checked = TRUE;
            }

          /* Go to the next bond-link in the linked list */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }
    }
}
/***********************************************************************

offshoot_splay  -  Check bonds springing off each ring to set them at
                   appropriate splay angles away from the ring

***********************************************************************/

void offshoot_splay(struct object *object_ptr,struct bond *ring_bond_ptr,
                    struct bond *other_ring_bond_ptr,
                    struct bond *connect_bond_ptr[CONNECTIONS],
                    int nstore,int end)
{
  int istore;

  float coord_store[3][3], matrix[4][3], splay_angle;
  float newx, newy, newz;
  float splay_start, splay_diff;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr, *common_atom_ptr;
  struct bond *bond_ptr;

/*    int counter, istore, nzeroz,
    stored_coords, ibond;
    int other_atom[CONNECTIONS + 1], other_bond[CONNECTIONS + 1],
    other_sign[CONNECTIONS + 1];
    struct link_bond *other_bond_ptr; */

  /* Get pointers to the two atoms of the ring-bond of interest */
  if (end == FIRST)
    {
      atom1_ptr = ring_bond_ptr->first_atom_ptr;
      atom2_ptr = ring_bond_ptr->second_atom_ptr;
    }
  else
    {
      atom1_ptr = ring_bond_ptr->second_atom_ptr;
      atom2_ptr = ring_bond_ptr->first_atom_ptr;
    }

  /* Store the pointer to the atom that is common to both ring bonds
     and hence is the atom about which all the transformations will be
     performed */
  common_atom_ptr = atom1_ptr;

  /* Store the atom coordinates, with the common atom between the
     two ring bonds as the second (middle) atom */
  coord_store[0][0] = atom2_ptr->x;
  coord_store[0][1] = atom2_ptr->y;
  coord_store[0][2] = atom2_ptr->z;
  coord_store[1][0] = atom1_ptr->x;
  coord_store[1][1] = atom1_ptr->y;
  coord_store[1][2] = atom1_ptr->z;

  /* Get the pointer to whichever of the other ring-bond's atoms is
     the third one of interest */
  if (other_ring_bond_ptr->first_atom_ptr == common_atom_ptr)
    atom_ptr = other_ring_bond_ptr->second_atom_ptr;
  else
    atom_ptr = other_ring_bond_ptr->first_atom_ptr;

  /* Store the coordinates of this third atom */
  coord_store[2][0] = atom_ptr->x;
  coord_store[2][1] = atom_ptr->y;
  coord_store[2][2] = atom_ptr->z;

  /* Calculate transformation that will place atoms 1 and 3 in the
     x-y plane, with their bisector pointing in the -ve x-direction */
  plane_transformation(coord_store,matrix,TRUE);

  /* Apply the transformation to all atoms in the object */
  transform_object(object_ptr,matrix);

  /* Calculate the splay-angle for the non-ring bonds */
  if (nstore == 1)
    {
      splay_start = 0.0;
      splay_diff = 0.0;
    }
  else if (nstore > 1)
    {
      splay_start = - PI / (nstore + 1);
      splay_diff = (float) (- 2.0 * splay_start);
    }

  /* Initialise values for splaying out the non-ring bonds */
  coord_store[0][0] = 0.0;
  coord_store[0][1] = 0.0;
  coord_store[0][2] = 0.0;

  /* Loop over all the stored non-ring bonds that are to be splayed */
  for (istore = 0; istore < nstore; istore++)
    {
      /* Get pointer to this bond */
      bond_ptr = connect_bond_ptr[istore];

      /* Get pointers to the bond's two atoms */
      atom1_ptr = bond_ptr->first_atom_ptr;
      atom2_ptr = bond_ptr->second_atom_ptr;

      /* Determine which is at the other end of the bond */
      if (bond_ptr->first_atom_ptr == common_atom_ptr)
        {
          atom_ptr = bond_ptr->second_atom_ptr;
          end = SECOND;
        }
      else
        {
          atom_ptr = bond_ptr->first_atom_ptr;
          end = FIRST;
        }

      /* Store ths atom's coords */
      coord_store[1][0] = atom_ptr->x;
      coord_store[1][1] = atom_ptr->y;
      coord_store[1][2] = atom_ptr->z;

      /* Calculate the transformation that will place the current
         bond along the +ve x-axis */
      splay_angle = splay_start + istore * splay_diff;

      line_transformation(coord_store,matrix,splay_angle);

      /* Apply the transformation to the current atom */
      apply_transformation(atom_ptr->x,atom_ptr->y,atom_ptr->z,
                           &newx,&newy,&newz,matrix);
      atom_ptr->x = newx;
      atom_ptr->y = newy;
      atom_ptr->z = newz;

      /* If the current bond is not an end-bond, then need to apply
         the transformation to all the atoms downstream of it */
      if (bond_ptr->end_bond == FALSE)
        {
          /* Transform all atoms downstream of the current bond */
          transform_downstream_atoms(bond_ptr,end,matrix);
        }
    }
}
/***********************************************************************

flatten_ring_offshoots  -  Check bonds springing off each ring in the
                           current object to set them at appropriate
                           splay angles away from the ring

***********************************************************************/

void flatten_ring_offshoots(struct object *object_ptr)
{
  int end, loop;
  int non_ring_bonds, nring_bonds, nstore;

  struct bond *bond_ptr, *other_bond_ptr, *other_ring_bond_ptr,
  *connect_bond_ptr[CONNECTIONS];
  struct bond_link *bond_link_ptr;
  struct object_bond *object_bond_ptr;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a ring bond, then process */
      if (bond_ptr->ring_bond == TRUE)
        {
          /* Loop over the bond's two atoms */
          for (loop = 0; loop < 2; loop++)
            {
              /* Get pointer to appropriate atom */
              if (loop == 0)
                end = FIRST;
              else
                end = SECOND;

              /* See how many bonds, and of what types, are attached to
                 this atom */

              /* Initalise counts and pointers */
              nstore = 0;
              non_ring_bonds = 0;
              nring_bonds = 0;
              other_ring_bond_ptr = NULL;
              
              /* Get the pointer to the first of the bonds coming off the 
                 current bond */
              bond_link_ptr = bond_ptr->first_bond_link_ptr;
              
              /* Check the bond-types of all the bonds coming off the first
                 atom*/
              while (bond_link_ptr != NULL)
                {
                  /* If it comes off the first atom the current bond, then
                     check its type */
                  if (bond_link_ptr->bond_end == end)
                    {
                      /* Get pointer to the other bond and store */
                      other_bond_ptr = bond_link_ptr->bond_ptr;
                      
                      /* Increment counts of ring- and non-ring bonds, and
                         store pointers to these bonds */
                      if (other_bond_ptr->ring_bond == TRUE)
                        {
                          other_ring_bond_ptr = other_bond_ptr;
                          nring_bonds++;
                        }
                      else
                        {
                          non_ring_bonds++;
                          connect_bond_ptr[nstore] = other_bond_ptr;
                          nstore++;
                        }
                    }
              
                  /* Go to the next bond-link in the linked list */
                  bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
                }

              /* If have both ring-bonds and non-ring bonds at this
                 junction, need to splay the non-ring bonds */
              if (nring_bonds == 1 && non_ring_bonds > 0)
                {
                  /* Perform the splay */
                  offshoot_splay(object_ptr,bond_ptr,other_ring_bond_ptr,
                                 connect_bond_ptr,nstore,end);
                }
              else if (nring_bonds > 1 && non_ring_bonds > 0)
                {
                  printf("*** WARNING. %d ring bonds and %d ",
                         nring_bonds + 1,non_ring_bonds);
                  printf(" non-ring bonds. May cause problems!\n");
                  Nwarnings++;
                }
            }
        }

      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }           
}
/***********************************************************************

flatten_bonds  -  Flatten out the bonds on the given side of the current
                  bond so that they lie in the x-y plane

***********************************************************************/

void flatten_bonds(struct bond *current_bond_ptr,int end)
{
  int have_ring, have_nonring, ibond, nstored, other_end;

  float coord_store[3][3], stored_coords[CONNECTIONS][3], matrix[4][3];
  float max_z, splay_angle;
  float new_x, new_y, new_z, x, y, z;

  struct bond *bond_ptr, *stored_bond_ptr[CONNECTIONS];
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr, *common_atom_ptr;

  /* Initialise variables */
  have_ring = FALSE;
  have_nonring = FALSE;
  nstored = 0;
  max_z = 0.0;

  /* Get pointer to the pivot atom (currently at the origin) */
  if (end == FIRST)
    common_atom_ptr = current_bond_ptr->first_atom_ptr;
  else
    common_atom_ptr = current_bond_ptr->second_atom_ptr;

  /* Get the pointer to the first of the bonds coming off the 
     test_bond */
  bond_link_ptr = current_bond_ptr->first_bond_link_ptr;

  /* Loop over all the bonds attached to the current one at this end */
  while (bond_link_ptr != NULL)
    {
      /* Check this bond comes off the appropriate end */
      if (bond_link_ptr->bond_end == end)
        {
          /* Get pointer to this bond */
          bond_ptr = bond_link_ptr->bond_ptr;

          /* Check if this is a ring bond. If it is, then don't
             want to apply the splaying */
          if (bond_ptr->ring_bond == TRUE)
            have_ring = TRUE;
          else
            have_nonring = TRUE;

          /* Get the pointer to the atom at the other end of the bond */
          if (bond_ptr->first_atom_ptr == common_atom_ptr)
            atom_ptr = bond_ptr->second_atom_ptr;
          else
            atom_ptr = bond_ptr->first_atom_ptr;

          /* Store ths atom's coords */
          if (nstored < CONNECTIONS)
            {
              stored_bond_ptr[nstored] = bond_ptr;
              stored_coords[nstored][0] = atom_ptr->x;
              stored_coords[nstored][1] = atom_ptr->y;
              stored_coords[nstored][2] = atom_ptr->z;
              nstored++;
            }

          /* Get maximum z-deviation */
          z = atom_ptr->z;
          if (fabs(z) > max_z)
            max_z = (float) (fabs(z));
        }

      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* If the bonds require flattening, then do it */
  if (max_z > 0.001)
    {
      /* If have both a ring, and other bonds springing off this end,
         then flatten as for one bond */
      if (have_ring == TRUE && have_nonring == TRUE)
        nstored = 1;

      /* If have just a ring, then flatten as for one bond */
      if (have_ring == TRUE)
        nstored = 1;
    
      /* If there is just the one bond coming off this one, swivel
         it, and all bonds beyond it, into the x-y plane */
      if (nstored == 1)
        {
          /* Create two dummy atoms, one at the origin, and the other
             having y- and z-coords opposite in sign to the first atom
             stored */
          coord_store[0][0] = stored_coords[0][0];
          coord_store[0][1] = stored_coords[0][1];
          coord_store[0][2] = stored_coords[0][2];
          coord_store[1][0] = 0.0;
          coord_store[1][1] = 0.0;
          coord_store[1][2] = 0.0;
          coord_store[2][0] = stored_coords[0][0];
          coord_store[2][1] = - stored_coords[0][1];
          coord_store[2][2] = - stored_coords[0][2];

          /* Calculate transformation that will this place bond (and its
             mirror image in the x-y plane, with their bisector pointing
             in the +ve x-direction */
          plane_transformation(coord_store,matrix,FALSE);
              
          /* Apply the transformation to all atoms springing off
             the current bond */
          transform_downstream_atoms(current_bond_ptr,end,matrix);
        }

  /* If there are two bonds coming off this one, then calculate the
     transformation that will swivel the plane defined by the two bonds
     into the x-y plane */
      else if (nstored == 2)
        {
          /* Move the first of the atoms to position 1, and replace
             position two with a dummy atom at the origin */
          coord_store[0][0] = stored_coords[0][0];
          coord_store[0][1] = stored_coords[0][1];
          coord_store[0][2] = stored_coords[0][2];
          coord_store[1][0] = 0.0;
          coord_store[1][1] = 0.0;
          coord_store[1][2] = 0.0;
          coord_store[2][0] = stored_coords[1][0];
          coord_store[2][1] = stored_coords[1][1];
          coord_store[2][2] = stored_coords[1][2];

          /* Calculate transformation that will this place atoms 1 and 3
             in the x-y plane, with their bisector pointing in the +ve
             x-direction */
          plane_transformation(coord_store,matrix,FALSE);
              
          /* Apply the transformation to all atoms springing off
             the current bond */
          transform_downstream_atoms(current_bond_ptr,end,matrix);
        }

      /* If there are more than two bonds coming off this one, then need
         to treat each one separately, splaying them out in the x-y
         plane at regularly-spaced angular intervals */
      else if (nstored > 2)
        {
          /* Loop through each of the bonds in turn */
          for (ibond = 0; ibond < nstored; ibond++)
            {
              /* Get the coordinates of this bond */
              coord_store[0][0] = 0.0;
              coord_store[0][1] = 0.0;
              coord_store[0][2] = 0.0;
              coord_store[1][0] = stored_coords[ibond][0];
              coord_store[1][1] = stored_coords[ibond][1];
              coord_store[1][2] = stored_coords[ibond][2];

              /* Calculate the splay-angle round from x-axis that want
                 this bond to be placed at in the x-y plane */
              splay_angle
                = (float) (PI - (ibond + 1) * (2.0 * PI / (nstored + 1)));
              
              /* Calculate the transformation that will place the current
                 bond at the correct splay angle */
              line_transformation(coord_store,matrix,splay_angle);

              /* Get the bond to apply the transformation to */
              bond_ptr = stored_bond_ptr[ibond];

              /* Get the pointer to the atom at the other end of the bond */
              if (bond_ptr->first_atom_ptr == common_atom_ptr)
                {
                  atom_ptr = bond_ptr->second_atom_ptr;
                  other_end = SECOND;
                }
              else
                {
                  atom_ptr = bond_ptr->first_atom_ptr;
                  other_end = FIRST;
                }

              /* Extract this atom's coordinates */
              x = atom_ptr->x;
              y = atom_ptr->y;
              z = atom_ptr->z;
        
              /* Apply the transformation */
              apply_transformation(x,y,z,&new_x,&new_y,&new_z,matrix);

              /* Save the new coordinates */
              atom_ptr->x = new_x;
              atom_ptr->y = new_y;
              atom_ptr->z = new_z;

              /* Apply the transformation to all atoms springing off
                 the current bond */
              transform_downstream_atoms(bond_ptr,other_end,matrix);
            }
        }
    }
}
/***********************************************************************

flip_object  -  Flip the current object about the x- or y-axis

***********************************************************************/

void flip_object(struct object *object_ptr,int xflip,int yflip)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         transformation */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Apply the transformation */
          atom_ptr->x = xflip * atom_ptr->x;
          atom_ptr->y = yflip * atom_ptr->y;

          /* Mark atom as being done */
          atom_ptr->checked = TRUE;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

write_pdb_file  -  Print out the co-ordinates in a PDB fashion

***********************************************************************/

void write_pdb_file(void)
{
/* v.4.0--> */
  int interface;
/* <--v.4.0 */

/* v.5.4.5 --> */
  float x, y, z;
/* <-- v.5.4.5 */

  struct coordinate *atom_ptr;

  /* Write out the title (if there is one) */
  printf("\nWriting out the PDB file...\n");
  if (Print_Title[0] != '\0')
    fprintf(ligplot_pdb,"TITLE  %s\n",Print_Title);

/* v.4.0--> */
  /* If this is a DIMPLOT, then write out identifier */
  if (Interface_Plot == TRUE)
    fprintf(ligplot_pdb,"DIMPLOT\n");

  /* Initialise interface identifier */
  interface = 0;
/* <--v.4.0 */

  /* Initialise pointer to start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms in the linked list */
  while (atom_ptr != NULL)
    {
/* v.4.0--> */
      /* If this is a DIMPLOT, then check which interface is being
         written out */
      if (Interface_Plot == TRUE)
        {
          if (atom_ptr->residue_ptr->object_ptr->interface != interface)
            {
              fprintf(ligplot_pdb,"SURFACE %d\n",
                      atom_ptr->residue_ptr->object_ptr->interface);
              interface = atom_ptr->residue_ptr->object_ptr->interface;
            }
        }
/* <--v.4.0 */

/* v.5.0--> */
      /* Write out if not from a dummy object */
      if (atom_ptr->residue_ptr->object_ptr->object_type != DUMMY)
        {
/* <--v.5.0 */

/* v.5.4.5 --> */
	  /* Save the coordinates */
	  x = atom_ptr->x;
	  y = atom_ptr->y;
	  z = atom_ptr->z;
	  if (PDBsum1 == TRUE)
	    {
	      x = atom_ptr->original_x;
	      y = atom_ptr->original_y;
	      z = atom_ptr->original_z;
	    }
/* <-- v.5.4.5 */

      /* Write out the atom record */
      fprintf(ligplot_pdb,"ATOM  %5s ",atom_ptr->atom_number);
      fprintf(ligplot_pdb,"%4s ",atom_ptr->atom_name);
/* v.4.3--> */
/*      fprintf(ligplot_pdb,"%3s %c%5s   %8.3f%8.3f%8.3f%6s%6s%8.3f\n", */
      fprintf(ligplot_pdb,"%3s %c%5s   %8.3f%8.3f%8.3f%6s%6s%8.3f  %2s\n",
/* <--v.4.3 */
              atom_ptr->residue_ptr->res_name,
              atom_ptr->residue_ptr->chain,
              atom_ptr->residue_ptr->res_num,
/* v.5.4.5 --> */
//              atom_ptr->x,
//              atom_ptr->y,
//              atom_ptr->z,
	      x, y, z,
/* <-- v.5.4.5 */
              atom_ptr->occupancy,
              atom_ptr->bvalue,
/* v.4.3--> */
/*            atom_ptr->accessibility); */
              atom_ptr->accessibility,
              atom_ptr->element);
/* <--v.4.3 */
/* v.5.0--> */
        }
/* <--v.5.0 */

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }

/* v.4.0--> */
  /* Write out the maximum accessibility/B-value */
  fprintf(ligplot_pdb,"ENDLIN                                        ");
  fprintf(ligplot_pdb,"                    %8.3f\n",Maximum_accessibility);
 /* <--v.4.0 */
}
/***********************************************************************

flatten_object  -  Transform each rotatable bond in the current object
                   to the origin and flatten bonds either side of it

***********************************************************************/

void flatten_object(struct object *object_ptr)
{
  float coord_store[3][3], matrix[4][3];

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object_bond *object_bond_ptr;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;

  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* If this is a rotatable bond perform unfolding of the object
         around this bond */
      if (bond_ptr->rotatable_bond == TRUE)
        {
          /* Get pointers to the bond's two atoms */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Store the two atoms' coordinates */
          coord_store[0][0] = atom1_ptr->x;
          coord_store[0][1] = atom1_ptr->y;
          coord_store[0][2] = atom1_ptr->z;
          coord_store[1][0] = atom2_ptr->x;
          coord_store[1][1] = atom2_ptr->y;
          coord_store[1][2] = atom2_ptr->z;

          /* Calculate the transformation that will put this bond along
             the x-axis, with atom 1 at the origin and atom 2 pointing
             in the -ve x-direction */
          line_transformation(coord_store,matrix,-PI);

          /* Apply the transformation to all atoms in the current
             object */
          transform_object(object_ptr,matrix);

          /* Perform flattening of bonds attached to atom 1 */
          flatten_bonds(bond_ptr,FIRST);
          
          /* Shift the entire structure onto second of this bond's
             two atoms and flip right over about the y-axis */
          move_object(object_ptr,atom2_ptr->x,atom2_ptr->y,
                                  atom2_ptr->z);
          flip_object(object_ptr,-1,1);
          
          /* Perform flattening of bonds attached to atom 2 */
          flatten_bonds(bond_ptr,SECOND);
        }
      
      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }
}
/***********************************************************************

unroll_structure  -  Unroll the structure of corrent object until
                     z-coords have been flattened completely

***********************************************************************/

void unroll_structure(struct object *object_ptr)
{
  int flat, flat_loop, loop;
  int natoms;
  float centre_x, centre_y, centre_z, max_off, x, y, z;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise variables */
  flat = FALSE;

  /* Loop until object is reasonably flat */
  for (flat_loop = 0; flat_loop < MAX_FLAT_LOOP && flat == FALSE;
       flat_loop++)
    {
      /* Initialise variables for this loop */
      centre_x = 0.0;
      centre_y = 0.0;
      centre_z = 0.0;
      natoms = 0;
      flat = TRUE;
      max_off = 0.0;

      /* Initialise all the atom flags */
      initialise_atoms();

      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;

      /* Loop through all object's bonds */
      while (object_bond_ptr != NULL)
        {
          /* Get the bond pointer */
          bond_ptr = object_bond_ptr->bond_ptr;

          /* Loop over the bond's two atoms */
          for (loop = 0; loop < 2; loop++)
            {
              /* Get pointer to appropriate atom */
              if (loop == 0)
                atom_ptr = bond_ptr->first_atom_ptr;
              else
                atom_ptr = bond_ptr->second_atom_ptr;

              /* Process only if this atom hasn't already been checked */
              if (atom_ptr->checked == FALSE)
                {
                  /* Get the coordinates of this atom */
                  x = atom_ptr->x;
                  y = atom_ptr->y;
                  z = atom_ptr->z;

                  /* Accumulate centre-of-mass data */
                  centre_x = centre_x + x;
                  centre_y = centre_y + y;
                  centre_z = centre_z + z;
                  natoms++;

                  /* If the z-coordinate is close to zero, set it
                     to exactly zero */
                  if (fabs(z) < 0.001)
                    atom_ptr->z = 0.0;

                  /* Otherwise, object is not flat, so check if this is
                     the highest deviant from the x-y plane so far */
                  else
                    {
                      flat = FALSE;
                      if (fabs(z) > max_off)
                        max_off = (float) (fabs(z));
                    }

                  /* Set this atom's flag as checked */
                  atom_ptr->checked = TRUE;
                }
            }
      
          /* Get pointer to this object's next bond entry */
          object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
        }

      /* If the object needs flattening, then do it */
      if (flat == FALSE)
        {
          /* If object contains rotatable bonds, then can unroll it */
          if (object_ptr->nrot_bonds > 0)
            flatten_object(object_ptr);

          /* Otherwise, just move object to origin */
          else
            if (natoms > 0)
              {
                centre_x = centre_x / natoms;
                centre_y = centre_y / natoms;
                centre_z = centre_z / natoms;
                move_object(object_ptr,centre_x,centre_y,centre_z);
              }
        }
    }
}
/***********************************************************************

adjust_to_origin  -  Translate all the atoms by the given amount

***********************************************************************/

void adjust_to_origin(float centre_x,float centre_y,float centre_z)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Apply the transformation */
      atom_ptr->x = atom_ptr->x - centre_x;
      atom_ptr->y = atom_ptr->y - centre_y;
      atom_ptr->z = atom_ptr->z - centre_z;
      
      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }          
}
/* v.5.3.4--> */
/***********************************************************************

adjust_object_to_origin  -  Translate all the atoms of the given object
                            by the given amount

***********************************************************************/

void adjust_object_to_origin(struct object *object_ptr,float centre_x,
                             float centre_y,float centre_z)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         transformation */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Apply the transformation */
          atom_ptr->x = atom_ptr->x - centre_x;
          atom_ptr->y = atom_ptr->y - centre_y;
          atom_ptr->z = atom_ptr->z - centre_z;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/* <--v.5.3.4 */
/***********************************************************************

align_object  -  Align object's principal axes with the x-y plane

***********************************************************************/

void align_object(struct object *object_ptr)
{
/* v.4.0--> */
/*  int iobject, iatom, iresid, natoms, nresid, nstored;
  int wanted; */
  int iatom, iresid, natoms, nresid, nstored;
/* <--v.4.0 */

  float coord_store[MAXATOMS][3], mass_x, mass_y, mass_z;
/* v.4.0--> */
/*  float new_x, new_y, new_z, x, y, z;
  float x, y, z; */
/* <--v.4.0 */
  float transformation[3][3];

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  nstored = 0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to first residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Store this atom's fitted coords */
          if (nstored < MAXATOMS)
            {
              coord_store[nstored][0] = atom_ptr->x;
              coord_store[nstored][1] = atom_ptr->y;
              coord_store[nstored][2] = atom_ptr->z;
              nstored++;
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Calculate the centre of mass of these atoms */
  centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

  /* Place the centre of mass at the origin */
/* v.5.3.4--> */
//  adjust_to_origin(mass_x,mass_y,mass_z);
  adjust_object_to_origin(object_ptr,mass_x,mass_y,mass_z);
/* <--v.5.3.4 */

  /* Apply the same translation to the coordinates stored in
     coord_store */
  adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

  /* Calculate transformation to orient the principal axis along
     the x-direction */
  principal_components(nstored,transformation,coord_store);

  /* Apply this transformation to the entire object */
  rotate_object(object_ptr,transformation);
}
/***********************************************************************

force_flat  -  Zero all the z-coordinates (whether or not the structure
               has been properly flattened) so that the flipping
               operations aren't thrown out by any non-planarity in
               the structure

***********************************************************************/

void force_flat(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         zeroing the z-coords */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Set the z-coord to zero */
          atom_ptr->z = 0.0;
      
          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

create_simple_hgroup  -  Adjust all the non-ligand residues for the
                         schematic plot, marking any unwanted atoms
                         for deletion

***********************************************************************/

void create_simple_hgroup(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  float hgroup_size;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Calculate size for the simple H-group */
  hgroup_size
    = (float) (5.0 * Text_Size_Val->Nonligand_Residue_Names * CHAR_ASPECT);

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is an H-bond */
/* v.5.4--> */
//      if (bond_ptr->bond_type == HBOND)
      if (bond_ptr->bond_type == HBOND ||
          bond_ptr->bond_type == SALT_BRIDGE ||
          bond_ptr->bond_type == DISULPHIDE)
/* <--v.5.4 */
        {
          /* Mark both atoms as being involved in the bond */
          bond_ptr->first_atom_ptr->checked = TRUE;
          bond_ptr->second_atom_ptr->checked = TRUE;
        }

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Set object type to simplified residue */
  object_ptr->object_type = SIMPLE_HGROUP;

  /* Adjust object- and residue- maximum atom sizes */
  object_ptr->max_atom_size = hgroup_size;
  residue_ptr->max_atom_size = hgroup_size;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      
      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, marking any not involved
         in H-bonds for deletion */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* If atom not involved in any interactions with other
             residues, then can mark it for deletion */
          if (atom_ptr->checked == FALSE)
            atom_ptr->deleted = TRUE;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

remove_mainchains  -  Remove any non-ligand mainchain atoms if mainchain
                      not involved in any interactions

***********************************************************************/

void remove_mainchains(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  int ninteract;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is an H-bond or a contact bond */
/* v.5.4--> */
//      if (bond_ptr->bond_type == HBOND || bond_ptr->bond_type == CONTACT)
      if (bond_ptr->bond_type == HBOND ||
          bond_ptr->bond_type == SALT_BRIDGE ||
          bond_ptr->bond_type == DISULPHIDE ||
          bond_ptr->bond_type == CONTACT)
/* <--v.5.4 */
        {
          /* Mark both atoms as being involved in the bond */
          bond_ptr->first_atom_ptr->checked = TRUE;
          bond_ptr->second_atom_ptr->checked = TRUE;
        }

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Initialise count of mainchain atoms involved in interactions */
      ninteract = 0;

      /* Loop over all the residue's atoms counting the number of
         mainchain atoms involved in interactions */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* If this is a mainchain atom and it is involved in
             an interaction, then increment count */
          if (atom_ptr->side_chain == FALSE && atom_ptr->checked == TRUE)
            ninteract++;

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* If no interactions with any of the mainchain atoms, then
         delete them all */
      if (ninteract == 0)
        {
          /* Re-initialise to this residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;

          /* Loop over all the residue's atoms, marking mainchain
             atoms for deletion */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* If this is a mainchain atom, mark it for deletion */
              if (atom_ptr->side_chain == FALSE)
                atom_ptr->deleted = TRUE;

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

mark_sidechain_atoms  -  Mark all the current residue's sidechain atoms
                         for deletion

***********************************************************************/

void mark_sidechain_atoms(struct residue *residue_ptr)
{
  int iatom, natoms;

  struct coordinate *atom_ptr;

  /* Get pointer to this residue's first atom */
  atom_ptr = residue_ptr->first_atom_ptr;

  /* Initialise atom count */
  iatom = 0;
  natoms = residue_ptr->natoms;

  /* Loop over all the residue's atoms, to mark all the sidechain atoms  */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* Initialise the side-chain flag */
      if (atom_ptr->side_chain == TRUE)
        atom_ptr->deleted = TRUE;

      /* Get pointer to the next atom in this residue */
      atom_ptr = atom_ptr->next;
      iatom++;
    }
}
/***********************************************************************

transform_downlist_atoms  -  Apply the given transformation to all the
                             atoms in the list from the given CA atom
                             onwards by finding the first downlist bond

***********************************************************************/

void transform_downlist_atoms(struct object *object_ptr,
                              struct coordinate *current_atom_ptr,
                              float matrix[4][3])
{
  int end, iatom, iresid, natoms, nresid;
  int got_bond, transform, wanted;

  struct bond *bond_ptr, *use_bond_ptr;
  struct coordinate *atom_ptr;
  struct residue *current_residue_ptr, *residue_ptr;

  /* Initialise variables */
  got_bond = FALSE;
  transform = FALSE;
  use_bond_ptr = NULL;
  current_residue_ptr = current_atom_ptr->residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL && got_bond == FALSE)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the
         transformation where required */
      while (iatom < natoms && atom_ptr != NULL && got_bond == FALSE)
        {
          /* Check whether this is the starting atom */
          if (atom_ptr == current_atom_ptr)
            transform = TRUE;

          /* Check if this is a sidechain atom of the starting
             residue */
          wanted = TRUE;
          if (residue_ptr == current_residue_ptr &&
              atom_ptr->side_chain == TRUE)
            wanted = FALSE;

          /* If have a potential starting bond, then try to find it */
          if (transform == TRUE && wanted == TRUE)
            {
              /* Initialise to first bond */
              bond_ptr = first_bond_ptr;

              /* Loop through all the bonds to find a covalent bond
                 between the two atoms */
              while (bond_ptr != NULL && got_bond == FALSE)
                {
                  /* Check if this is a covalent bond between the 
                     two atoms */
                  if (bond_ptr->bond_type == COVALENT)
                    {
                      /* Atoms in bond, with the far end the one
                         we're after */
                      if (bond_ptr->first_atom_ptr == current_atom_ptr &&
                          bond_ptr->second_atom_ptr == atom_ptr)
                        {
                          got_bond = TRUE;
                          use_bond_ptr = bond_ptr;
                          end = SECOND;
                        }

                      /* Atoms in bond, with the far end the one
                         we're after */
                      else if (bond_ptr->first_atom_ptr == atom_ptr &&
                               bond_ptr->second_atom_ptr == current_atom_ptr)
                        {
                          got_bond = TRUE;
                          use_bond_ptr = bond_ptr;
                          end = FIRST;
                        }
                    }

                  /* Go to the next bond in the linked list */
                  bond_ptr = bond_ptr->next_bond_ptr;
                }
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* If have a bond that can be used for the transformation, use it */
  if (use_bond_ptr != NULL)
    transform_downstream_atoms(use_bond_ptr,end,matrix);
}
/***********************************************************************

stretch_mainchain  -  Stretch the mainchain of the current residue out
                      into a straight line for the simplified ligand
                      representation

***********************************************************************/

void stretch_mainchain(struct object *object_ptr,
                       struct coordinate *n_ptr,struct coordinate *ca_ptr,
                       struct coordinate *c_ptr,struct coordinate *o_ptr,
                       struct coordinate *previous_ca_ptr,
                       struct coordinate *last_n_ptr,
                       struct coordinate *last_ca_ptr,
                       struct coordinate *last_c_ptr,
                       struct coordinate *last_o_ptr)
{
  int ncovalent, negside;
  int lose_C, nhbonds_O;
  int iatom, natoms;

  float coord_store[3][3], matrix[4][3];
  float x, x1, x2, y, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  lose_C = FALSE;

  /* Adjust the atom-sizes of the main-chain atoms */
  ca_ptr->atom_size = Size_Val->Simple_Residues;
  n_ptr->atom_size = Size_Val->Ligand_Bonds;
  c_ptr->atom_size = Size_Val->Ligand_Bonds;
  ca_ptr->plot_atom = FALSE;
  n_ptr->plot_atom = FALSE;
  c_ptr->plot_atom = FALSE;

  /* Update residue and object maximum sizes */
  if (ca_ptr->atom_size > ca_ptr->residue_ptr->object_ptr->max_atom_size)
    ca_ptr->residue_ptr->object_ptr->max_atom_size = ca_ptr->atom_size;
  if (ca_ptr->atom_size > ca_ptr->residue_ptr->max_atom_size)
    ca_ptr->residue_ptr->max_atom_size = ca_ptr->atom_size;

  /* Increment count of object's CA atoms */
  object_ptr->n_ca++;

  /* If have the previous CA or last N atoms, then can stretch the mainchain
     by pivoting at the last CA and getting the previous CA (or last N) in
     line with the current CA */
  if (last_ca_ptr != NULL && (previous_ca_ptr != NULL || last_n_ptr != NULL))
    {
      /* Store the coords of the last CA */
      coord_store[0][0] = last_ca_ptr->x;
      coord_store[0][1] = last_ca_ptr->y;
      coord_store[0][2] = last_ca_ptr->z;

      /* Store the coords of the previous CA, or last N, depending on
         which are available */
      if (previous_ca_ptr != NULL)
        {
          coord_store[1][0] = previous_ca_ptr->x;
          coord_store[1][1] = previous_ca_ptr->y;
          coord_store[1][2] = previous_ca_ptr->z;
        }
      else
        {
          coord_store[1][0] = last_n_ptr->x;
          coord_store[1][1] = last_n_ptr->y;
          coord_store[1][2] = last_n_ptr->z;
        }

      /* Calculate the transformation that will place the last CA at the
         origin, with the previous CA (or last N) in the negative
         x-direction */
      line_transformation(coord_store,matrix,-PI);

      /* Apply the transformation to all atoms in the current object */
      transform_object(object_ptr,matrix);

      /* Store the coords of the current CA atom */
      coord_store[0][0] = 0.0;
      coord_store[0][1] = 0.0;
      coord_store[0][2] = 0.0;
      coord_store[1][0] = ca_ptr->x;
      coord_store[1][1] = ca_ptr->y;
      coord_store[1][2] = ca_ptr->z;

      /* Calculate transformation to rotate about the last CA
         to put the current CA atom along the +ve x-axis */
      line_transformation(coord_store,matrix,0.0);

      /* Apply the transformation to all the atoms in the object
         from the last CA onwards */
      transform_downlist_atoms(object_ptr,last_ca_ptr,matrix);
    }

  /* If have a previous CA and C atoms, then stretch the N and O residues
     equidistantly between the two */
  if (last_ca_ptr != NULL)
    {
      /* Get the coordinates of the 2 CA atoms */
      x1 = last_ca_ptr->x;
      y1 = last_ca_ptr->y;
      x2 = ca_ptr->x;
      y2 = ca_ptr->y;

      /* Get point one third of the way along */
      x = (float) (x1 + (x2 - x1) / 3.0);
      y = (float) (y1 + (y2 - y1) / 3.0);

      /* Place the previous carbon there */
      if (last_c_ptr != NULL)
        {
          last_c_ptr->x = x;
          last_c_ptr->y = y;

          /* If there was an associated carbonyl oxygen, then place it 
             on one side */
          if (last_o_ptr != NULL)
            {
              last_o_ptr->x = x;
              last_o_ptr->y = (float) (y + 1.23);
            }
        }

      /* Get point two thirds of the way along */
      x = (float) (x1 + 2.0 * (x2 - x1) / 3.0);
      y = (float) (y1 + 2.0 * (y2 - y1) / 3.0);

      /* Place the current residue's nitrogen there */
      n_ptr->x = x;
      n_ptr->y = y;

      /* If residue is a proline, then need to tidy it up a little */
      if (!strncmp(n_ptr->residue_ptr->res_name,"PRO",3))
        {
          /* Get the Pro's first atom */
          atom_ptr = n_ptr->residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = n_ptr->residue_ptr->natoms;
          negside = 0;

          /* Loop over the residue's atoms to adjust the sidechain
             coords */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* If this is a sidechain atom, then fix its coords */
              if (!strncmp(atom_ptr->atom_name," CB ",4))
                {
                  /* Check which side of the x-axis this atom lies on */
                  if (negside == 0)
                    {
                      if (atom_ptr->y < 1)
                        {
                          negside = -1;
                        }
                      else
                        {
                          negside = 1;
                        }
                    }
                  atom_ptr->x = (float) (ca_ptr->x + 0.400);
                  atom_ptr->y = (float) (negside * 1.500);
                  atom_ptr->z = 0.0;
                }

              if (!strncmp(atom_ptr->atom_name," CG ",4))
                {
                  /* Check which side of the x-axis this atom lies on */
                  if (negside == 0)
                    {
                      if (atom_ptr->y < 1)
                        {
                          negside = -1;
                        }
                      else
                        {
                          negside = 1;
                        }
                    }
                  atom_ptr->x = (float) ((n_ptr->x + ca_ptr->x) / 2.0);
                  atom_ptr->y = (float) (negside * 2.300);
                  atom_ptr->z = 0.0;
                }

              if (!strncmp(atom_ptr->atom_name," CD ",4))
                {
                  /* Check which side of the x-axis this atom lies on */
                  if (negside == 0)
                    {
                      if (atom_ptr->y < 1)
                        {
                          negside = -1;
                        }
                      else
                        {
                          negside = 1;
                        }
                    }
                  atom_ptr->x = (float) (n_ptr->x - 0.400);
                  atom_ptr->y = (float) (negside * 1.500);
                  atom_ptr->z = 0.0;
                }

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }
        }
    }

  /* Place the carbonyl oxygen at the carbonyl carbon position and 
     alter the bonds, if required */
  if (lose_C == TRUE && o_ptr != NULL)
    {
      /* Replace the carbonyl carbon by the oxygen and delete the
         oxygen */
      strncpy(c_ptr->atom_name," O  ",4);
      strncpy(c_ptr->print_name," O  ",4);
      c_ptr->atom_name[4] = '\0';
      c_ptr->print_name[4] = '\0';
      o_ptr->deleted = TRUE;
    }

  /* Blank out the name of the CA atom */
  ca_ptr->print_name[0] = '\0';
  ncovalent = 0;

  /* Loop through all bonds to transfer any H-bonds from the oxygen
     to its new position on the carbonyl carbon (if required), or delete
     the oxygen if not involved in any H-bonds, and to transfer any
     non-bonded contacts from the carbon to the CA */
  nhbonds_O = 0;

  /* Get pointer to first bond in the linked list */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* If this is an H-bond, check whether the O is involved */
/* v.5.4--> */
//      if (bond_ptr->bond_type == HBOND)
      if (bond_ptr->bond_type == HBOND ||
          bond_ptr->bond_type == SALT_BRIDGE ||
          bond_ptr->bond_type == DISULPHIDE)
/* <--v.5.4 */
        {
          /* Check for H-bonds to the current O atom */
          if (bond_ptr->first_atom_ptr == o_ptr ||
              bond_ptr->second_atom_ptr == o_ptr)
            {
              /* Increment count of H-bonds to the O */
              nhbonds_O++;

          /* If either atom of the bond is the oxygen which
             has just been transferred, then transfer pointer, too */
              if (lose_C == TRUE)
                {
                  if (bond_ptr->first_atom_ptr == o_ptr)
                    bond_ptr->first_atom_ptr = c_ptr;
                  if (bond_ptr->second_atom_ptr == o_ptr)
                    bond_ptr->second_atom_ptr = c_ptr;
                }
            }
        }

      /* Transfer any contacts with the C to the CA */
      else if (bond_ptr->bond_type == CONTACT && lose_C == TRUE)
        {
          /* If either atom of the bond is the carbonyl carbon
             then transfer to the CA */
          if (bond_ptr->first_atom_ptr == c_ptr)
            bond_ptr->first_atom_ptr = ca_ptr;
          else if (bond_ptr->second_atom_ptr == c_ptr)
            bond_ptr->second_atom_ptr = ca_ptr;
        }

      /* Count the number of covalent bonds to the C */
      else if (bond_ptr->bond_type == COVALENT)
        {
          /* If either atom of the bond is the carbon which
             has just been overwritten, then increment count of
             its covalent bonds */
          if (bond_ptr->first_atom_ptr == c_ptr ||
              bond_ptr->second_atom_ptr == c_ptr)
            ncovalent++;
        }

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* If oxygen has no hydrogen bonds, then can be deleted */
  if (lose_C == FALSE && nhbonds_O == 0 && o_ptr != NULL)
    o_ptr->deleted = TRUE;
}
/***********************************************************************

create_simple_ligand  -  Adjust all the ligand residues for the
                         schematic plot, marking any unwanted atoms
                         for deletion

***********************************************************************/

void create_simple_ligand(struct object *object_ptr)
{
  int iatom, iresid, natoms, nresid;
  int sidechain_interaction;

  struct bond *bond_ptr;
  struct coordinate *c_ptr, *ca_ptr, *n_ptr, *o_ptr;
  struct coordinate *last_ca_ptr, *last_n_ptr, *last_c_ptr, *last_o_ptr;
  struct coordinate *previous_ca_ptr;
  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise all atoms */
  initialise_atoms();

  /* Initialise bond pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds */
  while (bond_ptr != NULL)
    {
      /* Process only if this is an H-bond or hydrophobic contact */
/* v.5.4--> */
//      if (bond_ptr->bond_type == HBOND || bond_ptr->bond_type == CONTACT)
      if (bond_ptr->bond_type == HBOND ||
          bond_ptr->bond_type == SALT_BRIDGE ||
          bond_ptr->bond_type == DISULPHIDE ||
          bond_ptr->bond_type == CONTACT)
/* <--v.5.4 */
        {
          /* Mark both atoms as being involved in the bond */
          bond_ptr->first_atom_ptr->checked = TRUE;
          bond_ptr->second_atom_ptr->checked = TRUE;
        }

      /* Get pointer for next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Loop over the object's residues to mark any unwanted sidechain
     atoms for deletion */

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Initialise flag */
      sidechain_interaction = FALSE;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, to locate all the mainchain
         and sidechain atoms  */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Initialise sidechain marker */
          atom_ptr->side_chain = FALSE;

          /* Check if atom is one of the mainchain atoms */
          if (!strncmp(atom_ptr->atom_name," N  ",4))
            n_ptr = atom_ptr;
          else if (!strncmp(atom_ptr->atom_name," CA ",4))
            ca_ptr = atom_ptr;
          else if (!strncmp(atom_ptr->atom_name," C  ",4))
            c_ptr = atom_ptr;
          else if (!strncmp(atom_ptr->atom_name," O  ",4))
            o_ptr = atom_ptr;

          /* Otherwise, if it is a sidechain atom, then see whether
             it is involved in an interaction */
          else
            {
              atom_ptr->side_chain = TRUE;
              if (atom_ptr->checked == TRUE)
                sidechain_interaction = TRUE;
            }

          /* Get pointer to the next atom in this residue */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* If have all the mainchain atoms, then can stretch this
         residue out into a straight line */
      if (n_ptr != NULL && ca_ptr != NULL && c_ptr != NULL)
        {
          /* If none of the sidechain atoms are involved in interactions,
             then mark the whole sidechain for deletion */
          if (sidechain_interaction == FALSE)
            mark_sidechain_atoms(residue_ptr);
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Loop to stretch out the mainchain */

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Initialise pointer to previous residue's CA */
  previous_ca_ptr = NULL;
  last_ca_ptr = NULL;
  last_n_ptr = NULL;
  last_c_ptr = NULL;
  last_o_ptr = NULL;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Initialise flags */
      n_ptr = NULL;
      ca_ptr = NULL;
      c_ptr = NULL;
      o_ptr = NULL;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, to locate all the mainchain
         and sidechain atoms  */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Check if atom is one of the mainchain atoms */
          if (!strncmp(atom_ptr->atom_name," N  ",4))
            n_ptr = atom_ptr;
          else if (!strncmp(atom_ptr->atom_name," CA ",4))
            ca_ptr = atom_ptr;
          else if (!strncmp(atom_ptr->atom_name," C  ",4))
            c_ptr = atom_ptr;
          else if (!strncmp(atom_ptr->atom_name," O  ",4))
            o_ptr = atom_ptr;

          /* Get pointer to the next atom in this residue */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* If have all the mainchain atoms, then can stretch this
         residue out into a straight line */
      if (n_ptr != NULL && ca_ptr != NULL && c_ptr != NULL)
        {
          /* Stretch the mainchain out */
          stretch_mainchain(object_ptr,n_ptr,ca_ptr,c_ptr,o_ptr,
                            previous_ca_ptr,last_n_ptr,last_ca_ptr,
                            last_c_ptr,last_o_ptr);

          /* Mark the residue as being of the simplified type */
          residue_ptr->residue_type = SIMPLE_LIGAND;

          /* Save the pointers (note that O's position is saved as the
             current carbonyl C) */
          previous_ca_ptr = last_ca_ptr;
          last_ca_ptr = ca_ptr;
          last_n_ptr = n_ptr;
          last_c_ptr = c_ptr;
          last_o_ptr = o_ptr;
        }

      /* Otherwise, no stretching to be done */
      else
        {
          last_n_ptr = NULL;
          last_c_ptr = NULL;
          last_o_ptr = NULL;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

Andrew Martin's matfit and qikfit routines for least-squares fitting of
one set of coordinates onto another (both sets centred around the origin)

***********************************************************************/
/*>static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column)
   ---------------------------------------------------------------
   Input:   REAL  umat[3][3]     The U matrix
            BOOL  column         TRUE: Create a column-wise matrix
                                 (other way round from normal).
   Output:  REAL  rm[3][3]       The output rotation matrix
  
   Does the actual fitting for matfit().
   04.02.91 Original
   01.06.92 ANSIed & doc'd
   11.03.94 column changed to BOOL
*/
static void qikfit(REAL  umat[3][3],
                   REAL  rm[3][3],
                   BOOL  column)
{
   
   REAL  rot[3][3],
         turmat[3][3],
         c[3][3],
         coup[3],
         dir[3],
         step[3],
         v[3],
         rtsum,rtsump,
         stp,stcoup,
         ud,tr,ta,cs,sn,ac,
         delta,
         gfac,
         cle,clep;
   int   i,j,k,l,m,
         jmax,
         ncyc,
         nsteep,
         nrem;

   /* Rotate repeatedly to reduce couple about initial direction
      to zero.
      Clear the rotation matrix */
   for(l=0;l<3;l++)
   {
      for(m=0;m<3;m++)
         rot[l][m] = 0.0;
      rot[l][l] = 1.0;
   }

   /* Copy vmat[][] (sp) into umat[][] (dp) */
   jmax = 30;
   rtsum = umat[0][0] + umat[1][1] + umat[2][2];
   delta = 0.0;

   for(ncyc=0;ncyc<jmax;ncyc++)
   {
      /* Modified CG. For first and every NSTEEP cycles, set previous
         step as zero and do an SD step */
      nsteep = 3;
      nrem = ncyc-nsteep*(int)(ncyc/nsteep);

      if(!nrem)
      {
         for(i=0;i<3;i++) step[i]=0.0;
         clep = 1.0;
      }
      
      /* Couple */
      coup[0] = umat[1][2]-umat[2][1];
      coup[1] = umat[2][0]-umat[0][2];
      coup[2] = umat[0][1]-umat[1][0];
      cle = (float) sqrt(coup[0]*coup[0] + coup[1]*coup[1] + coup[2]*coup[2]);

      /* Gradient vector is now -coup */
      gfac = (cle/clep)*(cle/clep);

      /* Value of rtsum from previous step */
      rtsump = rtsum;
      clep   = cle;
      if(cle < SMALL) break;

      /* Step vector conjugate to  previous */
      stp = 0.0;
      for(i=0;i<3;i++)
      {
         step[i]=coup[i]+step[i]*gfac;
         stp   += (step[i] * step[i]);
      }
      stp = 1.0/sqrt(stp);
         
      /* Normalised step */
      for(i=0;i<3;i++) dir[i] = stp*step[i];

      /* Couple resolved along step direction */
      stcoup = coup[0]*dir[0] + coup[1]*dir[1] + coup[2]*dir[2];

      /* Component of UMAT along direction */
      ud = 0.0;
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            ud += umat[l][m]*dir[l]*dir[m];


      tr = umat[0][0]+umat[1][1]+umat[2][2]-ud;
      ta = (float) sqrt(tr*tr + stcoup*stcoup);
      cs=tr/ta;
      sn=stcoup/ta;
         
      /* If cs<0 then posiiton is unstable, so don't stop */
      if((cs>0.0) && (ABS(sn)<SMALSN)) break;
            
      /* Turn matrix for correcting rotation:

         Symmetric part */
      ac = 1.0-cs;
      for(l=0;l<3;l++)
      {
         v[l] = ac*dir[l];
         for(m=0;m<3;m++)
            turmat[l][m] = v[l]*dir[m];
         turmat[l][l] += cs;
         v[l]=dir[l]*sn;
      }

      /* Asymmetric part */
      turmat[0][1] -= v[2];
      turmat[1][2] -= v[0];
      turmat[2][0] -= v[1];
      turmat[1][0] += v[2];
      turmat[2][1] += v[0];
      turmat[0][2] += v[1];

      /* Update total rotation matrix */
      for(l=0;l<3;l++)
      {
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*rot[k][m];
         }
      }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            rot[l][m] = c[l][m];

      /* Update umat tensor */
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*umat[k][m];
         }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            umat[l][m] = c[l][m];

      rtsum = umat[0][0] + umat[1][1] + umat[2][2];
      delta = rtsum - rtsump;

      /* If no improvement in this cycle then stop */
      if(ABS(delta)<SMALL) break;

      /* Next cycle */
   }

   /* Copy rotation matrix for output */
   if(column)
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[j][i] = rot[i][j];
   }
   else
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[i][j] = rot[i][j];
   }
}
/**********************************************************************/
/**********************************************************************/
/*>BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n,
               REAL *wt1, BOOL column)
   -----------------------------------------------------
   Input:   COOR  *x1         First array of coordinates
            COOR  *x2         Second array of coordinates
            int   n           Number of coordinates
            REAL  *wt1        Weight array or NULL
            int   column      TRUE: Output a column-wise matrix (as used
                                 by FRODO)
                              FALSE: Output a standard row-wise matrix.
   Output:  REAL  rm[3][3]    Returned rotation matrix
   Returns: rmsd              RMS distance between the two sets of
                              coordinates after fitting

   Fit coordinate array x1 to x2 both centred around the origin and of 
   length n. Optionally weighted with the wt1 array if iw flag is set.
   If flag is set, error messages will be issued internally. 
   If column is set the matrix will be returned column-wise rather
   than row-wise.

   04.02.91 Original
   01.06.92 ANSIed & doc'd
   17.06.93 various changes for release (including parameters)
   11.03.94 column changed to BOOL
   04.03.97 Calc of rmsd added (RAL)
*/
float matfit(COOR    *x1,        /* First coord array    */
             COOR    *x2,        /* Second coord array   */
             REAL    rm[3][3],   /* Rotation matrix      */
             int     n,          /* Number of points     */
             REAL    *wt1,       /* Weight array         */
             BOOL    column)     /* Column-wise output   */
{
   int  i, ipoint, j;
   REAL dist, rmsd, x, y, z;

   REAL umat[3][3];

   if(wt1)
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += wt1[j] * x1[j].x * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].x * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += wt1[j] * x1[j].y * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].y * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += wt1[j] * x1[j].z * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].z * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   } 
   else
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += x1[j].x * x2[j].x;
                  umat[i][1] += x1[j].x * x2[j].y;
                  umat[i][2] += x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += x1[j].y * x2[j].x;
                  umat[i][1] += x1[j].y * x2[j].y;
                  umat[i][2] += x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += x1[j].z * x2[j].x;
                  umat[i][1] += x1[j].z * x2[j].y;
                  umat[i][2] += x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   }
   qikfit(umat,rm,column);

   /* Calculate rmsd between the two fitted sets of coordinates */
   rmsd = 0.0;

   /* Loop over all the pairs of coordinates */
   for (ipoint = 0; ipoint < n; ipoint++)
     {
       /* Calculate the transformed coordinates of the current point */
       if (column)
         {
           x = rm[0][0] * x1[ipoint].x + rm[1][0] * x1[ipoint].y
             + rm[2][0] * x1[ipoint].z;
           y = rm[0][1] * x1[ipoint].x + rm[1][1] * x1[ipoint].y
             + rm[2][1] * x1[ipoint].z;
           z = rm[0][2] * x1[ipoint].x + rm[1][2] * x1[ipoint].y
             + rm[2][2] * x1[ipoint].z;
         }
       else
         {
           x = rm[0][0] * x1[ipoint].x + rm[0][1] * x1[ipoint].y
             + rm[0][2] * x1[ipoint].z;
           y = rm[1][0] * x1[ipoint].x + rm[1][1] * x1[ipoint].y
             + rm[1][2] * x1[ipoint].z;
           z = rm[2][0] * x1[ipoint].x + rm[2][1] * x1[ipoint].y
             + rm[2][2] * x1[ipoint].z;
         }

       /* Calculate squared distance between these two points */
       dist = (x - x2[ipoint].x) * (x - x2[ipoint].x)
         + (y - x2[ipoint].y) * (y - x2[ipoint].y)
         + (z - x2[ipoint].z) * (z - x2[ipoint].z);
       rmsd = rmsd + dist;
     }

   /* Calculate the rmsd */
   if (n > 0)
     {
       rmsd = rmsd / n;
       if (rmsd > 0.0)
         rmsd = sqrt(rmsd);
       else
         rmsd = 0.0;
     }

   /* Return the calculated rmsd */
   return((float) rmsd);
}
/***********************************************************************

get_object_distance_match  -  Calculate the rms deviation between the
                              atom-atom distances within the current
                              object and the equivalent distances in the
                              original coordinates system

***********************************************************************/

float get_object_distance_match(struct object *object_ptr)
{
  int iatom, iresid, jatom, jresid, natoms1, nresid1, nresid2, natoms2;
  int nterms;
  int atom_interval;

  float distance, dist2, orig_distance, rmsd, total_rmsd;
  float x1, x2, y1, y2;
  float orig_x1, orig_x2, orig_y1, orig_y2, orig_z1, orig_z2;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  total_rmsd = 0.0;

  /* Set the pointer to the object's first residue */
  residue1_ptr = object_ptr->first_residue_ptr;
  nresid1 = object_ptr->nresidues;
  iresid = 0;
  nterms = 0;

  /* If have only one residue, use every atom, for two residues use every
     other atom, and so on */
  atom_interval = nresid1;

  /* Loop over all this object's residues */
  while (iresid < nresid1 && residue1_ptr != NULL)
    {
      /* Set second pointer to the object's first residue */
      residue2_ptr = object_ptr->first_residue_ptr;
      nresid2 = object_ptr->nresidues;
      jresid = 0;

      /* Loop over all the object's residues for the second atoms */
      while (jresid < nresid2 && residue2_ptr != NULL)
        {
          /* Get pointer to residue's first atom */
          atom1_ptr = residue1_ptr->first_atom_ptr;

          /* Initialise atom count */
          iatom = 0;
          natoms1 = residue1_ptr->natoms;

          /* Loop over all the residue's atoms */
          while (iatom < natoms1 && atom1_ptr != NULL)
            {
              /* If this is a wanted atom, then carry on */
              if (iatom % atom_interval == 0)
                {
                  /* Get the atom's coordinates */
                  x1 = atom1_ptr->x;
                  y1 = atom1_ptr->y;
                  orig_x1 = atom1_ptr->original_x;
                  orig_y1 = atom1_ptr->original_y;
                  orig_z1 = atom1_ptr->original_z;

                  /* Get pointer to first compare atom */
                  atom2_ptr = residue2_ptr->first_atom_ptr;

                  /* Initialise atom count */
                  jatom = 0;
                  natoms2 = residue2_ptr->natoms;

                  /* Loop over all the atoms making up the compare atom */
                  while (jatom < natoms2 && atom2_ptr != NULL)
                    {
                      /* If this is a wanted atom, then carry on */
                      if (jatom % atom_interval == 0)
                        {
                          /* Skip if dealing with the same atom */
                          if (atom1_ptr != atom2_ptr)
                            {
                              /* Get atom's coordinates */
                              x2 = atom2_ptr->x;
                              y2 = atom2_ptr->y;
                              orig_x2 = atom2_ptr->original_x;
                              orig_y2 = atom2_ptr->original_y;
                              orig_z2 = atom2_ptr->original_z;

                              /* Calculate the distance between the
                                 two atoms */
                              dist2 = (x1 - x2) * (x1 - x2)
                                + (y1 - y2) * (y1 - y2);
                              distance = (float) sqrt((double) dist2);

                              /* Calculate the original distance between
                                 these two atoms in the structure */
                              dist2 = (orig_x1 - orig_x2) * (orig_x1 - orig_x2)
                                + (orig_y1 - orig_y2) * (orig_y1 - orig_y2)
                                  + (orig_z1 - orig_z2) * (orig_z1 - orig_z2);
                              orig_distance = (float) sqrt((double) dist2);

                              /* Calculate the rmsd between the two
                                 distances */
                              rmsd = distance - orig_distance;
                              rmsd = rmsd * rmsd;

                              /* Add to total rmsd value */
                              total_rmsd = total_rmsd + rmsd;
                              nterms++;
                            }
                        }

                      /* Get pointer to other residue's next atom */
                      atom2_ptr = atom2_ptr->next;
                      jatom++;
                    }
                }

              /* Get pointer to the next atom */
              atom1_ptr = atom1_ptr->next;
              iatom++;
            }

          /* Get pointer to the next residue */
          residue2_ptr = residue2_ptr->next_residue_ptr;
          jresid++;
        }

      /* Get pointer to the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return the calculated energy */
  if (nterms > 0)
    total_rmsd = total_rmsd / nterms;
  if (total_rmsd > 0.0)
    total_rmsd = (float) sqrt(total_rmsd);
  else
    total_rmsd = 0.0;
  return(total_rmsd);
}
/***********************************************************************

fit_object  -  Fit the given object as best as possible to the 
               original version (to give approximate orientation and
               relative disposition of the component atoms)

***********************************************************************/

float fit_object(struct object *object_ptr,int no_move)
{
  int i, iatom, iresid, istored, j;
  int natoms, nresid, nstored;
  int best_flipped;
/* v.4.0--> */
  int internal_fit, wanted;
/* <--v.4.0 */

  float transform[3][3];
  float mean_x, mean_y, mean_fit_x, mean_fit_y;
  float flipped_rmsd, rmsd;
  float distance_match, flipped_distance_match;
/* v.4.0--> */
  float residue_mean_x, residue_mean_y;
/* <--v.4.0 */

  BOOL column;

  REAL  flip_matrix[3][3], matrix[3][3];

  COOR coords1[MAXATOMS], coords1_flip[MAXATOMS], coords2[MAXATOMS];

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
/* v.4.0--> */
  if (no_move == TRUE)
    internal_fit = TRUE;
  else
    internal_fit = FALSE;
/* <--v.4.0 */
  best_flipped = FALSE;
  mean_x = 0.0;
  mean_y = 0.0;
  mean_fit_x = 0.0;
  mean_fit_y = 0.0;
  nstored = 0;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
/* v.4.0--> */
      /* Initialise residue mean positions */
      residue_mean_x = 0.0;
      residue_mean_y = 0.0;
/* <--v.4.0 */

      /* Get pointer to first residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms */
      while (iatom < natoms && atom_ptr != NULL)
        {
/* v.4.0--> */
          /* Determine whether this atom's coords are to be included
             in the fitting */
          if (object_ptr->have_anchors == FALSE || internal_fit == TRUE)
            wanted = TRUE;

          /* If have any anchor points, then only use if this atom is
             one of them */
          else if (atom_ptr->have_anchor == TRUE)
            wanted = TRUE;
          else
            wanted = FALSE;
/* <--v.4.0 */

          /* Store the current coords and the fit coords */
/* v.4.0--> */
/*        if (nstored < MAXATOMS) */
          if (nstored < MAXATOMS && wanted == TRUE)
/* <--v.4.0 */
            {
/* v.4.0--> */
              /* If this is an anchor-atom, use its coordinates instead */
              if (atom_ptr->have_anchor == TRUE)
                {
                  atom_ptr->fit_x = atom_ptr->anchor_pstn_x;
                  atom_ptr->fit_y = atom_ptr->anchor_pstn_y;
                }
/* <--v.4.0 */

              /* Store the two sets of coordinates to be
                 fitted */
              coords1[nstored].x = atom_ptr->x;
              coords1[nstored].y = atom_ptr->y;
              coords1[nstored].z = 0.0;
              coords2[nstored].x = atom_ptr->fit_x;
              coords2[nstored].y = atom_ptr->fit_y;
              coords2[nstored].z = 0.0;

              /* Increment count of stored coords */
              nstored++;

              /* Accumulate mean values */
              mean_x = mean_x + atom_ptr->x;
              mean_y = mean_y + atom_ptr->y;
              mean_fit_x = mean_fit_x + atom_ptr->fit_x;
              mean_fit_y = mean_fit_y + atom_ptr->fit_y;
            }

/* v.4.0--> */
          /* Accumulate residue means */
          residue_mean_x = residue_mean_x + atom_ptr->x;
          residue_mean_y = residue_mean_y + atom_ptr->y;
/* <--v.4.0 */

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

/* v.4.0--> */
      /* If this is an anchor residue, calculate residue mean and
         use this as a fitting point */
      if (residue_ptr->have_anchor == TRUE && natoms > 0 &&
          internal_fit == FALSE)
        {
          /* Get means */
          residue_mean_x = residue_mean_x / natoms;
          residue_mean_y = residue_mean_y / natoms;

          /* Store residue mean position and anchor position */
          coords1[nstored].x = residue_mean_x;
          coords1[nstored].y = residue_mean_y;
          coords1[nstored].z = 0.0;
          coords2[nstored].x = residue_ptr->anchor_pstn_x;
          coords2[nstored].y = residue_ptr->anchor_pstn_y;
          coords2[nstored].z = 0.0;

          /* Increment count of stored coords */
          nstored++;

          /* Accumulate mean values */
          mean_x = mean_x + residue_mean_x;
          mean_y = mean_y + residue_mean_y;
          mean_fit_x = mean_fit_x + residue_ptr->anchor_pstn_x;
          mean_fit_y = mean_fit_y + residue_ptr->anchor_pstn_y;
        }
/* <--v.4.0 */

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* Check that have some atoms belonging to this object */
  if (nstored > 0)
    {
      /* Calculate the mean values of the two sets of coords */
      mean_x = mean_x / nstored;
      mean_y = mean_y / nstored;
      mean_fit_x = mean_fit_x / nstored;
      mean_fit_y = mean_fit_y / nstored;

      /* Loop through all the stored coordinates to centre them at
         the origin */
      for (istored = 0; istored < nstored; istored++)
        {
          /* Adjust both sets of coords */
          coords1[istored].x = coords1[istored].x - mean_x;
          coords1[istored].y = coords1[istored].y - mean_y;
          coords2[istored].x = coords2[istored].x - mean_fit_x;
          coords2[istored].y = coords2[istored].y - mean_fit_y;

          /* Store the flipped version of the object also,
             to see if this gives a better fit */
          coords1_flip[istored].x = - coords1[istored].x;
          coords1_flip[istored].y = coords1[istored].y;
          coords1_flip[istored].z = 0.0;
        }

      /* Calculate transformation that gives the best fit between
         the flattened object and the centred original structure */
      if (nstored > 1)
        {
          /* Set flag to return column-wise rotation matrix */
          column = TRUE;

          /* Try fitting the flipped version of the object first */
          flipped_rmsd = matfit(coords1_flip,coords2,flip_matrix,
                                nstored,NULL,column);

          /* Calculate agreement in terms of all atom-atom distances
             within the object */
          flipped_distance_match = get_object_distance_match(object_ptr);
          flipped_rmsd = flipped_rmsd + flipped_distance_match;

          /* Now fit object straight */
          rmsd = matfit(coords1,coords2,matrix,nstored,NULL,column);

          /* Repeat calculation of atom-atom distance matches within
             object */
          distance_match = get_object_distance_match(object_ptr);
          rmsd = rmsd + distance_match;

          /* Determine which gives the best fit */
          if (flipped_rmsd < rmsd)
            {
              best_flipped = TRUE;
              rmsd = flipped_rmsd;
            }

          /* Get required rotation matrix */
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
              if (best_flipped == TRUE)
                transform[i][j] = (float) (flip_matrix[i][j]);
              else
                transform[i][j] = (float) (matrix[i][j]);
        }

      /* If only a single atom, then make matrix the unit matrix */
      else
        {
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
              if (i == j)
                transform[i][j] = (float) (1.0);
              else
                transform[i][j] = 0.0;
        }

      /* Move object to its desired position, if required */
      if (no_move == FALSE)
        {
          /* Translate object to centre of mass */
          move_object(object_ptr,mean_x,mean_y,0.0);

          /* If object is best flipped, then do so */
          if (best_flipped == TRUE)
            flip_object(object_ptr,-1,1);

          /* Apply the rotation matrix to coords */
          rotate_object(object_ptr,transform);

          /* Move object back to where the fit coords came from */
          move_object(object_ptr,- mean_fit_x,- mean_fit_y,0.0);
        }
    }

  /* Return the calculated rms deviation */
  return(rmsd);
}
/***********************************************************************

calc_internal_energy  -  Calculate the object's total energy from any
                         close contacts between its atoms

***********************************************************************/

/* v.5.2.1 --> */
// float calc_internal_energy(struct object *object_ptr,int *clash)
float calc_internal_energy(struct object *object_ptr,int *clash,
                           int with_fit,int with_abort)
/* <-- v.5.2.1 */
{
  int conatoms, iatom, in_range, iresid, jatom, jresid, katom, natoms,
  nresid, other_natoms, wanted;
/* v.5.2.1 --> */
  int abort;
/* <-- v.5.2.1 */

  float atom_size1, atom_size2, distance, dist2, interact_dist;
  float clash_energy, energy, fit_energy, rmsd, total_energy;
  float x1, x2, y1, y2;

  struct atom_link *atom_link_ptr;
  struct coordinate *atom_ptr, *other_atom_ptr;
  struct residue *residue_ptr, *other_residue_ptr;

  /* Initialise variables */
/* v.5.2.1 --> */
  abort = FALSE;
  rmsd = 0;
/* <-- v.5.2.1 */
  *clash = FALSE;
  clash_energy = 0.0;
  fit_energy = 0.0;
  total_energy = 0.0;

  /* Calculate the maximum and minimum coords of atoms in each
     residue to speed up computation */
  update_residue_boundaries(object_ptr);

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
/* v.5.2.1 --> */
//  while (iresid < nresid && residue_ptr != NULL)
  while (iresid < nresid && residue_ptr != NULL && abort == FALSE)
/* <-- v.5.2.1 */
    {
      /* Get a second pointer to this residue for calculating
         inter-residue contacts */
      other_residue_ptr = residue_ptr;
      jresid = iresid;

      /* Loop over all the possible other residues */
/* v.5.2.1 --> */
//      while (jresid < nresid && other_residue_ptr != NULL)
      while (jresid < nresid && other_residue_ptr != NULL &&
             abort == FALSE)
/* <-- v.5.2.1 */
        {
          /* Get the maximum atom-size of each residue */
          atom_size1 = residue_ptr->max_atom_size;
          atom_size2 = other_residue_ptr->max_atom_size;
          interact_dist = atom_size1 + atom_size2 + INTERNAL_INTERACT_DIST;

          /* Check whether the residues are within range of one another */
          in_range = TRUE;
          if ((other_residue_ptr->minx >
               (residue_ptr->maxx + interact_dist)) ||
              (residue_ptr->minx >
               (other_residue_ptr->maxx + interact_dist)) ||
              (other_residue_ptr->miny >
               (residue_ptr->maxy + interact_dist)) ||
              (residue_ptr->miny >
               (other_residue_ptr->maxy + interact_dist)))
            in_range = FALSE;

          /* If residues are within range, then check all distances
             between them */
          if (in_range == TRUE)
            {
              /* Get pointer to first residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;

              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the first residue's atoms */
/* v.5.2.1 --> */
//              while (iatom < natoms && atom_ptr != NULL)
              while (iatom < natoms && atom_ptr != NULL && abort == FALSE)
/* <-- v.5.2.1 */
                {
                  /* Get atom's coordinates */
                  x1 = atom_ptr->x;
                  y1 = atom_ptr->y;

                  /* Get the atom's size */
                  atom_size1 = atom_ptr->atom_size;

                  /* Get pointer to other residue's first atom */
                  other_atom_ptr = other_residue_ptr->first_atom_ptr;
                  
                  /* Initialise atom count */
                  jatom = 0;
                  other_natoms = other_residue_ptr->natoms;
                  
                  /* Loop over all the other residue's atoms */
/* v.5.2.1 --> */
//                  while (jatom < other_natoms && other_atom_ptr != NULL)
                  while (jatom < other_natoms && other_atom_ptr != NULL &&
                         abort == FALSE)
/* <-- v.5.2.1 */
                    {
                      /* Initialise flag */
                      wanted = TRUE;

                      /* If both are the same atom, then don't want them */
                      if (other_atom_ptr == atom_ptr)
                        wanted = FALSE;
                      
                      /* Otherwise, check that they are not covalently
                         bonded */
                      else
                        {
                          /* Get pointer to the first atom covalently
                             bonded to current atom of interest */
                          atom_link_ptr = atom_ptr->first_atom_link_ptr;
                          katom = 0;
                          conatoms = atom_ptr->natom_links;

                          /* Loop over all the connected atoms */
                          while (katom < conatoms && atom_link_ptr != NULL
                                 && wanted == TRUE)
                            {
                              /* Check whether the connected atom is
                                 the same as the other atom currently
                                 under consideration */
                              if (atom_link_ptr->atom_ptr ==
                                  other_atom_ptr)
                                wanted = FALSE;
        
                              /* Get pointer to next connected atom */
                              atom_link_ptr
                                = atom_link_ptr->next_atom_link_ptr;
                              katom++;
                            }
                        }

                      /* If atoms are not bonded, calculate the distance
                         between them */
                        if (wanted == TRUE)
                        {
                          /* Get other atom's coordinates */
                          x2 = other_atom_ptr->x;
                          y2 = other_atom_ptr->y;

                          /* Get the other atom's size */
                          atom_size2 = other_atom_ptr->atom_size;

                          /* Calculate the distance between the
                             two atoms */
                          dist2 = (x1 - x2) * (x1 - x2)
                            + (y1 - y2) * (y1 - y2);
                          distance = (float) sqrt((double) dist2);

                          /* Adjust distance to allow for the atom sizes */
                          if (distance > (atom_size1 + atom_size2))
                            distance = distance - (atom_size1 + atom_size2);
                          else
                            distance = 0.0;
                          dist2 = distance * distance;

                          /* Add energy of this interaction to the
                             overall energy sum */
                          if (dist2 < 0.001)
                            dist2 = (float) (0.001);
                          if (dist2 < INTERNAL_DIST2)
                            energy = (float) (Atom_Atom_Clash / dist2);
                          else
                            energy = 0.0;
                          clash_energy = clash_energy + energy;

/* v.5.2.1 --> */
                          /* If energy is already above minimum allowed
                             and we need to bail out, then set flag */
                          if (with_abort == TRUE &&
                              clash_energy > MAX_INTERNAL_ENERGY)
                            abort = TRUE;
/* <-- v.5.2.1 */
                        }

                      /* Get pointer to other residue's next atom */
                      other_atom_ptr = other_atom_ptr->next;
                      jatom++;
                    }

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }
            }

          /* Get pointer to the next residue */
          other_residue_ptr = other_residue_ptr->next_residue_ptr;
          jresid++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }

  /* See how well the current object can be fitted onto the original */
/* v.5.2.1 --> */
  if (ligfit == FALSE && with_fit == TRUE && abort == FALSE)
/* <-- v.5.2.1 */
    rmsd = fit_object(object_ptr,TRUE);

  /* Calculate an equivalent energy according to the rms deviation of the
     best fit */
  fit_energy = rmsd;
  
  /* Get total object energy from its components */
  total_energy = clash_energy + fit_energy;
  if (fabs(clash_energy) < 0.001)
    *clash = FALSE;
  else
    *clash = TRUE;

  /* Scale the energy up by its weight, as read in from the run parameters */
  total_energy = Internal_Energy_Weight * total_energy;

/* v.5.2.1 --> */
  /* If calculation aborted, make sure energy is above maximum */
  if (abort == TRUE)
    total_energy = 2 * MAX_INTERNAL_ENERGY;
/* <-- v.5.2.1 */

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

flip_about_bond  -  Flip part of the current object about the given bond,
                    transforming all the appropriate atoms' coordinates
                    from +y to -y

***********************************************************************/

void flip_about_bond(struct bond *flip_bond_ptr,int end)
{
  int loop, nbonds;

  struct bond *bond_ptr, *next_free_stack_ptr, *other_bond_ptr,
  *next_stack_ptr;
  struct bond_link *bond_link_ptr;
  struct coordinate *atom_ptr;

  /* Initialise variables */
  next_stack_ptr = NULL;
  next_free_stack_ptr = NULL;

  /* Initialise all the bond flags and stack pointers */
  initialise_bonds();

  /* Initialise all the atom flags */
  initialise_atoms();

  /* If haven't been given which end to use for the flip, then use
     whichever has fewer downstream of it */
  if (end != FIRST && end != SECOND)
    {
      /* If first atom has fewer bonds downstream of it, then want to
         process the bonds coming off it */
      if (flip_bond_ptr->nbonds_from_first_atom
          < flip_bond_ptr->nbonds_from_second_atom)
        end = FIRST;
  
      /* Otherwise, use the other end of the bond */
      else
        end = SECOND;
    }

  /* Set pointer to the first of the bonds sprouting off the flip
     bond */
  bond_link_ptr = flip_bond_ptr->first_bond_link_ptr;

  /* Set the flip-bond as already checked */
  flip_bond_ptr->checked = TRUE;

  /* Initialise count of bonds encountered */
  nbonds = 0;

  /* Initialise stack pointers for bond-search by putting the bonds
     coming off the selected end of the flip-bond onto the stack
     of bonds to be searched */
  while (bond_link_ptr != NULL)
    {
      /* Add this bond to the stack only if it comes off the relevant
         atom of the test-bond */
      if (bond_link_ptr->bond_end == end)
        {
          /* Get the bond this link is pointing to and add to stack */
          bond_ptr = bond_link_ptr->bond_ptr;

          /* Add to stack provided it is not an elastic bond */
          if (bond_ptr->elastic == FALSE)
            {
              /* If this is the first to be added to the stack, set both
                 stack-pointers to point to it */
              if (next_free_stack_ptr == NULL)
                {
                  next_free_stack_ptr = bond_ptr;
                  next_stack_ptr = bond_ptr;
                }

              /* Otherwise, make last bond on stack point to this one
                 and set this one as nopw the last on the stack */
              else
                {
                  next_free_stack_ptr->next_stack_ptr = bond_ptr;
                  next_free_stack_ptr = bond_ptr;
                }
            }

          /* Mark bond as added to the stack */
          bond_ptr->checked = TRUE;
        }
      
      /* Go to the next bond-link in the linked list */
      bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
    }

  /* Loop until stack of pointers to be processed has been exhausted */
  while (next_stack_ptr != NULL)
    {
      /* Pick up the next bond from the stack, and move the stack
       pointer onto the next position */
      bond_ptr = next_stack_ptr;
      next_stack_ptr = bond_ptr->next_stack_ptr;

      /* Increment count of bonds encountered */
      nbonds++;

      /* Loop over the bond's two atoms */
      for (loop = 0; loop < 2; loop++)
        {
          /* Get pointer to appropriate atom */
          if (loop == 0)
            atom_ptr = bond_ptr->first_atom_ptr;
          else
            atom_ptr = bond_ptr->second_atom_ptr;

          /* If first atom has not been flipped, then flip it */
          if (atom_ptr->checked == FALSE)
            {
              atom_ptr->y = - atom_ptr->y;
              atom_ptr->checked = TRUE;
            }
        }

     /* Loop through all the bonds connected to the current one
        to add to stack */

      /* Get the pointer to the first of the bonds coming off the 
         current bond */
      bond_link_ptr = bond_ptr->first_bond_link_ptr;

      /* Process each of these bonds in turn */
      while (bond_link_ptr != NULL)
        {
          /* Get the pointer to this bond */
          other_bond_ptr = bond_link_ptr->bond_ptr;

          /* If this bond hasn't already been processed, then move
             its atoms (if they need moving) and add the bond to the
             stack */
          if (other_bond_ptr->checked == FALSE &&
              other_bond_ptr->elastic == FALSE)
            {
              next_free_stack_ptr->next_stack_ptr = other_bond_ptr;
              next_free_stack_ptr = other_bond_ptr;

              /* If the stack has run out, restart it at the bond
                 just entered */
              if (next_stack_ptr == NULL)
                next_stack_ptr = other_bond_ptr;

              /* Mark bond as added to the stack */
              other_bond_ptr->checked = TRUE;
            }

          /* Go to the next bond-link in the linked list */
          bond_link_ptr = bond_link_ptr->next_bond_link_ptr;
        }
    }
    
  /* Update the counter of bonds downstream on flipped bond */
  if (end == FIRST)
    flip_bond_ptr->nbonds_from_first_atom = nbonds;
  else
    flip_bond_ptr->nbonds_from_second_atom = nbonds;
}
/***********************************************************************

flip_back  -  Flip all the flagged atoms back the way they were

***********************************************************************/

void flip_back(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* If the atom is flagged as flipped, then flip it back */
      if (atom_ptr->checked == TRUE)
        atom_ptr->y = - atom_ptr->y;

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }          
}
/***********************************************************************

flip_attached_groups  -  Flip all the groups attached by H-bonds and
                         non-bonded contacts to the ligand atoms that
                         have just been flipped

***********************************************************************/

void flip_attached_groups(struct object *flipped_object_ptr,
                          float move_x,float move_y,float matrix[4][3])
{
  int i, j;
  int wanted;
  float inverse_matrix[4][3];

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Check if this is an elastic bond */
      if (bond_ptr->elastic == TRUE)
        {
          /* Get the pointers to the two atoms at either end of the bond */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Check whether one of the atoms has been flipped and the
             other not */
          wanted = TRUE;
          if ((atom1_ptr->checked == TRUE && atom2_ptr->checked == TRUE) ||
              (atom1_ptr->checked == FALSE && atom2_ptr->checked == FALSE))
            wanted = FALSE;

          /* Check that one atom belongs to the flipped object and
             the other one doesn't */
          else
            {
              /* Get the objects these two atoms belong to */
              object1_ptr = atom1_ptr->residue_ptr->object_ptr;
              object2_ptr = atom2_ptr->residue_ptr->object_ptr;

              /* Check that one is the flipped and the other isn't */
              if (object1_ptr == flipped_object_ptr &&
                  object2_ptr != flipped_object_ptr)
                {
                  /* Confirm that the atom on the flipped object is the
                     checked one */
                  if (atom1_ptr->checked == TRUE)
                    {
                      /* Save the pointer to the other object */
                      object_ptr = object2_ptr;
                    }
                  else
                    wanted = FALSE;
                }               
              else if (object2_ptr == flipped_object_ptr &&
                       object1_ptr != flipped_object_ptr)
                {
                  /* Confirm that the atom on the flipped object is the
                     checked one */
                  if (atom2_ptr->checked == TRUE)
                    {
                      /* Save the pointer to the other object */
                      object_ptr = object1_ptr;
                    }
                  else
                    wanted = FALSE;
                }               

              /* Otherwise, don't want to flip anything */
              else
                wanted = FALSE;
            }

          /* If the bond is a contact bond, then only want it if the
             group involved is a hydrophonic group */
          if (wanted == TRUE && bond_ptr->bond_type == CONTACT &&
              object_ptr->object_type != HYDROPHOBIC)
            wanted = FALSE;

          /* If have a valid bond, perform the flip of this object */
          if (wanted == TRUE)
            {
              /* Transform the current object into the reference
                 frame of the flipped object */
              move_object(object_ptr,move_x,move_y,0.0);
              transform_object(object_ptr,matrix);

              /* Perform the flip */
              flip_object(object_ptr,1,-1);
              
              /* Generate the inverse rotation matrix to that used
                 above and perform the inverse rotation */
              for (i = 0; i < 4; i++)
                for (j = 0; j < 3; j++)
                  inverse_matrix[i][j] = matrix[i][j];
              inverse_matrix[0][1] = - matrix[0][1];
              inverse_matrix[1][0] = - matrix[1][0];
              transform_object(object_ptr,inverse_matrix);

              /* Move object back to where it came from */
              move_object(object_ptr,- move_x,- move_y,0.0);
            }
        }

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************

bond_flip  -  Apply the given flip(s) about the given bond

***********************************************************************/

void bond_flip(struct object *object_ptr,struct bond *bond_ptr,
/* v.5.2.1 --> */
//               int use_end,int both_ends,int flip_all)
               int use_end,int both_ends,int flip_all,int calc_en)
/* <-- v.5.2.1 */
{
  int clash;
/* v.5.2.1 --> */
  int with_abort, with_fit;
/* <-- v.5.2.1 */

  float coord_store[MAXATOMS][3], matrix[4][3], move_x, move_y;
  float internal_energy;

  struct coordinate *atom_ptr;

  /* Get the coordinates of the atom at the required end */
  if (use_end == FIRST)
    {
      move_x = bond_ptr->first_atom_ptr->x;
      move_y = bond_ptr->first_atom_ptr->y;
      atom_ptr = bond_ptr->second_atom_ptr;
    }
  else
    {
      move_x = bond_ptr->second_atom_ptr->x;
      move_y = bond_ptr->second_atom_ptr->y;
      atom_ptr = bond_ptr->first_atom_ptr;
    }

  /* Move the current object to place this atom at the origin */
  move_object(object_ptr,move_x,move_y,0.0);

  /* Store the coordinates of the bond two ends */
  coord_store[0][0] = 0.0;
  coord_store[0][1] = 0.0;
  coord_store[0][2] = 0.0;
  coord_store[1][0] = atom_ptr->x;
  coord_store[1][1] = atom_ptr->y;
  coord_store[1][2] = 0.0;

  /* Calculate the transformation that will put this bond along the
     x-axis, with atom 1 at the origin and atom 2 pointing in the
     -ve x-direction */
  line_transformation(coord_store,matrix,0.0);

  /* Apply the transformation to all atoms in the current object */
  transform_object(object_ptr,matrix);

  /* Flip one side of the object around this rotatable bond
     (currently lying along the x-axis) by changing the signs of all
     the y-coords */
  flip_about_bond(bond_ptr,use_end);

  /* If flipping both ends, then flip the other end too */
  if (both_ends == TRUE)
    flip_about_bond(bond_ptr,SECOND);

  /* Calculate energy of the flipped configuration */
/* v.5.2.1 --> */
//  internal_energy = calc_internal_energy(object_ptr,&clash);
  if (calc_en == TRUE)
    {
      with_fit = TRUE;
      with_abort = FALSE;
      internal_energy
        = calc_internal_energy(object_ptr,&clash,with_fit,with_abort);
/* <-- v.5.2.1 */
      object_ptr->internal_energy = internal_energy;
/* v.5.2.1 --> */
    }
/* <-- v.5.2.1 */

  /* If are flipping all the attached groups, as well,
     then do those */
  if (flip_all == TRUE)
    {
      bond_ptr->first_atom_ptr->checked = FALSE;
      bond_ptr->second_atom_ptr->checked = FALSE;
      flip_attached_groups(object_ptr,move_x,move_y,matrix);
    }

  /* Generate the inverse rotation matrix to that used
     above and perform the inverse rotation */
  matrix[0][1] = - matrix[0][1];
  matrix[1][0] = - matrix[1][0];
  transform_object(object_ptr,matrix);

  /* Move object back to where it came from */
  move_object(object_ptr,- move_x,- move_y,0.0);
}
/***********************************************************************

get_random_number  -  Routine for generating a random number from the
                      supplied seed. Returns a uniform random deviate
                      between 0.0 and 1.0. On the first call, the
                      seed must be set to any negative value to
                      initialise or reinitialise the sequence.
                      (Algorithm is that described for routine ran2
                      in Numerical Recipes).

***********************************************************************/

float get_random_number(int *iseed,int random_start)
{
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */

  static int ir[97], iy, j;

  static float random_number, seed;

  static int first_pass = TRUE;
  static int m = 714025;
  static int ia = 1366;
  static int ic = 150889;
  static float rm = (float) (1.4005112E-6);

  FILE *file_ptr;

  /* Initialise the random number */
  random_number = 0.0;

  /* If this is the first pass, or sequence to be re-initialised, then
     pick up seed from file */
  if (*iseed < 0 || first_pass == TRUE)
    {
      /* Initialise the seed */
      seed = 0.0;

      /* If a completely random start is required, then read in the
         number written out to the file ranseed.dat last time this
         routine was run */
      if (random_start == TRUE)
        {
          /* Open the random seed file and read in seed */
/* v.5.0--> */
          /* if ((file_ptr = fopen("ranseed.dat","r")) !=  NULL) */
          filename[0] = '\0';
          if (wkdir[0] != '\0')
            strcpy(filename,wkdir);
          strcat(filename,"ranseed.dat");
          if ((file_ptr = fopen(filename,"r")) == NULL)
/* <--v.5.0 */
            {
              /* If file exists, then read in the seed */
              fscanf(file_ptr,"%f",&seed);

              /* Close the file */
              fclose(file_ptr);
            }

          /* Store integer version of the seed */
          *iseed = (int) seed;
        }

      /* Otherwise, use a standard starting point */
      else
        *iseed = 54813;

      /* Check that the seed is non-zero */
      if (*iseed == 0) 
        *iseed = 1;

      /* Generate random number from seed */
      *iseed = (ic - *iseed) % m;
      for (j = 0; j < 97; j++)
        {
          *iseed = (ia * (*iseed) + ic) % m;
          ir[j] = *iseed;
        }
      *iseed = (ia * (*iseed) + ic) % m;
      iy = *iseed;
    }

  /* Check for error */
  j = (97 * iy) / m;
  if (j >= 97 || j < 0)
    {
      printf("*** Random number error\n");
      exit(1);
    }

  /* Get the random number */
  iy = ir[j];
  random_number = iy * rm;
  *iseed = (ia * (*iseed) + ic) % m;
  ir[j] = *iseed;

  /* If this is the first call, use the generated random number to
     create a new seed and write to disk */
  if (first_pass == TRUE)
    {
      first_pass = FALSE;
      seed = iy * rm * 100000;

      /* Open the random seed file */
/* v.5.0--> */
      /* if ((file_ptr = fopen("ranseed.dat","w")) !=  NULL) */
      filename[0] = '\0';
/* v.5.3.8--> */
      // if (wkdir[0] != '\0')
      if (wkdir[0] != '\0' && highlight_res[0] == '\0')
/* <--v.5.3.8 */
        strcpy(filename,wkdir);
      strcat(filename,"ranseed.dat");
      if ((file_ptr = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
        {
          /* If no file error, then write out the seed */
          fprintf(file_ptr,"%f\n",seed);

          /* Close the file */
          fclose(file_ptr);
        }
    }

  /* Return the random number generated */
  return(random_number);
}
/***********************************************************************

save_restore_coords  -  Save the given object's atom positions

***********************************************************************/

void save_restore_coords(struct object *object_ptr,int action)
{
  int iatom, iresid, natoms, nresid;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Set the pointer to the first of this object's residues */
  residue_ptr = object_ptr->first_residue_ptr;
  nresid = object_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid && residue_ptr != NULL)
    {
      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;

      /* Initialise atom count */
      iatom = 0;
      natoms = residue_ptr->natoms;

      /* Loop over all the residue's atoms, performing the save */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Save the x- and y-coords */
          if (action == SAVE)
            {
              atom_ptr->save_x = atom_ptr->x;
              atom_ptr->save_y = atom_ptr->y;
            }

          /* Restore the coords */
          else
            {
              atom_ptr->x = atom_ptr->save_x;
              atom_ptr->y = atom_ptr->save_y;
            }

          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }

      /* Get pointer to the next residue */
      residue_ptr = residue_ptr->next_residue_ptr;
      iresid++;
    }
}
/***********************************************************************

save_restore_all_objects  -  Loop through all objects to save/restore
                             their coordinates

***********************************************************************/

void save_restore_all_objects(int restore)
{
  int iobject;

  struct object *object_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Save/restore this object's coords */
      save_restore_coords(object_ptr,restore);

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }
}
/***********************************************************************

calc_rotation_matrix  -  Calculate the rotation matrix for a given
                         rotation about the z-axis

***********************************************************************/

void calc_rotation_matrix(float angle,float matrix[3][3])
{
  int i, j;

  float cos_theta, sin_theta;

  /* Initialise the matrix */
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
        {
          if (i == j)
            matrix[i][j] = (float) (1.0);
          else
            matrix[i][j] = 0.0;
        }
    }

  /* Calculate sines and cosines of the angle */
  cos_theta = (float) cos(- angle / RADDEG);
  sin_theta = (float) sin(- angle / RADDEG);

  /* Form the appropriate rotation matrix */
  matrix[0][0] = cos_theta;
  matrix[0][1] = sin_theta;
  matrix[1][0] = - sin_theta;
  matrix[1][1] = cos_theta;
}
/***********************************************************************

get_matrix  -  Calculate the rotation matrix for a rotation about the
               z-axis through the given angle

***********************************************************************/

void get_matrix(float angle,float matrix[4][3])
{
  int i, j;

  float cos_theta, sin_theta;

  /* Initialise the matrix */
  for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 3; j++)
        {
          if (i == j && i < 3)
            matrix[i][j] = (float) (1.0);
          else
            matrix[i][j] = 0.0;
        }
    }

  /* Calculate sines and cosines of the angle */
  cos_theta = (float) cos(angle);
  sin_theta = (float) sin(angle);

  /* Form the appropriate rotation matrix */
  matrix[0][0] = cos_theta;
  matrix[0][1] = sin_theta;
  matrix[1][0] = - sin_theta;
  matrix[1][1] = cos_theta;
}
/***********************************************************************

untangle_object  -  Untangle atom and bond overlaps by flipping structure
                    about its rotatable bonds

***********************************************************************/

void untangle_object(struct object *object_ptr,int *rseed)
{
  int end, iseed;
  int cycle, nflips;
  int clash, last_clash, swivel;
/* v.5.2.1 --> */
  int with_abort, with_fit;
/* <-- v.5.2.1 */

  float coord_store[3][3], matrix[4][3];
  float energy, last_energy;
  float rotate_angle, x;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr, *swap_ptr;
  struct object_bond *object_bond_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Initialise */
  nflips = 1;
/* v.5.2.1 --> */
  with_fit = TRUE;
  with_abort = FALSE;
/* <-- v.5.2.1 */

  /* Perform the energy minimization process */
  for (cycle = 0; cycle < MAXCYCLE && nflips != 0; cycle++)
    {
      /* Initialise count of flips made to structure */
      nflips = 0;
      if (cycle > 0 && cycle % 10 == 0)
        {
          if (object_ptr->object_type == LIGAND)
            printf("Untangling ligand,           cycle %4d\n",cycle);
          else
            printf("   Untangling residue %s %s %c, cycle %4d\n",
                   object_ptr->first_residue_ptr->res_name,
                   object_ptr->first_residue_ptr->res_num,
                   object_ptr->first_residue_ptr->chain,cycle);
        }

      /* Calculate current energy of the object from atom clashes */
/* v.5.2.1 --> */
//      energy = calc_internal_energy(object_ptr,&clash);
      energy
        = calc_internal_energy(object_ptr,&clash,with_fit,with_abort);
/* <-- v.5.2.1 */
      last_clash = clash;
      last_energy = energy;

      /* Get pointer to the first of this object's bonds */
      object_bond_ptr = object_ptr->first_object_bond_ptr;
      
      /* Loop through all object's bonds */
      while (object_bond_ptr != NULL)
        {
          /* Get the bond pointer */
          bond_ptr = object_bond_ptr->bond_ptr;

          /* If this is a rotatable bond see if object will look
             better if flipped about this bond */
          if (bond_ptr->rotatable_bond == TRUE)
            {
              /* Initialise swivel flag */
              swivel = FALSE;

              /* Get pointers to the bond's two atoms */
              atom1_ptr = bond_ptr->first_atom_ptr;
              atom2_ptr = bond_ptr->second_atom_ptr;

              /* If have been flipping a long time, and this is a
                 side-chain bond, may want to try sivelling the
                 sidechain about this bond */
              if (cycle + 1 >= START_SWIVELS &&
                  (atom1_ptr->side_chain == TRUE ||
                  atom2_ptr->side_chain == TRUE))
                {
                  /* Determine whether to flip or swivel */
                  x = get_random_number(&iseed,Random_Start);
                  if (x <= 0.5)
                    swivel = TRUE;
                }

              /* If swivelling, then store all the object's coordinates
                 ane determine about which end of the bond to perform the
                 swivel */
              if (swivel == TRUE)
                {
                  /* Store all the object's atoms */
                  save_restore_coords(object_ptr,SAVE);

                  /* Determine which end to use for the swivel */
                  x = get_random_number(&iseed,Random_Start);
                  if (x <= 0.5)
                    end = FIRST;
                  else
                    {
                      /* To swivel about the second end, swap the
                         bond's two atoms around */
                      end = SECOND;
                      swap_ptr = atom1_ptr;
                      atom1_ptr = atom2_ptr;
                      atom2_ptr = swap_ptr;
                    }
                }

              /* Store the two atoms' coordinates */
              coord_store[0][0] = atom1_ptr->x;
              coord_store[0][1] = atom1_ptr->y;
              coord_store[0][2] = atom1_ptr->z;
              coord_store[1][0] = atom2_ptr->x;
              coord_store[1][1] = atom2_ptr->y;
              coord_store[1][2] = atom2_ptr->z;

              /* Calculate the transformation that will put this bond
                 along the x-axis, with atom 1 at the origin and atom 2
                 pointing in the -ve x-direction */
              line_transformation(coord_store,matrix,0.0);

              /* Apply the transformation to all atoms in the current
                 object */
              transform_object(object_ptr,matrix);

              /* If swivelling about the bonds, then do so */
              if (swivel == TRUE)
                {
                  /* Get swivel angle */
                  x = get_random_number(&iseed,Random_Start);
                  if (x <= 0.5)
                    rotate_angle = SWIVEL_ANGLE / RADDEG;
                  else
                    rotate_angle = - SWIVEL_ANGLE / RADDEG;

                  /* Form the transformation matrix to perform the
                     swivel */
                  get_matrix(rotate_angle,matrix);

                  /* Apply the transformation to all atoms springing off
                     the current bond */
                  transform_downstream_atoms(bond_ptr,end,matrix);
                }

              /* Otherwise, flip one side of the object around this
                 rotatable bond (currently lying along the x-axis) by
                 changing the signs of all the y-coords */
              else
                flip_about_bond(bond_ptr,EITHER);

              /* Calculate energy of the flipped configuration */
/* v.5.2.1 --> */
//              energy = calc_internal_energy(object_ptr,&clash);
              energy
                = calc_internal_energy(object_ptr,&clash,with_fit,
                                       with_abort);
/* <-- v.5.2.1 */

              /* If the energy is greater than that of the previous
                 configuration, then restore to previous conformation */
              if (energy > (last_energy - 0.000001))
                {
                  /* If performed a swivel, then restore all the
                     coordinates as they were before */
                  if (swivel == TRUE)
                    save_restore_coords(object_ptr,RESTORE);

                  /* Otherwise, just flip the y-values back */
                  else
                    flip_back();

                  /* Retrieve the previous energy */
                  energy = last_energy;
                  clash = last_clash;
                }

              /* Otherwise, increment count of flips made */
              else
                nflips++;

              /* Save the current energy */
              last_energy = energy;
              object_ptr->internal_energy = energy;
            }

          /* Get pointer to this object's next bond entry */
          object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
        }

      /* Keep going until have started swivelling */
      if (clash == TRUE && nflips == 0)
        nflips = 1;
    }

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

flatten_all_objects  -  Flatten all objects into 2D by squashing rings
                        and unravelling along all rotatable bonds

***********************************************************************/

/* v.4.1.2--> */
/* void flatten_all_objects(int *iseed) */
void flatten_all_objects(int *iseed,int debug)
/* <--v.4.1.2 */
{
  int iobject, n_flattened;
  int rseed;
/* v.5.0--> */
  int fixed;
/* <--v.5.0 */

  struct object *object_ptr;

  /* Initialise pointer to the first stored object */
  printf("\nFlattening and untangling objects ...\n");
  object_ptr = first_object_ptr;
  iobject = 0;
  rseed = *iseed;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If plotting schematic diagram, prepare the ligand or H-group */
      if (Include->Simple_Nonligand_Residues == TRUE)
        {
          /* If this is an H-group, remove all unwanted atoms */
          if (object_ptr->object_type == HGROUP)
            create_simple_hgroup(object_ptr);

          /* Delete any unwanted atoms */
/* v.4.1.2--> */
/*        delete_atoms_and_bonds(); */
          delete_atoms_and_bonds(debug,"flatten_all_objects ... 1");
/* <--v.4.1.2 */
        }

      /* Otherwise, if not plotting non-ligand mainchain atoms, then
         get rid of mainchains not involved in interactions */
      else if (Include->Mainchain_Atoms == FALSE &&
               object_ptr->object_type == HGROUP)
        {
          /* Remove unwanted mainchain atoms */
          remove_mainchains(object_ptr);

          /* Delete any unwanted atoms */
/* v.4.1.2--> */
/*        delete_atoms_and_bonds(); */
          delete_atoms_and_bonds(debug,"flatten_all_objects ... 2");
/* <--v.4.1.2 */
        }

      /* Update each residue's maximum and minimum coords */
      update_residue_boundaries(object_ptr);

      /* If object is just a hydrophobic group, (or a H-bonded group
         in the simplified representation) then place at origin */
      if (object_ptr->object_type == HYDROPHOBIC ||
          object_ptr->object_type == SIMPLE_HGROUP ||
          object_ptr->object_type == WATER)
/* v.5.3.4--> */
        {
/* <--v.5.3.4 */
          place_group_at_origin(object_ptr);
/* v.5.3.4--> */
          fixed = fix_object(object_ptr);
        }
/* <--v.5.3.4 */

      /* Otherwise, perform the full flattening of rings, unrolling
         and untangling procedure */
      else
        {
          /* Flatten any rings in the structure so that their non-planarity
             doesn't mess up the process of unrolling the structure */
          n_flattened = flatten_all_rings(object_ptr);

          /* If any rings have been flattened, then may need to correct
             the lengths of bonds attached to them as these may have
             become deformed */
          if (n_flattened > 0)
            correct_bond_lengths(object_ptr);

          /* Fix any bonds springing off the rings */
          flatten_ring_offshoots(object_ptr);

          /* Unroll the structure by flattening out the bonds either
             side of each rotatable bond */
          unroll_structure(object_ptr);

          /* Prior to forcing flat, align object as best as possible
             with x-y plane */
          align_object(object_ptr);

          /* Force the object into x-y plane, whether or not it's quite
             there yet, so that flipping processes don't mess it up
             completely */
          force_flat(object_ptr);

          /* If plotting schematic diagram, prepare this object if
             it is a ligand */
          if (Include->Simple_Ligand_Residues == TRUE &&
              object_ptr->object_type == LIGAND)
            {
              /* Stretch the ligand residues out along the mainchain
                 and delete any unwanted sidechain atoms */
              create_simple_ligand(object_ptr);

              /* Delete any unwanted atoms */
/* v.4.1.2--> */
/*            delete_atoms_and_bonds(); */
              delete_atoms_and_bonds(debug,"flatten_all_objects ... 3");
/* <--v.4.1.2 */
            }

/* v.5.0--> */
          /* Check if all object's atoms are anchored and if we can set
             its starting coords to these anchor positions */
          fixed = fix_object(object_ptr);

          /* Flip about rotatable bonds to minimize atom and
             bond overlaps as far as is possible for this object */
          if (fixed == FALSE)
/* <--v.5.0 */
            untangle_object(object_ptr,&rseed);
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Retrieve the current random-number seed */
  *iseed = rseed;
}
/***********************************************************************
 
score_objects - Score all objects according to their importance

***********************************************************************/

void score_objects(void)
{
  int iresid, natoms, nresid;
  int weight;

  struct bond *bond_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue_ptr, *residue1_ptr, *residue2_ptr;

  /* Get pointer to the first of the objects */
  object_ptr = first_object_ptr;

  /* Loop through all the objects */
  while (object_ptr != NULL)
    {
      /* Initialise weight */
      weight = 0;

      /* Determine weight according to object type */

      /* Highest weight goes to ligand itself */
      if (object_ptr->object_type == LIGAND)
        weight = 100;

      /* Single-atom groups weight the lowest */
      else if (object_ptr->object_type == HYDROPHOBIC ||
               object_ptr->object_type == SIMPLE_HGROUP ||
               object_ptr->object_type == WATER)
        weight = 1;

      /* Otherwise, for ordinary H-groups, weight according to their
         number of atoms */
      else
        {
          /* Set the pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get the number of atoms in this residue and add to weight */
              natoms = residue_ptr->natoms;
              weight = weight + natoms;

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }
        }

      /* Store object's weight */
      object_ptr->weight = weight;

      /* Get pointer to the next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through the stored bond-list to see if already have bond */
  while (bond_ptr != NULL)
    {
      /* Initialise the weight */
      weight = 0;

      /* Get the two objects at either end of the bond */
      residue1_ptr = bond_ptr->first_atom_ptr->residue_ptr;
      residue2_ptr = bond_ptr->second_atom_ptr->residue_ptr;
      object1_ptr = residue1_ptr->object_ptr;
      object2_ptr = residue2_ptr->object_ptr;

      /* If this is a covalent bond, then only want it if it is between
         ligand and non-ligand residues */
      if (bond_ptr->bond_type == COVALENT)
        {
          /* If not the same object, then check if one is the ligand */
          if (object1_ptr != object2_ptr)
            {
              if (residue1_ptr->inligand == TRUE ||
                  residue2_ptr->inligand == TRUE)
                weight = 20;
            }
        }

      /* Check for H-bond */
/* v.5.4--> */
//      else if (bond_ptr->bond_type == HBOND)
      else if (bond_ptr->bond_type == HBOND ||
               bond_ptr->bond_type == SALT_BRIDGE ||
               bond_ptr->bond_type == DISULPHIDE)
/* <--v.5.4 */
        {
          /* Assign H-bond weight */
          weight = 15;
        }

      /* Check for non-bonded contact */
      else if (bond_ptr->bond_type == CONTACT)
        {
          /* Assign weight for non-bonded contact */
          weight = 5;
        }

      /* Add the weight for the current bond to both objects involved */
      object1_ptr->weight = object1_ptr->weight + weight;
      if (object1_ptr != object2_ptr)
        object2_ptr->weight = object2_ptr->weight + weight;

      /* Go to the next bond in the linked list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/***********************************************************************
 
sort_objects - Sort all objects according to their importance using their
               computed weights

***********************************************************************/
void sort_objects(void)
{
  int flipped;
  int weight, other_weight;

  struct object *object_ptr, *last_object_ptr, *other_object_ptr;

  /* Initialise variables */
  flipped = TRUE;

  /* Loop until all the values have been sorted */
  while (flipped == TRUE)
    {
      /* Initialise flag */
      flipped = FALSE;

      /* Get pointer to the first of the objects */
      object_ptr = first_object_ptr;
      last_object_ptr = NULL;

      /* Loop through all the objects */
      while (object_ptr != NULL)
        {
          /* Get the object's weight */
          weight = object_ptr->weight;

          /* Get pointer to the next object */
          other_object_ptr = object_ptr->next_object_ptr;

          /* If have the next object, then proceed */
          if (other_object_ptr != NULL)
            {
              /* Get the other object's weight */
              other_weight = other_object_ptr->weight;

              /* If other object's weight is greater, then swap the
                 order of the objects */
              if (other_weight > weight)
                {
                  /* Set flag to indicate that swap occurred */
                  flipped = TRUE;

                  /* If this is the first object, get the first pointer
                     to point to the second object */
                  if (last_object_ptr == NULL)
                      first_object_ptr = other_object_ptr;

                  /* Otherwise, set the last object's pointer to the
                     other object */
                  else
                    last_object_ptr->next_object_ptr = other_object_ptr;

                  /* Swap the pointers */
                  object_ptr->next_object_ptr
                    = other_object_ptr->next_object_ptr;
                  other_object_ptr->next_object_ptr = object_ptr;

                  /* Set the other object as the last object */
                  last_object_ptr = other_object_ptr;
                }

              /* Otherwise, set the current object as the last object */
              else
                last_object_ptr = object_ptr;
            }
          else
            last_object_ptr = object_ptr;

          /* Get pointer to the next object in the stack */
          object_ptr = last_object_ptr->next_object_ptr;
        }
    }
}
/* v.4.0--> */
/***********************************************************************

calc_residue_means  -  Calculate mean positions of each residue in the
                       plot

***********************************************************************/

void calc_residue_means(void)
{
  int iatom, iresid, natoms, ncoords, nresid;

  float mean_x, mean_y, mean_z;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Initialise the residue's average position */
          mean_x = 0.0;
          mean_y = 0.0;
          mean_z = 0.0;
          ncoords = 0;

          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all this residue's atoms */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* Update the mean coordinates */
              mean_x = mean_x + atom_ptr->original_x;
              mean_y = mean_y + atom_ptr->original_y;
              mean_z = mean_z + atom_ptr->original_z;
              ncoords++;

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }

          /* Calculate this residue's mean coordinates */
          if (ncoords > 0)
            {
              /* Get mean coords */
              mean_x = mean_x / ncoords;
              mean_y = mean_y / ncoords;
              mean_z = mean_z / ncoords;

              /* Store the object's original and flattened mean coords */
              residue_ptr->original_mean_x = mean_x;
              residue_ptr->original_mean_y = mean_y;
              residue_ptr->original_mean_z = mean_z;
            }

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

read_rcm_file  -  Read in the .rcm file, if present and match up the
                  relevant objects

***********************************************************************/

/* v.5.3 --> */
// void read_rcm_file(char pdb_name[FILENAME_LEN])
int read_rcm_file(char pdb_name[FILENAME_LEN])
/* <--v.5.3 */
{
/* v.5.3.4--> */
//  int iatom, iresid, natoms, nresid;
  int iatom, iresid, len, natoms, nresid;
/* <--v.5.3.4 */

/* v.5.0--> */
  char atom_no[6];
/* <--v.5.0 */
  char atom_name[5];
  char input_line[LINELEN + 1], rcm_name[FILENAME_LEN];
  char chain, res_name[4], res_num[6];

  int done, have_file;
  int line, natom_anchors, nres_anchors;
/* v.5.0--> */
  int atom_number;
/* <--v.5.0 */

  float x, y;
/* v.5.0--> */
  float anchor_weight;
/* <--v.5.0 */

/* v.5.0--> */
  // struct coordinate *atom_ptr;
  // struct object *object_ptr;
  // struct residue *residue_ptr;
  struct coordinate *atom_ptr, *last_atom_ptr;
  struct object *object_ptr, *dummy_object_ptr, *last_object_ptr;
  struct residue *residue_ptr, *last_residue_ptr;
/* <--v.5.0 */

  FILE *file_ptr;

  /* Initialise variables */
  have_file = TRUE;
  line = 0;
  natom_anchors = 0;
  nres_anchors = 0;
/* v.5.0--> */
  atom_number = 0;

  /* Initialise pointers */
  dummy_object_ptr = NULL;

  /* Loop through all the objects to get the last one in case need
     to create any new ones in the routine */
  last_object_ptr = get_last_object(first_object_ptr);

  /* Ditto for residues */
  last_residue_ptr = get_last_residue(first_residue_ptr);

  /* And for atoms */
  last_atom_ptr = get_last_atom(first_atom_ptr);
/* <--v.5.0 */

  /* Get the name of the .rcm file */
  get_filename(pdb_name,rcm_name,".rcm");
/* debug */
  printf("RCM file: %s\n",rcm_name);
  fflush(stdout);
/* debug */

  /* Open .rcm file to see if it is present */
  if ((file_ptr = fopen(rcm_name,"r")) ==  NULL)
    have_file = FALSE;

  /* If file is present, read through all the records */
  if (have_file == TRUE)
    {
      /* Read in the data from the PDB file */
      printf("\nReading in the residue restraints file [%s] ...\n",
             rcm_name);

      /* Loop through the PDB file, picking up all the atoms that
         are required */
      while (fgets(input_line,LINELEN,file_ptr)!= NULL)
        {
          /* If not a header record, get the details */
          if (!strncmp(input_line,"ATOM  ",6) ||
              !strncmp(input_line,"HETATM",6) ||
              !strncmp(input_line,"RESDUE",6))
            {
/* v.5.3.4--> */
              /* Truncate the string to strip off the newline */
              len = string_truncate(input_line,LINELEN);
/* <--v.5.3.4 */

              /* Get the details from this line */
              strncpy(res_name,input_line+17,3);
              res_name[3] = '\0';
              strncpy(res_num,input_line+22,5);
              res_num[5] = '\0';
              chain = input_line[21];
              if (chain == '-')
                chain = ' ';

              /* Get the x- and y-coordinates */
              sscanf(input_line+30," %f %f ",&x,&y);

/* v.5.0--> */
              /* Get the anchor weight */
/* v.5.3.4--> */
              if (len > 60)
/* <--v.5.3.4 */
                sscanf(input_line+60," %f ",&anchor_weight);
/* v.5.3.4--> */
              else
                anchor_weight = 1.0;
/* <--v.5.3.4 */

              /* Initialise flag */
              done = FALSE;

              /* If a dummy atom, need to create it */
              if (!strcmp(res_name,"..."))
                {
                  /* If this is the first one, create a new dummy
                     object and residue */
                  if (dummy_object_ptr == NULL)
                    {
                      /* Create a new object */
                      dummy_object_ptr
                        = create_object_record(&last_object_ptr);

                      /* Create a dummy residue record for it */
                      residue_ptr
                        = create_residue_record(&last_residue_ptr,
                                                res_name,res_num,chain);

                      /* Set type and store pointer to residue's object */
                      residue_ptr->residue_type = DUMMY_RESIDUE;
                      residue_ptr->object_ptr = dummy_object_ptr;

                      /* Update the details of the object to which
                         this residue belongs */
                      dummy_object_ptr->first_residue_ptr = residue_ptr;

                      /* Increment count of this object's residues */
                      dummy_object_ptr->nresidues++;

                      /* Store the object type */
                      dummy_object_ptr->object_type = DUMMY;
                    }

                  /* Otherwise, get the corresponding dummy residue */
                  else
                    residue_ptr = dummy_object_ptr->first_residue_ptr;

                  /* Get the atom name and number */
                  strncpy(atom_name,input_line+12,4);
                  atom_name[4] = '\0';

                  /* Increment atom-number */
                  atom_number++;
                  sprintf(atom_no,"%5d",atom_number);

                  /* Create an atom record */
                  atom_ptr = create_atom_record(&last_atom_ptr,atom_name,
                                                atom_no,x,y,0.0);

                  /* If residue's first atom, save pointer */
                  if (residue_ptr->first_atom_ptr == NULL)
                    residue_ptr->first_atom_ptr = atom_ptr;

                  /* Increment count of atoms in this residue */
                  residue_ptr->natoms++;

                  /* Get atom record to point to its residue record */
                  atom_ptr->residue_ptr = residue_ptr;

                  /* Increment count of this residue's
                     atom-anchor positions */
                  residue_ptr->nanchored_atoms++;

                  /* Store anchor position */
                  atom_ptr->have_anchor = TRUE;
                  atom_ptr->anchor_pstn_x = x;
                  atom_ptr->anchor_pstn_y = y;
/* v.5.0--> */
                  atom_ptr->anchor_weight = anchor_weight;
/* <--v.5.0 */
                  dummy_object_ptr->have_anchors = TRUE;
                  natom_anchors++;

                  /* Set flag that we're done */
                  done = TRUE;
                }
/* <--v.5.0 */

              /* Initialise pointer to the first stored object */
              object_ptr = first_object_ptr;
/* v.5.0--> 
              done = FALSE;
 <--v.5.0 */

              /* Loop through all objects to find the residue corresponding
                 to the anchor-position just read in */
              while (object_ptr != NULL && done == FALSE)
                {
                  /* Get pointer to this object's first residue */
                  residue_ptr = object_ptr->first_residue_ptr;
                  nresid = object_ptr->nresidues;
                  iresid = 0;

                  /* Loop over all this object's residues */
                  while (iresid < nresid && residue_ptr != NULL &&
                         done == FALSE)
                    {
                      /* See if this is the relevant residue */
                      if (!strncmp(residue_ptr->res_name,res_name,3) &&
                          !strncmp(residue_ptr->res_num,res_num,5) &&
                          residue_ptr->chain == chain)
                        {
                          /* If the anchor position refers to the CofM
                             of the residue, then store it */
                          if (!strncmp(input_line+12,"CofM",4))
                            {
                              residue_ptr->have_anchor = TRUE;
                              residue_ptr->anchor_pstn_x = x;
                              residue_ptr->anchor_pstn_y = y;
                              object_ptr->have_anchors = TRUE;
                              done = TRUE;
                              nres_anchors++;
                            }

                          /* Otherwise, take this to be the anchor for one
                             of this residue's atom-positions */
                          else
                            {
                              /* Get the atom name */
                              strncpy(atom_name,input_line+12,4);
                              atom_name[4] = '\0';

                              /* Get pointer to first residue's first atom */
                              atom_ptr = residue_ptr->first_atom_ptr;
                              iatom = 0;
                              natoms = residue_ptr->natoms;

                              /* Loop over all this residue's atoms */
                              while (iatom < natoms && atom_ptr != NULL &&
                                     done == FALSE)
                                {
                                  /* If this is the right atom, then
                                     store its anchor position */
                                  if (!strncmp(atom_ptr->atom_name,
                                               atom_name,4))
                                    {
                                      /* Increment count of this residue's
                                         atom-anchor positions */
                                      if (atom_ptr->have_anchor == FALSE)
                                        residue_ptr->nanchored_atoms++;

                                      /* Store anchor position */
                                      atom_ptr->have_anchor = TRUE;
                                      atom_ptr->anchor_pstn_x = x;
                                      atom_ptr->anchor_pstn_y = y;
                                      atom_ptr->anchor_weight
                                        = anchor_weight;
                                      object_ptr->have_anchors = TRUE;
                                      done = TRUE;
                                      natom_anchors++;

/* v.5.1--> */
                                      /* If generating plot for ligfit,
                                         and this is a ring atom, increase
                                         its anchor weight */
                                      if (ligfit == TRUE &&
                                          atom_ptr->ring_atom == TRUE &&
                                          anchor_weight == 1.0)
                                        atom_ptr->anchor_weight
                                          = RING_ANCHOR_WEIGHT;
/* <--v.5.1 */
                                    }

                                  /* Get pointer to the next atom */
                                  atom_ptr = atom_ptr->next;
                                  iatom++;
                                }

                              /* If appropriate atom not found, then
                                 display error message */
                              if (done == FALSE)
                                {
                                  /* Print message */
                                  printf("*** WARNING. Atom [%s] for ",
                                         atom_name);
                                  printf("residue [%s %s %c] in .rcm file",
                                         res_name,res_num,chain);
                                  printf(" not found.\n");
                                  printf("***          Restraint ignored.\n");
                                  Nwarnings++;
                                  done = TRUE;
                                }
                            }
                        }

                      /* Get pointer to the next residue */
                      residue_ptr = residue_ptr->next_residue_ptr;
                      iresid++;
                    }

                  /* Get pointer to next object */
                  object_ptr = object_ptr->next_object_ptr;
                }

              /* If residue not found, then print warning message */
              if (done == FALSE)
                {
                  /* Print message */
                  printf("*** WARNING. Residue [%s %s %c] in .rcm file",
                         res_name,res_num,chain);
                  printf(" not found. Restraint ignored.\n");
                  Nwarnings++;
                }
            }

          /* Increment line-count */
          line++;
        }
    }

  /* If have at least two anchor positions, then use anchoring */
  if ((nres_anchors + natom_anchors) > 1)
    {
      /* Set flag and show number of anchor-points read in */
      Have_Anchors = TRUE;
      printf("   Number of anchoring residues identified:    %d\n",
             nres_anchors);
      printf("   Number of atom-anchor positions identified: %d\n",
             natom_anchors);
    }

/* v.4.2--> */
  /* Close the .rcm file */
/* v.4.3--> */
  if (have_file == TRUE)
/* <--v.4.3 */
    fclose(file_ptr);
/* <--v.4.2 */

/* v.5.3 --> */
  /* Return whether rcm file was found */
  return(have_file);
/* <--v.5.3 */
}
/***********************************************************************

rotate_centres_of_mass  -  Apply the given rotation to the residue centres
                           of mass

***********************************************************************/

void rotate_centres_of_mass(int natoms, float coord_store[MAXATOMS][3],
                            float v[3][3])
{
  int iatom, icoord;
  float x, y, z;
  float new_coords[3];

  /* Loop over all the centres of mass */
  for (iatom = 0; iatom < natoms; iatom++)
    {
      /* Get the coordinates of this centre of mass */
      x = coord_store[iatom][0];
      y = coord_store[iatom][1];
      z = coord_store[iatom][2];

      /* Apply the transformation */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Calculate transformed coordinates */
              new_coords[icoord]
                =  ((x * v[0][icoord])
                    + (y * v[1][icoord])
                    + (z * v[2][icoord]));
        }

      /* Save the new coordinates */
      coord_store[iatom][0] = new_coords[0];
      coord_store[iatom][1] = new_coords[1];
      coord_store[iatom][2] = new_coords[2];
    }
}
/***********************************************************************

get_object_means  -  Calculate mean positions of eaach object's original
                     and flattened residues

***********************************************************************/

void get_object_means(struct object *store_object_ptr[MAXATOMS],
                      float coord_store[MAXATOMS][3],
                      float flattened_mean[MAXATOMS][2],int *nstored)
{
  int iatom, iresid, natoms, ncoords, nobj_coords, nresid;

  float flattened_mean_x, flattened_mean_y, mean_x, mean_y, mean_z;
  float flattened_object_mean_x, flattened_object_mean_y;
  float object_mean_x, object_mean_y, object_mean_z;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Initialise the object's average position */
  object_mean_x = 0.0;
  object_mean_y = 0.0;
  object_mean_z = 0.0;
  flattened_object_mean_x = 0.0;
  flattened_object_mean_y = 0.0;
  nobj_coords = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Initialise the residue's average position */
          mean_x = 0.0;
          mean_y = 0.0;
          mean_z = 0.0;
          flattened_mean_x = 0.0;
          flattened_mean_y = 0.0;
          ncoords = 0;

          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all this residue's atoms */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* Update the residue's mean coordinates */
              mean_x = mean_x + atom_ptr->original_x;
              mean_y = mean_y + atom_ptr->original_y;
              mean_z = mean_z + atom_ptr->original_z;
              ncoords++;

              /* Update mean coords of the current x-y coords of the
                 atom in the residue's flattened state */
              flattened_mean_x = flattened_mean_x + atom_ptr->x;
              flattened_mean_y = flattened_mean_y + atom_ptr->y;

              /* Repeat for current object */
              object_mean_x = object_mean_x + atom_ptr->original_x;
              object_mean_y = object_mean_y + atom_ptr->original_y;
              object_mean_z = object_mean_z + atom_ptr->original_z;
              flattened_object_mean_x = flattened_object_mean_x + atom_ptr->x;
              flattened_object_mean_y = flattened_object_mean_y + atom_ptr->y;
              nobj_coords++;

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }

          /* Get this residue's mean coordinates */
          if (ncoords > 0)
            {
              /* Get mean coords */
              mean_x = mean_x / ncoords;
              mean_y = mean_y / ncoords;
              mean_z = mean_z / ncoords;

              /* Calculate the flattened-out residue's mean position, and
                 store the coords */
              flattened_mean_x = flattened_mean_x / ncoords;
              flattened_mean_y = flattened_mean_y / ncoords;

              /* Store the residue's original and flattened mean coords */
              residue_ptr->flattened_mean_x = flattened_mean_x;
              residue_ptr->flattened_mean_y = flattened_mean_y;
            }

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Get this object's mean coordinates */
      if (nobj_coords > 0)
        {
          /* Get mean coords */
          object_mean_x = object_mean_x / nobj_coords;
          object_mean_y = object_mean_y / nobj_coords;
          object_mean_z = object_mean_z / nobj_coords;

          /* Store the current object pointer */
          store_object_ptr[*nstored] = object_ptr;

          /* Store the mean coords */
          coord_store[*nstored][0] = object_mean_x;
          coord_store[*nstored][1] = object_mean_y;
          coord_store[*nstored][2] = object_mean_z;

          /* Calculate the flattened-out object's mean position, and
             store the coords */
          flattened_object_mean_x = flattened_object_mean_x / nobj_coords;
          flattened_object_mean_y = flattened_object_mean_y / nobj_coords;
          flattened_mean[*nstored][0] = flattened_object_mean_x;
          flattened_mean[*nstored][1] = flattened_object_mean_y;

          /* Increment count of stored sets of coords */
          (*nstored)++;
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

place_interacting_residues  -  Lay out the current object's residues on
                               the appropriate side of the interface

***********************************************************************/

void place_interacting_residues(float coord_store[MAXATOMS][3],
                                struct object *store_object_ptr[MAXATOMS],
                                int nstored)
{
  int iatom, iresid, istored, natoms, nresid;

  float dx, dy, flattened_mean_x, flattened_mean_y, spread_factor, x_pos;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;
  struct object *object_ptr;

  /* Calculate the spread-factor */
  spread_factor = 2.0;

  /* Loop through all the stored positions to lay out the objects
     according to the principal component analysis */
  for (istored = 0; istored < nstored; istored++)
    {
      /* Get the corresponding object */
      object_ptr = store_object_ptr[istored];

      /* Get the mean position of the flattened object's coords */
      flattened_mean_x = object_ptr->first_residue_ptr->flattened_mean_x;
      flattened_mean_y = object_ptr->first_residue_ptr->flattened_mean_y;

      /* Calculate the adjustment to be made to each of this objects'
         flattened coords */

      /* The x-coordinate is just the transformed z-coord, with the
         residues from one domain/chain on the left side of it, and those
         from the other on the right */
      x_pos = (float) (spread_factor * fabs(coord_store[istored][2]));
      if (object_ptr->object_type == WATER)
        x_pos = 0.0;
      else if (object_ptr->interface == 1)
        x_pos = - x_pos;
      dx = x_pos - flattened_mean_x;

      /* The y-coordinate is just the transformed x-coord */
      dy = spread_factor * coord_store[istored][0] - flattened_mean_y;

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Get pointer to first residue's first atom */
          atom_ptr = residue_ptr->first_atom_ptr;
          iatom = 0;
          natoms = residue_ptr->natoms;

          /* Loop over all this residue's atoms */
          while (iatom < natoms && atom_ptr != NULL)
            {
              /* Transform the atom's coordinates */
              atom_ptr->x = atom_ptr->x + dx;
              atom_ptr->y = atom_ptr->y + dy;
              atom_ptr->z = 0.0;

              /* Get pointer to the next atom */
              atom_ptr = atom_ptr->next;
              iatom++;
            }

          /* Shift the flattened means by the same distance that
             the individual atomic coords have been shifted */
          residue_ptr->flattened_mean_x = residue_ptr->flattened_mean_x + dx;
          residue_ptr->flattened_mean_y = residue_ptr->flattened_mean_y + dy;

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }
    }
}
/***********************************************************************

lay_interface_objects_out  -  Lay out the interacting objects on either
                              side of an interface

***********************************************************************/

void lay_interface_objects_out(void)
{
  int nstored;

  float coord_store[MAXATOMS][3], flattened_mean[MAXATOMS][2];
  float mass_x, mass_y, mass_z;
  float transformation[3][3];

  struct object *store_object_ptr[MAXATOMS];

  /* Initialise variables */
  nstored = 0;

  /* Get the mean coordinates of the flattened objects */
  get_object_means(store_object_ptr,coord_store,flattened_mean,
                   &nstored);

  /* Calculate the centre of mass of the stored coordinates */
  centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

  /* Translate to mean coordinate to the origin */
  adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

  /* Calculate transformation to orient the principal axis along
     the x-direction */
  principal_components(nstored,transformation,coord_store);

  /* Apply this transformation to the centres of mass */
  rotate_centres_of_mass(nstored,coord_store,transformation);

  /* Place all the residues of this object on the appropriate
     side of the interface */
  place_interacting_residues(coord_store,store_object_ptr,nstored);
}
/* <--v.4.0 */
/***********************************************************************

lay_objects_out  -  Lay the objects out on the page in a starting
                    arrangement that approximately reflects the
                    relative interactions

***********************************************************************/

void lay_objects_out(void)
{
  int iatom, iresid, natoms, nresid;

  float spread_factor;
/* v.4.0--> */
  float x, y;
/* <--v.4.0 */

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
/* v.5.3.4--> */
      /* If object not already fixed, find where to place it */
      if (object_ptr->fixed == FALSE)
        {
/* <--v.5.3.4 */
          /* Calculate the spread-factor according to the object type */
          spread_factor = SPREAD_FACTOR;
          if (object_ptr->object_type == HYDROPHOBIC)
            spread_factor = (float) (10.0 * spread_factor);

          /* If this if a hydrophobic group or simplified H-group, use
             only a single atom to represent it */
          if (object_ptr->object_type == HYDROPHOBIC ||
              object_ptr->object_type == SIMPLE_HGROUP ||
              object_ptr->object_type == WATER)
            {
              /* Set object's position using its placement position */
/* v.4.0--> */
              /* object_ptr->minx = spread_factor * object_ptr->place_x;
                 object_ptr->maxx = spread_factor * object_ptr->place_x;
                 object_ptr->miny = spread_factor * object_ptr->place_y;
                 object_ptr->maxy = spread_factor * object_ptr->place_y; */

              /* Use the residue's centre-of-mass coords in the original
                 structure (times the spread factor) as starting coords
                 for object */ 
              x = spread_factor * object_ptr->place_x;
              y = spread_factor * object_ptr->place_y;

              /* If have an anchor-position for this object, use that instead */
              if (object_ptr->have_anchors == TRUE)
                {
                  /* Extract the anchor position */
                  x = object_ptr->first_residue_ptr->anchor_pstn_x;
                  y = object_ptr->first_residue_ptr->anchor_pstn_y;
                }

              /* Place the object */
              object_ptr->minx = x;
              object_ptr->maxx = x;
              object_ptr->miny = y;
              object_ptr->maxy = y;
/* <--v.4.0 */
            }
          else
            {
              /* Perform the fit for the current object */
              fit_object(object_ptr,FALSE);
            }

          /* Set the pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
             /* Get pointer to first residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all this residue's atoms */
              while (iatom < natoms && atom_ptr != NULL)
                {
                 /* If this if a hydrophobic group or simplified H-group,
                     set all atom coords to the placement position */
                  if (object_ptr->object_type == HYDROPHOBIC ||
                      object_ptr->object_type == SIMPLE_HGROUP ||
                      object_ptr->object_type == WATER)
                    {
                      /* Update the atom's coords */
                      atom_ptr->x = object_ptr->minx;
                      atom_ptr->y = object_ptr->miny;
                    }

                  /* Otherwise, move the object to its placement position */
                  else
                    {
                      /* Move object to its placement position */
                      atom_ptr->x = atom_ptr->x
                        + SPREAD_FACTOR * object_ptr->place_x;
                      atom_ptr->y = atom_ptr->y
                        + SPREAD_FACTOR * object_ptr->place_y;
                    }

                  /* Force z-coord to be zero */
                  atom_ptr->z = 0.0;

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Update residue- and object-boundaries */
          update_residue_boundaries(object_ptr);
/* v.5.3.4--> */
        }
/* <--v.5.3.4 */

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/***********************************************************************

get_atom_clash_energy  -  Calculate the energy for atom-clashes between
                          the given objects

***********************************************************************/

float get_atom_clash_energy(struct object *object1_ptr,
                            struct object *object2_ptr)
{
  int iatom, iresid, jatom, jresid, natoms1, nresid1, nresid2, natoms2;
  int object_type1, object_type2;
  int in_range;

  float atom_size1, atom_size2, distance, dist2, interact_dist;
  float energy, total_energy;
  float x1, x2, y1, y2;

  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;

  /* Initialise variables */
  total_energy = 0.0;

  /* Get the two object types */
  object_type1 = object1_ptr->object_type;
  object_type2 = object2_ptr->object_type;

  /* Set the pointer to the first object's first residue */
  residue1_ptr = object1_ptr->first_residue_ptr;
  nresid1 = object1_ptr->nresidues;
  iresid = 0;

  /* Loop over all this object's residues */
  while (iresid < nresid1 && residue1_ptr != NULL)
    {
      /* Set the pointer to the second object's first residue */
      residue2_ptr = object2_ptr->first_residue_ptr;
      nresid2 = object2_ptr->nresidues;
      jresid = 0;

      /* Loop over other object's residues */
      while (jresid < nresid2 && residue2_ptr != NULL)
        {
          /* Get the maximum atom-size of each residue */
          atom_size1 = residue1_ptr->max_atom_size;
          atom_size2 = residue2_ptr->max_atom_size;
          interact_dist = atom_size1 + atom_size2 + ATOM_INTERACT_DIST;

          /* Check whether the residues are within range of one another */
          in_range = TRUE;
          if ((residue2_ptr->minx > (residue1_ptr->maxx + interact_dist)) ||
              (residue1_ptr->minx > (residue2_ptr->maxx + interact_dist)) ||
              (residue2_ptr->miny > (residue1_ptr->maxy + interact_dist)) ||
              (residue1_ptr->miny > (residue2_ptr->maxy + interact_dist)))
            in_range = FALSE;

          /* If either residue doesn't have any atoms, then skip */
          if (residue1_ptr->natoms == 0 || residue2_ptr->natoms == 0)
            in_range = FALSE;

          /* If residues are within range, then check all distances
             between them */
          if (in_range == TRUE)
            {
              /* Get pointer to first residue's first atom */
              atom1_ptr = residue1_ptr->first_atom_ptr;

              /* Initialise atom count */
              iatom = 0;
              natoms1 = residue1_ptr->natoms;

              /* If this object is a hydrophobic group or a simplified
                 residue, then need only consider a single position */
              if (object_type1 == HYDROPHOBIC ||
                  object_type1 == SIMPLE_HGROUP)
                {
                  x1 = object1_ptr->minx;
                  y1 = object1_ptr->miny;
                  atom_size1 = object1_ptr->max_atom_size;
                  natoms1 = 1;
                }

              /* Loop over all the first residue's atoms */
              while (iatom < natoms1 && atom1_ptr != NULL)
                {
                  /* Get atom's coordinates and size */
                  if (object_type1 != HYDROPHOBIC &&
                      object_type1 != SIMPLE_HGROUP)
                    {
                      /* Get atom's coordinates */
                      x1 = atom1_ptr->x;
                      y1 = atom1_ptr->y;

                      /* Get the atom's size */
                      atom_size1 = atom1_ptr->atom_size;
                    }

                  /* Get pointer to other residue's first atom */
                  atom2_ptr = residue2_ptr->first_atom_ptr;
                  
                  /* Initialise atom count */
                  jatom = 0;
                  natoms2 = residue2_ptr->natoms;
                  
                  /* If this object is a hydrophobic group or a simplified
                     residue, then need only consider a single position */
                  if (object_type2 == HYDROPHOBIC ||
                      object_type2 == SIMPLE_HGROUP)
                    {
                      x2 = object2_ptr->minx;
                      y2 = object2_ptr->miny;
                      atom_size2 = object2_ptr->max_atom_size;
                      natoms2 = 1;
                    }

                  /* Loop over all the other residue's atoms */
                  while (jatom < natoms2 && atom2_ptr != NULL)
                    {
                      /* Get atom's coordinates and size */
                      if (object_type2 != HYDROPHOBIC &&
                          object_type2 != SIMPLE_HGROUP)
                        {
                          /* Get other atom's coordinates */
                          x2 = atom2_ptr->x;
                          y2 = atom2_ptr->y;

                          /* Get the other atom's size */
                          atom_size2 = atom2_ptr->atom_size;
                        }

                      /* Calculate the distance between the
                         two atoms */
                      dist2 = (x1 - x2) * (x1 - x2)
                        + (y1 - y2) * (y1 - y2);
                      distance = (float) sqrt((double) dist2);

                      /* Adjust distance to allow for the atom sizes */
                      if (distance > (atom_size1 + atom_size2))
                        distance = distance - (atom_size1 + atom_size2);
                      else
                        distance = 0.0;
                      dist2 = distance * distance;

                      /* Add energy of this interaction to the
                         overall energy sum */
                      if (dist2 < 0.001)
                        dist2 = (float) (0.001);
                      if (dist2 < ATOM_DIST2)
                        energy = (float) (Atom_Atom_Clash / dist2);
                      else
                        energy = 0.0;
                      total_energy = total_energy + energy;

                      /* Get pointer to other residue's next atom */
                      atom2_ptr = atom2_ptr->next;
                      jatom++;
                    }

                  /* Get pointer to the next atom */
                  atom1_ptr = atom1_ptr->next;
                  iatom++;
                }
            }

          /* Get pointer to the next residue */
          residue2_ptr = residue2_ptr->next_residue_ptr;
          jresid++;
        }

      /* Get pointer to the next residue */
      residue1_ptr = residue1_ptr->next_residue_ptr;
      iresid++;
    }

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

get_distance_energy  -  Calculate the energy given the bond-atom distance
                        and the atom-size

***********************************************************************/

float get_distance_energy(float distance, float atom_size)
{
  float energy;
  
  /* Adjust distance to allow for the atom sizes */
  if (distance > atom_size)
    distance = distance - atom_size;
  else
    distance = 0.0;

  /* If the distance is within the interaction distance, then
     compute bond-atom clash energy */
  if (distance < 0.03)
    distance = (float) 0.03;
  if (distance < BOND_INTERACT_DIST)
    energy = (float) (Bond_Atom_Clash / (distance * distance));
  else
    energy = 0.0;

  /* Return the calculated energy */
  return(energy);
}
/***********************************************************************

get_bond_clash_energy  -  Calculate the energy for atom-bond clashes
                          between the given bond and all objects/atoms

***********************************************************************/

float get_bond_clash_energy(struct bond *current_bond_ptr,
                            float matrix[4][3],int external_bond,
                            struct object *current_object_ptr)
{
  int iatom, iobject, iresid, natoms, nresid;
  int internal_bond, object_type;
  int in_range, want_atom, wanted;

  float atom_size, distance, interact_dist, length;
  float energy, total_energy;
  float max_x, max_y, min_x, min_y, x, x1, x2, y, y1, y2;
  float conv_x, conv_x1, conv_x2, conv_y, conv_y1, conv_y2;

  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;
  struct residue *residue_ptr;

  /* Initialise variables */
  total_energy = 0.0;

  /* Get the pointers to the two atoms at either end of the bond */
  atom1_ptr = current_bond_ptr->first_atom_ptr;
  atom2_ptr = current_bond_ptr->second_atom_ptr;

  /* Get the corresponding objects */
  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
  object2_ptr = atom2_ptr->residue_ptr->object_ptr;

  /* Check whether this is an internal bond */
  internal_bond = FALSE;
  if (object1_ptr == object2_ptr)
    internal_bond = TRUE;

  /* Get the coordinates of the two atoms */
  x1 = atom1_ptr->x;
  y1 = atom1_ptr->y;
  x2 = atom2_ptr->x;
  y2 = atom2_ptr->y;

  /* Determine max and min extents of the bond */
  min_x = x1;
  max_x = x1;
  min_y = y1;
  max_y = y1;
  if (x2 < min_x)
    min_x = x2;
  if (x2 > max_x)
    max_x = x2;
  if (y2 < min_y)
    min_y = y2;
  if (y2 > max_y)
    max_y = y2;

  /* Apply transformation to the atom coords */
  apply_xy_transformation(x1,y1,&conv_x1,&conv_y1,matrix);
  apply_xy_transformation(x2,y2,&conv_x2,&conv_y2,matrix);

  /* Calculate the bond length */
  length = conv_x2;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Set the pointer to the object's first residue */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;
      wanted = TRUE;

      /* If this if a hydrophobic group or simplified H-group, then
         is represented by a single atom */
      object_type = object_ptr->object_type;
      if (object_type == HYDROPHOBIC || object_type == SIMPLE_HGROUP)
        {
          /* Get atom size and coordinates */
          x = object_ptr->minx;
          y = object_ptr->miny;
          nresid = 1;
          natoms = 1;

          /* If bond is attached to this object at either end, then
             don't want to consider it */
          if (object_ptr == object1_ptr || object_ptr == object2_ptr)
            wanted = FALSE;
        }

      /* If the current bond is an internal one to a single object,
         then don't need to consider its clashes with other atoms in
         this object */
      if (internal_bond == TRUE && object_ptr == object1_ptr)
        wanted = FALSE;

      /* If the bond we're looking at is an external bond, then are
         only interested in its clashes with the object of interest */
      if (external_bond == TRUE && object_ptr != current_object_ptr)
        wanted = FALSE;

      /* Loop over all the object's residues */
      while (iresid < nresid && residue_ptr != NULL && wanted == TRUE)
        {
          /* Get the maximum atom-size of this residue */
          atom_size = residue_ptr->max_atom_size;
          interact_dist = atom_size + BOND_INTERACT_DIST;

          /* Check whether the residues is within range of the bond */
          in_range = TRUE;
          if ((min_x > (residue_ptr->maxx + interact_dist)) ||
              (max_x < (residue_ptr->minx - interact_dist)) ||
              (min_y > (residue_ptr->maxy + interact_dist)) ||
              (max_y < (residue_ptr->miny - interact_dist)))
            in_range = FALSE;

          /* If residue within range of the bond, then check all its
             atoms against the bond */
          if (in_range == TRUE)
            {
              /* Get pointer to residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;

              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* If object is a hydrophobic group or a simplified
                 residue, then only have a single position */
              if (object_type == HYDROPHOBIC ||
                  object_type == SIMPLE_HGROUP)
                natoms = 1;

              /* Loop over all the residue's atoms */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Get atom's coordinates and size */
                  if (object_type != HYDROPHOBIC &&
                      object_type != SIMPLE_HGROUP)
                    {
                      /* Get atom's coordinates */
                      x = atom_ptr->x;
                      y = atom_ptr->y;

                      /* Get the atom's size */
                      atom_size = atom_ptr->atom_size;
                    }

                  /* If current atom is at either end of the bond, then
                     don't want it */
                  want_atom = TRUE;
                  if (atom_ptr == atom1_ptr || atom_ptr == atom2_ptr)
                    want_atom = FALSE;

                  /* If atom is wanted, then calculate bond-clash energy */
                  if (want_atom == TRUE)
                    {
                      /* Initialise energy */
                      energy = 0.0;

                      /* Apply transformation to atom's coordinates */
                      apply_xy_transformation(x,y,&conv_x,&conv_y,
                                              matrix);
                      x = conv_x;
                      y = conv_y;

                      /* If atom's x-value lies within the length of
                         the bond, take its y-value as the distance */
                      if (x >= 0.0 && x <= length)
                        distance = (float) (fabs(y));

                      /* If atom is to the left of the bond, then get
                         its distance to the origin (where one end of
                         the bond currently lies) */
                      else if (x < 0)
                        distance = (float) sqrt(x * x + y * y);

                      /* Otherwise, calculate distance to the other
                         end of the bond */
                      else
                        distance = (float) sqrt((x - length) * (x - length)
                                        + y * y);

                      /* Calculate the energy based on this distance 
                         and the atom size */
                      energy = get_distance_energy(distance,atom_size);

                      /* Add to total */
                      total_energy = total_energy + energy;
                    }

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }
            }
          
          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

get_bond_overlap_energy  -  Calculate the energy for bond overlaps
                            between the given bond and all other bonds

***********************************************************************/

float get_bond_overlap_energy(struct bond *current_bond_ptr,
                              float matrix[4][3],int external_bond,
                              struct object *current_object_ptr)
{
  int bond_type, current_bond_type;
  int in_range, wanted;
  int point_object1, point_object2;

  float length;
  float overlap_energy, total_energy;
  float x1, x2, y1, y2;
  float other_x1, other_x2, other_y1, other_y2;
  float conv_x, conv_y, cut_x;
  float energy;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct coordinate *other_atom1_ptr, *other_atom2_ptr;
  struct object *object1_ptr, *object2_ptr;
  struct object *other_object1_ptr, *other_object2_ptr;

  /* Initialise variables */
  overlap_energy = Overlap_Score;
  total_energy = 0.0;

  /* Get the current bond's type */
  current_bond_type = current_bond_ptr->bond_type;

  /* Get the pointers to the two atoms at either end of the bond */
  atom1_ptr = current_bond_ptr->first_atom_ptr;
  atom2_ptr = current_bond_ptr->second_atom_ptr;

  /* Get the corresponding objects */
  object1_ptr = atom1_ptr->residue_ptr->object_ptr;
  object2_ptr = atom2_ptr->residue_ptr->object_ptr;

  /* Check whether either object at the end of the bond is a point
     object */
  point_object1 = FALSE;
  if (object1_ptr->object_type == HYDROPHOBIC ||
      object1_ptr->object_type == WATER ||
      object1_ptr->object_type == SIMPLE_HGROUP)
    point_object1 = TRUE;
  point_object2 = FALSE;
  if (object2_ptr->object_type == HYDROPHOBIC ||
      object2_ptr->object_type == WATER ||
      object2_ptr->object_type == SIMPLE_HGROUP)
    point_object2 = TRUE;

  /* Get the transformed coordinates of the two atoms */
  x1 = atom1_ptr->x;
  y1 = atom1_ptr->y;
  x2 = atom2_ptr->x;
  y2 = atom2_ptr->y;

  /* Apply transformation to the atom coords */
  apply_xy_transformation(x1,y1,&conv_x,&conv_y,matrix);
  x1 = conv_x;
  y1 = conv_y;
  apply_xy_transformation(x2,y2,&conv_x,&conv_y,matrix);
  x2 = conv_x;
  y2 = conv_y;

  /* Calculate the bond length */
  length = x2;

  /* Get pointer to first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds, one by one */
  while (bond_ptr != NULL)
    {
      /* Initialise flag giving whether this bond is wanted */
      wanted = TRUE;

      /* Get the bond type */
      bond_type = bond_ptr->bond_type;

      /* Check that this isn't the current bond */
      if (bond_ptr == current_bond_ptr)
        wanted = FALSE;

      /* Not interested in internal H-bonds  */
      if (bond_type == INTERNAL || bond_type == DELETED)
        wanted = FALSE;

      /* Get the pointers to the two atoms at either end of the bond */
      other_atom1_ptr = bond_ptr->first_atom_ptr;
      other_atom2_ptr = bond_ptr->second_atom_ptr;
          
      /* Check that the two bonds don't have any atoms in common */
      if (atom1_ptr == other_atom1_ptr || atom1_ptr == other_atom2_ptr ||
          atom2_ptr == other_atom1_ptr || atom2_ptr == other_atom2_ptr)
        wanted = FALSE;

      /* Get the corresponding objects */
      other_object1_ptr = other_atom1_ptr->residue_ptr->object_ptr;
      other_object2_ptr = other_atom2_ptr->residue_ptr->object_ptr;

      /* If bonds both spring from a hydrophobic group, or
         a simplified residue, then don't want to check for
         an overlap */
      if (point_object1 == TRUE &&
          (object1_ptr == other_object1_ptr || 
           object1_ptr == other_object2_ptr))
        wanted = FALSE;
      if (point_object2 == TRUE &&
          (object2_ptr == other_object1_ptr || 
           object2_ptr == other_object2_ptr))
        wanted = FALSE;

      /* If either bond is a contact bond which doesn't spring from a
         hydrophobic group (eg is between ligand and an H-group), then
         not interested in it */
      if (current_bond_type == CONTACT &&
          object1_ptr->object_type != HYDROPHOBIC &&
          object2_ptr->object_type != HYDROPHOBIC)
        wanted = FALSE;
      if (bond_type == CONTACT &&
          other_object1_ptr->object_type != HYDROPHOBIC &&
          other_object2_ptr->object_type != HYDROPHOBIC)
        wanted = FALSE;

      /* If the bond we're looking at is an external bond,
         then are only interested in its overlaps with bonds
         involving the object of interest */
      if (external_bond == TRUE &&
          other_object1_ptr != current_object_ptr &&
          other_object2_ptr != current_object_ptr)
        wanted = FALSE;

      /* If bond is still wanted, then check if within range */
      if (wanted == TRUE)
        {
          /* Get the coordinates of the two atoms */
          other_x1 = other_atom1_ptr->x;
          other_y1 = other_atom1_ptr->y;
          other_x2 = other_atom2_ptr->x;
          other_y2 = other_atom2_ptr->y;

          /* Apply transformation to the atom coords */
          apply_xy_transformation(other_x1,other_y1,&conv_x,&conv_y,
                                  matrix);
          other_x1 = conv_x;
          other_y1 = conv_y;
          apply_xy_transformation(other_x2,other_y2,&conv_x,&conv_y,
                                  matrix);
          other_x2 = conv_x;
          other_y2 = conv_y;

          /* Check whether the bonds are within range of one another */
          in_range = TRUE;
          if (other_y1 < 0.0 && other_y2 < 0.0)
            in_range = FALSE;
          else if (other_y1 > 0.0 && other_y2 > 0.0)
            in_range = FALSE;
          else if (other_x1 < 0.0 && other_x2 < 0.0)
            in_range = FALSE;
          else if (other_x1 > x2 && other_x2 > x2)
            in_range = FALSE;

          /* If the two bonds are within range of one another, calculate
             their cross-over point */
          if (in_range == TRUE)
            {
              /* If second bond is not parallel to the x-axis calculate
                 where it cuts the x-axis */
              if (fabs(other_y2 - other_y1) > 0.01)
                {           
                  cut_x = other_x2 - other_y2 * (other_x2 - other_x1)
                      / (other_y2 - other_y1);

                  /* If the cut-point is within the length of the first
                     bond, then the bonds cross one another */
                  if (cut_x >= 0.0 && cut_x <= length)
                    {
                      if (current_bond_type == CONTACT ||
                          bond_type == CONTACT)
                        energy = (float) (0.0005);
                      else
                        energy = overlap_energy;
                      total_energy = total_energy + energy;
                    }
                }
            }
        }

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Return the calculated energy */
  return(total_energy);
}
/***********************************************************************

calculate_energy  -  Calculate the energy between the current object and
                     all the others

***********************************************************************/

float calculate_energy(struct object *current_object_ptr,
                       float *atom_clash,float *bond_length,
/* v.4.0--> */
/*                     float *bond_clash,float *bond_overlap) */
                       float *bond_clash,float *bond_overlap,
/* v.5.0--> */
/*                       float *boundary_energy,float *rel_pstn_energy) */
                       float *boundary_energy,float *rel_pstn_energy,
                       int rel_anchors)
/* <--v.5.0 */
/* <--v.4.0 */
{
  int bond_type;
  int external_bond, in_range, to_ligand, wanted;

  float coord_store[3][3], matrix[4][3];
  float atom_size1, atom_size2, energy, interact_dist;
  float atom_clash_energy, total_atom_clash_energy;
  float bond_length_energy, total_bond_length_energy;
  float bond_clash_energy, total_bond_clash_energy;
  float bond_overlap_energy, total_bond_overlap_energy;
  float internal_energy;
  float diff, length, x1, x2, y1, y2;

  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct object *object_ptr, *object1_ptr, *object2_ptr;

  /* Initialise variables */
  energy = 0.0;
  total_atom_clash_energy = 0.0;
  total_bond_length_energy = 0.0;
  total_bond_clash_energy = 0.0;
  total_bond_overlap_energy = 0.0;
/* v.4.0--> */
  *rel_pstn_energy = 0.0;
  *boundary_energy = 0.0;
/* <--v.4.0 */

  /* Get this object's internal energy */
  internal_energy = current_object_ptr->internal_energy;  

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects to calculate the atom-atom energies */
  while (object_ptr != NULL)
    {
      /* If this is not the same as the current object, calculate
         the energy between them */
/* v.5.0--> */
      // if (object_ptr != current_object_ptr)

      /* Assume wanted */
      wanted = TRUE;

      /* If same object, then not wanted */
      if (object_ptr == current_object_ptr)
        wanted = FALSE;

      /* If either is the dummy object, then don't want it if other
         object is the ligand */
      if ((current_object_ptr->object_type == DUMMY &&
           object_ptr->object_type == LIGAND) ||
          (object_ptr->object_type == DUMMY &&
           current_object_ptr->object_type == LIGAND))
        wanted = FALSE;

      /* If wanted, calculate energy between these two objects */
      if (wanted == TRUE)
/* <--v.5.0 */
        {
          /* Get the maximum atom-size of each object */
          atom_size1 = current_object_ptr->max_atom_size;
          atom_size2 = object_ptr->max_atom_size;
          interact_dist = atom_size1 + atom_size2 + ATOM_INTERACT_DIST;

          /* Check whether the objects are within range of one another */
          in_range = TRUE;
          if ((object_ptr->minx >
               (current_object_ptr->maxx + interact_dist)) ||
              (current_object_ptr->minx >
               (object_ptr->maxx + interact_dist)) ||
              (object_ptr->miny >
               (current_object_ptr->maxy + interact_dist)) ||
              (current_object_ptr->miny >
               (object_ptr->maxy + interact_dist)))
            in_range = FALSE;

          /* If objects are within range, then check all distances
             between their constituent residues */
          if (in_range == TRUE)
            {
              /* Get the energy between the current pair of objects */
              atom_clash_energy
                = get_atom_clash_energy(current_object_ptr,object_ptr);

              /* Add the energy to the current object's total energy */
              total_atom_clash_energy = total_atom_clash_energy
                + atom_clash_energy;
            }
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Loop through all the elastic bonds involving the current object
     to get their contributions to the energy */

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

/* v.5.0--> */
  /* If current object is the dummy object, don't need to check bonds */
  if (current_object_ptr->object_type == DUMMY)
    bond_ptr = NULL;
/* <--v.5.0 */

  /* Loop through all the bonds, one by one */
  while (bond_ptr != NULL)
    {
      /* Initialise flag giving whether this bond is wanted */
      to_ligand = FALSE;
      wanted = TRUE;

      /* Get the bond type */
      bond_type = bond_ptr->bond_type;

      /* Only interested in external elastic bonds */
/* v.4.1.2--> */
/*      if (bond_ptr->elastic == FALSE || bond_type == INTERNAL) */
      if (bond_ptr->elastic == FALSE ||
          (bond_type == INTERNAL && Include->Internal_Hbonds == FALSE))
/* <--v.4.1.2 */
        wanted = FALSE;

      /* Get the pointers to the two atoms at either end of the bond */
      atom1_ptr = bond_ptr->first_atom_ptr;
      atom2_ptr = bond_ptr->second_atom_ptr;

      /* Get the objects these two atoms belong to */
      object1_ptr = atom1_ptr->residue_ptr->object_ptr;
      object2_ptr = atom2_ptr->residue_ptr->object_ptr;

      /* If this bond does not involve the current object, then
         mark it as an external one */
      external_bond = FALSE;
      if (object1_ptr != current_object_ptr &&
          object2_ptr != current_object_ptr)
        external_bond = TRUE;

      /* If this a wanted bond, get coords of its two ends */
      if (wanted == TRUE)
        {
          /* Get the atom coordinates */
          x1 = atom1_ptr->x;
          y1 = atom1_ptr->y;
          x2 = atom2_ptr->x;
          y2 = atom2_ptr->y;

          /* If this is an external bond, check whether it is too far
             from the current object to ever interact with it */
          if (external_bond == TRUE)
            {
              /* Get maximum interaction distance */
              interact_dist = current_object_ptr->max_atom_size
                + BOND_INTERACT_DIST;

              /* Check whether bond might be within range of the current
                 object */
              if (x1 > current_object_ptr->maxx + interact_dist &&
                  x2 > current_object_ptr->maxx + interact_dist)
                wanted = FALSE;
              else if (x1 < current_object_ptr->minx - interact_dist &&
                       x2 < current_object_ptr->minx - interact_dist)
                wanted = FALSE;
              else if (y1 > current_object_ptr->maxy + interact_dist &&
                       y2 > current_object_ptr->maxy + interact_dist)
                wanted = FALSE;
              else if (y1 < current_object_ptr->miny - interact_dist &&
                       y2 < current_object_ptr->miny - interact_dist)
                wanted = FALSE;
            }
        }

      /* If bond still wanted, then process it */
      if (wanted == TRUE)
        {
          /* Store the coordinates of the two atoms defining this bond */
          coord_store[0][0] = x1;
          coord_store[0][1] = y1;
          coord_store[0][2] = 0.0;
          coord_store[1][0] = x2;
          coord_store[1][1] = y2;
          coord_store[1][2] = 0.0;
              
          /* Calculate the transformation that will put this bond
             along the x-axis, with atom 1 at the origin and atom 2
             pointing in the +ve x-direction */
          line_transformation(coord_store,matrix,0.0);

          /* If this is any bond other than an internal one, calculate
             its length and compare with actual length */
/* v.5.4--> */
//          if ((bond_type == HBOND || bond_type == CONTACT ||
          if ((bond_type == HBOND || bond_type == CONTACT ||
               bond_type == SALT_BRIDGE || bond_type == DISULPHIDE ||
/* <--v.5.4 */
              (bond_type == COVALENT && bond_ptr->elastic == TRUE)) &&
              external_bond == FALSE)
            {
              /* Calculate length of bond */
              length = (float) sqrt((x2 - x1) * (x2 - x1)
                                 + (y2 - y1) * (y2 - y1));

              /* Get the difference from the bond's original length (times
                 the required "stretch factor" which determines how
                 closely the H-bond groups are to be packed in around
                 the ligand) */
              if (bond_type == COVALENT)
                diff = length - bond_ptr->bond_length;
              else
                diff = length - Bond_Stretch * bond_ptr->bond_length;
              diff = diff * diff;

              /* Calculate the energy due to the discrepancy in the
                 bond's length */
              if (bond_type == COVALENT)
                bond_length_energy = (float) (10.0 * diff);
/* v.5.4--> */
//              else if (bond_type == HBOND)
              else if (bond_type == HBOND ||
                       bond_type == SALT_BRIDGE ||
                       bond_type == DISULPHIDE)
/* <--v.5.4 */
                bond_length_energy = HB_Weight * diff;
              else
              {
                bond_length_energy = Nonbond_Weight * diff;
                if (object1_ptr->object_type != HYDROPHOBIC &&
                    object2_ptr->object_type != HYDROPHOBIC)
                  bond_length_energy = (float) (bond_length_energy / 100.0);
              }

              /* Check if either atom belongs to the ligand */
              if (atom1_ptr->residue_ptr->inligand == TRUE ||
                  atom2_ptr->residue_ptr->inligand == TRUE)
                to_ligand = TRUE;

              /* If the bond involves the ligand, then double its
                 weight */
/* v.5.4--> */
//              if (bond_type == HBOND && to_ligand == TRUE)
              if ((bond_type == HBOND ||
                   bond_type == SALT_BRIDGE ||
                   bond_type == DISULPHIDE) && to_ligand == TRUE)
/* <--v.5.4 */
                bond_length_energy = 3 * bond_length_energy;

/* v.5.3.4--> */
              /* For DIMPLOT, reduce strength of the bond length term */
              if (dimfit == TRUE)
                bond_length_energy = bond_length_energy / 1000;
/* <--v.5.3.4 */

              /* Add to total bond-length energy for this object */
              total_bond_length_energy = total_bond_length_energy
                + bond_length_energy;
            }

          /* Calculate energy from any clashes between the bond and atoms */
          if (bond_type == CONTACT)
              bond_clash_energy = 0.0;
          else
            {
              bond_clash_energy
                = get_bond_clash_energy(bond_ptr,matrix,external_bond,
                                        current_object_ptr);
            }

/* v.5.3.4--> */
          /* For DIMPLOT, reduce strength of the bond length term */
          if (dimfit == TRUE)
            bond_clash_energy = bond_clash_energy / 1000;
/* <--v.5.3.4 */

          /* Calculate energy from bond overlaps */
          bond_overlap_energy
            = get_bond_overlap_energy(bond_ptr,matrix,external_bond,
                                      current_object_ptr);

/* v.5.3.4--> */
          /* For DIMPLOT, reduce strength of the bond length term */
          if (dimfit == TRUE)
            bond_overlap_energy = bond_overlap_energy / 1000;
/* <--v.5.3.4 */

          /* Update energy totals */
          total_bond_clash_energy = total_bond_clash_energy
            + bond_clash_energy;
          total_bond_overlap_energy = total_bond_overlap_energy
            + bond_overlap_energy;
        }

      /* Get pointer to the next bond */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

/* v.4.0--> */
  /* Calculate object's anchor energy or its relative position energy */
/* v.5.0--> */
  if (current_object_ptr->object_type != DUMMY)
    {
      *rel_pstn_energy = 0.0;
/* <--v.5.0 */
      if (Have_Anchors == TRUE && current_object_ptr->have_anchors == TRUE)
/* v.5.0--> */
/*        *rel_pstn_energy = get_anchor_energy(current_object_ptr); */
        {
          /* If calculating relative anchor energy, then get that */
          if (rel_anchors == TRUE)
            *rel_pstn_energy
              = get_relative_anchor_energy(current_object_ptr);

          /* Otherwise, calculate absolute anchor energy */
          else
            *rel_pstn_energy = get_anchor_energy(current_object_ptr);
        }
/* <--v.5.0 */
      else if (Rel_Pstn_Energy_Weight > 0.0)
          *rel_pstn_energy = get_rel_pstn_energy(current_object_ptr);
    }

  /* For interface plot, calculate boundary energy of object */
  if (Interface_Plot == TRUE)
    *boundary_energy = get_boundary_energy(current_object_ptr);
/* <--v.4.0 */

  /* Calculate the total energy from the sum of its contributions */
  energy = total_atom_clash_energy + total_bond_length_energy
/* v.4.0--> */
/*    + total_bond_clash_energy + total_bond_overlap_energy; */
    + total_bond_clash_energy + total_bond_overlap_energy
      + (*rel_pstn_energy);
/* <--v.4.0 */

  /* Store this object's total energy */
/* v.4.0--> */
/*  current_object_ptr->total_energy = internal_energy + energy / 2.0; */
  current_object_ptr->total_energy
    = (float) (internal_energy + (*boundary_energy) + energy / 2.0);
/* <--v.4.0 */
  energy = current_object_ptr->total_energy;

  /* Return the calculated energies */
  *atom_clash = total_atom_clash_energy;
  *bond_length = total_bond_length_energy;
  *bond_clash = total_bond_clash_energy;
  *bond_overlap = total_bond_overlap_energy;
  return(energy);
}
/***********************************************************************

make_random_move  -  Apply a random translation to the given object

***********************************************************************/

void make_random_move(struct object *object_ptr,float max_move_size,
                      int *rseed)
{
  int iseed;
  int got_dir;

  float length, length2, move_size, move_x, move_y;
  float x, y, z;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Get the random number direction for the move */
  got_dir = FALSE;
  while (got_dir == FALSE)
    {
      /* Get trial vector for direction of move */
      x = get_random_number(&iseed,Random_Start);
      y = get_random_number(&iseed,Random_Start);

      /* Convert so that values lie between -1 and +1 */
      x = (float) (2.0 * x - 1.0);
      y = (float) (2.0 * y - 1.0);

      /* Check that the point representing this direction lies within
         a unit circle (to ensure direction is genuinely random */
      length2 = x * x + y * y;
      if (length2 <= 1.0)
        got_dir = TRUE;
    }

  /* Calculate the move size */
  z = get_random_number(&iseed,Random_Start);
  move_size = z * max_move_size;

  /* Make the size of the move equal to the required size */
  length = (float) sqrt(length2);
  move_x = x * move_size / length;
  move_y = y * move_size / length;

  /* Apply the move to the current object */
  move_object(object_ptr,move_x,move_y,0.0);

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

make_random_rotation  -  Apply a random rotation or swing to the given
                         object

***********************************************************************/

void make_random_rotation(struct object *object_ptr,float max_angle,
                          int *rseed,int swing)
{
  int iatom, natoms;
  int i, iseed, j, nstored;
  int swing_about_origin, wanted;

  float angle, cos_theta, sin_theta, matrix[3][3], x, y;
  float random_no;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
  struct coordinate *stored_atom_ptr[MAXATOMS];
  struct object *object1_ptr, *object2_ptr;
  struct residue *residue_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;
  x = y = 0.0;
  nstored = 0;

  /* Initialise both matrices */
  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
        {
          if (i == j)
            matrix[i][j] = (float) (1.0);
          else
            matrix[i][j] = 0.0;
          /*      inverse_matrix[i][j] = matrix[i][j];  */
        }
    }

  /* Generate a random number */
  random_no = get_random_number(&iseed,Random_Start);

  /* Convert so that values lie between plus or minus of the 
     maximum rotation angle */
  angle = (float) (2.0 * max_angle * random_no - max_angle);

  /* Determine whether to swing about the origin */
  swing_about_origin = FALSE;
  if (swing == TRUE)
    {
      random_no = get_random_number(&iseed,Random_Start);
      if (random_no < 0.1)
        swing_about_origin = TRUE;
    }

  /* Calculate sines and cosines of the angle */
  cos_theta = (float) cos(angle);
  sin_theta = (float) sin(angle);

  /* Form the appropriate rotation matrix */
  matrix[0][0] = cos_theta;
  matrix[1][0] = sin_theta;
  matrix[0][1] = - sin_theta;
  matrix[1][1] = cos_theta;

  /* Form the inverse matrix */
  /*  inverse_matrix[0][0] = cos_theta;
  inverse_matrix[0][1] = - sin_theta;
  inverse_matrix[1][0] = sin_theta;
  inverse_matrix[1][1] = cos_theta; */

  /* If performing a swing, get all the bonds from this object to
     other objects */
  if (swing == TRUE && swing_about_origin == FALSE)
    {
      /* Initialise pointer to start of linked list of bonds */
      bond_ptr = first_bond_ptr;

      /* Loop through all the bonds in the linked list */
      while (bond_ptr != NULL)
        {
          /* Check if this is an elastic bond */
          if (bond_ptr->elastic == TRUE)
            {
              /* Get the pointers to the two atoms at either end
                 of the bond, and the objects they belong to */
              atom1_ptr = bond_ptr->first_atom_ptr;
              atom2_ptr = bond_ptr->second_atom_ptr;
              object1_ptr = atom1_ptr->residue_ptr->object_ptr;
              object2_ptr = atom2_ptr->residue_ptr->object_ptr;

              /* Check whether one atom is on the current object and
                 the other is not */
              wanted = FALSE;
              if ((object1_ptr == object_ptr && object2_ptr != object_ptr) ||
                  (object2_ptr == object_ptr && object1_ptr != object_ptr))
                wanted = TRUE;

              /* For H-group, prefer to swing about an H-bond */
              if (object_ptr->object_type == HGROUP &&
/* v.5.4--> */
//                  bond_ptr->bond_type != HBOND)
                  bond_ptr->bond_type != HBOND &&
                  bond_ptr->bond_type != SALT_BRIDGE &&
                  bond_ptr->bond_type != DISULPHIDE)
/* <--v.5.4 */
                wanted = FALSE;

              /* If this is a possible bond, then store pointer to the
                 atom at the other end of it */
              if (wanted == TRUE && nstored < MAXATOMS)
                {
                  if (object1_ptr == object_ptr)
                    stored_atom_ptr[nstored] = atom2_ptr;
                  else
                    stored_atom_ptr[nstored] = atom1_ptr;
                  nstored++;
                }
            }

          /* Get pointer to next atom in linked-list */
          bond_ptr = bond_ptr->next_bond_ptr;
        }

      /* If have some possible bonds, then pick one of them to be the
         pivotal one */
      if (nstored > 0)
        {
          /* Generate a random number */
          random_no = get_random_number(&iseed,Random_Start);

          /* Get the corresponding atom position */
          iatom = (int) (nstored * random_no);
          if (iatom < 0)
            iatom = 0;
          else if (iatom > nstored - 1)
            iatom = nstored - 1;

          /* Get the corresponding atom pointer and its coordinates */
          atom_ptr = stored_atom_ptr[iatom];
          x = atom_ptr->x;
          y = atom_ptr->y;
        }
    }

  /* If performing rotation about one of its atoms, transform the object
     to place that atom at the origin */
  else if (swing == FALSE)
    {
      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;

      /* Get the number of atoms and decide which is to be the pivot
         about which to perform the swivel */
      natoms = residue_ptr->natoms;
      random_no = get_random_number(&iseed,Random_Start);
      iatom = (int) (natoms * random_no + 1.0);
      if (iatom < 1)
        natoms = 1;
      else if (iatom < natoms)
        natoms = iatom;

      /* Get pointer to this residue's first atom */
      atom_ptr = residue_ptr->first_atom_ptr;
      iatom = 0;

      /* Loop through all the atoms until get to the one about
         which the swivel is to be performed */
      while (iatom < natoms && atom_ptr != NULL)
        {
          /* Store the coordinates of this atom */
          x = atom_ptr->x;
          y = atom_ptr->y;
        
          /* Get pointer to the next atom */
          atom_ptr = atom_ptr->next;
          iatom++;
        }
    }

  /* Move object to place the pivotal atom at the origin */
  move_object(object_ptr,x,y,0.0);

  /* Apply the rotation matrix to the current object */
  rotate_object(object_ptr,matrix);

  /* Move object back to where it came from */
  move_object(object_ptr,-x,-y,0.0);

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

make_random_flip  -  Make a random flip of the object about one of its
                     rotatable bonds

***********************************************************************/

void make_random_flip(struct object *object_ptr,int *rseed,int flip_all)
{
  int iseed;
  int flip_bond, ibond;
/* v.5.2.1 --> */
//  int both_ends, clash, done, use_end;
  int both_ends, calc_en, done, use_end;
/* <-- v.5.2.1 */

  float x;
  float internal_energy;

  struct bond *bond_ptr;
  struct object_bond *object_bond_ptr;

/* v.5.2.1 --> */
  /* Initialise variables */
  calc_en = TRUE;
/* <-- v.5.2.1 */

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Save object's current internal energy */
  internal_energy = object_ptr->internal_energy;

  /* Get a random number */
  x = get_random_number(&iseed,Random_Start);

  /* Convert into the sequential number of the rotatable bond to
     be flipped */
  flip_bond = (int) (object_ptr->nrot_bonds * x);

  /* Initialise bond count */
  ibond = 0;
  both_ends = FALSE;
  done = FALSE;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;
      
  /* Loop through all object's bonds */
  while (object_bond_ptr != NULL && done == FALSE)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* Check whether this bond is a rotatable one */
      if (bond_ptr->rotatable_bond == TRUE)
        {
          /* If this is the bond selected for the flip, then perform
             the flip */
          if (ibond == flip_bond)
            {
              /* Determine which end of this bond to flip */
              x = get_random_number(&iseed,Random_Start);
              if (x <= 0.5)
                use_end = FIRST;
              else
                use_end = SECOND;

              /* Determine whether to flip whole object */
              if (flip_all == FALSE)
                {
                  /* Use previously-generated random number to decide
                     whether to flip the whole object */
                  if (fabs(x - 0.5) < 0.05)
/* v.5.1--> */
                    {
/* <--v.5.1 */
                      both_ends = TRUE;
/* v.5.1--> */
                      use_end = FIRST;
                    }
/* <--v.5.1 */
                }

              /* Flip about this bond */
/* v.5.2.1 --> */
//              bond_flip(object_ptr,bond_ptr,use_end,both_ends,flip_all);
              bond_flip(object_ptr,bond_ptr,use_end,both_ends,flip_all,
                        calc_en);
/* <-- v.5.2.1 */

              /* Set flag that we are done */
              done = TRUE;
            }

          /* Increment count of rotatable bonds */
          ibond++;
        }
      
      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

make_random_atom_flip  -  Make a random flip of the object about one of
                          its atoms

***********************************************************************/

void make_random_atom_flip(struct object *object_ptr,int *rseed,
                           int flip_about_x)
{
  int flip_atom, iseed;
  int iatom, natoms;

  float x, y;
  float random_no;

  struct coordinate *atom_ptr;
  struct residue *residue_ptr;

  /* Get the current random number generator seed */
  iseed = *rseed;

  /* Set the pointer to the first of this object's residues and get the
     number of atoms */
  residue_ptr = object_ptr->first_residue_ptr;
  natoms = residue_ptr->natoms;

  /* Generate a random number */
  random_no = get_random_number(&iseed,Random_Start);

  /* Convert into the sequential number of the atom to be flipped */
  flip_atom = (int) (natoms * random_no + 1.0);
  if (flip_atom < 1)
    natoms = 1;
  else if (flip_atom < natoms)
    natoms = flip_atom;

  /* Initialise atom count */
  iatom = 0;

  /* Get pointer to the first of this residue's atoms */
  atom_ptr = residue_ptr->first_atom_ptr;

  /* Loop through all the atoms until get to the one about
     which the swivel is to be performed */
  while (iatom < natoms && atom_ptr != NULL)
    {
      /* Store the coordinates of this atom */
      x = atom_ptr->x;
      y = atom_ptr->y;
        
      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
      iatom++;
    }

  /* Move object to the origin */
  move_object(object_ptr,x,y,0.0);

  /* Apply the flip to the current object */
  if (flip_about_x == TRUE)
    flip_object(object_ptr,-1,1);
  else
    flip_object(object_ptr,1,-1);

  /* Move object back */
  move_object(object_ptr,-x,-y,0.0);

  /* Return the current random number seed */
  *rseed = iseed;
}
/***********************************************************************

get_total_energy  -  Calculate and print total energy for all objects

***********************************************************************/

/* v.4.3--> */
/* float get_total_energy(int loop,int print_headings,int print_data) */
float get_total_energy(int loop,int print_headings,int print_data,
                       int write_to_res,int rel_anchors)
/* <--v.4.3 */
{
  int i, iobject;
/* v.5.0--> */
//  int false = FALSE;
/* <--v.5.0 */

  float total_energy;
  float atom_clash_energy, bond_length_energy, bond_clash_energy;
  float bond_overlap_energy, internal_energy;
/* v.4.0--> */
  float boundary_energy, rel_pstn_energy;
/* <--v.4.0 */
  float energy_sum[NTERMS], energy_term[NTERMS];

  struct object *object_ptr;

  /* If headings to be printed, then do so */
  if (print_headings == TRUE)
    {
      printf("\n");
/* v.4.0--> */
/*      printf("                      Atom     Bond     Bond     Bond\n");
      printf(" Loop     Internal    Clash   Length    Clash  overlap");
      printf("       Total\n"); */
      printf("                 Atom    Bond    Bond    Bond ");
      printf("  Inter-   Resid           \n");
      printf("Loop  Internal   Clash  Length   Clash overlap");
      printf("   face    pstns      Total\n");
/* <--v.4.0 */
    }

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;
  for (i = 0; i < NTERMS; i++)
    energy_sum[i] = 0.0;
  total_energy = 0.0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Get the object's internal energy */
      internal_energy = object_ptr->internal_energy;

      /* Calculate the energy between this object and all
         the others */
      calculate_energy(object_ptr,
                       &atom_clash_energy,
                       &bond_length_energy,
                       &bond_clash_energy,
/* v.4.0--> */
/*                                     &bond_overlap_energy); */
                       &bond_overlap_energy,
                       &boundary_energy,
/* v.5.0--> */
/*                       &rel_pstn_energy); */
                       &rel_pstn_energy,rel_anchors);
/* <--v.5.0 */
/* <--v.4.0 */

      /* Store the individual energy terms */
      energy_term[0] = internal_energy;
      energy_term[1] = (float) (atom_clash_energy / 2.0);
      energy_term[2] = (float) (bond_length_energy / 2.0);
      energy_term[3] = (float) (bond_clash_energy / 2.0);
      energy_term[4] = (float) (bond_overlap_energy / 2.0);
/* v.4.0--> */
      energy_term[5] = boundary_energy;
      energy_term[6] = (float) (rel_pstn_energy / 2.0);
/* <--v.4.0 */

      /* Increment energy sums */
      for (i = 0; i < NTERMS; i++)
        {
          energy_sum[i] = energy_sum[i] + energy_term[i];
          total_energy = total_energy + energy_term[i];
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
      iobject++;
    }

  /* Print out summary of the total energy and its components */
  if (print_data == TRUE)
    {
/* v.4.0--> */
/*      printf("%5d.   ",loop); */
      printf("%4d. ",loop);
/* <--v.4.0 */
      for (i = 0; i < NTERMS; i++)
/* v.4.0--> */
/*      printf(" %8.2f",energy_sum[i]);
      printf("  %10.2f\n",total_energy); */
        printf("%8.2f",energy_sum[i]);
      printf(" %10.2f\n",total_energy);
/* <--v.4.0 */
/* v.5.0--> */
      fflush(stdout);
/* <--v.5.0 */
    }

/* v.4.3--> */
  /* If writing out to .res file, then do so */
  if (write_to_res == TRUE)
    {
      /* Write out just the required energy terms */
      fprintf(ligplot_res,"#Energy: ");
      for (i = 0; i < NTERMS; i++)
        fprintf(ligplot_res,"%8.2f",energy_sum[i]);
      fprintf(ligplot_res," %10.2f\n",total_energy);
    }
/* <--v.4.3 */

  /* Return the total energy */
  return(total_energy);
}
/* v.5.1--> */
/* v.5.2.1 --> */
/***********************************************************************

get_ligand_conformation  -  Use the current conformation number to flip
                            about the corresponding ligand bonds

***********************************************************************/

int get_ligand_conformation(struct object *object_ptr,int nbonds,
                            int start_bond,int end_bond,int iconf,
                            int calc_en)
{
  int clash, nflips;
  int both_ends, flip_all, ok, use_end, with_abort, with_fit;

  long bit_test, ibond;

  float internal_energy;

  struct bond *bond_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise bit-test pattern and bond number */
  bit_test = 1;
  ibond = 0;
  nflips = 0;
  ok = TRUE;
  with_abort = TRUE;
  with_fit = FALSE;
  both_ends = flip_all = FALSE;
  use_end = FIRST;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;
      
  /* Loop over the bonds, determining which to flip to
     obtain the current conformation */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* Check whether this bond is a rotatable one */
      if (bond_ptr->rotatable_bond == TRUE)
        {
          /* Only consider if bond is one of the ones we're
             interested in */
          if (ibond >= start_bond && ibond <= end_bond)
            {
              /* Determine whether this bond is to be flipped
                 or not */
              if (iconf & bit_test)
                {
                  /* Perform flip on this bond */
                  bond_flip(object_ptr,bond_ptr,use_end,both_ends,
                            flip_all,calc_en);
                  nflips++;
                }

              /* Double the bit-test parameter */
              bit_test = 2 * bit_test;
            }
              
          /* Increment count of rotatable bonds */
          ibond++;
        }

      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }

  /* Calculate this object's internal energy */
  internal_energy
    = calc_internal_energy(object_ptr,&clash,with_fit,with_abort);
  object_ptr->internal_energy = internal_energy;
  if (internal_energy > MAX_INTERNAL_ENERGY)
    ok = FALSE;

  /* Return whether this conformation is OK to use */
  return(ok);
}
/***********************************************************************

flip_ligand_bonds  -  Try all possible flips about ligand's rotatable
                      bonds to get the best match to the internal
                      anchor-point distances

***********************************************************************/

float flip_ligand_bonds(void)
{
  long best_conf, iconf, nbonds, nconf;

  int loop, nloops, start_bond, end_bond, nflip_bonds;
  int best_found, calc_en, no_move, ok, rel_anchors;

  float energy, best_energy, rmsd;
  float atom_clash_energy, bond_length_energy, bond_clash_energy;
  float bond_overlap_energy, boundary_energy, rel_pstn_energy;

  struct object *object_ptr;

  /* Initialise variables */
  energy = 0;
  no_move = FALSE;
  rel_anchors = TRUE;
  best_found = FALSE;
  calc_en = FALSE;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects, moving each one in turn */
  while (object_ptr != NULL)
    {
      /* If this is a ligand, then minimize */
      if (object_ptr->object_type == LIGAND)
        {
          /* Initialise variables */
          best_energy = energy = 0;
          best_conf = -1;
          nloops = 1;

          /* Get the number of rotatable bonds */
          nflip_bonds = nbonds = object_ptr->nrot_bonds;
          start_bond = 0;

          /* If there are too many rotatable bonds to try all,
             determine how many loops of subsets of the bonds we
             need to perform */
          if (nbonds > MAX_ROTATABLE)
            {
              /* Determine the number of loops to be done */
              nflip_bonds = MAX_ROTATABLE;
              nloops = (nbonds / (nflip_bonds - ROTBOND_OVERLAP)) + 1;
              if (nbonds % (nflip_bonds -ROTBOND_OVERLAP) == 0)
                nloops--;
            }

          /* Loop the required number of times */
          for (loop = 0; loop < nloops; loop++)
            {
              /* Initialise variables */
              best_energy = energy = 0;
              best_conf = -1;
              calc_en = FALSE;

              /* Save the starting coordinates of this objects */
              save_restore_coords(object_ptr,SAVE);

              /* Determine the end bond */
              end_bond = start_bond + nflip_bonds - 1;
              if (end_bond > nbonds - 1)
                {
                  end_bond = nbonds - 1;
                  nflip_bonds = end_bond - start_bond + 1;
                }

              /* Calculate the number of conformations we need to try */
              nconf = pow(2,nflip_bonds);

              
              /* Loop over all the conformations to be generated */
              for (iconf = 0; iconf < nconf; iconf++)
                {
                  /* Perform the bond-flips to get the current
                     conformation */
                  ok = get_ligand_conformation(object_ptr,nbonds,
                                               start_bond,end_bond,iconf,
                                               calc_en);

                  /* If the internal energy not too high, calculate total
                     energy, including the relative anchor positions
                     energy */
                  if (ok == TRUE)
                    {
                      /* Calculate the energy of this conformation */
                      energy
                        = calculate_energy(object_ptr,&atom_clash_energy,
                                           &bond_length_energy,
                                           &bond_clash_energy,
                                           &bond_overlap_energy,
                                           &boundary_energy,
                                           &rel_pstn_energy,
                                           rel_anchors);

                      /* If it is the best energy so far, then store
                         conformation */
                      if (iconf == 0 || energy < best_energy)
                        {
                          /* Store conformation and its energy as the
                             best so far */
                          best_energy = energy;
                          best_conf = iconf;
                        }
                    }

                  /* Restore ligand's starting conformation */
                  save_restore_coords(object_ptr,RESTORE);
                }

              /* Retrieve the lowest-energy conformation */
              if (best_conf > -1)
                {
                  calc_en = TRUE;
                  get_ligand_conformation(object_ptr,nbonds,start_bond,
                                          end_bond,best_conf,calc_en);
                }

              /* Set new start-bond for next loop (if there is one) */
              start_bond = start_bond + nflip_bonds - ROTBOND_OVERLAP;
            }

          /* Fit the ligand (minimised for relative anchor positions) to
             the absolute anchor positions */
          rmsd = fit_object(object_ptr,no_move);
          printf("FITTED with rmsd = %8.3f\n",rmsd);

          /* Calculate the energy of this conformation */
          energy
            = calculate_energy(object_ptr,&atom_clash_energy,
                               &bond_length_energy,
                               &bond_clash_energy,
                               &bond_overlap_energy,
                               &boundary_energy,
                               &rel_pstn_energy,
                               rel_anchors);
          printf("Total energy  %8.3f\n",energy);
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

/* v.5.3.4--> */
  /* Update all the residue boundaries */
  update_all_boundaries();
/* <--v.5.3.4 */

  /* Return whether a 'best' conformation was found */
  return(energy);
}
/***********************************************************************

old_get_ligand_conformation  -  Use the current conformation number to flip
                            about the corresponding ligand bonds

***********************************************************************/

void old_get_ligand_conformation(struct object *object_ptr,int nbonds,
                             int iconf)
{
  int bit_test, ibond;
  int both_ends, flip_all, use_end;

  struct bond *bond_ptr;
  struct object_bond *object_bond_ptr;

  /* Initialise bit-test pattern and bond number */
  bit_test = 1;
  ibond = 0;
  both_ends = flip_all = FALSE;
  use_end = FIRST;

  /* Get pointer to the first of this object's bonds */
  object_bond_ptr = object_ptr->first_object_bond_ptr;
      
  /* Loop over the bonds, determining which to flip to
     obtain the current conformation */
  while (object_bond_ptr != NULL)
    {
      /* Get the bond pointer */
      bond_ptr = object_bond_ptr->bond_ptr;

      /* Check whether this bond is a rotatable one */
      if (bond_ptr->rotatable_bond == TRUE)
        {
          /* Determine whether this bond is to be flipped
             or not */
          if (iconf & bit_test)
            {
              /* Perform flip on this bond */
              // OLD  bond_flip(object_ptr,bond_ptr,use_end,both_ends,
              // OLD            flip_all);
            }
              
          /* Increment count of rotatable bonds */
          ibond++;

          /* Double the bit-test parameter */
          bit_test = 2 * bit_test;
        }

      /* Get pointer to this object's next bond entry */
      object_bond_ptr = object_bond_ptr->next_object_bond_ptr;
    }
}
/***********************************************************************

old_flip_ligand_bonds  -  Try all possible flips about ligand's rotatable
                      bonds to get the best match to the internal
                      anchor-point distances

***********************************************************************/

void old_flip_ligand_bonds(void)
{
  int best_conf, iconf, nbonds, nconf;
  int rel_anchors;

  float energy, best_energy;
  float atom_clash_energy, bond_length_energy, bond_clash_energy;
  float bond_overlap_energy, boundary_energy, rel_pstn_energy;

  struct object *object_ptr;
  /* Initialise variables */
  rel_anchors = TRUE;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects, moving each one in turn */
  while (object_ptr != NULL)
    {
      /* If this is a ligand, then minimize */
      if (object_ptr->object_type == LIGAND)
        {
          /* Initialise variables */
          best_energy = energy = 0;
          best_conf = -1;

          /* Get the number of rotatable bonds */
          nbonds = object_ptr->nrot_bonds;

          /* Calculate the number of conformations we need to try */
          nconf = pow(2,nbonds);

          /* Save the coordinates of all objects */
          save_restore_coords(object_ptr,SAVE);

          /* Loop over all the conformations to be generated */
          for (iconf = 0; iconf < nconf && nbonds < MAX_ROTATABLE;
               iconf++)
            {
              /* Perform the bond-flips to get the current conformation */
              //OUT  get_ligand_conformation(object_ptr,nbonds,iconf);

              /* Calculate the energy of this conformation */
              energy
                = calculate_energy(object_ptr,&atom_clash_energy,
                                   &bond_length_energy,&bond_clash_energy,
                                   &bond_overlap_energy,
                                   &boundary_energy,&rel_pstn_energy,
                                   rel_anchors);

              /* If it is the best energy so far, then store
                 conformation */
              if (iconf == 0 || energy < best_energy)
                {
                  /* Store conformation and its energy as the best so
                     far */
                  best_energy = energy;
                  best_conf = iconf;
                }

              /* Restore ligand's starting conformation */
              save_restore_coords(object_ptr,RESTORE);
            }

          /* Retrieve the lowest-energy conformation */
          //OUT  if (best_conf > -1)
          //OUT    get_ligand_conformation(object_ptr,nbonds,best_conf);
        }

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }
}
/* <-- v.5.2.1 */
/* <--v.5.1 */
/***********************************************************************

energy_minimize  -  Minimize overall energy of the entire structure to
                    make numbers of atom and bond clashes a minimum

***********************************************************************/

/* v.4.3--> */
/* void energy_minimize(int *iseed,int nligands) */
/* v.5.4.4 --> */
//float energy_minimize(int *iseed,int nligands)
float energy_minimize(int *iseed,int nligands,int atom_counter)
/* <-- v.5.4.4 */
/* <--v.4.3 */
{
  int loop, nlines, object_type;
  int accept_move, flip_about_x, flip_all, got_move, move_type, swing;
  int rseed;
  int energy_limit_reached;
  int n_accept, naccept_move, ntrial_move, total_accept;
  int naccept_swivel, ntrial_swivel, total_swivel;
  int naccept_swing, ntrial_swing, total_swing;
  int naccept_flip, ntrial_flip;
/* v.5.2.1 --> */
  int naccept_atom_flip, ntrial_atom_flip;
  int tot_accept_move, tot_accept_atom_flip, tot_accept_flip;
  int tot_accept_swing, tot_accept_swivel, tot_accept;
  int tot_trial_move, tot_trial_atom_flip, tot_trial_flip;
  int tot_trial_swing, tot_trial_swivel;
  int tot_trials;
/* <-- v.5.2.1 */
  int nloops, nmove;
  int print_data, print_headings;
/* v.5.0--> */
  int no_flip, rel_anchors;
/* <--v.5.0 */

  float atom_clash_energy, bond_length_energy, bond_clash_energy;
  float bond_overlap_energy, internal_energy;
  float energy, previous_energy;
  float energy_drop, last_total_energy;
  float move_size, x;
  float max_swivel_angle;
  float max_swing_angle;
/* v.4.0--> */
  float boundary_energy;
  float rel_pstn_energy;
/* <--v.4.0 */
/* v.5.0--> */
  float accept, ediff, efactor, start_temp, temp, try;
/* <--v.5.0 */

  struct object *move_object_ptr;

  /* Initialise variables */
  rseed = *iseed;
  total_accept = 0;
  total_swing = 0;
  total_swivel = 0;
  nlines = 0;
  last_total_energy = 0.0;
  print_headings = TRUE;
  energy_limit_reached = FALSE;
/* v.5.0--> */
  rel_anchors = FALSE;
  no_flip = FALSE;
  temp = start_temp = 300;
/* <--v.5.0 */
/* v.5.0.1--> */
  ntrial_move = 0;
  ntrial_swivel = 0;
  ntrial_swing = 0;
  ntrial_flip = 0;
/* <--v.5.0.1 */
/* v.5.2.1 --> */
  naccept_atom_flip = 0;
  ntrial_atom_flip = 0;
  tot_accept = 0;
  tot_accept_move = 0;
  tot_accept_atom_flip = 0;
  tot_accept_flip = 0;
  tot_accept_swing = 0;
  tot_accept_swivel = 0;
  tot_trials = 0;
  tot_trial_move = 0;
  tot_trial_atom_flip = 0;
  tot_trial_flip = 0;
  tot_trial_swing = 0;
  tot_trial_swivel = 0;
/* <-- v.5.2.1 */

  /* Update all the residue boundaries */
  update_all_boundaries();

  /* Initialise minimization variables */
  loop = 0;
  max_swivel_angle = PI;
  max_swing_angle = PI;
  move_size = Max_HBdist;
  move_type = NMOVE_TYPES;
  n_accept = 0;
  nloops = Max_Loops;
/* v.5.0--> */
  try = accept = 0;
/* <--v.5.0 */

/* v.5.4.4 --> */
  /* If this is a very large structure, reduce number of loops to 
     prevent minimization taking for ever */
  if (atom_counter > LARGE_STRUCTURE)
    {
/* debug */
      efactor = (float) LARGE_STRUCTURE/ (float) atom_counter;
      nloops = efactor * efactor * Max_Loops;
/* debug */
    }
/* <-- v.5.4.4 */

  /* Show number of loops */
  printf("\n");
  printf("\n");
  printf("Energy minimization ... %d loops",nloops);
  if (Energy_Drop_Maximum > 0.0)
    printf(", terminating if energy-drop < %f\n",Energy_Drop_Maximum);
  else
    printf("\n");

  /* Loop until energy reaches an acceptable level, or minimization
     has gone on for the required number of loops */
  while (loop < nloops && energy_limit_reached == FALSE)
    {
      /* Initialise counts */
      naccept_move = 0;
      naccept_swivel = 0;
      naccept_swing = 0;
      naccept_flip = 0;
/* v.5.2.1 --> */
      naccept_atom_flip = 0;
/* <-- v.5.2.1 */
      if (loop % 4 == 0)
        n_accept = 0;
      nmove = 0;

/* v.5.0--> */
      /* If minimizing ligand on its own, determine fit-type to anchors */
      if (ligfit == TRUE)
        {
          /* In early stages, minimize relative anchor positions */
          if (loop < STAGE1)
/* v.5.2.1 --> */
//            rel_anchors = TRUE;
            {
              rel_anchors = FALSE;
              no_flip = TRUE;
            }
/* <-- v.5.2.1 */

          /* In medium stages, don't allow flips which might undo
              conformation of ligand */
          else if (loop < STAGE2)
            {
/* v.5.2.1 --> */
//              no_flip = TRUE;
              no_flip = FALSE;
/* <-- v.5.2.1 */
/* v.5.0.2--> */
              rel_anchors = FALSE;
/* <--v.5.0.2 */
            }

/* v.5.2.1 --> */
/* v.5.1--> */
//          /* Try stage 1 again */
//          else if (loop < STAGE3)
//            {
//              no_flip = FALSE;
//              rel_anchors = TRUE;
//            }
//
//          /* Go back to stage 2 */
//          else if (loop < STAGE4)
//            {
//              no_flip = TRUE;
//              rel_anchors = FALSE;
//            }
/* <--v.5.1 */
/* <-- v.5.2.1 */

          /* Subsequently, minimize absolute anchor positions */
          else
            {
              no_flip = FALSE;
              rel_anchors = FALSE;
            }

/* v.5.1--> */
          /* Define energy minimization flags */
/* v.5.2.1 --> */
//          rel_anchors = FALSE;
//          no_flip = TRUE;
/* <-- v.5.2.1 */
/* <--v.5.1 */
        }
/* <--v.5.0 */

      /* Determine the move_type */
      move_type++;
      if (move_type >= NMOVE_TYPES)
        move_type = 0;

      /* Initialise pointer to the first stored object */
      move_object_ptr = first_object_ptr;

      /* Loop through all objects, moving each one in turn */
      while (move_object_ptr != NULL)
        {
         /* Initialise flags */
          accept_move = FALSE;
          flip_all = FALSE;
          object_type = move_object_ptr->object_type;

          /* Make sure the move type is compatible with the object
             type */
          got_move = FALSE;
          if (object_type == LIGAND)
            {
/* v.5.0--> */
//              if (nligands == 1 && move_type == FLIP_OBJECT)
//                got_move = TRUE;
              if (ligfit == FALSE && nligands == 1 &&
                  move_type == FLIP_OBJECT)
                got_move = TRUE;
              if (ligfit == TRUE && nligands == 1)
                {
                  if (rel_anchors == TRUE)
                    {
                      if (move_type == FLIP_OBJECT)
                        got_move = TRUE;
                      else
                        got_move = FALSE;
                    }
                  else
                    {
                      if (no_flip == TRUE && (move_type == FLIP_OBJECT ||
                                              move_type == FLIP_ABOUT_ATOM))
                        got_move = FALSE;
                      else
                        got_move = TRUE;
                    }
                }
/* <--v.5.0 */
              else if (nligands > 1)
                got_move = TRUE;
            }
          else if (object_type == HYDROPHOBIC || object_type == WATER ||
                   object_type == SIMPLE_HGROUP)
            {
              if (move_type == MOVE_OBJECT)
                got_move = TRUE;
            }
          else if (object_type == HGROUP)
            got_move = TRUE;

          /* If flipping about a ligand bond, determine whether to take
             all connected group when performing the flip */
/* v.5.0--> */
//          if (object_type == LIGAND && move_type == FLIP_OBJECT)
          if (object_type == LIGAND && move_type == FLIP_OBJECT &&
              ligfit == FALSE)
/* <--v.5.0 */
            {
              /* Generate a random number and use to determine whether
                 to take all connected objects with the flipped part
                 of the ligand */
              x = get_random_number(&rseed,Random_Start);
              if (x < 0.5)
                flip_all = TRUE;
            }

          /* If flipping about a nonligand bond, decide whether to
             change flip to be about an atom instead */
          if (object_type != LIGAND && move_type == FLIP_OBJECT)
            {
              /* Generate a random number and determine whether to
                 flip about an atom instead */
              x = get_random_number(&rseed,Random_Start);
              if (x < 0.2)
                {
                  move_type = FLIP_ABOUT_ATOM;
                  flip_about_x = TRUE;
                  if (x < 0.1)
                    flip_about_x = FALSE;
                }
            }

          /* Switch off swings if near the end */
          if (loop > (nloops - 1) && move_type == SWING_OBJECT)
            got_move = FALSE;

/* v.5.0--> */
          /* If object's coords were fixed to anchor positions, then
             disallow flips about bonds and atoms */
          if (move_object_ptr->fixed == TRUE &&
              (move_type == FLIP_OBJECT || move_type == FLIP_ABOUT_ATOM))
            got_move = FALSE;
/* <--v.5.0 */

/* v.5.3.4--> */
          /* If all object's coords fixed, then disallow any type of
             movement */
          if (move_object_ptr->coords_fixed == TRUE)
            got_move = FALSE;

          /* If fitting to an initial DIMHTML layout, don't allow moves */
          if (dimfit == TRUE &&
              (move_type == MOVE_OBJECT || move_type == SWING_OBJECT))
            got_move = FALSE;

          /* If this is a 1-point object and it has been anchored, then
             don't need to move it at all */
          if (move_object_ptr->fixed == TRUE &&
              (object_type == HYDROPHOBIC || object_type == WATER ||
               object_type == SIMPLE_HGROUP))
            got_move = FALSE;
/* <--v.5.3.4 */

          /* If have a valid move for this object, then do it */
          if (got_move == TRUE)
            {
              /* Increment count of valid moves on this loop */
              nmove++;
/* v.5.2.1 --> */
              if (move_type == MOVE_OBJECT)
                tot_trial_move++;
              else if (move_type == SWIVEL_OBJECT)
                tot_trial_swivel++;
              else if (move_type == SWING_OBJECT)
                tot_trial_swing++;
              else if (move_type == FLIP_OBJECT)
                tot_trial_flip++;
              else if (move_type == FLIP_ABOUT_ATOM)
                tot_trial_atom_flip++;
              tot_trials++;
/* <-- v.5.2.1 */

              /* Get the object's internal energy */
              internal_energy = move_object_ptr->internal_energy;

              /* If move might affect more than one object, compute
                 the total energy */
              if (flip_all == TRUE)
                {
                  print_headings = FALSE;
                  print_data = FALSE;
                  previous_energy
                    = get_total_energy(loop,print_headings,
/* v.4.3--> */
/*                                                   print_data); */
/* v.5.0--> */
/*                                                     print_data,FALSE); */
                                       print_data,FALSE,rel_anchors);
/* <--v.5.0 */
/* <--v.4.3 */

                  /* Save the coordinates of all objects */
                  save_restore_all_objects(SAVE);
                }

              /* Otherwise, just calculate the energy between this
                 object and all the others */
              else
                {
                  /* Calculate this object's energy */
                  previous_energy = calculate_energy(move_object_ptr,
                                                     &atom_clash_energy,
                                                     &bond_length_energy,
                                                     &bond_clash_energy,
/* v.4.0--> */
/*                                                   &bond_overlap_energy); */
                                                     &bond_overlap_energy,
                                                     &boundary_energy,
/* v.5.0--> */
//                                                     &rel_pstn_energy);
                                                     &rel_pstn_energy,
                                                     rel_anchors);
/* <--v.5.0 */
/* <--v.4.0 */

                  /* Save its coordinates */
                  save_restore_coords(move_object_ptr,SAVE);
                }                 

              /* If moving object, make a random move to the current
                 object */
              if (move_type == MOVE_OBJECT)
                {
                  /* Make the move */
                  make_random_move(move_object_ptr,move_size,&rseed);

                  /* Increment count of trial moves */
                  ntrial_move++;
                }

              /* Or, make a random rotation or swing of the
                 current object */
              else if (move_type == SWIVEL_OBJECT ||
                       move_type == SWING_OBJECT)
                {
                  /* Set flag about whether this is a swing or
                     a rotation */
                  swing = FALSE;
                  if (move_type == SWING_OBJECT)
                    swing = TRUE;

                  /* Apply a trial rotation/swing */
                  make_random_rotation(move_object_ptr,max_swivel_angle,
                                       &rseed,swing);

                  /* Increment count of trial rotations */
                  if (move_type == SWIVEL_OBJECT)
                    ntrial_swivel++;
                  else
                    ntrial_swing++;
                }

              /* Or, make a random flip about one of the object's
                 rotatable bonds */
              else if (move_type == FLIP_OBJECT)
                {
                  /* Make the flip */
                  make_random_flip(move_object_ptr,&rseed,flip_all);

                  /* Increment count of trial flips */
                  ntrial_flip++;
                }

              /* Or, make a random flip about one of the object's
                 atoms */
              else if (move_type == FLIP_ABOUT_ATOM)
                {
                  /* Make the flip */
                  make_random_atom_flip(move_object_ptr,&rseed,
                                        flip_about_x);

                  /* Increment count of trial flips */
/* v.5.2.1 --> */
//                  ntrial_flip++;
                  ntrial_atom_flip++;
/* <-- v.5.2.1 */
                }

              /* If move affected more than one object, compute the new
                 total energy */
              if (flip_all == TRUE)
                {
                  /* Update all residue- and object-boundaries */
                  update_all_boundaries();

                  /* Get total energy */
                  print_headings = FALSE;
                  print_data = FALSE;
/* v.4.3--> */
/*                energy = get_total_energy(loop,print_headings,print_data); */
                  energy = get_total_energy(loop,print_headings,print_data,
/* v.5.0--> */
/*                                            FALSE); */
                                            FALSE,rel_anchors);
/* <--v.5.0 */
/* <--v.4.3 */
                }

              /* Otherwise, calculate new energy between object and all
                 the others */
              else
                {
                  /* Update residue- and object-boundaries */
                  update_residue_boundaries(move_object_ptr);

                  /* Calculate energy */
                  energy = calculate_energy(move_object_ptr,
                                            &atom_clash_energy,
                                            &bond_length_energy,
                                            &bond_clash_energy,
/* v.4.0--> */
/*                                          &bond_overlap_energy); */
                                            &bond_overlap_energy,
                                            &boundary_energy,
/* v.5.0--> */
//                                            &rel_pstn_energy);
                                            &rel_pstn_energy,rel_anchors);
/* <--v.5.0 */
/* <--v.4.0 */
                }

              /* Adjust annealing temperature */
              temp
                = start_temp * (float) (nloops - loop)
                / (2 * (float) (nloops));

              /* Determine whether to accept or reject this move */
/* v.5.0--> */
//              if (energy < previous_energy)
//                accept_move = TRUE;
              if (ligfit == TRUE && energy > previous_energy &&
                  loop < nloops - 100 && temp > 0)
                {
                  /* Calculate parameters for choosing whether to accept
                     increase in energy */
                  try = get_random_number(&rseed,Random_Start);
                  ediff = energy - previous_energy;
                  efactor = - ediff / (BOLTZMANN_CONST * temp);
                  accept = exp(efactor);
                }
              else
                try = accept = 0;
              if (energy < previous_energy ||
                  (energy > previous_energy && try < accept))
                accept_move = TRUE;
/* <--v.5.0 */

              /* Update total energy, depending on whether move has
                 been accepted or not */
              if (accept_move == TRUE)
                {
                  /* Increment count of accepted moves */
/* v.5.2.1 --> */
//                  if (move_type == MOVE_OBJECT)
//                    naccept_move++;
//                  else if (move_type == SWIVEL_OBJECT)
//                    naccept_swivel++;
//                  else if (move_type == SWING_OBJECT)
//                    naccept_swing++;
//                  else if (move_type == FLIP_OBJECT ||
//                           move_type == FLIP_ABOUT_ATOM)
//                    naccept_flip++;
                  if (move_type == MOVE_OBJECT)
                    {
                      naccept_move++;
                      tot_accept_move++;
                    }
                  else if (move_type == SWIVEL_OBJECT)
                    {
                      naccept_swivel++;
                      tot_accept_swivel++;
                    }
                  else if (move_type == SWING_OBJECT)
                    {
                      naccept_swing++;
                      tot_accept_swing++;
                    }
                  else if (move_type == FLIP_OBJECT)
                    {
                      naccept_flip++;
                      tot_accept_flip++;
                    }
                  else if (move_type == FLIP_ABOUT_ATOM)
                    {
                      naccept_atom_flip++;
                      tot_accept_atom_flip++;
                    }

                  /* Increment total of all accepted moves */
                  tot_accept++;
/* <-- v.5.2.1 */
                }

              /* Otherwise, if no move, revert to previous conformation */
              else
                {
                  /* If several objects might have been moved, restore
                     the coordinates of all the objects */
                  if (flip_all == TRUE)
                    {
                      /* Restore the coords */
                      save_restore_all_objects(RESTORE);

                      /* Update all residue- and object-boundaries */
                      update_all_boundaries();
                    }
                  
                  /* Otherwise, restore the moved object only */
                  else
                    {
                      /* Restore the coords */
                      save_restore_coords(move_object_ptr,RESTORE);

                      /* Update residue- and object-boundaries */
                      update_residue_boundaries(move_object_ptr);
                    }

                  /* Update object's internal energy */
                  move_object_ptr->internal_energy = internal_energy;
                }
            }

          /* Get pointer to next object */
          move_object_ptr = move_object_ptr->next_object_ptr;
        }

      /* Increment loop count */
      loop++;

      /* Count numbers of trial moves accepted */
      n_accept = naccept_move + naccept_swivel + naccept_swing
        + naccept_flip;
      total_accept = total_accept + n_accept;
      total_swing = total_swing + naccept_swing;
      total_swivel = total_swivel + naccept_swivel;
          
      /* Show progress of minimization */
      if (loop % 20 == 0 || loop == nloops)
        {
          if (nlines % 20 == 0)
            print_headings = TRUE;
          else
            print_headings = FALSE;
          print_data = TRUE;
/* v.4.3--> */
/*        energy = get_total_energy(loop,print_headings,print_data); */
          energy
/* v.5.0--> */
/*            = get_total_energy(loop,print_headings,print_data,FALSE); */
            = get_total_energy(loop,print_headings,print_data,FALSE,
                               rel_anchors);
/* <--v.5.0 */
/* <--v.4.3 */

          /* If this is the first loop, save the current total energy */
          if (nlines == 0)
            last_total_energy = energy;

          /* Otherwise, check whether the energy drop since last time is
             smaller than the stop-limit */
          else
            {
              /* Get the drop in energy since the last loop */
              energy_drop = last_total_energy - energy;

              /* If this drop is lower than the user-defined limit, then
                 set flag to stop minimization process */
              if (energy_drop > 0.0 && energy_drop < Energy_Drop_Maximum)
                energy_limit_reached = TRUE;
            }

          /* Save the current total energy */
          last_total_energy = energy;
          nlines++;
        }

      /* Adjust move size, if necessary */
      if (ntrial_move > 0)
        {
          /* accept_percen = 100.0 * (float) naccept_move
            / (float) ntrial_move;
          if (accept_percen < 50.0)
            move_size = 0.9 * move_size;  */

            /* Check that move distance doesn't fall below minimum value */
            if (move_size < MIN_MOVE)
              move_size = MIN_MOVE;
        }

      /* Adjust rotation angle, if necessary */
      if (ntrial_swivel > 0)
        {
          /* accept_percen
            = 100.0 * (float) naccept_swivel / (float) ntrial_swivel;
          if (accept_percen < 50.0)
            max_swivel_angle = 0.9 * max_swivel_angle;  */
            
            /* Check that angle doesn't fall below minimum value */
            if (max_swivel_angle * RADDEG < MIN_ANGLE)
              max_swivel_angle = MIN_ANGLE / RADDEG;
        }

      /* Adjust swing angle, if necessary */
      if (ntrial_swing > 0)
        {
          /* accept_percen
            = 100.0 * (float) naccept_swing / (float) ntrial_swing;
          if (accept_percen < 50.0)
            max_swing_angle = 0.9 * max_swing_angle;  */

            /* Check that angle doesn't fall below minimum value */
            if (max_swing_angle * RADDEG < MIN_ANGLE)
              max_swing_angle = MIN_ANGLE / RADDEG;
        }
    }

  /* Retrieve the current random-number seed */
  *iseed = rseed;

/* v.4.3--> */
  /* If writing a .res file, then write out final energies to that file */
  if (Write_Res_File == TRUE)
    {
      print_headings = FALSE;
      print_data = FALSE;
      energy
/* v.5.0--> */
/*         = get_total_energy(loop,print_headings,print_data,TRUE); */
        = get_total_energy(loop,print_headings,print_data,TRUE,
                           rel_anchors);
/* <--v.5.0 */
    }

/* v.5.2.1 --> */
  /* Show minimization statistics */
  printf("Minimization statistics:\n");
  printf("    Trial moves:      %8d     Accepted %8d\n",tot_trial_move,
         tot_accept_move);
  printf("    Trial swivels:    %8d     Accepted %8d\n",tot_trial_swivel,
         tot_accept_swivel);
  printf("    Trial swings:     %8d     Accepted %8d\n",tot_trial_swing,
         tot_accept_swing);
  printf("    Trial flips:      %8d     Accepted %8d\n",tot_trial_flip,
         tot_accept_flip);
  printf("    Trial atom flips: %8d     Accepted %8d\n",tot_trial_atom_flip,
         tot_accept_atom_flip);
  printf("\n");
  printf("    Total trials:     %8d     Accepted %8d\n",tot_trials,
         tot_accept);
/* <-- v.5.2.1 */

  /* Return the final energy */
  return(energy);
/* <--v.4.3 */
}
/* v.5.3.4--> */
/***********************************************************************
 
sort_iface_objects - Sort all objects in the current interface according
                     to mean y-value

***********************************************************************/

void sort_iface_objects(struct object **store_object_ptr,int nobject)
{
  int ipos;
  int swapped;

  struct object *object_ptr, *other_object_ptr;

  /* Initialise variables */
  swapped = TRUE;

  /* Loop until all the values have been sorted */
  while (swapped == TRUE)
    {
      /* Initialise flag */
      swapped = FALSE;

      /* Loop over all the objects */
      for (ipos = 0; ipos < nobject - 1; ipos++)
        {
          /* Get this object and the next one */
          object_ptr = store_object_ptr[ipos];
          other_object_ptr = store_object_ptr[ipos + 1];
          
          /* If objects are the wrong way round, then swap */
          if ((object_ptr->miny + object_ptr->maxy) >
              (other_object_ptr->miny + other_object_ptr->maxy))
            {
              /* Swap the two objects */
              store_object_ptr[ipos] = other_object_ptr;
              store_object_ptr[ipos + 1] = object_ptr;

              /* Set flag to indicate that a swap took place in
                 this loop */
              swapped = TRUE;
            }
        }
    }
}
/***********************************************************************
 
pair_up_objects  -  Pair up objects on opposite sides of the interface,
                     including any waters that might be between them

***********************************************************************/

int pair_up_objects(struct object *store_object_ptr[3][MAXATOMS],
                    int nobject[3],struct object ***obj_pair_ptr,
                    int max_pairs)
{
  int i, iface, ipos, next[3], npairs;
  int done, pair_up;

  float miny, overlap;

  struct object *object_ptr[3], **object_pair_ptr;

  /* Initialise variables */
  done = FALSE;
  next[0] = next[1] = next[2] = 0;
  npairs = 0;

  /* Create the object_pair_ptr array */
  object_pair_ptr
    = (struct object **) malloc(3 * max_pairs * sizeof(struct object *));
  if (object_pair_ptr == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for object_pair_ptr\n");
      exit(1);
    }

  /* Initialise the object_pair_ptr */
  for (i = 0; i < 3 * max_pairs; i++)
    object_pair_ptr[i] = NULL;

  /* Loop through the sorted lists of objects to pair up */
  while (done == FALSE)
    {
      /* Identify the object with the lowest y-values */
      iface = -1;
      miny = 0;

      /* Get the next potential pair */
      for (i = 0; i < 3; i++)
        {
          /* If valid, get this object */
          if (next[i] < nobject[i])
            {
              /* Get the object */
              object_ptr[i] = store_object_ptr[i][next[i]];

              /* If this is the first, or has a lower y-value,
                 store the y-value */
              if (iface == -1 ||
                  (object_ptr[i]->miny - object_ptr[i]->max_atom_size)
                  < miny)
                {
                  miny = object_ptr[i]->miny - object_ptr[i]->max_atom_size;
                  iface = i;
                }
            }
          else
            object_ptr[i] = NULL;
        }

      /* If have a seed object, see if the other two can be paired up
         with it */
      if (iface > -1)
        {

          /* Store the seed object */
          ipos = 3 * npairs + iface;
          object_pair_ptr[ipos] = object_ptr[iface];

          /* Increment this object's subscript */
          next[iface]++;

          /* Check to see which of the other two objects can be paired
             with this one because they overlap significantly in the
             y-direction */
          for (i = 0; i < 3; i++)
            {
              /* Process if valid and not the seed object */
              if (object_ptr[i] != NULL && i != iface)
                {
                  /* Initialise flag */
                  pair_up = FALSE;

                  /* If object overlaps with seed, decide whether to pair
                     it */
                  if ((object_ptr[i]->miny
                       - object_ptr[i]->max_atom_size)
                      < (object_ptr[iface]->maxy
                         + object_ptr[iface]->max_atom_size))
                    {
                      /* If wholly inside, then accept */
                      if ((object_ptr[i]->maxy
                           + object_ptr[i]->max_atom_size)
                          < (object_ptr[iface]->maxy
                             + object_ptr[iface]->max_atom_size))
                        pair_up = TRUE;

                      /* Otherwise, determine the extent of the overlap */
                      else
                        {
                          /* Calculate overlap and determine whether to
                             keep */
                          overlap
                            = object_ptr[iface]->maxy
                            + object_ptr[iface]->max_atom_size
                            - object_ptr[i]->miny
                            + object_ptr[i]->max_atom_size;

                          /* If overlap sufficient, then pair up */
                          if (overlap > MIN_OVERLAP)
                            pair_up = TRUE;
                        }
                    }

                  /* If object to be paired up, store it */
                  if (pair_up == TRUE)
                    {
                      /* Store this object */
                      ipos = 3 * npairs + i;
                      object_pair_ptr[ipos] = object_ptr[i];

                      /* Increment this object's subscript */
                      next[i]++;
                    }
                }
            }

          /* Increment position in pairs list */
          npairs++;
          if (npairs > max_pairs)
            done = TRUE;
        }

      /* Otherwise, we are done */
      else
        done = TRUE;
    }

  /* Return the array */
  *obj_pair_ptr = object_pair_ptr;

  /* Return the number of pairs stored */
  return(npairs);
}
/***********************************************************************
 
pack_paired_objects  -  Align and pack the pairs as best as possible

***********************************************************************/

void pack_paired_objects(struct object **object_pair_ptr,int npairs)
{
  int i, iatom, ipair, ipos, iresid, natoms, nobj, nresid;
  int first;

  float adjust[3], dist, lasty, mindist, miny, maxy, newmin[3], shift;

  struct coordinate *atom_ptr;
  struct object *object_ptr[3], *last_object_ptr[3];
  struct residue *residue_ptr;

  /* Initialise variables */
  lasty = 0;
  for (i = 0; i < 3; i++)
    last_object_ptr[i] = NULL;

  /* Loop over the paired objects */
  for (ipair = 0; ipair < npairs; ipair++)
    {
      /* Initialise variables */
      nobj = 0;
      miny = maxy = 0;
      shift = 0;

      /* Store the three objects at this pair position */
      for (i = 0; i < 3; i++)
        {
          /* Initialise the adjust */
          adjust[i] = newmin[i] = 0;

          /* Get this position */
          ipos = 3 * ipair + i;

          /* Get this object */
          object_ptr[i] = object_pair_ptr[ipos];
        }

      /* Loop over the objects in this pair group */
      for (i = 0; i < 3; i++)
        {
          /* If valid, update the maximum and minimum y-coords */
          if (object_ptr[i] != NULL)
            {
              /* If the first, then store the max and min y-coords */
              if (nobj == 0)
                {
                  miny
                    = object_ptr[i]->miny - object_ptr[i]->max_atom_size;
                  maxy
                    = object_ptr[i]->maxy + object_ptr[i]->max_atom_size;
                }

              /* Otherwise update limits */
              else
                {
                  if ((object_ptr[i]->miny - object_ptr[i]->max_atom_size)
                      < miny)
                    miny
                      = object_ptr[i]->miny - object_ptr[i]->max_atom_size;
                  if ((object_ptr[i]->maxy + object_ptr[i]->max_atom_size)
                      > maxy)
                    maxy
                      = object_ptr[i]->maxy + object_ptr[i]->max_atom_size;
                }

              /* Increment object count */
              nobj++;
            }
        }

      /* Compute the shift that can be made to each object */
      for (i = 0; i < 3; i++)
        {
          /* If valid, centre and shift this object's coords */
          if (object_ptr[i] != NULL)
            {
              /* Calculate the shift to centre this object */
              adjust[i] = lasty - miny + (miny + maxy) / 2
                - (object_ptr[i]->miny + object_ptr[i]->maxy) / 2;

              /* Calculate this object's new minimum y-value */
              newmin[i] = object_ptr[i]->miny
                - object_ptr[i]->max_atom_size + adjust[i];
            }
        }

      /* See if we can pack this pair a little closer to the previous
         pair */
      if (last_object_ptr[0] != NULL || last_object_ptr[1] != NULL)
        {
          /* Initialise the shift */
          mindist = -1;
          shift = 0;

          /* Calculate the distance from the current new minimum to
             the last object's maximum */
          if (last_object_ptr[0] != NULL)
            mindist = newmin[0] - last_object_ptr[0]->maxy
              - last_object_ptr[0]->max_atom_size;
          if (last_object_ptr[1] != NULL)
            {
              dist = newmin[1] - last_object_ptr[1]->maxy
                - last_object_ptr[1]->max_atom_size;
              if (mindist < 0 || dist < mindist)
                mindist = dist;
            }

          /* Subtract the ideal separation */
          mindist = mindist - IDEAL_SEPARATION;

          /* If have a positive value left, then it means we can bunch
             the current pair a little closer */
          if (mindist > 0)
            shift = - mindist;
        }

      /* Use the maximum and minimum values to centre-up the residues
         in each pair and shift to the next spaced-out position */
      for (i = 0; i < 3; i++)
        {
          /* If valid, centre and shift this object's coords */
          if (object_ptr[i] != NULL)
            {
              /* Get the first of this object's residues */
              residue_ptr = object_ptr[i]->first_residue_ptr;
              nresid = object_ptr[i]->nresidues;
              iresid = 0;

              /* Loop over all this object's residues */
              while (iresid < nresid && residue_ptr != NULL)
                {
                  /* Get pointer to first residue's first atom */
                  atom_ptr = residue_ptr->first_atom_ptr;
                  iatom = 0;
                  natoms = residue_ptr->natoms;

                  /* Loop over all this residue's atoms */
                  while (iatom < natoms && atom_ptr != NULL)
                    {
                      /* Adjust this atom's y-coordinate */
                      atom_ptr->y = atom_ptr->y + adjust[i] + shift;

                      /* Get pointer to the next atom */
                      atom_ptr = atom_ptr->next;
                      iatom++;
                    }

                  /* Get pointer to the next residue */
                  residue_ptr = residue_ptr->next_residue_ptr;
                  iresid++;
                }
            }
        }

      /* Initialise the last y-value */
      lasty = 0;
      first = TRUE;

      /* Update and save the current three objects */
      for (i = 0; i < 3; i++)
        {
          /* Update object's boundaries */
          if (object_ptr[i] != NULL)
            {
              /* Update the boundaries */
              update_residue_boundaries(object_ptr[i]);

              /* See if this has the highest y-value so far */
              if (first == TRUE ||
                  (object_ptr[i]->maxy + object_ptr[i]->max_atom_size)
                   > lasty)
                {
                  lasty = object_ptr[i]->maxy
                    + object_ptr[i]->max_atom_size;
                  first = FALSE;
                }
            }
          last_object_ptr[i] = object_ptr[i];
        }

      /* Update the current last y-position */
      lasty = lasty + IDEAL_SEPARATION;
    }
}
/***********************************************************************
 
place_interface_objects  -  Place the objects so that they are
                            equidistant from the interface line

***********************************************************************/

void place_interface_objects(struct object *store_object_ptr[3][MAXATOMS],
                             int nobject[3],int *face)
{
  int iatom, iface, iobject, iresid;
  int natoms, nresid;

  float adjust;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Loop over the two interfaces */
  for (iface = 0; iface < 3; iface++)
    {
      /* Loop over the sorted list of objects in this loop */
      for (iobject = 0; iobject < nobject[iface]; iobject++)
        {
          /* Get this object */
          object_ptr = store_object_ptr[iface][iobject];

          /* If object is the left-hand interface, adjust to its max
             x-coord */
          if (face[iface] == LEFT_FACE)
            adjust = - object_ptr->maxx - IDEAL_GAP / 2
              - object_ptr->max_atom_size;

          /* Otherwise, adjust to its min x-coord */
          else if (face[iface] == RIGHT_FACE)
            adjust = - object_ptr->minx + IDEAL_GAP / 2
              + object_ptr->max_atom_size;

          /* Otherwise, for water, place on the interface line */
          else
            adjust = - (object_ptr->minx + object_ptr->maxx) / 2;

          /* Set the pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get pointer to first residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all this residue's atoms */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Adjust this atom's x-coordinate */
                  atom_ptr->x = atom_ptr->x + adjust;

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
            }

          /* Update this object's maximum and minimum coords */
          update_residue_boundaries(object_ptr);
        }
    }
}
/***********************************************************************

tighten_dimplot  -  Tighten up the plot by shuffling residues closer
                    together

***********************************************************************/

void tighten_dimplot(void)
{

  int face[3], i, iface;
  int max_pairs, nobject[3], npairs;
  int have_waters;

  float minx[3], maxx[3];

  struct object *object_ptr;
  struct object *store_object_ptr[3][MAXATOMS];

  struct object **object_pair_ptr;

  /* Initialise variables */
  have_waters = FALSE;
  minx[0] = minx[1] = maxx[0] = maxx[1] = 0;
  nobject[0] = nobject[1] = nobject[2] = 0;

  /* Update all residue- and object-boundaries */
  update_all_boundaries();

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Store interface and water objects in appropriate list */
      if (object_ptr->interface == 1 || object_ptr->interface == 2)
        iface = object_ptr->interface - 1;

      /* For waters, set flag */
      else
        {
          iface = 2;
          have_waters = TRUE;
        }
 
      /* Store this object's pointer */
      store_object_ptr[iface][nobject[iface]] = object_ptr;

      /* Update the min and max x values of each interface */
      if (nobject[iface] == 0)
        {
          /* Store current values as the max and min values */
          minx[iface] = object_ptr->minx;
          maxx[iface] = object_ptr->maxx;
        }

      /* Otherwise, update max and min values */
      else
        {
          if (object_ptr->minx < minx[iface])
            minx[iface] = object_ptr->minx;
          if (object_ptr->maxx > maxx[iface])
            maxx[iface] = object_ptr->maxx;
        }
          
      /* Increment object count to store their mean y coords */
      nobject[iface]++;

      /* Get pointer to next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Determine which interface is on the left, and which on the right */
  if (minx[0] < minx[1])
    {
      face[0] = LEFT_FACE;
      face[1] = RIGHT_FACE;
    }
  else
    {
      face[0] = RIGHT_FACE;
      face[1] = LEFT_FACE;
    }

  /* Get the maximum number of pairs we can form */
  max_pairs = 0;
  for (i = 0; i < 3; i++)
    max_pairs = max_pairs + nobject[i];

  /* If don't have any objects, then abort */
  if (max_pairs == 0)
    {
      printf("*** ERROR. No residues in the interface!\n");

/* v.5.4.4 --> */
      /* Write message out to the ligplot.drw file */
      fprintf(ligplot_drw,
              "ERROR: No residues in the interface!\n");

      /* Close the .drw file */
      fclose(ligplot_drw);
/* <-- v.5.4.4 */
      exit(1);
    }

  /* Sort the objects in this interface */
  for (iface = 0; iface < 3; iface++)
    {
      if (nobject[iface] > 1)
        sort_iface_objects(store_object_ptr[iface],nobject[iface]);
    }

  /* Pair up residues on opposite sides of the interface, including
     any waters that might be between them */
  npairs = pair_up_objects(store_object_ptr,nobject,&object_pair_ptr,
                           max_pairs);

  /* Align and pack the pairs as best as possible */
  pack_paired_objects(object_pair_ptr,npairs);

  /* Update all residue- and object-boundaries */
  update_all_boundaries();

  /* Place the objects so that they are equidistant from the interface
     line */
  place_interface_objects(store_object_ptr,nobject,face);
}
/* <--v.5.3.4 */
/***********************************************************************

rotate_whole_structure  -  Apply the given rotation to all the atoms

***********************************************************************/

void rotate_whole_structure(float v[3][3])
{
  int icoord;
  float x, y, z;
  float new_coords[3];

  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms to initialise the flag indicating whether
     coords have already been stored */
  while (atom_ptr != NULL)
    {
      /* Extract this atom's coordinates */
      x = atom_ptr->x;
      y = atom_ptr->y;
      z = atom_ptr->z;
        
      /* Apply the transformation */
      for (icoord = 0; icoord < 3; icoord++)
        {
          /* Calculate transformed coordinates */
          new_coords[icoord]
            =  ((x * v[0][icoord])
                + (y * v[1][icoord])
                + (z * v[2][icoord]));
        }

      /* Save the new coordinates */
      atom_ptr->x = new_coords[0];
      atom_ptr->y = new_coords[1];
      atom_ptr->z = new_coords[2];

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

swap_x_and_y  -  Swap the x- and y- coordinates

***********************************************************************/

void swap_x_and_y(void)
{
  float swap;

  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms */
  while (atom_ptr != NULL)
    {
      /* Swap the x- and y- coordinates */
      swap = atom_ptr->x;
      atom_ptr->x = atom_ptr->y;
      atom_ptr->y = swap;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

flip_about_x  -  Flip the whole structure about the x-axis

***********************************************************************/

void flip_about_x(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms */
  while (atom_ptr != NULL)
    {
      /* Perform the flip */
      atom_ptr->y = - atom_ptr->y;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

flip_about_y  -  Flip the whole structure about the y-axis

***********************************************************************/

void flip_about_y(void)
{
  struct coordinate *atom_ptr;

  /* Get pointer to the start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms */
  while (atom_ptr != NULL)
    {
      /* Perform the flip */
      atom_ptr->x = - atom_ptr->x;

      /* Get pointer to the next atom */
      atom_ptr = atom_ptr->next;
    }
}
/***********************************************************************

check_ligand_direction  -  Test if plot needs to be flipped to maintain
                           residue order down the page, or left-to-right
                           according to the orientation

***********************************************************************/

void check_ligand_direction(void)
{
  int iatom, iobject, iresid, natoms, ngot, nresid;
  int last_resno, resno;

  char resno_string[5];

  float C_x, C_y, N_x, N_y;

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Initialise maximum and minimum coordinates */

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;
  ngot = 0;
  last_resno = -999;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Process only if this is the ligand */
      if (object_ptr->object_type == LIGAND)
        {
          /* Set the pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* If this is the very first atom, store its
                     coordinates as the N-terminus coords */
                  if (ngot == 0)
                    {
                      N_x = atom_ptr->x;
                      N_y = atom_ptr->y;

                      /* Get the residue number */
                      strncpy(resno_string,atom_ptr->residue_ptr->res_num,4);
                      resno_string[4] = '\0';
                      resno = atoi(resno_string);
                    }

                  /* Otherwise, store as the C-terminal coords */
                  else
                    {
                      if (resno > last_resno)
                        {
                          C_x = atom_ptr->x;
                          C_y = atom_ptr->y;
                        }
                    }

                  /* Increment number of atoms processed */
                  ngot++;
                  if (resno > last_resno)
                    last_resno = resno;

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
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

  /* Determine whether whole picture needs to be flipped */

  /* If portrait plot needs to be flipped, the flip about y */
  if (Picture->Portrait == TRUE && N_y < C_y)
    flip_about_x();

  /* If landscape plot needs to be flipped, the flip about x */
  else if (Picture->Portrait == FALSE && N_x > C_x)
    flip_about_y();

}
/***********************************************************************

ligand_orientate  -  Orientate the ligand on the y-axis on the output
                     page

***********************************************************************/

void ligand_orientate(void)
{
  int iatom, iobject, iresid, natoms, nresid, nstored;
  int store_ca, store_mainchain, wanted;

  float transformation[3][3];
  float mass_x, mass_y, mass_z, coord_store[MAXATOMS][3];

  struct coordinate *atom_ptr;
  struct object *object_ptr;
  struct residue *residue_ptr;


  /* Initialise variables */
  nstored = 0;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;
  iobject = 0;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* Process only if this is the ligand */
      if (object_ptr->object_type == LIGAND)
        {
          /* Determine which atoms are to be stored */
          store_ca = FALSE;
          store_mainchain = FALSE;

          /* For the simplified plot, then use CA, providing have at
             least 2 of them */
          if (Include->Simple_Ligand_Residues && object_ptr->n_ca > 1)
            store_ca = TRUE;

          /* Otherwise, use just the mainchain atoms, if have a
             sufficient number */
          else if (object_ptr->nmain_chain > 9)
            store_mainchain = TRUE;

          /* Get pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              
              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Initialise flag */
                  wanted = FALSE;

                  /* Determine if atom is of the type required */
                  if (store_ca == TRUE)
                    {
                      /* If this is a CA, then want it */
                      if (!strncmp(atom_ptr->atom_name," CA ",4))
                        wanted = TRUE;
                    }
                  else if (store_mainchain == TRUE)
                    {
                      /* If this is a mainchain, then want it */
                      if (atom_ptr->side_chain == FALSE)
                        wanted = TRUE;
                    }
                  else
                    wanted = TRUE;

                  /* If atom is wanted, then store its coordinates */
                  if (wanted == TRUE)
                    {
                      /* Store the coords */
                      coord_store[nstored][0] = atom_ptr->x;
                      coord_store[nstored][1] = atom_ptr->y;
                      coord_store[nstored][2] = atom_ptr->z;
                      nstored++;
                    }

                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
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

  /* If have stored any coordinates, perform a principal components
     analysis */
  if (nstored > 1)
    {
      /* If plotting the simplified representation, then all residues
         will be in a straight line, so just use the first and last
         positions */
      if (store_ca == TRUE)
        {
          /* Move coords of last CA position to the second position */
          coord_store[1][0] = coord_store[nstored - 1][0];
          coord_store[1][1] = coord_store[nstored - 1][1];
          coord_store[1][2] = coord_store[nstored - 1][2];

          /* Set the number of points to 2 */
          nstored = 2;
        }

      /* Calculate the centre of mass of these atoms */
      centre_of_mass(coord_store,nstored,&mass_x,&mass_y,&mass_z);

      /* Place the centre of mass at the origin */
      adjust_to_origin(mass_x,mass_y,mass_z);

      /* Apply the same translation to the coordinates stored in
         coord_store */
      adjust_stored_coords(nstored,coord_store,mass_x,mass_y,mass_z);

      /* Calculate transformation to orient the principal axis along
         the x-direction */
      principal_components(nstored,transformation,coord_store);

      /* Apply this transformation to the entire structure */
      rotate_whole_structure(transformation);

      /* If the plot is to be in Portrait, rather than Landscape,
         orientation swap the x coords with the y coords */
      if (Picture->Portrait == TRUE)
          swap_x_and_y();

      /* Test if plot needs to be flipped to maintain residue order
         down the page, or left-to-right according to the orientation */
      check_ligand_direction();

      /* Update all the residue- and object-boundaries */
      update_all_boundaries();
    }
}
/***********************************************************************

write_bonds  -  Write out all the bonds in the structure

***********************************************************************/

void write_bonds(void)
{
  int ibond;
  char bond_type_char, ebond_char, rbond_char, single_double, source[7];
  struct bond *bond_ptr;
  struct coordinate *atom_ptr1, *atom_ptr2;

  printf("\nWriting out the bonds file...\n");

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;
  ibond = 0;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Write out provided bond hasn't been marked for deletion */
      if (bond_ptr->bond_type != DELETED)
        {
          bond_type_char = '.';
          if (bond_ptr->bond_type == HBOND)
            bond_type_char = 'H';
          else if (bond_ptr->bond_type == COVALENT)
            bond_type_char = 'c';
          else if (bond_ptr->bond_type == CONTACT)
            bond_type_char = 'n';
          else if (bond_ptr->bond_type == INTERNAL)
            bond_type_char = 'i';
/* v.5.4--> */
          else if (bond_ptr->bond_type == SALT_BRIDGE)
            bond_type_char = 'S';
          else if (bond_ptr->bond_type == DISULPHIDE)
            bond_type_char = 'D';
/* <--v.5.4 */

          /* Determine whether a rotatable bond */
          rbond_char = ' ';
          if (bond_ptr->rotatable_bond == TRUE)
            rbond_char = 'r';
          else if (bond_ptr->end_bond == TRUE)
            rbond_char = 'e';

          /* Determine whether an elastic bond */
          ebond_char = ' ';
          if (bond_ptr->elastic == TRUE)
            ebond_char = 'E';

          /* Get the bond's two atoms */
          atom_ptr1 = bond_ptr->first_atom_ptr;
          atom_ptr2 = bond_ptr->second_atom_ptr;

          /* Determine bond-order, single/double/treble */
          single_double = '-';
          if (bond_ptr->bond_order == DOUBLE)
            single_double = 'd';
          else if (bond_ptr->bond_order == TRIPLE)
            single_double = 't';

          /* Determine where bond info came from */
          if (bond_ptr->bond_source == CONECT)
            strcpy(source,"CONECT");
          else if (bond_ptr->bond_source == CALCULATED)
            strcpy(source,"Calc  ");
          else if (bond_ptr->bond_source == HBPLUS)
            strcpy(source,"HBPLUS");
          else
            strcpy(source,"      ");
          
          /* Write the bond information out to the ligplot.bonds file */
          fprintf(ligplot_bonds_out,"Bond %3d. Atoms%s [%s %s %s %c]",ibond,
                  atom_ptr1->atom_number,atom_ptr1->atom_name,
                  atom_ptr1->residue_ptr->res_name,
                  atom_ptr1->residue_ptr->res_num,
                  atom_ptr1->residue_ptr->chain);
          fprintf(ligplot_bonds_out," ->%s [%s %s %s %c]",
                  atom_ptr2->atom_number,atom_ptr2->atom_name,
                  atom_ptr2->residue_ptr->res_name,
                  atom_ptr2->residue_ptr->res_num,
                  atom_ptr2->residue_ptr->chain);
          fprintf(ligplot_bonds_out,"%c%c%c%7.4f %c %s\n",bond_type_char,
                  rbond_char,ebond_char,bond_ptr->bond_length,
                  single_double,source);

          /* If this is a hydrogen bond, write out to the ligplot.hhb
             file */
/* v.5.4--> */
//          if (bond_ptr->bond_type == HBOND ||
          if (bond_ptr->bond_type == HBOND ||
              bond_ptr->bond_type == SALT_BRIDGE ||
              bond_ptr->bond_type == DISULPHIDE ||
/* <--v.5.4 */
              bond_ptr->bond_type == INTERNAL)
            {
              fprintf(ligplot_hhb_out,"%s %c %s %s     ",
                      atom_ptr1->residue_ptr->res_name,
                      atom_ptr1->residue_ptr->chain,
                      atom_ptr1->residue_ptr->res_num,
                      atom_ptr1->atom_name);
/* v.5.4--> */
//              fprintf(ligplot_hhb_out,"%s %c %s %s    %4.2f\n",
              fprintf(ligplot_hhb_out,"%s %c %s %s    %4.2f",
/* <--v.5.4 */
                      atom_ptr2->residue_ptr->res_name,
                      atom_ptr2->residue_ptr->chain,
                      atom_ptr2->residue_ptr->res_num,
                      atom_ptr2->atom_name,bond_ptr->bond_length);
/* v.5.4--> */
              if (bond_ptr->bond_type == SALT_BRIDGE)
                fprintf(ligplot_hhb_out,"  +/-");
              if (bond_ptr->bond_type == DISULPHIDE)
                fprintf(ligplot_hhb_out,"  S-S");
              fprintf(ligplot_hhb_out,"\n");
/* <--v.5.4 */
            }

          /* If this is a non-bonded contact, write out to the ligplot.nnb
             file */
          else if (bond_ptr->bond_type == CONTACT)
            {
              fprintf(ligplot_nnb_out,"%s %c %s %s     ",
                      atom_ptr1->residue_ptr->res_name,
                      atom_ptr1->residue_ptr->chain,
                      atom_ptr1->residue_ptr->res_num,
                      atom_ptr1->atom_name);
              fprintf(ligplot_nnb_out,"%s %c %s %s    %4.2f\n",
                      atom_ptr2->residue_ptr->res_name,
                      atom_ptr2->residue_ptr->chain,
                      atom_ptr2->residue_ptr->res_num,
                      atom_ptr2->atom_name,bond_ptr->bond_length);
            }

          /* Increment bond counter */
          ibond++;
        }

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }
}
/* v.4.3--> */
/***********************************************************************

write_ligand_fit_coords  -  Write out the ligand coordinates

***********************************************************************/

void write_ligand_fit_coords(void)
{
/* v.5.0--> */
  char filename[FILENAME_LEN];
/* <--v.5.0 */

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;

  FILE *file_ptr;

  /* Set flag that have position anchors */
  Have_Anchors = TRUE;

  /* Open ligplot.lig file */
/* v.5.0--> */
  /* if ((file_ptr = fopen("ligplot.lig","w")) ==  NULL) */
  filename[0] = '\0';
  if (wkdir[0] != '\0')
    strcpy(filename,wkdir);
  strcat(filename,"ligplot.lig");
  if ((file_ptr = fopen(filename,"w")) == NULL)
/* <--v.5.0 */
    {
/* v.5.0--> */
      /* printf("\n*** Unable to open ligplot.lig file\n"); */
      printf("\n*** Unable to open %s file\n",filename);
/* <--v.5.0 */
      exit(1);
    }

  /* Initialise pointer to start of linked list of atom coords */
  atom_ptr = first_atom_ptr;

  /* Loop through all the atoms in the linked list */
  while (atom_ptr != NULL)
    {
      /* Check whether atom belongs to the ligand */
      if (atom_ptr->residue_ptr->inligand == TRUE)
        {
          /* Write out the atom record */
          fprintf(file_ptr,"ATOM  %5s ",atom_ptr->atom_number);
          fprintf(file_ptr,"%4s ",atom_ptr->atom_name);
          fprintf(file_ptr,"%3s %c%5s   %8.3f%8.3f%8.3f%6s%6s",
                  atom_ptr->residue_ptr->res_name,
                  atom_ptr->residue_ptr->chain,
                  atom_ptr->residue_ptr->res_num,
                  atom_ptr->fit_x,
                  atom_ptr->fit_y,
                  atom_ptr->fit_z,
                  atom_ptr->occupancy,
                  atom_ptr->bvalue);
          fprintf(file_ptr,"          %2s\n",atom_ptr->element);

          /* Set atom's coordinates as anchors */
          atom_ptr->have_anchor = TRUE;
          atom_ptr->anchor_pstn_x = atom_ptr->fit_x;
          atom_ptr->anchor_pstn_y = atom_ptr->fit_y;

          /* Mark residue as having anchored atoms */
          atom_ptr->residue_ptr->nanchored_atoms
            = atom_ptr->residue_ptr->natoms;

          /* Mark object as having anchored atoms */
          atom_ptr->residue_ptr->object_ptr->have_anchors = TRUE;
        }

      /* Get pointer to next atom in linked-list */
      atom_ptr = atom_ptr->next;
    }

  /* Write out a TER record */
  fprintf(file_ptr,"TER\n");

  /* Get pointer to the first bond */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Write out covalent bonds as CONECT records */
      if (bond_ptr->bond_type == COVALENT &&
          bond_ptr->bond_type != DELETED)
        {
          /* Get the bond's two atoms */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Only interested if both atoms belong to the ligand */
          if (atom1_ptr->residue_ptr->inligand == TRUE &&
              atom2_ptr->residue_ptr->inligand == TRUE)
            {
              /* Write out as a CONECT record */
              fprintf(file_ptr,"CONECT%5s%5s\n",
                      atom1_ptr->atom_number,atom2_ptr->atom_number);
            }
        }

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Close the file */
  fclose(file_ptr);
}
/* <--v.4.3 */
/* v.4.0--> */
/***********************************************************************

write_rcm_file  -  Write out the co-ordinates of the residue centres
                   of mass to the ligplot.rcm file

***********************************************************************/

void write_rcm_file(void)
{
  int iresid, nresid;

  char chain;

  struct object *object_ptr;
  struct residue *residue_ptr;

  /* Write out the header records */
  fprintf(ligplot_rcm,"                 Res.  Res       Flattened    \n");
  fprintf(ligplot_rcm,"            Atom Name  Num      coordinates   \n");
  fprintf(ligplot_rcm,"            ---- --- -----    ----------------\n");

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
      /* If chain is blank, write out as a hyphen */
      chain = object_ptr->first_residue_ptr->chain;
      if (chain == ' ')
        chain = '-';

      /* Set the pointer to the first of this object's residues */
      residue_ptr = object_ptr->first_residue_ptr;
      nresid = object_ptr->nresidues;
      iresid = 0;

/* v.5.0--> */
      /* Write out if not from a dummy object */
      if (residue_ptr != NULL &&
          residue_ptr->object_ptr->object_type == DUMMY)
        residue_ptr = NULL;
/* <--v.5.0 */

      /* Loop over all this object's residues */
      while (iresid < nresid && residue_ptr != NULL)
        {
          /* Write out the residue's centre-of-mass coords */
          fprintf(ligplot_rcm,"RESDUE      CofM %s %c%s   %8.3f%8.3f\n",
                  residue_ptr->res_name,
                  chain,
                  residue_ptr->res_num,
                  residue_ptr->flattened_mean_x,
                  residue_ptr->flattened_mean_y);

          /* Get pointer to the next residue */
          residue_ptr = residue_ptr->next_residue_ptr;
          iresid++;
        }

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Close the output file */
  fclose(ligplot_rcm);
}
/* <--v.4.0 */
/* v.4.1.1--> */
/***********************************************************************

write_sum_file  -  Calculate summaries of residue-pair interaction data
                   and write to the ligplot.sum file

***********************************************************************/

void write_sum_file(void)
{
  char chain1, chain2;
  char res_name1[4], res_name2[4];
  char res_num1[6], res_num2[6];
  int itype;
  int duplicate;
  struct bond *bond_ptr;
  struct coordinate *atom1_ptr, *atom2_ptr;
  struct residue *residue1_ptr, *residue2_ptr;
  struct pair *first_pair_ptr, *last_pair_ptr, *match_pair_ptr, *pair_ptr;

  printf("\nWriting out the summary interactions file...\n");

  /* Initialise variables */
  first_pair_ptr = NULL;
  last_pair_ptr = NULL;
  match_pair_ptr = NULL;

  /* Initialise pointer to start of linked list of bonds */
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Check for H-bonds and non-bonded contacts */
/* v.5.4--> */
//      if (bond_ptr->bond_type == HBOND ||
      if (bond_ptr->bond_type == HBOND ||
          bond_ptr->bond_type == SALT_BRIDGE ||
          bond_ptr->bond_type == DISULPHIDE ||
/* <--v.5.4 */
          bond_ptr->bond_type == CONTACT)
        {
          /* Get the bond's two atoms */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

          /* Get interaction type: mainchain/mainchain, sidechain/sidechain
             or mainchain/sidechain */
          if (atom1_ptr->side_chain == FALSE &&
              atom2_ptr->side_chain == FALSE)
            itype = MAIN_MAIN;
          else if (atom1_ptr->side_chain == TRUE &&
              atom2_ptr->side_chain == TRUE)
            itype = SIDE_SIDE;
          else
            itype = MAIN_SIDE;

          /* Get the corresponding residues */
          residue1_ptr = atom1_ptr->residue_ptr;
          residue2_ptr = atom2_ptr->residue_ptr;

          /* Set pointer to start of pairs-list */
          duplicate = FALSE;
          match_pair_ptr = NULL;
          pair_ptr = first_pair_ptr;

          /* Loop through all the stored pairs to check whether we already
             have this residue pair */
          while (pair_ptr != NULL)
            {
              /* Check whether both residues match */
              if ((residue1_ptr == pair_ptr->residue1_ptr &&
                   residue2_ptr == pair_ptr->residue2_ptr) ||
                  (residue2_ptr == pair_ptr->residue1_ptr &&
                   residue1_ptr == pair_ptr->residue2_ptr))

                /* Have match */
                match_pair_ptr = pair_ptr;

              /* If residues don't match exactly, then check whether
                 this is a mirror image of an interaction we already have */
              else
                {
                  /* Get the residue names and numbers */
                  chain1 = pair_ptr->residue1_ptr->chain;
                  strcpy(res_name1,pair_ptr->residue1_ptr->res_name);
                  strcpy(res_name2,pair_ptr->residue2_ptr->res_name);
                  chain2 = pair_ptr->residue2_ptr->chain;
                  strcpy(res_num1,pair_ptr->residue1_ptr->res_num);
                  strcpy(res_num2,pair_ptr->residue2_ptr->res_num);

                  /* Check if have a match */
                  if (!strcmp(residue1_ptr->res_name,res_name1) &&
                      !strcmp(residue1_ptr->res_num,res_num1) &&
                      residue1_ptr->chain == chain2 &&
                      !strcmp(residue2_ptr->res_name,res_name2) &&
                      !strcmp(residue2_ptr->res_num,res_num2) &&
                      residue2_ptr->chain == chain1)
                    duplicate = TRUE;
                  if (!strcmp(residue2_ptr->res_name,res_name1) &&
                      !strcmp(residue2_ptr->res_num,res_num1) &&
                      residue2_ptr->chain == chain2 &&
                      !strcmp(residue1_ptr->res_name,res_name2) &&
                      !strcmp(residue1_ptr->res_num,res_num2) &&
                      residue1_ptr->chain == chain1)
                    duplicate = TRUE;
                }

              /* Go to next pair in the list */
              pair_ptr = pair_ptr->next_pair_ptr;
            }

          /* If no match found, then this is a new pair, so need to create
             a new pair record for it */
          if (match_pair_ptr == NULL)
            {
              /* Allocate memory for structure to hold this pair */
              pair_ptr = (struct pair*)malloc(sizeof(struct pair));
              if (pair_ptr == NULL)
                {
                  printf("*** ERROR. Unable to allocate memory for");
                  printf(" struct pair\n");
                  printf("***        Program ligplot terminated with");
                  printf(" error.\n");
                  exit (1);
                }

              /* Initialise this record's fields */
              pair_ptr->residue1_ptr = residue1_ptr;
              pair_ptr->residue2_ptr = residue2_ptr;
              pair_ptr->nhbonds[0] = 0;
              pair_ptr->nhbonds[1] = 0;
              pair_ptr->nhbonds[2] = 0;
              pair_ptr->ncontacts = 0;
              pair_ptr->duplicate = duplicate;
              pair_ptr->next_pair_ptr = NULL;

              /* If this is not the first pair record, make previous pair
                 point here */
              if (last_pair_ptr != NULL)
                last_pair_ptr->next_pair_ptr = pair_ptr;

              /* Otherwise, set it as the first pair */
              else
                first_pair_ptr = pair_ptr;

              /* Save this pointer as the last pair created */
              last_pair_ptr = pair_ptr;
            }

          /* Otherwise, save the matched pair pointer */
          else
            pair_ptr = match_pair_ptr;

          /* Update the statistics on this pair */
/* v.5.4--> */
//          if (bond_ptr->bond_type == HBOND)
          if (bond_ptr->bond_type == HBOND ||
              bond_ptr->bond_type == SALT_BRIDGE ||
              bond_ptr->bond_type == DISULPHIDE)
/* <--v.5.4 */
            pair_ptr->nhbonds[itype]++;
          else
            pair_ptr->ncontacts++;
        }

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

  /* Set pointer to start of pair interactions list */
  pair_ptr = first_pair_ptr;

  /* Loop through all the stored summary pair interaction data and write
     out */
  while (pair_ptr != NULL)
    {
      /* Get the two interacting residues */
      residue1_ptr = pair_ptr->residue1_ptr;
      residue2_ptr = pair_ptr->residue2_ptr;

      /* Write out the residue pair */
      fprintf(ligplot_sum,"Residues: [%s %s %c] -> [%s %s %c]",
              residue1_ptr->res_name,
              residue1_ptr->res_num,
              residue1_ptr->chain,
              residue2_ptr->res_name,
              residue2_ptr->res_num,
              residue2_ptr->chain);
      fprintf(ligplot_sum," %6d%6d%6d   %6d",
              pair_ptr->nhbonds[0],
              pair_ptr->nhbonds[1],
              pair_ptr->nhbonds[2],
              pair_ptr->ncontacts);

      /* If a duplicate, then mark it as such */
      if (pair_ptr->duplicate == TRUE)
        fprintf(ligplot_sum,"   MIRROR");
      fprintf(ligplot_sum,"\n");

      /* Get pointer to next pair */
      pair_ptr = pair_ptr->next_pair_ptr;
    }

  /* Close the output file */
/* v.4.4.3--> */
/*  fclose(ligplot_sum); */
/* <--v.4.4.3 */
}
/* <--v.4.1.1 */
/* v.4.1--> */
/***********************************************************************

write_drw_file  -  Write out all the objects on the plot to the
                   ligplot.drw file

***********************************************************************/

void write_drw_file(void)
{
  int iatom, iobject, iresid;
  int atom_count, natoms, nresid, res_count;
  int npoints;
/* v.5.0.4--> */
  int wanted;
/* <--v.5.0.4 */

  char chain;
  char single_double;

  float max_x, max_y, min_x, min_y, title_x, title_y;

  struct bond *bond_ptr;
  struct coordinate *atom_ptr, *atom1_ptr, *atom2_ptr;
/* v.5.0.4--> */
/*  struct object *object_ptr; */
  struct object *object_ptr, *object1_ptr, *object2_ptr;
/* <--v.5.0.4 */
  struct residue *residue_ptr;

  printf("\nWriting out the drw file...\n");

  /* Initialise variables */
  atom_count = 0;
  iobject = 0;
  res_count = 0;
  title_x = 0.0;
  title_y = 0.0;
  npoints = 0;

  /* Initialise pointer to the first stored object */
  object_ptr = first_object_ptr;

  /* Loop through all objects */
  while (object_ptr != NULL)
    {
/* v.5.0--> */
      /* Write out if not the dummy object */
      if (object_ptr->object_type != DUMMY)
        {
/* <--v.5.0 */
          /* Write out #O record as start of object */
          fprintf(ligplot_drw,"#O%d\n",object_ptr->interface);

          /* If chain is blank, write out as a hyphen */
          chain = object_ptr->first_residue_ptr->chain;
          if (chain == ' ')
            chain = '-';

          /* Set the pointer to the first of this object's residues */
          residue_ptr = object_ptr->first_residue_ptr;
          nresid = object_ptr->nresidues;
          iresid = 0;

          /* Loop over all this object's residues */
          while (iresid < nresid && residue_ptr != NULL)
            {
              /* Write out #R record as start of residue */
              fprintf(ligplot_drw,"#R\n");

              /* Write out the residue's type, name, number and chain-id */
              fprintf(ligplot_drw,"%d%s%s%c",
                      residue_ptr->object_ptr->object_type,
                      residue_ptr->res_name,
                      residue_ptr->res_num,
                      chain);

              /* If residue has a name label, then add its location to line */
              if (residue_ptr->label_x != 0.0 && residue_ptr->label_y != 0.0)
                fprintf(ligplot_drw," %8.3f %8.3f",
                        residue_ptr->label_x,
                        residue_ptr->label_y);

              /* Line throw */
              fprintf(ligplot_drw,"\n");

              /* Write out #A record as start of atom-list for this residue */
              fprintf(ligplot_drw,"#A\n");

              /* Get pointer to this residue's first atom */
              atom_ptr = residue_ptr->first_atom_ptr;
              
              /* Initialise atom count */
              iatom = 0;
              natoms = residue_ptr->natoms;

              /* Loop over all the residue's atoms, checking whether
                 they are involved in interactions */
              while (iatom < natoms && atom_ptr != NULL)
                {
                  /* Write out this atom's name and coords */
                  fprintf(ligplot_drw,"%4s ",atom_ptr->atom_name);
                  fprintf(ligplot_drw,"%8.3f %8.3f",
                          atom_ptr->x,atom_ptr->y);

                  /* If atom has a name label, then add to line */
/* v.5.0-->
                     if (atom_ptr->print_name[0] != '\0' &&
                     object_ptr->object_type != WATER)
   <--v.5.0 */
                  fprintf(ligplot_drw," %8.3f %8.3f",
                          atom_ptr->label_x,atom_ptr->label_y);
/* v.5.0--> */
                  /* Add the atom's details from the PDB file */
                  fprintf(ligplot_drw," # %8.3f %8.3f %8.3f",
                          atom_ptr->original_x,atom_ptr->original_y,
                          atom_ptr->original_z);
/* <--v.5.0 */

                  /* Line throw */
                  fprintf(ligplot_drw,"\n");

                  /* Store the atom's reference number for use when writing 
                     out the bonds */
                  atom_ptr->atom_count = atom_count;

                  /* If this is the first atom, store its coords as the max
                     and min values */
                  if (npoints == 0)
                    {
                      /* Store atom's coords */
                      min_x = atom_ptr->x;
                      max_x = atom_ptr->x;
                      min_y = atom_ptr->y;
                      max_y = atom_ptr->y;
                    }

                  /* Otherwise, update the maximum and minimum values */
                  else
                    {
                      /* Update if outside existing values */
                      if (atom_ptr->x < min_x)
                        min_x = atom_ptr->x;
                      if (atom_ptr->x > max_x)
                        max_x = atom_ptr->x;
                      if (atom_ptr->y < min_y)
                        min_y = atom_ptr->y;
                      if (atom_ptr->y > max_y)
                        max_y = atom_ptr->y;
                    }

                  /* Increment count of points processed */
                  npoints++;
                  /* Get pointer to the next atom */
                  atom_ptr = atom_ptr->next;
                  iatom++;
                  atom_count++;
                }

              /* Get pointer to the next residue */
              residue_ptr = residue_ptr->next_residue_ptr;
              iresid++;
              res_count++;
            }

          /* Increment count of objects processed */
          iobject++;
/* v.5.0--> */
        }
/* <--v.5.0 */

      /* Get pointer for next object */
      object_ptr = object_ptr->next_object_ptr;
    }

  /* Initialise pointer to start of linked list of bonds */
  fprintf(ligplot_drw,"#B\n");
  bond_ptr = first_bond_ptr;

  /* Loop through all the bonds in the linked list */
  while (bond_ptr != NULL)
    {
      /* Write out provided bond hasn't been marked for deletion */
      if (bond_ptr->bond_type != DELETED)
        {
          /* If this is not a hydrogen bond, the store bond-order in
             the label */
          if (bond_ptr->bond_type != HBOND &&
/* v.5.4--> */
//              bond_ptr->bond_type != INTERNAL)
              bond_ptr->bond_type != INTERNAL &&
              bond_ptr->bond_type != SALT_BRIDGE &&
              bond_ptr->bond_type != DISULPHIDE)
/* <--v.5.4 */
            {
              /* Determine bond-order, single/double/treble */
              single_double = 's';
              if (bond_ptr->bond_order == DOUBLE)
                single_double = 'd';
              else if (bond_ptr->bond_order == TRIPLE)
                single_double = 't';
/* v.4.5--> */
              else if (bond_ptr->bond_order == AROMATIC)
                single_double = 'a';
/* <--v.4.5 */
/* v.4.4--> */
/* v.4.5.4--> */
              // if (Include->Double_Bonds == FALSE)
              if (Include->Double_Bonds == FALSE &&
                  Write_Bond_Types == FALSE)
/* <--v.4.5.4 */
                single_double = 's';
/* <--v.4.4 */

              /* Store in label, if not default (s) */
              if (single_double != 's')
                {
                  bond_ptr->label[0] = single_double;
                  bond_ptr->label[1] = '\0';
                }
            }

          /* Get the bond's two atoms */
          atom1_ptr = bond_ptr->first_atom_ptr;
          atom2_ptr = bond_ptr->second_atom_ptr;

/* v.5.0.4--> */
          /* If both atoms belong to the same hydrophobic group, then
             we don't want this bond */
          wanted = TRUE;
          object1_ptr = atom1_ptr->residue_ptr->object_ptr;
          object2_ptr = atom2_ptr->residue_ptr->object_ptr;
          if (object1_ptr == object2_ptr &&
              object1_ptr->object_type == HYDROPHOBIC)
            wanted = FALSE;

          /* If bond is wanted, then write out */
          if (wanted == TRUE)
            {
/* <--v.5.0.4 */
              /* Write the bond information out to the ligplot.drw file */
              fprintf(ligplot_drw,"%d %d %d",bond_ptr->bond_type,
                      atom1_ptr->atom_count,atom2_ptr->atom_count);

              /* Add the bond label, if there is one */
              if (bond_ptr->label[0] != '\0')
                fprintf(ligplot_drw," %s",bond_ptr->label);

              /* Line-throw */
              fprintf(ligplot_drw,"\n");
/* v.5.0.4--> */
            }
/* <--v.5.0.4 */
        }

      /* Get pointer to next atom in linked-list */
      bond_ptr = bond_ptr->next_bond_ptr;
    }

/* v.5.4.5 --> */
  /* Get the .drw title - if there is one */
  if (Drw_Title[0] != '\0')
    strcpy(Print_Title,Drw_Title);
/* <-- v.5.4.5 */

  /* Write out the plot title */
  if (Print_Title[0] != '\0')
    {
      /* Place the title such that it is centered left and right and
         below the lowest atom on the plot */
      if (npoints > 0)
        {
          title_x = (float) ((min_x + max_x) / 2.0);
          title_y = (float) (min_y - 3.0);
        }

      fprintf(ligplot_drw,"#T\n");
      fprintf(ligplot_drw,"%s %8.3f %8.3f\n",Print_Title,title_x,
              title_y);
    }

  /* Close the output file */
  fclose(ligplot_drw);
}
/* <--v.4.1 */
/* v.4.2--> */
/***********************************************************************

close_output_files  -  Close any output files that are open

***********************************************************************/

void close_output_files(void)
{
  /* Close all the relevant files */
  fclose(ligplot_pdb);
/* v.4.4.1--> */
/*  fclose(ligplot_frm); */
/*  fclose(ligplot_rcm); */
/* <--v.4.4.1 */
  if (Write_Res_File == TRUE)
    fclose(ligplot_res);
  fclose(ligplot_bonds_out);
  fclose(ligplot_hhb_out);
  fclose(ligplot_nnb_out);
  fclose(ligplot_sum);
}
/* <--v.4.2 */




/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
/* v.5.4--> */
//  int atom_counter, nhbonds, ncontacts, nlinks;
  int atom_counter, nlinks;
  int i, nbonds[S_NTYPES];
/* <--v.5.4 */
  int ndeleted, nligands, nobjects;
  int iseed;
  int nligand_objects;
  int show_deletions;
/* v.4.1.2--> */
  int debug;
/* <--v.4.1.2 */
/* v.4.2--> */
  int have_range;
  int iligand, nlig_molecs;
/* <--v.4.2 */
/* v.4.4.4--> */
  int bounds[4];
/* <--v.4.4.4 */
/* v.5.3 --> */
  int have_rcm;
/* <--v.5.3 */
/* v.5.4.2--> */
  int no_abort;
/* <--v.5.4.2 */

  float scale_factor;
  float centre_x, centre_y, min_x, max_x, min_y, max_y;
  float matrix[3][3];
/* v.4.3--> */
  float energy;
/* <--v.4.3 */
/* v.4.4.4--> */
  float adjust, height, width;
/* <--v.4.4.4 */

  char bonds_name[FILENAME_LEN], nnb_name[FILENAME_LEN],
    hhb_name[FILENAME_LEN], pdb_name[FILENAME_LEN];
  char chain_identi, res_num1[6], res_num2[6];
  char special_res[MAXSPECIAL_RES][15];
/* v.4.0--> */
/* v.4.2--> */
  /*  char res_name1[4], res_name2[4]; */
/* <--v.4.2 */
  char pdb_code[5], res_name1[4], res_name2[4];
/* <--v.4.0 */
/* debug 
  char chain, res_name[4], res_num[6];
 debug */

/* v.4.2--> */
  XMAS *xmas;
/* <--v.4.2 */

  /* Initialise global variables */
  scale_factor = (float) (1.0);
/* v.4.0--> */
  Interface_Plot = FALSE;
/* <--v.4.0 */
/* v.5.3.4--> */
  From_Dimer = FALSE;
/* <--v.5.3.4 */
/* v.4.3--> */
  Plot_from_drw = FALSE;
/* <--v.4.3 */
  Nwarnings = 0;
  Water_as_Ligand = FALSE;
  Metal_as_Ligand = FALSE;
  show_deletions = FALSE;
  Write_Res_File = FALSE;
/* v.4.3--> */
  Write_Lig_File = FALSE;
/* <--v.4.3 */
/* v.4.4.4--> */
  Scale_Plot = FALSE;
/* <--v.4.4.4 */
/* v.4.5.4--> */
  Write_Bond_Types = FALSE;
/* <--v.4.5.4 */
/* v.4.1.2--> */
  debug = FALSE;

/* <--v.4.1.2 */
/* v.4.4--> */
  ndeleted = 0;
/* <--v.4.4 */

/* v.5.4--> */
  /* Initialise counts */
  for (i = 0; i < S_NTYPES; i++)
    nbonds[i] = 0;
/* <--v.5.4 */

  /* Initialise the plot parameters and global variables */
  energy = 0.0;
  initialise_parameters();
/* v.5.3 --> */
  have_rcm = FALSE;
/* <--v.5.3 */

  /* Sort out the command line arguments*/
/* v.4.0--> */
/*  get_command_arguments(argv,argc - 1,pdb_name,res_num1,res_num2,
                        &chain_identi,hhb_name,nnb_name,bonds_name,
                        Print_Title); */
/* v.4.1.2--> */
/*  get_command_arguments(argv,argc - 1,pdb_name,res_name1,res_num1,
                        res_name2,res_num2,&chain_identi,hhb_name,
                        nnb_name,bonds_name,Print_Title); */
  get_command_arguments(argv,argc - 1,pdb_name,res_name1,res_num1,
                        res_name2,res_num2,&chain_identi,hhb_name,
/* v.5.4.2--> */
//                        nnb_name,bonds_name,Print_Title,&debug);
                        nnb_name,bonds_name,Print_Title,&debug,
                        &no_abort);
/* <--v.5.4.2 */
/* <--v.4.1.2 */
/* <--v.4.0 */

  /* Open the parameter file and read user-definable information in */
/* v.5.0--> */
//  if (run_type == LIGPLUS_PLOT)
//    read_ligplus_params();
//  else
/* <--v.5.0 */
    read_in_parameters(special_res);

/* v.5.0.3--> */
    /* Override certain of the parameter just read in, based on command-
       line parameters entered by user */

    /* If no title entered, and option to use filename as the default
       title, then store the filename as the title */
    if (Include->Filename_for_Title == TRUE && Print_Title[0] == '\0' &&
        Print_as_is == FALSE && Interface_Plot == FALSE)
      strcpy(Print_Title,root_name);

    /* If the ligand is a water molecule, then make sure that water
       molecules are allowed */
    if (Water_as_Ligand == TRUE)
      Include->Waters = TRUE;
/* <--v.5.0.3 */

  /* Check and adjust colour names in colour definitions */
/* v.4.2--> */
  match_colour_names();
/* <--v.4.2 */

  /* Initialise random number generator */
  iseed = -1;
  get_random_number(&iseed,Random_Start);

  /* Determine the maximum object and text sizes in the plot */
/* v.4.2--> */
/*  define_object_sizes(scale_factor,FALSE);
  define_text_sizes(scale_factor,FALSE); */
/* <--v.4.2 */

/* v.4.1--> */
  /* Write the plot parameters to the ligplot.drw file */
/* v.4.2--> */
/*  write_drw_headers(); */
/* <--v.4.2 */
/* <--v.4.1 */

/* v.4.2--> */
  /* Initialise number of ligand molecules to be plotted */
  nlig_molecs = 1;

  /* If reading from a XMAS file, open the file and cache all the data */
  if (Read_Xmas == TRUE)
    {
      /* Cache all the data */
      xmas = read_in_xmas_file(pdb_name);

      /* If looping through all the ligands, then get how many there are
         and which is the first */
      if (Plot_All_Xmas_Ligands == TRUE)
        nlig_molecs = get_nlig_xmas_molecs(xmas,&Ligand_Molecule_id);
    }

  /* Loop over all the ligand molecules to be processed */
  for (iligand = 0; iligand < nlig_molecs; iligand++)
    {
      /* Re-initialise the variables for the current plot */
      initialise_variables();
      scale_factor = (float) (1.0);

      /* Determine the maximum object and text sizes in the plot */
      define_object_sizes(scale_factor,FALSE);
      define_text_sizes(scale_factor,FALSE);

      /* If plotting all ligands, show which one we're on now */
      if (nlig_molecs > 1)
        {
          printf("\n");
          printf("Plotting ligand %3d of %3d (molecule-id = %3d)\n",
/* v.4.2.2--> */
                 /* iligand + 1,Ligand_Molecule_id,nlig_molecs); */
                 iligand + 1,nlig_molecs,Ligand_Molecule_id);
/* <--v.4.2.2 */
          printf("----------------------------------------------\n");
          printf("\n");
        }

      /* If have a molecule identifier for the ligand, convert into
         a residue range */
      if (Ligand_Molecule_id > -1)
        {
          have_range
            = get_xmaslig_start_end(xmas,res_name1,res_num1,
                                    res_name2,res_num2,&chain_identi);
          if (have_range == FALSE)
            {
              printf("*** ERROR. Unable to find ligand (molecule no. ");
              printf("%d) in XMAS file\n",Ligand_Molecule_id);
              printf("***        Program ligplot terminated with");
              printf(" error.\n");
              exit (1);
            }

          /* If plotting all ligands, form the plot title for this one */
          if (Plot_All_Xmas_Ligands == TRUE)
            form_plot_title(xmas,pdb_code,res_name1,res_num1,res_name2,
                            res_num2,chain_identi);
        }
/* <--v.4.2 */

/* v.4.2--> */
      /* Write the plot parameters to the ligplot.drw file */
/* v.4.3--> */
      if (Plot_from_drw == FALSE)
/* <--v.4.3 */
        write_drw_headers(pdb_code);
/* <--v.4.2 */

      /* Locate the start- and end-atoms and residues of the ligand
         in the PDB file */
/* v.4.0--> */
/* v.4.3--> */
/*      if (Interface_Plot == FALSE) */
      if (Interface_Plot == FALSE && Plot_from_drw == FALSE)
/* <--v.4.3 */
/* <--v.4.0 */
/* v.4.0--> */
        /*    get_ligand_start_end(pdb_name,res_num1,res_num2,chain_identi); */
/* v.4.2--> */
        /*    get_ligand_start_end(pdb_name,res_name1,res_num1,res_name2,res_num
2,
              chain_identi); */
        get_ligand_start_end(pdb_name,xmas,res_name1,res_num1,res_name2,
                             res_num2,chain_identi);
/* <--v.4.2 */
/* <--v.4.0 */

      /* If printing molecule as it stand, read in all the bonds and bond-
         lengths from the ligplot.bonds file */
      if (Print_as_is == TRUE)
        {
/* v.4.3--> */
          /* If reading in from the ligplot.drw file, get all the
             objects, residues, atoms and bonds */
          if (Plot_from_drw == TRUE)
            nobjects = read_drw_file(pdb_name);

          /* Otherwise, read in from the appropriate PDB-format file */
          else
            {
/* <--v.4.3 */
              /* Read in all the atom coordinates from the PDB file */
/* v.4.2--> */
              /*      atom_counter = read_pdb_file(pdb_name); */
              atom_counter = read_pdb_file(pdb_name,xmas);
/* <--v.4.2 */

              /* Read in the bonds from ligplot.bonds */
              get_bonds(bonds_name);

              /* Define the objects that will appear on the final 
                 LIGPLOT diagram */
              define_objects(&nobjects,&nligand_objects);
            }

          /* Update all the residue boundaries */
          update_all_boundaries();
        }

      /* Otherwise, define all the bonds using the CONEC records, distances,
         and the HBPLUS output */
      else
        {
          /* Use the CONEC records to find which protein residues the ligand
             is covalently bonded to */
          nlinks = 0;
/* v.4.2--> */
          /*      read_conec_records(pdb_name,&nlinks); */
          read_conec_records(pdb_name,&nlinks,xmas);
/* <--v.4.2 */

          /* Check the stored CONEC records by calculating distances between
             the atoms involved */
          if (nlinks > 0)
/* v.4.2--> */
            /*      verify_conec_records(pdb_name); */
            verify_conec_records(pdb_name,xmas);
/* <--v.4.2 */

          /* Read in the ligand's hydrogen bonds from the .hhb file */
          if (Include->Hbonds == TRUE)
            {
/* v.4.1.2--> */
              /*        read_hbplus_file(hhb_name,HHB_FILE,&nhbonds); */
/* v.4.2--> */
              /*        read_hbplus_file(hhb_name,HHB_FILE,&nhbonds,debug); */
/* v.5.4--> */
//              read_hbplus_file(hhb_name,HHB_FILE,&nhbonds,debug,xmas);
              read_hbplus_file(hhb_name,HHB_FILE,nbonds,debug,xmas);
/* <--v.5.4 */
/* <--v.4.2 */
/* <--v.4.1.2 */

              /* Include special case residues into the H_bond list - ie 
                 residues that are connected to other residues connected to 
                 the ligand */
              if (Include->Linked_Residues == TRUE)
/* v.4.2--> */
                /*          special_resinc(hhb_name,HHB_FILE,special_res); */
                special_resinc(hhb_name,HHB_FILE,special_res,xmas);
/* <--v.4.2 */
            }

          /* Read in the hydrophobic contacts from the .nnb file */
          if (Include->Hydrophobics == TRUE)
/* v.4.1.2--> */
            /*      read_hbplus_file(nnb_name,NNB_FILE,&ncontacts); */
/* v.4.2--> */
            /*      read_hbplus_file(nnb_name,NNB_FILE,&ncontacts,debug); */
/* v.4.2.1--> */
            /*          read_hbplus_file(hhb_name,NNB_FILE,&nhbonds,debug,xmas);
 */
/* v.5.4--> */
//            read_hbplus_file(nnb_name,NNB_FILE,&ncontacts,debug,xmas);
            read_hbplus_file(nnb_name,NNB_FILE,nbonds,debug,xmas);
/* <--v.5.4 */
/* <--v.4.2.1 */
/* <--v.4.2 */
/* <--v.4.1.2 */

          /* Read in all the atom coordinates from the PDB file */
/* v.4.2--> */
          /*      atom_counter = read_pdb_file(pdb_name); */
          atom_counter = read_pdb_file(pdb_name,xmas);
/* <--v.4.2 */
/* v.5.4--> */
//          printf("   Number of hydrogen bonds stored   = %7d\n",nhbonds);
//          printf("   Number of contacts stored         = %7d\n",ncontacts);
          printf("   Number of hydrogen bonds stored   = %7d\n",
                 nbonds[S_HBONDS]);
          printf("   Number of contacts stored         = %7d\n",
                 nbonds[S_CONTACTS]);
          printf("   Number of salt bridges stored     = %7d\n",
                 nbonds[S_SALT_BRIDGES]);
          printf("   Number of disulphides stored      = %7d\n",
                 nbonds[S_DISULPHIDES]);
/* <--v.5.4 */

          /* If there have been no atoms read in, then abort */
          if (atom_counter == 0)
            {
              printf("*** No atoms stored. Program terminated\n");
              return(-1);
            }

          /* Define the objects that will appear on the final
             LIGPLOT diagram */
          define_objects(&nobjects,&nligand_objects);

          /* Calculate covalent connectivity */
          covalent_connectivity();

/* v.4.4--> */
          /* Check whether any residues interpenetrate */
          ndeleted = interpenetration_check(ndeleted);
/* <--v.4.4 */

          /* Combine objects where necessary */
          combine_objects();

          /* Get rid of any objects that claim to be ligands, but aren't
             (eg protein residues with the same residue numbers as the
             ligand residues) */
          if (nligand_objects > 1)
            delete_false_ligands();

          /* Pick up H-bonds from .hhb file HBPLUS output, and add H-bond
             connectivity */
          add_hbplus_bonds();

          /* Open the output PDB file */
/* v.4.0--> */
          /*      open_pdb_output(res_num1,res_num2,chain_identi); */
          open_pdb_output(res_name1,res_num1,res_name2,res_num2,chain_identi);
/* <--v.4.0 */

          /* Open the output bonds files: ligplot.bonds, ligplot.hhb
             and ligplot.nnb */
          open_bonds_output();
        }

/* v.4.3--> */
      /* Read in the hbadd.bonds file, if present, to pick up the bond
         order information written out by hbadd.c from the components.cif
         file */
      get_bond_orders();

      /* Get the bond orders for any standard amino acid residues and DNA
         bases */
      get_standard_bond_orders();
/* <--v.4.3 */

/* v.4.4--> */
      /* Correct for "covalent bonds" to metal ions */
      check_metal_bonds();
/* <--v.4.4 */

      /* For each bond, determine which other bonds sprout off its two
         ends */
      bond_linkage();

      /* Remove any atoms that are not connected to other atoms in the
         same residue */
/* v.4.0--> */
      if (Print_as_is == FALSE)
/* <--v.4.0 */
/* v.4.1.2--> */
        /*    check_connections(); */
        check_connections(debug);
/* <--v.4.1.2 */

      /* Check whether any objects are cyclic */
      check_for_cyclics();

      /* Split objects if there are no covalent bonds spanning adjacent
         residues */
      split_objects(&nobjects,&nligands);

      /* Determine which bonds relate to each object */
      assign_bonds_to_objects(&ndeleted,show_deletions);

      /* If any cyclic peptides identified, need to delete any bond
         connections involving covalent bonds that have now been made
         elastic */
/* v.4.3--> */
      if (Print_as_is == FALSE)
/* <--v.4.3 */
        delete_unwanted_bond_connections();

      /* Mark all atoms that are reachable from the ligand residue(s) via
         one or more bonds */
/* v.4.0--> */
      if (Interface_Plot == FALSE && Print_as_is == FALSE)
        {
/* <--v.4.0 */
/* v.4.1.2--> */
          /* If in debug mode, list all the bonds */
          if (debug == TRUE)
            list_bonds();
/* <--v.4.1.2 */

/* v.5.1--> */
          /* If filtering out unwanted waters, mark them for deletion */
          if (Filter_Waters == TRUE)
            mark_unwanted_waters();
/* <--v.5.1 */

          /* Mark all atoms that are reachable from the ligand */
          mark_reachable_atoms();

          /* Delete any residues that are unreachable */
          delete_unreachable_residues(&ndeleted,show_deletions);
/* v.4.0--> */
        }
/* <--v.4.0 */

      /* If any residues have been deleted, then show explanatory message */
      if (show_deletions == TRUE && ndeleted > 0)
        show_deleted_residues_message(ndeleted);

      /* Determine each atom's size */
      get_atom_sizes();

      /* For each atom, create a linked list of the atoms covalently bonded
         to it */
      atom_linkage();

      /* Determine which bonds in the structure are "rotatable" */
      identify_rotatable_bonds();

      /* If structure not to be plotted just as it stands, then perform
         all the necessary processing */
      if (Print_as_is == FALSE)
        {
/* v.4.0--> */
          /* Write out any covalent bonds between the ligand and the protein
             to the .res file */
          if (Write_Res_File == TRUE)
            write_lig_prot_bonds();
/* <--v.4.0 */

          /* Print out the ligand and H-bonded sidechains into ligplot.frm */
/* v.4.2--> */
          /*      write_frm_file(pdb_name); */
          write_frm_file(pdb_name,xmas);
/* <--v.4.2 */

          /* Mark for deletion unwanted atoms in hydrophobic groups */
/* v.5.0.4--> */
/*          pare_down_hydrophobics(); */
/* <--v.5.0.4 */

          /* Delete any unwanted atoms */
/* v.4.1.2--> */
          /*      delete_atoms_and_bonds(); */
          delete_atoms_and_bonds(debug,"main");
/* <--v.4.1.2 */

          /* Transform the original coordinates, with each object oriented
             along its principal axes, for fitting of objects during flattening
             to be as close to original conformation as possible, and for
             laying the objects out in their approximate relative
             arrangements */
          calc_fit_coordinates(nligands);

/* v.5.0--> */
          /* Check for an input .rcm file containing residue centres of mass
             to which the current plot's residues are to be restrained */
/* v.5.3 --> */
//          read_rcm_file(pdb_name);
          have_rcm = read_rcm_file(pdb_name);
/* <--v.5.3 */
/* <--v.5.0 */

          /* Flatten all the objects, one by one */
/* v.4.1.2--> */
          /*      flatten_all_objects(&iseed); */
          flatten_all_objects(&iseed,debug);
/* <--v.4.1.2 */
          /* Score all the objects and then sort them by their relative
             importance */
          score_objects();

          /* Sort objects by their weight scores */
          sort_objects();

/* v.4.0--> */
          /* Calculate the mean coordinates of each residue, for use in the
             energy minimization process */
          calc_residue_means();

/* v.5.0--> 
          read_rcm_file(pdb_name);
 <--v.5.0 */

          /* For interface plot, work out how to place the residues on either
             side of an interface */
          if (Interface_Plot == TRUE)
/* v.5.3 --> */
            {
              if (have_rcm == FALSE)
/* <--v.5.3 */
                lay_interface_objects_out();
/* v.5.3 --> */
            }
/* <--v.5.3 */
          else
/* <--v.4.0 */

            /* Otherwise, lay objects out on the page around the ligand */
            lay_objects_out();

/* v.4.3--> */
          /* If writing for PDBsum generation, then write out the
             ligand coordinates (as previously transformed s.t. principal
             axes lie along the x-, y- and z-axes */
          if (Write_Lig_File == TRUE)
            write_ligand_fit_coords();
/* <--v.4.3 */

/* v.5.1--> */
          /* If only fitting ligand molecule to restraints written out
             by ligfit, then try all possible bond flips */
/* v.5.3.4--> */
//          if (ligfit == TRUE)
          if (ligfit == TRUE && dimfit == FALSE)
/* <--v.5.3.4 */
             energy = flip_ligand_bonds();
/* <--v.5.1 */

          /* Minimise the energy of the whole plot to reduce the numbers
             of clashes and overlaps to a minimum */
/* v.4.3--> */
          /*        energy_minimize(&iseed,nligands); */
/* v.5.2.1 --> */
          else
/* <-- v.5.2.1 */
/* v.5.4.4 --> */
//            energy = energy_minimize(&iseed,nligands);
            energy = energy_minimize(&iseed,nligands,atom_counter);
/* <-- v.5.4.4 */
/* <--v.4.3 */

/* v.4.4--> */
          /* If picture is a disaster, then abort */
/* v.5.4.2--> */
//          if (energy > ABORT_ENERGY)
          if (energy > ABORT_ENERGY && no_abort == FALSE)
/* <--v.5.4.2 */
            {
              /* Show error message and give up */
              printf("*** ERROR. Can't do this one - too difficult\n");
              printf("***        Minimization failed. Aborting ...\n");

/* v.5.4.4 --> */
              /* Write message out to the ligplot.drw file */
              fprintf(ligplot_drw,
                      "ERROR: Minimization failed. Problem too difficult\n");

              /* Close the .drw file */
              fclose(ligplot_drw);
/* <-- v.5.4.4 */

              exit(-2);
            }
/* <--v.4.4 */

          /* Reorientate ligand so its principal axis lies up the y-axis*/
/* v.4.0--> */
/* v.4.3--> */
          /*        if (Interface_Plot == FALSE && Have_Anchors == FALSE) */
          if (Interface_Plot == FALSE && Have_Anchors == FALSE &&
              Write_Lig_File == FALSE)
/* <--v.4.3 */
/* <--v.4.0 */
            ligand_orientate();
        }

/* v.5.3.4--> */
      /* If have an unconstrained DIMPLOT diagram, shuffle the residues
         up to give a more compact plot */
      if (Interface_Plot == TRUE && From_Dimer == TRUE && ligfit == FALSE)
        { 
         /* Tighten up the plot */
          tighten_dimplot();
        }
/* <--v.5.3.4 */

      /* If a further rotation is required, apply it now */
/* v.5.3.4--> */
//      if (Picture->Rotation_Angle != 0.0)
      if (Picture->Rotation_Angle != 0.0 && Picture->Allow_Rotation == TRUE)
/* <--v.5.3.4 */
        {
          /* Calculate the rotation matrix */
          calc_rotation_matrix(Picture->Rotation_Angle,matrix);
      
          /* Apply the rotation to the whole structure */
          rotate_whole_structure(matrix);

          /* Update all residue boundaries */
          update_all_boundaries();
        }

      if (Print_as_is == FALSE)
        {
          /* Print out the new atom coordinates as a PDB file */
          write_pdb_file();

          /* Write out the bond data to the ligplot.bonds, ligplot.hhb
             and ligplot.nnb files */
          write_bonds();

/* v.4.0--> */
          /* Write the residue centres of mass to file ligplot.rcm */
          write_rcm_file();
/* <--v.4.0 */

/* v.4.1.1--> */
          /* Calculate residue-pair interaction summaries */
          write_sum_file();
/* <--v.4.1.1 */
        }

      /* Plot molecule as a PostScript file */
  
/* v.4.4--> */
      /* If plotting in chemical notation, switch off all the atoms */
      if (Include->Chemical_Notation == TRUE)
        {
          /* Switch off "sphere" representations */
          Include->Ligand_Atoms = FALSE;
          Include->Nonligand_Atoms = FALSE;
          Include->Water_Atoms = FALSE;
        }
/* <--v.4.4 */

      /* Open the PostScript file */
      psopen_();

      /* If picture to be plotted in landscape mode, then rotate page */
      if (Picture->Portrait == FALSE)
        {
          adjust_for_landscape();
          psrot_(BBOXX2 + BBOXX1,0.0);
        }

      /* Print the title, if there is one */
      if (Print_Title[0] != '\0')
        print_title(Print_Title);

      /* Write out the key to symbols used, if required */
      if (Include->Key == TRUE)
        print_symbol_key();

      /* Calculate the minimum and maximum x-y extent of the atoms to
         be plotted */
      extract_minmax(&min_x,&max_x,&min_y,&max_y);

/* v.4.4.4--> */
      /* If plotting chemical LIGPLOTs for the Enzyme Structures Database,
         adjust max and min if have a small molecule */
      if (Scale_Plot == TRUE)
        {
          /* Get height and width of molecule */
          width = max_x - min_x;
          height = max_y - min_y;

          /* Check whether molecule is small compared to the "standard
             size" */
          if (width < CHEM_WIDTH && height < CHEM_HEIGHT)
            {
              /* Calculate adjustment to maximum and minumum parameters */
              adjust = (CHEM_WIDTH - width) / 2.0;
              min_x = min_x - adjust;
              max_x = max_x + adjust;
              adjust = (CHEM_HEIGHT - height) / 2.0;
              min_y = min_y - adjust;
              max_y = max_y + adjust;
            }
        }
/* <--v.4.4.4 */

      /* Calculate the scaling factor to fit the coords on the
         PostScript page */
      scale(&scale_factor,min_x,max_x,min_y,max_y);

      /* Get the equivalent transformation of the centre-point of the
         picture to be plotted, to the centre-point of the PostScript page */
      centre_point(&centre_x,&centre_y,scale_factor,
/* v.4.4.4--> */
/*                   min_x,max_x,min_y,max_y); */
/* v.5.1.4--> */
/*                   min_x,max_x,min_y,max_y,bounds); */
                   min_x,max_x,min_y,max_y,bounds,Clip_Diagram);
/* <--v.5.1.4 */
/* <--v.4.4.4 */

      /* Calculate scaled sizes of the different items in the plot */
      pscomm_(" ");
      pscomm_("MAIN PLOT");
      pscomm_("Scaled object and text sizes");
      define_object_sizes(scale_factor,TRUE);
      define_text_sizes(scale_factor,TRUE);

      /* Plot the PostScript picture */
      plot_postscript_picture(centre_x,centre_y,scale_factor);

      /* Close the PostScript file */
/* v.4.4.4--> */
/*      psclos_(); */
      psclos_(bounds);
/* <--v.4.4.4 */

/* v.4.1--> */
/* v.4.3--> */
      /* Write out .drw file and close all output files */
      if (Plot_from_drw == FALSE)
        {
/* <--v.4.3 */
          /* Write out all the plot objects to the ligplot.drw file */
          write_drw_file();
/* <--v.4.1 */

/* v.4.2--> */
          /* Close any open output files */
/* v.4.4.4--> */
          if (Print_as_is == FALSE)
/* <--v.4.4.4 */
            close_output_files();

          /* Free up memory used by all the structures */
/* v.4.5.3--> */
          // free_up_memory();
/* <--v.4.5.3 */
/* v.4.3--> */
        }
/* <--v.4.3 */

      /* If any more ligand molecules to pick up, then get identifier
         from XMAS file */
      if (iligand < nlig_molecs - 1)
        Ligand_Molecule_id = get_next_xmas_ligand(xmas,Ligand_Molecule_id);
    }
/* <--v.4.2 */

  printf("\nProgram complete\n\n");

  if (Nwarnings == 1)
    printf("+++ There was 1 warning or advisory message\n\n");
  else if (Nwarnings > 1)
    printf("+++ There were %d warning or advisory messages\n\n",Nwarnings);

  return(0);
}

