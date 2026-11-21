C**************************************************************************
C
C  PDBHTML.INC   -   Include file for PDBHTML.F
C
C----------------------------------------------------------------------+---

CVAX      IMPLICIT NONE

      INCLUDE 'fcommon.inc'

      CHARACTER*1   DIRSEP
      PARAMETER    (DIRSEP = '/')

      INTEGER       CHN, LIG, NUC
      PARAMETER    (CHN = 1, LIG = 2, NUC = 3)

      INTEGER       DBASE
      PARAMETER    (DBASE = 21)

      INTEGER       FILEN, MAXCOD, MAXDOM, MAXLEN, MAXPAT, MXCATH,
     -              MXCOBJ, MXENTR, MXEXTR, MXITEM, MXLINE, MXLINK,
     -              MXLNAM, MXMETP, MXNONS, MXNTRY, MXPCOD, MXPRIN,
     -              MXSCRP, MXSITE, MXSPEC, MXSWCH, MXSWIS, MXTITL,
     -              NAPPEN, NLEVEL, NUMLEV, PDILEN, RECLEN
      PARAMETER    (
     -              FILEN  =    512,
     -              MAXCOD =     10,
     -              MAXDOM =    100,
     -              MAXLEN =   1000,
     -              MAXPAT =    100,
     -              MXCATH = 1000000,
     -              MXCOBJ =      8,
     -              MXENTR =     10,
     -              MXEXTR =     25,
     -              MXITEM =    700,
     -              MXLINE =     60,
     -              MXLINK =     60,
     -              MXLNAM =     78,
     -              MXMETP =     10,
     -              MXNONS =    500,
     -              MXNTRY =     99,
     -              MXPCOD =   5000,
     -              MXPRIN =    500,
     -              MXSCRP =  20000,
     -              MXSPEC =      5,
     -              MXSITE =    500,
     -              MXSWCH =    100,
     -              MXSWIS =    100,
     -              MXTITL =  10000,
     -              NAPPEN =      3,
     -              NLEVEL =      4,
     -              NUMLEV =      4,
     -              PDILEN =     70,
     -              RECLEN =     80
     -             )

      CHARACTER*1   CHAIN(MXNTRY,3), SWICHN(MXSWCH,MXSWIS)
      CHARACTER*2   CALPHA(MXNTRY)
      CHARACTER*3   ATTNAM(MXNONS,MXITEM,NAPPEN,3), PATRES(2,MAXPAT)
      CHARACTER*4   BRCODE, EXTRA(MXEXTR), LNKCOD(MXLINK),
     -              METYPE(MXMETP), PRNPDB(MXPCOD),
     -              SPRSDE(MXITEM)
      CHARACTER*5   DOMPDB(MXCATH), PATRNO(2,MAXPAT)
      CHARACTER*7   CATPDB(MXCATH)
      CHARACTER*8   CNUMB, NDBREF
      CHARACTER*10  CHNCOL(MXNTRY), NUCCOL(MXNTRY), SPEGIF(MXSPEC)
      CHARACTER*11  DATE, RLDATE
      CHARACTER*12  ECTMPL, EENTRY(MXENTR), SWNAME(MXSWIS)
      CHARACTER*20  CATNEW, PDBPRE, PRCODE(MXPRIN)
      CHARACTER*30  ECNUMB
      CHARACTER*40  NODESC(MXNONS)
      CHARACTER*50  CATDSC(MXTITL)
      CHARACTER*51  NOFORM(MXNONS)
      CHARACTER*60  WWWGENFIGS
      CHARACTER*110 DOHEAD, DOINTR, DOLRES, DOMOLS, DONUCP, DOPROS,
     -              DOSAS, DOSPEE, PATCOL(MAXPAT), PATSEQ(MAXPAT)
      CHARACTER*120 CATFIL, DOMFIL, ENZPTH, EXPENZ, FNAME, IMBJEN,
     -              KYOENZ, LNKFIL, OBSFIL, PATFIL, PDB3DB, PDBENT,
     -              PDBLNK, PDBNDB, PRNFIL, RUNLIG, RUNWIR, RVLFIL,
     -              SPEFIL, SWADDR, SWIFIL, TITFIL, TOFFIL, UNPLST,
     -              UNPUBL, WITENZ, WWWCAT, WWWCDM, WWWDHS, WWWDOM,
     -              WWWENZ, WWWFSP, WWWGET, WWWHET, WWWMAC, WWWMMD,
     -              WWWNDB, WWWPDB, WWWPFA, WWWPRN, WWWPRO, WWWREP,
     -              WWWRSC, WWWSCP, WWWSPE, WWWTOP, WWWYEA
      CHARACTER*(FNAMLN) MCSLST, NMRDIR, PDBFIL, PDBSUM
      CHARACTER*(MXLNAM) LNAME, SEQUEN(MXLINE,MXNTRY,3)
      CHARACTER*(PDILEN) AUTHOR(MXITEM), COMPND(MXITEM), EXPDTA(MXITEM),
     -              HEADER, OUTREC, RELATE(MXITEM), REM5(MXITEM),
     -              SOURCE(MXITEM), TITLE(MXITEM)

      COMMON/CHARS/ CHAIN, SWICHN,

     -              CALPHA,

     -              ATTNAM, PATRES,

     -              BRCODE, EXTRA, LNKCOD,
     -              METYPE, PRNPDB,
     -              SPRSDE,

     -              DOMPDB, PATRNO,

     -              CATPDB,

     -              CNUMB, NDBREF,

     -              CHNCOL, NUCCOL, SPEGIF,

     -              DATE, RLDATE,

     -              ECTMPL, EENTRY, SWNAME,

     -              CATNEW,

     -              PDBPRE, PRCODE,

     -              ECNUMB,

     -              NODESC,

     -              CATDSC,

     -              NOFORM,

     -              WWWGENFIGS,

     -              DOHEAD, DOINTR, DOLRES, DOMOLS, DONUCP, DOPROS,
     -              DOSAS, DOSPEE, PATCOL, PATSEQ,

     -              CATFIL, DOMFIL, ENZPTH, EXPENZ, FNAME, IMBJEN,
     -              KYOENZ, LNKFIL, OBSFIL, PATFIL, PDB3DB, PDBENT,
     -              PDBLNK, PDBNDB, PRNFIL, RUNLIG, RUNWIR, RVLFIL,
     -              SPEFIL, SWADDR, SWIFIL, TITFIL, TOFFIL, UNPLST,
     -              UNPUBL, WITENZ, WWWCAT, WWWCDM, WWWDHS, WWWDOM,
     -              WWWENZ, WWWFSP, WWWGET, WWWHET, WWWMAC, WWWMMD,
     -              WWWNDB, WWWPDB, WWWPFA, WWWPRN, WWWPRO, WWWREP,
     -              WWWRSC, WWWSCP, WWWSPE, WWWTOP, WWWYEA,

     -              MCSLST, NMRDIR, PDBFIL, PDBSUM,

     -              LNAME, SEQUEN,

     -              AUTHOR, COMPND, EXPDTA,
     -              HEADER, OUTREC, RELATE, REM5,
     -              SOURCE, TITLE

      INTEGER       ACOUNT(MXNONS,MXITEM,NAPPEN,3),
     -              CATLEV(NUMLEV,MXTITL), CATNUM(NUMLEV,MXCATH),
     -              CLENTH(3), DOMNUM(MXCATH), DUPLIC(MXNTRY),
     -              METCOL(MXMETP), NALINE, NAMLEN(MXENTR), NATOMS,
     -              NBONDS(MXNTRY), NCATH, NCHAIN, NCLINE, NCODES,
     -              NDIFF, NDOMS, NEC, NELINE, NEXTRA, NHELIX(MXNTRY),
     -              NHYDRO, NLIGAN, NLINKS, NMETAL, NMETYP(MXMETP),
     -              NMODEL, NNUCL, NONUMB(MXNONS), NOTYPE, NPCODE,
     -              NPRLNK, NRES, NRESID(MXNTRY,3), NRLINE, NSLINE,
     -              NSPECI, NSTRND(MXNTRY), NSUPER, NSWCHN(MXSWCH),
     -              NSWISS, NTIMES(MXNTRY,3), NTITLE, NTLINE,
     -              NUMEXT(MXEXTR), NUMLEN(4), NUNIQ, NWATER, NXLINE,
     -              NXTCHN(MXNTRY), NXTCOL, PATOCC(MAXPAT),
     -              PRNUM(MXPCOD), RERROR

      COMMON/INTS/  ACOUNT,
     -              CATLEV, CATNUM,
     -              CLENTH, DOMNUM, DUPLIC,
     -              METCOL, NALINE, NAMLEN, NATOMS,
     -              NBONDS, NCATH, NCHAIN, NCLINE, NCODES,
     -              NDIFF, NDOMS, NEC, NELINE, NEXTRA, NHELIX,
     -              NHYDRO, NLIGAN, NLINKS, NMETAL, NMETYP,
     -              NMODEL, NNUCL, NONUMB, NOTYPE, NPCODE,
     -              NPRLNK, NRES, NRESID, NRLINE, NSLINE,
     -              NSPECI, NSTRND, NSUPER, NSWCHN,
     -              NSWISS, NTIMES, NTITLE, NTLINE,
     -              NUMEXT, NUMLEN, NUNIQ, NWATER, NXLINE,
     -              NXTCHN, NXTCOL, PATOCC,
     -              PRNUM, RERROR

      LOGICAL       ECEXPT, EMEXPT, ENSLNK, HAVCHK, HAVCLF, HAVIEW,
     -              HAVJPG, HAVLIG, HAVPRO, HAVRSC, HAVSRF, HAVUNP,
     -              HAVWRL, IFAIL, INCSAS, JAVASC, LAYER1, MCSG, NOLINK,
     -              NMR, NMRLNK, OBSOL, PDBSM1, PLAIN, RELINK, SFOPEN,
     -              TMODEL, TOOBIG, UPLOAD, XRAY, YEAST
      COMMON/LOGS/  ECEXPT, EMEXPT, ENSLNK, HAVCHK, HAVCLF, HAVIEW,
     -              HAVJPG, HAVLIG, HAVPRO, HAVRSC, HAVSRF, HAVUNP,
     -              HAVWRL, IFAIL, INCSAS, JAVASC, LAYER1, MCSG, NOLINK,
     -              NMR, NMRLNK, OBSOL, PDBSM1, PLAIN, RELINK, SFOPEN,
     -              TMODEL, TOOBIG, UPLOAD, XRAY, YEAST

      REAL          RESOL, RFACT, RFREE
      COMMON/REALS/ RESOL, RFACT, RFREE

C----------------------------------------------------------------------+---
C
C VARIABLE DESCRIPTIONS:-
C ---------------------
C
C ACOUNT(MXNONS,MXITEM,NAPPEN,3) - Count of attached HET groups /
C          non-standard residues
C ATTNAM(MXNONS,MXITEM,NAPPEN,3) - Het groups and non-standard residues
C AUTHOR(MXITEM) - Author line(s) in PDB file
C BRCODE - 4-character Brookhaven code
C CALPHA(MXNTRY) - Indicator of whether the chain is CA-only
C CATDSC(MXTITL) - CATH title
C CATFIL - Full path to file holding CATH number for each chain, cath.numbers
C CATH   - CATH number for given chain
C CATLEV(NUNLEV,MXTITL) - Entry number for current CATH level
C CATNEW - Compressed CATH number
C CATNUM(NUMLEV,MXCATH) - CATH number for each chain and domain of
C          each PDB entry
C CATPDB(MXCATH) - PDB code, chain and domain with equivalent CATH number
C          in CATNUM.
C CDESC  - CATH Class description
C CHAIN(MXNTRY,3) - Chain IDs for chains and ligands
C CHN    - Indicator for chain data
C CHNCOL(MXNTRY) - Colour for each chain in chain-trace
C CLENTH(3) - Length of string corresponding to E.C. numbers within the
C          same class or subclass, from the first class (length 1) downwards
C COMPND(MXITEM) - Molecule description from COMPND records in the PDB file
C DATE   - Date of PDB entry
C DBASE  - Unit name for database.dat output file
C DLEN   - Length of description
C DOHEAD - Command to run perl script for showing PDB header information
C DOINTR - Command to run perl script for showing intron splice sites
C DOLRES - Command to run perl script for loading RasMol with the given
C          ligand and its interacting residues
C DOMFIL - Full path to file holding CATH domains for each chain, cath.domains
C DOMOLS - Command to run perl script for generating a MolScript picture
C          of the given chain
C DOMPDB(MXCATH) - PDB code and chain for each entry in cath.domains
C DONUCP - Command to run perl script for displaying gzipped NUCPLOT
C          PostScript file
C DOPROS - Command to run perl script for generating the RasMol script
C          for a given PROSITE pattern
C DOSAS  - Command to run perl script for generating the SAS annotated
C          alignments
C DOSPEE - Command to run perl script for generating the binding site
C          analysis
C DUPLIC(MXNTRY) - Unique chain, set up by routine CHUNIQ. Duplicate
C          chains have a zero here and a chain-list number in NXTCHN
C ECEXPT - Flag indicating that expt is electron crystallography
C EMEXPT - Flag indicating that expt is electron microscopy
C ECNUMB - Uncompressed E.C. number
C ECOMPR - Compressed E.C. number
C ECTMPL - Template for blank E.C. numbers
C EENTRY(MXENTR) - Standardised set of E.C. numbers for current PDB file
C ENSLNK - Flag indicating that linked NMR PDB code holds the NMR ensemble
C ENZPTH - URL path to appropriate entry in Enzyme Structures Database
C EXPDTA(MXITEM) - Experimental methods
C EXPENZ - Link to ExPaSy ENZYME database entry for given enzyme code
C EXTRA(MXEXTR) - Unidentified single-atom residue name
C HAVCHK - Flag indicating that have PROCHECK results
C HAVCLF - Flag indicating that have clefts image file
C HAVIEW - Flag indicating that have both side-views of our structure
C HAVJPG - Flag indicating that have JPEG versions of trace files
C HAVLIG - Flag indicating that have ligands
C HAVPRO - Flag indicating that have PROMOTIF analyses
C HAVSRF - Flag indicating that have a SURFPLOT image
C HAVRSC - Flag indicating that can generate RasMol script
C HAVUNP - Flag indicating that have link to unpublished references
C          database
C HAVWRL - Flag indicating that have a VRML file
C IFAIL  - Error flag
C IMBJEN - Link to IMB Jena Image Library
C JAVASC - Flag indicating that JavaScript lines to be included
C KYOENZ - Link to Kyoto University database entry for given enzyme code
C LAYER1 - Flag indicating that this is a "Layer 1" release from the PDB
C          and, as such, may contain severe irregularities.
C LIG    - Indicator for ligand data
C LINE   - Line-count for records read from input file
C LNKCOD(MXLINK) - PDB codes of related NMR entry containing ensemble
C          or minimized structure
C MCSG   - Flag indicating whether structure comes from the Midwest
C          Center for Structural Genomics (MCSG)
C MCSGDR - Path to PDB files from Mid-West Center for Structural Genomics
C MCSLST - File listing all the MCSG structures deposited in the PDB
C METCOL(MXMETP) - Colour assigned to a given metal type
C METYPE(MXMETP) - Metal types
C MAXCOD - Maximum number of E.C. codes associated with each PDB file
C MAXLEN - Maximum length of output line
C MAXPAT - Maximum number of PROSITE patterns that can match a given chain
C MXCATH - Maximum number of CATH numbers that can be stored
C MXCOBJ - Maximum number of different colours for the different objects
C MXENTR - Maximum number of E.C. entries for a given PDB file
C MXEXTR - Maximum number of different types of unidantified single-atom
C          residues
C MXITEM - Maximum number of lines per given entry type
C MXLINE - Maximum number of sequence-info lines that can be stored for
C          each entry
C MXLINK - Maximum number of linked PDB codes
C MXLNAM - Maximum length of ligand name
C MXMETP - Maximum number of metal types
C MXNONS - Maximum number of non-standard residue-types per chain
C MXNTRY - Maximum number of chains and ligands
C MXSCRP - Maximum number of atom records for which link to Rasmol script
C          will be written
C MXSITE - Maximum number of SITE details that can be picked up from the
C          sites.html file
C MXSPEC - Maximum number of species entries per PDB file
C MXSWIS - Maximum number of SWISS-PROT codes per PDB entry
C MXTITL - Maximum number of CATH title entries
C NALINE - Number of author lines read in from the PDB header details
C NAMLEN(MXNTRY) - Name length of each E.C. entry for current PDB file
C NATOMS - Total number of atoms and hetatoms in the original PDB file
C NBONDS(MXNTRY) - Number of bonds in the ligand. If less than two, then
C          no LIGPLOT will have been generated for this ligand
C NCATH  - Number of CATH numbers read in
C NCHAIN - Number of chains in the file
C NCLINE - Number of lines of entries stored in COMPND array
C NCODES - Number of E.C. numbers associated with given PDB file
C NDBREF - NDB reference identifier for current structure
C NDIFF  - Number of different PROSITE patterns
C NEC    - Number of E.C. numbers for the given PDB file
C NELINE - Number of related entry lines read in from the PDB header details
C NEXTRA - Number of unidentified single-atom residues
C NHELIX(MXNTRY) - Number of helices in chain
C NLEVEL - Number of top CATH levels for which titles are required when
C          showing a given chain's CATH number
C NLIGAN - Number of ligands
C NLINKS - Number of linked PDB structures to current one
C NMETAL - Number of metal ions
C NMETYP(MXMETP) - Number of occurrences of each metal type in the protein
C NMODEL - Number of NMR models in the PDB file
C NMR    - Flag indicating that the structure is an NMR structure
C NMRLNK - Flag indicating that the structure has a related PDB entry
C          giving either the whole ensemble or the minimized structure.
C NNUCL  - Number of nucleotide chains in the PDB file
C NODESC(MXNONS) - Descriptions of non-standard residue names found in
C          PDB file
C NOFORM(MXNONS) - Formulae for non-standard residue names found in
C          PDB file
C NOLINK - Flag indicating that no links to external databases to be
C          written out as PDB file is local only
C NONUMB(MXNONS) - Numbers of occurrences of each of the different types
C          of non-standard residue names found in PDB file
C NOTYPE - Number of types of non-standard residue names found in PDB file
C NPCODE - Number of PDB entries which have PRINTS codes
C NPRLNK - Number of links to PRINTS entries
C NRES   - Number of residues read in
C NRESID(MXNTRY,3) - Number of residues for each chain/ligand
C NSPECI - Number of images associated with species found
C NSTRND(MXNTRY) - Number of strands in chain
C NSUPER - Number of superseded/ing entries
C NSWCHN(MXSWCH) - Number of chain entries associated with each SWISS-PROT
C          code
C NSWISS - Number of SWISS-PROT entries for given protein
C NTIMES(MXNTRY,3) - Number of occurrences of a given ligand
C NTITLE - Number of CATH title entries read in
C NTLINE - Number of lines of entries stored in TITLE array
C NUC    - Indicator for nucleotide chain
C NUCCOL(MXNTRY) - Colour for backbone of nucleotide sequence
C NUMEXT(MXEXTR) - Number of each unidentified single-atom residue type
C NUMLEN(4) - Length of each number in the 4 positions of the standard E.C.
C          number
C NUMLEV - Number of levels of the CATH number to be shown
C NUNIQ  - Number of unique protein chains
C NXTCHN(MXNTRY) - Next duplicate chain, set up by routine CHUNIQ
C NXLINE - Number of lines on experimental methods
C NXTCOL - Next available colour for metals
C OBSFIL - Name of file containing all the supersede data
C OBSOL  - Flag indicating that current entry is the obsolete one
C PATCOL(MAXPAT) - Residue colours for the given PROSITE pattern
C PATFIL - Full path to file holding PROSITE pattern-matches database
C PATOCC(MAXPAT) - Pattern occurrence number for the given PROSITE pattern
C PATRES(2,MAXPAT) - Start- and end-residues for given PROSITE pattern
C PATRNO(2,MAXPAT) - Start- and end-residue numbers for given PROSITE pattern
C PATSEQ(MAXPAT) - Residue sequence matching given PROSITE pattern
C PDBENT - Full path and filename for pdbent.dat file containing E.C.
C          number for each PDB file
C PDBFIL - Full name of PDB file
C PDBLNK - Link to entry for given PDB code on PDB server
C PDBNDB - File containing NDB reference nos for each PDB entry
C PDBSM1 - Flag indicating running for PDBsum1
C PDBSUM - Entry's PDBsum directory
C PDB3DB - Link to 3DB entry at PDB
C PDILEN - Maximum length of text information held on a PDB record
C PLAIN  - Flag indicating that output gif file requires no link and
C          no number
C PRCODE(MXPRIN) - PRINTS codes
C PRNPDB(MXPCOD) - PDB codes for which there are PRINTS entries
C PRNUM(MXPCOD) - Link to relevant PRINTS code from this PDB entry
C PRNFIL - File holding PRINTS codes for known PDB files
C PRODIR - Path to directory containing PDBsum files
C RECLEN - Maximum input length of a PDB record
C RELATE(MXITEM) - Lines of text dealing with related PDB entries
C RELINK - Flag indicating that linked PDB codes are "related" entries,
C          rather than links between NMR ensembles and minimized structures
C RERROR - Read-error indicator
C RFACT  - R-factor, where given
C RFREE  - R-free, where given
C RLDATE - PDB entry's release date
C RUNLIG - Address of CGI-bin script for running LIGPLOT to get a
C          PostScript version of the LIGPLOT diagram
C RUNWIR - Address of CGI-bin script for running WIRPLOT to get a
C          PostScript version of the secondary structure wiring diagram
C RVLFIL - R-factors file
C SEQUEN(MXLINE,MXNTRY,3) - Sequence for each chain and ligand
C SFOPEN - Flag indicating whether the superseded entries file is open
C SOURCE - Source of protein
C SPEFIL - Name of file containing all the species GIF file names
C SPEGIF(MXSPEC) - Names of GIF images holding images of species
C SPRSDE(MXITEM) - PDB code that supersedes current entry, or list of PDB
C          codes superseded by it
C SWADDR - Server address for Swiss Prot database
C SWICHN(MXSCHN,MXSWIS) - Chains relating to each Swiss-prot code
C SWIFIL - Filename for file holding the SwissProt links
C SWNAME(MXSWIS) - Swiss-prot codes for given entry
C TITFIL - Full path to file holding CATH titles
C TITLE(MXITEM) - Title of PDB entry (in TITLE record)
C TOFFIL - Full path to file holding PROSITE true/false matches file
C TMODEL - Flag indicating whether PDB structure is a theoretical model
C TOOBIG - Flag indicating that structure is too big for PDB viewers
C UNPLST - List containing entries in Stuart Moodie's Unpublished
C          references database
C UNPUBL - Link to Stuart Moodie's Unpublished references database
C WITENZ - Link to the What Is There database entry for given enzyme code
C WWWCAT - WWW reference to the CATH database
C WWWCDM - WWW reference to given CATH domaion
C WWWDHS - WWW reference to DHS database for a given CATH code
C WWWDOM - WWW reference to Perl script for generating Rasmol output
C          to shown protein/chain coloured by CATH domain definitions
C WWWENZ - WWW reference to the Enzyme database
C WWWFSP - WWW reference to the FSSP database of structurally aligned
C          protein fold families
C WWWGET - WWW reference to Perl script for getting PDB coords and
C          displaying using Rasmol
C WWWHET - WWW reference to the HET groups listing
C WWWMAC - WWW reference to the Macromolecule database at the EBI
C WWWMMD - WWW reference to the MMDB database (Entrez)
C WWWNDB - WWW reference to the NDB entry search page
C WWWPDB - WWW reference to the PDB entry search page
C WWWPFA - WWW reference to the PFam database
C WWWPRN - WWW reference to the PRINTS code Perl script
C WWWPRO - WWW reference to the PROSITE database
C WWWREP - WWW reference to PDBREPORT
C WWWRSC - WWW reference to Perl script for generating Rasmol output with
C          script defined in rasmol.scr
C WWWSCP - WWW reference to the SCOP classification database
C WWWSPE - WWW reference to species GIF images
C WWWTOP - WWW reference to TOPs diagrams at the EBI
C WWWYEA - WWW reference to the Yeast Structures database at Stanford
C XRAY   - Flag indicating that expt is X-ray experiment
C YEAST  - Flag indicating whether this is a yeast structure
C
C--------------------------------------------------------------------------
