C**************************************************************************
C
C  FCOMMON.INC   -   Include file for FCOMMON.F
C
C----------------------------------------------------------------------+---

      CHARACTER*4   TMPLCH
      PARAMETER    ( 
     -              TMPLCH =           "abcd"
     -             )

      INTEGER       FNAMLN, NOT_FOUND
      PARAMETER    ( 
     -              FNAMLN =           512,
     -              NOT_FOUND =         -1
     -             )

      INTEGER       NPDB_DIR, PDB_CURRENT, PDB_OBSOLETE, PDB_MODELS,
     -              PDB_MODELS_OBSOLETE, PDB_MISC, PDB_UPLOAD, PDB_MCSG,
     -              AF_PDB
      PARAMETER    ( 
     -              NPDB_DIR =            8, 
     -              PDB_CURRENT =         1,
     -              PDB_OBSOLETE =        2,
     -              PDB_MODELS =          3, 
     -              PDB_MODELS_OBSOLETE = 4,
     -              PDB_MISC =            5, 
     -              PDB_UPLOAD =          6, 
     -              PDB_MCSG =            7,
     -              AF_PDB =              8
     -             )

      INTEGER       NPDBSUM_DIR, PDBSUM_CURRENT, PDBSUM_UPLOAD
      PARAMETER    ( 
     -              NPDBSUM_DIR =         2, 
     -              PDBSUM_CURRENT =      1, 
     -              PDBSUM_UPLOAD =       2
     -             )

