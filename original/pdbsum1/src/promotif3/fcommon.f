C**************************************************************************
C
C  FCOMMON.F  -  Common FORTRAN subroutines
C
C                Written by Roman Laskowski, EBI,
C                November 2012.
C
C--------------------------------------------------------------------------
C
C Compilation (on unix)
C ---------------------
C
C gfortran -fbounds-check -c fcommon.f
C
C--------------------------------------------------------------------------
C
C Subroutines
C -----------
C
C   --> GETPDR  --> GETPAR (see below)
C   --> GETPSM  --> GETPAR (see below)
C   --> GETPAR  --> GETDET  --> INTERP
C                           --> EXPTOK  --> LOCTOK
C               --> GETVAL
C   --> PDBNAM  --> FRMFNM
C   --> SUMNAM  --> FRMFNM
C   --> FRMNUM
C   --> LOWCAS
C   --> UPPCAS
C   --> LENSTR
C
C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETPDR  -  Get the set of PDB directory templates
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETPDR(PDBDIR)

      INCLUDE 'fcommon.inc'

      CHARACTER*(FNAMLN) PDBDIR(NPDB_DIR)

      INTEGER       IPDB, LENSTR

      CHARACTER*(FNAMLN) TEMPLT(NPDB_DIR)

      DATA TEMPLT /
     -     'FTP_PDB_TEMPLATE', 'FTP_OBSOLETE_TEMPLATE',
     -     'FTP_MODELS_TEMPLATE', 'FTP_OBSOLETE_MODELS_TEMPLATE',
     -     'MISC_PDB_TEMPLATE', 'UPLOAD_PDB_TEMPLATE',
     -     'MCSG_PDB_TEMPLATE', 'ALPHAFOLD_PDB_TEMPLATE' /

C---- Loop over the PDB file templates
      DO 800, IPDB = 1, NPDB_DIR

C----     Get this PDB directory
          CALL GETPAR(TEMPLT(IPDB),PDBDIR(IPDB),FNAMLN)

 800  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETPSM  -  Get the set of PDBsum directory templates
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETPSM(PSMDIR)

      INCLUDE 'fcommon.inc'

      CHARACTER*(FNAMLN) PSMDIR(NPDBSUM_DIR)

      INTEGER       ISUM, LENSTR

      CHARACTER*(FNAMLN) TEMPLT(NPDBSUM_DIR)

      DATA TEMPLT /
     -     'PDBSUM_REAL_TEMPLATE', 'UPLOAD_PDBSUM_TEMPLATE'
     -     /

C---- Loop over the PDBsum templates
      DO 800, ISUM = 1, NPDBSUM_DIR

C----     Get this PDBsum directory
          CALL GETPAR(TEMPLT(ISUM),PSMDIR(ISUM),FNAMLN)

 800  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETPAR  -  Get the required file parameter
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETPAR(KEY,FILNAM,FILEN)

      SAVE

      INTEGER       MAXPAR, STRLEN
      PARAMETER    (MAXPAR = 5000, STRLEN = 256)

      CHARACTER*(*) KEY, FILNAM
      CHARACTER*(STRLEN) TOKNAM(MAXPAR), TOKVAL(MAXPAR)
      INTEGER       FILEN, LENSTR, NTOKEN
      LOGICAL       FIRST, HAVCAT, IFAIL

      DATA FIRST   / .TRUE. /
      DATA HAVCAT  / .FALSE. /

C---- Initialise variables
      FILNAM = ' '
      IFAIL = .FALSE.

C---- If this is the first call, then read in all the parameters from
C     the CATHPARAM file
      IF (FIRST) THEN

C----     Reset the flag
          FIRST = .FALSE.

C----     Read in all the parameters from the CATHPARAM file
          CALL GETDET(NTOKEN,TOKNAM,TOKVAL,STRLEN,MAXPAR,IFAIL)
          IF (IFAIL) GO TO 999

C----     Set flag indicating that we have the CATHPARAM file
          HAVCAT = .TRUE.
      ENDIF

C---- Get the value corresponding to the given token
      IF (HAVCAT) THEN
          CALL GETVAL(KEY,FILNAM,FILEN,NTOKEN,TOKNAM,TOKVAL,STRLEN,
     -        MAXPAR)

      ENDIF

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETDET  -  Read in all the parameters from the CATHPARAM file
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETDET(NTOKEN,TOKNAM,TOKVAL,STRLEN,MAXPAR,IFAIL)

      INCLUDE 'fcommon.inc'

      INTEGER       LINLEN, MAXNAM
      PARAMETER    (LINLEN = 1000,MAXNAM = 8)

      INTEGER       MAXPAR, STRLEN

      CHARACTER*1   TOKNAM(STRLEN,MAXPAR), TOKVAL(STRLEN,MAXPAR)
      CHARACTER*100 PARNAM(MAXNAM)
      CHARACTER*512 CTHPRM, FILNAM
      CHARACTER*(LINLEN) IREC
      INTEGER       IFILE, LEN, LENSTR, LINE, NTOKEN, STATUS
      LOGICAL       HAVFIL, IFAIL
CDEBUG
C      CHARACTER*256 STRNAM, STRVAL
C      INTEGER       IPOS
C      LOGICAL       FOUND
CDEBUG

      DATA PARNAM /
     -    'CATHPARAM',
     -    '../../param/CATHPARAM',
     -    './param/CATHPARAM',
     -    '../cath/param/CATHPARAM',
     -    '../cath/param/CATHPARAM',
     -    '/nfs/httpd/cgi-bin/cath/param/CATHPARAM',
     -    '/usr/local/httpd/cgi-bin/cath/param/CATHPARAM',
     -    '/ebi/research/thornton/www/databases/cgi-bin/param/CATHPARAM'
     -    /

C---- Initialise variables
      HAVFIL = .FALSE.
      LINE = 0
      NTOKEN = 0

C---- Get the location of the CATHPARAM file from the CATHPARAM
C     environment variable
      CALL GET_ENVIRONMENT_VARIABLE('CATHPARAM',CTHPRM,LEN,STATUS,
     -     .TRUE.)

C---- Loop while trying to locate the CATHPARAM file
      DO 300, IFILE = 0, MAXNAM

         IF (IFILE.EQ.0) THEN
            FILNAM = CTHPRM
         ELSE
            FILNAM = PARNAM(IFILE)
         ENDIF
C----     Open the parameter file
CDEBUG
C      INQUIRE(FILE=FILNAM,EXIST=FOUND)
C      PRINT*, 'Trying CATHPARAM: [',
C     -    FILNAM(1:LENSTR(FILNAM)), '] FOUND=', FOUND
CDEBUG
          OPEN(UNIT=1, FILE=FILNAM, STATUS='OLD',
     -        FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -        READONLY,
     -        ERR=300)

C----     File opened
          HAVFIL = .TRUE.
          GO TO 500
 300  CONTINUE

C---- If file not found, display error message
      IF (.NOT.HAVFIL) THEN
C          PRINT*, '*** ERROR. CATHPARAM file not found'
          IFAIL = .TRUE.
          GO TO 999
      ENDIF

C---- Read through the file, storing all the tokens and their translations
 500  CONTINUE

C----     Read in the next line from the file
          READ(1,520,END=800,ERR=800) IREC
 520      FORMAT(A)

C----     Get the line-length
          LEN = LENSTR(IREC)

C----     Process if line is not blank and is not a comment line
          IF (LEN.GT.1 .AND. IREC(1:1).NE.'#') THEN

C----         Interpret and expand the tokens on this line
              CALL INTERP(IREC,LEN,TOKNAM,TOKVAL,MAXPAR,STRLEN,NTOKEN)

C----         Expand token value if it is a compound name involving a
C             prior token value
              CALL EXPTOK(TOKNAM,TOKVAL,MAXPAR,STRLEN,NTOKEN)
CDEBUG
C              STRNAM = ' '
C              STRVAL = ' '
C              DO 201, IPOS = 1, 256
C                  IF (TOKNAM(IPOS,NTOKEN).NE.' ')
C     -                STRNAM(IPOS:IPOS) = TOKNAM(IPOS,NTOKEN)
C                  IF (TOKVAL(IPOS,NTOKEN).NE.' ')
C     -                STRVAL(IPOS:IPOS) = TOKVAL(IPOS,NTOKEN)
C 201          CONTINUE
C              PRINT*, 'Token ', STRNAM(1:LENSTR(STRNAM)), ' - ',
C     -            STRVAL(1:LENSTR(STRVAL))
CDEBUG

C----     If end of data, can stop here
          ELSE IF (IREC(1:10).EQ.'#---------') THEN
              GO TO 800
          ENDIF

C---- Loop back for next record in file
      GO TO 500

C---- End of file reached
 800  CONTINUE

C---- Close the file
      CLOSE(1)

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE INTERP  -  Interpret and expand the tokens on this line
C
C----------------------------------------------------------------------+---

      SUBROUTINE INTERP(IREC,LEN,TOKNAM,TOKVAL,MAXPAR,STRLEN,NTOKEN)

      INTEGER       MAXPAR, STRLEN

      CHARACTER*1   CH, TABCH, TOKNAM(STRLEN,MAXPAR),
     -              TOKVAL(STRLEN,MAXPAR)
      CHARACTER*(*) IREC
      INTEGER       I, IPOS, LEN, LENNAM, LENVAL, NTOKEN, OUTPOS, STATE
      LOGICAL       DONE, HAVTOK, INNAME, SAVECH

C---- Initialise variables
      DONE = .FALSE.
      HAVTOK = .FALSE.
      INNAME = .TRUE.
      LENNAM = 0
      LENVAL = 0
      OUTPOS = 0
      STATE = 0
      TABCH = CHAR(9)

C---- Loop through all the characters in the line
      DO 1000, IPOS = 1, LEN

C----     Get the character at this position
          CH = IREC(IPOS:IPOS)
          SAVECH = .FALSE.

C----     Process according to current state

C----     State 0: Looking for first non-blank character
          IF (STATE.EQ.0) THEN

C----         If non-blank, then initialise and switch state
              IF (CH.NE.' ' .AND. CH.NE.TABCH) THEN

C----             Initialise for token or token-string
                  OUTPOS = 0
                  SAVECH = .TRUE.
                  STATE = 1

              ENDIF

C----     State 1: Adding characters to current string
          ELSE IF (STATE.EQ.1) THEN

C----         If non-blank, then add to current string
              IF (CH.NE.' ' .AND. CH.NE.TABCH) THEN
                  SAVECH = .TRUE.

C----         Otherwise, have reached end of current token
              ELSE

C----             Terminate string and switch state
                  CH = ' '
                  IF (HAVTOK) THEN
                      DONE = .TRUE.
                  ELSE
                      STATE = 2
                  ENDIF
                  SAVECH = .FALSE.
              ENDIF

C----     State 2: Looking for equals sign
          ELSE IF (STATE.EQ.2) THEN

C----         If this is an equals sign, the switch state
              IF (CH.EQ.'=') THEN
                  STATE = 0
                  INNAME = .FALSE.
              ENDIF
          ENDIF

C----     If saving character, then do so
          IF (SAVECH) THEN

C----         If this is the first character of the name, then
C             increment token-count
              IF (INNAME .AND. OUTPOS.EQ.0) THEN

C----             Increment token count
                  NTOKEN = NTOKEN + 1

C----             Check that maximum not exceeded
                  IF (NTOKEN.GT.MAXPAR) THEN
                      NTOKEN = MAXPAR
                  ENDIF

C----             Initialise token name and value
                  DO 300, I = 1, STRLEN
                      TOKNAM(I,NTOKEN) = ' '
                      TOKVAL(I,NTOKEN) = ' '
 300              CONTINUE
              ENDIF

C----         Increment character position
              OUTPOS = OUTPOS + 1

C----         Check that not off end of string
              IF (OUTPOS.GT.STRLEN) THEN
                  OUTPOS = STRLEN
              ENDIF

C----         Save in token name or token value
              IF (INNAME) THEN
                  TOKNAM(OUTPOS,NTOKEN) = CH
                  LENNAM = OUTPOS
              ELSE
                  TOKVAL(OUTPOS,NTOKEN) = CH
                  LENVAL = OUTPOS
                  HAVTOK = .TRUE.
              ENDIF
          ENDIF

C----     If done, then end here
          IF (DONE) GO TO 999
 1000 CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE EXPTOK  -  Expand token value if it is a compound name
C                        involving a prior token value
C
C----------------------------------------------------------------------+---

      SUBROUTINE EXPTOK(TOKNAM,TOKVAL,MAXPAR,STRLEN,NTOKEN)

      INTEGER       MAXLEN
      PARAMETER    (MAXLEN = 1000)

      INTEGER       MAXPAR, STRLEN

      CHARACTER*1   CH, TOKNAM(STRLEN,MAXPAR), TOKVAL(STRLEN,MAXPAR)
      CHARACTER*(MAXLEN) STRING, TMPSTR, TOKSUB
      INTEGER       I, IEND, IPOS, ISTART, LEN, LENSTR, LENSUB, LENTOK,
     -              NTOKEN

C---- Extract the current token
      STRING = ' '
      LEN = STRLEN
      DO 100, IPOS = 1, STRLEN

C----     Get current character
          CH = TOKVAL(IPOS,NTOKEN)

C----     If not a space, copy across to string
          IF (CH.NE.' ') THEN
              STRING(IPOS:IPOS) = CH
              LEN = IPOS

C----     Otherwise, stop here
          ELSE
              GO TO 200
          ENDIF
 100  CONTINUE

C---- Loop until all token expansion complete
 200  CONTINUE

C----     Check whether expansion is required
          IPOS = INDEX(STRING,'$(')

C----     If no more expansion required, jump out of loop
          IF (IPOS.EQ.0) GO TO 999

C----     Find the token-string to be substituted
          ISTART = IPOS
          TMPSTR = STRING(IPOS:)
          IEND = IPOS + INDEX(TMPSTR,')') - 1

C----     If have valid start- and end-values, perform the substitution
          IF (ISTART.GT.0 .AND. IEND.GT.0) THEN

C----         Get the length of token-name to be looked up
              LENTOK = IEND - ISTART - 2

C----         Splice out the token-name to be substituted
              TOKSUB = TMPSTR(3:LENTOK + 2)

C----         Locate the token and replace its name by its value
              CALL LOCTOK(TOKSUB,LENTOK,TOKNAM,TOKVAL,MAXPAR,STRLEN,
     -            NTOKEN,TMPSTR,LENSUB)

C----         Splice together the new token value
              IF (ISTART.GT.1) THEN
                  TOKSUB = STRING(1:ISTART - 1) // TMPSTR(1:LENSUB) //
     -                STRING(IEND + 1:)
              ELSE
                  TOKSUB = TMPSTR(1:LENSUB) // STRING(IEND + 1:)
              ENDIF
              LEN = LENSTR(TOKSUB)

C----         Initialise token value
              DO 300, I = 1, STRLEN
                  TOKVAL(I,NTOKEN) = ' '
 300          CONTINUE

C----         Store the updated token value
              DO 400, I = 1, LEN
                  TOKVAL(I,NTOKEN) = TOKSUB(I:I)
 400          CONTINUE

C----         Store new string version
              STRING = TOKSUB
          ENDIF

C---- Loop back
      GO TO 200

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE LOCTOK  -  Locate the token and replace its name by its value
C
C----------------------------------------------------------------------+---

      SUBROUTINE LOCTOK(TOKSUB,LENTOK,TOKNAM,TOKVAL,MAXPAR,STRLEN,
     -    NTOKEN,TMPSTR,LEN)

      INTEGER       MAXPAR, STRLEN

      CHARACTER*(*) TMPSTR, TOKSUB
      CHARACTER*1   CH, TOKNAM(STRLEN,MAXPAR), TOKVAL(STRLEN,MAXPAR)
      INTEGER       IPOS, ITOKEN, LEN, LENTOK, NTOKEN
      LOGICAL       DONE

C---- Initialise variables
      DONE = .FALSE.
      LEN = 0
      TMPSTR = ' '

C---- Loop through all the tokens until get a match
      DO 800, ITOKEN = 1, NTOKEN

C----     Extract the current token name
          TMPSTR = ' '
          LEN = STRLEN
          DO 100, IPOS = 1, STRLEN

C----         Get current character
              CH = TOKNAM(IPOS,ITOKEN)

C----         If not a space, copy across to string
              IF (CH.NE.' ') THEN
                  TMPSTR(IPOS:IPOS) = CH
                  LEN = IPOS

C----         Otherwise, stop here
              ELSE
                  GO TO 200
              ENDIF
 100      CONTINUE

C----     Token name extracted
 200      CONTINUE

C----     If the token to be substituted equals the current name, then
C         have found what we're looking for
          IF (TMPSTR(1:LEN).EQ.TOKSUB(1:LENTOK)) THEN

C----         Extract the corresponding token value
              TMPSTR = ' '
              LEN = STRLEN
              DO 400, IPOS = 1, STRLEN

C----             Get current character
                  CH = TOKVAL(IPOS,ITOKEN)

C----             If not a space, copy across to string
                  IF (CH.NE.' ') THEN
                      TMPSTR(IPOS:IPOS) = CH
                      LEN = IPOS

C----             Otherwise, stop here
                  ELSE
                      GO TO 999
                  ENDIF
 400          CONTINUE

C----         Jump out of loop
              GO TO 999
          ENDIF
 800  CONTINUE

C---- If have found a match, then blank out return value
      TMPSTR = ' '
      LEN = 0

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETVAL  -  Get the value corresponding to the given token
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETVAL(KEY,FILNAM,FILEN,NTOKEN,TOKNAM,TOKVAL,STRLEN,
     -    MAXPAR)

      INTEGER       MAXLEN
      PARAMETER    (MAXLEN = 1000)

      INTEGER       MAXPAR, STRLEN

      CHARACTER*1   CH, TOKNAM(STRLEN,MAXPAR), TOKVAL(STRLEN,MAXPAR)
      CHARACTER*(MAXLEN) STRING
      CHARACTER*(*) KEY, FILNAM
      INTEGER       FILEN, IPOS, ITOKEN, LEN, LENSTR, LENTOK, NTOKEN

C---- Initialise variables
      DO 50, IPOS = 1, FILEN
          FILNAM(IPOS:IPOS) = ' '
 50   CONTINUE

C---- Get the length of the search key
      LEN = LENSTR(KEY)

C---- Loop through all the tokens until get a match
      DO 800, ITOKEN = 1, NTOKEN

C----     Extract the current token name
          STRING = ' '
          LENTOK = STRLEN
          DO 100, IPOS = 1, STRLEN

C----         Get current character
              CH = TOKNAM(IPOS,ITOKEN)

C----         If not a space, copy across to string
              IF (CH.NE.' ') THEN
                  STRING(IPOS:IPOS) = CH
                  LENTOK = IPOS

C----         Otherwise, stop here
              ELSE
                  GO TO 200
              ENDIF
 100      CONTINUE

C----     Token name extracted
 200      CONTINUE

C----     If the search key matches the token name, then have found
C         what we're after
          IF (STRING(1:LENTOK).EQ.KEY(1:LEN)) THEN

C----         Determine length to be transferred
              LEN = MIN(FILEN,STRLEN)

C----         Transfer token value across to answer-field
              DO 400, IPOS = 1, LEN

C----             If non-blank, then copy across
                  IF (TOKVAL(IPOS,ITOKEN).NE.' ') THEN
                      FILNAM(IPOS:IPOS) = TOKVAL(IPOS,ITOKEN)

C----             Otherwise, jump out
                  ELSE
                      GO TO 999
                  ENDIF
 400          CONTINUE

C----         Jump out of loop
              GO TO 999
          ENDIF
 800  CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PDBNAM  -  Locate the PDB by trying each of the directory
C                        templates in turn
C
C----------------------------------------------------------------------+---

      SUBROUTINE PDBNAM(PDBCOD,TEMPLT,FILNAM,MCSG,UPLOAD)

      INCLUDE 'fcommon.inc'

      CHARACTER*4   PDBCOD
      CHARACTER*(FNAMLN) FILNAM
      CHARACTER*(FNAMLN) TEMPLT(NPDB_DIR)

      INTEGER       IDIR, IPDB, ISTART, LENSTR

      LOGICAL       FOUND, MCSG, UPLOAD

C---- Initialise
      FILNAM = ' '
      IDIR = NOT_FOUND
      ISTART = 1
      IF (UPLOAD) THEN
          ISTART = PDB_UPLOAD
      ELSE IF (MCSG) THEN
          ISTART = PDB_MCSG
      ENDIF

C---- Loop over the PDB file templates
      DO 800, IPDB = ISTART, NPDB_DIR

C----     Form the file name for this template
          CALL FRMFNM(PDBCOD,TEMPLT(IPDB),FILNAM,FNAMLN)

C----     See if the file exists
          INQUIRE(FILE=FILNAM,EXIST=FOUND)

C----     If it exists, then can stop here
          IF (FOUND) THEN

C----         Set flag
              IDIR = IPDB
              GO TO 900
          ENDIF
 800  CONTINUE

C---- Rule out chance hits in the wrong directory
 900  CONTINUE
      IF (IDIR.EQ.PDB_UPLOAD) THEN
          IF (.NOT.UPLOAD) IDIR = NOT_FOUND
      ELSE IF (IDIR.EQ.PDB_MCSG) THEN
          IF (.NOT.MCSG) IDIR = NOT_FOUND
      ELSE
          IF (UPLOAD .OR. MCSG) IDIR = NOT_FOUND
      ENDIF

C---- If file not fouhnd, blank it out
      IF (IDIR.EQ.NOT_FOUND) FILNAM = ' '

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE SUMNAM  -  Identify the PDBsum directory
C
C----------------------------------------------------------------------+---

      SUBROUTINE SUMNAM(PDBCOD,TEMPLT,FILNAM,MCSG,UPLOAD)

      INCLUDE 'fcommon.inc'

      CHARACTER*4   PDBCOD
      CHARACTER*(FNAMLN) FILNAM
      CHARACTER*(FNAMLN) TEMPLT(NPDBSUM_DIR)

      INTEGER       IDIR, LENSTR

      LOGICAL       FOUND, MCSG, UPLOAD

C---- Initialise
      FILNAM = ' '
      FOUND = .FALSE.

C---- Determine the directory type
      IF (.NOT.UPLOAD .AND. .NOT.MCSG) THEN
          IDIR = PDBSUM_CURRENT
      ELSE
          IDIR = PDBSUM_UPLOAD
      ENDIF

C---- Form the name of the PDBsum directory
      CALL FRMFNM(PDBCOD,TEMPLT(IDIR),FILNAM,FNAMLN)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE FRMFNM  -  Form the filename using the given PDB code and
C                        filename template
C
C----------------------------------------------------------------------+---

      SUBROUTINE FRMFNM(PDBCOD,TEMPLT,FNAME,FILEN)

      INCLUDE 'fcommon.inc'

      INTEGER       FILEN

      CHARACTER*1   CH, COPYCH
      CHARACTER*4   PDBCOD
      CHARACTER*(FNAMLN) TEMPLT, FNAME

      INTEGER       I, ICHAR, IPOS, ITO, LEN, LENCOD, LENSTR

      LOGICAL       INCODE

C---- Check that string lengths passed to this subroutine aren't too long
      IF (FILEN.GT.FNAMLN) THEN
        PRINT*, '*** PROGRAM ERROR. String too long in FRMFNM!'
        PRINT*, '***    Submitted: ', FILEN, '    Allowed: ', FNAMLN
        STOP
      ENDIF

C---- Initialise variables
      INCODE = .FALSE.
      IPOS = 0

C---- Get length of template and PDB code
      LEN = LENSTR(TEMPLT)
      LENCOD = LENSTR(PDBCOD)

C---- Initialise the file name
      FNAME = " "

C---- Find any strings between square brackets in this template
      DO 800, I = 1, LEN

C----     Get this character
          CH = TEMPLT(I:I);
          COPYCH = " ";

C----     If this an open-bracket, set flag
          IF (CH.EQ.'[') THEN
            INCODE = .TRUE.

C----     If close-bracket, unset flag
          ELSE IF (CH.EQ.']') THEN
            INCODE = .FALSE.

          ELSE

C----         If in code, transfer appropriate part of PDB code to file name
              IF (INCODE) THEN

C----             See if this character is a template character
                  ICHAR = INDEX(TMPLCH,CH)

C----             If present, transfer corresponding character from PDB code
                  IF (ICHAR.GT.0 .AND. ICHAR.LE.LENCOD) THEN
                      COPYCH = PDBCOD(ICHAR:ICHAR)
                  ENDIF

C----         Otherwise, just transfer character from template
              ELSE
                COPYCH = CH
              ENDIF

C----         If have a character to copy across, perform the transfer
              IF (COPYCH.NE." ") THEN

C----             Transfer the copy character to the file name
                  IPOS = IPOS + 1
                  FNAME(IPOS:IPOS) = COPYCH
              ENDIF
          ENDIF
 800  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE FRMNUM  -  Format the given number
C
C----------------------------------------------------------------------+---

      SUBROUTINE FRMNUM(NUMBER,STRING,LEN)

      CHARACTER*(*) STRING
      INTEGER       I, ISTART, LEN, NUMBER


C---- Convert number to a string
      WRITE(STRING,10) NUMBER
 10   FORMAT(I14)

C---- Count the leading blanks to work out the length of the string
      ISTART = 1
      LEN = 14
      DO 100, I = 1, 14
          IF (STRING(I:I).EQ.' ') THEN
              LEN = LEN - 1
              ISTART = ISTART + 1
          ENDIF
 100  CONTINUE

C---- Reposition the number
      STRING = STRING(ISTART:)

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE LOWCAS  -  Convert text string of length NCHARS to lower-case
C
C----------------------------------------------------------------------+---

      SUBROUTINE LOWCAS(NCHARS,STRING)

      INTEGER       NCHARS

      CHARACTER*1   CH
      CHARACTER*(*) STRING
      INTEGER       IPOS

C---- Loop over all the characters
      DO 600, IPOS = 1, NCHARS

C----     Retrieve this character
          CH = STRING(IPOS:IPOS)

C----     If this is a letter convert into lower-case
          IF (LGE(CH,'A') .AND. LLE(CH,'Z')) THEN
              STRING(IPOS:IPOS) = CHAR(ICHAR('a')
     -            + ICHAR(CH) - ICHAR('A'))
          ENDIF
 600  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE UPPCAS  -  Convert text string of length NCHARS to upper-case
C
C----------------------------------------------------------------------+---

      SUBROUTINE UPPCAS(NCHARS,STRING)

      INTEGER       NCHARS

      CHARACTER*1   CH
      CHARACTER*(*) STRING
      INTEGER       IPOS

C---- Loop over all the characters
      DO 600, IPOS = 1, NCHARS

C----     Retrieve this character
          CH = STRING(IPOS:IPOS)

C----     If this is a letter convert into lower-case
          IF (LGE(CH,'a') .AND. LLE(CH,'z')) THEN
              STRING(IPOS:IPOS) = CHAR(ICHAR('A')
     -            + ICHAR(CH) - ICHAR('a'))
          ENDIF
 600  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  FUNCTION LENSTR  -  Return length of string 'CHARS' excluding
C                      trailing blanks
C
C----------------------------------------------------------------------+---

      INTEGER FUNCTION LENSTR(CHARS)

      CHARACTER*(*) CHARS
      INTEGER J, MVAL

      MVAL = LEN(CHARS)
      DO 10, J = MVAL, 1, -1
          IF (CHARS(J:J).NE.' ') GO TO 20
   10 CONTINUE
      J = 0
   20 CONTINUE
      LENSTR = J

      RETURN
      END

C--------------------------------------------------------------------------
