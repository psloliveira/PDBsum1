      PROGRAM disulphide
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: p_disulph_ensem.f                                        *
C*  FUNCTION: Program to identify and classify disulphide bridges in   * 
C*            proteins                                                 *
C* COPYRIGHT: (1995,1996) University College London                    *
C*---------------------------------------------------------------------*
C*    AUTHOR: GAIL HUTCHINSON                                          *
C*      DATE: 1.2.95                                                   *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
C* 14/02/95  GH     modified to allow running for more than 1 protein  *
C* 03/05/95  GH     to either write data all to one file or to write   *
C*                  separate files for each protein                    *
C* 27/03/96  GH     reorganize multiple/single file handling           *
C* 29/04/96  GH     this version created to allow command line options *
c  01/11/00  GH     ierror changed from integer to logical
c  21/11/00  GH     modified to read sstruc3 format
C***********************************************************************
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -o p_disulph3 p_disulph3.f
C
c***********************************************************************

      INCLUDE 'p_disulph.par'

      INTEGER I,NDSF,SEQNO(MAXDSF),PARTNR(MAXDSF),IOCODE,ISTART,IEND,
     ~     OPTION,NPS,NAMLEN,FINISH,FISTAT,ENDNFI,OPNERB,OK,
     ~     SEQLEN
CPDBSUM RAL 30/3/05 -->
      INTEGER   LENSTR, OPNERG
CPDBSUM RAL 30/3/05 <--
      REAL CHI(3,MAXDSF),CADIST(MAXDSF)
      CHARACTER*(FNAMLN) PDBNAM,SSTNAM,PDBCOD,FINAM,DSFNAM
      CHARACTER DSFDEF(MAXDSF)*1,AAID(MAXDSF)*6,coption
CPDBSUM RAL 30/3/05 -->
      CHARACTER*512 OUTNAM
CPDBSUM RAL 30/3/05 <--
      LOGICAL ENSEMB,FINUSE,IERROR

CPDBSUM RAL 30/3/05 -->
C      PARAMETER(FINISH=2,ENDNFI=6,OPNERB=3,OK=0)
      PARAMETER(FINISH=2,ENDNFI=6,OPNERB=3,OPNERG=4,OK=0)
CPDBSUM RAL 30/3/05 <--
      DATA NDSF/0/,FINUSE/.FALSE./,OPTION/1/,fistat/0/

      READ(*,'(A)') PDBNAM
      IF(PDBNAM(1:1).EQ.'%') THEN
         READ(*,'(A1)') COPTION
         ENSEMB=.TRUE.
         IF(COPTION.EQ.'L'.OR.COPTION.EQ.'l') then
            DSFNAM='DISULPHIDES'
            OPEN(UNIT=OUTFIL,FILE=DSFNAM,STATUS='UNKNOWN',
     ~           IOSTAT=IOCODE)
            IF(IOCODE.NE.0) FISTAT=OPNERB
         ENDIF
      ELSE
         ENSEMB=.FALSE.
         coption='p'
      ENDIF

C     if ensemble remove % from file name
      IF(ENSEMB) THEN
         FINAM=PDBNAM
         FINUSE=.TRUE.
         NAMLEN=INDEX(PDBNAM,SPACE)-1
         DO 180 I=1,NAMLEN-1
            FINAM(I:I)=FINAM(I+1:I+1)
 180     CONTINUE
         FINAM(NAMLEN:NAMLEN)=SPACE
         OPEN(NAMFI,FILE=FINAM,STATUS='OLD',IOSTAT=IOCODE)
c      ENDIF
      IF(NAMLEN.LE.0) THEN
         FISTAT=FINISH
      ENDIF
      endif
 200  CONTINUE
      IF(FISTAT.EQ.OK) THEN
         IF(ENSEMB) THEN
            READ(NAMFI,'(A)',IOSTAT=IOCODE) PDBNAM
            NAMLEN=INDEX(PDBNAM,SPACE)-1
            IF(IOCODE.LT.0) THEN
               CLOSE(NAMFI)
               FINUSE=.FALSE.
               FISTAT=ENDNFI
            ELSE
               WRITE(*,*) PDBNAM(1:NAMLEN)
            ENDIF
c         ELSE
c            PDBNAM=FINAM
         ENDIF
      ENDIF

      IF(FISTAT.EQ.OK) THEN
         CALL GETNAM(PDBNAM,ISTART,IEND,IERROR)
         PDBCOD=PDBNAM(ISTART:IEND)
         NPS=INDEX(PDBCOD,' ')
         SSTNAM=PDBCOD(1:NPS-1)//'.sst'
         OPEN(UNIT=STRKFI,FILE=SSTNAM,STATUS='OLD',IOSTAT=IOCODE)
         IF(IOCODE.NE.0) FISTAT=OPNERB
      ENDIF
      
c     read in secondary structure file
      IF(FISTAT.EQ.OK) THEN
         NDSF=0
         CALL RDSTRK(NDSF,SEQNO,AAID,CHI,PARTNR,DSFDEF,CADIST,SEQLEN)
         CLOSE(UNIT=STRKFI)
         
c     write out results
         IF(NDSF.GT.0) THEN
c     IF(OPTION.ne.2) THEN
            IF(COPTION.NE.'L'.AND.COPTION.NE.'l') then
               DSFNAM=PDBCOD(1:NPS-1)//'.dsf'
CPDBSUM RAL 14/4/05 -->
C               OPEN(UNIT=OUTFIL,FILE=DSFNAM,STATUS='UNKNOWN',
C     ~              IOSTAT=IOCODE)
               OPEN(UNIT=OUTFIL,STATUS='SCRATCH',IOSTAT=IOCODE)
CPDBSUM RAL 14/4/05 <--
               IF(IOCODE.NE.0) FISTAT=OPNERB
c               WRITE(OUTFIL,5) PDBCOD(NPS-4:NPS-1)
c 5             FORMAT('Disulphide Bridges - ',A4)
               write(outfil,5) pdbcod(1:nps-1)
 5             format('Disulphide Bridges - ',A)
CPDBSUM RAL 30/3/05 -->
C----          Open disulph.promotif data file
               outnam = "disulph.promotif"
               open(unit=outdis, file=outnam, status='unknown',
     -              iostat=iocode)
               if (iocode.ne.0)  fistat = opnerg
               call dished(outdis)
CPDBSUM RAL 30/3/05 <--
            ENDIF
            CALL RESULT(NDSF,PARTNR,SEQNO,CHI,AAID,CADIST,SEQLEN,
     ~           COPTION,PDBNAM)
C            IF(OPTION.ne.2) CLOSE(UNIT=OUTFIL)
            IF(COPTION.NE.'L'.AND.COPTION.NE.'l') CLOSE(UNIT=OUTFIL)
         ENDIF
      ELSE IF(FISTAT.EQ.OPNERB) THEN
         write(*,*)'Unable to open the 2D structure file ',
     +              sstnam(1:LENSTR(sstnam))
      ENDIF


      IF(ENSEMB) THEN
         IF(FISTAT.NE.FINISH) THEN
            FISTAT=OK
            IF(FINUSE) THEN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF

      END

c*********************************************************
c
C     SUBROUTINE RDSTRK - subroutine to read in .sst file
c
c**********************************************************
      SUBROUTINE RDSTRK(NDSF,SEQNO,AAID,CHI,PARTNR,DSFDEF,CADIST,
     ~     SEQLEN)

      INCLUDE 'p_disulph.par'

      INTEGER SEQLEN,I,NDSF,SEQNO(MAXDSF),PARTNR(MAXDSF),NRES
      REAL CHI(3,MAXDSF),CADIST(MAXDSF)
      CHARACTER LINE*133,AAID(MAXDSF)*6,DSFDEF(MAXDSF)*1

      DO 10 I=1,3
         READ(STRKFI,*)
 10   CONTINUE
      READ(STRKFI,15) SEQLEN
 15   FORMAT(19X,I4)
      DO 20 I=1,SEQLEN+11
         READ(STRKFI,*)
 20   CONTINUE

C     now at place where disulphides might be
      DO 40 NRES=1,SEQLEN
         READ(STRKFI,'(A)') LINE
         IF((LINE(14:16).EQ.'CYS'.OR.LINE(14:16).EQ.'MPR').
     ~        AND.LINE(63:63).NE.SPACE) THEN
            NDSF=NDSF+1
            READ(LINE,25) SEQNO(NDSF),AAID(NDSF),(CHI(I,NDSF),I=1,3),
     ~           PARTNR(NDSF),DSFDEF(NDSF),CADIST(NDSF)
 25         FORMAT(I5,1X,A6,16X,3(1X,F6.1),7X,I5,1X,A1,2X,F3.1)
         ENDIF
 40   CONTINUE
      
      RETURN
      END
C*****************************************************************
C
C     SUBROUTINE RESULT - to write out disulphide data
c
c******************************************************************
      SUBROUTINE RESULT(NDSF,PARTNR,SEQNO,CHI,AAID,CADIST,SEQLEN,
     ~     COPTION,PDBNAM)

      INCLUDE 'p_disulph.par'

      INTEGER NDSF,SEQNO(MAXDSF),PARTNR(MAXDSF),I,OPTION,J,K,SEQLEN
      REAL CADIST(MAXDSF),CHI(3,MAXDSF),CHI1P(MAXDSF),CHI2P(MAXDSF)
      CHARACTER AAID(MAXDSF)*6,AAIDP(MAXDSF)*6,DTYPE*3,PDBCOD*4,
     ~     COPTION
      CHARACTER*(FNAMLN) PDBNAM

      I=INDEX(PDBNAM,'.')
      IF(I.GT.4) then
         PDBCOD=PDBNAM(I-4:I-1)
      ELSE
         if(i.ne.0) then
            PDBCOD=PDBNAM(1:4)
         else
            k=index(pdbnam,' ')
            pdbcod=pdbnam(k-4:k-1)
         endif
      ENDIF

      DO 100 I=1,NDSF
C     write out each disulphide only once
         IF(PARTNR(I).GT.SEQNO(I)) THEN
c     search for partner residue
            DO 20 J=I+1,NDSF
               IF(SEQNO(J).EQ.PARTNR(I)) THEN
                  CHI2P(I)=CHI(2,J)
                  CHI1P(I)=CHI(1,J)
                  AAIDP(I)=AAID(J)
               ENDIF
 20         CONTINUE
c     test for different classes of disulphides
            DTYPE='   '
            IF(CHI(2,I).LT.0.AND.CHI(3,I).LT.0.AND.CHI2P(I).LT.0) THEN
c     left handed spiral
               DTYPE='LHS'
            ELSE IF(CHI(2,I).GT.0.AND.CHI(3,I).GT.0.AND.CHI2P(I).GT.0)
     ~              THEN
c     right handed spiral
               DTYPE='RHS'
            ELSE IF((CHI(2,I).GT.0.AND.CHI(3,I).GT.0.AND.CHI2P(I).LT.0).
     ~              OR.(CHI(2,I).LT.0.AND.CHI(3,I).GT.0.AND.
     ~              CHI2P(I).GT.0)) THEN
c     right handed hook
               DTYPE='RHH'
            ELSE IF(CHI(2,I).LT.0.AND.CHI(3,I).GT.0.AND.CHI2P(I).LT.0) 
     ~              THEN
c     short right handed hook
               DTYPE='SRH'
            ENDIF
C            IF(OPTION.EQ.2) THEN
            IF(COPTION.eq.'L'.or.COPTION.eq.'l') then
               WRITE(OUTFIL,16) PDBCOD,AAID(I),AAIDP(I),
     ~              (CHI(K,I),K=1,3),CHI2P(I),CHI1P(I),CADIST(I),DTYPE
 16            FORMAT(A4,1X,A6,1X,A6,5(1X,F6.1),1X,F3.1,1X,A3)
            ELSE
               WRITE(OUTFIL,15) AAID(I),AAIDP(I),(CHI(K,I),K=1,3),
     ~              CHI2P(I),CHI1P(I),CADIST(I),DTYPE,SEQNO(I),PARTNR(I)
 15            FORMAT(A6,1X,A6,5(1X,F6.1),1X,F3.1,1X,A3,1X,I4,1X,I4)
CPDBSUM RAL 30/3/05 -->
C----          Write out data to disulph.promotif data file
               WRITE(OUTDIS,640) AAID(I)(1:1), AAID(I)(1:1),
     -             AAID(I)(2:6), AAIDP(I)(1:1), AAIDP(I)(2:6),
     -             (CHI(K,I),K=1,3), CHI2P(I), CHI1P(I), CADIST(I),
     -             DTYPE, SEQNO(I), PARTNR(I)
 640           FORMAT(A1,'::',2(A1,'::',A5,'::'),5(F6.1,'::'),F3.1,'::',
     -             A3,'::',I4,'::',I4,'::')

C----          If disulphide bond spans two chains, write out the same
C              info for the other chains
               IF (AAID(I)(1:1).NE.AAIDP(I)(1:1)) THEN
                   WRITE(OUTDIS,640) AAIDP(I)(1:1), AAID(I)(1:1),
     -                 AAID(I)(2:6), AAIDP(I)(1:1), AAIDP(I)(2:6),
     -                 (CHI(K,I),K=1,3), CHI2P(I), CHI1P(I), CADIST(I),
     -                 DTYPE, SEQNO(I), PARTNR(I)
               ENDIF
CPDBSUM RAL 30/3/05 <--
            ENDIF
         ENDIF
 100  CONTINUE

      RETURN
      END
C**********************************************************************
      SUBROUTINE GETNAM(PDBFIL,ISTART,IEND,IERROR)

      CHARACTER*1   BSLASH, CLOSEB, OPENB, PCHAR, SLASH
      CHARACTER*78  PDBFIL
      INTEGER       IEND, IPOS, ISTART, ISTATE
      LOGICAL       FINISH, GOTDOT, IERROR

      OPENB = '['
      CLOSEB = ']'
      BSLASH = '\\'
      SLASH = '/'

C---- Initialise variables
      FINISH = .FALSE.
      IEND = 0
      IERROR = .FALSE.
      ISTART = 1
      ISTATE = 1
      IPOS = 78
      GOTDOT = .FALSE.

C---- Check through the filename from right to left
100   CONTINUE

C----     Pick off next character
          PCHAR = PDBFIL(IPOS:IPOS)

C----     State 1: Searching for first non-blank character
          IF (ISTATE.EQ.1) THEN
              IF (PCHAR.EQ.SLASH .OR. PCHAR.EQ.BSLASH .OR.
     -            PCHAR.EQ.CLOSEB) THEN
                  GO TO 900
              ENDIF
              IF (PCHAR.NE.' ' .AND. PCHAR.NE.'.') THEN
                  IEND = IPOS
                  ISTATE = 2
              ENDIF

C----     State 2: Searching for end of extension, or end of directory path
          ELSE IF (ISTATE.EQ.2) THEN

C----         If character is a dot, and is the first dot, then note position
              IF (PCHAR.EQ.'.' .AND. .NOT.GOTDOT) THEN
                  IEND = IPOS - 1
                  GOTDOT = .TRUE.

C----         If character signifies the end of a directory path, note pstn
              ELSE IF (PCHAR.EQ.SLASH .OR. PCHAR.EQ.BSLASH
     -            .OR. PCHAR.EQ.CLOSEB) THEN
                  ISTART = IPOS + 1
                  FINISH = .TRUE.
              ENDIF
          ENDIF

C----     Step back a character
          IPOS = IPOS - 1

C---- Loop back for next character
      IF (.NOT.FINISH .AND. IPOS.GT.0) GO TO 100

C---- Check whether file name is sensible
      IF (ISTART.GT.IEND) GO TO 900

      GO TO 999

C---- Error in file name
900   CONTINUE
      IEND = 40
      IF (PDBFIL(41:78).NE.' ') IEND = 78
      PRINT*,' *** ERROR in supplied name of file: [', PDBFIL(1:IEND),
     -    ']'
      IERROR = .TRUE.

999   CONTINUE
      RETURN
      END
c*************************************************************************
CPDBSUM RAL 30/3/05 -->
C**************************************************************************
C
C  SUBROUTINE DISHED  -  Write header records to disulph.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE DISHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 14)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN',
     -     'CHAIN1', 'RES_NUM1',
     -     'CHAIN2', 'RES_NUM2',
     -     'CHI1', 'CHI2', 'CHI3', 'CHIP2', 'CHIP1',
     -     'CA_DIST', 'TYPE', 'SEQ_NO', 'PARTNER' /

C---- Define field id tag
      IDTAG = 'DS_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:DISULPHIDES'
 10   FORMAT(500A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENSTR(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENSTR(FLDNAM(NFIELD)))

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE AA2NAM  -  Convert single-letter a.a. code to 3-letter residue
C                        name
C
C----------------------------------------------------------------------+---

      SUBROUTINE AA2NAM(AACODE,RESNAM)

      INTEGER       NAMINO

      PARAMETER    (NAMINO = 20)

      CHARACTER*1   AACODE
      CHARACTER*3   CODE(NAMINO), RESNAM
      CHARACTER*(NAMINO)  AMINOS

      INTEGER       IPOS

      DATA  AMINOS / 'ACDEFGHIKLMNPQRSTVWY'/
      DATA  CODE   / 'ALA','CYS','ASP','GLU','PHE','GLY','HIS',
     -               'ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG',
     -               'SER','THR','VAL','TRP','TYR' /

C---- Initialise residue name
      RESNAM = 'UNK'

C---- Get identifier for this single-letter code
      IPOS = INDEX(AMINOS,AACODE)

C---- If valid code, then store corresponding name
      IF (IPOS.GT.0 .AND. IPOS.LE.NAMINO) THEN
          RESNAM = CODE(IPOS)
      ENDIF

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
CPDBSUM RAL 30/3/05 <--
