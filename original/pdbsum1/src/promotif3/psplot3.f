      PROGRAM psplot
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: psplot.f                                                  *
C*  FUNCTION: Program to generate postscript files for motifs in promotif*
C* COPYRIGHT: University College London                                *
C*---------------------------------------------------------------------*
C*    AUTHOR: GAIL HUTCHINSON                                          *
C*      DATE: 1994                                                     *
C*    DEVELOPED AT UNIVERSITY COLLEGE IN THE BSM GROUP LED BY          *
C*                 JANET THORNTON                                      *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: postscript routines in psrout.f                          *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
C* 11/01/95  GH     change hairpin plotting routines                   *
C* 07/02/95  GH     add in summary page and disulphide page            *
C*                  add in subroutine GETNAM to cope with general      *
C*                  names for input files                              *
C* 11/11/95  GH     add sequence to summary page. Put helix info in    *
c*                  one table                                          *
C* 11/03/96  GH     put helix info back into 2 tables - not enough room*
C* 18/03/96  GH     added new routine NEWPAG to create new postscript  *
c*                  page in SUMMARY routine                            *
c* 18/03/96  GH     output file names changed to allow more than 9     *
C*                  pages of output per motif                          *
C* 20/03/96  GH     allowed for more than one output file of bab and psi
C* 25/03/96  GH     small correction to summary page                   *
c* 03/04/96  GH     small correction to SCHBUL 
c* 17/04/96  GH     new program to run multiple proteins without large *
c*                  arrays
c* 24/04/96  GH     small correction to hairpin drawing program        *
c* 24/04/96  GH     modify psi-loop drawing routine                    *
c* 24/04/96  GH     small correction to summary plot (helix data)      *
c* 21/05/96  GH     allow loop up to 125 residues in hairpin           *
c* 22/05/96  GH     correct drawing of table for helix interactions    *
c* 14/10/96  GH     correct drawing of bulges 
c* 14/10/96  GH     plot bulges if Phi,psi not = -999.9 (rather than 999.9)
C* 14/10/96  GH     improve plotting of loop residues in hairpins
c* 18/12/96  GH     allow hairpins to have more than 100 residues in loop
c  02/11/00  GH     some alterations to avoid exceeding array sizes (sseq1)
C  27/11/00  GH     modify to read in alternative format for sheet data
c                   (numbers to denote sheets rather than letters)
c  04/12/00  GH     further changes to correct reading in hairpin data
c  04/12/00  GH     remove '90 call hinttb' which seems unnecessary?
C***********************************************************************
C
C f77 -C -c psplot3.f
C f77 -C -c psrout.f
C f77 -C -c fcommon.f
C f77 -o psplot3 psplot3.o psrout.o fcommon.o
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -c psplot3.f
C f77 -Wimplicit -fbounds-check -c psrout.f
C f77 -Wimplicit -fbounds-check -c fcommon.f
C f77 -o psplot3 psplot3.o psrout.o fcommon.o
C
C***********************************************************************

      INCLUDE 'psplot.par'

CPDBSUM RAL 15/3/05 -->
C      CHARACTER*78 FINAM,INNAM,PDBNAM,PDBCOD
      CHARACTER*1   CHAIN, DIRSEP
CPDBSUM1 RAL 4/3/22 -->
C      CHARACTER*4   PDBCOD, POPT
      CHARACTER*4   POPT
      CHARACTER*78  PDBCOD
CPDBSUM1 RAL 4/3/22 <--
      CHARACTER*20  CHNAME, PLTYPE, SPCIAL
      CHARACTER*(FNAMLN) PDBSUM_TEMPLATE(NPDBSUM_DIR)

      INTEGER       MSTART, IUNIT

      LOGICAL       IFAIL, ONEPAG
CPDBSUM RAL 15/3/05 <--
CPDBSUM1 RAL 4/3/22 -->
      CHARACTER*(FNAMLN) ADDRES
CPDBSUM1 RAL 4/3/22 <--
      INTEGER ENDNFI,FINISH,FISTAT,ISTART,IEND,IERROR,IOCODE,NAMFI,
     ~     NAMLEN,OK,I
      CHARACTER SPACE
      LOGICAL ENSEMB,FINUSE

      PARAMETER (ENDNFI=6,FINISH=2,NAMFI=10,OK=0)
      DATA FINUSE/.false./,SPACE/' '/,FISTAT/0/

CPDBSUM1 RAL 4/3/22 -->
C---- Initialise PDBsum1 flag
      PDBSM1 = .FALSE.
CPDBSUM1 RAL 4/3/22 <--

c     read in colour parameters
      CALL PARAMS

c     read in name of input file
c      READ(*,'(A)') PDBNAM
c      CALL GETNAM(PDBNAM,ISTART,IEND,IERROR)
c      PDBCOD=PDBNAM(ISTART:IEND)
      
CPDBSUM RAL 15/3/05 -->
C      READ(*,'(A)') INNAM
C      IF(INNAM(1:1).EQ.'%') THEN
C         ENSEMB=.TRUE.
C      ELSE
C         ENSEMB=.FALSE.
C      ENDIF
C      IF(ENSEMB) THEN
C         FINAM=INNAM
C         FINUSE=.TRUE.
C         NAMLEN=INDEX(INNAM,SPACE)-1
C         DO 180 I=1,NAMLEN-1
C            FINAM(I:I)=FINAM(I+1:I+1)
C 180     CONTINUE
C         FINAM(NAMLEN:NAMLEN)=SPACE
C         OPEN(NAMFI,FILE=FINAM,STATUS='OLD',IOSTAT=IOCODE)
C      ENDIF
C 200  CONTINUE
C      IF(FISTAT.EQ.OK) THEN
C         IF(ENSEMB) THEN
C            READ(NAMFI,'(A)',IOSTAT=IOCODE) PDBNAM
C            namlen=index(pdbnam,space)-1
C            IF(IOCODE.LT.0) THEN
C               CLOSE(NAMFI)
C               FINUSE=.FALSE.
C               FISTAT=ENDNFI
C            ELSE
C               WRITE(*,*) PDBNAM(1:namlen)
C            ENDIF
C         ELSE
C            PDBNAM=INNAM
C         ENDIF
C      ENDIF

C---- Initialise variables
      ENSEMB = .FALSE.
      IFAIL = .FALSE.

C---- Read in PDB code and chain identifier and which plot is required
      READ(*,'(A)') PDBCOD
      READ(*,'(A)') CHNAME
      IF (CHNAME.EQ.'NONE') THEN
          CHAIN = ' '
      ELSE
          CHAIN = CHNAME(1:1)
      ENDIF
      READ(*,'(A)') PLTYPE
      READ(*,'(A)') SPCIAL
      READ(*,'(I10)') MSTART
      READ(*,'(A)') POPT
      IF (POPT.EQ.'ONE') THEN
          ONEPAG = .TRUE.
      ELSE
          ONEPAG = .FALSE.
      ENDIF
CPDBSUM1 RAL 4/3/22 -->
      READ(*,'(A)') ADDRES
      if (ADDRES.EQ.'PDBsum1') THEN
         PDBSM1 = .TRUE.
      ENDIF
CPDBSUM1 RAL 4/3/22 <--

C---- Get the paths to the files that may be required
C      CALL GETPTH(PRODIR,UPLDIR,DIRSEP,FILEN)
      DIRSEP = '/'

C---- Get the PDBsum directory templates
      CALL GETPSM(PDBSUM_TEMPLATE)

C---- Open the promotif.dat file and read through to the required section
      CALL OPNDAT(PDBCOD,CHAIN,PLTYPE,SPCIAL,PDBSUM_TEMPLATE,DIRSEP,
CPDBSUM1 RAL 4/3/22 -->
C     -    FILDAT,MSTART,IFAIL)
     -    FILDAT,MSTART,PDBSM1,IFAIL)
CPDBSUM1 RAL 4/3/22 <--
      IF (IFAIL) GO TO 999
CPDBSUM RAL 15/3/05 <--

      IF(FISTAT.EQ.OK) THEN
CPDBSUM RAL 15/3/05 -->
C         CALL GETNAM(PDBNAM,ISTART,IEND,IERROR)
C         PDBCOD=PDBNAM(ISTART:IEND)
C         IF(BETA) CALL BTURN(PDBCOD)
C         IF(GAMMA) CALL GAMTRN(PDBCOD)
C         IF(BULG) CALL BULGE(PDBCOD)
C         IF(HPIN) CALL HAIRPN(PDBCOD)
C         IF(SHT) CALL SHEET(PDBCOD)
C         IF(HLX) CALL HELIX(PDBCOD)
C         IF(DSF) CALL DISULPH(PDBCOD)
C         IF(SUM) CALL SUMMRY(PDBCOD)
         IF (PLTYPE.EQ.'BETA_TURNS') THEN
             CALL BTURN(PDBCOD,CHAIN,ONEPAG,MSTART)
         ELSE IF (PLTYPE.EQ.'GAMMA_TURNS') THEN
             CALL GAMTRN(PDBCOD,CHAIN,ONEPAG,MSTART)
         ELSE IF (PLTYPE.EQ.'BETA_BULGES') THEN
             CALL BULGE(PDBCOD,CHAIN,ONEPAG,MSTART)
         ELSE IF (PLTYPE.EQ.'BETA_HAIRPINS') THEN
             CALL HAIRPN(PDBCOD,CHAIN,ONEPAG,MSTART)
         ELSE IF(PLTYPE.EQ.'HELICES') THEN
             CALL HELIX(PDBCOD,CHAIN,ONEPAG,MSTART)
C         ELSE IF(PLTYPE.EQ.'SHEETS') THEN
C             CALL SHEET(PDBCOD,CHAIN)
C         ELSE IF(DSF) THEN
C             CALL DISULPH(PDBCOD,CHAIN)
C         ELSE IF(SUM) THEN
C             CALL SUMMRY(PDBCOD,CHAIN)
         ENDIF
CPDBSUM RAL 15/3/05 <--
      ENDIF

CPDBSUM RAL 15/3/05 -->
C      IF(ENSEMB) THEN
C         IF(FISTAT.NE.FINISH) THEN
C            FISTAT=OK
C            IF(FINUSE) GO TO 200
C         ENDIF
C      ENDIF
C
C      CLOSE(NAMFI)
CPDBSUM RAL 15/3/05 <--

C---- Close the file
      CLOSE(FILDAT)

 999  CONTINUE
      STOP
      END
C**********************************************************
C     SUBROUTINE PARAMS - to read in parameter file for promotif
C
C**********************************************************
      SUBROUTINE PARAMS
      
      INCLUDE 'psplot.par'
      
      INTEGER I,ICOL,ITRN,IBUL,IHPIN,IHLX
      CHARACTER*12 TCOLOR(7)*12,LINE*79,YESNO
      LOGICAL LOGICS

      OPEN(UNIT=10,FILE='promotif2.prm')

C     determine which plots to draw
      READ(10,'(A1)') YESNO
      BETA=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      GAMMA=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      BULG=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      HPIN=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      SHT=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      HLX=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      DSF=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      SUM=LOGICS(YESNO)
      READ(10,'(A1)') YESNO
      BW=LOGICS(YESNO)

c     read in colours
      IF(.NOT.BW) THEN
         REWIND(UNIT=10)
 100     READ(10,'(A)') LINE
         IF(LINE(1:8).NE.'Colours:') GO TO 100
         DO 10 ICOL=1,MAXCOL
            READ(10,*,END=15)(RGB(I,ICOL),I=1,3),COLNAM(ICOL)
 10      CONTINUE
 15      CONTINUE
         
C     now read in colour options for turns and match to rgb nums
         REWIND(UNIT=10)
 200     READ(10,'(A)') LINE
         IF(LINE(1:16).NE.'Turn Parameters:') GO TO 200
         DO 30 ITRN=1,7
            READ(10,'(A12)') TCOLOR(ITRN)
            DO 20 ICOL=1,20
               IF(TCOLOR(ITRN).EQ.COLNAM(ICOL)) THEN
                  TRNCOL(ITRN)=ICOL
               ENDIF
 20         CONTINUE
 30      CONTINUE
         
c     now read in colour options for gamma turns and match as above
         REWIND(UNIT=10)
 300     READ(10,'(A)') LINE
         IF(LINE(1:22).NE.'Gamma Turn Parameters:') GO TO 300
         DO 50 ITRN=1,6
            READ(10,'(A12)') TCOLOR(ITRN)
            DO 40 ICOL=1,20
               IF(TCOLOR(ITRN).EQ.COLNAM(ICOL)) THEN
                  GAMCOL(ITRN)=ICOL
               ENDIF
 40         CONTINUE
 50      CONTINUE
         
C     read in colour options for bulges
         REWIND(UNIT=10)
 400     READ(10,'(A)') LINE
         IF(LINE(1:17).NE.'Bulge Parameters:') GO TO 400
         DO 70 IBUL=1,6
            READ(10,'(A12)') TCOLOR(IBUL)
            DO 60 ICOL=1,20
               IF(TCOLOR(IBUL).EQ.COLNAM(ICOL)) THEN
                  BULCOL(IBUL)=ICOL
               ENDIF
 60         CONTINUE
 70      CONTINUE
         
c     read in colour options for hairpins
         REWIND(UNIT=10)
 500     READ(10,'(A)') LINE
         IF(LINE(1:19).NE.'Hairpin Parameters:') GO TO 500
         DO 90 IHPIN=1,6
            READ(10,'(A12)') TCOLOR(IHPIN)
            DO 80 ICOL=1,20
               IF(TCOLOR(IHPIN).EQ.COLNAM(ICOL)) THEN
                  HPCOL(IHPIN)=ICOL
               ENDIF
 80         CONTINUE
 90      CONTINUE
         
c     read in colour options for helical diagrams
         REWIND(UNIT=10)
 600     READ(10,'(A)') LINE
         IF(LINE(1:19).NE.'Helix Parameters:') GO TO 600
         DO 110 IHLX=1,4
            READ(10,'(A12)') TCOLOR(IHLX)
            DO 120 ICOL=1,20
               IF(TCOLOR(IHLX).EQ.COLNAM(ICOL)) THEN
                  HXCOL(IHLX)=ICOL
               ENDIF
 120        CONTINUE
 110     CONTINUE
      ENDIF
      
      CLOSE(UNIT=10)
      
      RETURN
      END

C**********************************************************
C      SUBROUTINE BTURN - To plot out postscript for beta turns
C
C**********************************************************
CPDBSUM RAL 15/3/05 -->
C      SUBROUTINE BTURN(PDBCOD)
      SUBROUTINE BTURN(PDBCOD,CHAIN,ONEPAG,MSTART)
CPDBSUM RAL 15/3/05 <--

      INCLUDE 'psplot.par'

      INTEGER I,J,NUMFI,NUMTFI,NPS,JPS,NUM
      REAL X0,Y0,PHI(2),PSI(2),CADIST,XDEL,YDEL,RSIZ,XPLOT1,
     ~     YPLOT1,XPLOT2,YPLOT2,FSIZ,XSTART,YSTART,CRAD,SSIZ,
     ~     AMINT,AMAXT,AHDLEN,SMAXT,SHDLEN,Y,CHI1(2)
CPDBSUM RAL 15/3/05 -->
      CHARACTER*1 CHAIN
      CHARACTER*3 C3, RSNAM(4)
      CHARACTER*4 SEQUEN
      CHARACTER*5 HBFLAG, RSNUM(4)
      CHARACTER*(RECLEN) IREC
CPDBSUM RAL 15/3/05 <--
      CHARACTER RES(4)*6,SEQ(4),TURNTY*4,EFIMOV*2,PDBCOD*78,
     ~     CA*5,C,HBPRES*1,TEXT3*3,TEXT6*6,TEXT8*8,TTL*4,C2*2
      CHARACTER*78 BNAME,BNAMTB,ITRNNM
      LOGICAL BT
CPDBSUM RAL 15/3/05 -->
      INTEGER MOTIF, MSTART
      LOGICAL EOF, ONEPAG
CPDBSUM RAL 15/3/05 <--

CPDBSUM RAL 5/4/05 -->
C      DATA XSTART/60.0/,YSTART/660.0/,XDEL/150.0/,YDEL/180.0/,
C     ~     RSIZ/50.0/,FSIZ/12.0/,CRAD/2.0/,SSIZ/4.0/,
      DATA XSTART/60.0/,YSTART/660.0/,XDEL/130.0/,YDEL/180.0/,
     ~     RSIZ/75.0/,FSIZ/12.0/,CRAD/2.0/,SSIZ/4.0/,
CPDBSUM RAL 5/4/05 <--
     ~     AMINT/0.5/,AMAXT/4.0/,AHDLEN/4.0/,
     ~     SMAXT/3.0/,SHDLEN/3.0/,
     ~     CA(4:5)/' A'/,TTL/'BETA'/

      BT=.FALSE.
      NUMFI=1
      NUMTFI=1
CPDBSUM RAL 15/3/05 -->
C      NPS=INDEX(PDBCOD,' ')
C      ITRNNM=PDBCOD(1:NPS-1)//'.bturns'
C      OPEN(UNIT=2,FILE=ITRNNM,STATUS='OLD',ERR=900)
C      READ(2,*,end=890)
C      READ(2,*,end=890)
CPDBSUM RAL 15/3/05 <--

C     plot a ramachandran and schematic diagram for each turn
      X0=XSTART
      Y0=YSTART
CPDBSUM RAL 15/3/05 -->
C      DO 20 I=1,NNN
C         READ(2,15,END=30) RES(1),SEQ(1),RES(2),SEQ(2),RES(3),SEQ(3),
C     ~        RES(4),SEQ(4),TURNTY,EFIMOV,PHI(1),PSI(1),PHI(2),PSI(2),
C     ~        CADIST,HBPRES,CHI1(1),CHI1(2)
C
C 15      FORMAT(4(A6,A1),1X,A4,1X,A2,1X,4F6.1,2X,F3.1,1X,A1,2(F6.1))
C         IF(RES(1).EQ.'      ') GO TO 30
      I = 0
      IF (ONEPAG) THEN
          MOTIF = MSTART
      ELSE
          MOTIF = 0
      ENDIF
 20   CONTINUE

C----     Get the next record from the promotif.dat file
          CALL GETREC(FILDAT,CHAIN,IREC,EOF)
          IF (EOF) GO TO 30

C----     Extract the fields from the line just read in
          READ(IREC,25) (RSNUM(J), RSNAM(J), J = 1, 4), SEQUEN,
     -        TURNTY,EFIMOV,PHI(1),PSI(1),PHI(2),PSI(2),
     -        CADIST,HBFLAG,CHI1(1),CHI1(2)
 25       FORMAT(1X,4(2X,A5,2X,A3),2X,A4,2X,A4,2X,A2,2X,4(F6.1,2X),F4.1,
     -        2X,A5,2X,2(F6.1,2X))
          I = I + 1

C----     Form incomplete fields
          DO 400, J = 1, 4
              RES(J) = CHAIN // RSNUM(J)
              SEQ(J) = SEQUEN(J:J)
 400      CONTINUE
          IF (HBFLAG.EQ.'TRUE ') THEN
              HBPRES = 'Y'
          ELSE
              HBPRES = 'N'
          ENDIF
CPDBSUM RAL 15/3/05 <--
         IF(I.EQ.1) THEN
c            NPS=INDEX(PDBCOD,' ')
CPDBSUM RAL 16/3/05 -->
C            BNAME=PDBCOD(1:NPS-1)//'_bturns_01.ps'
            BNAME = 'promotif001.ps'
CPDBSUM RAL 16/3/05 <--
            CALL PSOPEN(BNAME)
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'BETA',NUMFI,BNAME,1)
CPDBSUM RAL 16/3/05 <--
            BT=.TRUE.
CPDBSUM RAL 1/4/05 -->
C         ELSE IF(MOD(I,16).EQ.1) THEN
C            NUMFI=NUMFI+1
C            JPS=INDEX(BNAME,'.ps')
CC*****CHANGED HERE 18.3.96
C            if(numfi.lt.10) then
C               WRITE(C,'(I1)') NUMFI
C               BNAME(JPS-1:JPS-1)=C
C            else
C               WRITE(C2,'(I2)') NUMFI
C               BNAME(JPS-2:JPS-1)=C2
C            ENDIF
C            CALL PSOPEN(BNAME)
C            CALL TITLE(PDBCOD,'BETA',NUMFI,BNAME,1)
C----    If plotting all motifs (ie one-per-page, as opposed to just
C        a single page of several motifs) start a new page
         ELSE IF (.NOT.ONEPAG) THEN
            NUMFI=NUMFI+1
            JPS=INDEX(BNAME,'.ps')
            WRITE(C3,'(I3.3)') NUMFI
            BNAME(JPS-3:JPS-1) = C3
            CALL PSOPEN(BNAME)
CPDBSUM RAL 1/4/05 <--
            X0=XSTART
            Y0=YSTART
         ENDIF
C     draw ramachandran
CPDBSUM RAL 5/4/05 -->
C         IF(.NOT.BW) CALL COLOUR(TRNCOL(1))
C         CALL RAMACH(X0, Y0, X0+RSIZ, Y0+RSIZ)
         CALL RAMACH(X0, Y0, X0+RSIZ, Y0+RSIZ, TRNCOL(1))
CPDBSUM RAL 5/4/05 <--

c     plot phi,psi of i+1, i+2
CPDBSUM RAL 5/4/05 -->
C         IF(.NOT.BW) CALL COLOUR(TRNCOL(2))
         IF(.NOT.BW) CALL COLOUR(7)
CPDBSUM RAL 5/4/05 <--
         CALL COORD(PHI(1),PSI(1),RSIZ,X0,Y0,XPLOT1,YPLOT1)
         CALL FPCIRC(XPLOT1,YPLOT1,CRAD)
CPDBSUM RAL 5/4/05 -->
         CALL COLOUR(4)
CPDBSUM RAL 5/4/05 <--
         CALL COORD(PHI(2),PSI(2),RSIZ,X0,Y0,XPLOT2,YPLOT2)
         CALL FPSQUAR(XPLOT2,YPLOT2,SSIZ)
         IF(.NOT.BW) CALL COLOUR(TRNCOL(7))
         CALL FARROW(XPLOT1,YPLOT1,XPLOT2,YPLOT2
     ~                    ,AMINT,SMAXT,SHDLEN)
CPDBSUM RAL 5/4/05 -->
C         IF(.NOT.BW) CALL COLOUR(TRNCOL(3))
C         CALL PLTNUM(RES(1),RES(4),TURNTY,X0,Y0,RSIZ)
         CALL PLTNM2(RSNAM(1),RSNUM(1),RSNAM(4),RSNUM(4),TURNTY,X0,Y0,
     -       RSIZ,MOTIF)
CPDBSUM RAL 5/4/05 <--

c     draw schematic turn
         IF(.NOT.BW) CALL COLOUR(TRNCOL(4))
CPDBSUM RAL 5/4/05 -->
C         CALL PBTURN(X0, Y0-50.0, 50.0, 30.0,CRAD,SSIZ)
C         IF(.NOT.BW) CALL COLOUR(TRNCOL(5))
C         CALL TPRINT(SEQ(1),FSIZ,X0-FSIZ/2.0,Y0-50.0-FSIZ/2.0)
C         CALL TPRINT(SEQ(2),FSIZ,X0-FSIZ/2.0,Y0-20.0-FSIZ/2.0)
C         CALL TPRINT(SEQ(3),FSIZ,X0+50.0+FSIZ/2.0,Y0-20.0-FSIZ/2.0)
C         CALL TPRINT(SEQ(4),FSIZ,X0+50.0+FSIZ/2.0,Y0-50.0-FSIZ/2.0)
         CALL PBTURN(X0, Y0-50.0, RSIZ, 30.0,CRAD,SSIZ)
         CALL COLOUR(3)
         CALL BPRINT(SEQ(1),FSIZ,X0-FSIZ/2.0,Y0-50.0-FSIZ/2.0)
         CALL COLOUR(7)
         CALL BPRINT(SEQ(2),FSIZ,X0-FSIZ/2.0,Y0-20.0-FSIZ/2.0)
         CALL COLOUR(4)
         CALL BPRINT(SEQ(3),FSIZ,X0+RSIZ+FSIZ/2.0,Y0-20.0-FSIZ/2.0)
         CALL COLOUR(5)
         CALL BPRINT(SEQ(4),FSIZ,X0+RSIZ+FSIZ/2.0,Y0-50.0-FSIZ/2.0)
CPDBSUM RAL 5/4/05 <--

c     plot ca(i)-ca(i+3) distance
         IF(.NOT.BW) CALL COLOUR(TRNCOL(6))
         WRITE(CA(1:3),'(F3.1)') CADIST
CPDBSUM RAL 5/4/05 -->
C         CALL TPRINT(CA,FSIZ,X0+25.0,Y0-58.0)
         CALL TPRINT(CA,FSIZ,X0+RSIZ/2.0,Y0-58.0)
CPDBSUM RAL 5/4/05 <--
         CALL DASH(3)
CPDBSUM RAL 5/4/05 -->
C         CALL PLINE(X0+2.0,Y0-50.0,X0+48.0,Y0-50.0)
         CALL PLINE(X0+2.0,Y0-50.0,X0+RSIZ-2.0,Y0-50.0)
CPDBSUM RAL 5/4/05 <--
         CALL DASH(0)
         IF(.NOT.BW) CALL COLOUR(TRNCOL(7))
         
CPDBSUM RAL 5/4/05 -->
C         IF(HBPRES.EQ.'Y') CALL ARROW(X0+48.0,Y0-40.0,
C     ~        X0+2.0,Y0-40.0,AMINT,AMAXT,AHDLEN)
         IF(HBPRES.EQ.'Y') THEN
             CALL ARROW(X0+RSIZ-2.0,Y0-40.0,X0+2.0,Y0-40.0,AMINT,AMAXT,
     -           AHDLEN)
         ENDIF
CPDBSUM RAL 5/4/05 <--
C     check if end of page -close file and start again
CPDBSUM RAL 1/4/05 -->
C         IF(MOD(I,16).EQ.0.AND.I.GT.0) THEN
         IF(ONEPAG .AND. MOD(I,16).EQ.0) THEN
             GO TO 30
         ELSE IF (.NOT.ONEPAG) THEN
CPDBSUM RAL 1/4/05 <--
            CALL PSCLOSE
c     check if end of line
         ELSE IF(MOD(I,4).EQ.0) THEN
            Y0=Y0-YDEL
            X0=XSTART
c     otherwise just increment x
         ELSE
            X0=X0+XDEL
         ENDIF
CPDBSUM RAL 15/3/05 -->
C 20   CONTINUE
      IF (ONEPAG) MOTIF = MOTIF + 1
      GO TO 20
CPDBSUM RAL 15/3/05 <--
 30   IF(BT) CALL PSCLOSE

C     write out table of data
CPDBSUM RAL 15/3/05 -->
C      IF(BT) THEN
C         REWIND(UNIT=2)
C         BNAMTB=PDBCOD(1:NPS-1)//'_bturn_tab01.ps'
C         CALL PSOPEN(BNAMTB)
C         CALL TITLE(PDBCOD,'BETA',NUMTFI,BNAMTB,2)
C         READ(2,*,end=890)
C         READ(2,*,end=890)
C         Y=640.0
C         DO 40  I=1,NNN
C            NUM=I
C            READ(2,15,END=50) RES(1),SEQ(1),RES(2),SEQ(2),RES(3),
C     ~           SEQ(3),RES(4),SEQ(4),TURNTY,EFIMOV,PHI(1),PSI(1),
C     ~           PHI(2),PSI(2),CADIST,HBPRES,CHI1(1),CHI1(2)
C            IF(RES(1).EQ.'      ') GO TO 50
C            IF(Y.LT.50.0) THEN
C               CALL BTRNTB(Y)
C               CALL NEWPAG(NUMTFI,BNAMTB,PDBCOD,TTL)
C               Y=640.0
C            ENDIF
C            WRITE(TEXT6,'(A6)') RES(1)
C            CALL TPRINT(TEXT6,12.0,70.0,Y)
C            CALL TPRINT('-',12.0,90.0,Y)
C            WRITE(TEXT6,'(A6)') RES(4)
C            CALL TPRINT(TEXT6,12.0,110.0,Y)
C            WRITE(TEXT8,'(4(A1,1X))') SEQ(1),SEQ(2),SEQ(3),SEQ(4)
C            CALL TPRINT(TEXT8,12.0,155.0,Y)
C            CALL TPRINT(TURNTY,12.0,200.0,Y)
C            IF(PHI(1).NE.999.9) THEN
C               WRITE(TEXT6,'(F6.1)') PHI(1)
C               CALL TPRINT(TEXT6,12.0,240.0,Y)
C            ELSE
C               CALL TPRINT('-',12.0,240.0,Y)
C            ENDIF
C            IF(PSI(1).NE.999.9) THEN
C               WRITE(TEXT6,'(F6.1)') PSI(1)
C               CALL TPRINT(TEXT6,12.0,280.0,Y)
C            ELSE
C               CALL TPRINT('-',12.0,280.0,Y)
C            ENDIF
C            IF(PHI(2).NE.999.9) THEN
C               WRITE(TEXT6,'(F6.1)') PHI(2)
C               CALL TPRINT(TEXT6,12.0,320.0,Y)
C            ELSE
C               CALL TPRINT('-',12.0,320.0,Y)
C            ENDIF
C            IF(PSI(2).NE.999.9) THEN
C               WRITE(TEXT6,'(F6.1)') PSI(2)
C               CALL TPRINT(TEXT6,12.0,360.0,Y)
C            ELSE
C               CALL TPRINT('-',12.0,360.0,Y)
C            ENDIF
C            CALL TPRINT(EFIMOV,12.0,400.0,Y)
C            IF(CHI1(1).NE.999.9) THEN
C               WRITE(TEXT6,'(F6.1)') CHI1(1)
C               CALL TPRINT(TEXT6,12.0,440.0,Y)
C            ELSE
C               CALL TPRINT('-',12.0,440.0,Y)
C            ENDIF
C            IF(CHI1(2).NE.999.9) THEN
C               WRITE(TEXT6,'(F6.1)') CHI1(2)
C               CALL TPRINT(TEXT6,12.0,480.0,Y)
C            ELSE
C               CALL TPRINT('-',12.0,480.0,Y)
C            ENDIF
C            WRITE(TEXT3,'(F3.1)') CADIST
C            CALL TPRINT(TEXT3,12.0,525.0,Y)
C            Y=Y-20.0
C 40      CONTINUE      
C 50      IF(NUM.GT.1) CALL BTRNTB(Y)
C         IF(NUM.GT.1) CALL PSCLOSE
C      ENDIF
C      CLOSE(UNIT=2)
C      GO TO 900
C 890  WRITE(*,*) 'No beta turn data'
CPDBSUM RAL 15/3/05 -->

 900  RETURN
      END
C******************************************************************
C      SUBROUTINE GAMTRN - Plots out data for gamma turns
C
C******************************************************************
CPDBSUM RAL 16/3/05 -->
C      SUBROUTINE GAMTRN(PDBCOD)
      SUBROUTINE GAMTRN(PDBCOD,CHAIN,ONEPAG,MSTART)
CPDBSUM RAL 16/3/05 <--
      
      INCLUDE 'psplot.par'

      INTEGER I,J,NUMFI,NUMTFI,NPS,JPS,NUM
      REAL PHI,PSI,CADIST,XSTART,YSTART,XDEL,YDEL,X0,Y0,
     ~     RSIZ,XPLOT,YPLOT,CRAD,FSIZ,Y
CPDBSUM RAL 15/3/05 -->
      CHARACTER*1 CHAIN
      CHARACTER*3 C3, RSNAM(3), SEQUEN
      CHARACTER*5 RSNUM(3)
      CHARACTER*(RECLEN) IREC
CPDBSUM RAL 15/3/05 <--
      CHARACTER LINE*80
      CHARACTER RES(3)*6,SEQ(3),TURNTY*7,C,CA*5,
     ~     TEXT6*6,TEXT3*3,TTL*4,c2*2
      CHARACTER*78  GNAME,GNAMTB,ITRNNM,PDBCOD
      LOGICAL GT
CPDBSUM RAL 15/3/05 -->
      INTEGER MOTIF, MSTART
      LOGICAL EOF, ONEPAG
CPDBSUM RAL 15/3/05 <--

CPDBSUM RAL 5/4/05 -->
C      DATA XSTART/60.0/,YSTART/660.0/,XDEL/150.0/,YDEL/180.0/,
C     ~     RSIZ/50.0/,CRAD/2.0/,FSIZ/12.0/,
      DATA XSTART/60.0/,YSTART/660.0/,XDEL/130.0/,YDEL/180.0/,
     ~     RSIZ/75.0/,CRAD/2.0/,FSIZ/12.0/,
CPDBSUM RAL 5/4/05 <--
     ~     CA(4:5)/' A'/,TTL/'GAMM'/

      GT=.FALSE.
      NUMFI=1
      NUMTFI=1
CPDBSUM RAL 16/3/05 -->
C      NPS=INDEX(PDBCOD,' ')
C      ITRNNM=PDBCOD(1:NPS-1)//'.gturns'
C      OPEN(UNIT=2,FILE=ITRNNM,STATUS='OLD',ERR=900)
C      READ(2,*)
C      
C 10   READ(2,'(A80)',end=890) LINE
C      IF(LINE(2:6).NE.'Gamma') GO TO 10
CPDBSUM RAL 16/3/05 <--

      X0=XSTART
      Y0=YSTART
CPDBSUM RAL 16/3/05 -->
C      DO 30 I=1,NNN
C         READ(2,15,END=40) RES(1),SEQ(1),RES(2),SEQ(2),RES(3),SEQ(3),
C     ~        TURNTY,PHI,PSI,CADIST
C 15      FORMAT(3(A6,A1),1X,A7,2X,2(F6.1),1X,F3.1)
C         IF(RES(1).EQ.'      ') GO TO 40
       I = 0
       IF (ONEPAG) THEN
           MOTIF = MSTART
       ELSE
           MOTIF = 0
       ENDIF
 30    CONTINUE

C----     Get the next record from the promotif.dat file
          CALL GETREC(FILDAT,CHAIN,IREC,EOF)
          IF (EOF) GO TO 40

C----     Extract the fields from the line just read in
          READ(IREC,25) (RSNUM(J), RSNAM(J), J = 1, 3), SEQUEN, TURNTY,
     -        PHI, PSI, CADIST
 25       FORMAT(1X,3(2X,A5,2X,A3),2X,A3,2X,A7,5X,2(F6.1,2X),F4.1)
          I = I + 1

C----     Form incomplete fields
          DO 400, J = 1, 3
              RES(J) = CHAIN // RSNUM(J)
              SEQ(J) = SEQUEN(J:J)
 400      CONTINUE
CPDBSUM RAL 16/3/05 <--
         IF(I.EQ.1) THEN
CPDBSUM RAL 16/3/05 -->
C            NPS=INDEX(PDBCOD,' ')
C            GNAME=PDBCOD(1:NPS-1)//'_gturns_01.ps'
            GNAME = 'promotif001.ps'
CPDBSUM RAL 16/3/05 <--
            CALL PSOPEN(GNAME)
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'GAMM',NUMFI,GNAME,1)
CPDBSUM RAL 16/3/05 <--
            GT=.TRUE.
CPDBSUM RAL 4/4/05 -->
C         ELSE IF(MOD(I,16).EQ.1) THEN
C            JPS=INDEX(GNAME,'.ps')
C            NUMFI=NUMFI+1
C            if(numfi.lt.10) then
C               WRITE(C,'(I1)') NUMFI
C               GNAME(JPS-1:JPS-1)=C
C            else
C               WRITE(C2,'(I2)') NUMFI
C               GNAME(JPS-2:JPS-1)=C2
C            ENDIF
C            CALL PSOPEN(GNAME)
C            CALL TITLE(PDBCOD,'GAMM',NUMFI,GNAME,1)
C----    If plotting all motifs (ie one-per-page, as opposed to just
C        a single page of several motifs) start a new page
         ELSE IF (.NOT.ONEPAG) THEN
            NUMFI=NUMFI+1
            JPS=INDEX(GNAME,'.ps')
            WRITE(C3,'(I3.3)') NUMFI
            GNAME(JPS-3:JPS-1) = C3
            CALL PSOPEN(GNAME)
CPDBSUM RAL 4/4/05 <--
            X0=XSTART
            Y0=YSTART
         ENDIF
 
C     draw ramachandran 
CPDBSUM RAL 5/4/05 -->
C         IF(.NOT.BW) CALL COLOUR(GAMCOL(1))
C         CALL RAMACH(X0,Y0,X0+RSIZ,Y0+RSIZ)
         CALL RAMACH(X0,Y0,X0+RSIZ,Y0+RSIZ,GAMCOL(1))
CPDBSUM RAL 5/4/05 <--

C     plot phi,psi of i+1
CPDBSUM RAL 5/4/05 -->
C         IF(.NOT.BW) CALL COLOUR(GAMCOL(2))
         IF(.NOT.BW) CALL COLOUR(4)
CPDBSUM RAL 5/4/05 <--
         CALL COORD(PHI,PSI,RSIZ,X0,Y0,XPLOT,YPLOT)
         CALL FPCIRC(XPLOT,YPLOT,CRAD)
CPDBSUM RAL 6/4/05 -->
C         IF(.NOT.BW) CALL COLOUR(GAMCOL(3))
C         CALL PLTNUM(RES(1),RES(3),TURNTY,X0,Y0,RSIZ)
         CALL PLTNM2(RSNAM(1),RSNUM(1),RSNAM(3),RSNUM(3),TURNTY,X0,Y0,
     -       RSIZ,MOTIF)
CPDBSUM RAL 6/4/05 <--

C     plot  schematic diagram of gamma turn
CPDBSUM RAL 6/4/05 -->
         IF(.NOT.BW) CALL COLOUR(GAMCOL(4))
C         CALL PGTURN(X0, Y0-50.0, 50.0, 30.0,CRAD)
C         IF(.NOT.BW) CALL COLOUR(GAMCOL(5))
C         CALL TPRINT(SEQ(1),FSIZ,X0-FSIZ/2.0,Y0-50.0-FSIZ/2.0)
C         CALL TPRINT(SEQ(2),FSIZ,X0+25.0-FSIZ/2.0,Y0-18.0)
C         CALL TPRINT(SEQ(3),FSIZ,X0+50.0+FSIZ/2.0,Y0-50.0-FSIZ/2.0)
         CALL PGTURN(X0, Y0-50.0, RSIZ, 30.0,CRAD)
         CALL COLOUR(3)
         CALL BPRINT(SEQ(1),FSIZ,X0-FSIZ/2.0,Y0-50.0-FSIZ/2.0)
         CALL COLOUR(4)
         CALL BPRINT(SEQ(2),FSIZ,X0+RSIZ/2.0-FSIZ/2.0,Y0-18.0)
         CALL COLOUR(5)
         CALL BPRINT(SEQ(3),FSIZ,X0+RSIZ+FSIZ/2.0,Y0-50.0-FSIZ/2.0)
CPDBSUM RAL 6/4/05 <--

C     plot ca(i)-ca(i+3) distance
         IF(.NOT.BW) CALL COLOUR(GAMCOL(6))
         WRITE(CA(1:3),'(F3.1)') CADIST
CPDBSUM RAL 6/4/05 -->
C         CALL TPRINT(CA,FSIZ,X0+25.0,Y0-44.0)
         CALL TPRINT(CA,FSIZ,X0+RSIZ/2.0,Y0-44.0)
CPDBSUM RAL 6/4/05 <--
         CALL DASH(3)
CPDBSUM RAL 6/4/05 -->
C         CALL PLINE(X0+2.0,Y0-50.0,X0+48.0,Y0-50.0)
         CALL PLINE(X0+2.0,Y0-50.0,X0+RSIZ-2.0,Y0-50.0)
CPDBSUM RAL 6/4/05 <--
         CALL DASH(0)

c     check if at end of line or page and increment coordinates
CPDBSUM RAL 4/4/05 -->
C         IF(MOD(I,16).EQ.0) THEN
         IF(ONEPAG .AND. MOD(I,16).EQ.0) THEN
             GO TO 40
         ELSE IF (.NOT.ONEPAG) THEN
CPDBSUM RAL 4/4/05 <--
            CALL PSCLOSE
         ELSE IF(MOD(I,4).EQ.0) THEN
            Y0=Y0-YDEL
            X0=XSTART
         ELSE
            X0=X0+XDEL
         ENDIF
CPDBSUM RAL 16/3/05 -->
C 30   CONTINUE
         IF (ONEPAG) MOTIF = MOTIF + 1
      GO TO 30
CPDBSUM RAL 16/3/05 <--
 40   IF(GT) CALL PSCLOSE

c     write out table of data
CPDBSUM RAL 16/3/05 -->
C      IF(GT) THEN
C         REWIND(UNIT=2)
C         GNAMTB=PDBCOD(1:NPS-1)//'_gturn_tab01.ps'
C         CALL PSOPEN(GNAMTB)
C         CALL TITLE(PDBCOD,'GAMM',NUMTFI,GNAMTB,2)
C         
C         Y=640.0
C 110     READ(2,'(A80)',end=890) LINE
C         IF(LINE(2:6).NE.'Gamma') GO TO 110
C         
C         DO 50  I=1,NNN
C            NUM=I
C            READ(2,15,END=60) RES(1),SEQ(1),RES(2),SEQ(2),RES(3),SEQ(3),
C     ~           TURNTY,PHI,PSI,CADIST
C            IF(RES(1).EQ.'      ') GO TO 50
C            IF(Y.LT.50.0) THEN
C               CALL GTRNTB(Y)
C               CALL NEWPAG(NUMTFI,GNAMTB,PDBCOD,TTL)
C               Y=640.0
C            ENDIF
C            WRITE(TEXT6,'(A6)') RES(1)
C            CALL TPRINT(TEXT6,12.0,75.0,Y)
C            CALL TPRINT('-',12.0,100.0,Y)
C            WRITE(TEXT6,'(A6)') RES(3)
C            CALL TPRINT(TEXT6,12.0,125.0,Y)
C            WRITE(TEXT6,'(3(A1,1X))') SEQ(1),SEQ(2),SEQ(3)
C            CALL TPRINT(TEXT6,12.0,175.0,Y)
C            CALL TPRINT(TURNTY,12.0,250.0,Y)
C            WRITE(TEXT6,'(F6.1)') PHI
C            CALL TPRINT(TEXT6,12.0,325.0,Y)
C            WRITE(TEXT6,'(F6.1)') PSI
C            CALL TPRINT(TEXT6,12.0,375.0,Y)
C            WRITE(TEXT3,'(F3.1)') CADIST
C            CALL TPRINT(TEXT3,12.0,425.0,Y)
C            Y=Y-20.0
C 50      CONTINUE      
C 60      IF(NUM.GT.1) CALL GTRNTB(Y)
C         IF(NUM.GT.1) CALL PSCLOSE
C      ENDIF
C         
C      CLOSE(UNIT=2)
C      GO TO 900
C 890  WRITE(*,*) ' No gamma turn data'
CPDBSUM RAL 16/3/05 <--

 900  RETURN
      END
C***********************************************************************
C      SUBROUTINE BULGE - To plot out postscript for bulges
C
C**********************************************************
CPDBSUM RAL 16/3/05 -->
C      SUBROUTINE BULGE(PDBCOD)
      SUBROUTINE BULGE(PDBCOD,CHAIN,ONEPAG,MSTART)
CPDBSUM RAL 16/3/05 <--
     
      INCLUDE 'psplot.par'

      INTEGER I,J,K,NUMFI,NUMTFI,NPS,JPS,NUM
      REAL X0,Y0,PHI(5),PSI(5),XDEL,YDEL,RSIZ,XPLOT(5),
     ~     YPLOT(5),FSIZ,XSTART,YSTART,CRAD,SSIZ,AMINT,
     ~     SMAXT,SHDLEN,Y
CPDBSUM RAL 15/3/05 -->
      CHARACTER*1 CHAIN, CHAINS(5)
      CHARACTER*3 C3, RSNAM(5), RSNAME(8), RSNAMP(5)
      CHARACTER*5 RSNUM(5)
      CHARACTER*6 NUMBER
      CHARACTER*(RECLEN) IREC
CPDBSUM RAL 15/3/05 <--
      CHARACTER RES(5)*6,SEQ(5),AP,BTYP,C,
     ~     STYP(2)*12,CLASS(5)*7,CIND*5,SIND*2,TITL1*12,
     ~     TITL2*7,LINE*123,JSEQNO(8)*6,JSEQNOP(5)*6,
     ~     HBP(8,4)*6,SST(8),SSTP(5),TEXT3*3,TEXT6*6,TTL*4,c2*2
      CHARACTER*78 BLGNM2,BNAME,BNAMTB,IBLGNM,PDBCOD
      LOGICAL BG
CPDBSUM RAL 15/3/05 -->
      INTEGER MOTIF, MSTART
      LOGICAL EOF, ONEPAG
      REAL X
CPDBSUM RAL 15/3/05 <--

CPDBSUM RAL 6/4/05 -->
C      DATA XSTART/60.0/,YSTART/660.0/,YDEL/140.0/,RSIZ/50.0/,
      DATA XSTART/60.0/,YSTART/620.0/,YDEL/180.0/,RSIZ/75.0/,
CPDBSUM RAL 6/4/05 <--
     ~     FSIZ/12.0/,CRAD/2.0/,SSIZ/4.0/,AMINT/0.5/,
     ~     SMAXT/3.0/,SHDLEN/3.0/,XDEL/300.0/,
     ~     TTL/'BULG'/
      DATA (STYP(I),I=1,2)/'Antiparallel','  Parallel  '/
      DATA (CLASS(I),I=1,5)/'Classic','Special','  Bent ',
     ~     '   G1  ',' Wide  '/
      DATA CIND/'CSBGW'/,SIND/'AP'/

      BG=.FALSE.
      NUMFI=1
      NUMTFI=1
CPDBSUM RAL 17/3/05 -->
C      NPS=INDEX(PDBCOD,' ')
C      IBLGNM=PDBCOD(1:NPS-1)//'.blg'
C      OPEN(UNIT=2,FILE=IBLGNM,STATUS='OLD',ERR=900)
C      BLGNM2=PDBCOD(1:NPS-1)//'.blgplt'
C      OPEN(UNIT=3,FILE=BLGNM2)
CPDBSUM RAL 17/3/05 <--

c      READ(2,*)

C     plot a ramachandran and schematic diagram for each bulge

      X0=XSTART
      Y0=YSTART
CPDBSUM RAL 17/3/05 -->
C     DO 60 I=1,NNN
C         READ(2,'(A123)',END=70) LINE
      I = 0
      IF (ONEPAG) THEN
          MOTIF = MSTART
      ELSE
          MOTIF = 0
      ENDIF
 60   CONTINUE

C----     Get the next record from the promotif.dat file
          CALL GETREC(FILDAT,CHAIN,IREC,EOF)
          IF (EOF) GO TO 70

C----     Extract the fields from the line just read in
          READ(IREC,110) (CHAINS(J), RSNUM(J), RSNAM(J), J = 1, 5),
     -         AP, BTYP, (PHI(J), PSI(J), J=1,5)
 110      FORMAT(3X,5(A1,2X,A5,2X,A3,2X),A1,2X,A1,37X,10(F6.1,2X))
          I = I + 1

C----     Form incomplete fields
          DO 400, J = 1, 5
              RES(J) = CHAINS(J) // RSNUM(J)
              CALL NAM2AA(RSNAM(J),SEQ(J))
 400      CONTINUE

C----     Read in and interpret the plot data giving the bridge partners
C         around the bulge
          CALL GETPLT(IREC(234:),JSEQNO,JSEQNOP,RSNAME,RSNAMP,SST,SSTP,
     -        HBP)
CPDBSUM RAL 17/3/05 <--
         IF(I.EQ.1) THEN
CPDBSUM RAL 17/3/05 -->
C            NPS=INDEX(PDBCOD,' ')
C            BNAME=PDBCOD(1:NPS-1)//'_bulges_01.ps'
            BNAME = 'promotif001.ps'
CPDBSUM RAL 17/3/05 <--
            CALL PSOPEN(BNAME)
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'BULG',NUMFI,BNAME,1)
CPDBSUM RAL 16/3/05 <--
            BG=.TRUE.
CPDBSUM RAL 4/4/05 -->
C         ELSE IF(MOD(I,10).EQ.1) THEN
C            JPS=INDEX(BNAME,'.ps')
C            NUMFI=NUMFI+1
C            if(numfi.lt.10) then
C               WRITE(C,'(I1)') NUMFI
C               BNAME(JPS-1:JPS-1)=C
C            else
C               WRITE(C2,'(I2)') NUMFI
C               BNAME(JPS-2:JPS-1)=C2
C            ENDIF
C            CALL PSOPEN(BNAME)
C            CALL TITLE(PDBCOD,'BULG',NUMFI,BNAME,1)
C----    If plotting all motifs (ie one-per-page, as opposed to just
C        a single page of several motifs) start a new page
         ELSE IF (.NOT.ONEPAG) THEN
            NUMFI=NUMFI+1
            JPS=INDEX(BNAME,'.ps')
            WRITE(C3,'(I3.3)') NUMFI
            BNAME(JPS-3:JPS-1) = C3
            CALL PSOPEN(BNAME)
CPDBSUM RAL 4/4/05 <--
            X0=XSTART
            Y0=YSTART
         ENDIF
CPDBSUM RAL 17/3/05 -->
C         IF(LINE(1:1).EQ.'A'.OR.LINE(1:1).EQ.'P') THEN
C            READ(LINE,15,END=70) AP,BTYP,
C     ~           (RES(J),SEQ(J),PHI(J),PSI(J),J=1,5)
C 15         FORMAT(A1,1X,A1,5(1X,A6,1X,A1,2X,F6.1,1X,F6.1))
         IF(AP.EQ.'A' .OR. AP.EQ.'P') THEN
CPDBSUM RAL 17/3/05 <--
            IF(RES(2).EQ.'      ') GO TO 70

C     read in bridge partners around the bulge
CPDBSUM RAL 18/3/05 -->
C            DO 16 K=1,8
C               JSEQNO(K)='      '
C               SST(K)=' '
C               DO 116 J=1,4
C                  HBP(K,J)='     '
C 116           CONTINUE
C 16         CONTINUE
C            DO 17 K=1,5
C               JSEQNOP(K)='      '
C               SSTP(K)=' '
C 17         CONTINUE
C            READ(3,*)
C            DO 20 K=1,8
C               READ(3,25,end=221) JSEQNO(K),SST(K),(HBP(K,J),J=1,4)
C 25            FORMAT(7X,A6,1X,A1,1X,4(A6,1X))
C               IF(JSEQNO(K).EQ.'      ') GO TO 21
C               IF(K.EQ.8) READ(3,*)
C 20         CONTINUE
C 21         DO 22 K=1,5
C               READ(3,26,end=221) JSEQNOP(K),SSTP(K)
C 26            FORMAT(7X,A6,1X,A1)
C               IF(JSEQNOP(K).EQ.'      ') GO TO 221
C 22         CONTINUE
C            read(3,*)
CPDBSUM RAL 18/3/05 <--

C     draw ramachandran
 221        J=INDEX(SIND,AP)
            TITL1=STYP(J)
            K=INDEX(CIND,BTYP)
            TITL2=CLASS(K)
CPDBSUM RAL 5/4/05 -->
C            IF(.NOT.BW) CALL COLOUR(BULCOL(1))
C            CALL RAMACH(X0, Y0, X0+RSIZ, Y0+RSIZ)
            CALL RAMACH(X0, Y0, X0+RSIZ, Y0+RSIZ,BULCOL(1))
CPDBSUM RAL 5/4/05 <--
c     plot phi,psi of x
CPDBSUM RAL 6/4/05 -->
C            IF(.NOT.BW) CALL COLOUR(BULCOL(2))
            IF(.NOT.BW) CALL COLOUR(3)
CPDBSUM RAL 6/4/05 <--
            IF(PHI(1).NE.-999.9.AND.PSI(1).NE.-999.9) THEN
               CALL COORD(PHI(1),PSI(1),RSIZ,X0,Y0,XPLOT(1),YPLOT(1))
               CALL FPCIRC(XPLOT(1),YPLOT(1),CRAD)
            ENDIF
c     plot phi,psi of 1 and 2
            DO 30 J=2,3
               IF(PHI(J).NE.-999.9.AND.PSI(J).NE.-999.9) THEN
CPDBSUM RAL 6/4/05 -->
                  IF (J.EQ.2) THEN
                      CALL COLOUR(4)
                  ELSE
                      CALL COLOUR(5)
                  ENDIF
CPDBSUM RAL 6/4/05 <--
                  CALL COORD(PHI(J),PSI(J),RSIZ,X0,Y0,XPLOT(J),YPLOT(J))
                  CALL FPSQUAR(XPLOT(J),YPLOT(J),SSIZ)
               ENDIF
 30         CONTINUE
            IF(.NOT.BW) CALL COLOUR(BULCOL(5))
            IF(PHI(2).NE.-999.9.AND.PSI(2).NE.-999.9.AND.
     ~           PHI(3).NE.-999.9.AND.PSI(3).NE.-999.9)
     ~           CALL FARROW(XPLOT(2),YPLOT(2),XPLOT(3),YPLOT(3)
     ~           ,AMINT,SMAXT,SHDLEN)

C     plot phi,psi of 3 and 4 if they exist
            IF(.NOT.BW) CALL COLOUR(BULCOL(3))
            DO 40 J=4,5
               IF(PHI(J).NE.-999.9.AND.PSI(J).NE.-999.9) THEN
                  CALL COORD(PHI(J),PSI(J),RSIZ,X0,Y0,XPLOT(J),YPLOT(J))
                  CALL FPSQUAR(XPLOT(J),YPLOT(J),SSIZ)
               ENDIF
 40         CONTINUE

            IF(.NOT.BW) CALL COLOUR(BULCOL(4))
            CALL TPRINT(TITL1,FSIZ,X0+RSIZ/2,Y0+RSIZ+20.0)
            CALL TPRINT(TITL2,FSIZ,X0+RSIZ/2,Y0+RSIZ+8.0)

C----       Plot motif number, if required
            IF (MOTIF.GT.0) THEN

C----           Format the motif number
                WRITE(NUMBER,10) MOTIF
 10             FORMAT(I5,'.')

C----           Remove leading spaces
                CALL STRIM(NUMBER)

C----           Plot the motif number
                X = X0 - FSIZ / 2.0
                Y = Y0 + RSIZ + 20.0
                CALL COLOUR(1)
                CALL BPRINT(NUMBER,FSIZ,X,Y)
            ENDIF

c     plot a schematic diagram of bulge
            IF(.NOT.BW) CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 -->
C            CALL SCHBUL
C     ~           (X0,Y0,RES,AP,BTYP,HBP,JSEQNO,JSEQNOP,SST,SSTP)
            CALL SCHBUL(X0,Y0,RES,AP,BTYP,HBP,JSEQNO,JSEQNOP,RSNAME,
     -          RSNAMP,SST,SSTP)
CPDBSUM RAL 6/4/05 <--
C     check if end of page -close file and start again
CPDBSUM RAL 4/4/05 -->
C            IF(MOD(I,10).EQ.0.AND.I.GT.0) THEN
            IF(ONEPAG .AND. MOD(I,8).EQ.0) THEN
                GO TO 70
            ELSE IF (.NOT.ONEPAG) THEN
CPDBSUM RAL 4/4/05 <--
               CALL PSCLOSE
C     check if end of line
            ELSE IF(MOD(I,2).EQ.0) THEN
               Y0=Y0-YDEL
               X0=XSTART
C     otherwise just increment X
            ELSE
               X0=X0+XDEL
            ENDIF
         ENDIF
         NUM=I
CPDBSUM RAL 17/3/05 -->
C 60   CONTINUE
         IF (ONEPAG) MOTIF = MOTIF + 1
      GO TO 60
CPDBSUM RAL 17/3/05 <--
 70   IF(BG.AND.MOD(NUM,8).NE.0) CALL PSCLOSE

c     write out table of data
CPDBSUM RAL 17/3/05 -->
C      IF(BG) THEN
C         REWIND(UNIT=2)
C         BNAMTB=PDBCOD(1:NPS-1)//'_bulg_tab01.ps'
C         CALL PSOPEN(BNAMTB)
C         CALL TITLE(PDBCOD,'BULG',NUMTFI,BNAMTB,2)
Cc         READ(2,*)
C         Y=640.0
C         DO 80  I=1,NNN
C            NUM=I
C            READ(2,'(A123)',END=90) LINE
C            IF(LINE(1:1).EQ.'A'.OR.LINE(1:1).EQ.'P') THEN
C               READ(LINE,15,END=70) AP,BTYP,
C     ~              (RES(J),SEQ(J),PHI(J),PSI(J),J=1,5)
C               IF(RES(2).EQ.'      ') GO TO 90
C            ENDIF
C            IF(Y.LT.50.0) THEN
C               CALL BULGTB(Y)
C               CALL NEWPAG(NUMTFI,BNAMTB,PDBCOD,TTL)
C               Y=640.0
C            ENDIF
C            WRITE(TEXT6,'(A6)') RES(1)
C            CALL TPRINT(TEXT6,12.0,75.0,Y)
C            WRITE(TEXT6,'(A6)') RES(2)
C            CALL TPRINT(TEXT6,12.0,125.0,Y)
C            WRITE(TEXT6,'(A6)') RES(3)
C            CALL TPRINT(TEXT6,12.0,175.0,Y)
C            CALL TPRINT(SEQ(1),12.0,210.0,Y)
C            CALL TPRINT(SEQ(2),12.0,230.0,Y)
C            CALL TPRINT(SEQ(3),12.0,250.0,Y)
C            WRITE(TEXT3,'(A1,1X,A1)') AP,BTYP
C            CALL TPRINT(TEXT3,12.0,285.0,Y)
C            WRITE(TEXT6,'(F6.1)') PHI(1)
C            CALL TPRINT(TEXT6,12.0,330.0,Y)
C            WRITE(TEXT6,'(F6.1)') PSI(1)
C            CALL TPRINT(TEXT6,12.0,370.0,Y)
C            WRITE(TEXT6,'(F6.1)') PHI(2)
C            CALL TPRINT(TEXT6,12.0,410.0,Y)
C            WRITE(TEXT6,'(F6.1)') PSI(2)
C            CALL TPRINT(TEXT6,12.0,450.0,Y)
C            WRITE(TEXT6,'(F6.1)') PHI(3)
C            CALL TPRINT(TEXT6,12.0,490.0,Y)
C            WRITE(TEXT6,'(F6.1)') PSI(3)
C            CALL TPRINT(TEXT6,12.0,530.0,Y)
C            Y=Y-20.0
C 80      CONTINUE      
C 90      IF(NUM.GT.1) CALL BULGTB(Y)
Cc         IF(NUM.GT.1) CALL PSCLOSE
C      ENDIF
C         
C      CLOSE(UNIT=2)
CPDBSUM RAL 17/3/05 <--
      
 900  RETURN
      END
C******************************************************************
c     SUBROUTINE SCHBUL - Plots schematic diagram of bulge
c
c******************************************************************
CPDBSUM RAL 6/4/05 -->
C      SUBROUTINE SCHBUL
C     ~ (X0,Y0,RES,AP,BTYP,HBP,JSEQNO,JSEQNOP,SST,SSTP)
      SUBROUTINE SCHBUL(X0,Y0,RES,AP,BTYP,HBP,JSEQNO,JSEQNOP,RSNAME,
     -    RSNAMP,SST,SSTP)
CPDBSUM RAL 6/4/05 <--

      include 'psplot.par'

      INTEGER I,K,NRES,APART(6),PPART(6),NP,rbox,lbox,ijk,loopk,
     ~     loopl,l,twost,twostp,j
CPDBSUM RAL 6/4/05 -->
      INTEGER ICOL, JUST
      LOGICAL BOLD
CPDBSUM RAL 6/4/05 <--
      REAL XCOL1,XCOL2,XCOL3,Y0,YDEL,XSIZ,YSIZ,
     ~     Y1,Y2,FSIZSM,YSCH,X1,X2,AMINT,AMAXT,AHDLEN,XCOORD(8),
     ~     YCOORD(8),XCORDP(5),YCORDP(5),rbegx,rbegy,lbegx,lbegy,
     ~     X0
CPDBSUM RAL 6/4/05 -->
      REAL FSIZL, LSEP
      CHARACTER*3 RSNAME(8), RSNAMP(5)
CPDBSUM RAL 6/4/05 <--
      CHARACTER RES(5)*6,AP,BTYP,TEXT*6,JSEQNO(8)*6,JSEQNOP(5)*6,
     ~     HBP(8,4)*6,SST(8),SSTP(5),SHEET*2,NULL*6,BPOS(5)*1

CPDBSUM RAL 6/4/05 -->
C      DATA YDEL/15.0/,XSIZ/15.0/,YSIZ/8.0/,FSIZSM/6.0/,AMINT/0.5/,
      DATA YDEL/22.0/,XSIZ/30.0/,YSIZ/12.0/,FSIZSM/9.0/,AMINT/0.5/,
     -     FSIZL/12.0/,LSEP/22.0/,
CPDBSUM RAL 6/4/05 <--
     ~     AMAXT/4.0/,AHDLEN/5.0/,SHEET/'Ee'/,NULL/'      '/,
     ~     BPOS/'X','1','2','3','4'/

      DATA (APART(I),I=1,6)/5,4,3,0,2,1/
      DATA (PPART(I),I=1,6)/1,2,3,0,4,5/
      
CPDBSUM RAL 6/4/05 -->
C      XCOL1=X0+110.0
C      XCOL2=X0+150.0
C      XCOL3=X0+170.0
      BOLD = .FALSE.
      JUST = 0
      XCOL1=X0+120.0
      XCOL2=X0+170.0
      XCOL3=X0+200.0
CPDBSUM RAL 6/4/05 <--
c     first zero coordinate arrays
      DO 5 K=1,8
         XCOORD(K)=0.0
         YCOORD(K)=0.0
 5    CONTINUE
      DO 6 K=1,5
         XCORDP(K)=0.0
         YCORDP(K)=0.0
 6    CONTINUE
      rbox=0
      lbox=0
      YSCH=Y0-2.0*YDEL

c14/10/96 -search for any non-strand residues - don't plot any residues in 
c     strand before these as they may be in a different strand
c     1. check for more than one strand
      twost=0
      twostp=0
      do 750 k=3,5
         if(index(sheet,sstp(k)).ne.0) then
            do 749 j=k-1,1,-1
               if(index(sheet,sstp(k)).eq.0) then
                  do 748 l=j-1,1,-1
                     if(index(sheet,sstp(l)).ne.0) twostp=1
 748              continue
               endif
 749        continue
         endif
 750  continue
      do 760 k=3,8
         if(index(sheet,sst(k)).ne.0) then
            do 759 j=k-1,1,-1
               if(index(sheet,sst(j)).eq.0) then
                  do 758 l=j-1,1,-1
                     if(index(sheet,sst(l)).ne.0) twost=1
 758              continue
               endif
 759        continue
         endif
 760  continue
      if(twost.eq.1) then
         loopk=0
         do 555 k=1,3
            if(index(sheet,sst(k)).eq.0) then
               loopk=k
               go to 556
            endif
 555     continue
 556     if(loopk.gt.1) then
            do 557 k=1,loopk
               sst(k)=' '
 557        continue
         endif
         loopk=8
         do 565 k=8,5,-1
            if(index(sheet,sst(k)).eq.0.and.jseqno(k).ne.'      ')
     ~           then
               loopk=k
               go to 566
            endif
 565     continue
 566     if(loopk.lt.8) then
            do 567 k=loopk,8
               sst(k)=' '
 567        continue
         endif
      endif
      if(twostp.eq.1) then
         loopl=5
         do 558 l=5,3,-1
            if(index(sheet,sstp(l)).eq.0) then
               loopk=l
               go to 559
            endif
 558     continue
 559     if(loopl.lt.5) then
            do 560 l=loopl,5
               sstp(l)=' '
 560        continue
         endif
         loopl=0
         do 568 l=1,3
            if(index(sheet,sstp(l)).eq.0) then
               loopl=l
               go to 569
            endif
 568     continue
 569     if(loopl.gt.1) then
            do 570 l=1,loopl
               sstp(l)=' '
 570        continue
         endif
      endif
      IF(BTYP.EQ.'C'.OR.BTYP.EQ.'W') THEN
         DO 10 K=1,6
            YSCH=YSCH+YDEL
            X1=XCOL1
            Y2=YSCH
            IF(K.LT.3.OR.K.GT.4) THEN
               Y1=YSCH
               X2=XCOL2
            ELSE
               X2=XCOL3
               IF(K.EQ.3) THEN
                  Y1=YSCH+YDEL/2.0
               ELSE
                  Y1=YSCH
               ENDIF
            ENDIF
            TEXT=NULL
            IF(INDEX(SHEET,SST(K)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C               CALL PRECT(X2,Y2,XSIZ,YSIZ)
C               rbox=rbox+1
C               WRITE(TEXT,'(A6)') JSEQNO(K)
C               CALL TPRINT(TEXT,FSIZSM,X2,Y2)
               ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
               DO 7 I=2,5
                  IF(JSEQNO(K).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                     IF(.NOT.BW) THEN
C                        IF(I.LT.4) THEN
C                           CALL COLOUR(BULCOL(4))
C                        ELSE
C                           CALL COLOUR(BULCOL(3))
C                        ENDIF
C                     ENDIF
C                     CALL TPRINT(BPOS(I),FSIZSM,X2+17.0,Y2)
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                     CALL BLABEL(I,FSIZL,X2+LSEP,Y2,ICOL)
CPDBSUM RAL 6/4/05 <--
                  ENDIF
 7             CONTINUE
CPDBSUM RAL 6/4/05 -->
               CALL PRESNM(JSEQNO(K),RSNAME(K),FSIZSM,X2,Y2,ICOL,BOLD,
     -             JUST)
               CALL PRECT(X2,Y2,XSIZ,YSIZ)
               rbox=rbox+1
               CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
               XCOORD(K)=X2
               YCOORD(K)=Y2
               if(rbox.gt.1) call 
     ~              pline(x2,y2-ysiz/2.0,rbegx,rbegy)
               rbegx=x2
               rbegy=y2+ysiz/2.0
               if(k.eq.6) then
                  call FARROW(XCOL2,Y2+YSIZ/2.0,X2,
     ~                 Y2+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               else if(index(sheet,sst(k+1)).eq.0) then
CPDBSUM RAL 14/4/05 -->
C                  call FARROW(XCOL2,Y2+YSIZ/2.0,X2,
                  call FARROW(X2,Y2+YSIZ/2.0,XCOL2,
CPDBSUM RAL 14/4/05 <--
     ~                 Y2+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               endif
            ENDIF
            TEXT=NULL
            IF(K.NE.4) THEN
               IF(AP.EQ.'A') THEN
                  IF(INDEX(SHEET,SSTP(APART(K))).NE.0) THEN
                     XCORDP(APART(K))=X1
                     YCORDP(APART(K))=Y1
CPDBSUM RAL 6/4/05 -->
C                     WRITE(TEXT,'(A6)') JSEQNOP(APART(K))
C                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
                     ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
c                     IF(JSEQNOP(APART(K)).EQ.RES(1)) THEN
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(4))
c                        CALL TPRINT(BPOS(1),FSIZSM,X1-17.0,Y1)
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(5))
c                     ENDIF
                     DO 117 I=1,3
                        IF(JSEQNOP(APART(K)).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                           IF(.NOT.BW) THEN
C                              IF(I.EQ.1) THEN
C                                 CALL COLOUR(BULCOL(2))
C                              ELSE
C                                 CALL COLOUR(BULCOL(4))
C                              ENDIF
C                           ENDIF
C                           CALL TPRINT(BPOS(I),FSIZSM,X1-17.0,Y1)
c                           IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                           CALL BLABEL(I,FSIZL,X1-LSEP,Y1,ICOL)
CPDBSUM RAL 6/4/05 <--
                        ENDIF
 117                 CONTINUE
CPDBSUM RAL 6/4/05 -->
                     CALL PRESNM(JSEQNOP(APART(K)),RSNAMP(APART(K)),
     -                   FSIZSM,X1,Y1,ICOL,BOLD,JUST)
                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
                     CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
                     lbox=lbox+1
                     if(lbox.gt.1) call 
     ~                    pline(x1,y1-ysiz/2.0,lbegx,lbegy)
                     lbegx=x1
                     lbegy=y1+ysiz/2.0
                     if(k.eq.1) then
                        call FARROW(X1,Y1-YSIZ/2.0,X1,
     ~                       Y1-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                        else if(k.lt.5) then
                          IF(INDEX(SHEET,SSTP(APART(K-1)))
     ~                       .EQ.0) 
     ~                          CALL FARROW(X1,Y1-YSIZ/2.0,X1,
     ~                          Y1-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                        endif
                  ENDIF
               ELSE
                  IF(INDEX(SHEET,SSTP(PPART(K))).NE.0) THEN
                     XCORDP(PPART(K))=X1
                     YCORDP(PPART(K))=Y1
CPDBSUM RAL 6/4/05 -->
C                     WRITE(TEXT,'(A6)') JSEQNOP(PPART(K))
C                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
                     ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
c                     IF(JSEQNOP(PPART(K)).EQ.RES(1)) THEN
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(4))
c                        CALL TPRINT(BPOS(1),FSIZSM,X1-17.0,Y1)
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(5))
c                     ENDIF
                      DO 217 I=1,3
                        IF(JSEQNOP(PPART(K)).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                           IF(.NOT.BW) THEN
C                              IF(I.EQ.1) THEN
C                                 CALL COLOUR(BULCOL(2))
C                              ELSE
C                                 CALL COLOUR(BULCOL(4))
C                              ENDIF
C                           ENDIF
C                           CALL TPRINT(BPOS(I),FSIZSM,X1-17.0,Y1)
C                           IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                           CALL BLABEL(I,FSIZL,X1-LSEP,Y1,ICOL)
CPDBSUM RAL 6/4/05 <--
                        ENDIF
 217                 CONTINUE
CPDBSUM RAL 6/4/05 -->
                     CALL PRESNM(JSEQNOP(PPART(K)),RSNAMP(PPART(K)),
     -                   FSIZSM,X1,Y1,ICOL,BOLD,JUST)
                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
                     CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
                     lbox=lbox+1
                     if(lbox.gt.1) call 
     ~                    pline(x1,y1-ysiz/2.0,lbegx,lbegy)
                     lbegx=x1
                     lbegy=y1+ysiz/2.0
                     if(k.eq.6) then
                        call FARROW(X1,Y1+YSIZ/2.0,X1,
     ~                       Y1+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     else if(k.ne.3) then
                        IF(INDEX(SHEET,SSTP(PPART(K+1))).
     ~                       EQ.0) THEN 
                           CALL FARROW(X1,Y1+YSIZ/2.0,X1,
     ~                          Y1+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                          AHDLEN)
                        endif
                     ENDIF
                  ENDIF
               ENDIF
               IF(TEXT.NE.NULL) CALL TPRINT(TEXT,FSIZSM,X1,Y1)
            ENDIF
 10      CONTINUE
C     now look for hydrogen bonds
         CALL BHBOND(HBP,JSEQNOP,XCORDP,YCORDP,XCOORD,YCOORD,AP)
      ELSE IF(BTYP.EQ.'B') THEN
          DO 20 K=1,5
            YSCH=YSCH+YDEL
            IF(K.EQ.3) THEN
               X1=XCOL1-20.0
               X2=XCOL3
            ELSE
               X1=XCOL1
               X2=XCOL2
            ENDIF
            TEXT=NULL
            IF(INDEX(SHEET,SST(K)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C               CALL PRECT(X2,YSCH,XSIZ,YSIZ)
C               rbox=rbox+1
C               WRITE(TEXT,'(A6)') JSEQNO(K)
C               CALL TPRINT(TEXT,FSIZSM,X2,YSCH)
               ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
               DO 17 I=2,5
                  IF(JSEQNO(K).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(4))
C                     CALL TPRINT(BPOS(I),FSIZSM,X2+17.0,YSCH)
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                     CALL BLABEL(I,FSIZL,X2+LSEP,YSCH,ICOL)
CPDBSUM RAL 6/4/05 <--
                  ENDIF
 17            CONTINUE
CPDBSUM RAL 6/4/05 -->
               CALL PRESNM(JSEQNO(K),RSNAME(K),FSIZSM,X2,YSCH,ICOL,BOLD,
     -             JUST)
               CALL PRECT(X2,YSCH,XSIZ,YSIZ)
               rbox=rbox+1
               CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
               XCOORD(K)=X2
               YCOORD(K)=YSCH
               if(rbox.gt.1) call
     ~              pline(x2,ysch-ysiz/2.0,rbegx,rbegy)
               rbegx=x2
               rbegy=ysch+ysiz/2.0
               if(k.eq.6) then
                  call FARROW(X2,YSCH+YSIZ/2.0,X2,
     ~                 YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               else if(index(sheet,sst(k+1)).eq.0) then
                  call FARROW(X2,YSCH+YSIZ/2.0,X2,
     ~                 YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               endif
            ENDIF
            TEXT=NULL
            IF(AP.EQ.'A') THEN
               IF(INDEX(SHEET,SSTP(6-K)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C                  CALL PRECT(X1,YSCH,XSIZ,YSIZ)
C                  lbox=lbox+1
C                  WRITE(TEXT,'(A6)') JSEQNOP(6-K)
                  ICOL = BULCOL(5)
                  IF(JSEQNOP(6-K).EQ.RES(2)) THEN
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(4))
C                     CALL TPRINT(BPOS(2),FSIZSM,X1-17.0,YSCH)
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                     CALL BLABEL(2,FSIZL,X1-LSEP,YSCH,ICOL)
CPDBSUM RAL 6/4/05 <--
                  ENDIF
CPDBSUM RAL 6/4/05 -->
                  CALL PRESNM(JSEQNOP(6-K),RSNAMP(6-K),FSIZSM,X1,YSCH,
     -                ICOL,BOLD,JUST)
                  CALL PRECT(X1,YSCH,XSIZ,YSIZ)
                  lbox=lbox+1
                  CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
                  XCORDP(6-K)=X1
                  YCORDP(6-K)=YSCH
                  if(lbox.gt.1) call
     ~                 pline(x1,ysch-ysiz/2.0,lbegx,lbegy)
                  lbegx=x1
                  lbegy=ysch+ysiz/2.0
                  IF(k.eq.1) then
                     call FARROW(X1,YSCH-YSIZ/2.0,X1,
     ~                    YSCH-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                    AHDLEN)
                     else if(k.lt.5) then
                  IF(INDEX(SHEET,SSTP(APART(K-1)))
     ~                    .EQ.0) 
     ~                CALL FARROW(X1,YSCH-YSIZ/2.0,X1,
     ~                    YSCH-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                    AHDLEN)
                  endif
               ENDIF
            ELSE
c     parallel
               IF(INDEX(SHEET,SSTP(K)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C                  CALL PRECT(X1,YSCH,XSIZ,YSIZ)
C                  lbox=lbox+1
C                  WRITE(TEXT,'(A6)') JSEQNOP(K)
                  ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
                  DO 817 I=1,3
                  IF(JSEQNOP(K).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(4))
C                     CALL TPRINT(BPOS(I),FSIZSM,X1-17.0,YSCH)
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                     CALL BLABEL(I,FSIZL,X1-LSEP,YSCH,ICOL)
CPDBSUM RAL 6/4/05 <--
                  ENDIF
 817              CONTINUE
CPDBSUM RAL 6/4/05 -->
                  CALL PRESNM(JSEQNOP(K),RSNAMP(K),FSIZSM,X1,YSCH,ICOL,
     -                BOLD,JUST)
                  CALL PRECT(X1,YSCH,XSIZ,YSIZ)
                  lbox=lbox+1
                  CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
                  XCORDP(K)=X1
                  YCORDP(K)=YSCH
                  if(lbox.gt.1) call
     ~                 pline(x1,ysch-ysiz/2.0,lbegx,lbegy)
                  lbegx=x1
                  lbegy=ysch+ysiz/2.0
                  IF(k.eq.6) then
                     call FARROW(X1,YSCH+YSIZ/2.0,X1,
     ~                    YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                    AHDLEN)
                  ELSE 
                     if(k.ne.3) then
                        IF(INDEX(SHEET,SSTP(PPART(K+1))).
     ~                       EQ.0) 
     ~                      CALL FARROW(X1,YSCH+YSIZ/2.0,X1,
     ~                          YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                          AHDLEN)
                     endif
                  endif
               ENDIF
            ENDIF
            IF(TEXT.NE.NULL) CALL TPRINT(TEXT,FSIZSM,X1,YSCH)
 20      CONTINUE
         CALL BHBOND(HBP,JSEQNOP,XCORDP,YCORDP,XCOORD,YCOORD,AP)
      ELSE IF(BTYP.EQ.'G') THEN
         DO 30 K=1,6
            YSCH=YSCH+YDEL
            X1=XCOL1
            IF(K.LT.3.OR.K.GT.4) THEN
               X2=XCOL2
            ELSE
               X2=XCOL3
            ENDIF
            IF(K.EQ.3) THEN
               Y1=YSCH+YDEL/2.0
            ELSE
               Y1=YSCH
            ENDIF
            TEXT=NULL
            IF(INDEX(SHEET,SST(K)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C               IF(K.LE.3) THEN
C                  CALL PCIRC(X2,YSCH,2.0*YSIZ/3.0)
C               ELSE
C                  CALL PRECT(X2,YSCH,XSIZ,YSIZ)
C               ENDIF
C               rbox=rbox+1
C               WRITE(TEXT,'(A6)') JSEQNO(K)
C               CALL TPRINT(TEXT,FSIZSM,X2,YSCH)
               ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
               DO  27 I=2,5
                  IF(JSEQNO(K).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                     IF(.NOT.BW) THEN
C                        IF(I.LT.4) THEN
C                           CALL COLOUR(BULCOL(4))
C                        ELSE
C                           CALL COLOUR(BULCOl(3))
C                        ENDIF
C                     ENDIF
C                     CALL TPRINT(BPOS(I),FSIZSM,X2+17.0,YSCH)
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                     CALL BLABEL(I,FSIZL,X2+LSEP,YSCH,ICOL)
CPDBSUM RAL 6/4/05 <--
                  ENDIF
 27            CONTINUE
CPDBSUM RAL 6/4/05 -->
               CALL PRESNM(JSEQNO(K),RSNAME(K),FSIZSM,X2,YSCH,ICOL,BOLD,
     -             JUST)
               IF(K.LE.3) THEN
                  CALL PCIRC(X2,YSCH,2.0*YSIZ/3.0)
               ELSE
                  CALL PRECT(X2,YSCH,XSIZ,YSIZ)
               ENDIF
               rbox=rbox+1
               CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
               XCOORD(K)=X2
               YCOORD(K)=YSCH
               if(rbox.gt.1) call
     ~              pline(x2,ysch-ysiz/2.0,rbegx,rbegy)
               rbegx=x2
               rbegy=ysch+ysiz/2.0
               IF(k.eq.6) then
                  call FARROW(X2,YSCH+YSIZ/2.0,X2,
     ~                 YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               else if(index(sheet,sst(k+1)).eq.0) then
                  call FARROW(X2,YSCH+YSIZ/2.0,X2,
     ~                 YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               endif
            ENDIF
            TEXT=NULL
            IF(K.NE.4) THEN
               IF(AP.EQ.'A') THEN
                  IF(INDEX(SHEET,SSTP(APART(K))).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
C                     lbox=lbox+1
C                     WRITE(TEXT,'(A6)') JSEQNOP(APART(K))
c                     IF(JSEQNOP(APART(K)).EQ.RES(1)) THEN
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(4))
c                        CALL TPRINT(BPOS(1),FSIZSM,X1-17.0,Y1)
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(5))
c                     ENDIF
                     ICOL = BULCOL(5)
                     DO 317 I=1,3
                        IF(JSEQNOP(APART(K)).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                           IF(.NOT.BW) THEN
C                              IF(I.EQ.1) THEN
C                                 CALL COLOUR(BULCOL(2))
C                              ELSE
C                                 CALL COLOUR(BULCOL(4))
C                              ENDIF
C                           ENDIF
C                           CALL TPRINT(BPOS(I),FSIZSM,X1-17.0,Y1)
C                           IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                           CALL BLABEL(I,FSIZL,X1-LSEP,Y1,ICOL)
CPDBSUM RAL 6/4/05 <--
                        ENDIF
 317                 CONTINUE
CPDBSUM RAL 6/4/05 -->
                     CALL PRESNM(JSEQNOP(APART(K)),RSNAMP(APART(K)),
     -                   FSIZSM,X1,Y1,ICOL,BOLD,JUST)
                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
                     lbox=lbox+1
                     CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--

                     XCORDP(APART(K))=X1
                     YCORDP(APART(K))=Y1
                     if(lbox.gt.1) call
     ~                    pline(x1,y1-ysiz/2.0,lbegx,lbegy)
                     lbegx=x1
                     lbegy=y1+ysiz/2.0
                     IF(k.eq.1) then
                        call FARROW(X1,Y1-YSIZ/2.0,X1,
     ~                       Y1-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                        else if(k.lt.5) then
                      IF(INDEX(SHEET,SSTP(APART(K-1)))
     ~                       .EQ.0) 
     ~                        CALL FARROW(X1,Y1-YSIZ/2.0,X1,
     ~                       Y1-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     endif
                  ENDIF
               ELSE
                  IF(INDEX(SHEET,SSTP(PPART(K))).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
C                     lbox=lbox+1
C                     WRITE(TEXT,'(A6)') JSEQNOP(PPART(K))
                     ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
c                     IF(JSEQNOP(PPART(K)).EQ.RES(1)) THEN
c                        IF(.NOT.BW) CALL COLOUR(BULCOL(4))
c                        CALL TPRINT(BPOS(1),FSIZSM,X1-17.0,Y1)
cc                        IF(.NOT.BW) CALL COLOUR(BULCOL(5))
c                     ENDIF
                     DO 417 I=1,3
                        IF(JSEQNOP(PPART(K)).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                           IF(.NOT.BW) THEN
C                              IF(I.EQ.1) THEN
C                                 CALL COLOUR(BULCOL(2))
C                              ELSE
C                                 CALL COLOUR(BULCOL(4))
C                              ENDIF
C                           ENDIF
C                           CALL TPRINT(BPOS(I),FSIZSM,X1-17.0,Y1)
C                           IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                           CALL BLABEL(I,FSIZL,X1-LSEP,Y1,ICOL)
CPDBSUM RAL 6/4/05 <--
                        ENDIF
 417                 CONTINUE
CPDBSUM RAL 6/4/05 -->
                     CALL PRESNM(JSEQNOP(PPART(K)),RSNAMP(PPART(K)),
     -                   FSIZSM,X1,Y1,ICOL,BOLD,JUST)
                     CALL PRECT(X1,Y1,XSIZ,YSIZ)
                     lbox=lbox+1
                     CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--

                     XCORDP(PPART(K))=X1
                     YCORDP(PPART(K))=Y1
                     if(lbox.gt.1) call
     ~                    pline(x1,y1-ysiz/2.0,lbegx,lbegy)
                     lbegx=x1
                     lbegy=y1+ysiz/2.0
                     IF(k.eq.6) then
                        call FARROW(X1,Y1+YSIZ/2.0,X1,
     ~                       Y1+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     ELSE IF(INDEX(SHEET,SSTP(PPART(K+1))).
     ~                       EQ.0) THEN 
                        CALL FARROW(X1,Y1+YSIZ/2.0,X1,
     ~                       Y1+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     ENDIF
                  ENDIF
               ENDIF
               IF(TEXT.NE.NULL) CALL TPRINT(TEXT,FSIZSM,X1,Y1)
            ENDIF
 30      CONTINUE
         CALL BHBOND(HBP,JSEQNOP,XCORDP,YCORDP,XCOORD,YCOORD,AP)
C     special bulges
      ELSE IF(BTYP.EQ.'S') THEN
c     number of residues in  special bulges
         IF(RES(5).EQ.'      ') THEN
            NRES=7
         ELSE
            NRES=8
         ENDIF
         IF(AP.EQ.'A') THEN
            NP=6
         ELSE
            NP=0
         ENDIF
         DO 40 K=1,NRES
            YSCH=YSCH+YDEL
            X1=XCOL1
            IF(K.GE.3.AND.K.LE.NRES-2) THEN
               X2=XCOL3
               IF(K.EQ.4) Y1=YSCH+(MOD(NRES-1,2))*YDEL/2.0
            ELSE
               X2=XCOL2
               Y1=YSCH
            ENDIF
            TEXT=NULL
            IF(INDEX(SHEET,SST(K)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C               CALL PRECT(X2,YSCH,XSIZ,YSIZ)
C               rbox=rbox+1
C               WRITE(TEXT,'(A6)') JSEQNO(K)
C               CALL TPRINT(TEXT,FSIZSM,X2,YSCH)
               ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
               DO 37 I=2,5
                  IF(JSEQNO(K).EQ.RES(I)) THEN
                     IF(.NOT.BW) CALL COLOUR(BULCOL(4))
CPDBSUM RAL 6/4/05 -->
C                     CALL TPRINT(BPOS(I),FSIZSM,X2+17.0,YSCH)
C                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                     CALL BLABEL(I,FSIZL,X2+LSEP,YSCH,ICOL)
CPDBSUM RAL 6/4/05 <--
                  ENDIF
 37            CONTINUE
CPDBSUM RAL 6/4/05 -->
               CALL PRESNM(JSEQNO(K),RSNAME(K),FSIZSM,X2,YSCH,ICOL,BOLD,
     -             JUST)
               CALL PRECT(X2,YSCH,XSIZ,YSIZ)
               rbox=rbox+1
               CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
               XCOORD(K)=X2
               YCOORD(K)=YSCH
               if(rbox.gt.1) call
     ~              pline(x2,ysch-ysiz/2.0,rbegx,rbegy)
               rbegx=x2
               rbegy=ysch+ysiz/2.0
               IF (k.eq.NRES) then
                  call FARROW(X2,YSCH+YSIZ/2.0,X2,
     ~                 YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               else if(index(sheet,sst(k+1)).eq.0) then
                  call FARROW(X2,YSCH+YSIZ/2.0,X2,
     ~                 YSCH+YSIZ/2.0+YDEL,AMINT,AMAXT,AHDLEN)
               endif
            ENDIF
            TEXT=NULL
            IF(K.LT.3.OR.K.GT.NRES-2.OR.K.EQ.4) THEN
               IF(AP.EQ.'A') THEN
                  NP=NP-1
               ELSE
                  NP=NP+1
               ENDIF
               IF(INDEX(SHEET,SSTP(NP)).NE.0) THEN
CPDBSUM RAL 6/4/05 -->
C                  CALL PRECT(X1,Y1,XSIZ,YSIZ)
C                  lbox=lbox+1
C                  WRITE(TEXT,'(A6)') JSEQNOP(NP)
C                  CALL TPRINT(TEXT,FSIZSM,X1,Y1)
                  ICOL = BULCOL(5)
CPDBSUM RAL 6/4/05 <--
c                  IF(JSEQNOP(NP).EQ.RES(1)) THEN
c                     IF(.NOT.BW) CALL COLOUR(BULCOL(4))
c                     CALL TPRINT(BPOS(1),FSIZSM,X1-17.0,Y1)
c                     IF(.NOT.BW) CALL COLOUR(BULCOL(5))
c                  ENDIF
                   DO 517 I=1,3
                     IF(JSEQNOP(NP).EQ.RES(I)) THEN
CPDBSUM RAL 6/4/05 -->
C                        IF(.NOT.BW) THEN
C                           IF(I.EQ.1) THEN
C                              CALL COLOUR(BULCOL(2))
C                           ELSE
C                              CALL COLOUR(BULCOL(4))
C                           ENDIF
C                        ENDIF
C                        CALL TPRINT(BPOS(I),FSIZSM,X1-17.0,Y1)
C                        IF(.NOT.BW) CALL COLOUR(BULCOL(5))
                        CALL BLABEL(I,FSIZL,X1-LSEP,Y1,ICOL)
CPDBSUM RAL 6/4/05 <--
                     ENDIF
 517              CONTINUE
CPDBSUM RAL 6/4/05 -->
                  CALL PRESNM(JSEQNOP(NP),RSNAMP(NP),FSIZSM,X1,Y1,ICOL,
     -                BOLD,JUST)
                  CALL PRECT(X1,Y1,XSIZ,YSIZ)
                  lbox=lbox+1
                  CALL COLOUR(BULCOL(5))
CPDBSUM RAL 6/4/05 <--
                  XCORDP(NP)=X1
                  YCORDP(NP)=Y1
                  if(lbox.gt.1) call
     ~                 pline(x1,y1-ysiz/2.0,lbegx,lbegy)
                  lbegx=x1
                  lbegy=y1+ysiz/2.0
                  IF(AP.EQ.'A') THEN
                     IF(k.eq.1) then
                        call FARROW(X1,Y1-YSIZ/2.0,X1,
     ~                       Y1-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                        else if(k.lt.5) then
                      IF(INDEX(SHEET,SSTP(k-1))
     ~                       .EQ.0) 
     ~                   CALL FARROW(X1,Y1-YSIZ/2.0,X1,
     ~                       Y1-YSIZ/2.0-YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     endif
                  ELSE
                     IF(k.eq.NRES) then
                        call FARROW(X1,Y1+YSIZ/2.0,X1,
     ~                       Y1+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     ELSE IF(INDEX(SHEET,SSTP(np+1)).
     ~                       EQ.0) THEN 
                        CALL FARROW(X1,Y1+YSIZ/2.0,X1,
     ~                       Y1+YSIZ/2.0+YDEL,AMINT,AMAXT,
     ~                       AHDLEN)
                     ENDIF
                  ENDIF
               endif
            endif
 40      CONTINUE
         CALL BHBOND(HBP,JSEQNOP,XCORDP,YCORDP,XCOORD,YCOORD,AP)
      ENDIF
c      CALL FARROW(XCOL2,YSCH+YSIZ/2.0,XCOL2,YSCH+YSIZ/2.0+YDEL,
c     ~     AMINT,AMAXT,AHDLEN)
c      IF(AP.EQ.'P') THEN
c         CALL FARROW(XCOL1,YSCH+YSIZ/2.0,XCOL1,YSCH+YSIZ/2.0+YDEL,
c     ~        AMINT,AMAXT,AHDLEN)
c      ELSE
c         YSCH1=Y0-YDEL
c         CALL FARROW(XCOL1,YSCH1-YSIZ/2.0,XCOL1,YSCH1-YSIZ/2.0-YDEL,
c     ~        AMINT,AMAXT,AHDLEN)
c      ENDIF

      RETURN
      END
C***********************************************************************
C      SUBROUTINE HAIRPN - To plot out postscript for hairpins
C
C**********************************************************
CPDBSUM RAL 24/3/05 -->
      SUBROUTINE HAIRPN(PDBCOD,CHAIN,ONEPAG,MSTART)
CPDBSUM RAL 24/3/05 <--
      INCLUDE 'psplot.par'

CPDBSUM RAL 7/4/05 -->
      INTEGER LEFT, RIGHT, MXLLEN, MXSLEN
      PARAMETER (LEFT = -1,RIGHT = 1,MXLLEN=125,MXSLEN=30)
CPDBSUM RAL 7/4/05 <--

      INTEGER I,J,K,LPLEN,NUMFI,CLASS1,CLASS2,SLEN1,SLEN2,
     ~     HALFLEN,JNAR(6,2),LHB(3,3),NUMTFI,JSNG(3),NPS,JPS,NUM
      REAL X0,Y0,X1,X2,X3,X4,Y3,Y4,YEND1,YEND2,YDEL,
     ~     FSIZ,XSTART,YSTART,MSIZ,HMNWID,HMXWID,
     ~     HHDLN,XLOOP,YLOOP,YWID,AMINT,AMAXT,AHDLEN,XNAR1,XNAR2,
     ~     XLP1,XLP2,YLP1,YLP2,Y,PYLOOP
CPDBSUM RAL 24/3/05 -->
      CHARACTER*1 CHAIN
      CHARACTER*3 C3, RSNAM(4), RSNAME(4)
      CHARACTER*5 RSNUM(4)
      CHARACTER*6 NUMBER, RES(4)
      CHARACTER*(RECLEN) IREC
CPDBSUM RAL 24/3/05 <--
      CHARACTER RES1I*6,RES1E*6,RES2I*6,RES2E*6,C,
CPDBSUM RAL 7/4/05 -->
C     ~     LSEQ*125,TYPE*5,SSEQ1*30,SSEQ2*30,ECLASS*3,TTL*4,C2*2
     -     LSEQ*(MXLLEN),TYPE*5,SSEQ1*(MXSLEN),SSEQ2*(MXSLEN),ECLASS*3,
     -     TTL*4,C2*2
CPDBSUM RAL 7/4/05 <--
      CHARACTER*78 HNAME,HNAMTB,HPNM,PDBCOD
      CHARACTER TEXT2*2,TEXT6*6,TEXT9*9
      LOGICAL HP
CPDBSUM RAL 29/3/05 -->
      LOGICAL BOLD, EOF, EVEN, ONEPAG, UP
      INTEGER IPOS, JPOS, LENSTR, MOTIF, MSTART
      REAL    DY, FMSIZ, HSHIFT, RSHIFT, WIDTHS, WIDTHL, X, X1L, X1R,
     -        X2L, X2R, XCENTR, XMAX, XORIGN, YWDMIN, YWDMAX, YDROP,
     -        YDROPE, YDROPO, Y1, Y2, YSHIFT, YSTORE(5,2)
CPDBSUM RAL 29/3/05 <--

      DATA XSTART/60.0/,YSTART/660.0/,YDEL/150.0/,FSIZ/12.0/,
CPDBSUM RAL 6/4/05 -->
C     ~     MSIZ/8.0/,HMNWID/5.0/,HMXWID/10.0/HHDLN/5.0/,
C     ~     TYPE(3:3)/':'/,YWID/70.0/,AMINT/0.5/,AMAXT/4.0/,
     ~     MSIZ/10.0/,HMNWID/15.0/,HMXWID/22.0/HHDLN/10.0/,
     ~     TYPE(3:3)/':'/,YWID/45.0/,AMINT/0.5/,AMAXT/4.0/,FMSIZ/15.0/,
CPDBSUM RAL 6/4/05 <--
     ~     AHDLEN/5.0/,TEXT9(3:3)/':'/,
     ~     TTL/'HPIN'/
CPDBSUM RAL 6/4/05 -->
      DATA HSHIFT, RSHIFT, XCENTR, YSHIFT / 1.5, 10.0, 300.0, 18.0 /
      DATA YDROPE, YDROPO / 15.0, 11.0 /
CPDBSUM RAL 6/4/05 <--
      
      HP=.FALSE.
      NUMFI=1
      NUMTFI=1
CPDBSUM RAL 24/3/05 -->
C      NPS=INDEX(PDBCOD,' ')
C      HPNM=PDBCOD(1:NPS-1)//'.hpin'
C      OPEN(UNIT=2,FILE=HPNM,STATUS='OLD',ERR=900)
C      READ(2,*,end=890)
C      NPS = 5
C      HPNM = 'pdb1b6c.hpin'
      BOLD = .TRUE.
      DO 5, I = 1, 5
          YSTORE(I,1) = 0.0
          YSTORE(I,2) = 0.0
 5    CONTINUE
CPDBSUM RAL 24/3/05 <--

      X0=XSTART
      X1=X0
      X2=X0+50.0
      Y0=YSTART

c     read in hairpins one by one and draw schematic diagram
CPDBSUM RAL 29/3/05 -->
C      DO 60 I=1,NNN
C         READ(2,15,END=70) RES1I,RES1E,RES2I,RES2E,CLASS1,CLASS2,
C     ~        ECLASS,SLEN1,SLEN2
C         READ(2,16,END=70) SSEQ1
C         READ(2,17,END=70)LSEQ
C         READ(2,18,END=70)SSEQ2
C         READ(2,*)
C         READ(2,19,END=70)((JNAR(K,J),J=1,2),K=1,6),
C     ~        ((LHB(K,J),J=1,3),K=1,3),(JSNG(K),K=1,3)
      I = 0
      IF (ONEPAG) THEN
          MOTIF = MSTART
      ELSE
          MOTIF = 0
      ENDIF
 60   CONTINUE

C----     Get the next record from the promotif.dat file
          CALL GETREC(FILDAT,CHAIN,IREC,EOF)
          IF (EOF) GO TO 70

C----     Extract the fields from the line just read in
          READ(IREC,110) (RSNUM(J), RSNAM(J), J = 1, 4),
     -         CLASS1, CLASS2, ECLASS, SLEN1, SLEN2
 110      FORMAT(3X,4(A5,2X,A3,2X),2(I3,2X),A3,2X,2(I2,2X))
          I = I + 1

C----     Form incomplete fields
          RES1I = CHAIN // RSNUM(1)
          RES1E = CHAIN // RSNUM(2)
          RES2I = CHAIN // RSNUM(3)
          RES2E = CHAIN // RSNUM(4)

C----     Get the sequence of strand 1
          IPOS = 75
          JPOS = INDEX(IREC(IPOS:),'::')
          IF (JPOS.GT.1) THEN
              SSEQ1 = IREC(IPOS:IPOS + JPOS - 2)
              IPOS = IPOS + JPOS + 1
          ELSE
              SSEQ1 = ' '
          ENDIF

C----     Get the sequence of loop
          JPOS = INDEX(IREC(IPOS:),'::')
          IF (JPOS.GT.1) THEN
              LSEQ = IREC(IPOS:IPOS + JPOS - 2)
              IPOS = IPOS + JPOS + 1
          ELSE
              LSEQ = ' '
          ENDIF

C----     Get the sequence of strand 2
          JPOS = INDEX(IREC(IPOS:),'::')
          IF (JPOS.GT.1) THEN
              SSEQ2 = IREC(IPOS:IPOS + JPOS - 2)
              IPOS = IPOS + JPOS + 1
          ELSE
              SSEQ2 = ' '
          ENDIF

C----     Get the plot data
          JPOS = INDEX(IREC(IPOS:),'::')
          IF (JPOS.EQ.79) THEN
              READ(IREC(IPOS:),19,ERR=99) ((JNAR(K,J),J=1,2),K=1,6),
     -            ((LHB(K,J),J=1,3),K=1,3),(JSNG(K),K=1,3)
          ELSE
              PRINT*, '*** ERROR. Length error in plot data at line', I
              PRINT 98, IPOS, JPOS
 98           FORMAT('   IPOS =',I3,'   JPOS =',I3)
              GO TO 60
          ENDIF
          GO TO 100

 99       CONTINUE
          PRINT*, '*** ERROR. Format error in plot data at line', I
          GO TO 60

 100      CONTINUE
CPDBSUM RAL 29/3/05 <--
 15      FORMAT(8X,4(A6,1X),I3,1X,I3,1X,A3,1X,2(I2,1X))
 16      FORMAT(A30)
 17      FORMAT(A80)
 18      FORMAT(A30)
 19      FORMAT(12(I2,1X),9(I1,1X),3(I2,1X))
         IF(I.EQ.1) THEN
CPDBSUM RAL 29/3/05 -->
C            NPS=INDEX(PDBCOD,' ')
C            HNAME=PDBCOD(1:NPS-1)//'_hairpin_01.ps'
            HNAME = 'promotif001.ps'
CPDBSUM RAL 29/3/05 <--
            CALL PSOPEN(HNAME)
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'HPIN',NUMFI,HNAME,1)
CPDBSUM RAL 16/3/05 <--
            HP=.TRUE.
CPDBSUM RAL 4/4/05 -->
C         ELSE IF(MOD(I,5).EQ.1.AND.I.GT.1) THEN
C            JPS=INDEX(HNAME,'.ps')
C            NUMFI=NUMFI+1
C            if(numfi.lt.10) then
C               WRITE(C,'(I1)') NUMFI
C               HNAME(JPS-1:JPS-1)=C
C            else
C               WRITE(C2,'(I2)') NUMFI
C               HNAME(JPS-2:JPS-1)=C2
C            ENDIF
C            CALL PSOPEN(HNAME)
C            CALL TITLE(PDBCOD,'HPIN',NUMFI,HNAME,1)
C----    If plotting all motifs (ie one-per-page, as opposed to just
C        a single page of several motifs) start a new page
         ELSE IF (.NOT.ONEPAG) THEN
            NUMFI=NUMFI+1
            JPS=INDEX(HNAME,'.ps')
            WRITE(C3,'(I3.3)') NUMFI
            HNAME(JPS-3:JPS-1) = C3
            CALL PSOPEN(HNAME)
CPDBSUM RAL 4/4/05 <--
            X0=XSTART
            Y0=YSTART
         ENDIF
         X0=XSTART
CPDBSUM RAL 7/4/05 -->
C         X1=X0
C         X2=X0+50.0
CC     calculate the end points of the arrows - strand 1
C         YEND1=Y0+SLEN1*5.0
CC     calculate end point of arrows - strand 2
C         YEND2=YEND1-SLEN2*5.0
CC     draw arrows
C         IF(.NOT.BW) CALL COLOUR(HPCOL(1))
C         CALL FARROW(X1,Y0,X1,YEND1,HMNWID,HMXWID,HHDLN)
C         CALL FARROW(X2,YEND1,X2,YEND2,HMNWID,HMXWID,HHDLN)
C
CC     write on beginning and end residues
C         IF(.NOT.BW) CALL COLOUR(HPCOL(2))
C         CALL ARROW(X1,Y0,X1,YEND1,HMNWID,HMXWID,HHDLN)
C         CALL ARROW(X2,YEND1,X2,YEND2,HMNWID,HMXWID,HHDLN)
C         CALL PHCIRC((X1+X2)/2.0,YEND1,25.0)
C         CALL TPRINT(RES1I,FSIZ,X1-20.0,Y0)
C         CALL TPRINT(RES1E,FSIZ,X1-20.0,YEND1)
C         CALL TPRINT(RES2I,FSIZ,X2+23.0,YEND1)
C         CALL TPRINT(RES2E,FSIZ,X2+23.0,YEND2)
C        
Cc     add hairpin type
C         IF(.NOT.BW) CALL COLOUR(HPCOL(3))
C         WRITE(TYPE(1:2),'(I2)') CLASS1
C         WRITE(TYPE(4:5),'(I2)') CLASS2
C         CALL TPRINT(TYPE,FSIZ,X0+25.0,YEND1+35.0)
C         IF(ECLASS.NE.'   ') 
C     ~        CALL TPRINT(ECLASS,FSIZ,X0+45.0,YEND1+35.0)
C
c     print out strand and loop sequences
C         IF(.NOT.BW) CALL COLOUR(HPCOL(4))
C         X3=X0+300.0
C         X4=X0+300.0
C         Y3=Y0+YWID
C         Y4=Y0
C         DO 20 J=1,30
C            IF(SSEQ1(31-J:31-J).NE.' ') THEN
C               CALL TPRINT(SSEQ1(31-J:31-J),MSIZ,X3,Y3)
C               X3=X3-10.0
C            ENDIF
C            IF(SSEQ2(J:J).NE.' ') THEN
C               CALL TPRINT(SSEQ2(J:J),MSIZ,X4,Y4)
C               X4=X4-10.0
C            ENDIF
C 20      CONTINUE
C         CALL FARROW(X0+250.0,Y3+20.0,X0+300.0,Y3+20.0,
C     ~        AMINT,AMAXT,AHDLEN)
C         CALL FARROW(X0+300.0,Y4-20.0,X0+250.0,Y4-20.0,
C     ~        AMINT,AMAXT,AHDLEN)
C
Cc     draw in hydrogen bonds
C         CALL DASH(4)
C         IF(.NOT.BW) CALL COLOUR(HPCOL(5))
C         DO 25 J=1,6
C            IF(JNAR(J,1).NE.0.AND.JNAR(J,2).NE.0) THEN
C               XNAR1=X0+300.0-(JNAR(J,1)-1)*10.0
C               XNAR2=X0+300.0-(JNAR(J,2)-1)*10.0
C               CALL FARROW(XNAR1+2.0,Y3-4.0,XNAR2+2.0,Y4+4.0,
C     ~              AMINT,AMAXT,AHDLEN)
C               CALL FARROW(XNAR2-2.0,Y4+4.0,XNAR1-2.0,Y3-4.0,
C     ~              AMINT,AMAXT,AHDLEN)
C            ENDIF
C 25      CONTINUE
C         IF(JSNG(1).NE.0.AND.JSNG(2).NE.0) THEN
C            XNAR1=X0+300.0-(JSNG(1)-1)*10.0
C            XNAR2=X0+300.0-(JSNG(2)-1)*10.0
C            IF(JSNG(3).EQ.1) THEN
C               CALL FARROW(XNAR1+2.0,Y3-4.0,XNAR2+2.0,Y4+4.0,
C     ~              AMINT,AMAXT,AHDLEN)
C            ELSE
C               CALL FARROW(XNAR2-2.0,Y4+4.0,XNAR1-2.0,Y3-4.0,
C     ~              AMINT,AMAXT,AHDLEN)
C            ENDIF
C         ENDIF
C         CALL DASH(0)
C
Cc     draw in loop residues
C         IF(.NOT.BW) CALL COLOUR(HPCOL(6))
C         XLOOP=X0+300.0
C         YLOOP=Y3
Cc         HALFLEN=INT((CLASS2+1)/2.0) 11.1.95
Cc         HALFLEN=INT((CLASS1+1)/2.0) 20.3.95
C         DO 331 J=1,125
C            IF(LSEQ(J:J).EQ.' ') THEN
C               LPLEN=J-1
C               HALFLEN=INT(J/2.0)
C               GO TO 332
C            ENDIF
C 331     CONTINUE
CC         DO 30 J=1,CLASS2 11.1.95
CC         DO 30 J=1,CLASS1 20.3.95
C 332     DO 30 J=1,LPLEN
CC            YLOOP=YLOOP-YWID/(CLASS2+1) 11.1.95
Cc            YLOOP=YLOOP-YWID/(CLASS1+1) 20.3.95
C            YLOOP=YLOOP-YWID/(LPLEN+1) 
C               IF(LPLEN.GT.10) THEN
C                  IF(YLOOP.GE.Y3-YWID/2.0) THEN
C                     IF(YLOOP-Y3+YWID/2.0.LT.6.0) THEN
C                        PYLOOP=Y3-YWID/2.0+6.0
C                     ELSE
C                        PYLOOP=YLOOP
C                     ENDIF
C                  ELSE
C                     IF(Y3-YWID/2.0-YLOOP.LT.6.0) THEN
C                        PYLOOP=Y3-YWID/2.0+6.0
C                     ELSE
C                        PYLOOP=YLOOP
C                     ENDIF
C                  ENDIF
C               ELSE
C                  PYLOOP=YLOOP
C               ENDIF
C               IF(J.LE.HALFLEN) THEN
C                  XLOOP=XLOOP+10.0
CC     ELSE IF(J.EQ.HALFLEN+1.AND.MOD(CLASS2,2).EQ.0) THEN 11.1.95
Cc     ELSE IF(J.EQ.HALFLEN+1.AND.MOD(CLASS1,2).EQ.0) THEN 20.3.95
C               ELSE IF(J.EQ.HALFLEN+1.AND.MOD(LPLEN,2).EQ.0) THEN
C                  XLOOP=XLOOP
CC     ELSE IF(J.NE.HALFLEN+1.OR.MOD(CLASS2,2).NE.0) THEN 11.1.95
Cc     ELSE IF(J.NE.HALFLEN+1.OR.MOD(CLASS1,2).NE.0) THEN 20.3.95
C               ELSE IF(J.NE.HALFLEN+1.OR.MOD(LPLEN,2).NE.0) THEN
C                  XLOOP=XLOOP-10.0
C               ENDIF
C               CALL TPRINT(LSEQ(J:J),MSIZ,XLOOP,PYLOOP)
C 30      CONTINUE
C
CC     draw loop hydrogen bonds
C         CALL DASH(4)
C         IF(.NOT.BW) CALL COLOUR(HPCOL(5))
C         DO 45 J=1,3
C            IF(LHB(J,1).NE.0.or.LHB(J,2).NE.0) THEN
C               XLP1=X0+300.0+LHB(J,1)*10.0
CC               YLP1=Y3-(LHB(J,1)*YWID)/(CLASS2+1) 11.1.95
CC               YLP1=Y3-(LHB(J,1)*YWID)/(CLASS1+1) 20.3.95
C               YLP1=Y3-(LHB(J,1)*YWID)/(LPLEN+1)
C               XLP2=X0+300.0+LHB(J,2)*10.0
CC               YLP2=Y4+(LHB(J,2)*YWID)/(CLASS2+1) 11.1.95
CC               YLP2=Y4+(LHB(J,2)*YWID)/(CLASS1+1) 20.3.95
C               YLP2=Y4+(LHB(J,2)*YWID)/(LPLEN+1)
C               IF(LHB(J,3).EQ.1) THEN
C                  CALL FARROW(XLP1-2.0,YLP1-4.0,XLP2-2.0,
C     ~                 YLP2+4.0,AMINT,AMAXT,AHDLEN)
C               ELSE
C                  CALL FARROW(XLP2+2.0,YLP2+4.0,XLP1+2.0,
C     ~                 YLP1-4.0,AMINT,AMAXT,AHDLEN)
C               ENDIF
C            ENDIF
C 45      CONTINUE
C         CALL DASH(0)

CPDBSUM RAL 6/4/05 -->

C----    Calculate approx width of hairpin strands plus loop
         IF (SLEN1.GE.SLEN2) THEN
             WIDTHS = (SLEN1 + 2) * RSHIFT
         ELSE
             WIDTHS = (SLEN2 + 2) * RSHIFT
         ENDIF
         WIDTHL = (LENSTR(LSEQ) / 2.0 + 2) * RSHIFT

C----    Determine x-position of hairpin origin (ie at loop start- and
C        end-points)
         XORIGN = XCENTR - (WIDTHS + WIDTHL) / 2.0 + WIDTHS

C----    Print out the hairpin motif, including sequences
         CALL COLOUR(HPCOL(4))
         X3 = XORIGN
         X4 = XORIGN
         Y3 = Y0+YWID
         Y4 = Y0

C----    Initialise arrow end-points
         X1L = X3
         X1R = X3
         X2L = X4
         X2R = X4

C----    Print sequences of strands 1 and 2
         DO 20 J = 1, MXSLEN

C----       Strand 1
            IF (SSEQ1(31-J:31-J).NE.' ') THEN
               CALL COLOUR(5)
               CALL TPRINT(SSEQ1(31-J:31-J),MSIZ,X3,Y3)
               X1L = X3
               X3 = X3 - RSHIFT
            ENDIF

C----       Strand 2
            IF (SSEQ2(J:J).NE.' ') THEN
               CALL COLOUR(4)
               CALL TPRINT(SSEQ2(J:J),MSIZ,X4,Y4)
               X2L = X4
               X4 = X4 - RSHIFT
            ENDIF
 20      CONTINUE

C----    Draw the arrows representing the strands
         CALL PSLWID(1.0)
         X1 = X1L - 3.0 * RSHIFT / 4.0
         X2 = X1R + RSHIFT
         CALL COLOUR(5)
         CALL ARROW(X1,Y3,X2,Y3,HMNWID,HMXWID,HHDLN)
         X2 = X2L - RSHIFT
         X1 = X2R + 3.0 * RSHIFT / 4.0
         CALL COLOUR(4)
         CALL ARROW(X1,Y4,X2,Y4,HMNWID,HMXWID,HHDLN)
         CALL PSLWID(0.1)

C----    Label arrows with start and end residues
         X1 = X1L + RSHIFT / 2.0
         X2 = X1R - RSHIFT / 2.0
         CALL PRESNM(RES1I,RSNAM(1),FSIZ,X1,Y3+YSHIFT,5,BOLD,RIGHT)
         CALL PRESNM(RES1E,RSNAM(2),FSIZ,X2,Y3+YSHIFT,5,BOLD,LEFT)
         X1 = X2R - RSHIFT / 2.0
         X2 = X2L + RSHIFT / 2.0
         CALL PRESNM(RES2I,RSNAM(3),FSIZ,X1,Y4-YSHIFT,4,BOLD,LEFT)
         CALL PRESNM(RES2E,RSNAM(4),FSIZ,X2,Y4-YSHIFT,4,BOLD,RIGHT)

C----    Draw in hydrogen bonds spanning between the two strands
         CALL DASH(4)
         DY = HMNWID / 2.0
         IF(.NOT.BW) CALL COLOUR(HPCOL(5))
         DO 25 J=1,6
            IF(JNAR(J,1).NE.0.AND.JNAR(J,2).NE.0) THEN
               XNAR1=XORIGN-(JNAR(J,1)-1)*RSHIFT
               XNAR2=XORIGN-(JNAR(J,2)-1)*RSHIFT
               CALL FARROW(XNAR1+HSHIFT,Y3-DY,XNAR2+HSHIFT,Y4+DY,
     ~              AMINT,AMAXT,AHDLEN)
               CALL FARROW(XNAR2-HSHIFT,Y4+DY,XNAR1-HSHIFT,Y3-DY,
     ~              AMINT,AMAXT,AHDLEN)
            ENDIF
 25      CONTINUE
         IF(JSNG(1).NE.0.AND.JSNG(2).NE.0) THEN
            XNAR1=XORIGN-(JSNG(1)-1)*RSHIFT
            XNAR2=XORIGN-(JSNG(2)-1)*RSHIFT
            IF(JSNG(3).EQ.1) THEN
               CALL FARROW(XNAR1+HSHIFT,Y3-DY,XNAR2+HSHIFT,Y4+DY,
     ~              AMINT,AMAXT,AHDLEN)
            ELSE
               CALL FARROW(XNAR2-HSHIFT,Y4+DY,XNAR1-HSHIFT,Y3-DY,
     ~              AMINT,AMAXT,AHDLEN)
            ENDIF
         ENDIF
         CALL DASH(0)

C----    Plot the loop residues
         CALL COLOUR(3)
         XLOOP = XORIGN + RSHIFT / 2.0
         YLOOP = Y3

C----    Determine the length and half-length of the loop
         DO 330 J = 1, MXLLEN
            IF (LSEQ(J:J).EQ.' ') THEN
               LPLEN = J - 1
               HALFLEN = INT(J / 2.0)
               GO TO 340
            ENDIF
 330     CONTINUE

C----    Determine whether length is odd or even
 340     CONTINUE
         IF (MOD(LPLEN,2).EQ.0) THEN
             EVEN = .TRUE.
             YDROP = YDROPE
         ELSE
             EVEN = .FALSE.
             YDROP = YDROPO
         ENDIF
         XMAX = XLOOP

C----    Loop over the residues
         DO 30 J = 1, LPLEN

C----       Determine this residue's x-position
            IF (J.LE.HALFLEN) THEN
               XLOOP = XLOOP + RSHIFT
               XMAX = XLOOP
            ELSE IF (J.EQ.HALFLEN+1 .AND. EVEN) THEN
               XLOOP = XLOOP
            ELSE IF (J.NE.HALFLEN+1 .OR. .NOT.EVEN) THEN
               XLOOP = XLOOP - RSHIFT
            ENDIF

C----       Determine the y-position for this residue
            IF (EVEN .AND. J.EQ.HALFLEN) THEN
                PYLOOP = Y3 - YDROP
                YLOOP = Y4
            ELSE IF (EVEN .AND. J.EQ.HALFLEN+1) THEN
                PYLOOP = Y4 + YDROP
            ELSE IF (.NOT.EVEN .AND. J.EQ.HALFLEN - 1) THEN
                PYLOOP = Y3 - YDROP
                YLOOP = Y4
            ELSE IF (.NOT.EVEN .AND. J.EQ.HALFLEN) THEN
                PYLOOP = Y3 - 2 * YDROP
            ELSE IF (.NOT.EVEN .AND. J.EQ.HALFLEN+1) THEN
                PYLOOP = Y4 + YDROP
            ELSE
                PYLOOP = YLOOP 
            ENDIF

C----       Store y-position if residue close to loop start/end
            IF (J.LE.5) YSTORE(J,1) = PYLOOP
            IF (LPLEN - J + 1.LE.5) YSTORE(LPLEN - J + 1,2) = PYLOOP

C----       Print this residue
            CALL TPRINT(LSEQ(J:J),MSIZ,XLOOP,PYLOOP)
 30      CONTINUE

C----    Draw any loop hydrogen bonds
         XLOOP = XORIGN + 3.0 * RSHIFT / 2.0
         CALL DASH(4)
         IF(.NOT.BW) CALL COLOUR(HPCOL(5))
         DO 45 J = 1, 3
            IF (LHB(J,1).NE.0 .OR. LHB(J,2).NE.0) THEN
               XLP1 = XLOOP + (LHB(J,1) - 1) * RSHIFT
               IF (LHB(J,1).EQ.0) XLP1 = XLP1 - RSHIFT / 2.0
               IF (LHB(J,1).EQ.0 .OR. LHB(J,1).GT.5) THEN
                   YLP1 = Y3
               ELSE
                   YLP1 = YSTORE(LHB(J,1),1)
               ENDIF
               XLP2 = XLOOP + (LHB(J,2) - 1) * RSHIFT
               IF (LHB(J,2).EQ.0) XLP2 = XLP2 - RSHIFT / 2.0
               IF (LHB(J,2).EQ.0 .OR. LHB(J,2).GT.5) THEN
                   YLP2 = Y4
               ELSE
                   YLP2 = YSTORE(LHB(J,2),2)
               ENDIF
               IF(LHB(J,3).EQ.1) THEN
                  CALL FARROW(XLP1-2.0,YLP1-4.0,XLP2-2.0,
     ~                 YLP2+4.0,AMINT,AMAXT,AHDLEN)
               ELSE
                  CALL FARROW(XLP2+2.0,YLP2+4.0,XLP1+2.0,
     ~                 YLP1-4.0,AMINT,AMAXT,AHDLEN)
               ENDIF
            ENDIF
 45      CONTINUE
         CALL DASH(0)

C----    Print the hairpin type
         CALL COLOUR(3)
         WRITE(TYPE(1:2),'(I2)') CLASS1
         WRITE(TYPE(4:5),'(I2)') CLASS2
         X1 = XMAX + RSHIFT
         Y1 = (Y3 + Y4) / 2.0
         CALL PSLTXT(X1,Y1,FSIZ,TYPE(1:LENSTR(TYPE)))

C----    Plot motif number, if required
         IF (MOTIF.GT.0) THEN

C----        Format the motif number
             WRITE(NUMBER,10) MOTIF
 10          FORMAT(I5,'.')

C----        Remove leading spaces
             CALL STRIM(NUMBER)

C----        Plot the motif number
             X = XCENTR
             Y = Y4 - 2.5 * YSHIFT
             CALL COLOUR(1)
             CALL BPRINT(NUMBER(1:LENSTR(NUMBER)),FMSIZ,X,Y)
         ENDIF

CPDBSUM RAL 6/4/05 <--

C     check if end of page -close file and start again
CPDBSUM RAL 4/4/05 -->
C        IF(MOD(I,5).EQ.0) THEN
        IF(ONEPAG .AND. MOD(I,5).EQ.0) THEN
            GO TO 70
        ELSE IF (.NOT.ONEPAG) THEN
CPDBSUM RAL 4/4/05 <--
           CALL PSCLOSE
C     otherwise just increment y
        ELSE
           Y0=Y0-YDEL
        ENDIF
        NUM=I
CPDBSUM RAL 29/3/05 -->
C 60   CONTINUE
        IF (ONEPAG) MOTIF = MOTIF + 1
        GO TO 60
CPDBSUM RAL 29/3/05 <--
 70   IF(HP.AND.MOD(NUM,5).NE.0) CALL PSCLOSE

c     write out table of data
CPDBSUM RAL 29/3/05 -->
C      REWIND(UNIT=2)
C      READ(2,*,end=890)
C      Y=640.0
C      DO 80  I=1,NNN
C         READ(2,15,END=90) RES1I,RES1E,RES2I,RES2E,CLASS1,CLASS2,
C     ~        ECLASS,SLEN1,SLEN2
C         READ(2,16,END=90)SSEQ1
C         READ(2,17,END=90)LSEQ
C         READ(2,18,END=90)SSEQ2
C         READ(2,*)
C         READ(2,19,END=90)((JNAR(K,J),J=1,2),K=1,6),
C     ~        ((LHB(K,J),J=1,3),K=1,3)
CC     if at least one hairpin open table file
C         IF(I.EQ.1) THEN
C            HNAMTB=PDBCOD(1:NPS-1)//'_hairpin_tab01.ps'
C            CALL PSOPEN(HNAMTB)
C            CALL TITLE(PDBCOD,'HPIN',NUMTFI,HNAMTB,2)
C         ENDIF
C         IF(Y.LT.50.0) THEN
C            CALL HPINTB(Y)
C            CALL NEWPAG(NUMTFI,HNAMTB,PDBCOD,TTL)
C            Y=640.0
C         ENDIF
C         WRITE(TEXT6,'(A6)') RES1I
C         CALL TPRINT(TEXT6,12.0,75.0,Y)
C         WRITE(TEXT6,'(A6)') RES1E
C         CALL TPRINT(TEXT6,12.0,125.0,Y)
C         WRITE(TEXT6,'(A6)') RES2I
C         CALL TPRINT(TEXT6,12.0,175.0,Y)
C         WRITE(TEXT6,'(A6)') RES2E
C         CALL TPRINT(TEXT6,12.0,225.0,Y)
C         WRITE(TEXT2,'(I2)') SLEN1
C         CALL TPRINT(TEXT2,12.0,275.0,Y)
C         WRITE(TEXT2,'(I2)') SLEN2
C         CALL TPRINT(TEXT2,12.0,325.0,Y)
C         WRITE(TEXT9,'(I2,1X,I2,1X,A3)') CLASS1,CLASS2,
C     ~        ECLASS
C         TEXT9(3:3)=':'
C         CALL TPRINT(TEXT9,12.0,375.0,Y)
C         Y=Y-20.0
C 80   CONTINUE      
C 90   IF(HP) THEN
C         CALL HPINTB(Y)
C         CALL PSCLOSE
C      ENDIF
C      CLOSE(UNIT=2)
C      GO TO 900
C 890  WRITE(*,*) 'No hairpin data'
CPDBSUM RAL 29/3/05 <--

 900  RETURN
      END
c*******************************************************************
c     SUBROUTINE SHEET - plots out data for sheets and strands
c
c*******************************************************************
CPDBSUM RAL 30/3/05 -->
C      SUBROUTINE SHEET(PDBCOD)
      SUBROUTINE SHEET(PDBCOD,CHAIN)
CPDBSUM RAL 30/3/05 <--

      INCLUDE 'psplot.par'
      
      INTEGER NUMSTR,NUMTFI,NUMSFI,NUMMFI,I,J,STNUM,NUMRES,HLXLEN,
     ~     SLEN1,SLEN2,LPLEN,NSTR,NUMIFI,NPS,IOCODE,SHTLT
      REAL Y,AMINT,AMAXT,AHDLEN,Y0,Y1
      CHARACTER SHETYP,
     ~     BARREL,LINE*80,TEXT2*2,TYPES*3,CLASS(3)*12,TPLGY*120,
     ~     STBEG*6,STEND*6,SEQUEN*25,TEXT3*3,STBEG1*6,STEND1*6,
     ~     STBEG2*6,STEND2*6,TTL*4
      CHARACTER*78 IBABNM,IPSNM,ISHTNM,ISTRNM,PDBCOD,SHTFI,STRFI,
     ~     MOTFI
CPDBSUM RAL 30/3/05 -->
      CHARACTER*1 CHAIN
      CHARACTER*(RECLEN) IREC
CPDBSUM RAL 30/3/05 <--
c     CHARACTER SINTFI
      LOGICAL SH,BAB,PSI

      DATA TYPES/'APM'/,TTL/'SHET'/
      DATA AMINT/0.5/,AMAXT/4.0/,AHDLEN/5.0/
      DATA (CLASS(J),J=1,3)
     ~     /'Antiparallel','  Parallel  ','   Mixed    '/
      
      SH=.FALSE.
      BAB=.FALSE.
      NSTR=0
      NUMIFI=1
      NUMTFI=1
      NUMSFI=1
      NUMMFI=1
CPDBSUM RAL 30/3/05 -->
C      NPS=INDEX(PDBCOD,' ')
C      ISHTNM=PDBCOD(1:NPS-1)//'.sht'
      NPS = 5
      ISHTNM='pdb'//PDBCOD(1:NPS-1)//'.sht'
CPDBSUM RAL 30/3/05 <--
      OPEN(UNIT=2,FILE=ISHTNM,STATUS='OLD',ERR=900)
C      READ(2,*)

C 10   READ(2,'(A)') LINE
C      IF(LINE(26:36).NE.'beta sheets') GO TO 10
      DO 20 I=1,NNN
         READ(2,'(A80)',END=30) LINE
         write(*,*) line
         IF(I.EQ.1) THEN
CPDBSUM RAL 30/3/05 -->
C            NPS=INDEX(PDBCOD,' ')
CPDBSUM RAL 30/3/05 <--
            SHTFI=PDBCOD(1:NPS-1)//'_sheets_tab01.ps'
            CALL PSOPEN(SHTFI)
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'SHET',NUMTFI,SHTFI,2)
CPDBSUM RAL 16/3/05 <--
            Y=640.0
            SH=.TRUE.
         ENDIF
         IF(LINE(3:3).EQ.' ') GO TO 30 
         READ(LINE,15) SHTLT,NUMSTR,SHETYP,BARREL
c 15      FORMAT(A1,1X,I2,1X,A1,1X,A1)
 15      FORMAT(I3,1X,I2,1X,A1,1X,A1)
         READ(2,'(10X,A120)',END=30) TPLGY
c         READ(2,*)
c         READ(2,*)
c         READ(2,*)
         IF(Y.LT.50.0) THEN
            CALL SHTTB(Y)
            CALL NEWPAG(NUMTFI,SHTFI,PDBCOD,TTL)
            Y=640.0
         ENDIF
         WRITE(TEXT3,'(I3)') SHTLT
         CALL TPRINT(TEXT3,12.0,75.0,Y)
         WRITE(TEXT2,'(I2)') NUMSTR
         CALL TPRINT(TEXT2,12.0,125.0,Y)
         J=INDEX(TYPES,SHETYP)
         CALL TPRINT(CLASS(J),12.0,190.0,Y)
         IF(BARREL.EQ.'Y') CALL TPRINT('Barrel',12.0,265.0,Y)
         CALL TPRINT(TPLGY,7.0,425.0,Y)
         Y=Y-20.0
 20   CONTINUE
 30   IF(SH) THEN
         CALL SHTTB(Y)
         CALL PSCLOSE
         CLOSE(2)
      ENDIF

c     now produce output for strands
      IF(SH) THEN
CPDBSUM RAL 30/3/05 -->
C         ISTRNM=PDBCOD(1:NPS-1)//'.str'
         ISTRNM='pdb'//PDBCOD(1:NPS-1)//'.str'
CPDBSUM RAL 30/3/05 <--
         OPEN(UNIT=7,FILE=ISTRNM,STATUS='OLD')
         TTL='STRN'
         STRFI=PDBCOD(1:NPS-1)//'_strands_tab01.ps'
         CALL PSOPEN(STRFI)
CPDBSUM RAL 16/3/05 -->
C         CALL TITLE(PDBCOD,'STRN',NUMSFI,STRFI,2)
CPDBSUM RAL 16/3/05 <--
         Y=640.0
C         REWIND(2)
C 40      READ(2,'(A)') LINE
C         IF(LINE(26:37).NE.'beta strands') GO TO 40
C         READ(2,*)
         DO 50 I=1,NNN
            READ(7,45,END=60) STNUM,STBEG,STEND,SHTLT,NUMRES,
     ~           SEQUEN
c 45         FORMAT(I4,1X,2(A6,1X),A1,2X,I2,2X,A25)
 45         FORMAT(I4,1X,2(A6,1X),I3,6X,I2,2X,A25)
            IF(STBEG(5:5).EQ.' ') GO TO 60 
            IF(Y.LT.50.0) THEN
               CALL STRTB(Y)
               CALL NEWPAG(NUMSFI,STRFI,PDBCOD,TTL)
               Y=640.0
            ENDIF
            WRITE(TEXT3,'(I3)') STNUM
            CALL TPRINT(TEXT3,12.0,75.0,Y)
            CALL TPRINT(STBEG,12.0,125.0,Y)
            CALL TPRINT(STEND,12.0,175.0,Y)
            WRITE(TEXT3,'(I3)') SHTLT
            CALL TPRINT(TEXT3,12.0,225.0,Y)
            WRITE(TEXT2,'(I2)') NUMRES
            CALL TPRINT(TEXT2,12.0,275.0,Y)
            CALL TPRINT(SEQUEN,12.0,400.0,Y)
            Y=Y-20.0
            NSTR=NSTR+1
 50      CONTINUE
 60      CALL STRTB(Y)
         CALL PSCLOSE
         CLOSE(7)
      ENDIF

c     now produce output for other motifs: beta-alpha-beta units etc

c     first beta-alpha-beta motifs
      IBABNM=PDBCOD(1:NPS-1)//'.bab'
      OPEN(UNIT=8,FILE=IBABNM,STATUS='OLD',IOSTAT=IOCODE)
      IF(IOCODE.EQ.0) THEN
C 70   READ(2,'(A)') LINE
C      IF(LINE(2:16).NE.'BETA ALPHA BETA') GO TO 70
      TTL='MOTF'
      DO 80 I=1,NNN
         READ(8,75,END=90) STBEG1,STEND1,STBEG2,STEND2,SLEN1,SLEN2,
     ~        LPLEN,HLXLEN
 75      FORMAT(A6,1X,A6,3X,A6,1X,A6,4(1X,I3))
         IF(STBEG1(5:5).EQ.' ') GO TO 90
         IF(I.EQ.1) THEN
            MOTFI=PDBCOD(1:NPS-1)//'_motifs_tab01.ps'
            CALL PSOPEN(MOTFI)            
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'MOTF',NUMMFI,MOTFI,2)
CPDBSUM RAL 16/3/05 <--
            Y=640.0
            BAB=.TRUE.
         ENDIF
         IF(Y.LT.50.0) THEN
            CALL BABTB(Y)
            CALL NEWPAG(NUMMFI,MOTFI,PDBCOD,TTL)
            Y=640.0
         ENDIF
         CALL TPRINT(STBEG1,12.0,75.0,Y)
         CALL TPRINT(STEND1,12.0,125.0,Y)
         CALL TPRINT(STBEG2,12.0,175.0,Y)
         CALL TPRINT(STEND2,12.0,225.0,Y)
         WRITE(TEXT3,'(I3)') SLEN1
         CALL TPRINT(TEXT3,12.0,275.0,Y)
         WRITE(TEXT3,'(I3)') SLEN2
         CALL TPRINT(TEXT3,12.0,325.0,Y)
         WRITE(TEXT3,'(I3)') LPLEN
         CALL TPRINT(TEXT3,12.0,375.0,Y)
         WRITE(TEXT3,'(I3)') HLXLEN
         CALL TPRINT(TEXT3,12.0,425.0,Y)
         Y=Y-20.0
 80   CONTINUE
 90   IF(BAB) THEN
         CALL BABTB(Y)
      ENDIF
      CLOSE(8)
      ENDIF
c     now do psi-loops if any
      PSI=.FALSE.
      IPSNM=PDBCOD(1:NPS-1)//'.psi'
      OPEN(UNIT=9,FILE=IPSNM,STATUS='OLD',IOSTAT=IOCODE)
      IF(IOCODE.EQ.0) THEN
C 100  READ(2,'(A)') LINE
C      IF(LINE(2:10).NE.'PSI LOOPS') GO TO 100
         TTL='MOTF'
         DO 110 I=1,NNN
         READ(9,95,END=140) STBEG1,STEND1,STBEG2,STEND2
 95      FORMAT(A6,1X,A6,3X,A6,1X,A6)
         IF(STBEG1(5:5).EQ.' ') GO TO 140
         IF(I.EQ.1) THEN
            PSI=.TRUE.
            IF(.NOT.BAB) THEN
               MOTFI=PDBCOD(1:NPS-1)//'_motifs_tab01.ps'
               CALL PSOPEN(MOTFI)
CPDBSUM RAL 16/3/05 -->
C               CALL TITLE(PDBCOD,'MOTF',NUMMFI,MOTFI,2)
CPDBSUM RAL 16/3/05 <--
               Y=640.0
               PSI=.TRUE.
               Y0=Y+20.0
               Y1=Y-60.0
            ELSE
               Y0=Y-50.0
               Y1=Y-130.0
            ENDIF
         ENDIF
         IF(Y1.LT.50.0) THEN
            CALL PSITB(Y1,Y0)
            CALL NEWPAG(NUMMFI,MOTFI,PDBCOD,TTL)
            Y=640.0
            Y0=Y+20.0
            Y1=Y-60.0
         ENDIF
         CALL TPRINT(STBEG1,12.0,175.0,Y1)
         CALL TPRINT(STEND1,12.0,225.0,Y1)
         CALL TPRINT(STBEG2,12.0,275.0,Y1)
         CALL TPRINT(STEND2,12.0,325.0,Y1)
         Y1=Y1-20.0
 110  CONTINUE
 140  IF(PSI) THEN
         CALL PSITB(Y1,Y0)
      ENDIF
      CLOSE(UNIT=9)
      ENDIF
      IF(BAB.OR.PSI) CALL PSCLOSE
 900  RETURN
      END
C*******************************************************************
c     SUBROUTINE COORD - CONVERTS PHI,PSI VALUES TO POSTSCRIPT COORDS
C
C******************************************************************
      SUBROUTINE COORD(X,Y,RSIZ,X0,Y0,XPLOT,YPLOT)

      REAL X,Y,XPLOT,YPLOT,RSIZ,X0,Y0

      XPLOT=X+180.0
      XPLOT=XPLOT*RSIZ/360.0
      XPLOT=XPLOT+X0

      YPLOT=Y+180.0
      YPLOT=YPLOT*RSIZ/360.0
      YPLOT=YPLOT+Y0


      RETURN
      END
C********************************************************************
C     SUBROUTINE PLTNUM - Plots beginning and end residues on ramach 
C                         and turntype
C
C******************************************************************
      SUBROUTINE PLTNUM(R1,R2,TURNTY,X0,Y0,RSIZ)
      
      REAL X0,Y0,RSIZ,XMAX,XMID,YMAX,FSIZ,YDEL,YDEL2
      CHARACTER R1*6,R2*6,RTITLE*13
      CHARACTER TURNTY*(*)
      
      DATA FSIZ/12.0/
      
      RTITLE(1:6)=R1
      RTITLE(7:7)='-'
      RTITLE(8:13)=R2
      
      XMAX=X0+RSIZ
      XMID=(X0+XMAX)/2.0
      YMAX=Y0+RSIZ
      
      YDEL=(FSIZ+1.0)/2.0
      YDEL2=FSIZ+YDEL
      
      CALL TPRINT(TURNTY,FSIZ,XMID,YMAX+YDEL)
      CALL TPRINT(RTITLE,FSIZ,XMID,YMAX+YDEL2)
      
      RETURN
      END
CPDBSUM RAL 5/4/05 -->
C**************************************************************************
C
C  SUBROUTINE PLTNM2  -  Plot residue range above Ramachandran
C
C----------------------------------------------------------------------+---

      SUBROUTINE PLTNM2(RSNAM1,RSNUM1,RSNAM2,RSNUM2,TURNTY,X0,Y0,RSIZ,
     -    MOTIF)
      
      CHARACTER*3   RSNAM1, RSNAM2
      CHARACTER*5   RSNUM1, RSNUM2
      CHARACTER*6   NUMBER
      CHARACTER*17  RTITLE
      CHARACTER     TURNTY*(*)

      INTEGER       LENSTR, MOTIF
      
      REAL FSIZ, RSIZ, X, XDEL, XMAX, XMID, X0, Y, YDEL, YDEL2, Y0, YMAX

      DATA FSIZ/11.0/

C---- Initialise variables
      XMAX = X0 + RSIZ 
      XMID = (X0 + XMAX) / 2.0
      YMAX = Y0 + RSIZ
      XDEL = FSIZ / 2.0
      YDEL = FSIZ
      YDEL2 = FSIZ + YDEL

C---- Plot the motif number if required
      IF (MOTIF.GT.0) THEN

C----     Format the motif number
          WRITE(NUMBER,10) MOTIF
 10       FORMAT(I5,'.')

C----     Remove leading spaces
          CALL STRIM(NUMBER)

C----     Plot the motif number
          X = X0 - XDEL
          Y = YMAX + YDEL2
          CALL COLOUR(1)
          CALL PSRTXT(X,Y,FSIZ,NUMBER(1:LENSTR(NUMBER)))
      ENDIF

C---- Fix case for residue names
      CALL FRMRES(RSNAM1)
      CALL FRMRES(RSNAM2)

C---- Remove leading spaces from residue number
      CALL STRIM(RSNUM1)
      CALL STRIM(RSNUM2)

C---- Form name of first residue
      RTITLE = RSNAM1(1:LENSTR(RSNAM1)) // RSNUM1

C---- Print range-start using right-justified
      X = XMID - XDEL
      Y = YMAX + YDEL2
      CALL COLOUR(3)
      CALL PSRTXT(X,Y,FSIZ,RTITLE(1:LENSTR(RTITLE)))

C---- Plot hyphen
      X = XMID
      CALL COLOUR(1)
      CALL BPRINT('-',FSIZ,X,Y)

C---- Form name of last residue
      RTITLE = RSNAM2(1:LENSTR(RSNAM2)) // RSNUM2

C---- Print range-end using left-justified
      X = XMID + XDEL
      CALL COLOUR(5)
      CALL PSLTXT(X,Y,FSIZ,RTITLE(1:LENSTR(RTITLE)))
      
C---- Plot turn type
      X = XMID
      Y = YMAX + YDEL
      CALL COLOUR(9)
      CALL BPRINT(TURNTY,FSIZ,X,Y)
      
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PRESNM  -  Plot given residue number
C
C----------------------------------------------------------------------+---

      SUBROUTINE PRESNM(SEQNO,RSNAME,FNTSIZ,X,Y,ICOL,BOLD,JUST)
      
      CHARACTER*3   RESNAM, RSNAME
      CHARACTER*5   RESNUM
      CHARACTER*6   SEQNO
      CHARACTER*8   TITLE

      INTEGER       ICOL, JUST, LENSTR
      LOGICAL       BOLD

      REAL          FNTSIZ, X, Y

C---- Fix case for residue name
      RESNAM = RSNAME
      CALL FRMRES(RESNAM)

C---- Remove leading spaces from residue number
      RESNUM = SEQNO(2:)
      CALL STRIM(RESNUM)

C---- Form name of residue
      TITLE = RESNAM(1:LENSTR(RESNAM)) // RESNUM

C---- Print residue name
      CALL COLOUR(ICOL)
      IF (JUST.EQ.-1) THEN
          CALL PSLTXT(X,Y,FNTSIZ,TITLE(1:LENSTR(TITLE)))
      ELSE IF (JUST.EQ.1) THEN
          CALL PSRTXT(X,Y,FNTSIZ,TITLE(1:LENSTR(TITLE)))
      ELSE IF (BOLD) THEN
          CALL BPRINT(TITLE(1:LENSTR(TITLE)),FNTSIZ,X,Y)
      ELSE
          CALL TPRINT(TITLE(1:LENSTR(TITLE)),FNTSIZ,X,Y)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE FRMRES  -  Format residue name
C
C----------------------------------------------------------------------+---

      SUBROUTINE FRMRES(RESNAM)

      INTEGER       IPOS

      CHARACTER*1   CH
      CHARACTER*3   RESNAM

C---- Loop over all the characters to be converted to lower case
      DO 600, IPOS = 2, 3

C----     Retrieve this character
          CH = RESNAM(IPOS:IPOS)

C----     Check if this is a letter
          IF (LGE(CH,'A') .AND. LLE(CH,'Z')) THEN

C----         Convert to lower-case
              CH = CHAR(ICHAR('a') + ICHAR(CH) - ICHAR('A'))

C----         Put back into residue name
              RESNAM(IPOS:IPOS) = CH
          ENDIF
 600  CONTINUE

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE STRIM  -  Remove leading space from the given string
C
C----------------------------------------------------------------------+---

      SUBROUTINE STRIM(STRING)

      INTEGER       ILEN, IPOS, NOBLNK

      CHARACTER*1   CH
      CHARACTER*(*) STRING

C---- Initialise variables
      NOBLNK = 0
      ILEN = LEN(STRING)

C---- Loop over all the leading characters
      DO 600, IPOS = 1, ILEN

C----     Retrieve this character
          CH = STRING(IPOS:IPOS)

C----     If first non-blank character, store its position
          IF (CH.NE.' ' .AND. NOBLNK.EQ.0) THEN

C----         Store position
              NOBLNK = IPOS
          ENDIF
 600  CONTINUE

C---- Shift string to left-justify
      IF (NOBLNK.GT.0) THEN
          STRING = STRING(NOBLNK:)
      ENDIF

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE BLABEL  -  Print the label identifying the bulge residue
C
C----------------------------------------------------------------------+---

      SUBROUTINE BLABEL(IRES,FNTSIZ,X,Y,ICOL)

      CHARACTER*1   BPOS(5)
      INTEGER       ICOL, ICOLRS(5), IRES
      REAL          FNTSIZ, X, Y

      DATA BPOS   / 'X', '1', '2', '3', '4' /
      DATA ICOLRS /   3,   4,   5,   9,   9 /

C---- Call appropriate colour
      ICOL = ICOLRS(IRES)
      CALL COLOUR(ICOL)

C---- Print the label
      CALL BPRINT(BPOS(IRES),FNTSIZ,X,Y)

      RETURN
      END

C--------------------------------------------------------------------------
CPDBSUM RAL 5/4/05 <--
C******************************************************************
C     SUBROUTINE TITLE - Prints out title at top of page
C
C******************************************************************
      SUBROUTINE TITLE(PDBCOD,TYPE,PAGNUM,FNAME,NCOL)

      INCLUDE 'psplot.par'

      INTEGER PAGNUM,ILEN,NCOL
      REAL TFSIZ,FSIZ
      CHARACTER PDBCOD*78,TITL*40,TYPE*4,PAGE*7
      CHARACTER*(*) FNAME
      
      DATA TFSIZ/20.0/,FSIZ/9.0/
C      TITL='                          Page  '
      PAGE='Page   '

      IF(NCOL.EQ.1) THEN
         IF(.NOT.BW) CALL COLOUR(15)
      ENDIF
CPDBSUM RAL 16/3/05 -->
C      ILEN=INDEX(PDBCOD,' ')
      ILEN = 5
CPDBSUM RAL 16/3/05 <--
      TITL='                                        '

      IF(TYPE.EQ.'BETA') THEN
         TITL(1:12)= 'Beta turns: '
         TITL(13:13+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'GAMM') THEN
         TITL(1:13)= 'Gamma turns: '
         TITL(14:14+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'BULG') THEN
         TITL(1:8)= 'Bulges: '
         TITL(9:9+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'HPIN') THEN
         TITL(1:10)='Hairpins: '
         TITL(11:11+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'SHET') THEN
         TITL(1:8)='Sheets: '
         TITL(9:9+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'STRN') THEN
         TITL(1:9)='Strands: '
         TITL(10:10+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'SINT') THEN
         TITL(1:21)='Strand Interactions: '
         TITL(22:22+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'MOTF') THEN
         TITL(1:8)='Motifs: '
         TITL(9:9+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'HELX') THEN
         TITL(1:9)='Helices: '
         TITL(10:10+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'HINT') THEN
         TITL(1:20)='Helix Interactions: '
         TITL(21:21+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'HGEM') THEN
         TITL(1:16)='Helix Geometry: '
         TITL(17:17+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'HWHL') THEN
         TITL(1:14)='Helical Nets: '
         TITL(15:15+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'DSLF') THEN
         TITL(1:20)='Disulphide Bridges: '
         TITL(21:21+ILEN-1)=PDBCOD(1:ILEN-1)
      ELSE IF(TYPE.EQ.'SUMM') THEN
         TITL(1:9)='Summary: '
         TITL(10:10+ILEN-1)=PDBCOD(1:ILEN-1)
      ENDIF

      WRITE(PAGE(6:7),'(I2)') PAGNUM
         
      ILEN=INDEX(TITL,'   ')      
      IF(ILEN.NE.0) THEN
         CALL TPRINT(TITL(1:ILEN-1),TFSIZ,300.0,760.0)
      ELSE
         CALL TPRINT(TITL,TFSIZ,300.0,760.0)
      ENDIF
      ILEN=INDEX(FNAME,' ')
      IF(ILEN.NE.0) THEN
         CALL TPRINT(FNAME(1:ILEN-1),FSIZ,150.0,760.0)
      ELSE
         CALL TPRINT(FNAME,FSIZ,150.0,760.0)
      ENDIF
      CALL TPRINT(PAGE,TFSIZ,500.0,760.0)
      
      RETURN
      END
C**********************************************************************
      SUBROUTINE DOTTED(X1,Y1,X2,Y2)

      REAL X1,Y1,X2,Y2

      WRITE(1,5) X1,Y1,X2,Y2
 5    FORMAT(4(F5.1,1X),'Dline')

      RETURN
      END
c************************************************************************
C     SUBROUTINE CONVERT - converts a yes/no answer into logical true 
c                          and false
c*************************************************************************
      LOGICAL FUNCTION LOGICS(YESNO)
      
      CHARACTER YESNO

      IF(YESNO.EQ.'Y') THEN
         LOGICS=.TRUE.
      ELSE
         LOGICS=.FALSE.
      ENDIF
      
      END
c******************************************************************
      SUBROUTINE BHBOND(HBP,JSEQNOP,XCORDP,YCORDP,XCOORD,YCOORD,AP)
c     plots hydrogen bonds on schematic diagrams of bulges
      INCLUDE 'psplot.par'

      INTEGER I,J,K
      REAL XCOORD(8),YCOORD(8),XCORDP(5),YCORDP(5),YHB,XPART,
     ~     YPART,AMINT,AMAXT,AHDLEN,XSIZ
      CHARACTER HBP(8,4)*6,JSEQNOP(5)*6,AP

CPDBSUM RAL 6/4/05 -->
C      DATA AMINT/0.5/,AMAXT/4.0/,AHDLEN/5.0/,XSIZ/15.0/,
      DATA AMINT/0.5/,AMAXT/4.0/,AHDLEN/5.0/,XSIZ/30.0/,
CPDBSUM RAL 6/4/05 <--
     ~     YHB/2.0/

      IF(.NOT.BW) CALL COLOUR(BULCOL(6))
      DO 110 K=1,8
         DO 109 J=1,4
            IF(HBP(K,J).NE.'000000') THEN
               DO 108 I=1,5
                  IF(HBP(K,J).EQ.JSEQNOP(I)) THEN
c     there is a hydrogen bond to be plotted - find coordinates
                     XPART=XCORDP(I)
                     YPART=YCORDP(I)
                     IF(XCOORD(K).NE.0.0.AND.XPART.NE.0.0) THEN
                        IF(AP.EQ.'A') THEN
                           IF(MOD(J,2).EQ.1) THEN
                              CALL FARROW(XCOORD(K)-XSIZ/2.0,
     ~                             YCOORD(K)+YHB,XPART+XSIZ/2.0,
     ~                             YPART+YHB,AMINT,AMAXT,AHDLEN)
                           ELSE
                              CALL FARROW(XPART+XSIZ/2.0,YPART-YHB,
     ~                             XCOORD(K)-XSIZ/2.0,YCOORD(K)-YHB,
     ~                             AMINT,AMAXT,AHDLEN)
                           ENDIF
                        ELSE
                           IF(MOD(J,2).EQ.1) THEN
                              CALL FARROW(XCOORD(K)-XSIZ/2.0,
     ~                             YCOORD(K)-YHB,XPART+XSIZ/2.0,
     ~                             YPART+YHB,AMINT,AMAXT,AHDLEN)
                           ELSE
                              CALL FARROW(XPART+XSIZ/2.0,YPART-YHB,
     ~                             XCOORD(K)-XSIZ/2.0,YCOORD(K)+YHB,
     ~                             AMINT,AMAXT,AHDLEN)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
 108           CONTINUE
            ENDIF
 109     CONTINUE
 110  CONTINUE
      IF(.NOT.BW) CALL COLOUR(BULCOL(5))
      
      RETURN
      END
C**********************************************************************
C      SUBROUTINE BTRNTB - DRAWS OUT A TABLE FORMAT FOR BTURN DATA
C**********************************************************************
      SUBROUTINE BTRNTB(Y)

      REAL Y
      
      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(130.0,Y,130.0,700.0)
      CALL PLINE(180.0,Y,180.0,700.0)
      CALL PLINE(220.0,Y,220.0,700.0)
      CALL PLINE(260.0,Y,260.0,675.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(340.0,Y,340.0,675.0)
      CALL PLINE(380.0,Y,380.0,700.0)
      CALL PLINE(420.0,Y,420.0,700.0)
      CALL PLINE(460.0,Y,460.0,700.0)
      CALL PLINE(500.0,Y,500.0,700.0)
      CALL PLINE(550.0,Y,550.0,700.0)
      CALL PLINE(50.0,700.0,550.0,700.0)
      CALL PLINE(50.0,650.0,550.0,650.0)
      CALL PLINE(50.0,Y,550.0,Y)
      CALL PLINE(220.0,675.0,380.0,675.0)
      CALL TPRINT('Residue',12.0,90.0,685.0)
      CALL TPRINT('numbers',12.0,90.0,667.0)
      CALL TPRINT('Sequence',11.0,155.0,675.0)
      CALL TPRINT('Turn',12.0,200.0,685.0)
      CALL TPRINT('Type',12.0,200.0,667.0)
      CALL SYMB('f',12.0,240.0,667.0)
      CALL SYMB('y',12.0,280.0,667.0)
      CALL TPRINT('(i+1)',12.0,260.0,685.0)
      CALL SYMB('f',12.0,320.0,667.0)
      CALL SYMB('y',12.0,360.0,667.0)
      CALL TPRINT('(i+2)',12.0,340.0,685.0)
      CALL SYMB('f',12.0,395.0,685.0)
      CALL SYMB('y',12.0,405.0,685.0)
      CALL TPRINT('Region',12.0,400.0,667.0)
      CALL TPRINT('Chi1',12.0,440.0,685.0)
      CALL TPRINT('(i+1)',12.0,440.0,667.0)
      CALL TPRINT('Chi1',12.0,480.0,685.0)
      CALL TPRINT('(i+2)',12.0,480.0,667.0)
      CALL TPRINT('i to i+3',12.0,525.0,685.0)
      CALL TPRINT('Distance',12.0,525.0,667.0)
      RETURN
      END
c**********************************************************************
C      SUBROUTINE GTRNTB - DRAWS OUT A TABLE FORMAT FOR GTURN DATA
C**********************************************************************
      SUBROUTINE GTRNTB(Y)
      
      REAL Y 

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(350.0,Y,350.0,675.0)
      CALL PLINE(400.0,Y,400.0,700.0)
      CALL PLINE(450.0,Y,450.0,700.0)
      CALL PLINE(50.0,700.0,450.0,700.0)
      CALL PLINE(50.0,650.0,450.0,650.0)
      CALL PLINE(50.0,Y,450.0,Y)
      CALL PLINE(300.0,675.0,400.0,675.0)
      CALL TPRINT('Residue numbers',12.0,100.0,675.0)
      CALL TPRINT('Sequence',12.0,175.0,675.0)
      CALL TPRINT('Turn Type',12.0,250.0,675.0)
      CALL SYMB('f',12.0,325.0,667.0)
      CALL SYMB('y',12.0,375.0,667.0)
      CALL TPRINT('(i+1)',12.0,350.0,685.0)
      CALL TPRINT('i to i+2',12.0,425.0,685.0)
      CALL TPRINT('Distance',12.0,425.0,667.0)
      RETURN
      END
c**********************************************************************
C      SUBROUTINE BULGTB - DRAWS OUT A TABLE FORMAT FOR BULGE DATA
C**********************************************************************
      SUBROUTINE BULGTB(Y)

      REAL Y
      
      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(260.0,Y,260.0,700.0)
      CALL PLINE(310.0,Y,310.0,700.0)
      CALL PLINE(350.0,Y,350.0,675.0)
      CALL PLINE(390.0,Y,390.0,700.0)
      CALL PLINE(430.0,Y,430.0,675.0)
      CALL PLINE(470.0,Y,470.0,700.0)
      CALL PLINE(510.0,Y,510.0,675.0)
      CALL PLINE(550.0,Y,550.0,700.0)

      CALL PLINE(50.0,700.0,550.0,700.0)
      CALL PLINE(50.0,650.0,550.0,650.0)
      CALL PLINE(50.0,Y,550.0,Y)
      CALL PLINE(310.0,675.0,550.0,675.0)
      CALL TPRINT('Residue numbers',12.0,125.0,685.0)
      CALL TPRINT('X',12.0,75.0,667.0)
      CALL TPRINT('1',12.0,125.0,667.0)
      CALL TPRINT('2',12.0,175.0,667.0)
      CALL TPRINT('Sequence',12.0,230.0,685.0)
      CALL TPRINT('X',12.0,210.0,667.0)
      CALL TPRINT('1',12.0,230.0,667.0)
      CALL TPRINT('2',12.0,250.0,667.0)
      CALL TPRINT('Bulge',12.0,285.0,685.0)
      CALL TPRINT('Type',12.0,285.0,667.0)
      CALL SYMB('f',12.0,330.0,667.0)
      CALL SYMB('y',12.0,370.0,667.0)
      CALL TPRINT('X',12.0,350.0,685.0)
      CALL SYMB('f',12.0,410.0,667.0)
      CALL SYMB('y',12.0,450.0,667.0)
      CALL TPRINT('1',12.0,430.0,685.0)
      CALL SYMB('f',12.0,490.0,667.0)
      CALL SYMB('y',12.0,530.0,667.0)
      CALL TPRINT('2',12.0,510.0,685.0)
      RETURN
      END
c************************************************************************
C      SUBROUTINE HPINTB - DRAWS OUT A TABLE FORMAT FOR HAIRPIN DATA
C**********************************************************************
      SUBROUTINE HPINTB(Y)

      REAL Y
      
      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,675.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,675.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,675.0)
      CALL PLINE(350.0,Y,350.0,700.0)
      CALL PLINE(400.0,Y,400.0,700.0)

      CALL PLINE(50.0,700.0,400.0,700.0)
      CALL PLINE(50.0,650.0,400.0,650.0)
      CALL PLINE(50.0,675.0,350.0,675.0)
      CALL PLINE(50.0,Y,400.0,Y)
      CALL TPRINT('Strand 1',12.0,100.0,685.0)
      CALL TPRINT('Start',12.0,75.0,667.0)
      CALL TPRINT('End',12.0,125.0,667.0)
      CALL TPRINT('Strand 2',12.0,200.0,685.0)
      CALL TPRINT('Start',12.0,175.0,667.0)
      CALL TPRINT('End',12.0,225.0,667.0)
      CALL TPRINT('Number of residues',12.0,300.0,685.0)
      CALL TPRINT('Strand 1',12.0,275.0,667.0)
      CALL TPRINT('Strand 2',12.0,325.0,667.0)
      CALL TPRINT('Hairpin',12.0,375.0,685.0)
      CALL TPRINT('Class',12.0,375.0,667.0)
      RETURN
      END
c************************************************************************
C      SUBROUTINE SHTTB - DRAWS OUT A TABLE FORMAT FOR SHEET DATA
C
C**********************************************************************
      SUBROUTINE SHTTB(Y)

      REAL Y
      
      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(550.0,Y,550.0,700.0)

      CALL PLINE(50.0,700.0,550.0,700.0)
      CALL PLINE(50.0,650.0,550.0,650.0)
      CALL PLINE(50.0,Y,550.0,Y)
      CALL TPRINT('Sheet',12.0,75.0,685.0)
      CALL TPRINT('Letter',12.0,75.0,667.0)
      CALL TPRINT('No. of',12.0,125.0,685.0)
      CALL TPRINT('Strands',12.0,125.0,667.0)
      CALL TPRINT('Sheet type',12.0,225.0,675.0)
      CALL TPRINT('Topology',12.0,425.0,675.0)
      RETURN
      END
c************************************************************************
C      SUBROUTINE STRTB - DRAWS OUT A TABLE FORMAT FOR STRAND DATA
C
C**********************************************************************
      SUBROUTINE STRTB(Y)
      
      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(500.0,Y,500.0,700.0)

      CALL PLINE(50.0,700.0,500.0,700.0)
      CALL PLINE(50.0,650.0,500.0,650.0)
      CALL PLINE(50.0,Y,500.0,Y)
      CALL TPRINT('Strand',12.0,75.0,685.0)
      CALL TPRINT('Number',12.0,75.0,667.0)
      CALL TPRINT('Start',12.0,125.0,685.0)
      CALL TPRINT('Residue',12.0,125.0,667.0)
      CALL TPRINT('End',12.0,175.0,685.0)
      CALL TPRINT('Residue',12.0,175.0,667.0)
      CALL TPRINT('Sheet',12.0,225.0,685.0)
      CALL TPRINT('Label',12.0,225.0,667.0)
      CALL TPRINT('No. of',12.0,275.0,685.0)
      CALL TPRINT('Residues',12.0,275.0,667.0)
      CALL TPRINT('Sequence',12.0,400.0,675.0)
      RETURN
      END
C************************************************************************
C      SUBROUTINE SINTTB - DRAWS OUT A TABLE FORMAT FOR STRAND INTERACTION DATA
C
C**********************************************************************
      SUBROUTINE SINTTB(Y)
      
      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(220.0,Y,220.0,675.0)
      CALL PLINE(290.0,Y,290.0,675.0)
      CALL PLINE(360.0,Y,360.0,675.0)
      CALL PLINE(430.0,Y,430.0,700.0)

      CALL PLINE(50.0,700.0,430.0,700.0)
      CALL PLINE(50.0,650.0,430.0,650.0)
      CALL PLINE(150.0,675.0,430.0,675.0)
      CALL PLINE(50.0,Y,430.0,Y)
      CALL TPRINT('Strand',12.0,75.0,685.0)
      CALL TPRINT('Number',12.0,75.0,667.0)
      CALL TPRINT('No. of',12.0,125.0,685.0)
      CALL TPRINT('Partners',12.0,125.0,667.0)
      CALL TPRINT('Partners',12.0,290.0,685.0)
      RETURN
      END
c************************************************************************
C     SUBROUTINE BABTB - DRAWS OUT A TABLE FORMAT FOR BETA-ALPHA-BETA DATA
C
C************************************************************************
      SUBROUTINE BABTB(Y)

      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      
      CALL TPRINT('Beta-Alpha-Beta Units',14.0,250.0,720.0)

      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(350.0,Y,350.0,700.0)
      CALL PLINE(400.0,Y,400.0,700.0)
      CALL PLINE(450.0,Y,450.0,700.0)

      CALL PLINE(50.0,700.0,450.0,700.0)
      CALL PLINE(50.0,650.0,450.0,650.0)
      CALL PLINE(50.0,Y,450.0,Y)
      CALL TPRINT('Strand 1',12.0,75.0,685.0)
      CALL TPRINT('Start',12.0,75.0,667.0)
      CALL TPRINT('Strand 1',12.0,125.0,685.0)
      CALL TPRINT('End',12.0,125.0,667.0)
      CALL TPRINT('Strand 2',12.0,175.0,685.0)
      CALL TPRINT('Start',12.0,175.0,667.0)
      CALL TPRINT('Strand 2',12.0,225.0,685.0)
      CALL TPRINT('End',12.0,225.0,667.0)
      CALL TPRINT('Strand 1',12.0,275.0,685.0)
      CALL TPRINT('Length',12.0,275.0,667.0)
      CALL TPRINT('Strand 2',12.0,325.0,685.0)
      CALL TPRINT('Length',12.0,325.0,667.0)
      CALL TPRINT('Loop',12.0,375.0,685.0)
      CALL TPRINT('Length',12.0,375.0,667.0)
      CALL TPRINT('Helix',12.0,425.0,685.0)
      CALL TPRINT('Length',12.0,425.0,667.0)
      RETURN
      END
c****************************************************************
C     SUBROUTINE PSITB - DRAWS OUT A TABLE FORMAT FOR PSI-LOOPS 
C
C************************************************************************
      SUBROUTINE PSITB(Y1,Y0)

      REAL Y1,Y0

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      
      CALL TPRINT('Psi-Loops ',14.0,250.0,Y0+10.0)

      CALL PLINE(150.0,Y1,150.0,Y0-10.0)
      CALL PLINE(200.0,Y1,200.0,Y0-10.0)
      CALL PLINE(250.0,Y1,250.0,Y0-10.0)
      CALL PLINE(300.0,Y1,300.0,Y0-10.0)
      CALL PLINE(350.0,Y1,350.0,Y0-10.0)

      CALL PLINE(150.0,Y0-10.0,350.0,Y0-10.0)
      CALL PLINE(150.0,Y0-60.0,350.0,Y0-60.0)
      CALL PLINE(150.0,Y1,350.0,Y1)
      CALL TPRINT('Strand 1',12.0,175.0,Y0-25.0)
      CALL TPRINT('Start',12.0,175.0,Y0-45.0)
      CALL TPRINT('Strand 1',12.0,225.0,Y0-25.0)
      CALL TPRINT('End',12.0,225.0,Y0-45.0)
      CALL TPRINT('Strand 2',12.0,275.0,Y0-25.0)
      CALL TPRINT('Start',12.0,275.0,Y0-45.0)
      CALL TPRINT('Strand 2',12.0,325.0,Y0-25.0)
      CALL TPRINT('End',12.0,325.0,Y0-45.0)
      RETURN
      END
c****************************************************************
c     SUBROUTINE HELIX - Draws schematic diagrams for helices
c
c*****************************************************************
CPDBSUM RAL 29/3/05 -->
C      SUBROUTINE HELIX(PDBCOD)
      SUBROUTINE HELIX(PDBCOD,CHAIN,ONEPAG,MSTART)
CPDBSUM RAL 29/3/05 <--

      INCLUDE 'psplot.par'
      
      INTEGER NUMTFI,HLXNUM,NUMRES,I,NUMGFI,HLX1,HLX2,ILEN,NINT,
     ~     NINT1,NINT2,NUMIFI,NUMWFI,NPS,NUM
      REAL Y,LEN,RISE,RPT,PITCH,DEV,DIST,ANG,X
CPDBSUM RAL 29/3/05 -->
      CHARACTER*1 CHAIN
      CHARACTER*3 RSNAMB, RSNAME
      CHARACTER*(RECLEN) IREC
CPDBSUM RAL 29/3/05 <--
      CHARACTER DASH,HLXBEG*6,HLXEND*6,HTYPE,
     ~     TEXT3*3,TEXT2*2,LINE*125,TEXT6*6,TEXT5*5,
     ~     TEXT4*4,HSEQ*55,TTL*4
      CHARACTER*78 PDBCOD,HLXFI,HINTFI,HGEFI,HDIAGF,IHXNM
      LOGICAL HL
CPDBSUM RAL 29/3/05 -->
      LOGICAL EOF, ONEPAG
      INTEGER IPOS, MOTIF, MSTART
CPDBSUM RAL 29/3/05 <--

      DATA DASH/'-'/
      
CPDBSUM RAL 29/3/05 -->
C      HL=.FALSE.
      HL=.TRUE.
CPDBSUM RAL 29/3/05 <--
      NUMTFI=1
      NUMGFI=1
      NUMIFI=1
      NUMWFI=1
      TTL='HELX'

CPDBSUM RAL 29/3/05 -->
C      NPS=INDEX(PDBCOD,' ')
C      IHXNM=PDBCOD(1:NPS-1)//'.hlx'
C      OPEN(UNIT=2,FILE=IHXNM,STATUS='OLD',ERR=900)
C      READ(2,*,END=890)
C
c     first do standard helix data
C 10   READ(2,'(A)',END=890) LINE
C      IF(LINE(2:11).NE.'Helix data') GO TO 10
C      
C      DO 20 I=1,100
C         READ(2,'(A125)',END=30) LINE
C         IF(LINE(3:3).EQ.' ') GO TO 30
C         HSEQ=
C     ~    '                                                       '
c     modification for all helix info in 1 table 11/11/95
c         READ(LINE,15) HLXNUM,HLXBEG,HLXEND,HTYPE,NUMRES,len,rise,rpt,
c     ~        pitch,dev,HSEQ
c 15      FORMAT(I3,3X,A6,1X,A6,1X,A1,1X,I3,10X,4(1x,f6.2),1x,
c     ~        f5.1,1X,A55)
c     remove this 11/3/96
C         READ(LINE,15) HLXNUM,HLXBEG,HLXEND,HTYPE,NUMRES,HSEQ
C 15      FORMAT(I3,3X,A6,1X,A6,1X,A1,1X,I3,45X,A55)
C         IF(I.EQ.1) THEN
C            NPS=INDEX(PDBCOD,' ')
C            HLXFI=PDBCOD(1:NPS-1)//'_helix_tab01.ps'
C            CALL PSOPEN(HLXFI)
C            CALL TITLE(PDBCOD,'HELX',NUMTFI,HLXFI,2)
C            Y=640.0
C            HL=.TRUE.
C         ENDIF
C         IF(Y.LT.50.0) THEN
C            CALL HLXTB(Y)
C            CALL NEWPAG(NUMTFI,HLXFI,PDBCOD,TTL)
C            Y=640.0
C         ENDIF
C         IF(HLXBEG(6:6).EQ.'-') HLXBEG(6:6)=' '
C         IF(HLXEND(6:6).EQ.'-') HLXEND(6:6)=' '
C         WRITE(TEXT3,'(I3)') HLXNUM
CC     11/11/95 all helix info in 1 table
CC         CALL TPRINT(TEXT3,12.0,70.0,Y)
CC         CALL TPRINT(HLXBEG,12.0,120.0,Y)
CC         CALL TPRINT(HLXEND,12.0,170.0,Y)
C         CALL TPRINT(TEXT3,12.0,75.0,Y)
C         CALL TPRINT(HLXBEG,12.0,125.0,Y)
C         CALL TPRINT(HLXEND,12.0,175.0,Y)
C         IF(HTYPE.EQ.'H') THEN
Cc11/11/95            CALL SYMB('a',12.0,210.0,Y)
C            CALL SYMB('a',12.0,225.0,Y)
C         ELSE IF(HTYPE.EQ.'G') THEN
Cc11/11/95            CALL TPRINT('3',12.0,205.0,Y)
Cc11/11/95            CALL TPRINT('10',10.0,214.0,Y-4.0)
C            CALL TPRINT('3',12.0,221.0,Y)
C            CALL TPRINT('10',10.0,229.0,Y-4.0)
C         ENDIF
C         WRITE(TEXT2,'(I2)') NUMRES
Cc11/11/95         CALL TPRINT(TEXT2,12.0,250.0,Y)
C         CALL TPRINT(TEXT2,12.0,275.0,Y)
Cc11/11/95         IF(LEN.NE.0.0) THEN
Cc            WRITE(TEXT6,'(F6.2)') LEN
Cc            CALL TPRINT(TEXT6,12.0,290.0,Y)
Cc         ELSE
Cc            CALL TPRINT(DASH,12.0,290.0,Y)
Cc         ENDIF
Cc         IF(RISE.NE.0.0) THEN
Cc               WRITE(TEXT6,'(F6.2)') RISE
Cc               CALL TPRINT(TEXT6,12.0,330.0,Y)
Cc            ELSE
Cc               CALL TPRINT(DASH,12.0,330.0,Y)
Cc            ENDIF
Cc            IF(RPT.NE.0.0) THEN
Cc               WRITE(TEXT6,'(F6.2)') RPT
Cc               CALL TPRINT(TEXT6,12.0,377.0,Y)
Cc            ELSE
Cc               CALL TPRINT(DASH,12.0,377.0,Y)
Cc            ENDIF
Cc            IF(PITCH.NE.0.0) THEN
Cc               WRITE(TEXT6,'(F6.2)') PITCH
Cc               CALL TPRINT(TEXT6,12.0,423.0,Y)
Cc            ELSE
Cc               CALL TPRINT(DASH,12.0,423.0,Y)
Cc            ENDIF
Cc            IF(DEV.NE.0.0) THEN
Cc               WRITE(TEXT5,'(F5.1)') DEV
Cc               CALL TPRINT(TEXT5,12.0,468.0,Y)
Cc            ELSE
Cc               CALL TPRINT(DASH,12.0,468.0,Y)
Cc            ENDIF
C         ILEN=INDEX(HSEQ,' ')
C         IF(ILEN.NE.0) THEN
Cc11/11/95            CALL TPRINT(HSEQ(1:ILEN-1),7.0,530.0,Y)
C            CALL TPRINT(HSEQ(1:ILEN-1),7.0,425.0,Y)
C         ELSE
Cc11/11/95            CALL TPRINT(HSEQ,7.0,530.0,Y)
C            CALL TPRINT(HSEQ,7.0,425.0,Y)
C         ENDIF
C         Y=Y-20.0
C 20   CONTINUE
C 30   IF(HL) THEN
C         CALL HLXTB(Y)
C         CALL PSCLOSE
C      ENDIF
C
Cc     now do helix geometry
C      IF(HL) THEN
C         HGEFI=PDBCOD(1:NPS-1)//'_helix_geom01.ps'
C         CALL PSOPEN(HGEFI)
C         CALL TITLE(PDBCOD,'HGEM',NUMGFI,HGEFI,2)
C         Y=640.0
C         REWIND(2)
C 40      READ(2,'(A)',END=890) LINE
C         IF(LINE(2:11).NE.'Helix data') GO TO 40
C         
C         TTL='HGEM'
C         DO 50 I=1,100
C            READ(2,'(A80)',END=60) LINE
C            IF(LINE(3:3).EQ.' ') GO TO 60
C            READ(LINE,45) HLXNUM,LEN,RISE,RPT,PITCH,DEV
C 45         FORMAT(I3,32X,4(1X,F6.2),1X,F5.1)
C            IF(Y.LT.50.0) THEN
C               CALL HGETB(Y)
C               CALL NEWPAG(NUMGFI,HGEFI,PDBCOD,TTL)
C               Y=640.0
C            ENDIF
C            WRITE(TEXT3,'(I3)') HLXNUM
C            CALL TPRINT(TEXT3,12.0,75.0,Y)
C            IF(LEN.NE.0.0) THEN
C               WRITE(TEXT6,'(F6.2)') LEN
C               CALL TPRINT(TEXT6,12.0,125.0,Y)
C            ELSE
C               CALL TPRINT(DASH,12.0,125.0,Y)
C            ENDIF
C            IF(RISE.NE.0.0) THEN
C               WRITE(TEXT6,'(F6.2)') RISE
C               CALL TPRINT(TEXT6,12.0,175.0,Y)
C            ELSE
C               CALL TPRINT(DASH,12.0,175.0,Y)
C            ENDIF
C            IF(RPT.NE.0.0) THEN
C               WRITE(TEXT6,'(F6.2)') RPT
C               CALL TPRINT(TEXT6,12.0,225.0,Y)
C            ELSE
C               CALL TPRINT(DASH,12.0,225.0,Y)
C            ENDIF
C            IF(PITCH.NE.0.0) THEN
C               WRITE(TEXT6,'(F6.2)') PITCH
C               CALL TPRINT(TEXT6,12.0,275.0,Y)
C            ELSE
C               CALL TPRINT(DASH,12.0,275.0,Y)
C            ENDIF
C            IF(DEV.NE.0.0) THEN
C               WRITE(TEXT5,'(F5.1)') DEV
C               CALL TPRINT(TEXT5,12.0,325.0,Y)
C            ELSE
C               CALL TPRINT(DASH,12.0,325.0,Y)
C            ENDIF
C            Y=Y-20.0
C 50      CONTINUE
C 60      CALL HGETB(Y)
C         CALL PSCLOSE
C      ENDIF
CPDBSUM RAL 29/3/05 <--

C     now draw out helical wheels and nets
      IF(HL) THEN
CPDBSUM RAL 29/3/05 -->
C         HDIAGF=PDBCOD(1:NPS-1)//'_helix_diag01.ps'
         HDIAGF = 'promotif001.ps'
CPDBSUM RAL 29/3/05 <--
         CALL PSOPEN(HDIAGF)
         IF(.NOT.BW) CALL COLOUR(HXCOL(2))
CPDBSUM RAL 29/3/05 -->
C         CALL TITLE(PDBCOD,'HWHL',NUMWFI,HDIAGF,1)
CPDBSUM RAL 29/3/05 <--
c     comment this out 13.3.95 - no room for it
c         CALL PRECT(500.0,760.0,100.0,50.0)
c         CALL COLOUR(HXCOL(1))
c         CALL TPRINT('Green: Hydrophobic',10.0,500.0,780.0)
c         CALL COLOUR(HXCOL(3))
c         CALL TPRINT('Blue : Polar      ',10.0,489.0,760.0)
c         CALL COLOUR(HXCOL(2))
c         CALL TPRINT('Red  : Charged    ',10.0,492.0,740.0)
         Y=640.0
         X=50.0
CPDBSUM RAL 29/3/05 -->
C         REWIND(2)
C 140      READ(2,'(A)',END=890) LINE
C         IF(LINE(2:11).NE.'Helix data') GO TO 140
CPDBSUM RAL 29/3/05 <--
         TTL='HWHL'
CPDBSUM RAL 29/3/05 -->
C         DO 150 I=1,100
C            READ(2,'(A125)',END=160) LINE
C            IF(LINE(3:3).EQ.' ') GO TO 160
C            READ(LINE,'(I3,3X,A6,1X,A6,3X,I3,45X,A55)') 
C     ~           HLXNUM,HLXBEG,HLXEND,NUMRES,HSEQ
          I = 0
          IF (ONEPAG) THEN
             MOTIF = MSTART
          ELSE
             MOTIF = 0
          ENDIF
 150      CONTINUE
              I = I + 1

C----         Get the next record from the promotif.dat file
              CALL GETREC(FILDAT,CHAIN,IREC,EOF)
              IF (EOF) GO TO 160

C----         Extract the fields from the line just read in
              HLXBEG(1:1) = IREC(1:1)
              HLXEND(1:1) = IREC(1:1)
              READ(IREC,25) HLXNUM, HLXBEG(2:6), RSNAMB, HLXEND(2:6),
     -            RSNAME, HTYPE, NUMRES
 25           FORMAT(3X,I3,2X,2(A5,2X,A3,2X),A1,2X,I3)
              IPOS = INDEX(IREC(91:),'::')
              IF (IPOS.GT.55) THEN
                  HSEQ = IREC(91:91 + 54)
              ELSE
                  HSEQ = IREC(91:91 + IPOS - 1)
              ENDIF
CPDBSUM RAL 29/3/05 <--
              IF(HLXBEG(6:6).EQ.'-') HLXBEG(6:6)=' '
              IF(HLXEND(6:6).EQ.'-') HLXEND(6:6)=' '
CPDBSUM RAL 4/4/05 -->
C              IF(Y.LT.50.0) THEN
C----         If plotting all motifs (ie one-per-page, as opposed to just
C             a single page of several motifs) start a new page
              IF (.NOT.ONEPAG .AND. I.GT.1) THEN
CPDBSUM RAL 4/4/05 <--
                  CALL NEWPAG(NUMWFI,HDIAGF,PDBCOD,TTL)
CPDBSUM RAL 7/4/05 -->
Cc               CALL PRECT(500.0,760.0,100.0,50.0)
C                  IF(.NOT.BW) CALL COLOUR(HXCOL(1))
Cc               CALL TPRINT('Green: Hydrophobic',10.0,500.0,780.0)
C                  IF(.NOT.BW) CALL COLOUR(HXCOL(3))
Cc               CALL TPRINT('Blue : Polar      ',10.0,489.0,760.0)
C                  IF(.NOT.BW) CALL COLOUR(HXCOL(2))
Cc               CALL TPRINT('Red  : Charged    ',10.0,492.0,740.0)
CPDBSUM RAL 7/4/05 <--
                  Y=640.0
                  X=50.0
              ENDIF
CPDBSUM RAL 7/4/05 -->
C              CALL HNET(HLXNUM,HLXBEG,HLXEND,NUMRES,HSEQ,Y,X)
              CALL HNET(HLXNUM,HLXBEG,RSNAMB,HLXEND,RSNAME,NUMRES,HSEQ,
     -            Y,X,MOTIF)
CPDBSUM RAL 7/4/05 <--
              CALL HWHEEL(NUMRES,HSEQ,Y,X)
              IF(MOD(I,2).EQ.0) THEN
                  X=X-300.0
                  Y=Y-140.0
              ELSE
                  X=X+300.0
              ENDIF
CPDBSUM RAL 29/3/05 -->
C 150     CONTINUE
          IF (ONEPAG .AND. MOD(I,8).EQ.0) THEN
              GO TO 160
          ELSE
              IF (ONEPAG) MOTIF = MOTIF + 1
              GO TO 150
          ENDIF
CPDBSUM RAL 29/3/05 <--
 160      CONTINUE
          CALL PSCLOSE
      ENDIF

c     now do helix interactions
CPDBSUM RAL 29/3/05 -->
C      IF(HL.and.hlxnum.gt.1) THEN
C 70      READ(2,'(A)',end=90) LINE
C         IF(LINE(2:19).NE.'Helix Interactions') GO TO 70
C 
C         TTL='HINT'
C         DO 80 I=1,100
C            NUM=I
C            READ(2,55,END=90) HLX1,HLX2,DIST,ANG,NINT,NINT1,NINT2
C 55         FORMAT(I3,1X,I3,6X,F4.1,5X,F6.1,2X,I3,1X,I2,1X,I2)
C            IF(I.EQ.1) THEN
C               HINTFI=PDBCOD(1:NPS-1)//'_helix_int01.ps'
C               CALL PSOPEN(HINTFI)
C               CALL TITLE(PDBCOD,'HINT',NUMIFI,HINTFI,2)
C               Y=640.0
C            ENDIF
C            IF(Y.LT.50.0) THEN
C               CALL HINTTB(Y)
C               CALL NEWPAG(NUMIFI,HINTFI,PDBCOD,TTL)
C               Y=640.0
C            ENDIF
C            WRITE(TEXT3,'(I3)') HLX1
C            CALL TPRINT(TEXT3,12.0,75.0,Y)
C            WRITE(TEXT3,'(I3)') HLX2
C            CALL TPRINT(TEXT3,12.0,125.0,Y)
C            WRITE(TEXT4,'(F4.1)') DIST
C            CALL TPRINT(TEXT4,12.0,175.0,Y)
C            WRITE(TEXT6,'(F6.1)') ANG
C            CALL TPRINT(TEXT6,12.0,225.0,Y)
C            WRITE(TEXT3,'(I3)') NINT
C            CALL TPRINT(TEXT3,12.0,275.0,Y)
C            WRITE(TEXT3,'(I3)') NINT1
C            CALL TPRINT(TEXT3,12.0,325.0,Y)
C            WRITE(TEXT3,'(I3)') NINT2
C            CALL TPRINT(TEXT3,12.0,375.0,Y)
C            Y=Y-20.0
C 80      CONTINUE
C 90      IF(NUM.GT.1) CALL HINTTB(Y)
C 101     IF(NUM.GT.0) CALL PSCLOSE
C      ENDIF
C      GO TO 900
C 890  WRITE(*,*) ' No helix data found'
CPDBSUM RAL 29/3/05 <--
 900  RETURN
      END
C********************************************************************
C     SUBROUTINE HLXTB - DRAWS OUT A TABLE FORMAT FOR HELIX DATA
C     
C********************************************************************
      SUBROUTINE HLXTB(Y)

      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(550.0,Y,550.0,700.0)

      CALL PLINE(50.0,700.0,550.0,700.0)
      CALL PLINE(50.0,650.0,550.0,650.0)
      CALL PLINE(50.0,Y,550.0,Y)

      CALL TPRINT('Helix',12.0,75.0,685.0)
      CALL TPRINT('Number',12.0,75.0,667.0) 
      CALL TPRINT('Start',12.0,125.0,685.0)
      CALL TPRINT('Residue',12.0,125.0,667.0)
      CALL TPRINT('End',12.0,175.0,685.0)
      CALL TPRINT('Residue',12.0,175.0,667.0)
      CALL TPRINT('Helix',12.0,225.0,685.0)
      CALL TPRINT('Type',12.0,225.0,667.0)
      CALL TPRINT('No. of',12.0,275.0,685.0)
      CALL TPRINT('Residues',12.0,275.0,667.0)
      CALL TPRINT('Sequence',12.0,425.0,675.0)

      RETURN
      END
C*******************************************************************
C     SUBROUTINE HGETB - DRAWS OUT A TABLE FORMAT FOR HELIX GEOMETRY DATA
C     
C********************************************************************
      SUBROUTINE HGETB(Y)

      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(350.0,Y,350.0,700.0)
      
      CALL PLINE(50.0,700.0,350.0,700.0)
      CALL PLINE(50.0,650.0,350.0,650.0)
      CALL PLINE(50.0,Y,350.0,Y)
      CALL TPRINT('Helix',12.0,75.0,685.0)
      CALL TPRINT('Number',12.0,75.0,667.0) 
      CALL TPRINT('Length',12.0,125.0,685.0)
      CALL TPRINT('(A)',12.0,125.0,667.0)
      CALL TPRINT('Unit',12.0,175.0,685.0)
      CALL TPRINT('Rise (A)',12.0,175.0,667.0)
      CALL TPRINT('Residues',12.0,225.0,685.0)
      CALL TPRINT('per turn',12.0,225.0,667.0)
      CALL TPRINT('Pitch (A)',12.0,275.0,675.0)
      CALL TPRINT('Deviation',10.0,325.0,690.0)
      CALL TPRINT('from ideal',10.0,325.0,672.0)
      CALL TPRINT('(degrees)',10.0,325.0,660.0)

      RETURN
      END
C*******************************************************************
C     SUBROUTINE HINTTB - DRAWS OUT A TABLE FORMAT FOR HELIX INTERACTION DATA
C     
C********************************************************************
      SUBROUTINE HINTTB(Y)

      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,675.0)
      CALL PLINE(350.0,Y,350.0,675.0)
      CALL PLINE(400.0,Y,400.0,700.0)
      
      CALL PLINE(50.0,700.0,400.0,700.0)
      CALL PLINE(50.0,650.0,400.0,650.0)
      CALL PLINE(250.0,675.0,400.0,675.0)
      CALL PLINE(50.0,Y,400.0,Y)
      CALL TPRINT('Helix 1',12.0,75.0,675.0)
      CALL TPRINT('Helix 2',12.0,125.0,675.0)
      CALL TPRINT('Distance',12.0,175.0,685.0)
      CALL TPRINT('(A)',12.0,175.0,667.0)
      CALL TPRINT('Omega',12.0,225.0,685.0)
      CALL TPRINT('(degrees)',12.0,225.0,667.0)
      CALL TPRINT('No. of interacting residues',12.0,325.0,685.0)
      CALL TPRINT('Total',12.0,275.0,667.0)
      CALL TPRINT('Helix 1',12.0,325.0,667.0)
      CALL TPRINT('Helix 2',12.0,375.0,667.0)
      RETURN
      END
C*******************************************************************
CPDBSUM RAL 7/4/05 -->
C      SUBROUTINE HNET(HLXNUM,HLXBEG,HLXEND,NUMRES,HSEQ,Y,X)
      SUBROUTINE HNET(HLXNUM,HLXBEG,RSNAMB,HLXEND,RSNAME,NUMRES,HSEQ,
     -    YSTART,XSTART,MOTIF)
CPDBSUM RAL 7/4/05 <--

      INCLUDE 'psplot.par'

CPDBSUM RAL 7/4/05 -->
C      INTEGER NUMRES,I,T,RPT,P,HLXNUM,K
C      REAL Y,XN,YN,Y0,X,LEN
C      CHARACTER HSEQ*55,HLXBEG*6,HLXEND*6,STRING*26,
C     ~     HYPHOB*20
      CHARACTER    HSEQ*55, HLXBEG*6, HLXEND*6, STRING*26, HYPHOB*20
      CHARACTER*1  CH
      CHARACTER*3  NUMBER, RSNAMB, RSNAME
      CHARACTER*5  RSNUMB, RSNUME
      INTEGER      I, IPOS, IRES, LENSTR, MOTIF, MXRES, NUMRES, T, RPT,
     -             RPTLEN, P, HLXNUM, K
      REAL         FNTNET, FNTSEQ, FNTSIZ, LEN, RSHIFT, SCALE, WIDNET,
     -             XN, X, XSTART, Y, YN, Y0, YSTART
CPDBSUM RAL 7/4/05 <--

c      DATA HYPHOB/'WYMCFVILGAPTSQKNEHDR'/
      DATA HYPHOB/'AVFPMILGDEKRSTYHCNQW'/
      DATA P/6/,RPT/36/,T/50/
CPDBSUM RAL 7/4/05 -->
      DATA FNTNET, FNTSIZ, FNTSEQ / 12.0, 14.0, 12.0 /
      DATA RSHIFT, SCALE, WIDNET / 10.0, 1.5, 75.0 /
      DATA MXRES /18/,RPTLEN/10/
CPDBSUM RAL 7/4/05 <--

CPDBSUM RAL 7/4/05 -->
C      IF(.NOT.BW) CALL COLOUR(HXCOL(4))
C      WRITE(STRING,10) HLXNUM,HLXBEG,HLXEND
C 10   FORMAT('Helix ',I3,2X,'(',A6,'-',A6,')')
C      CALL TPRINT(STRING,12.0,X+90.0,Y-15.0)
C      CALL TPRINT(HSEQ(1:numres),8.0,X+90.0,Y-25.0)

C---- Initialise coordinates
      X = XSTART
      Y = YSTART - 15.0

C---- Print helix number and range, if required
      CALL COLOUR(1)
      IF (MOTIF.GT.0) THEN

C----     Get the motif number
          WRITE(NUMBER,10) MOTIF
 10       FORMAT(I3)
          CALL STRIM(NUMBER)

C----     Format the helix's start and end residues
          CALL FRMRES(RSNAMB)
          CALL FRMRES(RSNAME)
          RSNUMB = HLXBEG(2:6)
          RSNUME = HLXEND(2:6)
          CALL STRIM(RSNUMB)
          CALL STRIM(RSNUME)

C----     Form string holding helix info
          STRING = 'Helix ' // NUMBER(1:LENSTR(NUMBER)) // '. ' //
     -        RSNAMB(1:LENSTR(RSNAMB)) // RSNUMB(1:LENSTR(RSNUMB)) //
     -        '-' //
     -        RSNAME(1:LENSTR(RSNAME)) // RSNUME(1:LENSTR(RSNUME))

C----     Print the string
          CALL PSLTXT(X,Y,FNTSIZ,STRING(1:LENSTR(STRING)))

C----     Drop the y-value
          Y = Y - FNTSIZ
      ENDIF

C---- Show the sequence
      CALL PSLTXT(X,Y,FNTSEQ,'Sequence: ')
      X = X + 62.0
      IPOS = 0
      DO 100, IRES = 1, NUMRES

C----     Get this residue
          CH = HSEQ(IRES:IRES)

C----     Determine its colour
          K = INDEX(HYPHOB,CH)
          IF(K.LE.8) THEN
              CALL COLOUR(4)
          ELSE IF(K.LE.12) THEN
              CALL COLOUR(3)
          ELSE
              CALL COLOUR(5)
          ENDIF

C----     If off end of line, start a new line
          IPOS = IPOS + 1
          IF (IPOS.GT.MXRES) THEN
              X = X + 62.0
              IPOS = 1
              Y = Y - FNTSEQ
          ENDIF

C----     Write out this residue
          CALL BPRINT(CH,FNTSEQ,X,Y)

C----     Increment x-position
          X = X + RSHIFT
 100  CONTINUE

C---- Show sequence length
      X = XSTART
      Y = Y - FNTSEQ
      WRITE(STRING,120) NUMRES
 120  FORMAT('Length: ',I3,' residues')
      CALL COLOUR(1)
      CALL PSLTXT(X,Y,FNTSEQ,STRING(1:LENSTR(STRING)))

C---- Reset y-value
      Y = YSTART
      X = XSTART + 5.0

CPDBSUM RAL 7/4/05 <--
C---  Calculate length of net
      LEN = REAL(P * NUMRES * RPTLEN / RPT)
      Y0 = Y + 40.0 - LEN / 2.0

C---- Plot the residues
      DO 500 I = 1, NUMRES
          XN = X + SCALE * REAL(MOD(I*T*RPTLEN/RPT,T))
          YN = Y0 + SCALE * REAL(P*I*RPTLEN/RPT)
          K = INDEX(HYPHOB,HSEQ(I:I))
          IF (K.LE.8) THEN
              IF(.NOT.BW) CALL COLOUR(HXCOL(1))
          ELSE IF(K.LE.12) THEN
              IF(.NOT.BW) CALL COLOUR(HXCOL(2))
          ELSE
              IF(.NOT.BW) CALL COLOUR(HXCOL(3))
          ENDIF
          CALL BPRINT(HSEQ(I:I),FNTNET,XN,YN)
 500  CONTINUE
      CALL COLOUR(1)
      CALL PSLWID(1.0)
      CALL PLINE(X,Y0-6.0,X,YN+6.0)
      CALL PLINE(X+WIDNET,Y0-6.0,X+WIDNET,YN+6.0)

      RETURN
      END
C*********************************************************************
      SUBROUTINE HWHEEL(NUMRES,HSEQ,Y,X)

      INCLUDE 'psplot.par'

      INTEGER NUMRES,I,K,L
      REAL Y,XN,YN,Y0,X0,THETA,DTHETA,RPT,PI,R1,R2,R3,R4,R,
     ~     THETA1,X
      CHARACTER HSEQ*55,HYPHOB*20
CPDBSUM RAL 7/4/05 -->
      REAL FNTSIZ, RADIUS
      DATA FNTSIZ, RADIUS / 12.0, 25.0 /
CPDBSUM RAL 7/4/05 <--

      DATA RPT/3.6/,PI/3.1416/,R1/7.0/,R2/15.0/,R3/23.0/,
     ~     R4/49.0/
c      DATA HYPHOB/'WYMCFVILGAPTSQKNEHDR'/
      DATA HYPHOB/'AVFPMILGDEKRSTYHCNQW'/
c     green - hydrophobic; red - charged; blue - polar (Branden +Tooze)
      DTHETA=2*PI/RPT

      Y0=Y+40.0
      X0=X+150.0
      CALL PCIRC(X0,Y0,RADIUS)
c      CALL FARROW(X0+10.0,Y0,X0+18.0,Y0,HMNWID,HMXWID,HHDLN)
      CALL TPRINT('*',9.0,X0+RADIUS-2.0,Y0-2.0)
      THETA=DTHETA
      DO 100 I = 1, NUMRES
         THETA = THETA - DTHETA
         K = MOD(I,18)
         IF(K.EQ.15.OR.K.EQ.8.OR.K.EQ.1.OR.K.EQ.12.OR.K.EQ.5)
     ~        THEN
            THETA1=THETA-3.*PI/180.0
         ELSE IF(K.EQ.14.OR.K.EQ.3.OR.K.EQ.10.OR.K.EQ.17.OR.
     ~           K.EQ.6) THEN
            THETA1=THETA+3.*PI/180.0
         ELSE
            THETA1=THETA
         ENDIF
         IF(I.LE.18) THEN
            R = RADIUS + R1
         ELSE IF(I.LE.36) THEN
            R = RADIUS + R2
         ELSE IF(I.LE.54) THEN
            R = RADIUS + R3
         ELSE
            R = RADIUS + R4
         ENDIF
         XN=X0+R*COS(THETA1)
         YN=Y0+R*SIN(THETA1)
         L=INDEX(HYPHOB,HSEQ(I:I))
         IF(L.LE.8) THEN
            IF(.NOT.BW) CALL COLOUR(HXCOL(1))
         ELSE IF(L.LE.12) THEN
            IF(.NOT.BW) CALL COLOUR(HXCOL(2))
         ELSE
            IF(.NOT.BW) CALL COLOUR(HXCOL(3))
         ENDIF
         CALL BPRINT(HSEQ(I:I),FNTSIZ,XN,YN)

 100  CONTINUE

      RETURN
      END
C*********************************************************************
      SUBROUTINE DISULPH(PDBCOD)
      
      include 'psplot.par'

      INTEGER I,J,NUMDFI,NPS
      REAL CADIST,CHI1,CHI2,CHI3,CHI2P,CHI1P,Y
      CHARACTER DTYPE*3,LINE*80,NAMES(4)*21,RES1*6,RES2*6,
     ~     TEXT3*3,TEXT6*6,TYPES(4)*3,TTL*4
      CHARACTER*78 PDBCOD,DSFFI,DSFNM
      LOGICAL DS

      DATA DS/.FALSE./,TTL/'DSLF'/
      DATA (TYPES(I),I=1,4)/'LHS','RHS','RHH','SRH'/
      DATA (NAMES(I),I=1,4)/'  Left Hand Spiral   ',
     ~     '  Right hand spiral  ','   Right hand hook   ',
     ~     'Short right hand hook'/

      NUMDFI=1
      NPS=INDEX(PDBCOD,' ')
      DSFNM=PDBCOD(1:NPS-1)//'.dsf'
      OPEN(UNIT=2,FILE=DSFNM,STATUS='OLD',ERR=100)
      READ(2,*,end=890)

      
      DO 20 I=1,NNN
         READ(2,'(A80)',END=30) LINE
         READ(LINE,15) RES1,RES2,CHI1,CHI2,CHI3,CHI2P,CHI1P,CADIST,
     ~        DTYPE
 15      FORMAT(A6,1X,A6,5(1X,F6.1),1X,F3.1,1X,A3)
         IF(I.EQ.1) THEN
            NPS=INDEX(PDBCOD,' ')
            DSFFI=PDBCOD(1:NPS-1)//'_disulph_01.ps'
            CALL PSOPEN(DSFFI)
CPDBSUM RAL 16/3/05 -->
C            CALL TITLE(PDBCOD,'DSLF',NUMDFI,DSFFI,2)
CPDBSUM RAL 16/3/05 <--
            Y=640.0
            DS=.TRUE.
         ENDIF
         IF(Y.LT.50.0) THEN
            CALL DSFTB(Y)
            CALL NEWPAG(NUMDFI,DSFFI,PDBCOD,TTL)
            Y=640.0
         ENDIF
         CALL TPRINT(RES1,12.0,75.0,Y)
         CALL TPRINT(RES2,12.0,125.0,Y)
         WRITE(TEXT6,'(F6.1)') CHI1
         CALL TPRINT(TEXT6,12.0,175.0,Y)
         WRITE(TEXT6,'(F6.1)') CHI2
         CALL TPRINT(TEXT6,12.0,225.0,Y)
         WRITE(TEXT6,'(F6.1)') CHI3
         CALL TPRINT(TEXT6,12.0,275.0,Y)
         WRITE(TEXT6,'(F6.1)') CHI2P
         CALL TPRINT(TEXT6,12.0,325.0,Y)
         WRITE(TEXT6,'(F6.1)') CHI1P
         CALL TPRINT(TEXT6,12.0,375.0,Y)
         WRITE(TEXT3,'(F3.1)') CADIST
         CALL TPRINT(TEXT3,12.0,425.0,Y)
         DO 25 J=1,4
            IF(TYPES(J).EQ.DTYPE) THEN
               CALL TPRINT(NAMES(J),10.0,500.0,Y)
               GO TO 26
            ENDIF
 25      CONTINUE
 26      Y=Y-20.0
 20   CONTINUE
 30   IF(DS) THEN
         CALL DSFTB(Y)
         CALL PSCLOSE
      ENDIF
      GO TO 100
 890  WRITE(*,*) ' No disulphide data found'
 100  RETURN
      END
c*********************************************************************
C     SUBROUTINE DSFTB - DRAWS OUT A TABLE FORMAT FOR DISULPHIDE DATA
C     
C********************************************************************
      SUBROUTINE DSFTB(Y)

      REAL Y

      WRITE(1,5)
 5    FORMAT(' 2 setlinewidth')
      
      CALL PLINE(50.0,Y,50.0,700.0)
      CALL PLINE(100.0,Y,100.0,700.0)
      CALL PLINE(150.0,Y,150.0,700.0)
      CALL PLINE(200.0,Y,200.0,700.0)
      CALL PLINE(250.0,Y,250.0,700.0)
      CALL PLINE(300.0,Y,300.0,700.0)
      CALL PLINE(350.0,Y,350.0,700.0)
      CALL PLINE(400.0,Y,400.0,700.0)
      CALL PLINE(450.0,Y,450.0,700.0)
      CALL PLINE(550.0,Y,550.0,700.0)
      
      CALL PLINE(50.0,700.0,550.0,700.0)
      CALL PLINE(50.0,650.0,550.0,650.0)
      CALL PLINE(50.0,Y,550.0,Y)
      CALL TPRINT('Cysteine',12.0,75.0,685.0)
      CALL TPRINT('1',12.0,75.0,667.0)
      CALL TPRINT('Cysteine',12.0,125.0,685.0)
      CALL TPRINT('2',12.0,125.0,667.0)
      CALL TPRINT('Chi1',12.0,175.0,675.0)
      CALL TPRINT('Chi2',12.0,225.0,675.0)
      CALL TPRINT('Chi3',12.0,275.0,675.0)
      CALL TPRINT('Chi2''',12.0,325.0,675.0)
      CALL TPRINT('Chi1''',12.0,375.0,675.0)
      CALL TPRINT('Calpha',12.0,425.0,685.0)
      CALL TPRINT('Distance',12.0,425.0,667.0)
      CALL TPRINT('Type',12.0,500.0,675.0)
      RETURN
      END
c***********************************************************************
      SUBROUTINE SUMMRY(PDBCOD)

      INCLUDE 'psplot.par'

      INTEGER CLASS1,CLASS2,I,ILEN,NBAB,NBL,NBT,NDS,NGT,NHP,NHX,
     ~     NPSI,NPS,NSHT,NST,NUMFI,NUMSTR,TNEED,IOCODE,SHTLT
      REAL X0,Y,Y0,YB,YBT,YDS,YG,YH,YHX,YS,YT
      CHARACTER RES(4)*6,SEQ(4),TURNTY*4,TEXT13*13,
     ~     TEXT4*4,GTRNTY*7,TEXT3*3,LINE*80,DASH,LINE2*125,AP,BTYP,
     ~     TEXT6*6,HLXBEG*6,HLXEND*6,HTYPE,SHETYP,BARREL,TEXT2*2,
     ~     STBEG*6,STEND*6,RES1I*6,RES2E*6,ECLASS*2,TEXT8*8,
     ~     COLON,SSSUM*80,TPLGY*120,TITL*80,TTL*4,TEXT9*9
      CHARACTER*78 PDBCOD,SUMNAM,SNAM

      DATA DASH/'-'/,COLON/':'/,TTL/'SUMM'/
      NUMFI=1
      X0=80.0

c     read in  secondary structure summary information first
c      IF(.NOT.BW) CALL COLOUR(4)
      NPS=INDEX(PDBCOD,' ')
      SNAM=PDBCOD(1:NPS-1)//'.sum'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=1900)
c     open output summary file
      SUMNAM=PDBCOD(1:NPS-1)//'_summary_01.ps'
      CALL PSOPEN(SUMNAM)
      IF(.NOT.BW) CALL COLOUR(3)
CPDBSUM RAL 16/3/05 -->
C      CALL TITLE(PDBCOD,'SUMM',NUMFI,SUMNAM,2)
CPDBSUM RAL 16/3/05 <--
      Y0=720.0

      READ(2,*,END=890)
      READ(2,*,END=890)
C     read title of protein
      READ(2,'(A80)',END=890) TITL
      CALL TPRINT(TITL,11.0,300.0,Y0)
      Y0=Y0-17.0
c     read sequence
      DO 10 I=1,75
         READ(2,'(A80)',END=11) SSSUM
         IF(I.EQ.1) THEN
            CALL TPRINT('Sequence: ',12.0,100.0,Y0)
            Y0=Y0-17.0
            CALL TPRINT(SSSUM,9.0,300.0,Y0)
            Y0=Y0-12.0
         ELSE
            IF(SSSUM(1:6).NE.'Number') THEN
               if(y0.lt.50.0) then 
                  call newpag(numfi,sumnam,pdbcod,ttl)
                  y0=720.0
               endif
               CALL TPRINT(SSSUM,9.0,300.0,Y0)
               Y0=Y0-12.0
            ELSE
               Y0=Y0-5.0
               CALL PLINE(100.0,Y0+8.0,500.0,Y0+8.0)
               GO TO 11
            ENDIF
         ENDIF
 10   CONTINUE
 11   CALL TPRINT(SSSUM(1:62),9.0,300.0,Y0)
      read(2,*,end=890)
      Y0=Y0-17.0
      if(y0.lt.50.0) then
         call newpag(numfi,sumnam,pdbcod,ttl)
         y0=720.0
      endif
      READ(2,'(A80)',end=890) SSSUM
      CALL TPRINT(SSSUM(1:25),9.0,150.0,Y0)
      READ(2,'(A80)',end=890) SSSUM
      CALL TPRINT(SSSUM(1:31),9.0,300.0,Y0)
      READ(2,'(A80)',end=890) SSSUM
      CALL TPRINT(SSSUM(1:30),9.0,450.0,Y0)
      READ(2,'(A80)',end=890) SSSUM
      Y0=Y0-17.0
      if(y0.lt.50.0) then
         call newpag(numfi,sumnam,pdbcod,ttl)
         y0=720.0
      endif
      CALL TPRINT(SSSUM(1:30),9.0,150.0,Y0)
      READ(2,'(A80)',end=890) SSSUM
      CALL TPRINT(SSSUM(1:35),9.0,300.0,Y0)
      READ(2,'(A80)',end=890) SSSUM
      CALL TPRINT(SSSUM(1:34),9.0,450.0,Y0)
      IF(.NOT.BW) CALL COLOUR(5)
      Y0=Y0-17.0
      if(y0.lt.50.0) then
         call newpag(numfi,sumnam,pdbcod,ttl)
         y0=720.0
      endif
      CALL PLINE(100.0,Y0+8.0,500.0,Y0+8.0)
      CALL TPRINT('Secondary Structure Summary: ',12.0,145.0,Y0)
      Y0=Y0-17.0
      if(y0.lt.50.0) then
         call newpag(numfi,sumnam,pdbcod,ttl)
         y0=720.0
      endif
      DO 16 I=1,20
         READ(2,'(A80)',END=17) SSSUM
         ILEN=INDEX(SSSUM,' ')
         IF(ILEN.EQ.1) THEN
            GO TO 17
         ELSE IF(ILEN.EQ.0) THEN
            CALL TPRINT(SSSUM,9.0,300.0,Y0)
            Y0=Y0-15.0
         ELSE
            CALL TPRINT(SSSUM(1:ILEN-1),9.0,300.0,Y0)
            Y0=Y0-15.0
            if(y0.lt.50.0) then
               call newpag(numfi,sumnam,pdbcod,ttl)
               y0=720.0
            endif
         ENDIF
 16   CONTINUE

 17   YS=Y0
      YT=YS
      GO TO 1900
 890  WRITE(*,*) ' No secondary structure summary data'
C     now do sheets and strands
 1900 IF(.NOT.BW) CALL COLOUR(12)
      Y0=YS-30.0
      Y=Y0-25.0
      SNAM=PDBCOD(1:NPS-1)//'.sht'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=400)
      NSHT=0
C 70   READ(2,'(A)') LINE
C      IF(LINE(1:3).NE.'The') GO TO 70
      DO 80 I=1,NNN
         READ(2,'(A)',END=85) LINE
         IF(LINE(3:3).EQ.' ') GO TO 85
C         IF(LINE(8:8).EQ.'Y'.OR.LINE(8:8).EQ.'N') THEN
            READ(LINE,75) SHTLT,NUMSTR,SHETYP,BARREL
 75         FORMAT(I3,1X,I2,1X,A1,1X,A1)
            READ(2,'(10x,A120)') TPLGY
            NSHT=NSHT+1
            WRITE(TEXT2,'(I2)') NUMSTR
            CALL TPRINT(TEXT2,9.0,150.0,Y)
            WRITE(TEXT3,'(I3)') SHTLT 
            CALL TPRINT(TEXT3,9.0,175.0,Y)
            IF(SHETYP.EQ.'A') THEN
               CALL TPRINT('ANTI',9.0,205.0,Y)
            ELSE IF(SHETYP.EQ.'P') THEN
               CALL TPRINT('PAR',9.0,205.0,Y)
            ELSE IF(SHETYP.EQ.'M') THEN
                CALL TPRINT('MIX',9.0,205.0,Y)
            ENDIF
            CALL TPRINT(BARREL,9.0,240.0,Y)
            ILEN=INDEX(TPLGY,'      ')
            IF(ILEN.EQ.0) THEN
               CALL TPRINT(TPLGY,7.0,365.0,Y)
            ELSE
               CALL TPRINT(TPLGY(1:ILEN-1),7.0,365.0,Y)
            ENDIF
            Y=Y-15.0
            IF(Y.LT.50.0) THEN
               IF(NSHT.GT.0) CALL SSMTB(Y,Y0)
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               Y=695.0
            ENDIF
C         ENDIF
 80   CONTINUE
 85   IF(NSHT.GT.0) CALL SSMTB(Y,Y0)
      YT=Y

c     now do strands
      TNEED=0
      IF(.NOT.BW) CALL COLOUR(15)
      NST=0
      SNAM=PDBCOD(1:NPS-1)//'.str'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD')
C      READ(2,'(A)') LINE
      DO 180 I=1,NNN
         READ(2,'(A)',END=190) LINE
         IF(LINE(4:4).EQ.' ') GO TO 190
         READ(LINE,95) STBEG,STEND,SHTLT
c         TNEED=1
         NST=NST+1
c         IF(NST.EQ.1) THEN
         if(i.eq.1) then
            Y0=YT-30.0
            Y=Y0-20.0
         ENDIF
 95      FORMAT(5X,A6,1X,A6,1X,I3)
         IF(Y.LT.50.0) THEN
c            IF(NST.GT.1) CALL STSMTB(Y0,Y,X0)
            if(tneed.eq.1) then
               call stsmtb(y0,y,x0)
               TNEED=0
            endif
            Y0=YT-30.0
            Y=Y0-20.0
            X0=X0+140.0
            IF(X0.GT.500.0) THEN
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               X0=80.0
               Y=Y0-20.0
               YT=Y0+30.0
               IF(.NOT.BW) CALL COLOUR(15)
            ENDIF
         ENDIF
C         ELSE
            WRITE(TEXT13,'(A6,A1,A6)') STBEG,DASH,STEND
            CALL TPRINT(TEXT13,9.0,X0-10.0,Y)
            WRITE(TEXT3,'(I3)') SHTLT
            CALL TPRINT(TEXT3,9.0,X0+25.0,Y)
            Y=Y-15.0
            tneed=1
C         ENDIF
 180  CONTINUE
 190  IF(TNEED.EQ.1) CALL STSMTB(Y0,Y,X0)
      YS=Y
      CLOSE(2)

c     now do helices
 400  IF(.NOT.BW) CALL COLOUR(9)
      SNAM=PDBCOD(1:NPS-1)//'.hlx'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=500)
      NHX=0
      TNEED=0
 110  READ(2,'(A125)',end=490) LINE2
      IF(LINE2(2:11).NE.'Helix data') GO TO 110
      DO 130 I=1,NNN
         READ(2,'(A125)',END=140) LINE2
         IF(LINE2(3:3).EQ.' ') GO TO 140
         READ(LINE2,115) HLXBEG,HLXEND,HTYPE
C         NHX=NHX+1
C         IF(NHX.EQ.1) THEN
C            IF(TNEED.EQ.0) THEN
C               Y0=YT-30.0
C            ELSE
C               Y0=YS-30.0
C            ENDIF
C            Y=Y0-20.0
C         ENDIF
C         TNEED=1
         IF(I.EQ.1) THEN
            Y0=Y-30.0
            Y=Y0-20.0
         ENDIF
 115     FORMAT(6X,A6,1X,A6,1X,A1)
         IF(Y.LT.50.0) THEN
            IF(TNEED.EQ.1) THEN 
               CALL HSMTB(Y,Y0,X0)
               TNEED=0
            ENDIF
            X0=X0+140.0
            Y0=YT-30.0
            Y=Y0-20.0
c24/4/96            NHX=0
            IF(X0.GT.500.0) THEN
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               X0=80.0
               Y=Y0-20.0
               YT=Y0+30.0
               IF(.NOT.BW) CALL COLOUR(9)
            ENDIF
         ENDIF
C         ELSE
            WRITE(TEXT13,'(A6,A1,A6)') HLXBEG,DASH,HLXEND      
            CALL TPRINT(TEXT13,9.0,X0-10.0,Y)
            IF(HTYPE.EQ.'G') HTYPE='3'
            CALL TPRINT(HTYPE,9.0,X0+25.0,Y)
            Y=Y-15.0
            TNEED=1
C         ENDIF
 130  CONTINUE
c24/4/96 140  IF(TNEED.EQ.1.AND.NHX.GT.0) CALL HSMTB(Y,Y0,X0)
 140  IF(TNEED.EQ.1) CALL HSMTB(Y,Y0,X0)
      YHX=Y   
 490  CLOSE(UNIT=2)

C     now disulphide bridges
 500  IF(.NOT.BW) CALL COLOUR(17)
      SNAM=PDBCOD(1:NPS-1)//'.dsf'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=600)
      NDS=0
      tneed=0
      READ(2,*,end=590)
      DO 300 I=1,NNN
         READ(2,'(A6,1X,A6)',END=310) RES1I,RES2E
         IF(RES1I.EQ.'      ') GO TO 310
         NDS=NDS+1
c         IF(NDS.EQ.1) THEN
c            IF(TNEED.EQ.0) THEN
c               Y0=YT-30.0
c            ELSE
c               Y0=YHX-30.0
c            ENDIF
         if(i.eq.1) then
            y0=y-30.0
            Y=Y0-20.0
         ENDIF
c         TNEED=1
         IF(Y.LT.50.0) THEN
c            IF(NDS.GT.1) CALL DSSMTB(Y,Y0,X0)
            if(tneed.eq.1) then
               call dssmtb(y,y0,x0)
               TNEED=0
            endif
            X0=X0+140.0
            Y0=YT-30.0
            Y=Y0-20.0
            IF(X0.GT.500.0) THEN
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               X0=80.0
               Y=Y0-20.0
               YT=Y0+30.0
               IF(.NOT.BW) CALL COLOUR(17)
            ENDIF
C           WRITE(TEXT6,'(A6)') RES1I
C            CALL TPRINT(TEXT6,9.0,X0-20.0,Y)
C            WRITE(TEXT6,'(A6)') RES2E
C            CALL TPRINT(TEXT6,9.0,X0+20.0,Y)
C            Y=Y-15.0
         ENDIF
C         ELSE
            WRITE(TEXT6,'(A6)') RES1I
            CALL TPRINT(TEXT6,9.0,X0-20.0,Y)
            WRITE(TEXT6,'(A6)') RES2E
            CALL TPRINT(TEXT6,9.0,X0+20.0,Y)
            Y=Y-15.0
            tneed=1
C         ENDIF
 300  CONTINUE
c 310  IF(TNEED.EQ.1.AND.NDS.GT.0) CALL DSSMTB(Y,Y0,X0)
 310  IF(TNEED.EQ.1) CALL DSSMTB(Y,Y0,X0)
      YDS=Y
 590  CLOSE(2)

c     beta turns first
      
 600  IF(.NOT.BW) CALL COLOUR(4)
      SNAM=PDBCOD(1:NPS-1)//'.bturns'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=700)
      NBT=0
      tneed=0
      READ(2,*,end=690)
      READ(2,*,end=690)
      DO 20 I=1,NNN
         READ(2,15,END=30) RES(1),SEQ(1),SEQ(2),SEQ(3),RES(4),SEQ(4),
     ~        TURNTY
 15      FORMAT(A6,A1,6X,A1,6X,A1,A6,A1,1X,A4)
         IF(RES(1).EQ.'      ') GO TO 30
         NBT=NBT+1
c         IF(NBT.EQ.1) THEN
c            IF(X0.EQ.80.0) THEN
c               X0=X0+140.0
c               Y0=YT-30.0
c            ELSE
c               IF(TNEED.EQ.0) THEN
c                  Y0=YT-30.0
c               ELSE
c                  Y0=YDS-30.0
c               ENDIF
c            ENDIF
         if(i.eq.1) then
            y0=y-30.0
            Y=Y0-20.0
         ENDIF
c         TNEED=1
         IF(Y.LT.50.0) THEN
c            IF(NBT.GT.1) CALL BSUMTB(Y,Y0,X0)
            if(tneed.eq.1) then
               call bsumtb(y,y0,x0)
               TNEED=0
            endif
            X0=X0+140.0
            Y0=YT-30.0
            Y=Y0-20.0
            IF(X0.GT.500.0) THEN
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               X0=80.0
               Y=Y0-20.0
               YT=Y0+30.0
               IF(.NOT.BW) CALL COLOUR(4)
            ENDIF
         ENDIF
C         ELSE
            WRITE(TEXT13,'(A6,A1,A6)') RES(1),DASH,RES(4)
            WRITE(TEXT4,'(4A1)') SEQ(1),SEQ(2),SEQ(3),SEQ(4)
            CALL TPRINT(TEXT13,9.0,X0-30.0,Y)
            CALL TPRINT(TEXT4,9.0,X0+10.0,Y)
            CALL TPRINT(TURNTY,9.0,X0+40.0,Y)
            Y=Y-15.0
            tneed=1
C         ENDIF
 20   CONTINUE
c 30   IF(TNEED.EQ.1.AND.NBT.GT.0) CALL BSUMTB(Y,Y0,X0)
 30   IF(TNEED.EQ.1) CALL BSUMTB(Y,Y0,X0)
      YBT=Y
 690  CLOSE(2)

c     now gamma turns
      IF(.NOT.BW) CALL COLOUR(4)
      SNAM=PDBCOD(1:NPS-1)//'.gturns'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD')
 40   READ(2,'(A80)',end=790) LINE
      IF(LINE(2:6).NE.'Gamma') GO TO 40
      NGT=0
      tneed=0
      DO 50 I=1,NNN
         READ(2,25,END=60) RES(1),SEQ(1),SEQ(2),RES(3),SEQ(3),GTRNTY
 25      FORMAT(A6,A1,6X,A1,A6,A1,1X,A7)
         IF(RES(1).EQ.'      ') GO TO 50
         NGT=NGT+1
c         IF(NGT.EQ.1) THEN
c            IF(X0.EQ.220.0) THEN
c               X0=X0+140.0
c               Y0=YT-30.0
c            ELSE
c               IF(TNEED.EQ.0) THEN
c                  Y0=YT-30.0
c               ELSE
c                  Y0=YBT-30.0
c               ENDIF
c            ENDIF
c            Y=Y0-20.0
c         ENDIF
c         TNEED=1
         if(i.eq.1) then
            y0=y-30.0
            y=y0-20.0
         endif
         IF(Y.LT.50.0) THEN
c     IF(NGT.GT.1) CALL GSUMTB(Y,Y0,X0)
            if(tneed.eq.1) then
               call gsumtb(y,y0,x0)
               TNEED=0
            endif
            X0=X0+140.0
            Y0=YT-30.0
            Y=Y0-20.0
            IF(X0.GT.500.0) THEN
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               X0=80.0
               Y=Y0-20.0
               YT=Y0+30.0
               IF(.NOT.BW) CALL COLOUR(4)
            ENDIF
         ENDIF
C         ELSE
            WRITE(TEXT13,'(A6,A1,A6)') RES(1),DASH,RES(3)
            WRITE(TEXT3,'(3A1)') SEQ(1),SEQ(2),SEQ(3)
            CALL TPRINT(TEXT13,9.0,X0-30.0,Y)
            CALL TPRINT(TEXT3,9.0,X0+10.0,Y)
            CALL TPRINT(GTRNTY(1:3),9.0,X0+40.0,Y)
            Y=Y-15.0
            tneed=1
C         ENDIF
 50   CONTINUE
c 60   IF(TNEED.EQ.1.AND.NGT.GT.0) CALL GSUMTB(Y,Y0,X0)
 60   IF(TNEED.EQ.1) CALL GSUMTB(Y,Y0,X0)
      YG=Y
 790  CLOSE(UNIT=2)

c     now do bulges
 700  IF(.NOT.BW) CALL COLOUR(14)
      SNAM=PDBCOD(1:NPS-1)//'.blg'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=800)
      NBL=0
      tneed=0
C      READ(2,*)
      DO 90 I=1,NNN
         READ(2,'(A123)',END=100) LINE2
         IF(LINE2(1:1).EQ.'A'.OR.LINE2(1:1).EQ.'P') THEN
            READ(LINE2,35,END=100) AP,BTYP,RES(1),SEQ(1),RES(2),SEQ(2),
     ~           RES(3),SEQ(3)
 35         FORMAT(A1,1X,A1,3(1X,A6,1X,A1,15X))
            IF(RES(2).EQ.'      ') GO TO 90
            NBL=NBL+1
c            IF(NBL.EQ.1) THEN
c               IF(TNEED.EQ.0) THEN
c                  Y0=YT-30.0
c               ELSE
c                  Y0=YG-30.0
c               ENDIF
c               Y=Y0-20.0
c            ENDIF
c            TNEED=1
            if(i.eq.1) then
               y0=y-30.0
               y=y0-20.0
            endif
            IF(Y.LT.50.0) THEN
c               IF(NBL.GT.1) CALL BLSMTB(Y,Y0,X0)
c               TNEED=0
               if(tneed.eq.1) then
                  call blsmtb(y,y0,x0)
                  tneed=0
               endif
               X0=X0+140.0
               Y0=YT-30.0
               Y=Y0-20.0
               IF(X0.GT.500.0) THEN
                  CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
                  Y0=720.0
                  X0=80.0
                  Y=Y0-20.0
                  YT=Y0+30.0
                  IF(.NOT.BW) CALL COLOUR(14)
               ENDIF
            ENDIF
C            ELSE
               WRITE(TEXT6,'(A6)') RES(1)
               CALL TPRINT(TEXT6,9.0,X0-52.0,Y)
               WRITE(TEXT6,'(A6)') RES(2)
               CALL TPRINT(TEXT6,9.0,X0-25.0,Y)
               WRITE(TEXT6,'(A6)') RES(3)
               CALL TPRINT(TEXT6,9.0,X0+2.0,Y)
               WRITE(TEXT3,'(A1,1X,A1)') AP,BTYP
               CALL TPRINT(TEXT3,9.0,X0+27.0,Y)
               WRITE(TEXT3,'(3A1)') SEQ(1),SEQ(2),SEQ(3)
               CALL TPRINT(TEXT3,9.0,X0+52.0,Y)
               Y=Y-15.0
               tneed=1
C            ENDIF
         ENDIF
 90   CONTINUE
c 100  IF(TNEED.EQ.1.AND.NBL.GT.0) CALL BLSMTB(Y,Y0,X0)
 100  IF(TNEED.EQ.1) CALL BLSMTB(Y,Y0,X0)
      YB=Y
      CLOSE(UNIT=2)
      
c     now do hairpins
 800  IF(.NOT.BW) CALL COLOUR(7)
      SNAM=PDBCOD(1:NPS-1)//'.hpin'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',ERR=850)
      NHP=0
      tneed=0
      READ(2,*,end=990)
      DO 150 I=1,NNN
         READ(2,125,END=160) RES1I,RES2E,CLASS1,CLASS2,ECLASS
 125     FORMAT(8X,A6,15X,A6,1X,I3,1X,I3,1X,A3)
         READ(2,*)
         READ(2,*)
         READ(2,*)
         READ(2,*)
         READ(2,*)
         NHP=NHP+1
c         IF(NHP.EQ.1) THEN
c            IF(TNEED.EQ.0) THEN
c               Y0=YT-30.0
c            ELSE
c               Y0=YB-30.0
c            ENDIF
c            Y=Y0-20.0
c         ENDIF
c         TNEED=1
         if(i.eq.1) then
            y0=y-30.0
            y=y0-20.0
         endif
         IF(Y.LT.50.0) THEN
c            IF(NHP.GT.1) CALL HPSMTB(Y,Y0,X0)
c            TNEED=0
            if(tneed.eq.1) then
               call hpsmtb(y,y0,x0)
               tneed=0
            endif
            X0=X0+140.0
            Y0=YT-30.0
            Y=Y0-20.0
            IF(X0.GT.500.0) THEN
               CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
               Y0=720.0
               X0=80.0
               Y=Y0-20.0
               YT=Y0+30.0
               IF(.NOT.BW) CALL COLOUR(7)
            ENDIF
         ENDIF
C         ELSE
         TNEED=1
            WRITE(TEXT13,'(A6,A1,A6)') RES1I,DASH,RES2E
            CALL TPRINT(TEXT13,9.0,X0-25.0,Y)
            WRITE(TEXT9,'(I2,A1,I2,1X,A3)') CLASS1,COLON,CLASS2,ECLASS
            CALL TPRINT(TEXT9,9.0,X0+25.0,Y)
            Y=Y-15.0
            tneed=1
C         ENDIF
 150  CONTINUE
c     160        IF(TNEED.EQ.1.AND.NHP.GT.0) CALL HPSMTB(Y,Y0,X0)
 160  IF(TNEED.EQ.1) CALL HPSMTB(Y,Y0,X0)
C      IF(NHP.GT.0) CALL HPSMTB(Y,Y0,X0)
 990  CLOSE(UNIT=2)
      YH=Y

c     now do bab and psi loops
 850  IF(.NOT.BW) CALL COLOUR(8)
      SNAM=PDBCOD(1:NPS-1)//'.bab'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',IOSTAT=IOCODE)
      IF(IOCODE.EQ.0) THEN
         NBAB=0
         tneed=0
C     265  READ(2,'(A)') LINE
C     IF(LINE(2:5).NE.'BETA') GO TO 265
         DO 280 I=1,NNN
            READ(2,275,END=290) RES1I,RES2E
 275        FORMAT(A6,17X,A6)
            IF(RES1I.EQ.'      ') GO TO 290
            NBAB=NBAB+1
            IF(NBAB.EQ.1) THEN
c               IF(TNEED.EQ.0) THEN
c                  Y0=YT-30.0
c               ELSE
c                  Y0=YH-30.0
c               ENDIF
               y0=y-30.0
               Y=Y0-20.0
            ENDIF
c            TNEED=1
            IF(Y.LT.50.0) THEN
c               IF(NBAB.GT.1) CALL BBSMTB(Y,Y0,X0)
c               TNEED=0
               if(tneed.eq.1) then
                  call bbsmtb(y,y0,x0)
                  tneed=0
               endif
               X0=X0+140.0
               Y0=YT-30.0
               Y=Y0-20.0
               IF(X0.GT.500.0) THEN
                  CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
                  Y0=720.0
                  X0=80.0
                  Y=Y0-20.0
                  YT=Y0+30.0
                  IF(.NOT.BW) CALL COLOUR(8)
               ENDIF
            ENDIF
C     ELSE
            WRITE(TEXT6,'(A6)') RES1I
            CALL TPRINT(TEXT6,9.0,X0-25.0,Y)
            WRITE(TEXT6,'(A6)') RES2E
            CALL TPRINT(TEXT6,9.0,X0+25.0,Y)
            Y=Y-15.0
            tneed=1
C     ENDIF
 280     CONTINUE
c 290     IF(TNEED.EQ.1.AND.NBAB.GT.0) CALL BBSMTB(Y,Y0,X0)
 290     IF(TNEED.EQ.1) CALL BBSMTB(Y,Y0,X0)
         CLOSE(2)
      ENDIF
         
C     psi loops
 900     SNAM=PDBCOD(1:NPS-1)//'.psi'
      OPEN(UNIT=2,FILE=SNAM,STATUS='OLD',IOSTAT=IOCODE)
      IF(IOCODE.EQ.0) THEN
         IF(.NOT.BW) CALL COLOUR(16)
         NPSI=0
         tneed=0
C 195  READ(2,'(A)') LINE
C      IF(LINE(2:4).NE.'PSI') GO TO 195
         DO 210 I=1,NNN
            READ(2,275,END=220) RES1I,RES2E
            IF(RES1I.EQ.'      ') GO TO 220
            NPSI=NPSI+1
            IF(NPSI.EQ.1) THEN
c               IF(TNEED.EQ.0) THEN
c                  Y0=YT-30.0
c               ELSE
c                  Y0=Y-30.0
c               ENDIF
               y0=y-30.0
               Y=Y0-20.0
            ENDIF
c            TNEED=1
            IF(Y.LT.50.0) THEN
c               IF(NPSI.GT.1) CALL PSSMTB(Y,Y0,X0)
c               TNEED=0
               if(tneed.eq.1) then
                  call pssmtb(y,y0,x0)
                  tneed=0
               endif
               X0=X0+140.0
               Y0=YT-30.0
               Y=Y0-20.0
               IF(X0.GT.500.0) THEN
                  CALL NEWPAG(NUMFI,SUMNAM,PDBCOD,TTL)
                  Y0=720.0
                  X0=80.0
                  Y=Y0-20.0
                  YT=Y0+30.0
                  IF(.NOT.BW) CALL COLOUR(16)
               ENDIF
            ENDIF
C     ELSE
            WRITE(TEXT6,'(A6)') RES1I
            CALL TPRINT(TEXT6,9.0,X0-25.0,Y)
            WRITE(TEXT6,'(A6)') RES2E
            CALL TPRINT(TEXT6,9.0,X0+25.0,Y)
            Y=Y-15.0
            tneed=1
C     ENDIF
 210     CONTINUE
c 220     IF(TNEED.EQ.1.AND.NPSI.GT.0) CALL PSSMTB(Y,Y0,X0)
 220     IF(TNEED.EQ.1) CALL PSSMTB(Y,Y0,X0)
         CLOSE(2)
      endif
      CALL PSCLOSE
      RETURN
      END
C********************************************************************
      SUBROUTINE newpag(numfi,sumnam,pdbcod,ttl)
      
      INCLUDE 'psplot.par'

      INTEGER J,NUMFI
      REAL X0,Y0
      CHARACTER C,PDBCOD*78,SUMNAM*78,TTL*4,C2*2
CPDBSUM RAL 4/4/05 -->
      CHARACTER*3 C3
CPDBSUM RAL 4/4/05 <--

C     to close existing ps output file and open another

      CALL PSCLOSE
CPDBSUM RAL 4/4/05 -->
C      J=INDEX(SUMNAM,'.')
C      NUMFI=NUMFI+1
C      if(numfi.lt.10) then
C         WRITE(C,'(I1)') NUMFI
C         SUMNAM(J-1:J-1)=C
C      else
C         WRITE(C2,'(I2)') NUMFI
C         SUMNAM(J-2:J-1)=C2
C      ENDIF
      NUMFI=NUMFI+1
      J=INDEX(SUMNAM,'.ps')
      WRITE(C3,'(I3.3)') NUMFI
      SUMNAM(J-3:J-1) = C3
CPDBSUM RAL 4/4/05 <--
      CALL PSOPEN(SUMNAM)
      IF(.NOT.BW) CALL COLOUR(3)
CPDBSUM RAL 16/3/05 -->
C      CALL TITLE(PDBCOD,TTL,NUMFI,SUMNAM,2)
CPDBSUM RAL 16/3/05 <--
      Y0=720.0
      X0=80.0


      RETURN
      END
C********************************************************************
      SUBROUTINE BSUMTB(Y,Y0,X0)
      
      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-55.0,Y+YDEL,X0-55.0,Y0)
      CALL PLINE(X0-5.0,Y+YDEL,X0-5.0,Y0)
      CALL PLINE(X0+25.0,Y+YDEL,X0+25.0,Y0)
      CALL PLINE(X0+55.0,Y+YDEL,X0+55.0,Y0)
      CALL PLINE(X0-55.0,Y+YDEL,X0+55.0,Y+YDEL)
      CALL PLINE(X0-55.0,Y0-15.0,X0+55.0,Y0-15.0)
      CALL PLINE(X0-55.0,Y0,X0+55.0,Y0)
      CALL TPRINT('Location',9.0,X0-30.0,Y0-10.0)
      CALL TPRINT('Seq',9.0,X0+10.0,Y0-10.0)
      CALL TPRINT('Type',9.0,X0+40.0,Y0-10.0)
      CALL TPRINT('B-TURNS',12.0,X0,Y0+10.0)
      RETURN
      END
C*************************************************************************
      SUBROUTINE GSUMTB(Y,Y0,X0)
      
      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-55.0,Y+YDEL,X0-55.0,Y0)
      CALL PLINE(X0-5.0,Y+YDEL,X0-5.0,Y0)
      CALL PLINE(X0+25.0,Y+YDEL,X0+25.0,Y0)
      CALL PLINE(X0+55.0,Y+YDEL,X0+55.0,Y0)
      CALL PLINE(X0-55.0,Y+YDEL,X0+55.0,Y+YDEL)
      CALL PLINE(X0-55.0,Y0-15.0,X0+55.0,Y0-15.0)
      CALL PLINE(X0-55.0,Y0,X0+55.0,Y0)
      CALL TPRINT('Location',9.0,X0-30.0,Y0-10.0)
      CALL TPRINT('Seq',9.0,X0+10.0,Y0-10.0)
      CALL TPRINT('Type',9.0,X0+40.0,Y0-10.0)
      CALL TPRINT('G-TURNS',12.0,X0,Y0+10.0)
      RETURN
      END
C*************************************************************************
      SUBROUTINE BLSMTB(Y,Y0,X0)
      
      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-65.0,Y+YDEL,X0-65.0,Y0)
      CALL PLINE(X0+15.0,Y+YDEL,X0+15.0,Y0)
      CALL PLINE(X0+40.0,Y+YDEL,X0+40.0,Y0)
      CALL PLINE(X0+65.0,Y+YDEL,X0+65.0,Y0)
      CALL PLINE(X0-65.0,Y+YDEL,X0+65.0,Y+YDEL)
      CALL PLINE(X0-65.0,Y0-15.0,X0+65.0,Y0-15.0)
      CALL PLINE(X0-65.0,Y0,X0+65.0,Y0)
      CALL TPRINT('Location',9.0,X0-25.0,Y0-10.0)
      CALL TPRINT('Type',9.0,X0+27.0,Y0-10.0)
      CALL TPRINT('Seq',9.0,X0+52.0,Y0-10.0)
      CALL TPRINT('BULGES',12.0,X0,Y0+10.0)
      RETURN
      END
C*************************************************************************
      SUBROUTINE HSMTB(Y,Y0,X0)
      
      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-35.0,Y+YDEL,X0-35.0,Y0)
      CALL PLINE(X0+15.0,Y+YDEL,X0+15.0,Y0)
      CALL PLINE(X0+35.0,Y+YDEL,X0+35.0,Y0)
      CALL PLINE(X0-35.0,Y0-15.0,X0+35.0,Y0-15.0)
      CALL PLINE(X0-35.0,Y0,X0+35.0,Y0)
      CALL PLINE(X0-35.0,Y+YDEL,X0+35.0,Y+YDEL)
      CALL TPRINT('Location',9.0,X0-10.0,Y0-10.0)
      CALL TPRINT('Type',9.0,X0+25.0,Y0-10.0)
      CALL TPRINT('HELICES',12.0,X0,Y0+10.0)
      RETURN
      END
C**************************************************************************
      SUBROUTINE SSMTB(Y,Y0)
      
      REAL Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(125.0,Y+YDEL,125.0,Y0)
      CALL PLINE(165.0,Y+YDEL,165.0,Y0)
      CALL PLINE(185.0,Y+YDEL,185.0,Y0)
      CALL PLINE(225.0,Y+YDEL,225.0,Y0)
      CALL PLINE(255.0,Y+YDEL,255.0,Y0)
      CALL PLINE(475.0,Y+YDEL,475.0,Y0)
      CALL PLINE(125.0,Y0-20.0,475.0,Y0-20.0)
      CALL PLINE(125.0,Y0,475.0,Y0)
      CALL PLINE(125.0,Y+YDEL,475.0,Y+YDEL)
      CALL TPRINT('Number',9.0,150.0,Y0-5.0)
      CALL TPRINT('strands',9.0,150.0,Y0-15.0)
      CALL TPRINT('Name',9.0,175.0,Y0-15.0)
      CALL TPRINT('Type',9.0,205.0,Y0-15.0)
      CALL TPRINT('Barrel',9.0,240.0,Y0-15.0)
      CALL TPRINT('Topology',9.0,365.0,Y0-15.0)
      CALL TPRINT('SHEETS',12.0,300.0,Y0+10.0)
      call pline(100.0,y0+22.0,500.0,y0+22.0)
      RETURN
      END
C**************************************************************************
      SUBROUTINE STSMTB(Y0,Y,X0)

      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-35.0,Y+YDEL,X0-35.0,Y0)
      CALL PLINE(X0+15.0,Y+YDEL,X0+15.0,Y0)
      CALL PLINE(X0+35.0,Y+YDEL,X0+35.0,Y0)
      CALL PLINE(X0-35.0,Y+YDEL,X0+35.0,Y+YDEL)
      CALL PLINE(X0-35.0,Y0-15.0,X0+35.0,Y0-15.0)
      CALL PLINE(X0-35.0,Y0,X0+35.0,Y0)
      CALL TPRINT('Location',9.0,X0-10.0,Y0-10.0)
      CALL TPRINT('Sheet',9.0,X0+25.0,Y0-10.0)
      CALL TPRINT('STRANDS',12.0,X0,Y0+10.0)
      RETURN
      END
C**************************************************************************
      SUBROUTINE HPSMTB(Y,Y0,X0)

      REAL Y,Y0,X0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-50.0,Y+YDEL,X0-50.0,Y0)
      CALL PLINE(X0,Y+YDEL,X0,Y0)
      CALL PLINE(X0+50.0,Y+YDEL,X0+50.0,Y0)
      CALL PLINE(X0-50.0,Y+YDEL,X0+50.0,Y+YDEL)
      CALL PLINE(X0-50.0,Y0-15.0,X0+50.0,Y0-15.0)
      CALL PLINE(X0-50.0,Y0,X0+50.0,Y0)
      CALL TPRINT('Location',9.0,X0-25.0,Y0-10.0)
      CALL TPRINT('Class',9.0,X0+25.0,Y0-10.0)
      CALL TPRINT('HAIRPINS',12.0,X0,Y0+10.0)
      RETURN
      END
C**************************************************************************
      SUBROUTINE BBSMTB(Y,Y0,X0)

      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-50.0,Y+YDEL,X0-50.0,Y0)
      CALL PLINE(X0,Y+YDEL,X0,Y0)
      CALL PLINE(X0+50.0,Y+YDEL,X0+50.0,Y0)
      CALL PLINE(X0-50.0,Y+YDEL,X0+50.0,Y+YDEL)
      CALL PLINE(X0-50.0,Y0-15.0,X0+50.0,Y0-15.0)
      CALL PLINE(X0-50.0,Y0,X0+50.0,Y0)
      CALL TPRINT('Start',9.0,X0-25.0,Y0-10.0)
      CALL TPRINT('End',9.0,X0+25.0,Y0-10.0)
      CALL TPRINT('BETA ALPHA BETA UNITS',12.0,X0,Y0+10.0)
      RETURN
      END
C**************************************************************************
      SUBROUTINE PSSMTB(Y,Y0,X0)

      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-50.0,Y+YDEL,X0-50.0,Y0)
      CALL PLINE(X0,Y+YDEL,X0,Y0)
      CALL PLINE(X0+50.0,Y+YDEL,X0+50.0,Y0)
      CALL PLINE(X0-50.0,Y+YDEL,X0+50.0,Y+YDEL)
      CALL PLINE(X0-50.0,Y0-15.0,X0+50.0,Y0-15.0)
      CALL PLINE(X0-50.0,Y0,X0+50.0,Y0)
      CALL TPRINT('Start',9.0,X0-25.0,Y0-10.0)
      CALL TPRINT('End',9.0,X0+25.0,Y0-10.0)
      CALL TPRINT('PSI LOOPS',12.0,X0,Y0+10.0)
      RETURN
      END
C**************************************************************************
      SUBROUTINE DSSMTB(Y,Y0,X0)
      
      REAL X0,Y,Y0,YDEL
      DATA YDEL/8.0/

      CALL PLINE(X0-40.0,Y+YDEL,X0-40.0,Y0)
      CALL PLINE(X0+5.0,Y+YDEL,X0+5.0,Y0)
      CALL PLINE(X0+40.0,Y+YDEL,X0+40.0,Y0)
      CALL PLINE(X0-40.0,Y0-15.0,X0+40.0,Y0-15.0)
      CALL PLINE(X0-40.0,Y0,X0+40.0,Y0)
      CALL PLINE(X0-40.0,Y+YDEL,X0+40.0,Y+YDEL)
      CALL TPRINT('Cys 1',9.0,X0-20.0,Y0-10.0)
      CALL TPRINT('Cys 2',9.0,X0+20.0,Y0-10.0)
      CALL TPRINT('DISULPHIDES',12.0,X0,Y0+10.0)
      RETURN
      END
C**************************************************************************
C**************************************************************************
C
C  SUBROUTINE OPNDAT  -  Open the promotif.dat file and read through to
C                        the appropriate section
C
C----------------------------------------------------------------------+---

      SUBROUTINE OPNDAT(PDBCOD,CHAIN,PLTYPE,SPCIAL,PDBSUM_TEMPLATE,
CPDBSUM1 RAL 4/3/22 -->
C     -    DIRSEP,FILDAT,MSTART,IFAIL)
     -    DIRSEP,FILDAT,MSTART,PDBSM1,IFAIL)
CPDBSUM1 RAL 4/3/22 <--

      INCLUDE 'fcommon.inc'

      CHARACTER*1   CHAIN, DIRSEP
      CHARACTER*4   PDBCOD
      CHARACTER*20  PLTYPE, SPCIAL
      CHARACTER*512 FNAME, IREC, PDBSUM
      CHARACTER*(FNAMLN) PDBSUM_TEMPLATE(NPDBSUM_DIR)

      INTEGER       FILDAT, MSTART, LEN, LENSTR, LINE, NDONE

CPDBSUM1 RAL 4/3/22 -->
C      LOGICAL       GOTDAT, HAVFIL, IFAIL, MCSG, UPLOAD
      LOGICAL       GOTDAT, HAVFIL, IFAIL, MCSG, PDBSM1, UPLOAD
CPDBSUM1 RAL 4/3/22 <--

C---- Initialise variables
      IFAIL = .FALSE.
      GOTDAT = .FALSE.
      HAVFIL = .FALSE.
      LINE = 0
      NDONE = 0
      UPLOAD = .FALSE.
      MCSG = .FALSE.
      IF (SPCIAL.EQ.'UPLOAD') THEN
          UPLOAD = .TRUE.
      ELSE IF (SPCIAL.EQ.'MCSG') THEN
          MCSG = .TRUE.
      ENDIF

C---- Get length of search key
      LEN = LENSTR(PLTYPE)

C---- Find the name of the PDBsum directory
      CALL SUMNAM(PDBCOD,PDBSUM_TEMPLATE,PDBSUM,MCSG,UPLOAD)

C---- Form full name of the promotif.dat file
      FNAME = PDBSUM(1:LENSTR(PDBSUM)) // DIRSEP // 'promotif.dat'
CPDBSUM1 RAL 4/3/22 -->
      IF (PDBSM1) THEN
         FNAME = PDBSUM(1:LENSTR(PDBSUM)) // DIRSEP // 
     -        '/promotif/promotif.dat'
      ENDIF
CPDBSUM1 RAL 4/3/22 <--

C---- Open the promotif.dat file
      OPEN(UNIT=FILDAT, FILE=FNAME, STATUS='OLD',
     -        FORM='FORMATTED', ACCESS='SEQUENTIAL',
CVAX     -        READONLY,
     -        ERR=900)

C---- Read through the file until we get to the point we're after
 500  CONTINUE

C----     Read in the next line from the file
          READ(FILDAT,520,END=904,ERR=902) IREC
 520      FORMAT(A)

C----     If have reached the right section, determine whether this
C         is the right place to stop
          IF (GOTDAT) THEN

C----         If this is the next START record, then we've gone too
C             far and have to end here
              IF (IREC(1:6).EQ.'START:') GO TO 904

C----         If this is the right chain, increment count of records
C             of correct chain read in
              IF (IREC(1:1).EQ.CHAIN .AND. IREC(2:3).EQ.'::') THEN

C----             Increment counter
                  NDONE = NDONE + 1

C----             If this is one record prior to the start point, then
C                 end here
                  IF (NDONE.EQ.MSTART - 1) GO TO 800
              ENDIF

C----     Check whether line signifies the start of the data we are
C         looking for
          ELSE IF (IREC(1:6).EQ.'START:' .AND.
     -        IREC(7:LEN + 6).EQ.PLTYPE(1:LEN)) THEN

C----         Set flag that we have reached the right section
              GOTDAT = .TRUE.

C----         If starting at the first motif, then can stop reading
C             here
              IF (MSTART.EQ.1) GO TO 800
          ENDIF


C---- Loop back for next record in file
      GO TO 500

C---- Start-point for data reached
 800  CONTINUE

      GO TO 999

C---- Error messages
 900  CONTINUE
      PRINT*, '*** ERROR. Failed to locate promotif.dat file'
      PRINT*, '*** ', FNAME(1:LENSTR(FNAME))
      GO TO 990

 902  CONTINUE
      PRINT*, '*** ERROR. Read error on promotif.dat file'
      GO TO 990

 904  CONTINUE
      PRINT*, '*** ERROR. Required data not found in promotif.dat file'
      GO TO 990

 990  CONTINUE
      IFAIL = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETREC  -  Get next record from the promotif.dat file
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETREC(FILDAT,CHAIN,IREC,EOF)

      CHARACTER*1   CHAIN
      CHARACTER*3   KEY
      CHARACTER*(*) IREC

      INTEGER       FILDAT

      LOGICAL       EOF, HAVREC

C---- Initialise variables
      EOF = .FALSE.
      HAVREC = .FALSE.
      KEY = CHAIN // '::'

C---- Read through the file, storing all the tokens and their translations
 500  CONTINUE

C----     Read in the next line from the file
          READ(FILDAT,520,END=990,ERR=990) IREC
 520      FORMAT(A)

C----     Check whether have hit next set of data
          IF (IREC(1:6).EQ.'START:') GO TO 990

C----     Check whether this belongs to our chain of interest
          IF (IREC(1:3).EQ.KEY) HAVREC = .TRUE.

C---- If didn;t find a suitable record, loop back for next record in file
      IF (.NOT.HAVREC) GO TO 500
      GO TO 999

 990  CONTINUE
      EOF = .TRUE.

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GETPLT  -  Retrieve the plot data from the end of the line
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETPLT(IREC,JSEQNO,JSEQNOP,RSNAME,RSNAMP,SST,SSTP,HBP)

      CHARACTER*1   AACODE, SST(8), SSTP(5)
      CHARACTER*3   RSNAME(8), RSNAMP(5)
      CHARACTER*6   HBP(8,4), JSEQNO(8), JSEQNOP(5)
      CHARACTER*(*) IREC

      INTEGER       IHBOND, IPOS, IRES, J, K, LENSTR, NHBOND, NRES1,
     -              NRES2

C---- Initialise variables
      IPOS = 0
      DO 20 K = 1, 8
          JSEQNO(K) = ' '
          SST(K) = ' '
          DO 10 J = 1, 4
              HBP(K,J) = '000000'
 10       CONTINUE
 20   CONTINUE
      DO 40 K = 1, 5
          JSEQNOP(K) = ' '
          SSTP(K) = ' '
 40   CONTINUE

C---- Retrieve the number of residues in the 1, 2 strand
      READ(IREC,60) NRES1
 60   FORMAT(I1)
      IPOS = 2

C---- Get the residues and H-bond partners of the 1, 2 strand
      DO 500, IRES = 1, NRES1

C----     Read in this residue's details
          READ(IREC(IPOS:),80) JSEQNO(IRES), AACODE, SST(IRES), NHBOND
 80       FORMAT(A6,A1,A1,I1)
          IPOS = IPOS + 9

C----     Convert a.a. code to residue name
          CALL AA2NAM(AACODE,RSNAME(IRES))

C----     Loop over the residue's H-bond partners
          DO 400, IHBOND = 1, NHBOND

C----         Read in the H-bond partners
              READ(IREC(IPOS:),120) HBP(IRES,IHBOND), AACODE
 120          FORMAT(A6,A1)
              IPOS = IPOS + 7
 400      CONTINUE
 500  CONTINUE

C---- Read in number of residues in the X strand
      READ(IREC(IPOS:),60) NRES2
      IPOS = IPOS + 1

C---- Get the residues of the X strand
      DO 800, IRES = 1, NRES2

C----     Read in this residue's details
          READ(IREC(IPOS:),520) JSEQNOP(IRES), AACODE, SSTP(IRES)
 520      FORMAT(A6,A1,A1)
          IPOS = IPOS + 8

C----     Convert a.a. code to residue name
          CALL AA2NAM(AACODE,RSNAMP(IRES))
 800  CONTINUE

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
C  SUBROUTINE NAM2AA  -  Convert 3-letter residue name to single-letter
C                        a.a. code
C
C----------------------------------------------------------------------+---

      SUBROUTINE NAM2AA(RESNAM,AACODE)

      INTEGER       NAMINO

      PARAMETER    (NAMINO = 20)

      CHARACTER*1   AACODE
      CHARACTER*3   CODE(NAMINO), RESNAM
      CHARACTER*(NAMINO)  AMINOS

      INTEGER       IAMINO

      DATA  AMINOS / 'ACDEFGHIKLMNPQRSTVWY'/
      DATA  CODE   / 'ALA','CYS','ASP','GLU','PHE','GLY','HIS',
     -               'ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG',
     -               'SER','THR','VAL','TRP','TYR' /

C---- Initialise a.a. code
      AACODE = 'X'

C---- Get single-letter code for this residue
      DO 100, IAMINO = 1, NAMINO

C----     If correct residue, then store corresponding code
          IF (CODE(IAMINO).EQ.RESNAM) THEN
              AACODE = AMINOS(IAMINO:IAMINO)
              GO TO 999
          ENDIF
 100  CONTINUE

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
CPDBSUM RAL 15/3/05 <--
