C***********************************************************************
C*      NAME: psrout.f                                                 *
C*  FUNCTION: postscript subroutines for promotif                      *
C* COPYRIGHT: (1994-1996) University College London                    *
C*---------------------------------------------------------------------*
C*    AUTHOR: GAIL HUTCHINSON                                          *
C*      DATE: 5/7/94                                                         *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *

C**************************************************************
C     SUBROUTINE PSOPEN - Open postscript file and write out header
C
C***************************************************************
      SUBROUTINE PSOPEN(FNAME)
      
      INCLUDE 'psplot.par'

      INTEGER I,ICOL
      CHARACTER*(*) FNAME
      CHARACTER COLNO*2
c      OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED',
c     ~     ACCESS='SEQUENTIAL',CARRIAGECONTROL='LIST')
      OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED',
     ~     ACCESS='SEQUENTIAL')
      WRITE(1,5)
 5    FORMAT(
     ~     '%!PS-Adobe-3.0'/,
     ~     '%%Creator: Promotif',/,
CPDBSUM RAL 19/10/11 -->
     ~     '%%DocumentNeededResources: font Times-Roman Symbol',/,
     ~     '%%BoundingBox: (atend)',/,
     ~     '%%Pages: 1',/,
CPDBSUM RAL 19/10/11 <--
     ~     '%%EndComments',/,
     ~     '%%BeginProlog')
      WRITE(1,15)
 15   FORMAT(
     ~   '/tr {/Times-Roman findfont exch scalefont setfont } bind def',/,
     ~   '/sym {/Symbol findfont exch scalefont setfont} bind def',/,
     ~   '/Centre { ',/,
     ~   ' 2 copy dup tr exch stringwidth pop -2 div exch -3 div',
     ~   ' rmoveto tr show } bind def',/,
CPDBSUM RAL 5/4/05 -->
     ~   '/Bcentre { ',/,
     ~   '  dup /Times-Roman findfont exch scalefont setfont',/,
     ~   '  exch stringwidth pop -2 div exch -3 div rmoveto',/,
     ~   ' } bind def',/,
CPDBSUM RAL 5/4/05 <--
     ~   '/Scentre { ',/,
     ~   ' 2 copy dup sym exch stringwidth pop -2 div exch -3 div',
     ~   ' rmoveto sym show } bind def')
      WRITE(1,25)
 25   FORMAT(
     ~     '/Circle { gsave newpath 0 360 arc gsave ',
     ~     ' grestore stroke grestore }',/
     ~     ' bind def',/,
     ~     '/FCircle { gsave newpath 0 360 arc gsave fill',
     ~     ' grestore stroke grestore }',/
     ~     ' bind def',/,
     ~     '/HCircle { gsave newpath 0 180 arc gsave ',
     ~     ' grestore stroke grestore }',/
     ~     ' bind def',/,
     ~     '/FRectangle { gsave moveto lineto lineto lineto closepath',
     ~     ' fill stroke grestore} bind def',/, 
     ~     '/Rectangle { moveto lineto lineto lineto closepath',
     ~     ' stroke } bind def',/, 
     ~     '/Line { moveto lineto stroke } bind def')
      WRITE(1,35)
 35   FORMAT('/Col   { setrgbcolor } bind def')
CPDBSUM RAL 5/4/05 -->
      WRITE(1,30)
 30   FORMAT(
     -    '/Bprint { /Times-Bold findfont exch scalefont setfont show ',
     -    '} bind def',/,
     -    '/Lalign {',/,
     -    '  dup /Times-Bold findfont exch scalefont setfont',/,
     -    '  exch stringwidth pop pop -3 div 0.0 exch rmoveto',/,
     -    ' } bind def',/,
     -    '/Ralign {',/,
     -    '  dup /Times-Bold findfont exch scalefont setfont',/,
     -    '  exch stringwidth pop -1 mul exch -3 div rmoveto',/,
     -    ' } bind def')
CPDBSUM RAL 5/4/05 <--
      DO 40 ICOL=1,MAXCOL
         WRITE(COLNO,'(I2)') ICOL
         IF(COLNO(1:1).EQ.' ') COLNO(1:1)='0'
         WRITE(1,45) COLNO,(RGB(I,ICOL),I=1,3)
 45      FORMAT('/Col',A2,'{ ', 3F8.4,' setrgbcolor } def')
 40   CONTINUE
      WRITE(1,55)
 55   FORMAT(
     ~     '/arrowdict 14 dict def',/,'arrowdict begin',/,
     ~     '/mtrx matrix def',/,'end',/,
     ~     '/arrow',/,'{arrowdict begin',/,
     ~     '/headlength exch def',/,'/halfheadthickness exch 2 div def',
     ~     '/halfthickness exch 2 div def',/,
     ~     /,'/tipy exch def /tipx exch def',/,
     ~     '/taily exch def /tailx exch def',/,
     ~     '/dx tipx tailx sub def ',/,
     ~     '/dy tipy taily sub def ')
      WRITE(1,65)
 65   FORMAT(
     ~     '/arrowlength dx dx mul dy dy mul add sqrt def',/,
     ~     '/angle dy dx atan def',/,
     ~     '/base arrowlength headlength sub def',/,
     ~     '/savematrix mtrx currentmatrix def',/,
     ~     'tailx taily translate angle rotate',/,
     ~     '0 halfthickness neg moveto base halfthickness neg lineto',/,
     ~     'base halfheadthickness neg lineto arrowlength 0 lineto',/,
     ~     'base halfheadthickness lineto base halfthickness lineto',/,
     ~     '0 halfthickness lineto closepath',/,
     ~     'savematrix setmatrix',/,'end',/,' } def')
      WRITE(1,75)
 75   FORMAT(
     ~     '/vshowdict 4 dict def',/,
     ~     '/vshow',/,
     ~     '{vshowdict begin',/,
     ~     '/thestring exch def',/,
     ~     '/lineskip exch def',/,
     ~     'thestring',/,
     ~     '{',/,
     ~     '/charcode exch def',/,
     ~     '/thechar () dup 0 charcode put def',/,
     ~     '0 lineskip neg rmoveto',/,
     ~     'gsave',/,
     ~     'thechar stringwidth pop 2 div neg 0 rmoveto',/,
     ~     'thechar show',/,
     ~     'grestore',/,
     ~     '}forall',/,
     ~     'end',/,
     ~     '}def',/,     
     ~     '%%EndProlog')
CPDBSUM RAL 19/10/11 -->
      WRITE(1,78)
 78   FORMAT(
     ~     '%%BeginSetup',/,
     ~     '1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray',
     ~     ' [ ] 0 setdash newpath',/,
     ~     '%%EndSetup',/,
     ~     '%%Page: 1 1')
CPDBSUM RAL 19/10/11 <--
      WRITE(1,85)
 85   FORMAT(' 0.1 setlinewidth')
      RETURN
      END
C**************************************************************
C     SUBROUTINE PSCLOSE - Write out final lines to file and close
C
C**************************************************************
      SUBROUTINE PSCLOSE

      WRITE(1,5)
 5    FORMAT
CPDBSUM RAL 19/10/11 -->
C     ~     ('grestore stroke',/,
C     ~     'showpage',/,
C     ~     '%%EOF')
     ~     (
     ~     'showpage',/,
     ~     '%%Trailer',/,
     ~     '%%BoundingBox: 0 0 700 900',/,
     ~     '%%EOF')
CPDBSUM RAL 19/10/11 <--

      CLOSE(1)
      RETURN
      END
C**************************************************************
C     SUBROUTINE RAMACH - To plot out Ramachandran plot
C
C**************************************************************
CPDBSUM RAL 5/4/05 -->
C      SUBROUTINE RAMACH(XMIN,YMIN,XMAX,YMAX)
      SUBROUTINE RAMACH(XMIN,YMIN,XMAX,YMAX,ICOL)
CPDBSUM RAL 5/4/05 <--

      INCLUDE 'psplot.par'

      REAL XMIN,YMIN,XMAX,YMAX,XMID,YMID,FSIZE,FSIZE2,
     ~     XSIZ,SHADE
CPDBSUM RAL 5/4/05 -->
      INTEGER ICOL
      REAL DELTA
CPDBSUM RAL 5/4/05 <--

      DATA FSIZE/12.0/,FSIZE2/8.0/,SHADE/0.9/
      XMID=(XMIN+XMAX)/2.0
      YMID=(YMIN+YMAX)/2.0
      XSIZ=XMAX-XMIN
CPDBSUM RAL 5/4/05 -->
      DELTA = 0.04 * XSIZ
CPDBSUM RAL 5/4/05 <--

C     draw surrounding box
      IF(.NOT.BW) THEN
CPDBSUM RAL 5/4/05 -->
         CALL COLOUR(19)
         CALL FPSQUAR(XMID+DELTA,YMID+DELTA,XSIZ)
         CALL COLOUR(ICOL)
CPDBSUM RAL 5/4/05 <--
         CALL FPSQUAR(XMID,YMID,XSIZ)
      ELSE
         WRITE(1,100) SHADE
 100     FORMAT('gsave ',F3.1,' setgray')
         CALL FPSQUAR(XMID,YMID,XSIZ)
         WRITE(1,101)
 101     FORMAT('grestore')
      ENDIF

C     draw axes
      IF(.NOT.BW) CALL COLOUR(5)
      CALL PSQUAR(XMID,YMID,XSIZ)
      CALL PLINE(XMIN,YMID,XMAX,YMID)
      CALL PLINE(XMID,YMIN,XMID,YMAX)
      CALL SYMB('f',FSIZE,XMID,YMIN-FSIZE/2.0)
      CALL SYMB('y',FSIZE,XMIN-FSIZE/2.0,YMID)
      CALL TPRINT('-180',FSIZE2,XMIN,YMIN-FSIZE/2.0)
      CALL TPRINT('180',FSIZE2,XMAX,YMIN-FSIZE/2.0)
      CALL TPRINT('-180',FSIZE2,XMIN-FSIZE,YMIN)
      CALL TPRINT('180',FSIZE2,XMIN-FSIZE,YMAX)
      RETURN
      END
C***************************************************************
C     SUBROUTINE PLINE - Writes out postscript for drawing line
C
C***************************************************************
      SUBROUTINE PLINE(X1,Y1,X2,Y2)
      
      REAL X1,Y1,X2,Y2

      WRITE(1,5) X1,Y1,X2,Y2
 5    FORMAT(4(F5.1,1X),'Line')
      
      RETURN
      END
C***************************************************************
C     SUBROUTINE ARCT - Writes out postscript for drawing arc at top
C
C***************************************************************
      SUBROUTINE ARCT(X,Y,R)
      
      REAL X,Y,R

      WRITE(1,5) X,Y,R
 5    FORMAT('3 setlinewidth ',3(F5.1,1X),
     ~     '180 0 arcn stroke 1 setlinewidth')
      
      RETURN
      END
C***************************************************************
C     SUBROUTINE ARCB - Writes out postscript for drawing arc at bottom
C
C***************************************************************
      SUBROUTINE ARCB(X,Y,R)
      
      REAL X,Y,R

      WRITE(1,5) X,Y,R
 5    FORMAT('3 setlinewidth ',3(F5.1,1X),
     ~     '180 0 arc stroke 1 setlinewidth')
      
      RETURN
      END
C***************************************************************
C     SUBROUTINE TPRINT - Writes out text
C
C***************************************************************
      SUBROUTINE TPRINT(TEXT,SIZE,X,Y)
      REAL SIZE,X,Y
      CHARACTER*(*) TEXT

      WRITE(1,5) X,Y
 5    FORMAT(2(F5.1,1X),'moveto')
      WRITE(1,15) TEXT,SIZE
 15   FORMAT('(',A,')',1x,F4.1,' Centre')

      RETURN
      END
C***************************************************************
C     SUBROUTINE SYMB - Plots a symbol (e.g. phi,psi) on graph
C
C***************************************************************
      SUBROUTINE SYMB(LETTER,SIZE,X,Y)

      REAL SIZE,X,Y
      CHARACTER LETTER

      WRITE(1,5) X,Y
 5    FORMAT(2(F5.1,1X),'moveto')
      WRITE(1,15) LETTER, SIZE
 15   FORMAT('(',A1,')',1x,F4.1,' Scentre')

      RETURN
      END
c**************************************************************
c     SUBROUTINE CROSS - PLOTS A CROSS OF DEFINED SIZE
C
C**************************************************************
      SUBROUTINE CROSS(SIZE,X,Y)

      REAL SHADE,SIZE,X,Y
      
      DATA SHADE/0.5/

      WRITE(1,100) SHADE
 100  FORMAT('gsave ',F3.1,' setgray')

      WRITE(1,5) X,Y
 5    FORMAT(2(F5.1,1X),'moveto')
      WRITE(1,15) -SIZE/2.0, SIZE/2.0
 15   FORMAT(2(F5.1,1X),'rmoveto')
      WRITE(1,25) SIZE,-SIZE
 25   FORMAT(2(F5.1,1X),'rlineto stroke')
      WRITE(1,5) X,Y
      WRITE(1,15) -SIZE/2.0,-SIZE/2.0
      WRITE(1,25) SIZE,SIZE

      WRITE(1,101)
 101  FORMAT('grestore')

      RETURN
      END
C***************************************************************
C     SUBROUTINE PCIRC - Draws out a circle at given point
C
C***************************************************************
      SUBROUTINE PCIRC(X,Y,RADIUS)

      REAL X,Y,RADIUS

      WRITE(1,5) X,Y,RADIUS
 5    FORMAT(3(F5.1,1X),' Circle')

      RETURN
      END
C***************************************************************
C     SUBROUTINE PHCIRC - Draws out a circle at given point
C
C***************************************************************
      SUBROUTINE PHCIRC(X,Y,RADIUS)

      REAL X,Y,RADIUS

      WRITE(1,5) X,Y,RADIUS
 5    FORMAT(3(F5.1,1X),' HCircle')

      RETURN
      END
c***************************************************************
C     SUBROUTINE FPCIRC - Draws out a filled circle at given point
C
C***************************************************************
      SUBROUTINE FPCIRC(X,Y,RADIUS)

      REAL X,Y,RADIUS

      WRITE(1,5) X,Y,RADIUS
 5    FORMAT(3(F5.1,1X),' FCircle')

      RETURN
      END
c***************************************************************
C     SUBROUTINE PSQUAR - Draws a square centred on a given point
C
C****************************************************************
      SUBROUTINE PSQUAR(X,Y,SIZ)

      REAL X,Y,SIZ,XMIN,YMIN,XMAX,YMAX

      XMIN=X-SIZ/2.0
      YMIN=Y-SIZ/2.0
      XMAX=X+SIZ/2.0
      YMAX=Y+SIZ/2.0
      
      WRITE(1,5) XMIN,YMIN,XMIN,YMAX,XMAX,YMAX,XMAX,YMIN
 5    FORMAT(8(F5.1,1X),'Rectangle')

      RETURN
      END
C***************************************************************
C     SUBROUTINE FPSQUAR - Draws a filled square centred on a given point
C
C****************************************************************
      SUBROUTINE FPSQUAR(X,Y,SIZ)

      REAL X,Y,SIZ,XMIN,YMIN,XMAX,YMAX

      XMIN=X-SIZ/2.0
      YMIN=Y-SIZ/2.0
      XMAX=X+SIZ/2.0
      YMAX=Y+SIZ/2.0
      
      WRITE(1,5) XMIN,YMIN,XMIN,YMAX,XMAX,YMAX,XMAX,YMIN
 5    FORMAT(8(F5.1,1X),'FRectangle')

      RETURN
      END
C***************************************************************
C     SUBROUTINE PRECT - Draws a rectangle centred on a given point
C
C****************************************************************
      SUBROUTINE PRECT(X,Y,XSIZ,YSIZ)

      REAL X,Y,XSIZ,YSIZ,XMIN,YMIN,XMAX,YMAX

      XMIN=X-XSIZ/2.0
      YMIN=Y-YSIZ/2.0
      XMAX=X+XSIZ/2.0
      YMAX=Y+YSIZ/2.0
      
      WRITE(1,5) XMIN,YMIN,XMIN,YMAX,XMAX,YMAX,XMAX,YMIN
 5    FORMAT(8(F5.1,1X),'Rectangle')

      RETURN
      END
C***************************************************************
c     SUBROUTINE PBTURN - Draws schematic diagram for beta turn
C
C***************************************************************
      SUBROUTINE PBTURN(XMIN,YMIN,XSIZ,YSIZ,CRAD,SSIZ)
      
      INCLUDE 'psplot.par'

      REAL XMIN,YMIN,XSIZ,YSIZ,XMAX,YMAX,CRAD,SSIZ

      XMAX=XMIN+XSIZ
      YMAX=YMIN+YSIZ

      CALL PLINE(XMIN,YMIN,XMIN,YMAX)
      CALL PLINE(XMIN,YMAX,XMAX,YMAX)
      CALL PLINE(XMAX,YMAX,XMAX,YMIN)
CPDBSUM RAL 5/4/05 -->
C      IF(.NOT.BW) CALL COLOUR(11)
      IF(.NOT.BW) CALL COLOUR(7)
CPDBSUM RAL 5/4/05 <--
      CALL FPCIRC(XMIN,YMAX,CRAD)
CPDBSUM RAL 5/4/05 -->
      CALL COLOUR(4)
CPDBSUM RAL 5/4/05 <--
      CALL FPSQUAR(XMAX,YMAX,SSIZ)

      RETURN
      END
C**************************************************************
c     SUBROUTINE PGTURN - Draws schematic diagram for gamma turn
C
C***************************************************************
      SUBROUTINE PGTURN(XMIN,YMIN,XSIZ,YSIZ,CRAD)

      INCLUDE 'psplot.par'

      REAL XMIN,YMIN,XSIZ,YSIZ,XMAX,YMAX,XMID,CRAD

      XMAX=XMIN+XSIZ
      YMAX=YMIN+YSIZ
      XMID=(XMIN+XMAX)/2.0

      CALL PLINE(XMIN,YMIN,XMID,YMAX)
      CALL PLINE(XMID,YMAX,XMAX,YMIN)
      IF(.NOT.BW) CALL COLOUR(4)
      CALL FPCIRC(XMID,YMAX,CRAD)

      RETURN
      END
C**************************************************************
C     SUBROUTINE DASH - Turns on and off dashed lines
C
C***************************************************************
      SUBROUTINE DASH(ONOFF)

      INTEGER ONOFF

      IF(ONOFF.EQ.0) THEN
         WRITE(1,5)
 5       FORMAT(' [] 0 setdash')
      ELSE
         WRITE(1,15) ONOFF,ONOFF
 15      FORMAT(' [',I2,' ',I2,'] 0 setdash')
      ENDIF

      RETURN
      END
c***************************************************************
C     SUBROUTINE ARROW - Draws arrows
C
C****************************************************************
      SUBROUTINE ARROW(X1,Y1,X2,Y2,MINWID,MAXWID,HEADLN)
      
      REAL X1,Y1,X2,Y2,MINWID,MAXWID,HEADLN

      WRITE(1,5) X1,Y1,X2,Y2,MINWID,MAXWID,HEADLN
 5    FORMAT(7(F5.1,1X),'arrow stroke')
      
      RETURN
      END
C****************************************************************
C     SUBROUTINE FARROW - Draws filled arrows
C
C****************************************************************
      SUBROUTINE FARROW(X1,Y1,X2,Y2,MINWID,MAXWID,HEADLN)
      
      REAL X1,Y1,X2,Y2,MINWID,MAXWID,HEADLN

      WRITE(1,5) X1,Y1,X2,Y2,MINWID,MAXWID,HEADLN
 5    FORMAT(7(F5.1,1X),'arrow fill stroke')
      
      RETURN
      END
C****************************************************************
c     SUBROUTINE COLOUR - Calls colour routine
C
C****************************************************************
      SUBROUTINE COLOUR(IC)
      
      INTEGER IC

      IF(IC.LT.10) THEN
         WRITE(1,5) IC
 5       FORMAT('Col0',I1)
      ELSE
         WRITE(1,15) IC
 15      FORMAT('Col',I2)
      ENDIF

      RETURN
      END
C******************************************************************
C     SUBROUTINE VPRINT - to print text vertically
c
C*******************************************************************
      SUBROUTINE VPRINT(TEXT,X,Y)

      REAL X,Y,YSPACE
      CHARACTER TEXT*30

      DATA YSPACE/16.0/

      WRITE(1,5)
 5    FORMAT('/Helvetica findfont 16 scalefont setfont')
      WRITE(1,15) X,Y
 15   FORMAT(2(F5.1,1X),'moveto')
      WRITE(1,25) YSPACE,TEXT
 25   FORMAT(F5.1,1X,'(',A30,') vshow')

      RETURN
      END
C*******************************************************************
CPDBSUM RAL 5/4/05 -->
C**************************************************************************
C
C  SUBROUTINE BPRINT  -  Write out centred, bold text
C
C----------------------------------------------------------------------+---

      SUBROUTINE BPRINT(TEXT,SIZE,X,Y)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(1,100) X, Y
 100  FORMAT(2F7.2,' moveto')
      WRITE(1,200) TEXT, SIZE
 200  FORMAT('(',A,') ',F7.2,' Bcentre')
      WRITE(1,300) TEXT, SIZE
 300  FORMAT('(',A,') ',F7.2,' Bprint')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSRTXT  -  Write out right-justified text
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSRTXT(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(1,100) X, Y
 100  FORMAT(2F7.2,' moveto')
      WRITE(1,200) TEXT, SIZE
 200  FORMAT('(',A,') ',F7.2,' Ralign')
      WRITE(1,300) TEXT, SIZE
 300  FORMAT('(',A,') ',F7.2,' Bprint')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSLTXT  -  Write out left-justified text
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSLTXT(X,Y,SIZE,TEXT)

      CHARACTER*(*) TEXT
      REAL          SIZE, X, Y

      WRITE(1,100) X, Y
 100  FORMAT(2F7.2,' moveto')
      WRITE(1,200) TEXT, SIZE
 200  FORMAT('(',A,') ',F7.2,' Lalign')
      WRITE(1,300) TEXT, SIZE
 300  FORMAT('(',A,') ',F7.2,' Bprint')

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSLWID  -  Write line-width out to PostScript file
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSLWID(LWIDTH)

      REAL          LWIDTH

      WRITE(1,100) LWIDTH
 100  FORMAT(F5.2,' setlinewidth')

      RETURN
      END

C--------------------------------------------------------------------------
CPDBSUM RAL 5/4/05 <--



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

