      PROGRAM hera
C***********************************************************************
C*      NAME: p_hera3.f                                                *
C*  FUNCTION: Program to draw schematic diagrams of protein structures-*
C*            hydrogen bonding diagram                                 *
c             (modified version for use with PROMOTIF)                 *
C* COPYRIGHT: (1990-1996) University College London                    *
C----------------------------------------------------------------------+---
C*    AUTHOR: GAIL HUTCHINSON                                          *
C*      DATE: 1990                                                         *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C
C Compilation
C -----------
C
C f77 -o p_hera3 p_hera3.f
C
C or
C
C f77 -Wimplicit -fbounds-check -o p_hera3 p_hera3.f
C
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
c  25/09/91  GH     corrected direction of hydrogen bond arrows
c  26/09/91  GH     changed font size of text and automatically
c                   included insertion codes
c  10/09/91  GH     corrected ishape(2,2) to ishape(1,2) for extra
c                   integer file
c                   changed the interpretation of extra file data
c                   so that the first value in the file is attached to 
C                   the first residue selected to be plotted rather than 
C                   the first residue in the protein
c  02/12/91  GH     subroutine conect modified to cope better with
c                   bifurcated sheets e.g. in aspartic proteinases
c  30/03/92  GH     subroutine ldraw modified to cope with N terminii
c                   at very top of page
c  27/07/92  GH     correct plotting for large helical proteins -
c                   increased size of paginc for plotting helices only; 
C                   also increased array sizes to allow all helices from 
C                   1PRC to be plotted. this involved allowing for a 
C                   fifth page to be drawn
c  06/04/94  GH     declare all variables, increase sheet array sizes,
c                   change routines for fitting helices around sheets
c  12/09/94  GH     add extra option for choosing the starting orientation
c                   for helices where just helices are being plotted
c  14/09/94  GH     debugging for DSSP input files correcting number of 
C                   residues to be used
c  10/10/94  GH     remove variable LTCODE; correct NUMRES for DSSP files
c  17/11/94  GH     correct DSSP input code
C   8/05/12  RAL    Increase sizex from 200 to 400 todeal with very long
C                   chains
C  19/08/13  RAL    Introduced sizey to allow for very long helices
C  14/03/22  RAL    Fix to stop crashing when subscript out of bounds
C
c     **************************************************************

      include 'p_hera.par'
      integer access,arr1(lmax,7),arr2a(lmax,7),arr2b(lmax,7),
     ~     arr3(lmax,7),arr4(lmax,7),arr5(lmax,7),barr,barrel(msht),
     ~     bot(sizex,10),clmn(sizex,2),clmni(sizex,2),
     ~     conn(nsmax,msht*3),dir(nsmax,msht*3),eint(mres,3),finish,
     ~     endnfi,fistat,hbond(mres*4,2),hbp(mres,4),hbplot(mres*4),
     ~     hbpsn(mres*4,4),i,ia,iarr(sizex,2),
     ~     iarri(sizex,2),ibp(mres,2),ibuff(sizex,sizey),icol,icolp,
     ~     icolq,iconn(nsmax,3),idir(nsmax,3),iflip(msht),ifull,ihb,
     ~     ij,ins(mhel),iocode,option,iplot,
     ~     iplt1,iplt2,iseqno(mres),ishape(3,2),ist,ien
      integer ishord(msht),isht,istand,istph(mhel,2),
     ~     istpl(lmax,2),istps(mstr,7),
     ~     j,jbot(sizex,3),jflip(msht),jj,jjoin(3*nsmax,3),
     ~     jtop(sizex,3),jseqno(mres),jx,k,kloop,ks,lp(lmax,4),m,
     ~     maxlen,maxres,mcol,minres,n,n1,n2a,n2b,n3,n4,n5,
     ~     narray(sizex,sizey),ncol,nhrow,ninp,nj,nlrow,nr,nres,nrow,
     ~     nsrow,numres,numsht,numstr(msht),ok,opnera,paginc(npage,3),
     ~     paglim,posn(mres,2),shcode,sheet(300,6),shnum,smax,smin,
     ~     sstchk,str(nsmax,2),top(sizex,10),xcoord,ntop,
     ~     opnerb,nps,namlen,m1,m2,delxl,delyt,delyb,sht(mres)
      real efp(mres,3),energy,fshape(3,2),hbe(mres,4)
      character aa(mres),cif*2,cshape(3),echar(mres,3),
     ~     hlxchr*4,hxchr1*2,insert(mres),jsigni(sizex,2),
     ~     jsign(sizex,2),lsheet(msht),ofile*80,percnt,shape(mres),
     ~     shchr1*2,shchr2*2,shtref*62,space*1,sstr(mres),
     ~     title*4,chain,chn(mres),yesno
      character*3 resnam(mres)
      character*80 finam,pdbcod,pdbnam,sstnam,wkrec
      logical ierror,logics,ENSEMB,finuse
C     ***************
      parameter(opnera=1,ok=0,finish=2,space=' ',endnfi=6,opnerb=3)
      data energy/-0.5/,cif/' '/,hlxchr/'hHgG'/,hxchr1/'HG'/,
     ~     percnt/'%'/,shchr1/'eE'/,shchr2/'eE'/
      data(paginc(1,i),i=1,3)/0,0,0/,
     ~    (paginc(2,i),i=1,3)/29,43,43/,
     ~    (paginc(3,i),i=1,3)/56,84,84/,
     ~    (paginc(4,i),i=1,3)/84,123,123/,
     ~    (paginc(5,i),i=1,3)/124,165,165/,
     ~    (paginc(6,i),i=1,3)/153,197,197/

      data shtref(1:26)/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/,
     ~     shtref(27:52)/'abcdefghijklmnopqrstuvwxyz'/,
     ~     shtref(53:62)/'1234567890'/
      
c     check whether colour or black and white
CDEBUG
      FISTAT = OK
CDEBUG
      open(unit=10,file='promotif2.prm',err=1005)
      do 1001 i=1,8
         read(10,*)
 1001 continue
      read(10,'(a1)') yesno
      bw=logics(yesno)
      close(10)
      go to 1006
 1005 bw=.false.
 1006 continue

c     read in colours
      if(.not.bw) then
         open(unit=1,file='hera_colours')
         do 1010 icol=1,maxcol
            read(1,*) (rgb(i,icol),i=1,3),colnam(icol)
 1010    continue
         close(1)
      endif

C     ***************     
      do 2 i=1,msht
         ishord(i)=i
 2    continue
      do 3 i=1,100
         ins(i)=0
 3    continue
c     READ IN OPTIONS
      call rdopt(option,iplot,kloop,iplt1,iplt2,shcode,istand,
     ~     cshape,ishape,fshape,energy,minres,maxres,ntop,chain)
c      if(option.eq.1) call shtopt(smin,smax,shnum,ishord,jflip,iflip)
c     CHANGE FOR ENSEMBLE VERSION 21.5.96
      READ(*,'(A)') PDBNAM
      IF(PDBNAM(1:1).EQ.'%') THEN
         ENSEMB=.TRUE.
      ELSE
         ENSEMB=.FALSE.
      ENDIF
      IF(ENSEMB) THEN
         FINAM=PDBNAM
         FINUSE=.TRUE.
         NAMLEN=INDEX(PDBNAM,SPACE)-1
         DO 180 I=1,NAMLEN-1
            FINAM(I:I)=FINAM(I+1:I+1)
 180     CONTINUE
         FINAM(NAMLEN:NAMLEN)=SPACE
         OPEN(NAMFI,FILE=FINAM,STATUS='OLD',IOSTAT=IOCODE)
      ENDIF
 
 2000 CONTINUE
      IF(FISTAT.EQ.OK) THEN
         IF(ENSEMB) THEN
            READ(NAMFI,'(A)',IOSTAT=IOCODE,END=501) PDBNAM
            NAMLEN=INDEX(PDBNAM,SPACE)-1
            IF(IOCODE.LT.0) THEN
               CLOSE(NAMFI)
               FINUSE=.FALSE.
               FISTAT=ENDNFI
            ELSE
               WRITE(*,*) PDBNAM(1:NAMLEN)
            ENDIF
         ENDIF
      ENDIF

      IF(FISTAT.EQ.OK) THEN
         CALL GETNAM(PDBNAM,IST,IEN,IERROR)
         PDBCOD=PDBNAM(IST:IEN)
         NPS=INDEX(PDBCOD,' ')
         SSTNAM=PDBCOD(1:NPS-1)//'.sst'
         OPEN(UNIT=7,FILE=SSTNAM,STATUS='OLD',IOSTAT=IOCODE)
         IF(IOCODE.NE.0) FISTAT=OPNERB
      ENDIF

C     if .sst file exists then process further
      IF(FISTAT.EQ.OK) THEN
c     OPEN OUTPUT FILE
      ofile=pdbcod(1:nps-1)//'.her'
      open(unit=hbfile,file=ofile)
      call psinit(hbfile)

c     ZERO ARRAYS FOR EACH NEW PROTEIN
      DO 2050 I=1,SIZEX
         DO 2040 J=1,sizey
            NARRAY(I,J)=0
 2040    CONTINUE
         DO 2045 J=1,2
            CLMN(I,J)=0
            CLMNI(I,J)=0
            IARR(I,J)=0
 2045    CONTINUE
 2050 CONTINUE
      DO 2100 I=1,MRES*4
         DO 2070 J=1,4
            HBPSN(I,J)=0
 2070    CONTINUE
         DO 2080 J=1,2
            HBOND(I,J)=0
 2080    CONTINUE
 2100 CONTINUE
      DO 2200 I=1,MRES
         DO 2180 J=1,2
            POSN(I,J)=0
 2180    CONTINUE
 2200 CONTINUE
      DO 2300 I=1,LMAX
         DO 2280 J=1,2
            ISTPL(I,J)=0
 2280    CONTINUE
         DO 2290 J=1,4
            LP(I,J)=0
 2290    CONTINUE
 2300 CONTINUE
      DO 2400 I=1,MHEL
         INS(I)=0
         DO 2380 J=1,2
            ISTPH(I,J)=00
 2380    CONTINUE
 2400 CONTINUE
      IFULL=0
      IHB=0
      MAXLEN=0
      N1=0
      N2A=0
      N2B=0
      N3=0
      N4=0
      N5=0
      NHROW=1
      NINP=1
      JX=0
      M1=0
      M2=0
      DELXL=0
      DELYT=0
      DELYB=0
      SHNUM=0
c     
c     READ IN DATA FROM SECONDARY STRUCTURE FILE
         call readst(numres,iseqno,insert,jseqno,aa,resnam,sstr,sht,ibp,
     ~        hbp,hbe,nres,title,minres,maxres,option,access,chn,
     ~     chain,sstchk)
         if(sstchk.eq.1) go to 499
c     NOW READ IN EXTRA DATA FROM ANOTHER FILE IF REQUIRED
         if(option.eq.1) then
            if(shcode.gt.2.or.iplt1.eq.4.or.iplt2.eq.2)
     ~           call newdat(eint,efp,echar,cif,icolp,icolq,iplt1,
     ~           iplt2,maxres,ninp,minres)
           call shapes(shcode,istand,minres,maxres,sstr,cshape,shape,
     ~           hlxchr,shchr2,eint,efp,echar,ishape,fshape,ninp)
        endif
         smax=0
         smin=10
         if((iplot.ne.2).and.option.ne.3) then
c     FOR CALCULATIONS INVOLVING SHEETS AND HELICES
            do 110 m=minres,maxres
               if(index(shchr1,sstr(m)).ne.0) then
                  ia = sht(m)
                  if(smax.lt.ia.and.ia.ne.0) smax=ia
                  if(smin.gt.ia.and.ia.ne.0) smin=ia
                  go to 110
               else
                  if(index(hxchr1,sstr(m)).ne.0) 
     ~                call hlxrec(sstr,istph,m,nhrow,
     ~                 minres,maxres,maxlen)
               endif
 110        continue
            nhrow=nhrow-1
         else 
c     CASE WHERE ONLY WANT TO DEAL WITH SHEETS
            do 120 m=minres,maxres
               if(index(shchr1,sstr(m)).ne.0) then
                  ia = sht(m)
                  if(smax.lt.ia.and.ia.ne.0) smax=ia
                  if(smin.gt.ia.and.ia.ne.0) smin=ia
               endif
 120        continue
         endif
C     
C*******************************************************************
C     
C     FIRST DEAL WITH SHEETS
C     
C********************************************************************
         if(((option.eq.1.and.iplot.ne.3).or.option.eq.3).and.
     ~        smax.ne.0) then
            do 130 i=1,mstr
               do 125 j=1,7
                  istps(i,j)=0
 125           continue
 130        continue
            isht=0
            nsrow=0
            numsht=0
            do 131 i=1,msht
               barrel(i)=0
 131        continue
            do 133 i=1,nsmax
               do 132 j=1,msht*3
                  conn(i,j)=0
                  dir(i,j)=0
 132           continue
 133        continue
c            if(option.eq.3) shnum=smax
c     ALSO MODIFIED FOR MULTIPLE PROTEIN VERSION
            shnum=smax
            do 2600 ij=1,shnum
               ishord(ij)=ij
               jflip(i)=0
               iflip(i)=0
 2600       continue
c     END OF CHANGE
            do 200 ij=1,shnum
               jj=ishord(ij)
               call adist(ibp,sstr,jj,sht,sheet,ks,nr,str,nrow,
     ~              minres,maxres,posn,istps,
     ~              nsrow,shtref,numstr,numsht,lsheet)
               if(nrow.gt.1.and.nr.gt.2) then
                  call combin(sheet,str,iconn,ks,nr,jjoin,barr)
                  call shtins(iconn,idir,sheet,ibuff,str,
     ~                 ncol,ks,nr,posn,clmni,iarri,
     ~                 hbp,hbe,jsigni,jjoin,jtop,jbot)
                  if(option.eq.1) call orient(narray,ibuff,clmn,
     ~                 clmni,iarr,iarri,posn,ncol,jx,jsigni,jsign,
     ~                 iflip(ij),jflip(ij),iconn,idir,isht,jtop,
     ~                 top,jbot,bot,ifull)
                  if(ifull.eq.1) go to 501
                  do 170 j=1,3
                     do 160 i=1,nsmax
                        conn(i,isht+j)=iconn(i,j)
                        dir(i,isht+j)=idir(i,j)
 160                 continue
 170              continue
c     RECORD BARRELS
                  if(barr.eq.1) barrel(ij)=1
                  isht=isht+3
               else
                  write(*,171) shtref(jj:jj),ij
 171              format('Sheet ',a1,' consists of 2 residues ',
     ~                 '- will not be plotted ',i3)
               endif
 200        continue
         endif
C     
C*******************************************************************
C     
C     NOW DEAL WITH HELICES
C     
C********************************************************************
         if(iplot.eq.3.or.(iplot.eq.1.and.smax.eq.0)) then
            do 320 nj=1,nhrow
c     JUST DRAWS HELICES
               call dpnhel(clmn,istph,narray,posn,nj,numres,iseqno,iarr,
     ~              top,bot,jx,ntop)
 320        continue
         else
c     FIT HELICES AROUND SHEETS
            if(iplot.eq.1.and.nhrow.ge.1) 
     ~           call psnhel(clmn,istph,narray,posn,nhrow,iseqno,iarr,
     ~           ins,jsign,top,bot,jx)
         endif
         if(ifull.eq.1) go to 501
         do 330 i=1,sizex
            if(iarr(i,1).ne.0) mcol=i
 330     continue
         
C*****************************************************************
C     
C     NOW SORT OUT LOOPS - DEFINED AS ANYTHING NOT ASSIGNED A
C     POSITION ALREADY - AND HYDROGEN BONDS
C     
C*******************************************************************
         if(option.eq.1) then
c            if(nhrow.ne.0.or.shnum.gt.0.or.nsheet.gt.0) then
            if(nhrow.ne.0.or.smax.gt.0) then
            if(kloop.eq.1) then
               call dpsnlp(posn,lp,nlrow,istpl,minres,maxres)
               call plotlp(istpl,lp,nlrow,n1,n2a,n2b,n3,n4,n5,arr1,
     ~              arr2a,arr2b,arr3,arr4,arr5,clmn,mcol,top,bot)
            endif
c     OPEN POSTSCRIPT OUTPUT FILE AND WRITE PLOT DATA
            do 3351 k=scfile,scfil5
                open(unit=k,status='scratch')
 3351       continue
C     DRAW RESIDUES IN SHEETS AND HELICES
            call plotar(narray,aa,iplt1,jseqno,title,iseqno,iplot,
     ~           iplt2,shape,maxlen,cif,eint,efp,echar,icolp,
     ~           icolq,paginc,insert,resnam,chn,chain,numres,minres,
     -           maxres)
C     DRAW LOOPS
            if(kloop.eq.1) call ldraw(arr1,arr2a,arr2b,arr3,arr4,arr5,
     ~           n1,n2a,n2b,n3,n4,n5,narray,iarr,mcol,lp,clmn,nlrow,
     ~           maxlen,iplot,paginc,icol,m1,m2,delxl,delyt,delyb)
C     CALCULATE AND DRAW HYDROGEN BONDS
            call hbtran(hbond,energy,hbe,hbp,hbpsn,ihb,posn,minres,
     ~           maxres)
            call plothb(ihb,hbpsn,hbond,narray,iplt1,jseqno,
     ~           hbplot,clmn,iplot,maxlen,paginc,iseqno)
C     WRITE DATA TO OUTPUT FILE
            do 400 i=scfile,scfil5
               paglim=paginc(i-12,iplot)-1
               if(jx.ge.paglim.or.icol.ge.xcoord(paglim))
     ~              then
                  rewind(unit=i)
 350              read(i,'(A)',end=370) wkrec
                  n=index(wkrec,percnt)
                  write(hbfile,*) wkrec(1:n-1)
                  go to 350
 370              write(hbfile,340)
 340              format(' stroke showpage ')
                  if(iplot.eq.1.or.(iplot.eq.3.and.maxlen.gt.35)) then 
                     if(i.eq.13) then
                        write(hbfile,380)
                     else if(i.eq.14) then 
                        write(hbfile,382)
                     else
                        write(hbfile,384)
                     endif
                  endif
 380              format(' -520 0 translate')
 382              format(' -1020 0 translate')
 384              format(' -1560 0 translate')
               endif
                  close(unit=i,status='delete')
 400        continue
            write(hbfile,401)
 401        format(' %%Trailer')
            close(unit=hbfile)
         else
            write(*,*)' NO SECONDARY STRUCTURE FOUND FOR THIS PROTEIN'
            write(*,*)' NO DIAGRAM WILL BE DRAWN'
         endif
      endif
      else
         write(*,*) 'NO secondary structure file found'
         write(*,*) 'File: ',SSTNAM
      endif
 499  IF(ENSEMB) THEN
         IF(FISTAT.NE.FINISH) THEN
            FISTAT=OK
            IF(FINUSE) GO TO 2000
         ENDIF
      ENDIF
 501  stop   
      end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     
c     CALCULATION ROUTINES - INPUT
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine rdopt(option,iplot,kloop,iplt1,iplt2,shcode,
     ~     istand,cshape,ishape,fshape,energy,minres,maxres,ntop,
     ~     chain)

      integer iplot,kloop,iplt1,iplt2,shcode,istand,ishape(3,2),
     ~     minres,maxres,option,ntop
      real fshape(3,2),energy
      character cshape(3),chain
      chain=' '

c     SUBROUTINE TO TAKE IN INPUT OPTIONS FOR PROGRAM
C     
      read(*,*) option
      if(option.eq.1) then 
         call hbopt(iplot,kloop,iplt1,iplt2,shcode,istand,
     ~        cshape,ishape,fshape,energy,minres,
     ~        maxres,ntop,chain)
      endif
      return
      end
c*********************************************************************
      subroutine hbopt(iplot,kloop,iplt1,iplt2,shcode,istand,
     ~     cshape,ishape,fshape,energy,minres,maxres,ntop,chain)

C     READ IN OPTIONS FOR HYDROGEN BONDING DIAGRAM
      integer iplot,kloop,iplt1,iplt2,shcode,istand,nshape,
     ~     ishape(3,2),kshape,minres,maxres,num,opteng,ntop
      real fshape(3,2),energy
      character cshape(3),chain,ans

c      write(*,5)
c 5    format(' DO YOU WANT TO DRAW DIAGRAM FOR:'/
c     ~     '    HELICES AND SHEETS                      (1)'/
c     ~     '    JUST THE SHEETS                         (2)'/
c     ~     '    OR JUST THE HELICES                     (3)'/)
      read(*,'(i1)') iplot
c      write(*,15)
c 15   format(' PLOT LINES REPRESENTING LOOPS              (1)'/
c     ~       '               OR NOT                       (0)'/)
      read(*,'(i1)') kloop
c     choose what information to plot on the diagram
c      write(*,20)
c 20   format(' PLOT CONSECUTIVE NUMBERS                   (1)'/
c     ~       '      BROOKHAVEN SEQUENCE NUMBERS           (2)'/
c     ~       '      ACCESSIBILITIES                       (3)'/
c     ~       '      ADDED DATA                            (4)'/
c     ~       '      OR NONE OF THESE                      (5)'/)
      read(*,'(i2)') iplt1
c      write(*,30)
c 30   format(' PLOT ONE LETTER CODE                       (1)'/
c     ~       '      ADDED DATA                            (2)'/
c     ~       '   OR NONE OF THESE                         (3)'/)
      read(*,'(i2)') iplt2
c      write(*,35)
c 35   format(' CHOOSE ONE OF THE FOLLOWING TO CODE FOR'/
c     ~       '  THE SHAPE OF THE BOX:'/
c     ~       '      SECONDARY STRUCTURE                   (1)'/
c     ~       '      DSSP ACCESSIBILITY                    (2)'/
c     ~       '      INTEGER DATA FROM ANOTHER FILE        (3)'/
c     ~       '      CHARACTER DATA FROM ANOTHER FILE      (4)'/
c     ~       '      FLOATING POINT DATA FROM ANOTHER FILE (5)'/)
      read(*,'(i2)') shcode
      if(shcode.eq.1) then
c         write(*,40) 
c 40      format(' DO YOU WANT STANDARD CODING (box for sheet and'/
c     ~        '  circle for helix) (1 FOR YES, 0 FOR NO)'/)
         read(*,'(i2)') istand
         if(istand.eq.1) go to 200
      endif
c      write(*,50)
c 50   format(' INPUT NUMBER OF DIFFERENT REPRESENTATIONS'/
c     ~     '  (UP TO 2)'/)
      read(*,'(I2)') nshape
      num=0
      if(shcode.eq.1.or.shcode.eq.4) then
c         WRITE(*,55)
c 55      format(' TYPE THE CHARACTER TO BE REPRESENTED BY A CIRCLE'/)
         read(*,'(a1)') cshape(1)
         if(cshape(1).ne.' ') num=num+1
         if(num.eq.nshape) go to 190
c         write(*,60)
c 60      format(' TYPE THE CHARACTER TO BE REPRESENTED BY A SQUARE'/)
         read(*,'(a1)') cshape(2)
      else if(shcode.eq.2.or.shcode.eq.3) then
c         write(*,65)
c 65      format(' DO YOU WANT TO CODE FOR EACH SHAPE WITH:'/
c     ~        '  A SINGLE VALUE       (1)'/
c     ~        '  OR A RANGE OF VALUES (2)'/)
         read(*,*) kshape
         if(kshape.eq.1) then
c            write(*,70)
c 70         format(' TYPE THE INTEGER TO BE REPRESENTED BY A CIRCLE'/
c     ~           '  (FORMAT I4)')
            read(*,*) ishape(1,1)
            ishape(1,2)=ishape(1,1)
            if(ishape(1,1).ne.1000) num=num+1
            if(num.eq.nshape) go to 190
c            write(*,80)
c 80         format(' TYPE THE INTEGER TO BE REPRESENTED BY A SQUARE'/
c     ~           '  (FORMAT I4)'/)
            read(*,*) ishape(2,1)
            ishape(2,2)=ishape(2,1)
            if(ishape(2,1).ne.1000) num=num+1
         else
c            write(*,85)
c 85         format(' TYPE THE INTEGER RANGE TO BE REPRESENTED '/
c     ~           '  BY A CIRCLE'/)
            read(*,*) ishape(1,1),ishape(1,2)
            if(ishape(1,1).ne.1000) num=num+1
            if(num.eq.nshape) go to 190
c            write(*,95)
c 95         format(' TYPE THE INTEGER RANGE TO BE REPRESENTED '/
c     ~           '  BY A SQUARE'/)
            read(*,*) ishape(2,1),ishape(2,2)
         endif
      else
c         write(*,65)
         read(*,*) kshape
         if(kshape.eq.1) then
c            write(*,100)
c 100        format(' TYPE THE VALUE TO BE REPRESENTED BY A CIRCLE'/
c     ~           '  (FORMAT F5.1)'/)
            read(*,105) fshape(1,1)
 105        format(f5.1)
            fshape(1,2)=fshape(1,1)
            if(fshape(1,1).ne.1000.0) num=num+1
            if(num.eq.nshape) go to 190
c            write(*,110)
c 110        format(' TYPE THE VALUE TO BE REPRESENTED BY A SQUARE'/
c     ~           '  (FORMAT F6.1)'/)
            read(*,105) fshape(2,1)
            fshape(2,2)=fshape(2,1)
         else
c            write(*,115)
c 115        format(' TYPE THE RANGE TO BE REPRESENTED BY A CIRCLE'/
c     ~           '  (FORMAT F5.1,1X,F5.1)'/)
            read(*,120) fshape(1,1),fshape(1,2)
 120        format(f5.1,1x,f5.1)
            if(fshape(1,1).ne.1000.0) num=num+1
            if(num.eq.nshape) go to 190
c            write(*,125)
c 125        format(' TYPE THE RANGE TO BE REPRESENTED BY A SQUARE'/
c     ~           '  (FORMAT F5.1,1X,F5.1)'/)
            read(*,120) fshape(2,1),fshape(2,2)
         endif
      endif
 190  continue
 200  continue
c      write(*,405)
c 405  format
c     ~ (' THE DEFAULT THRESHOLD ENERGY FOR PLOTTING HYDROGEN BONDS'/
c     ~ ' IS -0.5. DO YOU WANT TO CHANGE THIS? (0 FOR NO,1 FOR YES)'/)
      read(*,*) opteng
      if(opteng.eq.1) then
c         write(*,4051)
c 4051    format(' GIVE VALUE FOR THRESHOLD ENERGY')
         read(*,406) energy
      endif
 406  format(f4.2)
c      write(*,415)
c 415  format(' INPUT INITIAL RESIDUE NUMBER (DEF 0)')
c      write(*,4151)
c 4151 format(' PLEASE GIVE SEQUENTIAL NUMBER, NOT BROOKHAVEN NUMBER'/)
      read(*,*) minres
c      write(*,420)
c 420  format(' INPUT FINAL RESIDUE NUMBER (DEF 0)')
c      write(*,4151)
      read(*,*) maxres
c      write(*,421)
c 421  format('Would you like to include just 1 chain instead (Y or N)')
      read(*,'(a1)') ans
      if(ans.eq.'Y'.or.ans.eq.'y') then
c         write(*,422)
c 422     format('Type in a chain letter')
         read(*,'(a1)') chain
      endif
      if(iplot.eq.3) then
c         write(*,430)
c 430     format(' WOULD YOU LIKE THE N TERMINUS OF THE PROTEIN TO BE:'/
c     ~        ' AT THE TOP OF THE PAGE       (TYPE 1)'/
c     ~        ' OR AT THE BOTTOM OF THE PAGE (TYPE 0)'/)
         read(*,*) ntop
      endif
      return
      end
c**********************************************************************
      subroutine shapes(shcode,istand,minres,maxres,sstr,cshape,shape,
     ~     hlxchr,shchr2,eint,efp,echar,ishape,fshape,ninp)

C     READ OPTIONS FOR SHAPES OF BOXES AROUND EACH RESIDUE (CIRCLES OR
C     RECTANGLES) FOR HYDROGEN BONDING DIAGRAM AND ASSIGN A SHAPE TO
C     EACH RESIDUE
      include 'p_hera.par'
      integer eint(mres,3),icol,ishape(3,2),
     ~     istand,maxres,minres,n,ninp,shcode
      real efp(mres,3),fshape(3,2)
      character sstr(mres),cshape(3),shape(mres),hlxchr*4,shchr2*2,
     ~     echar(mres,3)
      data icol/1/

      if(shcode.eq.1) then
         if(istand.eq.0) then
            do 5 n=minres,maxres
               if(sstr(n).eq.cshape(1).and.cshape(1).ne.' ') then
                  shape(n)='C'
               else if(sstr(n).eq.cshape(2).and.cshape(2).ne.' ') then
                  shape(n)='B'
               else
                  shape(n)=' '
               endif
 5          continue
         else
c            if(dssp.eq.0) then
               do 10 n=minres,maxres
                  if(index(hlxchr,sstr(n)).ne.0) then
                     shape(n)='C'
                  else
                     if(index(shchr2,sstr(n)).ne.0) shape(n)='B'
                  endif
 10            continue
c           else
c               if(dssp.eq.1) then
c                  do 11 n=minres+1,maxres
c                     if(index(hlxchr,sstr(n)).ne.0) then
c                        shape(n)='C'
c                     else
c                        if(index(shchr2,sstr(n)).ne.0) shape(n)='B'
c                     endif
c 11               continue
c               endif
c            endif
         endif
      else if (shcode.eq.2) then
c         do 15 n=minres,maxres
c            if(iacc(n).ge.ishape(1,1).and.iacc(n).le.ishape(1,2)) then
c               shape(n)='C'
c            else if(iacc(n).ge.ishape(2,1).and.iacc(n).le.ishape(2,2))
c    ~              then
c               shape(n)='B'
c            else 
c               shape(n)=' '
c            endif
c 15      continue
      else if(shcode.eq.3) then
         if(ninp.gt.1) then
c            write(*,20)
c 20         format(' WHICH OF THE INTEGERS IN THE EXTRA DATA FILE'/
c    ~           ' CONTAINS THE RELEVANT DATA FOR BOX SHAPE')
            read(*,*) icol
         endif
         do 30 n=minres,maxres
            if(eint(n,icol).ge.ishape(1,1).and.eint(n,icol).le.
     ~           ishape(1,2)) then
               shape(n)='C'
            else if(eint(n,icol).ge.ishape(2,1).and.eint(n,icol).le.
     ~              ishape(2,2)) then
               shape(n)='B'
            else
               shape(n)=' '
            endif
 30      continue
      else if(shcode.eq.4) then
         if(ninp.gt.1) then
c            write(*,35)
c 35         format(' WHICH OF THE CHARS IN THE EXTRA DATA FILE '/
c     ~           ' CONTAINS THE RELEVANT DATA FOR BOX SHAPE')
            read(*,*) icol
         endif
         do 40 n=minres,maxres
            if(echar(n,icol).eq.cshape(1).and.cshape(1).ne.' ') then
               shape(n)='C'
            else if(echar(n,icol).eq.cshape(2).and.cshape(2).ne.' ')
     ~              then
               shape(n)='B'
            else
               shape(n)=' '
            endif
 40      continue
      else
         if(ninp.gt.1) then
c            write(*,45)
c 45         format(' WHICH OF THE VALUES IN THE EXTRA DATA FILE '/
c     ~           ' CONTAINS THE RELEVANT DATA FOR BOX SHAPE')
            read(*,*) icol
         endif
         do 50 n=minres,maxres
            if(efp(n,icol).ge.fshape(1,1).and.efp(n,icol).le.
     ~           fshape(1,2)) then
               shape(n)='C'
            else if(efp(n,icol).ge.fshape(2,1).and.efp(n,icol).le.
     ~              fshape(2,2)) then
               shape(n)='B'
            else
               shape(n)=' '
            endif
 50      continue
      endif
      return
      end
c**********************************************************************
      subroutine shtopt(smin,smax,shnum,ishord,jflip,iflip)
      include 'p_hera.par'

      integer i,iflip(msht),ishord(msht),ishuf,ixshuf,iyshuf,
     ~     jflip(msht),jopt,shnum,smax,smin

c     RECORDS OPTIONS FOR SELECTING SHEETS AND FLIPPING THEM

c      write(*,5)
c 5    format(' PLOT ALL SHEETS             (0)'/
c     ~     ' OR SELECTED RANGE           (1)'/)
      read(*,*) jopt
      if(jopt.eq.1) then
c         write(*,*)' TYPE IN INITIAL AND FINAL SHEET
c     ~        (1 on each line)'
         read(*,*) smin
         read(*,*) smax
      endif
      shnum=smax-smin+1
C     OPTION FOR SHUFFLING SHEETS AROUND
c      write(*,*)' DO YOU WANT TO CHANGE THE ORDER OF THE SHEETS '
c      write(*,*)'  TYPE 1 IF YES, 0 IF NO'
      read(*,*) ishuf
c      if(ishuf.eq.1) then
c         write(*,*) ' TYPE IN THE ORDER OF SHEETS'
c         write(*,*)'  ( 1 NUMBER PER LINE)'
c         do 10 i=1,shnum
c            read(*,*) ishord(i)
c 10      continue
c      else
c         do 15 i=1,shnum
c            ishord(i)=i+smin-1
c 15      continue
c      endif
c      write(*,*) ' DO YOU WANT TO FLIP ANY SHEETS FROM LEFT TO RIGHT'
c      write(*,*)'  TYPE 1 IF YES, O IF NO'
      read(*,*) ixshuf
      if(ixshuf.eq.1) then
c         write(*,16)
c 16      format(' TYPE IN 1 FOR EACH SHEET TO BE FLIPPED'/
c     ~              '          0 FOR OTHER SHEETS')
         do 20 i=1,shnum
            read(*,*) jflip(i)
 20      continue
      else
         do 25 i=1,shnum
            jflip(i)=0
 25      continue
      endif
c      write(*,*) ' DO YOU WANT TO FLIP ANY SHEETS FROM TOP TO BOTTOM'
c      write(*,*) '  TYPE 1 IF YES, 0 IF NO'
      read(*,*) iyshuf
      if(iyshuf.eq.1) then
c         write(*,16)
         do 30 i=1,shnum
            read(*,*) iflip(i)
 30      continue
      else
         do 35 i=1,shnum
            iflip(i)=0
 35      continue
      endif
      return
      end
c**********************************************************************
      subroutine hbtran(hbond,energy,hbe,hbp,hbpsn,ihb,posn,
     ~     minres,maxres)
      include 'p_hera.par'
c     TRANSLATES THE INPUT HYDROGEN BOND DATA INTO DATA RELATING RESIDUES
c     AT DIFFERENT POSITIONS IN narray
      integer hbpsn(mres*4,4),posn(mres,2),hbp(mres,4),hbond(mres*4,2),
     ~     i,ihb,i1,i2,j,k,m,minres,maxres,n,np
      real hbe(mres,4),energy

      do 20 m=minres,maxres
         do 10 n=1,4
            if(hbe(m,n).le.energy) then
c               if(dssp.eq.0) then
                  np=hbp(m,n)
                  if(np.lt.m.and.np.ge.minres) go to 10
c               else
c                  np=hbp(m,n)+m
c               endif
               if(hbp(m,n).gt.0.or.np.lt.minres) then
                  ihb=ihb+1
                  if(n.eq.1.or.n.eq.3) then
                     hbond(ihb,1)=m
                     hbond(ihb,2)=np
                  else
                     hbond(ihb,1)=np
                     hbond(ihb,2)=m
                  endif
               endif
            endif
 10      continue
 20   continue

      do 30 j=1,4
         do 25 i=1,mres
            hbpsn(i,j)=0
 25      continue
 30   continue
c    
      do 40 k=1,ihb
         i1=hbond(k,1)
         i2=hbond(k,2)
         if(i1.gt.0) then
            hbpsn(k,1)=posn(i1,1)
            hbpsn(k,2)=posn(i1,2)
         else
            hbpsn(k,1)=0
            hbpsn(k,2)=0
         endif
         if(i2.gt.0) then
            hbpsn(k,3)=posn(i2,1)
            hbpsn(k,4)=posn(i2,2)
         else
            hbpsn(k,3)=0
            hbpsn(k,4)=0
         endif
 40   continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                 
      subroutine newdat(eint,efp,echar,cif,icolp,icolq,iplt1,iplt2,
     ~     maxres,ninp,minres)
      include 'p_hera.par'

C     READ IN DATA FROM ADDITIONAL FILES IF REQUIRED
      integer eint(mres,3),i,icolp,icolq,iplt1,iplt2,
     ~     j,ninp,minres,maxres
      real efp(mres,3)
      character cif*2,echar(mres,3),fmt*40,ifile*80
      data fmt/'(                                      )'/
C     
c      write(*,*)' GIVE THE NAME OF THE INPUT FILE FOR EXTRA DATA'
      read(*,'(A)') ifile
      open(unit=8,file=ifile,status='old',err=100)
c      write(*,2)
c 2    format(' CHOOSE THE TYPE OF INPUT DESIRED'/
c     1     '   CHARACTER ONLY                   (TYPE C)'/
c     1     '   INTEGER ONLY                 (TYPE I)'/
c     1     '   FLOATING POINT ONLY          (TYPE F)'/
c     1     '   CHARACTER AND INTEGER        (TYPE CI)'/
c     1     '   CHARACTER AND FLOATING POINT (TYPE CF)'/
c     1     '   INTEGER AND FLOATING POINT  (TYPE IF)')
      read(*,*) cif
      if(cif.eq.'C '.or.cif.eq.'I '.or.cif.eq.'F '.
     ~     or.cif.eq.'c '.or.cif.eq.'i '.or.cif.eq.'f ') then
c         write(*,*)' HOW MANY OF THESE DO YOU WANT TO READ IN(UP TO 3)?'
c         read(*,*) ninp
      endif
c      write(*,*)' GIVE INPUT FORMAT'
      read(*,'(a38)') fmt(2:39)
      
      if(cif.eq.'C '.or.cif.eq.'c ') then
         do 10 i=minres,maxres
            read(8,fmt,end=110)(echar(i,j),j=1,ninp)
 10      continue
      else if(cif.eq.'I '.or.cif.eq.'i ') then
         do 11 i=minres,maxres
            read(8,fmt,end=110)(eint(i,j),j=1,ninp)
 11      continue
      else if(cif.eq.'F '.or.cif.eq.'f ') then
         do 12 i=minres,maxres
            read(8,fmt,end=110)(efp(i,j),j=1,ninp)
 12      continue
      else if(cif.eq.'CI'.or.cif.eq.'ci') then
         if(index(fmt,'A').lt.index(fmt,'I').or.
     ~        index(fmt,'a').lt.index(fmt,'i')) then
            do 13 i=minres,maxres
               read(8,fmt,end=110)echar(i,1),eint(i,1)
 13         continue
         else 
            do 23 i=minres,maxres
               read(8,fmt,end=110) eint(i,1),echar(i,1)
 23         continue
         endif
      else if(cif.eq.'CF'.or.cif.eq.'cf') then
         if(index(fmt,'A').lt.index(fmt,'F').or.
     ~        index(fmt,'a').lt.index(fmt,'f')) then
            do 14 i=minres,maxres
               read(8,fmt,end=110)echar(i,1),efp(i,1)
 14         continue
         else 
            do 24 i=minres,maxres
               read(8,fmt,end=110) efp(i,1),echar(i,1)
 24         continue
         endif
      else
         if(cif.eq.'IF'.or.cif.eq.'if') then
            if(index(fmt,'i').lt.index(fmt,'f').or.
     ~           index(fmt,'I').lt.index(fmt,'F')) then
               do 15 i=minres,maxres
                  read(8,fmt,end=110)eint(i,1),efp(i,1)
 15            continue
            else
               do 25 i=minres,maxres
                  read(8,fmt,end=110)efp(i,1),eint(i,1)
 25            continue
            endif
         endif
      endif
      if(ninp.eq.1) then
         icolp=1
         icolq=1
      else
         if(iplt1.eq.4) then
            if(cif.eq.'I ') then
c               write(*,20)
c 20            format(' WHICH OF THE INTEGERS IN THE EXTRA DATA FILE'/
c     ~              ' CONTAINS THE RELEVANT DATA FOR DRAWING IN BOXES')
               read(*,*) icolp
            else
               if(cif.eq.'F ') then
c                  write(*,30)
c 30               format(' WHICH OF THE REALS IN THE EXTRA DATA FILE'/
c    ~               ' CONTAINS THE RELEVANT DATA FOR DRAWING IN BOXES')
                  read(*,*) icolp
               endif
            endif
         endif
         if(iplt2.eq.2) then
            if(cif.eq.'C ') then
c               write(*,35)
c 35            format(' WHICH OF THE CHARACTERS IN THE EXTRA DATA FILE'/
c     ~              ' CONTAINS THE RELEVANT DATA FOR DRAWING IN BOXES')
               read(*,*) icolq
            endif
         endif
      endif
 110  continue
 100  return
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     
c     CALCULATION ROUTINES - SHEET CALCULATIONS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine shtins(iconn,idir,sheet,ibuff,str,jxx
     ~     ,ks,nr,posn,clmni,iarri,hbp,hbe,jsigni,
     ~     jjoin,jtop,jbot)
      include 'p_hera.par'
C     ASSIGN RESIDUES IN THIS SHEET TO A POSITION IN A TEMPORARY
C     ARRAY (IBUFF)
      integer brp(40,2),clmni(sizex,2),col(nsmax),dst,hbp(mres,4),
     ~     hlen,i,iarri(sizex,2),ibuff(sizex,sizey),iconn(nsmax,3),
     ~     idir(nsmax,3),ie,ii,ik,ilen,iposn,iscen,iscol,
     ~     iydel,j,j1,j2,jbot(sizex,3),ji,jjoin(3*nsmax,3),jk,
     ~     jtop(sizex,3),jx,jxx,jy,jys,k,kk,ks,l,len,m,mindel,
     ~     mindst,mm,n,nbr,nr,posn(mres,2),rp(40),sgn(nsmax),sgni,sgnj,
     ~     sheet(300,6),str(nsmax,2),vac,y1(40),y2(40)
      real hbe(mres,4),x
      character jsigni(sizex,2)
      
      mindel=500
      do 5 i=1,sizex        
         do 3 j=1,sizey
            ibuff(i,j)=0
 3       continue
         do 4 j=1,2
            iarri(i,j)=0
            clmni(i,j)=0
 4       continue
         do 55 j=1,3
            jtop(i,j)=0
            jbot(i,j)=0
 55      continue
 5    continue
      do 8 j=1,2
         do 7 i=1,40
            brp(i,j)=0
 7       continue
 8    continue
      do 10 j=1,3
         do 9 i=1,nsmax
            idir(i,j)=0
 9       continue
 10   continue
      do 13 j=1,nsmax
         col(j)=0
 13   continue
C
      j=iconn(1,1)
      sgn(1)=1
      jx=4
      idir(1,1)=sgn(1)
C     INSERT FIRST STRAND IN CONNECTIVITY ARRAY
      jys=iscen(str(j,1),str(j,2))
      call sins(str(j,1),str(j,2),jx,jys,ibuff,posn,
     1     sgn(1),clmni,iarri,jsigni,jtop,jbot)
C     COL CONTAINS MARKER FOR THE COLUMN OF THE STRAND
      col(j)=jx
C     LOOP THROUGH REMAINDER OF THE STRANDS IN THE ARRAY
      do 100 i=1,ks
         nbr=0
         do 50 j=1,2
            do 49 k=1,40
               brp(k,j)=0
 49         continue
 50      continue
         j1=iconn(i,1)
         if(iconn(i+1,1).ne.0) then
            j2=iconn(i+1,1)
         else
            go to 101
         endif
C     COLLECT BRIDGES BETWEEN THESE TWO STRANDS
         do 15 m=1,nr
            if(sheet(m,4).eq.j1) then
               do 14 n=5,6
                  if(sheet(m,n).eq.j2) then
                     nbr=nbr+1
                     brp(nbr,1)=sheet(m,1)
                     brp(nbr,2)=sheet(m,n-3)
                  endif
 14            continue
            endif
 15      continue
C     FIND DIRECTION OF STRAND
         call direct(brp,sgn(i+1),sgn(i),str,hbe,hbp,j1,j2)
         if(sgn(i+1).eq.1) then
            do 30 m=1,nbr
               rp(m)=brp(m,2)-str(j2,1)+1
 30         continue
         else
            do 31 m=1,nbr
               rp(m)=str(j2,2)-brp(m,2)+1
 31         continue
         endif
C     POSITION OF 2ND STRAND CALCULATED TO MINIMIZE TOTAL BRIDGE LENGTHS
         len=str(j2,2)-str(j2,1)+1
         hlen=int(len/2)
         mindst=100
         do 33 m=1,nbr
            mm=brp(m,1)
            y1(m)=posn(mm,2)
 33      continue        
         do 40 l=1,sizey
            do 34 m=1,nbr
               y2(m)=rp(m)+l
 34         continue
            dst=0
            do 35 m=1,nbr
               dst=dst+(y2(m)-y1(m))**2
 35         continue
            if(dst.lt.mindst) then
               mindst=dst
               jys=l+1
            endif
 40      continue        
C     INSERT IT
         jx=jx+2
         idir(i+1,1)=sgn(i+1)
         call sins(str(j2,1),str(j2,2),jx,jys,ibuff,posn,
     1        sgn(i+1),clmni,iarri,jsigni,jtop,jbot)
         col(j2)=jx
 100  continue        
 101  continue
      jxx=jx
C     NOW INSERT EXTRA STRANDS FROM JJOIN ARRAY
      do 201 ik=1,3
         do 200 i=1,3*ks
            if(jjoin(i,ik).ne.0) then
c               j1=(i+1)/2
C*************************
c     NEW 17.11.94
               x=(i+1)/3.0
               j1=nint(x)
               j2=jjoin(i,ik)
C     END OF NEW BIT
C*************************
C     J1 IS ALREADY ATTACHED STRAND, J2 IS NEW STRAND
               nbr=0
               do 88 j=1,2
                  do 87 ii=1,40
                     brp(ii,j)=0
 87               continue
 88            continue
               
C     COLLECT BRIDGES
               do 105 m=1,nr
                  if(sheet(m,4).eq.j1) then
                     do 104 n=5,6
                        if(sheet(m,n).eq.j2) then
                           nbr=nbr+1
                           brp(nbr,1)=sheet(m,1)
                           brp(nbr,2)=sheet(m,n-3)
                        endif
 104                    continue
                     endif
 105              continue
C     CALCULATE DIRECTION OF EXISTING STRAND
                  ii=str(j1,1)
                  ie=str(j1,2)
        if(posn(ie,2).gt.posn(ii,2)) then
           sgni=1
        else
           sgni=-1
        endif
C     CALCULATE DIRECTION OF NEW STRAND
        call direct(brp,sgnj,sgni,str,hbe,hbp,j1,j2)
C     CALCULATE VERTICAL POSITION OF NEW STRAND RELATIVE TO OLD
        if(sgnj.eq.1) then
           do 110 m=1,nbr
              rp(m)=brp(m,2)-str(j2,1)+1
 110       continue
        else
           do 111 m=1,nbr
              rp(m)=str(j2,2)-brp(m,2)+1
 111       continue
        endif
        mindst=100
        do 113 m=1,nbr
           mm=brp(m,1)
           y1(m)=posn(mm,2)
 113    continue
        
        do 120 l=1,sizey
           do 114 m=1,nbr
              y2(m)=rp(m)+l
 114       continue
           dst=0
           do 115 m=1,nbr
              dst=dst+(y2(m)-y1(m))**2
 115       continue
           if(dst.lt.mindst) then
              mindst=dst
              jys=l+1
           endif
 120    continue
        
C     CALCULATE BEST COLUMN TO PUT STRAND IN + BEST POSITION
        ilen=str(j2,2)-str(j2,1)        
C     COLUMN CONTAINING EXISTING STRAND
        iscol=col(j1)
        jk=0
        jx=0
        jy=0
        do 130 k=0,ilen+5
           do 129 ji=iscol-2,iscol+2,4
              vac=0
            do 125 jy=jys-2+k,jys+ilen+k+1
CRAL 14/03/22 -->
C               if(jy.lt.sizey.and.jy.gt.0) then
               if(jy.lt.sizey.and.ji.gt.0.and.jy.gt.0) then
CRAL <-- 14/03/22 -->
                  if(ibuff(ji,jy).ne.0) vac=vac+1
               endif
 125        continue
            if(vac.eq.0) then
               iydel=k
               jx=ji
               go to 131
            endif
            vac=0
            do 127 jy=jys-2-k,jys+ilen-k+1
               if(jy.lt.sizey.and.jy.gt.0) then
                  if(ibuff(ji,jy).ne.0) vac=vac+1
               endif
 127        continue
            if(vac.eq.0) then
               iydel=-k
               jx=ji
               go to 131
            endif
 129     continue
 130  continue
 131  continue
c      IF(JX.EQ.0) THEN
C     STILL NEED TO FIND ANOTHER PATTERN

c
      call sins(str(j2,1),str(j2,2),jx,jys+iydel,ibuff,posn,
     ~     sgnj,clmni,iarri,jsigni,jtop,jbot)
      col(j2)=jx
c     ADD TO CONNECTIVITY ARRAY - SEARCH FOR EXISTING STRAND IN 'ICONN'
      do 140 kk=1,3
         do 139 k=1,ks
            if(iconn(k,kk).eq.j1) then
               iposn=k
               go to 141
            endif
 139     continue
 140  continue
 141  continue
      if(jx.eq.iscol-2) then
         if(iposn.gt.1) then
            if(iconn(iposn-1,2).eq.0) then
               iconn(iposn-1,2)=j2
               idir(iposn-1,2)=sgnj
            else
               iconn(iposn-1,3)=j2
               idir(iposn-1,3)=sgnj
            endif  
         else
c     CASE IF EXTEND BEYOND THE END OF THE ARRAY -SHIFT ALL STRANDS UP ONE
            do 150 k=2,1,-1
               do 149 j=nsmax,2,-1
                  iconn(j,k)=iconn(j-1,k)
                  idir(j,k)=idir(j-1,k)
 149           continue
               iconn(1,k)=0
               idir(i,k)=0
 150        continue
            if(iconn(1,2).eq.0) then
               iconn(1,2)=j2
               idir(1,2)=sgnj
            else
               iconn(1,3)=j2
               idir(1,3)=sgnj
            endif  
         endif                
      else 
         if(jx.eq.iscol+2) then
            if(iconn(iposn+1,2).eq.0) then
               iconn(iposn+1,2)=j2
               idir(iposn+1,2)=sgnj
            else
               iconn(iposn+1,3)=j2
               idir(iposn+1,3)=sgnj
            endif
         endif
      endif
      endif
 200  continue
 201  continue
      return
      end
C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c
      subroutine direct(brp,sgnj,sgni,str,hbe,hbp,j1,j2)
C     CALCULATE DIRECTION OF NEW STRAND
      include 'p_hera.par'
      integer brp(40,2),hb(20,2),hbp(mres,4),i,ihb,j,j1,j2,
     ~     m,n,np,sgni,sgnj,str(nsmax,2)
      real hbe(mres,4)

      ihb=0
      do 10 i=1,2
         do 9 j=1,20
            hb(j,i)=0
 9       continue
 10   continue
      if(brp(1,2).gt.brp(2,2).and.brp(2,2).ne.0) then
         sgnj=-sgni
      else if(brp(1,2).le.brp(2,2).and.brp(2,2).ne.0) then
         sgnj=sgni
      else
c         if(dssp.eq.1) then
c            do 20 n=str(j1,1)-1,str(j1,2)+1
c               do 19 m=1,4
c                  if(hbe(n,m).le.-0.5) then
c                     np=n+hbp(n,m)
c                     if(np.ge.str(j2,1)-1.and.np.le.str(j2,2)+1)
c     ~                    then
c                        ihb=ihb+1
c                        hb(ihb,1)=n
c                        hb(ihb,2)=np
c                     endif
c                  endif
c 19            continue
c 20         continue
c         else
            do 30 n=str(j1,1)-1,str(j1,2)+1
               do 29 m=1,4
                  if(hbe(n,m).le.-0.5.and.hbp(n,m).gt.0) then
                     np=hbp(n,m)
                     if(np.ge.str(j2,1)-1.and.np.le.str(j2,2)+1) then
                        ihb=ihb+1
                        hb(ihb,1)=n
                        hb(ihb,2)=np
                     endif
                  endif
 29            continue
 30         continue
c         endif
         if(hb(2,2).gt.hb(1,2).and.hb(2,1).gt.hb(1,1)) then
            sgnj=sgni
         else if(hb(2,2).lt.hb(1,2).and.hb(2,1).lt.hb(1,1)) then
            sgnj=sgni
         else
            sgnj=-sgni
         endif
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine sins(m,n,jx,jys,ibuff,posn,sgn,clmni,iarri,
     ~     jsigni,jtop,jbot)
C     INSERT STRAND INTO IBUFF
      include 'p_hera.par'

      integer clmni(sizex,2),iarri(sizex,2),ibuff(sizex,sizey),
     ~     i1,i2,j,jbot(sizex,3),jtop(sizex,3),jx,jy,jys,m,n,
     ~     posn(mres,2),sgn
      character jsigni(sizex,2)

      if(jx.gt.sizex.or.jx.le.0) then
         write(*,*)' no room left in output array'
         go to 100
      endif
      if(sgn.eq.1) then
         jsigni(jx,1)='a'
         jsigni(jx,2)='b'
         i1=min0(m,n)
         i2=max0(m,n)
         if(posn(i1,1).ne.0) i1=i1+1
         if(posn(i2,1).ne.0) i2=i2-1
      else
         jsigni(jx,1)='b'
         jsigni(jx,2)='a'
         i1=max0(m,n)
         i2=min0(m,n)
         if(posn(i1,1).ne.0) i1=i1-1
         if(posn(i2,1).ne.0) i2=i2+1
      endif
      jy=jys
      do 10 j=i1,i2,sgn
         if(posn(j,1).eq.0) then
            ibuff(jx,jy)=j
            posn(j,1)=jx
            posn(j,2)=jy
            jy=jy+1
         endif
 10   continue
      if(clmni(jx,2).eq.0) then
         clmni(jx,2)=i2
         iarri(jx,2)=posn(i2,2)
         jtop(jx,1)=i2
      else
         call ttop(jtop,jx,i2)
         j=clmni(jx,2)
         if(posn(j,2).lt.posn(i2,2)) then
            clmni(jx,2)=i2
            iarri(jx,2)=posn(i2,2)
         endif
      endif      
      if(clmni(jx,1).eq.0) then
         clmni(jx,1)=i1
         iarri(jx,1)=posn(i1,2)
         jbot(jx,1)=i1
      else
          call ttop(jbot,jx,i1)
          j=clmni(jx,1)
          if(posn(j,2).gt.posn(i1,2)) then
             clmni(jx,1)=i1
             iarri(jx,1)=posn(i1,2)
          endif
       endif
       
 100   return
       end        
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
       subroutine ttop(array,jx,i)
c       RECORD RESIDUES AT THE TOP AND BOTTOM OF EACH COLUMN
       include 'p_hera.par'
       integer array(sizex,3),jx,k,i
       do 5 k=1,3
          if(array(jx,k).eq.0) then
          array(jx,k)=i
          go to 6
       endif
 5    continue
 6    return
      end
c     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      subroutine orient(narray,ibuff,clmn,clmni,iarr,iarri,posn
     1     ,ncol,jx,jsigni,jsign,iflip,jflip,iconn,idir,isht
     1           ,jtop,top,jbot,bot,ifull)
C     CALCULATE ORIENTATION OF THE NEW SHEET RELATIVE TO ALL PREVIOUS 
C     SHEETS DRAWN AND WRITE FROM IBUFF ARRAY INTO FINAL ARRAY (NARRAY)
      include 'p_hera.par'

      integer bot(sizex,10),clmn(sizex,2),clmni(sizex,2),del(sizex,2),
     ~     i,iarr(sizex,2),iarri(sizex,2),ibuff(sizex,sizey),
     ~     iconn(nsmax,3),idir(nsmax,3),idiff,iflip,ifull,isame,isht,
     ~     istart,ixflip,iyflip,i1,i2,j,jbot(sizex,3),j1,j2,
     ~     jconn(nsmax,3),jdiff,jdir(nsmax,3),jflip,jsame,
     ~     jtop(sizex,3),jx,k,mdel1(sizex,2),mdel2(sizex,2),mi1,mj1,
     ~     mj2,n,narray(sizex,sizey),ncol,nj,posn(mres,2),top(sizex,10),
     ~     ycg,yshift,ytot
      character jsigni(sizex,2),jsign(sizex,2)
      isame=0
      idiff=0
      jsame=0
      jdiff=0
      ytot=0
      n=0
      istart=4
c     
C     CHECK IF THERE IS ANY DATA IN COLUMN 2 OF IBUFF ARRAY
      do 1000 i=1,sizey
         if(ibuff(2,i).ne.0) then 
            istart=2
            go to 1001
         endif
 1000 continue
 1001 continue
      if(ncol+jx.gt.sizex) then
         write(*,*)'*************************************************'
         WRITE(*,*)'                                                 '
         WRITE(*,*)'                WARNING !!!!                     '
         WRITE(*,*)' THE PROTEIN YOU ARE TRYING TO PLOT IS TOO BIG   '
         WRITE(*,*)' FOR THE PROGRAM. PLEASE RERUN USING PART OF THE '
         WRITE(*,*)'  PROTEIN'
         WRITE(*,*)'                                                 '
         WRITE(*,*)'*************************************************'
         ifull=1
         go to 56
      endif
      if(isht.eq.0) go to 21
      do 2 i=1,3
         do 1 j=1,nsmax
            jconn(j,i)=0
            jdir(j,i)=0
 1       continue
 2    continue
c     
      do 10 i1=2,ncol,2
         do 9 j1=1,2
            if(clmni(i1,j1).ne.0) then
               mdel1(i1,j1)=500
               do 5 i2=2,jx,2
                  do 4 j2=1,2
                     if(jsign(i2,j2).ne.jsigni(i1,j1)) then
                        del(i2,j2)=abs(clmni(i1,j1)-clmn(i2,j2))
                        if(del(i2,j2).lt.mdel1(i1,j1)) then
                           mdel1(i1,j1)=del(i2,j2)
                           mj1=j2
                           mi1=i2
                        endif
                     endif
 4                continue
 5             continue
               mdel2(i1,j1)=500
               do 25 i2=2,ncol,2
                  do 24 j2=1,2
                     if(jsigni(i1,j1).ne.jsigni(i2,j2)) then
                        del(i2,j2)=abs(clmni(i1,j1)-clmni(i2,j2))
                        if(del(i2,j2).lt.mdel2(i1,j1).and.i2.ne.i1) then
                           mdel2(i1,j1)=del(i2,j2)
                           mj2=j2
                        endif
                     endif
 24               continue
 25            continue
               if(mdel1(i1,j1).lt.mdel2(i1,j1)) then
                  jsame=jsame+i1
                  jdiff=jdiff+ncol-i1
                  if(mj1.eq.j1) then
                     isame=isame+1
                  else
                     if(mj1.eq.abs(3-j1)) idiff=idiff+1
                  endif
               endif
            endif
 9       continue
 10   continue
c     
 21   if(idiff.gt.isame.and.iflip.eq.0) then
         iyflip=1
      else if(idiff.le.isame.and.iflip.eq.1) then
         iyflip=1
      else
         iyflip=0
      endif
c     
      if(jdiff.gt.jsame.and.jflip.eq.0) then
         ixflip=1
      else if(jdiff.le.jsame.and.jflip.eq.1) then
         ixflip=1
      else
         ixflip=0
      endif
C     CALCULATE CENTRE OF GRAVITY OF THE SHEET TO ALLOW CORRECT Y POSITION
      do 23 i=istart,ncol,2
         do 22 j=1,sizey
            if(ibuff(i,j).ne.0) then
               ytot=ytot+j
               n=n+1
            endif
 22      continue
 23   continue
c     
      if(n.ne.0) then 
         ycg=ytot/n
      else
         ycg=30
      endif
c     
      if(iyflip.eq.1) then
         yshift=ycg-30
      else
         yshift=30-ycg
      endif            
c     
C     INSERT INTO MAIN ARRAY
        if(iyflip.eq.0) then
           do 15 i=istart,ncol,2
              do 14 j=1,sizey
                 if(ibuff(i,j).ne.0) then
                    n=ibuff(i,j)
                    if(ixflip.eq.0) then
                       narray(jx+i,j+yshift)=n
                       posn(n,1)=jx+i
                    else
                       narray(jx+ncol-i+4,j+yshift)=n
                       posn(n,1)=jx+ncol-i+istart
                    endif
                    posn(n,2)=j+yshift
                 endif
 14           continue
c     
              if(ixflip.eq.1) then
                 nj=jx+ncol-i+istart
              else
                 nj=jx+i
              endif
c     
              do 13 j=1,2
                 jsign(nj,j)=jsigni(i,j)
                 clmn(nj,j)=clmni(i,j)
                 iarr(nj,j)=iarri(i,j)+yshift
 13           continue
              do 16 j=1,3
                 top(nj,j)=jtop(i,j)
                 bot(nj,j)=jbot(i,j)
 16           continue
 15        continue
        else
           do 20 i=istart,ncol,2
              if(ixflip.eq.0) then
                 nj=jx+i
              else
                 nj=jx+ncol-i+istart
              endif
              do 19 j=1,sizey
                 if(ibuff(i,j).ne.0) then
                    n=ibuff(i,j)
                    narray(nj,sizey-j+1+yshift)=n
                    posn(n,1)=nj
                    posn(n,2)=sizey + 1 - j + yshift
                 endif
 19           continue
              do 18 j=1,2
                 jsign(nj,3-j)=jsigni(i,j)
                 clmn(nj,3-j)=clmni(i,j)
                 if(iarri(i,j).ne.0) iarr(nj,3-j)=sizey + 1 - iarri(i,j)
     -              + yshift
 18           continue
              do 211 j=1,3
                 top(nj,j)=jbot(i,j)
                 bot(nj,j)=jtop(i,j)
 211          continue
 20        continue
        endif
        jx=jx+ncol
        if(ixflip.eq.1) then
           do 35 i=1,nsmax
              if(iconn(i,1).eq.0) then
                 k=i
                 go to 40
              endif
 35        continue
 40        continue
           do 45 i=1,k-1
              do 44 j=1,3
                 jconn(i,j)=iconn(k-i,j)
                 jdir(i,j)=idir(k-i,j)
 44           continue
 45        continue
           do 50 i=1,k-1
              do 49 j=1,3
                 iconn(i,j)=jconn(i,j)
                 idir(i,j)=jdir(i,j)
 49           continue
 50        continue
        endif
c     
        if(iyflip.eq.1) then
           do 55 i=1,20
              do 54 j=1,3
                 idir(i,j)=-idir(i,j)
 54           continue
 55        continue
        endif
 56     return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
        subroutine adist(ibp,sstr,jj,sht,sheet,ks,nr,str,nrow,
     ~       minres,maxres,posn,istps,nsrow,
     ~       shtref,numstr,numsht,lsheet)
C     EXTRACT RESIDUES IN SHEET AND THEIR PARTNERS AND ASSIGN STRAND
C     NUMBERS
        include 'p_hera.par'

        integer del(4),delt,i,ibp(mres,2),
     ~       ij,il,ipos(mres),ires,istps(mstr,7),i1,i2,j,jj,jres,
     ~       j1,j2,k,kres,ks,l,maxres,minres,nr,nrow,nsrow,
     ~       numsht,posn(mres,2),sheet(300,6),str(nsmax,2),
     ~       numstr(msht),sht(mres)
        character sstr(mres),hlxchr*4,shchr1*2,shtref*62,
     ~       lsheet(msht)
        data hlxchr/'hHgG'/,shchr1/'eE'/
        nrow=1
        nr=0
        ks=1
        numstr(jj)=0
        do 147 i=1,nsmax
           do 146 j=1,2
              str(i,j)=0
 146       continue
 147    continue
        do 5 i=1,mres
           ipos(i)=0
 5      continue
        do 7 i=1,300
           do 6 j=1,6
              sheet(i,j)=0
 6         continue
 7      continue
c     FIRST RECORD THE RESIDUES IN THE SHEET AND THEIR BRIDGE PARTNERS
        do 10 i=minres,maxres
           if(index(hlxchr,sstr(i)).eq.0) then
C     EXCLUDE ALL HELICAL RESIDUES
              if(i.eq.1) then
                 if(index(shchr1,sstr(i)).ne.0) then
                    if(sht(i).eq.jj.or.sht(i+1).eq.jj) then
                       nr=nr+1
                       sheet(nr,1)=i
                       ipos(i)=nr
                       do 9 j=1,2
                          if(ibp(i,j).ge.minres.and.ibp(i,j).le.maxres) 
     &                         then
                             k=ibp(i,j)
                             sheet(nr,j+1)=k
                          endif
 9                     continue
                    endif
                 endif
              else
                 if(sstr(i).eq.'E') then
                    if(sht(i).eq.jj.or.(sht(i-1).eq.jj.and.
     -                   sht(i+1).eq.jj)) then
                       nr=nr+1
                       sheet(nr,1)=i
                       ipos(i)=nr
                       do 12 j=1,2
                          if(ibp(i,j).ge.minres.and.ibp(i,j).le.maxres)
     ~                         then
                             k=ibp(i,j)
                             sheet(nr,j+1)=k
                          endif
 12                    continue
                    endif
                 else if(sstr(i).eq.'e'.and.(sht(i-1).eq.jj.or.
     -                   sht(i+1).eq.jj).and.posn(i,1).eq.0) then
                    nr=nr+1
                    sheet(nr,1)=i
                    ipos(i)=nr
                    do 13 j=1,2
                       if(ibp(i,j).ge.minres.and.ibp(i,j).le.maxres) 
     ~                      then
                          k=ibp(i,j)
                          sheet(nr,j+1)=k
                       endif
 13                 continue
                 else
       endif
      endif
      endif
 10   continue
C     NOW ALLOCATE STRAND NUMBERS WITHIN SHEET
      sheet(1,4)=1
      str(1,1)=sheet(1,1)
      do 20 i=2,nr
         ires=sheet(i,1)
         jres=sheet(i-1,1)
         kres=sheet(i+1,1)
         if(ires.eq.jres+1.and.
     ~        (sstr(jres).eq.'E'.or.str(nrow,1).eq.jres)) then
            sheet(i,4)=ks
         else if(abs(ires-jres).lt.5.and.sstr(ires).ne.'e'.and.
     ~           sstr(jres).ne.'e') then
            i1=sheet(i,2)
            i2=sheet(i,3)
            j1=sheet(i-1,2)
            j2=sheet(i-1,3)
            del(1)=abs(i1-j1)
            del(2)=abs(i1-j2)
            del(3)=abs(i2-j1)
            del(4)=abs(i2-j2)
            delt=100
            do 11 ij=1,4
               if(del(ij).lt.delt.and.del(ij).ne.0) delt=del(ij)
 11         continue
            if(delt.le.5.and.delt.ne.0) then
               if(i1.ne.jres.and.i2.ne.jres) then
                  sheet(i,4)=ks
               else
                  ks=ks+1
                  sheet(i,4)=ks
               endif
            else
               ks=ks+1
               sheet(i,4)=ks
            endif            
         else
            ks=ks+1
            sheet(i,4)=ks
         endif
         if(sheet(i,4).ne.sheet(i-1,4)) then
            nrow=nrow+1
            str(nrow,1)=sheet(i,1)
         endif
         if(sheet(i,4).ne.sheet(i+1,4).or.i.eq.nr) 
     ~        str(nrow,2)=sheet(i,1)
 20   continue
C     NOW ALLOCATE STRAND NUMBERS TO PARTNERS
      do 30 i=1,nr
         do 25 j=2,3
            if(sheet(i,j).ne.0) then
               l=sheet(i,j)
               il=ipos(l)
               if(il.ne.0) then
                  if(j.eq.2) then
                     sheet(i,5)=sheet(il,4)
                  else
                     sheet(i,6)=sheet(il,4)
                  endif
               endif
            endif
 25      continue
 30   continue
      if(nrow.gt.1) then
         do 130 i=1,nrow
            nsrow=nsrow+1
            istps(nsrow,1)=str(i,1)
            istps(nsrow,2)=str(i,2)
            istps(nsrow,3)=i
            istps(nsrow,4)=jj
 130     continue
         numsht=numsht+1
         numstr(numsht)=ks
         lsheet (numsht)=shtref(jj:jj)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine combin(sheet,str,iconn,ks,nr,jjoin,barr)
      include 'p_hera.par'
C     CALCULATE STRAND CONNECTIVITY IN SHEET
      integer barr,carr(nsmax,3),freq(nsmax),i,ibnd,ibuff(2),
     ~     icol(nsmax),iconn(nsmax,3),id,id1,id2,ie,iend,il,
     ~     iord,islen,istr,iuse(nsmax),j,jcol(nsmax),
     ~     jconn(nsmax,3),jjoin(3*nsmax,3),jmc,jp(2),jr(2),
     ~     juse(nsmax),k,ks,l,lplen,m,mc,minlen,mk,n,nend,nr,
     ~     sheet(300,6),str(nsmax,2),strp

      mc=0
      jmc=0
      ie=0
      barr=0
      ibnd=0
      lplen=0
      do 4 i=1,nsmax
         freq(i)=0
         iuse(i)=0
         juse(i)=0
         icol(i)=0
         jcol(i)=0
         do 44 j=1,3
            iconn(i,j)=0
            jconn(i,j)=0
 44      continue
         do 45 j=1,3
            carr(i,j)=0
 45      continue
 4    continue
      do 8 i=1,2
         ibuff(i)=0
 8    continue
      do 391 j=1,3*nsmax
         do 390 i=1,3
            jjoin(j,i)=0
 390     continue
 391  continue
      
      do 40 i=1,nr
         do 35 j=5,6
            jmc=0
            ibuff(1)=sheet(i,4)
            ibuff(2)=sheet(i,j)
            if(ibuff(2).ne.0.and.ibuff(2).ne.ibuff(1)) then
               call compar(ibuff,mc,carr,jmc)
               if(jmc.eq.0) then
                  mc=mc+1
                  carr(mc,1)=min0(ibuff(1),ibuff(2))
                  carr(mc,2)=max0(ibuff(1),ibuff(2))
               else
                  carr(jmc,3)=carr(jmc,3)+1           
               endif
            endif
 35      continue
 40   continue
C     SHOULD PRODUCE A UNIQUE SET OF CONNECTED PAIRS-NOW FIND THE ENDS
      do 50 i=1,ks
         do 48 j=1,mc
            do 46 l=1,2
               if(carr(j,l).eq.i) freq(i)=freq(i)+1
 46         continue
 48      continue
 50   continue
C     NOW FIND THE STRAND ORDER
      islen=0
      iend=0
      nend=0
      do 70 i=1,ks
         if(freq(i).eq.1) then
            il=str(i,2)-str(i,1)
            nend=nend+1
            if(il.gt.islen.or.islen.eq.0) then
               islen=il
               iend=i
            endif
         endif
 70   continue
c     DEFINE BARRELS LT 2 END STRANDS        
      if(nend.lt.2) barr=1
      if(iend.ne.0) then
         iconn(1,1)=iend
c     the following line was moved up 3 lines 18.1.95 - TEST
         iord=1
         icol(iend)=iord
         istr=1
         iuse(iend)=1
         call conect(iord,carr,iconn,mc,iuse,istr,icol,ks,jjoin,freq,
     ~        str)
      else
c     BARREL
         do 78 i=1,ks
c            if(freq(i).eq.2) then
               minlen=100
               n=0
               do 72 j=1,mc
                  do 71 k=1,2
                     if(carr(j,k).eq.i) then
                        n=n+1
                        jp(n)=strp(carr,j,k)
                        jr(n)=j
                     endif
 71               continue
 72            continue
               do 74 j=1,2
                  do 200 m=1,nsmax
                     iuse(m)=0
                     juse(m)=0
                     jcol(m)=0
                     icol(m)=0
 200              continue
                  jconn(1,1)=i
                  jcol(i)=1
                  juse(i)=1
                  jconn(2,1)=jp(j)
                  jcol(jp(j))=2
                  juse(jp(j))=1
                  if(j.eq.1) then
                     ibnd=jr(2)
                  else
                     ibnd=jr(1)
                  endif   
                  iord=2
                  istr=2
                  call conect(iord,carr,jconn,mc,juse,istr
     1                 ,jcol,ks,jjoin,freq,str)
                  do 75 n=1,ks-1
                     id1=jcol(n)
                     id2=jcol(n+1)
                     id=abs(id2-id1)
                     lplen=lplen+id
 75               continue
                  lplen=lplen+ibnd
                  if(lplen.le.minlen) then
                     minlen=lplen
                     do 77 k=1,nsmax
                        do 67 mk=1,3
                           iconn(k,mk)=jconn(k,mk)
 67                     continue
 77                  continue
                     iconn(istr+1,1)=iconn(1,1)
                  endif                
 74            continue
 78         continue
         endif
         return
         end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        subroutine compar(ibuff,mc,carr,jmc)
        include 'p_hera.par'
        integer carr(nsmax,3),i,ibuff(2),jmc,mc
        jmc=0
        do 10 i=1,mc
           if(carr(i,1).eq.ibuff(1).and.carr(i,2)
     1          .eq.ibuff(2)) then
              jmc=i
              go to 11
           else
              if(carr(i,1).eq.ibuff(2).and.carr(i,2)
     1             .eq.ibuff(1)) then
                 jmc=i
                 go to 11
              endif
           endif
 10     continue
 11     return
        end
C     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     
        integer function strp(carr,m,n)
C     IDENTIFY STRAND PARTNER 
        include 'p_hera.par'
        integer carr(nsmax,3),m,n
        if(n.eq.2) then
           strp=carr(m,1)
        else
           strp=carr(m,2)
        endif
        return
        end
C     
CcCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     
c     CALCULATION ROUTINES - LOOP CALCULATIONS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine lsrt13(arr,n)
        include 'p_hera.par'
        integer arr(lmax,7),i,ix1,ix2,j,ja,jb,jj,jx1,jx2,jy1,jy2,k,n
        real deli,delj
c     
c     SORTS THE LOOPS INTO INCREASING LENGTH -LOOP3, LOOP1
        do 12 j=2,n
           jx1=arr(j,1)
           jx2=arr(j,3)
           jy1=arr(j,2)
           jy2=arr(j,4)
           jj=arr(j,5)
           ja=arr(j,6)
           jb=arr(j,7)
           delj=abs(jx1-jx2)
           do 11 i=j-1,1,-1
              ix1=arr(i,1)
              ix2=arr(i,3)
              deli=abs(ix1-ix2)
              if(deli.le.delj) go to 10
              do 9 k=1,7
                 arr(i+1,k)=arr(i,k)
 9            continue
 11        continue
           i=0
 10        arr(i+1,1)=jx1
           arr(i+1,2)=jy1
           arr(i+1,3)=jx2
           arr(i+1,4)=jy2
           arr(i+1,5)=jj
           arr(i+1,6)=ja
           arr(i+1,7)=jb
 12     continue
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine lsrt2a(arr,n)
        include 'p_hera.par'
        integer arr(lmax,7),deli,delj,i,ix1,ix2,j,ja,jb,jj,
     ~       jx1,jx2,jy1,jy2,k,n
c     
c     SORTS THE LOOPS INTO INCREASING LENGTH - LOOP2A
        do 12 j=2,n
           jx1=arr(j,1)
           jx2=arr(j,3)
           jy1=arr(j,2)
           jy2=arr(j,4)
           jj=arr(j,5)
           ja=arr(j,6)
           jb=arr(j,7)
           delj=jx1+jx2
           do 11 i=j-1,1,-1
              ix1=arr(i,1)
              ix2=arr(i,3)
              deli=ix1+ix2
              if(deli.le.delj) go to 10
              do 9 k=1,7
                 arr(i+1,k)=arr(i,k)
 9            continue
 11        continue
           i=0
 10        arr(i+1,1)=jx1
           arr(i+1,2)=jy1
           arr(i+1,3)=jx2
           arr(i+1,4)=jy2
           arr(i+1,5)=jj
           arr(i+1,6)=ja
           arr(i+1,7)=jb
 12     continue
        return
        end
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine lsrt2b(arr,n)
        include 'p_hera.par'
        integer arr(lmax,7),deli,delj,i,ix1,ix2,j,ja,jb,jj,jx1,
     ~       jx2,jy1,jy2,k,n
c     
c     SORTS THE LOOPS INTO decreasing LENGTH - LOOP2B
        do 12 j=2,n
           jx1=arr(j,1)
           jx2=arr(j,3)
           jy1=arr(j,2)
           jy2=arr(j,4)
           jj=arr(j,5)
           ja=arr(j,6)
           jb=arr(j,7)
           delj=jx1+jx2
           do 11 i=j-1,1,-1
              ix1=arr(i,1)
              ix2=arr(i,3)
              deli=ix1+ix2
              if(deli.ge.delj) go to 10
              do 9 k=1,7
                 arr(i+1,k)=arr(i,k)
 9            continue
 11        continue
           i=0
 10        arr(i+1,1)=jx1
           arr(i+1,2)=jy1
           arr(i+1,3)=jx2
           arr(i+1,4)=jy2
           arr(i+1,5)=jj
           arr(i+1,6)=ja
           arr(i+1,7)=jb
 12     continue
        return
        end
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine dpsnlp(posn,lp,nlrow,istpl,minres,maxres)
        include 'p_hera.par'

C     DEFINE THE BEGINNING AND END OF LOOP REGIONS -TO BE DRAWN AS
C     LINES
        integer i,istpl(lmax,2),lp(lmax,4),maxres,minres,n,nlrow,
     ~       posn(mres,2)

        nlrow=1
        lp(nlrow,1)=0
        lp(nlrow,2)=0
        if(posn(minres,1).eq.0) then
           istpl(nlrow,1)=minres
        else
           lp(nlrow,3)=posn(minres,1)
           lp(nlrow,4)=posn(minres,2)
           istpl(nlrow,1)=1
           istpl(nlrow,2)=1
        endif
        if(posn(minres+1,1).ne.0.and.posn(minres,1).eq.0) 
     ~       call endlp(lp,nlrow,minres+1,posn,istpl)
        do 10 n=minres+1,maxres-1
           if(posn(n-1,1).ne.0) then
              if(posn(n,1).eq.0) then
                 call addlp(lp,nlrow,n-1,posn,istpl)
              else if(posn(n,1).ne.posn(n-1,1)) then
                 call addlp(lp,nlrow,n-1,posn,istpl)
                 call endlp(lp,nlrow,n,posn,istpl)
                 go to 10
              else if(abs(posn(n,2)-posn(n-1,2)).ne.1) then
                 call addlp(lp,nlrow,n-1,posn,istpl)
                 call endlp(lp,nlrow,n,posn,istpl)
                 go to 10
              else
                 if(n.eq.maxres-1.and.posn(maxres,1).eq.0) 
     ~                call addlp(lp,nlrow,n,posn,istpl)
              endif
           endif
           if(posn(n+1,1).ne.0.and.posn(n,1).eq.0) 
     ~          call endlp(lp,nlrow,n+1,posn,istpl)
 10     continue
        if(posn(maxres,1).eq.0) then
           lp(nlrow,3)=1000
           lp(nlrow,4)=1000
           istpl(nlrow,2)=maxres
        else
           nlrow=nlrow+1
           lp(nlrow,1)=posn(maxres,1)
           lp(nlrow,2)=posn(maxres,2)
           do 11 i=3,4
              lp(nlrow,i)=1000
 11        continue
           istpl(nlrow,1)=maxres
           istpl(nlrow,2)=maxres
        endif        
        return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine addlp(lp,nlrow,m,posn,istpl)
        include 'p_hera.par'
C     START A NEW LOOP
        integer istpl(lmax,2),lp(lmax,4),m,nlrow,posn(mres,2)

        nlrow=nlrow+1
        lp(nlrow,1)=posn(m,1)
        lp(nlrow,2)=posn(m,2)
        istpl(nlrow,1)=m

        return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
        subroutine endlp(lp,nlrow,m,posn,istpl)
        include 'p_hera.par'
C     END OF EXISTING LOOP
        integer istpl(lmax,2),lp(lmax,4),m,nlrow,posn(mres,2)

        lp(nlrow,3)=posn(m,1)
        lp(nlrow,4)=posn(m,2)
        istpl(nlrow,2)=m

        return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        CALCULATION ROUTINES -HELIX CALCULATIONS
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine hlxrec(sstr,istph,n,nhrow,minres,maxres,maxlen)
        include 'p_hera.par'
C        TRANSFERS THE INITIAL AND FINAL RESIDUES OF HELIX TO ARRAY ISTPH
C        BEARING IN MIND THAT THE FIRST AND LAST RESIDUES ARE NOT ACTUALLY
C        CLASSIFIED AS BEING IN HELIX, THOUGH THEY ARE HYDROGEN BONDED

        integer istph(mhel,2),nhrow,n,maxlen,
     ~       maxres,minres,len
        character sstr(mres),hlxchr*2,hxchr2*2
        data hlxchr/'HG'/,hxchr2/'hg'/

        if(n.eq.minres.or.n.eq.1) then
           istph(nhrow,1)=n
        else 
           if(index(hlxchr,sstr(n-1)).eq.0) then
c              if(dssp.eq.0) then
                 if(index(hxchr2,sstr(n-1)).eq.0) then
                    istph(nhrow,1)=n
                 else if(n.le.2) then
                    istph(nhrow,1)=n-1
                 else
                    if(index(hlxchr,sstr(n-2)).ne.0) then
                       istph(nhrow,1)=n
                    else
                       istph(nhrow,1)=n-1
                    endif
                 endif
c              else
c                 istph(nhrow,1)=n-1
c              endif
           endif
        endif
        if(n.eq.maxres) then
           istph(nhrow,2)=n
           len=istph(nhrow,2)-istph(nhrow,1)
           if(len.gt.maxlen) maxlen=len
           nhrow=nhrow+1
        else
           if(index(hlxchr,sstr(n+1)).eq.0) then
c              if(dssp.eq.0) then
                 if(index(hxchr2,sstr(n+1)).eq.0) then
                    istph(nhrow,2)=n
                 else 
                    istph(nhrow,2)=n+1
                 endif
c              else
c                 istph(nhrow,2)=n+1
c              endif
              len=istph(nhrow,2)-istph(nhrow,1)
              if(len.gt.maxlen) maxlen=len
              nhrow=nhrow+1
           endif
        endif
        return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     
        subroutine hins(ifun,m,n,jx,jys,narray,posn,sgn)
        include 'p_hera.par'

C     INSERT HELIX INTO NARRAY (FINAL ARRAY)
        integer ifun(mres),i1,i2,j,jx,jy,jys,m,n,posn(mres,2),
     ~       narray(sizex,sizey),sgn
        
        if(jx.gt.sizex.or.jx.le.0) then
          write(*,*)' no room left in output array, jx = ', jx
          go to 100
        endif
        jy=jys
        if(jy.gt.sizey.or.jy.le.0) then
          write(*,*)' no room left in output array, jy = ', jy
          go to 100
        endif
        if(sgn.eq.1) then
           i1=min0(m,n)
           i2=max0(m,n)
        else
           i1=max0(m,n)
           i2=min0(m,n)
        endif
        do 10 j=i1,i2,sgn
           if(ifun(j).ne.0.and.posn(j,1).eq.0) then
              narray(jx,jy)=j
              posn(j,1)=jx
              posn(j,2)=jy
              jy=jy+1
           endif
 10     continue
 100    return
        end        
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        integer function iscen(m,n)
        integer m,n
        
        iscen=int(30-abs((m-n)/2))
        return
        end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        subroutine dpnhel(clmn,istph,narray,posn,nrow,numres,
     ~       iseqno,iarr,top,bot,jx,ntop)
        include 'p_hera.par'

C     THIS SUBROUTINE FINDS THE BEST COLUMN TO INSERT THE HELIX
C     AND INSERTS THE HELIX IN THE CORRECT ORIENTATION AND EITHER
C     CENTRED OR SHIFTED UP OR DOWN ACCORDING TO WHAT ELSE IS IN
C     THE COLUMN.
C     FOR DRAWING HELICES ONLY
        integer bot(sizex,10),clmn(sizex,2),del(sizex,2),
     ~       delt(sizex),hlxpsn(2),i,iarr(sizex,2),ic1,ic2,
     ~       ie,ii,iscen,iseqno(mres),istph(mhel,2),j,jx,jys,
     ~       mindel,mincol,narray(sizex,sizey),nrow,numres,posn(mres,2),
     ~       sgn,top(sizex,10),ntop

c        sgn=1
        mindel=numres
        do 10 j=1,2
          do 9 i=2,sizex,3
            del(i,j)=abs(clmn(i,j)-istph(nrow,1))
C        FINDS NEAREST RESIDUE TO X,X IS AT BEGINNING OF HELIX
            if(del(i,j).lt.mindel.and.clmn(i,j).ne.0) then
              mindel=del(i,j)
              hlxpsn(1)=i
              hlxpsn(2)=j
            endif
9          continue
10        continue
        mincol=30
        if(hlxpsn(1).eq.0) then
           hlxpsn(1)=1
           if(ntop.eq.1) then
              hlxpsn(2)=2
           else
              hlxpsn(2)=1
           endif
        endif
        do 20 i=2,sizex,3
          delt(i)=abs(i-hlxpsn(1))
          if(delt(i).lt.mincol.and.clmn(i,1).eq.0) then
            mincol=delt(i)
            jx=i
          endif
20        continue
        ii=istph(nrow,1)
        ie=istph(nrow,2)
        if(hlxpsn(2).eq.2) then
          sgn=-1
        else
          sgn=1
        endif
        jys=iscen(ii,ie)
        if(sgn.eq.1) then
           if(posn(ii,1).ne.0) ii=ii+1
           if(posn(ie,1).ne.0) ie=ie-1
           clmn(jx,1)=ii
           bot(jx,1)=ii
           clmn(jx,2)=ie
           top(jx,1)=ie
        else
           if(posn(ii,1).ne.0) ii=ii+1
           if(posn(ie,1).ne.0) ie=ie-1
           clmn(jx,2)=ii 
           top(jx,1)=ii
           clmn(jx,1)=ie
           bot(jx,1)=ie
        endif
        call hins(iseqno,ii,ie,jx,jys,narray,posn,sgn)
        ic1=clmn(jx,1)
        ic2=clmn(jx,2)
        iarr(jx,1)=posn(ic1,2)
        iarr(jx,2)=posn(ic2,2)
        return
        end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
        subroutine psnhel(clmn,istph,narray,posn,nhrow,
     ~       iseqno,iarr,ins,jsign,top,bot,jx)

C     FINDS BEST PLACE TO PUT EACH HELIX -FOR DIAGRAMS IN WHICH HELICES
C     AND SHEETS ARE TO BE DRAWN 
        include 'p_hera.par'
        integer bot(sizex,10),clmn(sizex,2),del,hlxpsn(mhel,3),
     ~       i,iarr(sizex,2),ichang,ins(mhel),irej,irow(sizex,2),
     ~       iseqno(mres),istph(mhel,2),iwait(mhel),j,jx,k,l,
     ~       mindel(mhel),n,narray(sizex,sizey),nhrow,posn(mres,2),
     ~       top(sizex,10)
        character jsign(sizex,2)

        do 100 i=1,2
           do 90 j=1,sizex
              irow(j,i)=0
 90        continue
 100    continue
 1      do 10 l=1,nhrow
           iwait(l)=0
           if(ins(l).eq.0) then
              do 5 i=1,3
                 hlxpsn(l,i)=0
 5            continue
              mindel(l)=40
              do 8 k=1,2
                 do 7 j=1,2
                    do 6 i=2,sizex,2
                       if(clmn(i,j).ne.0.and.irow(i,j).eq.0) then
                          if(jsign(i,j).eq.'b'.and.k.eq.1) then
                             del=istph(l,k)-clmn(i,j)
                          else if(jsign(i,j).eq.'a'.and.k.eq.2) 
     ~                            then
                             del=clmn(i,j)-istph(l,k)
                          else
                             del=1000
                          endif
                          if(del.lt.mindel(l).and.del.ge.0) then
                             mindel(l)=del
                             hlxpsn(l,1)=i
                             hlxpsn(l,2)=j
                             hlxpsn(l,3)=k
                          endif
                       endif
 6                  continue
 7               continue
 8            continue
              do 9 k=1,l-1
                 if(hlxpsn(k,1).eq.hlxpsn(l,1).and.hlxpsn(k,2).
     ~                eq.hlxpsn(l,2).and.ins(k).eq.0) then
                    if(mindel(k).lt.mindel(l)) then
                       iwait(l)=1
                    else
                       if(hlxpsn(l,1).ne.0) iwait(k)=1
                    endif
                 endif
 9              continue
           endif
 10     continue
        ichang=0
        do 30 l=1,nhrow
           if(ins(l).eq.0) then
              n=l
              if(hlxpsn(n,1).ne.0.and.iwait(n).eq.0) then
                 ins(n)=1
                 ichang=1
                 irej=0
                 call hcalc(clmn,narray,posn,n,hlxpsn,iseqno,
     ~                iarr,top,bot,jsign,irej,ins,jx,mindel,irow,istph,
     ~                nhrow)
              endif
              if(ichang.eq.1) then
                 go to 1
              endif
           endif
 30        continue
           do 40 l=1,nhrow
              if(ins(l).ne.1) then
                 irej=0
                 call hcalc(clmn,narray,posn,l,hlxpsn,iseqno,
     ~                iarr,top,bot,jsign,irej,ins,jx,mindel,irow,istph,
     ~                nhrow)
                 ins(l)=1
                 go to 1
              endif
 40        continue
           return
           end
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
        subroutine search(clmn,jx,jys,sgn,hlxpsn,ii,ie
     1       ,k,ins)
C     FOR INSERTING HELICES
        include 'p_hera.par'
        integer clmn(sizex,2),hlxpsn(mhel,3),idiff,ie,ii,
     ~       ins(mhel),j1,jjx,jx,jys,k,mindif,sgn
        jjx=0
        mindif=90
        do 10 j1=4,sizex,2
           if(clmn(j1,1).eq.0) then
              idiff=abs(j1-jx)
              if(idiff.lt.mindif) then
                 jjx=j1
                 mindif=idiff
              endif
           endif
 10     continue
        if(jjx.ne.0) then
           jx=jjx
           if(hlxpsn(k,2).eq.2) then
              if(hlxpsn(k,3).eq.1) then
                 sgn=-1
              else
                 sgn=1
              endif
           else
              if(hlxpsn(k,3).eq.1) then
                 sgn=1
              else
                 sgn=-1
              endif
           endif
           if(k.gt.1) then
              if(ins(k-1).eq.1) then
                 if(sgn.eq.1) then
                    jys=3
                 else
                    jys=sizey - 2 - (ie-ii+1)
                 endif
              else
                 if(sgn.eq.-1) then
                    jys=3
                 else
                    jys=sizey - 2- (ie-ii+1)
                 endif
              endif
           else
              if(sgn.eq.1) then
                 jys=3
              else
                 jys=sizey - 2 - (ie-ii+1)
              endif
           endif
        else
           write(*,*) ' no room for helix'
        endif
        return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        integer function sgn(x,y)
        integer x
        character*1 y

        if(y.eq.'+') sgn=x
        if(y.eq.'-') sgn=-x
        if(y.eq.' ') sgn=0
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ttop1(array,jx,i)
        include 'p_hera.par'
        integer array(sizex,10),i,jx,k

        do 5 k=1,10
         if(array(jx,k).eq.0) then
          array(jx,k)=i
          go to 6
         endif
5        continue
6        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hcalc(clmn,narray,posn,n,hlxpsn,iseqno,
     ~       iarr,top,bot,jsign,irej,ins,maxx,mindel,irow,istph,
     ~       nhrow)
C     FOR INSERTING HELICES
      
      include 'p_hera.par'
      integer narray(sizex,sizey),clmn(sizex,2),iarr(sizex,2),
     ~     top(sizex,10),bot(sizex,10),ib,ic1,ic2,ie,ien,ifind,
     ~     ii,ins(mhel),irb,irej,irt,it,j,ji,jys,k,maxx,
     ~     hlxpsn(mhel,3),n,sgn,iseqno(mres),nhrow,posn(mres,2),
     ~     mindel(mhel),irow(sizex,2),istph(mhel,2),itry,jx
        character jsign(sizex,2)

        itry=0
        ii=istph(n,1)
        ie=istph(n,2)
        if(posn(ie,1).ne.0) ie=ie-1
        if(posn(ii,1).ne.0) ii=ii+1
 1      jx=hlxpsn(n,1)
        if(hlxpsn(n,2).eq.2) then
C     INSERT HELIX AT TOP OF EXISTING HELIX OR STRAND
           ib=clmn(jx,2)
           irb=posn(ib,2)
           jys=irb+2
c     CALCULATE WHERE END OF HELIX WILL BE
 121       ien=jys+ie-ii+1
C     CHECK WHETHER AREA AROUND START OF HELIX IS SAFE
           if(ien.le.sizey) then
              do 125 ji=1,4
                 if(narray(jx-1,jys+ji-3).ne.0.or.narray(jx+1,jys+ji-3)
     1                .ne.0) then
                    jys=jys+1
                    go to 121
                 endif
 125          continue 
              if(hlxpsn(n,3).eq.1) then
                 sgn=1
                 clmn(jx,2)=ie
                 jsign(jx,2)='b'
                 call ttop1(top,jx,ie)
                 call ttop1(bot,jx,ii)
              else
                 sgn=-1
                 clmn(jx,2)=ii
                 jsign(jx,2)='a'
                 call ttop1(top,jx,ii)
                 call ttop1(bot,jx,ie)
              endif 
              call hins(iseqno,ii,ie,jx,jys,narray,posn,sgn)
              ic2=clmn(jx,2)
              iarr(jx,2)=posn(ic2,2)
              if(jx.gt.maxx) maxx=jx
           else
              irej=1
              k=hlxpsn(n,1)
              j=hlxpsn(n,2)
              irow(k,j)=1
              call srcha(n,clmn,hlxpsn,mindel,irow,istph,jsign,
     ~             ifind,nhrow)
              itry=itry+1
              if(ifind.eq.1.and.itry.lt.2) go to 1
              call search(clmn,jx,jys,sgn,hlxpsn,ii,ie,n,ins) 
              call hins(iseqno,ii,ie,jx,jys,narray,posn,sgn)
              if(jx.gt.maxx) maxx=jx
              if(sgn.eq.1) then
                 clmn(jx,1)=ii
                 jsign(jx,1)='a'
                 iarr(jx,1)=posn(ii,2)
                 clmn(jx,2)=ie
                 jsign(jx,2)='b'
                 iarr(jx,2)=posn(ie,2)
                 call ttop1(bot,jx,ii)
                 call ttop1(top,jx,ie)
              else
                 clmn(jx,1)=ie
                 jsign(jx,1)='b'
                 iarr(jx,1)=posn(ie,2)
                 clmn(jx,2)=ii
                 jsign(jx,2)='a'
                 iarr(jx,2)=posn(ii,2)
                 call ttop1(bot,jx,ie)
                 call ttop1(top,jx,ii)
              endif
           endif
        else if(hlxpsn(n,2).eq.1) then
C     INSERT AT BOTTOM OF EXISTING STRAND OR HELIX
           it=clmn(jx,1)
           irt=posn(it,2)
 131       jys=irt-1-(ie-ii+1)
           if(jys.ge.5) then
              do 135 ji=1,4
                 if(irt-3-ji.gt.0) then
                    if(narray(jx-1,irt-3-ji).ne.0.or.
     ~                   narray(jx+1,irt-3-ji).ne.0) then
                       irt=irt-1
                       go to 131
                    endif
                 endif
 135          continue 
           endif
           if(jys.ge.5) then
              if(hlxpsn(n,3).eq.1) then
                 sgn=-1
                 clmn(jx,1)=ie
                 jsign(jx,1)='b'
                 call ttop1(bot,jx,ie)
                 call ttop1(top,jx,ii)
              else
                 sgn=1
                 clmn(jx,1)=ii
                 jsign(jx,1)='a'
                 call ttop1(bot,jx,ii)
                 call ttop1(top,jx,ie)
              endif
              call hins(iseqno,ii,ie,jx,jys,narray,posn,sgn)
              if(jx.gt.maxx) maxx=jx
              ic1=clmn(jx,1)
              iarr(jx,1)=posn(ic1,2)
           else
              irej=1
              k=hlxpsn(n,1)
              j=hlxpsn(n,2)
              irow(k,j)=1
              itry=itry+1
              call srcha(n,clmn,hlxpsn,mindel,irow,istph,jsign,
     ~             ifind,nhrow)
              if(ifind.eq.1.and.itry.lt.2) go to 1
              call search(clmn,jx,jys,sgn,hlxpsn,ii,ie,n,ins)
              call hins(iseqno,ii,ie,jx,jys,narray,posn,sgn)
              if(jx.gt.maxx) maxx=jx
              if(sgn.eq.1) then
                 clmn(jx,1)=ii
                 jsign(jx,1)='a'
                 iarr(jx,1)=posn(ii,2)
                 clmn(jx,2)=ie
                 jsign(jx,2)='b'
                 iarr(jx,2)=posn(ie,2)
                 call ttop1(bot,jx,ii)
                 call ttop1(top,jx,ie)
              else
                 clmn(jx,1)=ie
                 jsign(jx,1)='b'
                 iarr(jx,1)=posn(ie,2)
                 clmn(jx,2)=ii
                 jsign(jx,2)='a'
                 iarr(jx,2)=posn(ii,2)
                 call ttop1(bot,jx,ie)
                 call ttop1(top,jx,ii)
              endif
           endif
        else
c     CASE IF HLXPSN IS ZERO
           call search(clmn,jx,jys,sgn,hlxpsn,ii,ie,n,ins)
              call hins(iseqno,ii,ie,jx,jys,narray,posn,sgn)
              if(jx.gt.maxx) maxx=jx
              if(sgn.eq.1) then
                 clmn(jx,1)=ii
                 jsign(jx,1)='a'
                 iarr(jx,1)=posn(ii,2)
                 clmn(jx,2)=ie
                 jsign(jx,2)='b'
                 iarr(jx,2)=posn(ie,2)
                 call ttop1(bot,jx,ii)
                 call ttop1(top,jx,ie)
              else
                 clmn(jx,1)=ie
                 jsign(jx,1)='b'
                 iarr(jx,1)=posn(ie,2)
                 clmn(jx,2)=ii
                 jsign(jx,2)='a'
                 iarr(jx,2)=posn(ii,2)
                 call ttop1(bot,jx,ie)
                 call ttop1(top,jx,ii)
              endif

        endif 
        return
        end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
c     
c     POSTSCRIPT PLOTTING ROUTINES
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine plotar(narray,aa,iplt1,jseqno,title,iseqno,
     ~       iplot,iplt2,shape,maxlen,cif,eint,efp,echar,
     ~       icolp,icolq,paginc,insert,resnam,chn,chain,nres,minres,
     -       maxres)
C     PLOT RESIDUES IN NARRAY - STRANDS AND/OR HELICES
        include 'p_hera.par'

        integer eint(mres,3),i,icolp,icolq,ij,iplot,
     ~       iplt1,iplt2,iseqno(mres),ix,iy,izero,j,jj,jseqno(mres),
     ~       maxlen,minres,maxres,n,narray(sizex,sizey),nres,
     ~       paginc(npage,3),xcoord,ycoord
        real efp(mres,3)
        character aa(mres),chain,chn(mres),cif*2,echar(mres,3),
     ~       insert(mres),plotch*6,shape(mres),title*4
        character*3 resnam(mres)
        data izero/0/

        do 10 i=scfile,scfil5
           write(i,1)
           write(i,2)
           write(i,4)
           write(i,5)
 10     continue
        if(iplot.eq.1.or.(iplot.eq.3.and.maxlen.gt.35)) then
           do 15 i=scfile,scfil5
              write(i,9)
              write(i,19) title
 15        continue
        else
           write(13,21)
           write(14,22)
           write(15,23)
           write(16,24)
           write(17,241)
           do 16 i=scfile,scfil5
              write(i,9)
              write(i,191) title
 16        continue
        endif
        do 17 i=scfile,scfil5
           write(i,1)
           write(i,311)
 17     continue

 1      format(' /Helvetica-Bold findfont %')
 2      format(' [15 0 0 24 0 0] makefont setfont %')
 4      format(' 0.5 setlinewidth %')
 5      format(' 10 10 translate %')
 21     format(' 300 428 translate 90 rotate -428 -300 translate %')
 22     format(' 300 428 translate 90 rotate -1217 -300 translate %')
 23     format(' 300 428 translate 90 rotate -2006 -300 translate %')
 24     format(' 300 428 translate 90 rotate -2776 -300 translate %')
 241    format(' 300 428 translate 90 rotate -3546 -300 translate %')
 19     format(' 20 20 moveto (',a4,')',1x,'show %')
 191    format(' 20 100 moveto (',a4,')',1x,'show %')
c 3      format(' [6 0 0 12 0 0] makefont setfont %')
 311     format(' [5 0 0 12 0 0] makefont setfont %')
 9      format(' 1 0.55 scale %')
C     
C----   Write out the layout info to the hera.dat file
        CALL WRIDAT(narray,sizex,sizey,chain,chn,jseqno,resnam,insert,
     -      shape,mres,nres,minres,maxres)

C----   Write out the residue numbers
        do 100 i=1,sizex
           do 90 j=1,sizey
              plotch='      '
              if(narray(i,j).gt.0) then
                 n=narray(i,j)
                 if(iplt1.eq.1) then
                    write(plotch,'(i4)') iseqno(n)
                 else if(iplt1.eq.2) then
                    write(plotch,'(i4)') jseqno(n)
                 else if(iplt1.eq.3) then
c                    write(plotch,'(i4)') iacc(n)
                 else if(iplt1.eq.4) then
                    if(cif.eq.'i '.or.cif.eq.'I '.or.cif.eq.
     ~                   'CI'.or.cif.eq.'ci') then
                       jj=eint(n,icolp)
                       write(plotch,'(i4)') jj
                    else
                       if(cif.eq.'f '.or.cif.eq.'F '.or.cif.eq.
     ~                      'CF'.or.cif.eq.'cf') then
                          jj=int(efp(n,icolp)+0.5)
                          write(plotch,'(i4)') jj
                       endif
                    endif
                 endif
c                 if(iplt2.eq.1) then
c                    write(plotch(5:5),'(a1)') aa(n)
c                 else if(iplt2.eq.2) then
c                    write(plotch(5:5),'(a1)') insert(n)
c                 else
c                    if(iplt2.eq.3) write(plotch(5:5),'(a1)') 
c     ~                   echar(n,icolq)
c                 endif
                 write(plotch(5:5),'(a1)') insert(n)
                 if(iplt2.eq.1) then
                    write(plotch(6:6),'(a1)') aa(n)
                 else 
                    if(iplt2.eq.2) then
                       write(plotch(6:6),'(a1)') echar(n,icolq)
                    endif
                 endif
                 ix=xcoord(i)
                 iy=ycoord(j,iplot,maxlen)
                 if(plotch.ne.'      ') then
                    do 32 ij=1,npage-1
                       if(i.ge.paginc(ij,iplot).and.i.le.
     ~                      paginc(ij+1,iplot)) then
                          if(.not.bw) call colour(3,ij+12)
                          write(ij+12,30) ix,plotch
                          write(ij+12,31) iy-3,izero,plotch,plotch
                       endif
 32                 continue
 30                 format(i5,1x,'(',A6,') stringwidth pop 2 div sub %')
 31                 format(i5,1x,'moveto ',i1,1x,'(',A6,
     ~                   ') stringwidth 2 div neg exch pop rmoveto (',
     ~                   a6,') show %')
                 endif
              endif
 90        continue
 100    continue
        
        do 120 i=1,sizex
           do 110 j=1,sizey - 1
              if(narray(i,j).ne.0.and.narray(i,j+1).ne.0.and.narray(i,j)
     ~             .ne.1000.and.narray(i,j+1).ne.1000) then
                 ix=xcoord(i)
                 iy=ycoord(j,iplot,maxlen)
                 do 109 ij=1,npage-1
                    if(i.ge.paginc(ij,iplot).and.i.le.paginc
     ~                   (ij+1,iplot)) write(ij+12,111) ix,iy+9,ix,iy+15
 109             continue
              endif
 110       continue
 120    continue
c     PLOTS THE LINES BETWEEN ADJACENT RESIDUES
c     NOW DRAW BOXES,CIRCLES ETC AROUND THE RESIDUES 
        call label(narray,shape,maxlen,iplot,paginc)
 111    format(2(i5,1x),'moveto ',2(i5,1x),'lineto stroke %')
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        integer function xcoord(i)
        integer i
        
        xcoord=i*19
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        integer function ycoord(i,iplot,maxlen)
        integer i,iplot,maxlen
            
        if(iplot.eq.1.or.(iplot.eq.3.and.maxlen.gt.35)) then      
            ycoord=(i*24)+23
         else
            ycoord=(i*26)-270
         endif
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine label(narray,shape,maxlen,iplot,paginc)
      
      include 'p_hera.par'
c     THIS SUBROUTINE PLOTS THE BOXES ETC AROUND THE RESIDUES
      integer i,ij,iplot,ix,iy,j,maxlen,n,narray(sizex,sizey),
     ~     xcoord,ycoord,paginc(npage,3)
      character shape(mres)
      
      do 6 i=1,sizex
         do 5 j=1,sizey
            if(narray(i,j).gt.0) then
               n=narray(i,j)
               ix=xcoord(i)
               iy=ycoord(j,iplot,maxlen)
               IF(SHAPE(N).EQ.'C') THEN
                  do 3 ij=1,npage-1
                     if(i.ge.paginc(ij,iplot).and.i.le.paginc
     ~                    (ij+1,iplot)) then
                        if(.not.bw) call colour(5,ij+12)
                        write(ij+12,10) ix,iy
                     endif
 3                continue
 10               format(2(i5,1x),'11 0 360 arc stroke %')
               ELSE 
                  IF(SHAPE(N).EQ.'B') THEN
                     do 13 ij=1,npage-1
                        if(i.ge.paginc(ij,iplot).and.i.le.paginc
     ~                       (ij+1,iplot)) then
                           if(.not.bw) call colour(4,ij+12)
                           write(ij+12,20) ix,iy
                        endif
 13                     continue
 20                     format(2(i5,1x),'moveto box %')
                     endif
                  ENDIF
               endif
 5          continue
 6       continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine plothb(ihb,hbpsn,hbond,narray,iplt,jseqno,
     ~           hbplot,clmn,iplot,maxlen,paginc,iseqno)
        include 'p_hera.par'
C        ACTUALLY PLOTS OUT THE ARROWS FROM HYDROGEN BOND DONORS TO ACCEPTORS

        integer clmn(sizex,2),dy1,hbond(mres*4,2),hbplot(mres*4),
     ~       hbpsn(mres*4,4),i,ihb,iplot,iplt,ir,iseqno(mres),ispace,
     ~       ix1,ix2,iy1,iy2,j,jseqno(mres),jx1,jx2,jy,jy1,jy2,
     ~       maxlen,narray(sizex,sizey),paginc(npage,3),xcoord,
     ~       ycoord
        data dy1/3/

        do 5 i=scfile,scfil5
           write(i,1)
           write(i,3)
           if(.not.bw) call colour(17,i)
 5      continue
 1      format(' /Helvetica-Bold findfont [5 0 0 9 0 0] makefont',
     ~   ' setfont %')
 3      format(' 0.1 setlinewidth %')
        do 100 i=1,ihb
           ix1=hbpsn(i,1)
           iy1=hbpsn(i,2)
           ix2=hbpsn(i,3)
           iy2=hbpsn(i,4)
           jx1=xcoord(ix1)
           jy1=ycoord(iy1,iplot,maxlen)
           jx2=xcoord(ix2)
           jy2=ycoord(iy2,iplot,maxlen)
C     COORDINATES OF THE TWO ENDS OF THE HYDROGEN BONDS
           if(ix1.gt.0) then
              if(ix2.le.0) then
                 if(iplt.eq.1) then
                    ir=iseqno(hbond(i,2))
                 else
                    ir=jseqno(hbond(i,2))
                 endif
                 j=hbond(i,1)
                 jy=jy1-dy1
                 if(ix1.eq.2) then
                    if(narray(ix1-1,iy1).ne.-1) then
                       call arr(jx1,jy,'-',ir,2,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix1-1,iy1)=-1
                    else
                       call arr(jx1,jy,'+',ir,2,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix1+1,iy1)=-1
                    endif
                 else if(clmn(ix1-2,1).eq.0.and.narray(ix1-1,iy1)
     ~                   .ne.-1) then
                    call arr(jx1,jy,'-',ir,2,iplot,paginc,ix1)
                    hbplot(i)=1
                    narray(ix1-1,iy1)=-1
                 else if(clmn(ix1+2,1).eq.0.and.narray(ix1+1,iy1).
     ~                   ne.-1) then
                    call arr(jx1,jy,'+',ir,2,iplot,paginc,ix1)
                    hbplot(i)=1
                    narray(ix1+1,iy1)=-1
                 else if(narray(ix1+1,iy1).ne.-1) then
                    call arr(jx1,jy,'+',ir,2,iplot,paginc,ix1)
                    hbplot(i)=1
                    narray(ix1+1,iy1)=-1
                 else
                    call arr(jx1,jy,'-',ir,2,iplot,paginc,ix1)
                    hbplot(i)=1
                    narray(ix1-1,iy1)=-1
                 endif
              else if(abs(ix2-ix1).gt.2.or.abs(iy1-iy2).gt.5) then
                 if(iplt.eq.1) then
                    ir=iseqno(hbond(i,2))
                 else
                    ir=jseqno(hbond(i,2))
                 endif
                 j=hbond(i,1)
                 jy=jy1-dy1
                 if(ix1.eq.2) then
                    if(narray(ix1-1,iy1).ne.-1) then
                       call arr(jx1,jy,'-',ir,2,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix1-1,iy1)=-1
                    else
                       call arr(jx1,jy,'+',ir,2,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix1+1,iy1)=-1
                    endif
                 else if(clmn(ix1-2,1).eq.0.and.narray(ix1-1,iy1).
     ~                   ne.-1) then
                    call arr(jx1,jy,'-',ir,2,iplot,paginc,ix1)
                    hbplot(i)=1
                    narray(ix1-1,iy1)=-1
                 else if(narray(ix1+1,iy1).ne.-1) then
                    call arr(jx1,jy,'+',ir,2,iplot,paginc,ix1)
                    hbplot(i)=1
                    narray(ix1+1,iy1)=-1
                 else
                    if(iplt.eq.1) then
                       ir=iseqno(hbond(i,1))
                    else
                       ir=jseqno(hbond(i,1))
                    endif
                    jy=jy2+dy1
                    if(ix2.eq.2) then
                       if(narray(ix2-1,iy2).ne.-2) then
                          call arr(jx2,jy,'-',ir,1,iplot,paginc,ix1)
                          hbplot(i)=1
                          narray(ix2-1,iy2)=-2
                       else
                          call arr(jx2,jy,'+',ir,1,iplot,paginc,ix1)
                          hbplot(i)=1
                          narray(ix2+1,iy2)=-2
                       endif
                    else if(clmn(ix2-2,1).eq.0.and.narray(ix2-1,iy2)
     ~                      .ne.-2) then
                       call arr(jx2,jy,'-',ir,1,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix2-1,iy2)=-2
                    else if(clmn(ix2+2,1).eq.0.and.narray(ix2+1,iy2)
     ~                      .ne.-2)  then
                       call arr(jx2,jy,'+',ir,1,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix2+1,iy2)=-2
                    else if(narray(ix2+1,iy2).ne.-2) then
                       call arr(jx2,jy,'+',ir,1,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix2+1,iy2)=-2
                    else
                       call arr(jx2,jy,'-',ir,1,iplot,paginc,ix1)
                       hbplot(i)=1
                       narray(ix2-1,iy2)=-2
                    endif                 
                 endif
              else if(ix2.eq.ix1.and.abs(iy1-iy2).lt.6) then
                 j=(iy1+iy2)/2
                 if(narray(ix1,j).eq.0) then
                    ispace=0
                 else
                    ispace=1
                 endif
                 call hlxhb(hbpsn,i,jx1,jy1,jy2
     1                ,ispace,paginc,iplot,ix1)
                 hbplot(i)=1
              else
                 if(abs(ix1-ix2).eq.2.and.abs(iy1-iy2).le.5) then
                    call strhb(hbpsn,i,ihb,hbond,ix1,iy1,ix2,iy2
     1                   ,iplot,maxlen)
                    hbplot(i)=1
                 endif
              endif        
           else
              if(ix2.gt.0) then
                 if(iplt.eq.1) then
                    ir=iseqno(hbond(i,1))
                 else
                    ir=jseqno(hbond(i,1))
                 endif
                 j=hbond(i,2)
                 jy=jy2+dy1
                 if(ix2.eq.2) then
                    if(narray(ix2-1,iy2).ne.-2) then
                       call arr(jx2,jy,'-',ir,1,iplot,paginc,ix2)
                       hbplot(i)=1
                       narray(ix2-1,iy2)=-2
                    else
                       call arr(jx2,jy,'+',ir,1,iplot,paginc,ix2)
                       hbplot(i)=1
                       narray(ix2+1,iy2)=-2
                    endif
                 else if(clmn(ix2-2,1).eq.0.and.narray(ix2-1,iy2).
     ~                   ne.-2) then
                    call arr(jx2,jy,'-',ir,1,iplot,paginc,ix2)
                    hbplot(i)=1
                    narray(ix2-1,iy2)=-2
                 else if(clmn(ix2+2,1).eq.0.and.narray(ix2+1,iy2).
     ~                    ne.-2) then
                    call arr(jx2,jy,'+',ir,1,iplot,paginc,ix2)
                    hbplot(i)=1
                    narray(ix2+1,iy2)=-2
                 else if(narray(ix2+1,iy2).ne.-2) then
                    call arr(jx2,jy,'+',ir,1,iplot,paginc,ix2)
                    hbplot(i)=1
                    narray(ix2+1,iy2)=-2
                 else
                    call arr(jx2,jy,'-',ir,1,iplot,paginc,ix2)
                    hbplot(i)=1
                    narray(ix2-1,iy2)=-2
                 endif
              endif
           endif
 100    continue
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine arr(ix,iy,c,ir,ift,iplot,paginc,i)
C     DRAWS ARROWS
        include 'p_hera.par'
        integer dx1,dx2,dx3,paginc(npage,3),dx4,dx5,dx6,i,ift,ij,
     ~       iplot,ir,ix,ix1,iy
        character c
        data dx1/12/,dx2/10/,dx3/13/,dx4/13/,dx5/15/,
     ~       dx6/23/
        
        if(c.eq.'+') then
           if(ift.eq.1) then
              call newarb(paginc,iplot,ix+dx1,iy,ix+dx2,iy,i)
           else
              call newarb(paginc,iplot,ix+dx4,iy,ix+dx5,iy,i)
           endif
           ix1=ix+dx3
        else
           if(ift.eq.1) then
              call newarb(paginc,iplot,ix-dx1,iy,ix-dx2,iy,i)
           else
              call newarb(paginc,iplot,ix-dx4,iy,ix-dx5,iy,i)
           endif
           ix1=ix-dx6
        endif
        do 5 ij=1,npage-1
           if(i.ge.paginc(ij,iplot).and.i.le.paginc(ij+1,iplot))
     ~          write(ij+12,10) ix1,iy-3,ir
 5      continue
 10     format(2(i5,1x),' moveto (',i4,') show %')
        
        return
        end
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine hlxhb(hbpsn,i,jx1,jy1,jy2
     1      ,ispace,paginc,iplot,ix)
C     DRAWS HYDROGEN BONDS FOR HELICES
        include 'p_hera.par'
        integer hbpsn(mres*4,4),i,idel,ii,ij,
     ~       iplot,ispace,ix,ix1,ix2,ix3,ix3a,ix3b,ix4,ix5,
     ~       iy1,iy2,iy3,iy3a,iy3b,iy4,iy5,j,jx1,jy1,
     ~       jy2,paginc(npage,3)

        j=abs(hbpsn(i,2)-hbpsn(i,4))
c        SHOULD BE 3 OR 4 ACCORDING TO WHETHER THE HELIX IS ALPHA OR 3,10
        if(j.eq.3) then
C     PI HELIX
           if(hbpsn(i,2).gt.hbpsn(i,4)) then
              ix1=jx1-18
              iy1=jy1-15
              ix2=jx1-22
              iy2=jy1-31
              ix3=jx1
              iy3=jy1-38
              ix5=jx1+18
              iy5=jy2+15
              ix4=jx1+22
              iy4=jy2+31
              do 1 ii=1,npage-1
                if(ix.ge.paginc(ii,iplot).and.ix.le.paginc(ii+1,iplot))
     ~                then
                   write(ii+12,6) jx1-10,jy1,ix1,iy1,ix2,iy2,ix3,iy3
                   write(ii+12,6) jx1+10,jy2,ix5,iy5,ix4,iy4,ix3,iy3
                endif
 1           continue
             call newarb(paginc,iplot,ix5,iy5,jx1+10,jy2,ix)
          else 
             if(hbpsn(i,2).lt.hbpsn(i,4)) then
C     CASE iv
                ix1=jx1-18
                iy1=jy1+15
                ix2=jx1-22
                iy2=jy1+31
                ix3=jx1
                iy3=jy1+38
                ix5=jx1+18
                iy5=jy2-15
                ix4=jx1+22
                iy4=jy2-31
             do 2 ii=1,npage-1
                if(ix.ge.paginc(ii,iplot).and.ix.le.paginc(ii+1,iplot))
     ~               then
                   write(ii+12,6) jx1-10,jy1,ix1,iy1,ix2,iy2,ix3,iy3
                   write(ii+12,6) jx1+10,jy2,ix5,iy5,ix4,iy4,ix3,iy3
                endif
 2           continue
             call newarb(paginc,iplot,ix5,iy5,jx1+10,jy2,ix)
          endif
       endif
      else if(j.eq.4) then
C     ALPHA HELIX
         if(hbpsn(i,2).gt.hbpsn(i,4)) then
C     CASE i
            ix1=jx1-22
            iy1=jy1
            ix2=jx1-19
            iy2=jy1-25
            if(ispace.eq.0) then
               ix3a=jx1
               ix3b=ix3a
               iy3a=(jy1+jy2)/2
               iy3b=iy3a
            else
               ix3a=jx1-10
               ix3b=jx1+8
               iy3a=6+(jy1+jy2)/2 
               iy3b=-7+(jy1+jy2)/2
            endif
            ix5=jx1+22
            iy5=jy2
            ix4=jx1+19
            iy4=jy2+25
             do 3 ii=1,npage-1
                if(ix.ge.paginc(ii,iplot).and.ix.le.paginc(ii+1,iplot)) 
     ~               then
                   write(ii+12,6) jx1-10,jy1,ix1,iy1,ix2,iy2,
     ~                  ix3a,iy3a
                   write(ii+12,6) jx1+10,jy2,ix5,iy5,ix4,iy4,
     ~                  ix3b,iy3b
                endif
 3           continue
             call newarb(paginc,iplot,ix5,iy5,jx1+10,jy2,ix)
          else if(hbpsn(i,2).lt.hbpsn(i,4)) then
C     CASE ii           
             ix1=jx1-22
             iy1=jy1
             ix2=jx1-19
             iy2=jy1+25
             if(ispace.eq.0) then
                ix3a=jx1
                ix3b=ix3a
                iy3a=(jy1+jy2)/2
                iy3b=iy3a
               else
                ix3a=jx1-10
                ix3b=jx1+8
                iy3a=-6+(jy1+jy2)/2 
                iy3b=7+(jy1+jy2)/2
             endif
             ix5=jx1+22
             iy5=jy2
             ix4=jx1+19
             iy4=jy2-25
             do 4 ii=1,npage-1
                if(ix.ge.paginc(ii,iplot).and.ix.le.paginc(ii+1,iplot)) 
     ~               then
                   write(ii+12,6) jx1-10,jy1,ix1,iy1,ix2,iy2,
     ~                  ix3a,iy3a
                   write(ii+12,6) jx1+10,jy2,ix5,iy5,ix4,iy4,
     ~                  ix3b,iy3b
                endif
 4           continue
             call newarb(paginc,iplot,ix5,iy5,jx1+10,jy2,ix)
          endif
       else
          if(jy1.gt.jy2) then
             idel=25
             ix1=jx1-15
             iy1=jy1-15
             ix2=jx1-idel
             iy2=jy1-30
             ix3=jx1-10
             iy3=jy2
             do 5 ij=1,npage-1
                if(ix.ge.paginc(ij,iplot).and.ix.le.paginc(ij+1,iplot))
     ~               write(ij+12,6) jx1-10,jy1,ix1,iy1,ix2,iy2,ix3,iy3
 5           continue
             call newarb(paginc,iplot,ix2,iy2,ix3,iy3,ix)
          else
             idel=25
             ix1=jx1-15
             iy1=jy2-15
             ix2=jx1-idel
             iy2=jy2-30
             ix3=jx1-10
             iy3=jy1
             do 7 ij=1,npage-1
                if(ix.ge.paginc(ij,iplot).and.ix.le.paginc(ij+1,iplot))
     ~               write(ij+12,6) jx1-10,jy2,ix1,iy1,ix2,iy2,ix3,iy3
 7           continue
             call newarb(paginc,iplot,ix2,iy2,ix3,iy3,ix)
          endif
       endif
       
 6     format(2(i5,1x),'moveto ',6(i5,1x),'curveto stroke %')
       return
       end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine strhb(hbpsn,i,ihb,hbond,ix1,iy1,ix2,iy2
     1      ,iplot,maxlen)
C     DRAWS HYDROGEN BONDS FOR STRANDS
        include 'p_hera.par'
C        
        integer hbond(mres*4,2),i,iplot,ix1,ix2,iy1,iy2,ihb,
     ~       hbpsn(mres*4,4),dx,dy,j,jx1,jx2,jy1,jy2,maxlen
     ~       ,nplot,xcoord,ycoord
        data dx/9/,dy/4/
        jx1=xcoord(ix1)
        jy1=ycoord(iy1,iplot,maxlen)
        jx2=xcoord(ix2)
        jy2=ycoord(iy2,iplot,maxlen)
        nplot=1
        do 10 j=1,ihb
           if(nplot.eq.0) go to 11
           if((hbond(j,1).eq.hbond(i,2).and.hbond(j,2).eq.
     ~          hbond(i,1)).or.jy1.eq.jy2) then
c     I.E. ANTIPARALLEL SHEET
              nplot=0
              if(hbpsn(i,3).gt.hbpsn(i,1)) then
                 call newarr(jx1+dx,jy1+dy,jx2-dx,jy2+dy)
              else
                 call newarr(jx1-dx,jy1-dy,jx2+dx,jy2-dy)
              endif
           endif
 10     continue
 11     continue
        if(nplot.eq.1) then
           if(hbpsn(i,3).gt.hbpsn(i,1)) then
              if(hbpsn(i,4).gt.hbpsn(i,2)) then
                 call newarr(jx1+dx,jy1+dy,jx2-dx,jy2-dy)
              else
                 call newarr(jx1+dx,jy1-dy,jx2-dx,jy2+dy)
              endif
           else
              if(hbpsn(i,4).gt.hbpsn(i,2)) then
                 call newarr(jx1-dx,jy1+dy,jx2+dx,jy2-dy)
              else
                 call newarr(jx1-dx,jy1-dy,jx2+dx,jy2+dy)
              endif
           endif
        endif
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine plotlp(istpl,lp,nlrow,n1,n2a,n2b,n3,n4,n5,
     ~       arr1,arr2a,arr2b,arr3,arr4,arr5,clmn,mcol,top,bot)
        include 'p_hera.par'
C     DECIDE WHICH KIND OF LOOP TO DRAW FOR EACH RESIDUE PAIR TO BE
C     JOINED UP
        integer arr1(lmax,7),arr2a(lmax,7),arr2b(lmax,7),arr3(lmax,7),
     ~       arr4(lmax,7),arr5(lmax,7),bot(sizex,10),clmn(sizex,2),i,
     ~       istpl(lmax,2),ixa,ixb,iya,iyb,j,k,lp(lmax,4),mcol,n,na,nb,
     ~       n1,n2a,n2b,n3,n4,n5,nlrow,top(sizex,10)
        logical tba,tbb,bba,bbb

cc    TBA==1ST RESIDUE AT TOP, TBB==LAST RESIDUE AT TOP
CC    BBA==1ST RESIDUE AT BOTTOM, BBB==LAST RESIDUE AT BOTTOM
C     
        do 1000 i=1,70
           do 999 j=1,7
              arr1(i,j)=0
              arr2a(i,j)=0
              arr2b(i,j)=0
              arr3(i,j)=0
              arr4(i,j)=0
              arr5(i,j)=0
 999       continue
 1000   continue
        do 100 n=2,nlrow-1
           na=istpl(n,1)
           nb=istpl(n,2)
           ixa=lp(n,1)
           iya=lp(n,2)
           ixb=lp(n,3)
           iyb=lp(n,4)
           tba=.false.
           tbb=.false.
           bba=.false.
           bbb=.false.
           if(ixa.eq.ixb) then
              call loop(ixa,iya,ixb,iyb,arr5,n5,1,na,nb)
           else
              do 5 k=1,10
                 if(top(ixa,k).eq.na) tba=.true.
                 if(top(ixb,k).eq.nb) tbb=.true.
                 if(bot(ixa,k).eq.na) bba=.true.
                 if(bot(ixb,k).eq.nb) bbb=.true.
 5            continue
             if(bba) then
                 if(bbb) then
                    call loop(ixa,iya,ixb,iyb,arr3,n3,1,na,nb)
                 else 
                    if(tbb) then
                       if(abs(ixa-ixb).eq.2.and.
     ~                      clmn((ixa+ixb)/2,2).eq.0) then
                          call loop(ixb,iyb,ixa,iya,arr4,n4,2,nb,na)
                       else
                          if(ixa+ixb.lt.mcol-1) then 
                             call loop(ixa,iya,ixb,iyb,arr2a,n2a,1,
     ~                            na,nb)
                          else
                             call loop(ixa,iya,ixb,iyb,arr2b,n2b,1,
     ~                            na,nb)
                          endif    
                       endif
                    endif
                 endif
              else
                 if(tba) then
                    if(tbb) then
                       call loop(ixa,iya,ixb,iyb,arr1,n1,1,na,nb)
                    else
                       if(bbb) then
                          if(abs(ixa-ixb).eq.2.and.
     ~                         clmn((ixa+ixb)/2,2).eq.0) then
                             call loop(ixa,iya,ixb,iyb,arr4,n4,1,
     ~                            na,nb)
                          else
                             if(ixa+ixb.lt.mcol-1) then 
                                call loop(ixb,iyb,ixa,iya,arr2a,n2a,2,
     ~                               nb,na)
                             else
                                call loop(ixb,iyb,ixa,iya,arr2b,n2b,2,
     ~                               nb,na)
                             endif    
                          endif
                       endif
                    endif
                 endif
              endif
           endif
 100    continue
        call lsrt13(arr3,n3)
        call lsrt13(arr1,n1)
        call lsrt2a(arr2a,n2a)
        call lsrt2b(arr2b,n2b)
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine ldraw(arr1,arr2a,arr2b,arr3,arr4,arr5,n1,n2a,n2b
     ~   ,n3,n4,n5,narray,iarr,mcol,lp,clmn,nlrow,maxlen,iplot,
     ~       paginc,icol,m1,m2,delxl,delyt,delyb)
        include 'p_hera.par'

C     DRAWS N AND C TERMINII AND CALLS THE OTHER LOOP DRAWING ROUTINES
        integer arr1(lmax,7),arr2a(lmax,7),arr2b(lmax,7),arr3(lmax,7),
     ~       arr4(lmax,7),arr5(lmax,7),clmn(sizex,2),delxl,delyb,
     ~       delyt,i,iarr(sizex,2),icol,ij,iplot,ix,ixn,iy,iya,iyn,
     ~       lower(50,5),lp(lmax,4),maxlen,mcol,m1,m2,n1,n2a,n2b,
     ~       n3,n4,n5,narray(sizex,sizey),nlrow,paginc(npage,3),
     ~       upper(50,5),xcoord,ycoord,ymax(3),ymin(3)
        character cap
        data ymax/60,50,60/,ymin/1,14,0/
c
        do 5 i=scfile,scfil5
           write(i,1)
           write(i,2)
           if(.not.bw) call colour(12,i)
 5      continue
1       format(' /Helvetica-Bold findfont %')
2       format(' [6 0 0 12 0 0] makefont setfont %')

c     DRAW N- AND C-TERMINI
        do 50 i=1,2
           if(i.eq.1) then
              ixn=lp(1,3)
              iyn=lp(1,4)
              cap='N'
           else
              ixn=lp(nlrow,1)
              iyn=lp(nlrow,2)
              cap='C'
           endif
           ix=xcoord(ixn)
           iy=ycoord(iyn,iplot,maxlen)
           if(narray(ixn,iyn+1).eq.0) then
              if(iyn.lt.sizey - 1) then
                 if(narray(ixn,iyn+2).eq.0) then
                    iya=ycoord(iyn+2,iplot,maxlen)
                    narray(ixn,iyn+2)=-101
                 else
                    iya=ycoord(iyn+1,iplot,maxlen)
                    narray(ixn,iyn+1)=-101
                 endif
              else
                 if(iyn.ne.sizey) then
                    iya=ycoord(sizey,iplot,maxlen)
                 else
                    iya=ycoord(sizey + 1,iplot,maxlen)
                 endif
              endif
              do 7 ij=1,npage-1
                 if(ixn.ge.paginc(ij,iplot).and.ixn.le.paginc
     ~                (ij+1,iplot)) write(ij+12,13) ix,iy+10,ix,iya-10
 7            continue
           else
              if(iyn.gt.2) then
                 if(narray(ixn,iyn-2).eq.0) then
                    iya=ycoord(iyn-2,iplot,maxlen)
                    narray(ixn,iyn-2)=-101
                 else
                    iya=ycoord(iyn-1,iplot,maxlen)
                    narray(ixn,iyn-1)=-101
                 endif
              else
                 iya=ycoord(0,1,maxlen)
              endif
              do 8 ij=1,npage-1
                 if(ixn.ge.paginc(ij,iplot).and.ixn.le.paginc
     ~                (ij+1,iplot)) write(ij+12,13) ix,iy-10,ix,iya+10
 8            continue
           endif
           do 11 ij=1,npage-1
              if(ixn.ge.paginc(ij,iplot).and.ixn.le.paginc
     ~             (ij+1,iplot)) write(ij+12,10) ix,cap,iya-3,cap
              if(ixn.ge.paginc(ij,iplot).and.ixn.le.paginc
     ~             (ij+1,iplot)) write(ij+12,14) ix,iya
 11        continue
           if(i.eq.1) then
              call newarb(paginc,iplot,ix,iya,ix,(iy+iya)/2,iyn)
           else
              call newarb(paginc,iplot,ix,iy,ix,(iy+iya)/2,iyn)
           endif
 50     continue
c     
        call lp5(arr5,n5,maxlen,iplot,paginc)
        call lp4(paginc,arr4,n4,maxlen,iplot)
        call lp3(paginc,arr3,n3,narray,m1,lower,clmn,maxlen,iplot,
     ~       delyb,ymin)
        call lp1(paginc,arr1,n1,narray,m2,upper,clmn,maxlen,iplot,
     ~       delyt,ymax)
        call lp2a(paginc,arr2a,n2a,narray,m2,m1,upper,lower,delxl,
     ~       clmn,maxlen,iplot,delyt,delyb,ymin,ymax)
        call lp2b(paginc,arr2b,n2b,narray,m2,m1,upper,lower,
     ~       mcol,clmn,maxlen,iplot,icol,delyt,delyb,ymin,ymax)

 10     format(i5,1x,'(',a1,') stringwidth pop 2 div sub ',i5,1x,
     ~       'moveto (',a1,') show %')
 13     format(2(i5,1x),'moveto ',2(i5,1x),'lineto stroke %')
 14     format(2(i5,1x),'moveto box % ')
        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine loop(ixa,iya,ixb,iyb,arr,n,k,na,nb)
        include 'p_hera.par'
        integer arr(lmax,7),ixa,ixb,iya,iyb,k,n,na,nb
        n=n+1
        arr(n,1)=ixa
        arr(n,2)=iya
        arr(n,3)=ixb
        arr(n,4)=iyb
        arr(n,5)=k
        arr(n,6)=na
        arr(n,7)=nb
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine lp5(arr,n,maxlen,iplot,paginc)
        include 'p_hera.par'
C        SINGLE LINES BETWEEN RESIDUES
        integer arr(lmax,7),dy1,i,ij,ixa,iya,iyb,iplot,maxlen,
     ~       n,xcoord,ycoord,paginc(npage,3)
            data dy1/9/
        do 10 i=1,n
           ixa=xcoord(arr(i,1))
           iya=ycoord(arr(i,2),iplot,maxlen)
           iyb=ycoord(arr(i,4),iplot,maxlen)
             if(iya.lt.iyb) then
              call newarb(paginc,iplot,ixa,iya+dy1,ixa,(iya+iyb)/2,0)
              do 4 ij=1,npage-1
                 if(arr(i,1).ge.paginc(ij,iplot).and.arr(i,1).le.paginc
     ~                (ij+1,iplot)) 
     ~                write(ij+12,13) ixa,iya+dy1,ixa,iyb-dy1
 4            continue
           else
              call newarb(paginc,iplot,ixa,iya-dy1,ixa,(iya+iyb)/2,0)
              do 5 ij=1,npage-1
                 if(arr(i,1).ge.paginc(ij,iplot).and.arr(i,1).le.paginc
     ~                (ij+1,iplot)) 
     ~             write(ij+12,13) ixa,iya-dy1,ixa,iyb+dy1
 5            continue
           endif
 10     continue
        
 13     format(2(i5,1x),'moveto ',2(i5,1x),'lineto stroke %')     
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine lp4(paginc,arr,n,maxlen,iplot)
        include 'p_hera.par'
C     LOOPS BETWEEN RESIDUES IN ADJACENT COLUMNS
        integer arr(lmax,7),i,ij,iplot,ixa,ixb,ixe,ix0,ix1,ix2,ix3,
     ~       ix4,ixm,iya,iyb,iye,iy0,iy1,iy2,iy3,iy4,jx1,jx2,jx3,
     ~       jx4,jy1,jy2,jy3,jy4,maxlen,n,xcoord,ycoord,dy1,
     ~       paginc(npage,3)
        data dy1/12/

        do 10 i=1,n
          ixa=arr(i,1)
          ixb=arr(i,3)
          iya=arr(i,2)
          iyb=arr(i,4)
          ix0=xcoord(ixa)
          iy0=ycoord(iya,iplot,maxlen)
          ixe=xcoord(ixb)
          iye=ycoord(iyb,iplot,maxlen)
          
          iy1=iy0+dy1
          iy2=iy0+20
          iy4=iy0+20
          iy3=iy0+dy1
          jy1=iye-dy1
          jy2=iye-20
          jy4=iye-20
          jy3=jy1
          ix1=ix0
          jx3=ixe
          ixm=int((ix0+ixe)/2)
          ix3=ixm
          jx1=ixm
          if(ixa.lt.ixb) then
            ix2=ix0+8
            ix4=ix0+16
            jx2=ixe-16
            jx4=ixe-8
          else
            ix2=ix0-8
            ix4=ix0-16
            jx2=ixe+16
            jx4=ixe+8
          endif
          if(arr(i,5).eq.2) then
             call newarb(paginc,iplot,ix3,jy1,ix3,(jy1+iy3)/2,0)
          else
             call newarb(paginc,iplot,ix3,iy3,ix3,(jy1+iy3)/2,0)
          endif        
          do 6 ij=1,npage-1
             write(ij+12,15) ix1,iy1,ix2,iy2,ix4,iy4,ix3,iy3
             write(ij+12,13) ix3,iy3
             write(ij+12,14) jx1,jy1
             write(ij+12,15) jx3,jy3,jx4,jy4,jx2,jy2,jx1,jy1
 6        continue
10        continue

13      format(2(i5,1x),'moveto %')
14        format(2(i5,1x),'lineto stroke %')
15        format(2(i5,1x),'moveto ',6(i5,1x),'curveto stroke %')
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        subroutine lp1(paginc,arr,n,narray,m2,upper,clmn,maxlen,iplot,
     ~       delyt,ymax)
        include 'p_hera.par'
C     LOOPS BETWEEN TWO RESIDUES AT TOPS OF COLUMNS
        integer arr(lmax,7),clmn(sizex,2),delyt,dy1,dy2,free,i,iplot,
     ~       iss,
     ~       ixa,ixb,ixe,ix0,ix1,ix2,ix3,iy,iya,iyb,iye,iy0,iy1,iy2,iy3,
     ~       j,jx,jx1,jx2,k,kn,maxlen,m2,n,na,nb,narray(sizex,sizey),
     ~       paginc(npage,3),upper(50,5),xcoord,ycoord,ymax(3)
        data dy1/9/,dy2/18/

        do 100 k=1,n
          na=arr(k,6)
          nb=arr(k,7)
          ixa=arr(k,1)
          ixb=arr(k,3)
          iya=arr(k,2)
          iyb=arr(k,4)
          ix0=xcoord(ixa)
          iy0=ycoord(iya,iplot,maxlen)
          ixe=xcoord(ixb)
          iye=ycoord(iyb,iplot,maxlen)
          iy=max0(iya,iyb)
          jx1=min0(ixa,ixb)
          jx2=max0(ixa,ixb)
          do 5 i=1,m2
             if(upper(i,4).eq.k.and.upper(i,5).eq.1) then
                iss=upper(i,1)
                go to 20
             endif
 5        continue
 6        continue
          do 10 j=iy,sizey
             free=1
             do 9 i=jx1,jx2
                kn=narray(i,j)
                if(kn.gt.0.or.(kn.lt.-5.and.mod(kn,100).ne.0)) free=0
 9           continue
         if(free.eq.1) then
            iss=j
            go to 11
         endif
 10   continue
 11   continue
c     
      do 12 i=1,m2
         if(upper(i,1).eq.iss) then
            if(jx2.lt.upper(i,2)) then
               go to 12
            else if(jx1.gt.upper(i,3)) then
               go to 12
            else
               iy=iss+1
               iss=iy
               if(iy.lt.sizey) go to 6
            endif
         endif
 12   continue
      m2=m2+1
      if(iss.gt.ymax(iplot)) iss=ymax(iplot)
      upper(m2,1)=iss
      upper(m2,2)=jx1
      upper(m2,3)=jx2
      upper(m2,4)=k
      upper(m2,5)=1
c        
 20   if(clmn(ixa,2).eq.na) then
         ix1=ix0
      else
         call lroute(narray,iss,iya+1,ixa,jx,1)
         ix1=xcoord(jx)
      endif
      iy1=iy0+dy1
      ix2=(ix0+ixe)/2
      if(iss.lt.ymax(iplot)) then
         iy2=ycoord(iss,iplot,maxlen)-2
      else
         do 22 i=1,m2-1
            if(upper(i,1).eq.iss.and.jx2.ge.upper(i,2).and.
     ~           jx1.le.upper(i,3)) then
               delyt=delyt+3
               go to 23
            endif
 22      continue
 23      iy2=ycoord(iss,iplot,maxlen)+delyt
      endif
      if(clmn(ixb,2).eq.nb) then
         ix3=ixe
      else
         call lroute(narray,iss,iyb+1,ixb,jx,1)
         ix3=xcoord(jx)
      endif  
      iy3=iye+dy2
      do 106 i=1,npage-1
         write(i+12,13) ix0,iy1
         if(ix1.ne.ix0) write(i+12,15) ix0,iy0+dy2,ix1,iy0+dy2
         write(i+12,15) ix1,iy2,ix3,iy2
         if(ix3.ne.ixe) write(i+12,15) ix3,iy3,ixe,iy3
         write(i+12,14) ixe,iye+dy1
 106     continue
      call newarb(paginc,iplot,ix1,iy1,ix1,(iy1+iy2)/2,0)
      call newarb(paginc,iplot,ix1,iy2,(ix1+ix3)/2,iy2,0)
      call newarb(paginc,iplot,ix3,iy2,ix3,(iy2+iye+dy1)/2,0)
 100  continue
      
 13   format(2(i5,1x),'moveto %')
 15   format(2(i5,1x),'lineto ',2(i5,1x),'lineto %')
 14   format(2(i5,1x),'lineto stroke %')
      
      return 
      end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine lp2a(paginc,arr,n,narray,m2,m1,upper,lower,
     ~     delxl,clmn,maxlen,iplot,delyt,delyb,ymin,ymax)
      include 'p_hera.par'
C     LOOPS TO LEFT OF DIAGRAM
      integer arr(lmax,7),clmn(sizex,2),delxl,delyb,delyt,dy1,dy2,
     ~     free,i,imax,imin,iplot,ixa,ixb,ixc,
     ~     ixd,ixe,ix0,iya,iyb,iye,iy0,j,jmax,jmin,jx,k,kn,
     ~     lower(50,5),maxlen,m1,m2,n,na,narray(sizex,sizey),nb,
     ~     paginc(npage,3),upper(50,5),xcoord,ycoord,ymax(3),ymin(3)
      data dy1/9/,dy2/18/
      
      do 100 k=1,n
         na=arr(k,6)
         nb=arr(k,7)
         ixa=arr(k,1)
         ixb=arr(k,3)
         iya=arr(k,2)
         iyb=arr(k,4)
         ix0=xcoord(ixa)
         iy0=ycoord(iya,iplot,maxlen)
         ixe=xcoord(ixb)
         iye=ycoord(iyb,iplot,maxlen)
         delxl=delxl+5
c     
         do 3 i=1,m2
            if(upper(i,4).eq.k.and.upper(i,5).eq.20) then
               jmax=upper(i,1)
               go to 10
            endif
 3       continue
c
         jmax=iyb
 6       do 101 j=jmax,sizey
            free=1
            do 99 i=1,ixb
               kn=narray(i,j)
               if(kn.gt.0.or.(kn.lt.-5.and.mod(kn,100).ne.0)) free=0
 99         continue
            if(free.eq.1) then
               jmax=j
               go to 102
            endif
 101     continue
 102     continue
         
         do 8 i=1,m2
            if(upper(i,1).eq.jmax) then
               if(upper(i,2).gt.ixb) then
                  go to 8
               else
                  jmax=jmax+1
                  go to 6
               endif
            endif
 8       continue   
         m2=m2+1
         if(jmax.gt.ymax(iplot)) jmax=ymax(iplot)
         upper(m2,1)=jmax
         upper(m2,2)=1
         upper(m2,3)=ixb
         upper(m2,4)=k
         upper(m2,5)=20
c     
 10      if(jmax.lt.ymax(iplot)) then
            imax=ycoord(jmax,iplot,maxlen)-2
         else
            do 22 i=1,m2-1
               if(upper(i,1).eq.jmax.and.upper(i,2).le.ixb) then
                  delyt=delyt+3
                  go to 23
               endif
 22         continue
 23         imax=ycoord(jmax,iplot,maxlen)+delyt
         endif
         do 12 i=1,m1
            if(lower(i,4).eq.k.and.lower(i,5).eq.20) then
               jmin=lower(i,1)
               go to 20
            endif
 12      continue

         jmin=iya
 16      do 111 j=jmin,1,-1
            free=1
            do 109 i=1,ixa
               if(narray(i,j).gt.1.or.narray(i,j).lt.-100) free=0
 109        continue
            if(free.eq.1) then
               jmin=j
               go to 112
            endif
 111     continue
 112     continue
         
         do 18 i=1,m1
            if(lower(i,1).eq.jmin) then
               if(lower(i,2).gt.ixa) then
                  go to 18
               else
                  jmin=jmin-1
                  if(jmin.ge.0) go to 16
               endif
            endif
 18      continue   
         m1=m1+1
         if(jmin.lt.ymin(iplot)) jmin=ymin(iplot)
         lower(m1,1)=jmin
         lower(m1,2)=1
         lower(m1,3)=ixa
         lower(m1,4)=k
         lower(m1,5)=20
c     
 20      if(jmin.gt.ymin(iplot)) then
            imin=ycoord(jmin,iplot,maxlen)        
         else
            do 32 i=1,m1-1
               if(lower(i,1).eq.jmin.and.lower(i,2).le.ixa) then
                  delyb=delyb+3
                  go to 33
               endif
 32         continue
 33         imin=ycoord(jmin,iplot,maxlen)-delyb
         endif
         if(clmn(ixa,1).eq.na) then
            ixc=ix0
         else
            call lroute(narray,jmin,iya-1,ixa,jx,-1)
            ixc=xcoord(jx)
         endif
         if(clmn(ixb,2).eq.nb) then
            ixd=ixe
         else
            call lroute(narray,jmax,iyb+1,ixb,jx,1)
            ixd=xcoord(jx)
         endif
         if(arr(k,5).eq.1) then
            call newarb(paginc,iplot,ix0,iy0,ix0,(iy0+imin)/2,0)
            call newarb(paginc,iplot,ix0,imin,(ix0+30)/2,imin,0)
            call newarb(paginc,iplot,30-delxl,imin,30-delxl,
     ~           (imin+imax)/2,0)
            call newarb(paginc,iplot,30-delxl,imax,(30+ixe)/2,imax,0
     ~        )
            call newarb(paginc,iplot,ixe,imax,ixe,(iye+imax)/2,0)
         else
            call newarb(paginc,iplot,ix0,imin,ix0,(iy0+imin)/2,0)
            call newarb(paginc,iplot,30,imin,(ix0+30)/2,imin,0)
            call newarb(paginc,iplot,30-delxl,imax,30-delxl,
     ~           (imin+imax)/2,0)
            call newarb(paginc,iplot,ixe,imax,(30+ixe)/2,imax,0)
            call newarb(paginc,iplot,ixe,iye,ixe,(iye+imax)/2,0)
         endif
         do 105 i=1,npage-1
            write(i+12,13) ix0,iy0-dy1
            if(ixc.ne.ix0) write(i+12,15) ix0,iy0-dy2,ixc,iy0-dy2
            write(i+12,15) ixc,imin,30-delxl,imin
            write(i+12,15) 30-delxl,imax,ixd,imax
            if(ixd.ne.ixe) write(i+12,15) ixd,iye+dy2,ixe,iye+dy2
            write(i+12,14) ixe,iye+dy1
 105     continue
 100  continue
      
 13   format(2(i5,1x),'moveto %')
 14   format(2(i5,1x),'lineto stroke %')
 15   format(2(i5,1x),'lineto ',2(i5,1x),'lineto %')
      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine lp2b(paginc,arr,n,narray,m2,m1,upper,lower,
     ~     mcol,clmn,maxlen,iplot,icol,delyt,delyb,ymin,ymax)
      include 'p_hera.par'
C     LOOPS ROUND TO RIGHT OF DIAGRAM
      integer arr(lmax,7),clmn(sizex,2),delxr,delyb,delyt,dy1,
     ~     dy2,free,i,icol,imax,imin,iplot,ixa,
     ~     ixb,ixc,ixd,ixe,ix0,iya,iyb,iye,iy0,j,
     ~     jmax,jmin,jx,k,kn,lower(50,5),maxlen,mcol,m1,m2,n,
     ~     na,narray(sizex,sizey),nb,paginc(npage,3),upper(50,5),
     ~     xcoord,ycoord,ymax(3),ymin(3)
      data delxr/0/,dy1/9/,dy2/18/
      
      do 100 k=1,n
         na=arr(k,6)
         nb=arr(k,7)
         ixa=arr(k,1)
         iya=arr(k,2)
         ixb=arr(k,3)
         iyb=arr(k,4)
         ix0=xcoord(ixa)
         iy0=ycoord(iya,iplot,maxlen)
         ixe=xcoord(ixb)
         iye=ycoord(iyb,iplot,maxlen)
         delxr=delxr+5
c     
         do 3 i=1,m2
            if(upper(i,4).eq.k.and.upper(i,5).eq.25) then
               jmax=upper(i,1)
               go to 10
            endif
 3       continue
         
         jmax=iyb
 6       do 101 j=jmax,sizey
            free=1
            do 99 i=ixb,mcol
               kn=narray(i,j)
               if(kn.gt.0.or.(kn.lt.-5.and.mod(kn,100).ne.0)) free=0
 99         continue
            if(free.eq.1) then
               jmax=j
               go to 102
            endif
 101     continue
 102     continue
c     
         do 8 i=1,m2
            if(upper(i,1).eq.jmax) then
               if(upper(i,3).lt.ixb) then
                  go to 8
               else
                  jmax=jmax+1
                  go to 6
               endif
            endif
 8       continue   
         m2=m2+1
         if(jmax.gt.ymax(iplot)) jmax=ymax(iplot)
         upper(m2,1)=jmax
         upper(m2,2)=ixb
         upper(m2,3)=mcol
         upper(m2,4)=k
         upper(m2,5)=25
 10      if(jmax.lt.ymax(iplot)) then
            imax=ycoord(jmax,iplot,maxlen)-2
         else
            do 32 i=1,m2-1
               if(upper(i,1).eq.jmax.and.upper(i,3).ge.ixb) then
                  delyt=delyt+3
                  go to 33
               endif
 32         continue
 33         imax=ycoord(jmax,iplot,maxlen)+delyt
         endif
         
         do 12 i=1,m1
            if(lower(i,4).eq.k.and.lower(i,5).eq.25) then
               jmin=lower(i,1)
               go to 20
            endif
 12      continue
         jmin=iya
 16      do 111 j=jmin,1,-1
            free=1
            do 109 i=ixa,mcol
               kn=narray(i,j)
               if(kn.gt.0.or.(kn.lt.-5.and.mod(kn,100).ne.0)) free=0
 109        continue
            if(free.eq.1) then
               jmin=j
               go to 112
            endif
 111     continue
 112     continue
c     
         do 18 i=1,m1
            if(lower(i,1).eq.jmin) then
               if(lower(i,3).lt.ixa) then
                  go to 18
               else
                  jmin=jmin-1
                  go to 16
               endif
            endif
 18      continue   
         m1=m1+1
         if(jmin.lt.ymin(iplot)) jmin=ymin(iplot)
         lower(m1,1)=jmin
         lower(m1,2)=ixa
         lower(m1,3)=mcol
         lower(m1,4)=k
         lower(m1,5)=25
c     
 20      if(jmin.gt.ymin(iplot)) then
            imin=ycoord(jmin,iplot,maxlen)
         else
            do 22 i=1,m1-1
               if(lower(i,1).eq.jmin.and.lower(i,3).ge.ixa) then
                  delyb=delyb+3
                  go to 23
               endif
 22         continue
 23         imin=ycoord(jmin,iplot,maxlen)-delyb
         endif
         icol=xcoord(mcol+2)+delxr
         if(clmn(ixa,1).eq.na) then
            ixc=ix0
         else
            call lroute(narray,jmin,iya-1,ixa,jx,-1)
            ixc=xcoord(jx)
         endif
         if(clmn(ixb,2).eq.nb) then
            ixd=ixe
         else
            call lroute(narray,jmax,iyb+1,ixb,jx,1)
            ixd=xcoord(jx)
         endif
         
         if(arr(k,5).eq.1) then
            call newarb(paginc,iplot,ix0,iy0,ix0,(iy0+imin)/2,0)
            call newarb(paginc,iplot,ix0,imin,(ix0+icol)/2,imin,0)
            call newarb(paginc,iplot,icol,imin,icol,(imax+imin)/2,0)
            call newarb(paginc,iplot,icol,imax,(icol+ixe)/2,imax,0)
            call newarb(paginc,iplot,ixe,imax,ixe,(iye+imax)/2,0)
         else
            call newarb(paginc,iplot,ix0,imin,ix0,(iy0+imin)/2,0)
            call newarb(paginc,iplot,icol,imin,(ix0+icol)/2,imin,0)
            call newarb(paginc,iplot,icol,imax,icol,(imax+imin)/2,0)
            call newarb(paginc,iplot,ixe,imax,(icol+ixe)/2,imax,0)
            call newarb(paginc,iplot,ixe,iye,ixe,(iye+imax)/2,0)
         endif
         do 105 i=1,npage-1
            write(i+12,13) ix0,iy0-dy1
            if(ixc.ne.ix0) write(i+12,15) ix0,iy0-dy2,ixc,iy0-dy2
            write(i+12,15) ixc,imin,icol,imin
            write(i+12,15) icol,imax,ixd,imax
            if(ixd.ne.ixe) write(i+12,15) ixd,iye+dy2,ixe,iye+dy2
            write(i+12,14) ixe,iye+dy1
 105     continue
 100  continue
      
 13   format(2(i5,1x),'moveto %')
 14   format(2(i5,1x),'lineto stroke %')
 15   format(2(i5,1x),'lineto ',2(i5,1x),'lineto %')
      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine lp3(paginc,arr,n,narray,m1,lower,clmn,maxlen,iplot,
     ~     delyb,ymin)
      include 'p_hera.par'
C     LOOPS BETWEEN TWO RESIDUES AT BOTTOM OF COLUMNS
      integer arr(lmax,7),clmn(sizex,2),delyb,dy1,dy2,free,i,iplot,
     ~     iss,ixa,ixb,ix0,ix1,ix2,ix3,ixe,iy,iya,iyb,iye,
     ~     iy0,iy1,iy2,iy3,j,jx,jx1,jx2,k,kn,lower(50,5),maxlen,m1,n,
     ~     na,narray(sizex,sizey),nb,paginc(npage,3),xcoord,ycoord,
     ~     ymin(3)
      data dy1/9/,dy2/18/

      do 100 k=1,n
         iss=0
         na=arr(k,6)
         nb=arr(k,7)
         ixa=arr(k,1)
         ixb=arr(k,3)
         iya=arr(k,2)
         iyb=arr(k,4)
         ix0=xcoord(ixa)
         iy0=ycoord(iya,iplot,maxlen)
         ixe=xcoord(ixb)
         iye=ycoord(iyb,iplot,maxlen)
         iy=min0(iya,iyb)
         jx1=min0(ixa,ixb)
         jx2=max0(ixa,ixb)
         do 5 i=1,m1
            if(lower(i,4).eq.k.and.lower(i,5).eq.3) then
               iss=lower(i,1)
               go to 20
            endif
 5       continue
 6       continue
c     IF NO VALUE CALCULATED ALREADY CALCULATE ONE!
         do 10 j=iy,1,-1
            free=1
            do 9 i=jx1,jx2
               kn=narray(i,j)
               if(kn.gt.0.or.(kn.lt.-5.and.mod(kn,100).ne.0)) free=0
 9          continue
            if(free.eq.1) then
               iss=j
               go to 11
            endif
 10      continue
 11      continue
         do 12 i=1,m1
            if(lower(i,1).eq.iss) then
               if(jx2.lt.lower(i,2)) then
                  go to 12
               else if(jx1.gt.lower(i,3)) then
                  go to 12
               else
                  iy=iss-1
                  iss=iy
                  if(iy.ge.0.and.free.eq.1) go to 6
               endif
            endif
 12      continue
 20      if(iss.lt.ymin(iplot)) iss=ymin(iplot)
         m1=m1+1
         lower(m1,1)=iss
         lower(m1,2)=jx1
         lower(m1,3)=jx2
         lower(m1,4)=k
         lower(m1,5)=3
c     
         if(clmn(ixa,1).eq.na) then
            ix1=ix0
         else
            call lroute(narray,iss,iya-1,ixa,jx,-1)
            ix1=xcoord(jx)
         endif
         iy1=iy0-dy1
         ix2=(ix0+ixe)/2
         if(iss.gt.ymin(iplot)) then
            iy2=ycoord(iss,iplot,maxlen)+dy1
         else
            do 32 i=1,m1-1
               if(lower(i,1).eq.iss.and.jx2.ge.lower(i,2).and.
     ~              jx1.le.lower(i,3)) then
                  delyb=delyb+3
                  go to 33
               endif
 32         continue
 33         iy2=ycoord(iss,iplot,maxlen)+dy1-delyb
         endif
         if(clmn(ixb,1).eq.nb) then
            ix3=ixe
         else
            call lroute(narray,iss,iyb-1,ixb,jx,-1)
            ix3=xcoord(jx)
         endif
         iy3=iye-dy1
         do 106 i=1,npage-1
            write(i+12,13) ix0,iy1
            if(ix1.ne.ix0) write(i+12,15) ix0,iy0-dy2,ix1,iy0-dy2
            write(i+12,15) ix1,iy2,ix3,iy2
            if(ix3.ne.ixe) write(i+12,15) ix3,iye-dy2,ixe,iye-dy2
            write(i+12,14) ixe,iye-dy1
 106     continue
         call newarb(paginc,iplot,ix1,iy1,ix1,(iy1+iy2)/2,0)
         call newarb(paginc,iplot,ix1,iy2,(ix1+ix3)/2,iy2,0)
         call newarb(paginc,iplot,ix3,iy2,ix3,(iy2+iy3)/2,0)   
 100  continue
      
 13   format(2(i5,1x),'moveto %')
 14   format(2(i5,1x),'lineto stroke %')
 15   format(2(i5,1x),'lineto ',2(i5,1x),'lineto %')
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lroute(narray,iss,iy,ix,jx,sgn)
C     CALCULATES A MORE COMPLICATED ROUTE FOR LOOP WHERE DIRECT ROUTE
C     WOULD DRAW OVER EXISTING PARTS OF DIAGRAM
      include 'p_hera.par'
      integer i,iblock,ipath,iss,ix,iy,j,jx,narray(sizex,sizey),sgn
      
      iblock=0
      ipath=0
      jx=0
      do 10 j=iy,iss,sgn
         if(narray(ix,j).gt.1.or.narray(ix,j).le.-100)
     ~        iblock=iblock+1
 10   continue
      if(iblock.eq.0) then
           ipath=1
           jx=ix
           go to 30
        endif
        do 20 i=ix-1,ix+1,2
           iblock=0
           ipath=0
           do 15 j=iy,iss,sgn
              if(narray(i,j).gt.1.or.narray(i,j).le.-100)
     ~             iblock=iblock+1
 15        continue
           if(iblock.eq.0) then
              ipath=1
              jx=i
              go to 30
           endif
 20     continue
 30     continue
        if(ipath.eq.1) then
           do 40 j=iy,iss,sgn
              narray(jx,j)=narray(jx,j)-100
 40        continue
        endif
        return
        end
cccccccccccccccccccccccccccccccccccccccccccc
        subroutine newarr(jx1,jy1,jx2,jy2)
C     ARROW 
        include 'p_hera.par'
        integer ij,jx1,jx2,jy1,jy2

        do 4 ij=1,npage-1
           write(ij+12,5) jx1,jy1,jx2,jy2
 4      continue
5        format(' newpath',4(i5,1x),'arrow stroke %')
        return
        end

cccccccccccccccccccccccccccccccccccccccccccc        
        subroutine newarb(paginc,iplot,jx1,jy1,jx2,jy2,ix)
C     ARROW HEAD
        include 'p_hera.par'
        integer ij,iplot,ix,jx1,jy1,jx2,jy2,paginc(npage,3)
        
        write(scfile,5) jx1,jy1,jx2,jy2
        if(ix.eq.0) then
           do 2 ij=1,npage-1
              write(ij+12,5) jx1,jy1,jx2,jy2
 2         continue
        else
           do 3 ij=1,npage-1
              if(ix.ge.paginc(ij,iplot).and.ix.le.paginc
     ~             (ij+1,iplot)) write(ij+12,5) jx1,jy1,jx2,jy2
 3         continue
        endif
5        format(' newpath',4(i5,1x),'arhead stroke %')
        return
        end

cccccccccccccccccccccccccccccccccccccccccccc        
        subroutine psinit(nfil)
C     INITIALISE POSTSCRIPT
        include 'p_hera.par'
        integer nfil,i,icol
        character colno*2

        write(nfil,1)
 1      format('%!PS-Adobe-1.0'/
     &         '%%Creator: Gail Hutchinson'/
     &         '%%Title: HERA Diagram'/
     &         '%%EndComments')
        write(nfil,5)
5        format(' /arrowdict 14 dict def'/
     &         ' arrowdict begin'/
     &         '  /mtrx matrix def end'/
     &         ' /arrow {arrowdict begin'/
     &         '  /tipy exch def /tipx exch def'/
     &         '  /taily exch def /tailx exch def'/
     &         '  /dx tipx tailx sub def'/
     &         '  /dy tipy taily sub def'/
     &         '  /arrowlength dx dx mul dy dy mul add sqrt def'/
     &         '  /angle dy dx atan def'/
     &         '  /savematrix mtrx currentmatrix def'/
     &         '  tailx taily translate'/
     &         '  angle rotate'/
     &         '  0 0 moveto'/
     &         '  arrowlength 0 lineto'/
     &         '  -5 3 rlineto'/
     &         '  arrowlength 0 moveto'/
     &         '  -5 -3 rlineto'/
     &         '   savematrix setmatrix end} def'/
     &         ' 0.1 setlinewidth')
        write(nfil,6)
6        format(' /arheaddict 14 dict def'/
     &         ' arheaddict begin'/
     &         '  /mtrx matrix def'/
     &         '  end'/
     &         ' /arhead'/
     &         ' {arheaddict begin'/
     &         '  /tipy exch def /tipx exch def'/
     &         '  /taily exch def /tailx exch def'/
     &         '  /dx tipx tailx sub def'/
     &         '  /dy tipy taily sub def'/
     &         '  /arrowlength dx dx mul dy dy mul add sqrt def'/
     &         '  /angle dy dx atan def'/
     &         '  /savematrix mtrx currentmatrix def'/
     &         '  tailx taily translate'/
     &         '  angle rotate'/
     &         '  arrowlength 0 moveto'/
     &         '  -5 3 rlineto'/
     &         '  arrowlength 0 moveto'/
     &         '  -5 -3 rlineto'/
     &         '   savematrix setmatrix end} def'/)
               write(nfil,7)
7        format(' /box'/
     ~               ' {-9 -9 rmoveto '/
     ~         ' 0 18 rlineto '/
     ~               ' 18 0 rlineto '/
     ~         ' 0 -18 rlineto '/
     ~               ' -18 0 rlineto stroke }def')
        write (nfil,8)
 8      format(' /Col { setrgbcolor} bind def')
        do 10 icol=1,maxcol
           write(colno,'(i2)') icol
           if(colno(1:1).eq.' ') colno(1:1)='0'
           write(nfil,15) colno,(rgb(i,icol),i=1,3)
 15        format(' /Col',A2,'{ ',3f8.4,' setrgbcolor } def')
 10     continue
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     MISCELLANEOUS ROUTINES
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         subroutine conect(iord,carr,iconn,mc,iuse,istr,icol,
     ~       numstr,jjoin,freq,str)
C     DEFINES CONNECTIONS BETWEEN STRANDS IN SHEETS
        include 'p_hera.par'
        integer carr(nsmax,3),freq(nsmax),i,icol(nsmax),
     ~       iconn(nsmax,3),idone(nsmax),ii,ilen,iord,ip,istr,
     ~       iuse(nsmax),j,jjoin(3*nsmax,3),jlen,k,kb(nsmax,3)
     ~       ,ki,kk,m,mc,mconn,n,nconn,nopt,nr,numstr,
     ~       str(nsmax,2),strp
        
1        ilen=0
        kk=0
        mconn=0
        do 10 i=1,mc
           do 9 j=1,2
              if(carr(i,j).eq.iconn(iord,1)) then
                 n=strp(carr,i,j)
                 if(iuse(n).eq.0) then
                    kk=kk+1
                    kb(kk,1)=n
                    kb(kk,2)=carr(i,3)
                    kb(kk,3)=freq(n)
                 endif
              endif
 9         continue
 10     continue
c     1 STRAND ATTACHED TO CURRENT STRAND
        if(kk.eq.1) then
           iord=iord+1
           istr=istr+1 
           iconn(iord,1)=kb(1,1)
           icol(kb(1,1))=iord
           iuse(kb(1,1))=1
           if(istr.lt.numstr) go to 1
C     MORE THAN ONE STRAND ATTACHED TO CURRENT STRAND
        else if(kk.gt.1) then
           ki=kb(1,1)
C     CHECK WHETHER ANY OF THESE PARTNER STRANDS IS AT THE END OF A BRANCH
           nopt=0 
           do 11 i=1,kk
              if(kb(i,3).gt.1) then
                 nopt=nopt+1
                 ki=kb(i,1)
              endif
 11        continue
c     IF >1 OF THESE NOT AT END OF BRANCH SEE WHICH ONE TO USE!
           if(nopt.gt.1) then
              do 15 i=1,kk
                 do 5 k=1,nsmax
                    idone(k)=0
 5               continue
                 nconn=0
                 do 116 nr=1,3
                    m=kb(i,1)
 12                 do 14 k=1,mc
                       do 13 j=1,2
                          if(carr(k,j).eq.m.and.idone(k).eq.0) then
                             ii=strp(carr,k,j) 
                             if(iuse(ii).eq.0) then
                                m=ii
                                nconn=nconn+1
                                idone(k)=1
                                go to 12
                             endif
                          endif
 13                    continue
 14                 continue
 116             continue
                 jlen=str(m,2)-str(m,1)
                 if(nconn.gt.mconn.or.mconn.eq.0.or.
     ~                (nconn.eq.mconn.and.jlen.gt.ilen)) then
                    mconn=nconn
                    ki=kb(i,1)
                    ilen=str(m,2)-str(m,1)
                 endif
 15           continue
           endif
           iord=iord+1
           istr=istr+1
           iconn(iord,1)=ki
           icol(ki)=iord
           iuse(ki)=1
           if(istr.lt.numstr) go to 1
c     NO STRANDS ATTACHED TO CURRENT STRAND-
        else
C     FIRST LEVEL OF BRANCHING
           do 20 i=1,numstr
              if(iuse(i).eq.0) then
c     FIND A STRAND NOT CONNECTED AND FIND ITS PARTNER
                 ip=0
                 do 17 m=1,mc
                    do 16 n=1,2
                       if(carr(m,n).eq.i) ip=strp(carr,m,n)
                       if(ip.ne.0) then
                          if(iuse(ip).eq.1) then
c                             if(jjoin(2*ip-1,1).eq.0) then
c                                jjoin(2*ip-1,1)=i
c                             else 
c                                if(jjoin(2*ip,1).eq.0) then
c                                   jjoin(2*ip,1)=i
c                                endif
c                             endif              
                             if(jjoin(3*ip-2,1).eq.0) then
                                jjoin(3*ip-2,1)=i
                             else if(jjoin(3*ip-1,1).eq.0) then
                                jjoin(3*ip-1,1)=i
                             else if(jjoin(3*ip,1).eq.0) then
                                jjoin(3*ip,1)=i
                             else
                                write(*,*)'no jjoin 1'
                             endif
                             iuse(i)=2
                             istr=istr+1
                             if(istr.eq.numstr) then
                                go to 50
                             else
                                go to 19
                             endif
                          endif
                       endif
 16                 continue
 17              continue
 19              continue
              endif
 20        continue
c     SECOND LEVEL OF BRANCHING
           do 30 i=1,numstr
              if(iuse(i).eq.0) then
c     FIND A STRAND NOT CONNECTED AND FIND ITS PARTNER
                 ip=0
                 do 27 m=1,mc
                    do 26 n=1,2
                       if(carr(m,n).eq.i) ip=strp(carr,m,n)
                       if(ip.ne.0) then
                          if(iuse(ip).eq.2) then
c                             if(jjoin(2*ip-1,2).eq.0) then
c                                jjoin(2*ip-1,2)=i
c                             else 
c                                if(jjoin(2*ip,2).eq.0)
c     ~                               jjoin(2*ip,2)=i
c                             endif              
                             if(jjoin(3*ip-2,2).eq.0) then
                                jjoin(3*ip-2,2)=i
                             else if(jjoin(3*ip-1,2).eq.0) then
                                jjoin(3*ip-1,2)=i
                             else if(jjoin(3*ip,2).eq.0) then
                                jjoin(3*ip,2)=i
                             else
                                write(*,*)'no jjoin 2'
                             endif
                             iuse(i)=3
                             istr=istr+1
                             go to 29
                          else
                             if(iuse(ip).eq.3) then
c                                if(jjoin(2*ip-1,3).eq.0) then
c                                   jjoin(2*ip-1,3)=i
c                                else 
c                                   if(jjoin(2*ip,3).eq.0) 
c     ~                                  jjoin(2*ip,3)=i
c                                endif              
                                if(jjoin(3*ip-2,3).eq.0) then
                                   jjoin(3*ip-2,3)=i
                                else if(jjoin(3*ip-1,3).eq.0) then
                                   jjoin(3*ip-1,3)=i
                                else if(jjoin(3*ip,3).eq.0) then
                                   jjoin(3*ip,3)=i
                                else
                                   write(*,*)'no jjoin 3'
                                endif
                                iuse(i)=3
                                istr=istr+1
                                go to 29
                             endif
                          endif
                       endif
 26                 continue
 27              continue
              endif
 29           continue
 30        continue
           
C     THIRD LEVEL OF BRANCHING
      do 40 i=1,numstr
         if(iuse(i).eq.0) then
c     FIND A STRAND NOT CONNECTED AND FIND ITS PARTNER
            ip=0
            do 37 m=1,mc
               do 36 n=1,2
                  if(carr(m,n).eq.i) ip=strp(carr,m,n)
                  if(ip.ne.0) then
                     if(iuse(ip).eq.3) then
c                        if(jjoin(2*ip-1,3).eq.0) then
c                           jjoin(2*ip-1,3)=i
c                        else 
c                           if(jjoin(2*ip,3).eq.0) 
c     ~                          jjoin(2*ip,3)=i
c                        endif      
                        if(jjoin(3*ip-2,3).eq.0) then
                           jjoin(3*ip-2,3)=i
                        else if(jjoin(3*ip-1,3).eq.0) then
                           jjoin(3*ip-1,3)=i
                        else if(jjoin(3*ip,3).eq.0) then
                           jjoin(3*ip,3)=i
                        else
                           write(*,*)'no jjoin 4'
                        endif
                        iuse(i)=4
                        istr=istr+1
                        go to 39
                     endif
                  endif
 36            continue
 37         continue
         endif
 39      continue
 40   continue
      endif
 50   return
      end
c***********************************************************************
      subroutine readst(numres,iseqno,insert,jseqno,aa,resnam,sstr,sht,
     ~     ibp,hbp,hbe,nres,title,minres,maxres,option,access,chn,
     ~     chain,sstchk)
C     READ SECONDARY STRUCTURE FILE
      include 'p_hera.par'
      
      integer access,hbp(mres,4),i,ibp(mres,2),
     ~     iseqno(mres),j,jseqno(mres),m,minres,maxres,nres,
     ~     numres,option,sstchk,sht(mres)
      real hbe(mres,4)
      character aa(mres),insert(mres),sstr(mres),title*4,chn(mres),chain
      CHARACTER*3 resnam(mres)

      sstchk=0
      access=0
      numres=0
      nres=0
      do 200 i=1,mres
         iseqno(i)=0
         jseqno(i)=0
         insert(i)=' '
         aa(i)=' '
         resnam(i)=' '
         sstr(i)=' '
         sht(i)=0
         do 110 j=1,2
              ibp(i,j)=0
 110       continue
           do 120 j=1,4
              hbp(i,j)=0
              hbe(i,j)=0.0
 120       continue
 200    continue

        read(7,*,end=250)
        read(7,*,end=250)
        read(7,5,end=250) title
        read(7,10,end=250) numres
        read(7,*,end=250)
        read(7,*,end=250)
        read(7,*,end=250)
        do 20 m=1,numres
           read(7,15,end=21) iseqno(m),chn(m),jseqno(m),insert(m),
     ~          resnam(m), aa(m),sstr(m),sht(m),(ibp(m,i),i=1,2),
     ~          (hbp(m,i),hbe(m,i),i=1,4)
 20     continue
c 21     if(maxres.eq.0.or.option.ne.1) maxres=numres
c        if(minres.eq.0.or.option.ne.1) minres=1
c     CHANGE TO RUN ON MULTIPLE FILES
 21     MAXRES=NUMRES
        MINRES=1
        if(chain.ne.' ') then
           do 221 m=1,numres
              if(chn(m).eq.chain) then
                 minres=m
                 go to 222
              endif
 221       continue
 222       do 223 m=numres,1,-1
              if(chn(m).eq.chain) then
                 maxres=m
                 go to 224
              endif
 223       continue
 224       continue
        endif
        close(unit=7)
        go to 300
        
 5      format(1x,a4)
 10     format(19x,i4)
 15     format(1x,i4,1x,a1,i4,a1,3x,a3,1x,a1,5x,a1,2x,i3,7x,2(1x,i4),
     ~       40x,4(1x,i4,1x,f4.1))
 250    write(*,*) 'No data in secondary structure file'
        sstchk=1
 300    return
        end
c**************************************************************************
        subroutine srcha(n,clmn,hlxpsn,mindel,irow,istph,jsign,
     ~       ifind,nhrow)
C     HELIX PLACING 
        include 'p_hera.par'
        integer clmn(sizex,2),del,hlxpsn(mhel,3),i,ifind,j,ji,jx,
     ~       k,m,mindel(mhel),n,nhrow,irow(sizex,2),istph(mhel,2)
        character jsign(sizex,2)

        jx=hlxpsn(n,1)
        ji=hlxpsn(n,2)
        ifind=0
        mindel(n)=40
        do 8 k=1,2
           do 7 j=1,2
              do 6 i=2,sizex,2
                 if(clmn(i,j).ne.0.and.irow(i,j).eq.0) then
                    if(jsign(i,j).eq.'b'.and.k.eq.1) then
                       del=istph(n,k)-clmn(i,j)
                    else if(jsign(i,j).eq.'a'.and.k.eq.2) then
                       del=clmn(i,j)-istph(n,k)
                    else
                       del=1000
                    endif
                    if(del.gt.0.and.del.lt.mindel(n).and.(i.ne.jx.
     ~                   or.j.ne.ji)) then
                       do 5 m=1,nhrow
                          if(hlxpsn(m,1).eq.i.and.hlxpsn(m,2).eq.j.
     ~                         and.mindel(m).lt.del) then
                             mindel(n)=40
                             hlxpsn(n,1)=jx
                             hlxpsn(n,2)=ji
                             go to 6
                          endif
                          mindel(n)=del
                          hlxpsn(n,1)=i
                          hlxpsn(n,2)=j
                          hlxpsn(n,3)=k
                          ifind=1
 5                     continue
                    endif
                 endif
 6            continue
 7         continue
 8      continue
        return
        end

C**************************************************************************
C
C  SUBROUTINE WRIDAT  -  Write out the helix and strand data
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRIDAT(NARRAY,SIZEX,SIZEY,CHAIN,CHN,JSEQNO,RESNAM,
     -     INSERT,SHAPE,MRES,NRES,MINRES,MAXRES)

      INTEGER       MRES, SIZEX, SIZEY

      CHARACTER*1   CHAIN, CHN(MRES), CSHAPE, INSERT(MRES), LSHAPE,
     -              SHAPE(MRES)
      CHARACTER*3   RESNAM(MRES)
      CHARACTER*6   TYPE
      INTEGER       I, ISEQ, ISTART, J, JSEQNO(MRES), LSTSEQ, MINRES,
     -              MAXRES, NARRAY(SIZEX,SIZEY), NRES

C---- Open the output file
      OPEN(UNIT=21, FILE='hera.dat', STATUS='UNKNOWN',
     -     FORM='FORMATTED',ACCESS='SEQUENTIAL',ERR=900)

C---- Loop through the residues to write out
      WRITE(21,10) CHAIN, NRES, MINRES, MAXRES
 10   FORMAT(A1,1X,3I6)
      DO 200, I = MINRES, MAXRES
          WRITE(21,20) CHN(I), RESNAM(I), JSEQNO(I), INSERT(I)
 20       FORMAT(A1,A3,I4,A1)
 200  CONTINUE

C---- Loop through the array to identify all the strands and helices
      DO 800, I = 1, SIZEX

C----     Initialise last residue number
          LSTSEQ = 0
          LSHAPE = ' '
          ISTART = 0

C----     Loop down the y-column
          DO 600, J = 1, SIZEY + 1

C----         Initialise current shape
              CSHAPE = ' '

C----         Get this sequence position
              IF (J.LE.SIZEY) THEN
                  ISEQ = NARRAY(I,J)
                  IF (ISEQ.GT.0 .AND. ISEQ.LE.MRES) THEN
                      CSHAPE = SHAPE(ISEQ)
                   ENDIF
              ELSE
                  ISEQ = 0
              ENDIF

C----         If shape has changed, end previous secondary structure
C             element and start new one
              IF (CSHAPE.NE.LSHAPE) THEN

C----             If have a previous SSE, write out its details
                  IF (LSHAPE.NE.' ' .AND. ISTART.GT.0) THEN

C----                 Define SSE type
                      IF (LSHAPE.EQ.'C') THEN
                          TYPE = 'HELIX '
                      ELSE
                          TYPE = 'STRAND'
                      ENDIF

C----                 Write out the SSE details
                      WRITE(21,110) TYPE, NARRAY(I,ISTART),
     -                    NARRAY(I,J - 1), I, ISTART, J - 1
 110                  FORMAT(A6,2I5,3I4)
                  ENDIF

C----             If this is the start of a new SSE, store start
C                 position
                  IF (CSHAPE.NE.' ') THEN
                      ISTART = J
                  ELSE
                      ISTART = 0
                  ENDIF
              ENDIF

C----         Store current shape
              LSHAPE = CSHAPE
 600      CONTINUE
 800  CONTINUE

C---- Close the output file
      CLOSE(21)
      GO TO 999

C---- Fatal error
 900  CONTINUE
      PRINT*, '*** Error opening output file hera.dat'

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C*******************************************************************

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

c************************************************************************
      subroutine colour(ic,nfil)

      integer ic,nfil

      if(ic.lt.10) then
         write(nfil,5) ic
 5       format(' Col0',i1,' %')
      else
         write(nfil,15) ic
 15      format(' Col',i2,' %')
      endif

      return
      end
c***************************************************************************
      logical function logics(yesno)

      character yesno
      
      if(yesno.eq.'Y'.or.yesno.eq.'y') then
         logics=.true.
      else
         logics=.false.
      endif
      
      end
c**************************************************************************

