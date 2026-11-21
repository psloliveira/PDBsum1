      PROGRAM bulge
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: bulge.f                                                  *
C*  FUNCTION: Program to extract bulges from .sst files and classify   * 
C*            them                                                     *
C* COPYRIGHT: (1993-1996) University College London                    *
C*---------------------------------------------------------------------*
C*    AUTHOR: EDITH CHAN and GAIL HUTCHINSON                           *
C*      DATE:                                                          *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
C* 01/02/93  GH     removed all 'iseqno' from criteria for identifying *
C*                  bulges-think they were superfluous - generates same*
C*                  results without                                    *
C* 23/02/93  GH     modified and checked against list of all bulges    *
C* 14/02/95  GH     modified to run on more than one protein           *
C*  9/05/95  GH     to either write data all to one file or to write   *
c*                   out separate files for each protein               *
c* 27/03/96  GH     reorganize multiple/single file handling           *
c* 24/04/96  GH     modify to avoid some duplication of classic bulges *
c* 29/04/96  GH     this version modified to allow command line options*
c* 11/01/99  GH     slight modification around line 200 (to avoid problems
c                   when hbp() is zero
c  1/11/00   GH     correct format statement in readst to remove 'x' format
c                   change ierror from 'integer' to 'logical' in main
c 27/11/00   GH     changed to read in new sstruc format
c=========================================================================
C
C f77 -o p_bulge3 p_bulge3.f
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -o p_bulge3 p_bulge3.f
C
c***********************************************************************

      include 'p_bulge.par'
      integer finish,ok,endnfi,opnera,opnerb,sstfi,namfi
      parameter (finish=2,opnera=1,opnerb=3,ok=0,endnfi=6,
     ~     sstfi=1,namfi=3)

      integer fistat,hbp(mres,4),ibp(mres,2),iseqno(mres),
     ~     numres,iocode,i,istart,iend,nps,NAMLEN,option
CPDBSUM RAL 16/3/05 -->
      integer   lenstr
CPDBSUM RAL 16/3/05 <--
      real*8 hbe(mres,4), phi(mres,2)
      character aa(mres),sht(mres)*5,sstr(mres),brcode*4,
     ~     jseqno(mres)*6,coption
      character*(fnamln) sstnam,pdbnam,pdbcod,FINAM,BLGNAM,BLNAM2
CPDBSUM RAL 16/3/05 -->
      character*(fnamln) outnam
CPDBSUM RAL 16/3/05 <--
      LOGICAL ENSEMB,FINUSE,IERROR

      data fistat/0/,FINUSE/.FALSE./,OPTION/1/

      READ(*,'(A)') PDBNAM
      IF(PDBNAM(1:1).EQ.'%') THEN
         READ(*,'(A1)') COPTION
         ENSEMB=.TRUE.
         IF(COPTION.EQ.'L'.OR.COPTION.EQ.'l') THEN
            BLGNAM='BULGES'
            OPEN(UNIT=7,FILE=BLGNAM,STATUS='UNKNOWN',IOSTAT=IOCODE)
            IF(IOCODE.NE.0) FISTAT=OPNERB
         ENDIF
      ELSE
         ENSEMB=.FALSE.
         coption='p'
      ENDIF

c     if ensemble remove % from file name
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
         CALL GETNAM(pdbnam,istart,iend,ierror)
         PDBCOD=PDBNAM(ISTART:IEND)
         NPS=INDEX(PDBCOD,' ')
         SSTNAM=PDBCOD(1:NPS-1)//'.sst'
         OPEN(unit=sstfi,file=sstnam,status='old',iostat=iocode)
         if(iocode.ne.0) fistat=opnerb
      ENDIF
      
      IF(FISTAT.EQ.OK) THEN
c         IF(OPTION.NE.2) THEN
         IF(COPTION.NE.'L'.AND.COPTION.NE.'l') then
            BLGNAM=PDBCOD(1:NPS-1)//'.blg'
CPDBSUM RAL 14/4/05 -->
C            open(unit=7,file=BLGNAM,status='unknown',IOSTAT=IOCODE)
            open(unit=7,status='scratch',IOSTAT=IOCODE)
CPDBSUM RAL 14/4/05 <--
            IF(IOCODE.NE.0) FISTAT=OPNERB
CPDBSUM RAL 16/3/05 -->
            outnam = "bulges.promotif"
            open(unit=outblg, file=outnam, status='unknown',
     -          IOSTAT=IOCODE)
            if(iocode.ne.0) then
                fistat=opnera
            else
                call blghed(outblg)
            endif
CPDBSUM RAL 16/3/05 <--
         ENDIF
         BLNAM2=PDBCOD(1:NPS-1)//'.blgplt'
CPDBSUM RAL 14/4/05 -->
C         open(unit=8,file=BLNAM2,status='unknown',IOSTAT=IOCODE)
         open(unit=8,status='scratch',IOSTAT=IOCODE)
CPDBSUM RAL 14/4/05 <--
c     READ IN SST FILE
         call readst(numres,iseqno,jseqno,aa,sstr,sht,ibp,hbp,hbe,
     e        brcode,phi)
C     EXTRACT BULGES
         call extract(ibp,numres,sht,hbp,sstr,jseqno,aa,phi,coption,
     ~        pdbnam)
c         IF(OPTION.NE.2) close (unit=7)
         IF(COPTION.NE.'L'.AND.COPTION.NE.'l') close(unit=7)
      ELSE IF(FISTAT.EQ.OPNERB) THEN
         write(*,*) 'Unable to open the 2D structure file ',
CPDBSUM RAL 16/3/05 -->
C     +        sstnam(1:namlen)
     +        sstnam(1:lenstr(sstnam))
CPDBSUM RAL 16/3/05 <--
      ENDIF
      
C     if another file then run again
      IF(ENSEMB) THEN
         IF(FISTAT.NE.FINISH) THEN
            FISTAT=OK
            IF(FINUSE) THEN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF
         
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      subroutine readst(numres,iseqno,jseqno,aa,sstr,sht,ibp,hbp,hbe,
     e     brcode,phi)
C     READ SECONDARY STRUCTURE FILE
      include 'p_bulge.par'
            
      integer hbp(mres,4),ibp(mres,2),iseqno(mres),numres,i,j,m
      real*8 hbe(mres,4), phi(mres,2)
      character aa(mres),sht(mres)*5,sstr(mres),brcode*4,jseqno(mres)*6
      
      numres=0
      do 200 i=1,mres
         iseqno(i)=0
         jseqno(i)='      '
         aa(i)=' '
         sstr(i)=' '
         sht(i)='     '
         do 110 j=1,2
              ibp(i,j)=0
              phi(i,j)=0.0
 110       continue
           do 120 j=1,4
              hbp(i,j)=0
              hbe(i,j)=0.0
 120       continue
 200    continue

        read(1,*)
        read(1,*)
        read(1,5) brcode
        read(1,10) numres
        read(1,*)
        read(1,*)
        read(1,*)
        do 20 m=1,numres
           read(1,15) iseqno(m),jseqno(m),aa(m),sstr(m),sht(m),
     e          (ibp(m,i),i=1,2),(phi(m,i),i=1,2),
     e           (hbp(m,i), hbe(m,i),i=1,4)
 20     continue
c
 5      format(1x,a4)
 10     format(19x,i4)
c 15     format(1x,i4,1x,a6,7x,a1,5x,a1,2x,a3,5x,2(1x,i4),
c     e        2(1x,f6.1),26x,4(1x,i4,1x,f4.1))
 15     format(i5,1x,a6,7x,a1,5x,a1,2x,a5,5x,2(i5),
     e        2(1x,f6.1),26x,4(i5,1x,f4.1))
         end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CPDBSUM RAL 17/3/05 -->
C      subroutine g1bulge(jseqno,aa,phi,hbp,sstr,M,coption,pdbcod,sht)
      subroutine g1bulge(jseqno,aa,phi,hbp,sstr,M,coption,pdbcod,sht,
     -     numres)
CPDBSUM RAL 17/3/05 <--
c     ROUTINE FOR EXTRACTING G1 BULGES
c     
      include 'p_bulge.par'
      integer hbp(mres,4),M,i,j,k,l,kk,option
CPDBSUM RAL 17/3/05 -->
      integer numres
CPDBSUM RAL 17/3/05 <--
      real*8 phi(mres,2)
      character aa(mres),sstr(mres),g1type(3)*3,
     ~     jseqno(mres)*6,null*24,pdbcod*4,sht(mres)*5,coption
c
c      data g1type/'G1C','G1 ','G1T'/
      data g1type/'A C','A G','A G'/
c      data null/'                        '/
      data null/'          -999.9 -999.9 '/

c     4.2.93 DON'T RESTRICT TO JUST GLYCINE
      do 19  i=1,3,2
         do 18 j=2,4,2
            if(hbp(m,i).eq.hbp(m,j)+1.and.hbp(m,j).gt.0.and.
     e           hbp(m,i).gt.0) then
            if((sstr(m).eq.'E'.or.sstr(m).eq.'e').and.
     e           (sstr(hbp(m,i)).eq.'E'.or.
     e           sstr(hbp(m,i)).eq.'e')) then
               k=hbp(m,i)
               l=hbp(m,j)
               if (hbp(k,1).ne.hbp(k,2).and.hbp(k,1).gt.0.0.and.
     e              hbp(k,2).gt.0.0) then
c                  if(option.eq.2) then
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,161)pdbcod,g1type(3),jseqno(m),aa(m)
     e                    ,(phi(m,kk),kk=1,2),jseqno(l),aa(l),
     e                    (phi(l,kk),kk=1,2),jseqno(k),aa(k),
     e                    (phi(k,kk),kk=1,2),null,null
                  else
                     write(7,16)g1type(3),jseqno(m),aa(m)
     e                    ,(phi(m,kk),kk=1,2),jseqno(l),aa(l),
     e                    (phi(l,kk),kk=1,2),jseqno(k),aa(k),
     e                    (phi(k,kk),kk=1,2),null,null,m,l,k
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),g1type(3),m,
     -                    l,k,0,0,jseqno,aa,phi,sstr,hbp,mres,
     -                    m,l-2,l+3,numres)
CPDBSUM RAL 16/3/05 <--
                  endif
                  write(8,116) m,l,k
                  call wr_plot(m,l-2,l+3,jseqno,sstr,hbp)
               elseif (phi(l,1).lt.0.0) then
c     24/4/96 added in this extra condition to avoid duplication of 
c     classic bulges
                  if(index(sht(l),'*').eq.0.and.
     ~                 index(sht(k),'*').eq.0) then
c                  if(option.eq.2) then
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,161)pdbcod,g1type(1),jseqno(m),aa(m),
     e                    (phi(m,kk),kk=1,2),jseqno(l),aa(l),
     e                    (phi(l,kk),kk=1,2),jseqno(k),aa(k),
     e                    (phi(k,kk),kk=1,2),null,null
                  else
                     write(7,16)g1type(1),jseqno(m),aa(m),
     e                    (phi(m,kk),kk=1,2),jseqno(l),aa(l),
     e                    (phi(l,kk),kk=1,2),jseqno(k),aa(k),
     e                    (phi(k,kk),kk=1,2),null,null,m,l,k
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),g1type(1),m,
     -                    l,k,0,0,jseqno,aa,phi,sstr,hbp,mres,
     -                    m,l-2,l+3,numres)
CPDBSUM RAL 16/3/05 <--
                  endif
                     write(8,116) m,l,k
                  call wr_plot(m,l-2,l+3,jseqno,sstr,hbp)
                  endif
               else
c                  if(option.eq.2) then
                  if(coption.eq.'L'.or.coption.eq.'l') then 
                    write(7,161)pdbcod,g1type(2),jseqno(m),aa(m),
     e                    (phi(m,kk),kk=1,2),jseqno(l),aa(l),
     e                    (phi(l,kk),kk=1,2),jseqno(k),aa(k),
     e                    (phi(k,kk),kk=1,2),null,null
                  else
                     write(7,16)g1type(2),jseqno(m),aa(m),
     e                    (phi(m,kk),kk=1,2),jseqno(l),aa(l),
     e                    (phi(l,kk),kk=1,2),jseqno(k),aa(k),
     e                    (phi(k,kk),kk=1,2),null,null,m,l,k
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),g1type(2),m,
     -                    l,k,0,0,jseqno,aa,phi,sstr,hbp,mres,
     -                    m,l-2,l+3,numres)
CPDBSUM RAL 16/3/05 <--
                  endif
                  write(8,116) m,l,k
                  call wr_plot(m,l-2,l+3,jseqno,sstr,hbp)
               endif
            endif
            endif
 18      continue
 19   continue

 16   format(a3,1x,3(a6,1x,a1,1x,2(1x,f6.1),1x),a24,a24,3(1x,i4))
 161  format(a4,1x,a3,1x,3(a6,1x,a1,1x,2(1x,f6.1),1x),a24,a24)
 116  format(4x,3(i6,18x))

      return
      end
c**************************************************************************
      subroutine extract(ibp,numres,sht,hbp,sstr,jseqno,aa,phi,cOPTION,
     ~     PDBNAM)

      include 'p_bulge.par'

      integer hbp(mres,4),ibp(mres,2),d(2),numres,m,k,j,i,option,l
      real*8 phi(mres,2)
      character aa(mres),sht(mres)*5,sstr(mres),jseqno(mres)*6,
     ~     type(17)*3,null*24,PDBCOD*4,coption
      character*(fnamln) pdbnam

c      data type/'AGX','PE2','PE1','AE2','AE1','PWI','AM1','AM2',
c     ~     'AM3','AWI','POU','ASW','ACL','AS3','AS4','PS3','PS4'/
      data type/'AGX','P B','P B','A B','A B','P W','A C','A S',
     ~     'AM3','A W','P C','A W','A C','A S','A S','P S','P S'/
c      data null/'                        '/
      data null/'          -999.9 -999.9 '/

      I=INDEX(PDBNAM,'.')
      IF(I.GT.4) THEN
         PDBCOD=PDBNAM(I-4:I-1)
c         PDBCOD=PDBNAM(I-5:I-2)
      ELSE
         if(i.ne.0) then
            PDBCOD=PDBNAM(1:4)
         else
            l=index(pdbnam,' ')
            pdbcod=pdbnam(l-4:l-1)
         endif
      ENDIF

      do 20 m=2,numres-1
         d(1)=ibp(m+1,1)-ibp(m-1,1)
         d(2)=ibp(m+1,2)-ibp(m-1,2)
         if(m.ge.10.and.m.le.50) then
         endif
c     
         IF(INDEX(SHT(M),'*').NE.0) THEN
c         do 109 k = 1,2
c            if (sht(m)(k+1:k+2).eq.'*') then
         do 109 k = 1,2
CPDBSUM RAL 16/3/05 -->
C            if (sht(m)(k+3:k+4).eq.'*') then
            if (sht(m)(k+3:k+3).eq.'*') then
CPDBSUM RAL 16/3/05 <--
c     extra (PE2 or P_B)
               if(d(k).eq.2) then  
                  if (hbp(m,3).eq.ibp(m-1,k).and.
     e                 hbp(m,1).eq.ibp(m+1,k)-1) then
c                     IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,171) pdbcod,type(2),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)-1),
     e                       aa(ibp(m+1,k)-1),(phi(ibp(m+1,k)-1,j)
     e                       ,j=1,2),null,null
                     ELSE
                        write(7,17) type(2),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)-1),
     e                       aa(ibp(m+1,k)-1),(phi(ibp(m+1,k)-1,j)
     e                       ,j=1,2),null,null,m,ibp(m+1,k)-1
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(2),0,m,
     -                       ibp(m+1,k)-1,0,0,jseqno,aa,phi,sstr,hbp,
     -                       mres,ibp(m+1,k)-1,m-2,m+3,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,117) m,ibp(m+1,k)-1
                     call wr_plot(ibp(m+1,k)-1,m-2,m+3,jseqno,sstr,
     ~                    hbp)
C     PE1 or P_B
                  elseif (hbp(m,1).eq.ibp(m-1,k)) then
c                     IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,171) pdbcod,type(3),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2)
     e                       ,jseqno(ibp(m+1,k)-1)
     e                       ,aa(ibp(m+1,k)-1),(phi(ibp(m+1,k)-1,j)
     e                       ,j=1,2),null,null
                     ELSE
                        write(7,17) type(3),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2)
     e                       ,jseqno(ibp(m+1,k)-1)
     e                       ,aa(ibp(m+1,k)-1),(phi(ibp(m+1,k)-1,j)
     e                       ,j=1,2),null,null,m,ibp(m+1,k)-1
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(3),0,m,
     -                       ibp(m+1,k)-1,0,0,jseqno,aa,phi,sstr,hbp,
     -                       mres,ibp(m+1,k)+1,m-2,m+3,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,117) m,ibp(m+1,k)-1
                     call wr_plot(ibp(m+1,k)+1,m-2,m+3,jseqno,sstr,
     ~                    hbp)
                  endif
               elseif (d(k).eq.-2) then
                  if (hbp(m,1).eq.ibp(m+1,k)+1.and.
     e                 hbp(m,3).eq.ibp(m+1,k)) then
c     AE2 OR A_B
c                     if(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,171) pdbcod,type(4),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),
     e                       aa(ibp(m+1,k)+1),(phi(ibp(m+1,k)+1,j)
     e                       ,j=1,2),null,null
                     ELSE
                        write(7,17) type(4),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),
     e                       aa(ibp(m+1,k)+1),(phi(ibp(m+1,k)+1,j)
     e                       ,j=1,2),null,null,m,ibp(m+1,k)+1
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(4),0,m,
     -                       ibp(m+1,k)+1,0,0,jseqno,aa,phi,sstr,hbp,
     -                       mres,ibp(m+1,k)+1,m-2,m+3,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,117) m,ibp(m+1,k)+1
                     call wr_plot(ibp(m+1,k)+1,m-2,m+3,jseqno,sstr,
     ~                    hbp)
                  elseif (hbp(m,1).eq.ibp(m+1,k)+1) then
C     AE1 or A_B
c                     IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,171) pdbcod,type(5),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),
     e                       aa(ibp(m+1,k)+1),(phi(ibp(m+1,k)+1,j),
     e                       j=1,2),null,null
                     ELSE
                        write(7,17) type(5),null,
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),
     e                       aa(ibp(m+1,k)+1),(phi(ibp(m+1,k)+1,j),
     e                       j=1,2),null,null,m,ibp(m+1,k)+1
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(5),0,m,
     -                       ibp(m+1,k)+1,0,0,jseqno,aa,phi,sstr,hbp,
     -                       mres,ibp(m+1,k)+1,m-2,m+3,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,17) m,ibp(m+1,k)+1
                     call wr_plot(ibp(m+1,k)+1,m-2,m+3,jseqno,sstr,
     ~                    hbp)
                  endif
c     P-wide (PWI or P_W)
               elseif (d(k).eq.3) then
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,161) PDBCOD,type(6),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+1),
     e                    aa(ibp(m-1,k)+1),(phi(ibp(m-1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m+1,k)-1),
     e                    aa(ibp(m+1,k)-1),(phi(ibp(m+1,k)-1,j),j=1,2),
     e                    null,null
                  ELSE
                     WRITE(7,16) type(6),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+1),
     e                    aa(ibp(m-1,k)+1),(phi(ibp(m-1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m+1,k)-1),
     e                    aa(ibp(m+1,k)-1),(phi(ibp(m+1,k)-1,j),j=1,2),
     e                    null,null,m,ibp(m-1,k)+1,ibp(m+1,k)-1
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),type(6),m,
     -                    ibp(m-1,k)+1,ibp(m+1,k)-1,0,0,jseqno,aa,phi,
     -                    sstr,hbp,mres,m,ibp(m-1,k)-1,ibp(m-1,k)+4,
     -                    numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                   write(8,116) m,ibp(m-1,k)+1,ibp(m+1,k)-1
                   call wr_plot(m,ibp(m-1,k)-1,ibp(m-1,k)+4,jseqno,
     ~                  sstr,hbp)
c     A-mix
               elseif (d(k).eq.-3) then
c     AM1 or A_C
                  if (hbp(m,1).eq.ibp(m+1,k)+2.and.hbp(m,2).eq.
     e                 ibp(m+1,k)+1) then
c                     IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,161) PDBCOD,type(7),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                       (phi(ibp(m+1,k)+1,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                       (phi(ibp(m+1,k)+2,j),j=1,2),null,null
                     ELSE
                        WRITE(7,16) type(7),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                       (phi(ibp(m+1,k)+1,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                       (phi(ibp(m+1,k)+2,j),j=1,2),null,null,
     e                       m,ibp(m+1,k)+1,ibp(m+1,k)+2
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(7),m,
     -                       ibp(m+1,k)+1,ibp(m+1,k)+2,0,0,jseqno,aa,
     -                       phi,sstr,hbp,mres,m,ibp(m+1,k)-1,
     -                       ibp(m+1,k)+4,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,116) m,ibp(m+1,k)+1,ibp(m+1,k)+2
                   call wr_plot(m,ibp(m+1,k)-1,ibp(m+1,k)+4,jseqno,
     ~                    sstr,hbp)
c     AM2 or A_S
                  elseif (hbp(m,1).eq.ibp(m+1,k)+2.and.hbp(m,2).eq.0) 
     e                    then
c                     IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,161) PDBCOD,type(8),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                       (phi(ibp(m+1,k)+1,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                       (phi(ibp(m+1,k)+2,j),j=1,2),null,null
                     ELSE
                        WRITE(7,16) type(8),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                       (phi(ibp(m+1,k)+1,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                       (phi(ibp(m+1,k)+2,j),j=1,2),null,null,
     e                       m,ibp(m+1,k)+1,ibp(m+1,k)+2
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(8),m,
     -                       ibp(m+1,k)+1,ibp(m+1,k)+2,0,0,jseqno,aa,
     -                       phi,sstr,hbp,mres,m,ibp(m+1,k)-1,
     -                       ibp(m+1,k)+4,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,116) m,ibp(m+1,k)+1,ibp(m+1,k)+2
                   call wr_plot(m,ibp(m+1,k)-1,ibp(m+1,k)+4,jseqno,
     ~                    sstr,hbp)
                  elseif (hbp(m,1).eq.ibp(m+1,k)+1.and.hbp(m,2).eq.0) 
     e                    then
c     SEEMS TO BE SUPERFLUOUS AT THE MOMENT - NONE OF THESE ARE IDENFITIED
c     AM3 
c                     write(7,16) brcode,type(9),
c     e                    brcode,jseqno(m),aa(m),(phi(m,j),j=1,2),
c     e                    brcode,jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
c     e                    (phi(ibp(m+1,k)+1,j),j=1,2),
c     e                    brcode,jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
c     e                    (phi(ibp(m+1,k)+2,j),j=1,2),null,null
c     A-wide
                  elseif (hbp(m,1).ne.ibp(m+1,k)+2.or.
     e                    hbp(m,2).ne.ibp(m+1,k)+2) then
c     AWI or A_W
c                     IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,161) PDBCOD,type(10),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                       (phi(ibp(m+1,k)+1,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                       (phi(ibp(m+1,k)+2,j),j=1,2),null,null
                     ELSE
                        WRITE(7,16) type(10),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                       (phi(ibp(m+1,k)+1,j),j=1,2),
     e                       jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                       (phi(ibp(m+1,k)+2,j),j=1,2),null,null,
     e                       m,ibp(m+1,k)+1,ibp(m+1,k)+2
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(m)(1:1),type(10),m,
     -                       ibp(m+1,k)+1,ibp(m+1,k)+2,0,0,jseqno,aa,
     -                       phi,sstr,hbp,mres,m,ibp(m+1,k)-1,
     -                       ibp(m+1,k)+4,numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                        write(8,116) m,ibp(m+1,k)+1,ibp(m+1,k)+2
                   call wr_plot(m,ibp(m+1,k)-1,ibp(m+1,k)+4,jseqno,
     ~                    sstr,hbp)
                  endif
c     classic POU or P_C
               elseif (d(k).eq.1) then
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,161) PDBCOD,type(11),
     e                    jseqno(ibp(m-1,k)),aa(ibp(m-1,k)),
     e                    (phi(ibp(m-1,k),j),j=1,2),
     e                    jseqno(m-1),aa(m-1),(phi(m-1,j),j=1,2),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),null,
     e                    null
                  ELSE
                     WRITE(7,16) type(11),
     e                    jseqno(ibp(m-1,k)),aa(ibp(m-1,k)),
     e                    (phi(ibp(m-1,k),j),j=1,2),
     e                    jseqno(m-1),aa(m-1),(phi(m-1,j),j=1,2),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),null,
     e                    null,ibp(m-1,k),m-1,m
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(ibp(m-1,k))(1:1),
     -                    type(11),ibp(m-1,k),m-1,m,0,0,jseqno,aa,phi,
     -                    sstr,hbp,mres,ibp(m-1,k),m-2,m+3,numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                  write(8,116) ibp(m-1,k),m-1,m
                  call wr_plot(ibp(m-1,k),m-2,m+3,jseqno,sstr,hbp)
               elseif (d(k).eq.-1) then
c     ASW or A_W
                  if (hbp(m,2).eq.ibp(m-1,k)) then
c                     IF(OPTION.EQ.2) THEN
                     if(coption.eq.'L'.or.coption.eq.'l') then
                        write(7,161) PDBCOD,type(12),
     e                       jseqno(ibp(m-1,k)-1),aa(ibp(m-1,k)-1),
     e                       (phi(ibp(m-1,k)-1,j),j=1,2),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(m+1),aa(m+1),(phi(m+1,j),j=1,2),
     e                       null,null
                     ELSE
                        WRITE(7,16) type(12),
     e                       jseqno(ibp(m-1,k)-1),aa(ibp(m-1,k)-1),
     e                       (phi(ibp(m-1,k)-1,j),j=1,2),
     e                       jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                       jseqno(m+1),aa(m+1),(phi(m+1,j),j=1,2),
     e                       null,null,ibp(m-1,k)-1,m,m+1
CPDBSUM RAL 16/3/05 -->
C----                   Write out data to .promotif data file
                        CALL WRIBLG(outblg,jseqno(ibp(m-1,k)-1)(1:1),
     -                       type(12),ibp(m-1,k)-1,m,m+1,0,0,jseqno,aa,
     -                       phi,sstr,hbp,mres,ibp(m-1,k),m-2,m+3,
     -                       numres)
CPDBSUM RAL 16/3/05 <--
                     ENDIF
                     write(8,116) ibp(m-1,k)-1,m,m+1
                     call wr_plot(ibp(m-1,k),m-2,m+3,jseqno,sstr,hbp)
                  else
c     ACL or A_C
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,161) PDBCOD,type(13),
     e                    jseqno(ibp(m+1,k)),aa(ibp(m+1,k)),
     e                    (phi(ibp(m+1,k),j),j=1,2),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(m+1),aa(m+1),
     e                    (phi(m+1,j),j=1,2),null,null
                  ELSE
                     WRITE(7,16) type(13),
     e                    jseqno(ibp(m+1,k)),aa(ibp(m+1,k)),
     e                    (phi(ibp(m+1,k),j),j=1,2),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(m+1),aa(m+1),
     e                    (phi(m+1,j),j=1,2),null,null,
     e                    ibp(m+1,k),m,m+1
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(ibp(m+1,k))(1:1),
     -                    type(13),ibp(m+1,k),m,m+1,0,0,jseqno,aa,
     -                    phi,sstr,hbp,mres,ibp(m+1,k),m-2,m+3,numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                     write(8,116) ibp(m+1,k),m,m+1
                  call wr_plot(ibp(m+1,k),m-2,m+3,jseqno,sstr,hbp)
                  endif
cspecial
c     This is for special and anti
               elseif (d(k).eq.-4) then
c     AS3 or A_S
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,241) PDBCOD,type(14),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                    (phi(ibp(m+1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                    (phi(ibp(m+1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+3),aa(ibp(m+1,k)+3),
     e                    (phi(ibp(m+1,k)+3,j),j=1,2),null
                  ELSE
                     WRITE(7,24) type(14),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                    (phi(ibp(m+1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                    (phi(ibp(m+1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+3),aa(ibp(m+1,k)+3),
     e                    (phi(ibp(m+1,k)+3,j),j=1,2),null,
     e                    m,ibp(m+1,k)+1,ibp(m+1,k)+2,
     e                    ibp(m+1,k)+3
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),
     -                    type(14),m,ibp(m+1,k)+1,ibp(m+1,k)+2,
     -                    ibp(m+1,k)+3,0,jseqno,aa,phi,sstr,hbp,mres,
     -                    m,ibp(m+1,k)-1,ibp(m+1,k)+5,numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                     write(8,124) m,ibp(m+1,k)+1,ibp(m+1,k)+2,
     ~                 ibp(m+1,k)+3
                  call wr_plot(m,ibp(m+1,k)-1,ibp(m+1,k)+5,jseqno,
     ~                 sstr,hbp)
c     AS4 or A_S
               elseif (d(k).eq.-5) then
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,251) PDBCOD,type(15),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                    (phi(ibp(m+1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                    (phi(ibp(m+1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+3),aa(ibp(m+1,k)+3),
     e                    (phi(ibp(m+1,k)+3,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+4),aa(ibp(m+1,k)+4),
     e                    (phi(ibp(m+1,k)+4,j),j=1,2)
                  ELSE
                     WRITE(7,25) type(15),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+1),aa(ibp(m+1,k)+1),
     e                    (phi(ibp(m+1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+2),aa(ibp(m+1,k)+2),
     e                    (phi(ibp(m+1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+3),aa(ibp(m+1,k)+3),
     e                    (phi(ibp(m+1,k)+3,j),j=1,2),
     e                    jseqno(ibp(m+1,k)+4),aa(ibp(m+1,k)+4),
     e                    (phi(ibp(m+1,k)+4,j),j=1,2),
     e                    m,ibp(m+1,k)+1,ibp(m+1,k)+2,
     e                    ibp(m+1,k)+3,ibp(m+1,k)+4
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),type(15),m,
     -                    ibp(m+1,k)+1,ibp(m+1,k)+2,ibp(m+1,k)+3,
     -                    ibp(m+1,k)+4,jseqno,aa,phi,sstr,hbp,mres,
     -                    m,ibp(m+1,k)-1,ibp(m+1,k)+6,numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                     write(8,125) m,ibp(m+1,k)+1,ibp(m+1,k)+2,
     ~                 ibp(m+1,k)+3,ibp(m+1,k)+4
                  call wr_plot(m,ibp(m+1,k)-1,ibp(m+1,k)+6,jseqno,
     ~                 sstr,hbp)
               elseif (d(k).eq.4) then
c     PS3 or P_S
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,241) PDBCOD,type(16),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+1),aa(ibp(m-1,k)+1),
     e                    (phi(ibp(m-1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+2),aa(ibp(m-1,k)+2),
     e                    (phi(ibp(m-1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+3),aa(ibp(m-1,k)+3),
     e                    (phi(ibp(m-1,k)+3,j),j=1,2),null
                  ELSE
                     WRITE(7,24) type(16),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+1),aa(ibp(m-1,k)+1),
     e                    (phi(ibp(m-1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+2),aa(ibp(m-1,k)+2),
     e                    (phi(ibp(m-1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+3),aa(ibp(m-1,k)+3),
     e                    (phi(ibp(m-1,k)+3,j),j=1,2),null,
     e                    m,ibp(m-1,k)+1,ibp(m-1,k)+2,
     e                    ibp(m-1,k)+3
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),type(16),m,
     -                    ibp(m-1,k)+1,ibp(m-1,k)+2,ibp(m-1,k)+3,0,
     -                    jseqno,aa,phi,sstr,hbp,mres,m,ibp(m-1,k)-1,
     -                    ibp(m-1,k)+5,numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                     write(8,124) m,ibp(m-1,k)+1,ibp(m-1,k)+2,
     ~                 ibp(m-1,k)+3
                  call wr_plot(m,ibp(m-1,k)-1,ibp(m-1,k)+5,jseqno,
     ~                 sstr,hbp)
               elseif (d(k).eq.5) then
c     PS4 or P_S
c                  IF(OPTION.EQ.2) THEN
                  if(coption.eq.'L'.or.coption.eq.'l') then
                     write(7,251) PDBCOD,type(17),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+1),aa(ibp(m-1,k)+1),
     e                    (phi(ibp(m-1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+2),aa(ibp(m-1,k)+2),
     e                    (phi(ibp(m-1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+3),aa(ibp(m-1,k)+3),
     e                    (phi(ibp(m-1,k)+3,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+4),aa(ibp(m-1,k)+4),
     e                    (phi(ibp(m-1,k)+4,j),j=1,2)
                  ELSE
                     WRITE(7,25) type(17),
     e                    jseqno(m),aa(m),(phi(m,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+1),aa(ibp(m-1,k)+1),
     e                    (phi(ibp(m-1,k)+1,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+2),aa(ibp(m-1,k)+2),
     e                    (phi(ibp(m-1,k)+2,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+3),aa(ibp(m-1,k)+3),
     e                    (phi(ibp(m-1,k)+3,j),j=1,2),
     e                    jseqno(ibp(m-1,k)+4),aa(ibp(m-1,k)+4),
     e                    (phi(ibp(m-1,k)+4,j),j=1,2),
     e                    m,ibp(m-1,k)+1,ibp(m-1,k)+2,
     e                    ibp(m-1,k)+3,ibp(m-1,k)+4
CPDBSUM RAL 16/3/05 -->
C----                Write out data to .promotif data file
                     CALL WRIBLG(outblg,jseqno(m)(1:1),type(17),m,
     -                    ibp(m-1,k)+1,ibp(m-1,k)+2,ibp(m-1,k)+3,
     -                    ibp(m-1,k)+4,jseqno,aa,phi,sstr,hbp,mres,
     -                    m,ibp(m-1,k)-1,ibp(m-1,k)+6,numres)
CPDBSUM RAL 16/3/05 <--
                  ENDIF
                  write(8,125) m,ibp(m-1,k)+1,ibp(m-1,k)+2,
     ~                 ibp(m-1,k)+3,ibp(m-1,k)+4
                  call wr_plot(m,ibp(m-1,k)-1,ibp(m-1,k)+6,jseqno,
     ~                 sstr,hbp)
               endif
            endif
 109     continue
      else
c     G1 BULGES!!
         IF(HBP(M,1).NE.0.AND.HBP(M,2).NE.0) THEN
            IF(INDEX(SHT(HBP(M,1)),'*').EQ.0.AND.
     ~         INDEX(SHT(HBP(M,2)),'*').EQ.0)
     ~          CALL G1BULGE(JSEQNO,AA,PHI,HBP,SSTR,M,cOPTION,PDBCOD,
CPDBSUM RAL 17/3/05 -->
C     ~           SHT)
     ~           SHT,NUMRES)
CPDBSUM RAL 17/3/05 <--
         ENDIF
      ENDIF
c     
 20   continue
 16   format(a3,1x,3(a6,1x,a1,1x,2(1x,f6.1),1x),a24,a24,3(1x,i4))
 161   format(A4,1X,a3,1x,3(a6,1x,a1,1x,2(1x,f6.1),1x),a24,a24)
 116  format(4x,3(i6,18x))
 17   format(a3,1x,a24,2(a6,1x,a1,1x,2(1x,f6.1),1x),
     ~     a24,a24,3(1x,i4))
 171  format(A4,1X,a3,1x,a24,2(a6,1x,a1,1x,2(1x,f6.1),1x),
     ~     a24,a24)
 117  format(28x,2(i6,18x))
 24   format(a3,1x,4(a6,1x,a1,1x,2(1x,f6.1),1x),a24,4(1x,i4))
 241  format(A4,1X,a3,1x,4(a6,1x,a1,1x,2(1x,f6.1),1x),a24)
 124  format(4x,4(i6,18x))
 25   format(a3,1x,5(a6,1x,a1,1x,2(1x,f6.1),1x),5(1x,i4))
 251  format(A4,1X,a3,1x,5(a6,1x,a1,1x,2(1x,f6.1),1x))
 125  format(4x,5(i6,18x))

      return
      end
c*******************************************************************************
      CHARACTER*6 FUNCTION Seq(ires,jseqno)
c     converts to author sequence number assigning null if residue number is
c     zero

      include 'p_bulge.par'
     
      integer ires
      character jseqno(mres)*6,NONE*6

      data NONE/'000000'/

      if(ires.ne.0) then
         Seq=jseqno(ires)
      else
         Seq=NONE
      endif
      
      return
      end
c***********************************************************************
      CHARACTER FUNCTION Sst(ires,sstr)
c     converts to author sequence number assigning NONE if residue number is
c     zero

      include 'p_bulge.par'
     
      integer ires
      character sstr(mres),NONE

      data NONE/'0'/

      if(ires.ne.0) then
         Sst=sstr(ires)
      else
         Sst=NONE
      endif
      
      return
      end
c**********************************************************************
      subroutine wr_plot(m,i1,i2,jseqno,sstr,hbp)

      include 'p_bulge.par'

      integer hbp(mres,4),i1,i2,j,m
      character jseqno(mres)*6,sstr(mres),seq1*6,seq2*6,
     ~     seq3*6,seq4*6,Seq*6

      do 20 j=i1,i2
         if(j.gt.0) then
            seq1=Seq(hbp(j,1),jseqno)
            seq2=Seq(hbp(j,2),jseqno)
            seq3=Seq(hbp(j,3),jseqno)
            seq4=Seq(hbp(j,4),jseqno)
            write(8,10) j,jseqno(j),sstr(j),
     ~           seq1,seq2,seq3,seq4
         else
            write(8,50)
         endif
 20   continue
      write(8,*)
      do 30 j=m-2,m+2
         if(j.gt.0) then
            write(8,40) j,jseqno(j),sstr(j)
         else
            write(8,50)
         endif
 30   continue
      write(8,*)
 10   format(i6,1x,a6,1x,a1,1x,4(a6,1x))
 40   format(i6,1x,a6,1x,a1)
 50   format('000000 000000')
      return
      end
c************************************************************

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

c**********************************************************************
CPDBSUM RAL 16/3/05 -->
C**************************************************************************
C
C  SUBROUTINE BLGHED  -  Write header records to bulges.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE BLGHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 40)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN',
     -     'CHAIN1', 'RESNUM1', 'RESNAM1',
     -     'CHAIN2', 'RESNUM2', 'RESNAM2',
     -     'CHAIN3', 'RESNUM3', 'RESNAM3',
     -     'CHAIN4', 'RESNUM4', 'RESNAM4',
     -     'CHAIN5', 'RESNUM5', 'RESNAM5',
     -     'TYPE', 'CLASS',
     -     'HAVE1', 'HAVE2', 'HAVE3', 'HAVE4', 'HAVE5', 
     -     'PHI1', 'PSI1', 'PHI2', 'PSI2', 'PHI3', 'PSI3',
     -     'PHI4', 'PSI4', 'PHI5', 'PSI5',
     -     'CHAIN_STRAND1', 'START_STRAND1', 'END_STRAND1',
     -     'CHAIN_STRAND2', 'START_STRAND2', 'END_STRAND2',
     -     'PLOT_DATA' /

C---- Define field id tag
      IDTAG = 'BL_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:BETA_BULGES'
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
C  SUBROUTINE WRIBLG  -  Write current beta bulge to bulges.promotif file
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRIBLG(IUNIT,CHAIN,TYPE,IPOS1,IPOS2,IPOS3,IPOS4,IPOS5,
     -     JSEQNO,AA,PHI,SSTR,HBP,MRES,MPOSX,MPOS1,MPOS2,NRESID)

      INTEGER       MRES

      CHARACTER*1   AA(MRES), CHAIN, CHAINS(5), SSTR(MRES), STRCHN(2)
      CHARACTER*3   RSNAM(5), TYPE
      CHARACTER*5   RSNUM(5), TRUFAL(5), STRNUM(2,2)
      CHARACTER*6   JSEQNO(MRES), HBOND(4)
      CHARACTER*20  TMPSTR
      CHARACTER*900 HBONDS, PLTDAT

      INTEGER       HBP(MRES,4), I, IHBOND, IPOINT, IPOS(5), IPOS1,
     -              IPOS2, IPOS3, IPOS4, IPOS5, IPT, IUNIT, J, JPOS(5),
     -              LENSTR, MPOSX, MPOS1, MPOS2, NHBOND, NPOS, NRESID,
     -              NRES1, NRES2

      REAL*8        PHI(MRES,2), PHIVAL(5,2)


C---- Save the residue positions
      IPOS(1) = IPOS1
      IPOS(2) = IPOS2
      IPOS(3) = IPOS3
      IPOS(4) = IPOS4
      IPOS(5) = IPOS5
      NPOS = 0
      DO 20, I = 1, 5
          JPOS(I) = 0
 20   CONTINUE
      DO 40, J = 1, 2
          STRCHN(J) = ' '
          DO 30, I = 1, 2
              STRNUM(I,J) = ' '
 30       CONTINUE
 40   CONTINUE

C---- Store the data for the 5 bulge residues
      DO 100, I = 1, 5
          IF (IPOS(I).GT.0) THEN
              CHAINS(I) = JSEQNO(IPOS(I))(1:1)
              RSNUM(I) = JSEQNO(IPOS(I))(2:6)
              CALL AA2NAM(AA(IPOS(I)),RSNAM(I))
              PHIVAL(I,1) = PHI(IPOS(I),1)
              PHIVAL(I,2) = PHI(IPOS(I),2)
              TRUFAL(I) = 'TRUE'
              NPOS = NPOS + 1
              JPOS(NPOS) = IPOS(I)
          ELSE
              CHAINS(I) = ' '
              RSNAM(I) = ' '
              RSNUM(I) = ' '
              PHIVAL(I,1) = -999.9
              PHIVAL(I,2) = -999.9
              TRUFAL(I) = 'FALSE'
          ENDIF
 100  CONTINUE

C---- Calculate number of residues in the 1,2 and X strands
      NRES1 = MPOS2 - MPOS1 + 1
      NRES2 = 5

C---- Initialise current plot-details line
      IPOINT = 0
      PLTDAT = ' '

C---- Add count of residues in the 1, 2 strand
      WRITE(PLTDAT,120) NRES1
 120  FORMAT(I1)
      IPOINT = 1

C---- Store the H-bond data required for the bulge plot
      DO 300, J = MPOS1, MPOS2

C----     Initialise count of this residue's H-bond partners
          NHBOND = 0

C----     If valid, get the chain and residue number
          IF (J.GT.0 .AND. J.LE.NRESID) THEN

C----         If this is the first residue of this strand, store
              IF (STRNUM(1,1).EQ.' ') THEN
                  STRCHN(1) = JSEQNO(J)(1:1)
                  STRNUM(1,1) = JSEQNO(J)(2:6)
              ENDIF

C----         Store residue as last of this strand
              STRNUM(2,1) = JSEQNO(J)(2:6)

C----         Initialise H-bonds line
              IPT = 0
              HBONDS = ' '

C----         Determine number of H-bond entries to be written out
              DO 150, IHBOND = 1, 4

C----             If a valid H-bond, then store position
                  IF (HBP(J,IHBOND).GT.0) NHBOND = iHBOND
 150          CONTINUE

C----         Loop over this residue's H-bond partners
              DO 200, IHBOND = 1, NHBOND

C----             If a valid H-bond, then store partner
                  IF (HBP(J,IHBOND).GT.0) THEN

C----                 Add H-bond to H-bonds line
                      HBONDS(IPT + 1:) = JSEQNO(HBP(J,IHBOND))
                      IPT = IPT + 6

C----                 Add single-letter a.a. code
                      HBONDS(IPT + 1:) = AA(HBP(J,IHBOND))
                      IPT = IPT + 1

C----             Otherwise, write out blank entry
                  ELSE

C----                 Write blank entry
                      HBONDS(IPT + 1:) = '0000000'
                      IPT = IPT + 7
                  ENDIF
 200          CONTINUE

C----         Write residue number and secondary structure details to
C             current line
              WRITE(TMPSTR,220) JSEQNO(J), AA(J), SSTR(J), NHBOND
 220          FORMAT(A6,A1,A1,I1)
              PLTDAT(IPOINT + 1:) = TMPSTR
              IPOINT = IPOINT + LENSTR(TMPSTR)

C----         Write H-bond details to current line
              PLTDAT(IPOINT + 1:) = HBONDS
              IPOINT = IPOINT + LENSTR(HBONDS)

C----     If residue not valid, add empty reference
          ELSE

C----         Write details for empty residue
              WRITE(TMPSTR,240) ' ', ' ', 0
 240          FORMAT(A7,A1,I1)

C----         Add to current line
              PLTDAT(IPOINT + 1:) = TMPSTR
              IPOINT = IPOINT + LENSTR(TMPSTR)
          ENDIF
 300  CONTINUE

C---- Add count of residues in the X strand
      WRITE(TMPSTR,120) NRES2
      PLTDAT(IPOINT + 1:) = TMPSTR
      IPOINT = IPOINT + 1

C---- Store the residues in this strand for the bulge plot
      DO 500, J = MPOSX - 2, MPOSX + 2
          
C----     If valid, write the residue details
          IF (J.GT.0 .AND. J.LE.NRESID) THEN
          
C----         If this is the first residue of this strand, store
              IF (STRNUM(1,2).EQ.' ') THEN
                  STRCHN(2) = JSEQNO(J)(1:1)
                  STRNUM(1,2) = JSEQNO(J)(2:6)
              ENDIF

C----         Store residue as last of this strand
              STRNUM(2,2) = JSEQNO(J)(2:6)

C----         Write residue number and secondary structure details to
C             current line
              WRITE(TMPSTR,320) JSEQNO(J), AA(J), SSTR(J)
 320          FORMAT(A6,A1,A1)
              PLTDAT(IPOINT + 1:) = TMPSTR
              IPOINT = IPOINT + LENSTR(TMPSTR)

C----     If residue not valid, add empty reference
          ELSE

C----         Write details for empty residue
              WRITE(TMPSTR,340) ' ', ' '
 340          FORMAT(A7,A1)

C----         Add to current line
              PLTDAT(IPOINT + 1:) = TMPSTR
              IPOINT = IPOINT + LENSTR(TMPSTR)
          ENDIF
 500  CONTINUE

C---- Write the bulge data
      WRITE(IUNIT,600) CHAIN, (CHAINS(I), RSNUM(I), RSNAM(I), I = 1, 5),
     -    TYPE(1:1), TYPE(3:3), (TRUFAL(I), I = 1, 5),
     -    ((PHIVAL(I,J), J = 1, 2), I = 1, 5),
     -    (STRCHN(J), (STRNUM(I,J), I = 1, 2), J = 1, 2),
     -    PLTDAT(1:IPOINT)
 600  FORMAT(A1,'::',5(A1,'::',A5,'::',A3,'::'),2(A1,'::'),5(A5,'::'),
     -     10(F6.1,'::'),2(A1,'::',2(A5,'::')),A)

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
CPDBSUM RAL 16/3/05 <--
