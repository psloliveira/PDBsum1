      PROGRAM hairpin
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: p_hairpin_ensem.f                                        *
C*  FUNCTION: Program to identify and classify hairpins in proteins    * 
C* COPYRIGHT: (1993-1996) University College London                    *
C*---------------------------------------------------------------------*
C*    AUTHOR: GAIL HUTCHINSON                                          *
C*      DATE: 17.11.93                                                 *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
C* 08/08/94  GH     new version started                                *
C* 11/01/95  GH     modified to change the data for the terminii of    *
C*                  strands and strand and loop sequences to be strict *
C*                  extended Kabsch and Sander                         *
c* 14/02/95  GH     add in possibility of analysing more than 1 protein*
c* 03/02/95  GH     adapted to produce either single file for each     *
C*                  protein or all data in a single file               *
c* 09/05/95  GH     corrected definition of 1st of 2 numbers in hairpin*
C*                  class (ixlen)                                      *
c* 27/03/96  GH     reorganise multiple/single file handling           *
c* 24/04/96  GH     small correction in subroutine narrow to correct   *
C*                  for crashing in 3ldh                               *
c* 29/04/96  GH     corrected loop hydrogen bond routine               *
c  29/04/96  GH     increased carr array size to allow up to 50 strands*
c  29/04/96  GH     this version created to allow command line options *
c  13/05/96  GH     corrected to allow up to 8500 residues
c                   small change to definitions of hairpin class to
c                   avoid problems with 1asz
c  21/05/96  GH     allow up to 125 residues in loop
c  05/11/96  GH     allow hairpin lengths of >100 residues (1dot)
c  03/06/97  GH     modified to be strict Kabsch and Sander and to
c                   write out Kabsch and Sander strand terminii to
c                   output file, for consistency with rest of PROMOTIF
c  03/06/97  GH     also modify value written out for length of strands
c  27/08/97  GH     increase size of lhb array to (5,3); print out (4 x 3)
c                   values
c  11/01/99  GH     small modification to line 454 to avoid problems with
c                   zero values of ik1,ik2,jk1,jk2
c  18/10/00  GH     changed array size of sheet to ressht (currently 800)
c   1/11/00  GH     removed 'x' specifier from format statements;changed
c                   reading of map(i,j) to be integer only
c                   made associated changes to phipsi.mat
c                   changed ierror from integer to logical in main
c                   changed to write out 6 lines per hairpin - each <80 characters
c                    long
c  21/11/00  GH     modified to read in new secondary structure format
c  28/11/00  GH     increased maxsht to 150 (formerly 100)
c  28/11/00  GH     changed shpart and abpart 2nd dimension from 5 to 6
c                   to cope with 1cm4
c  28/11/00  GH     increased possible loop length to 300 from 125 residues
c  20/03/01  GH     changed pshtlt to integer in main program
c                   removed double declaration of 'k' in subr result
c  29/08/19  RAL    fixed various compiler errors
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -o p_hairpin3 p_hairpin3.f
C
c***********************************************************************
      INCLUDE 'p_hairpin.par'
      integer hbp(maxres,4),ibp(maxres,2),iseqno(maxres),
     ~     sheet(ressht,6),numres,isht,map(36,36),pnsht,nstr,
     ~     strks(maxstr,2),i,pshtlt(maxsht)
      integer abpart(maxstr,6),ok,fistat,finish,opnera,
     ~     iocode,endnfi,opnerb,ftt,ttt,nsrow,istart,iend,
     ~     nps,NAMLEN,option,sht(maxres),shtks(maxstr,2)
CPDBSUM RAL 24/3/05 -->
      integer   lenstr
CPDBSUM RAL 24/3/05 <--
      real hbe(maxres,4),phi(maxres),psi(maxres)
      character aa(maxres),sstrks(maxres),pname*4,
     ~     sstnam*78,bseqno(maxres)*6,brgl(maxres,2),
     ~     chnbrk(maxres),pdbnam*78,
     ~     hpfil*78,COPTION
      character null*11,hpin(maxstr),pdbcod*80,finam*80
CPDBSUM RAL 24/3/05 -->
      character*512 outnam
CPDBSUM RAL 24/3/05 <--
      LOGICAL ENSEMB,finuse,IERROR
      parameter(finish=2,opnera=1,endnfi=6,opnerb=3,ok=0)
C     ***************
      data fistat/0/,ftt/0/,ttt/0/,FINUSE/.FALSE./,OPTION/1/
C                        ***************
      null=' 999.9 99.9'

      CALL RDPHPS(MAP)
      FISTAT=OK

      READ(*,'(A)') PDBNAM
      IF(PDBNAM(1:1).EQ.'%') THEN
         READ(*,'(A1)') COPTION
         ENSEMB=.TRUE.
         IF(COPTION.EQ.'L'.OR.COPTION.EQ.'l') THEN
            HPFIL='HAIRPINS'
c            OPEN(UNIT=HPINFI,FILE=HPFIL,recl=512)
            OPEN(UNIT=HPINFI,FILE=HPFIL)
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
      ENDIF

 200  CONTINUE
      IF(FISTAT.EQ.OK) THEN
         IF(ENSEMB) THEN
            READ(NAMFI,'(A)',IOSTAT=IOCODE) PDBNAM
            IF(IOCODE.LT.0) THEN
               CLOSE(NAMFI)
               FINUSE=.FALSE.
               FISTAT=ENDNFI
            ELSE
               WRITE(*,*) PDBNAM
            ENDIF
c         ELSE
c            PDBNAM=FINAM
         ENDIF
      ENDIF

      IF(FISTAT.EQ.OK) then
         CALL GETNAM(pdbnam,istart,iend,ierror)
         PDBCOD=PDBNAM(istart:iend)
         NPS=INDEX(pdbcod,' ')
         SSTNAM=PDBCOD(1:nps-1)//'.sst'
         OPEN(UNIT=7,FILE=SSTNAM,STATUS='OLD',IOSTAT=IOCODE)
         IF(IOCODE.NE.0) FISTAT=OPNERB
      ENDIF
      
      IF(FISTAT.EQ.OK) THEN
C         IF(OPTION.ne.2) THEN
         IF(COPTION.NE.'L'.AND.COPTION.NE.'l') then
            HPFIL=PDBCOD(1:NPS-1)//'.hpin'
c            OPEN(unit=hpinfi,file=hpfil,recl=512)
CPDBSUM RAL 14/4/05 -->
C            OPEN(unit=hpinfi,file=hpfil)
            OPEN(unit=hpinfi,status='scratch')
CPDBSUM RAL 14/4/05 <--
c            I=INDEX(pdbnam,'.')
c            if(i.gt.4) then
c               write(hpinfi,123) pdbnam(i-4:i-1)
c            else
c               write(hpinfi,123) pdbnam(1:4)
c            endif
            write(hpinfi,123) pdbcod(1:nps-1)
 123        format('Promotif output for protein :',a)
CPDBSUM RAL 24/3/05 -->
            outnam = "hairpins.promotif"
            open(unit=outpin, file=outnam, status='unknown', err=900)
            call pinhed(outpin)
CPDBSUM RAL 24/3/05 <--
         ENDIF
c     end of change 14.2.95
         call rdsst(pname,numres,iseqno,bseqno,aa,sstrks,sht,
     ~        ibp,hbp,hbe,brgl,phi,psi,chnbrk)
C     C    ***************
c     ZERO ARRAYS FOR WHOLE PROTEIN
         isht=0
         nsrow=0
         do 5 i=1,maxstr
            hpin(i)='N'
 5       continue
         call strands(sstrks,sht,numres,pnsht,pshtlt,strks,shtks,
     ~        nstr)
         call sheets(pnsht,pshtlt,strks,shtks,nstr,ibp,sheet,
     ~        abpart,hbe,hbp,hpin)
         call result(pname,nstr,strks,bseqno,hbp,phi,psi,
     ~        AA,map,sstrks,ibp,chnbrk,ftt,ttt,hpin,pdbnam,coption)
      ELSE IF(FISTAT.EQ.OPNERB) then
         write(*,*) 'Unable to open the 2D structure file ',
     ~        sstnam(1:namlen)
      ENDIF

c     1500 continue

      IF(ENSEMB) THEN
         IF(FISTAT.NE.FINISH) THEN
            FISTAT=OK
            IF(FINUSE) THEN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF
CPDBSUM RAL 24/3/05 -->
      goto 999

C---- Fatal errors
 900  write (*, 7000) outnam(1:lenstr(outnam))
 7000 format ('*** ERROR. Unable to open output file "',A,'"')
      go to 999
CPDBSUM RAL 14/3/05 <--

CPDBSUM RAL 24/3/05 -->
  999 continue
CPDBSUM RAL 24/3/05 <--
      stop
      end
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine strands(sstrks,sht,numres,pnsht,pshtlt,
     ~     strks,shtks,nstr)
c     PURPOSE - TO CALCULATE STRANDS

      INCLUDE 'p_hairpin.par'
      
      integer nchn,n,strks(maxstr,2),nstr,pnsht,
     ~     numres,j,i,sht(maxres),shtks(maxstr,2),
     ~     cursht(maxchn),pshtlt(maxsht)
      character sstrks(maxres),shchr*2
      logical pnew
      data shchr/'eE'/

      pnsht=0
      nstr=0
      do 6 j=1,2
         do 5 i=1,maxstr
            strks(i,j)=0
            shtks(i,j)=0
 5       continue
 6    continue
      nchn=1
      do 9 i=1,12
         cursht(i)=0
 9    continue
      if(index(shchr,sstrks(1)).ne.0) then
         nstr=nstr+1
         strks(nstr,1)=1
      endif
      do 10 i=2,numres
c     RECORD WHETHER NEW CHAIN - FROM BROOKHAVEN CODE
         if(sstrks(i).eq.'E') then
            pnew=.false.
            if(sht(i).ne.0.and.sht(i).ne.cursht(nchn)) then
               pnew=.true.
               do 12 n=1,pnsht
                  if(sht(i).eq.pshtlt(n)) pnew=.false.
 12            continue
               if(pnew) then
                  pnsht=pnsht+1
                  pshtlt(pnsht)=sht(i)
               endif
            endif
         endif
c     CALCULATE STRANDS
c     MODIFY 21.8.92 for CASE OF 'EeE' TO HAVE NO OVERLAP
c     RESIDUE WITH 'e' GOES IN FIRST STRAND
c         if(sstrks(i).eq.'e') then
c            if(sstrks(i-1).eq.'E') strks(nstr,2)=i
c            if(i.lt.numres) then
c               if(sstrks(i-1).eq.'E'.and.sstrks(i+1).eq.'E') then
c                  nstr=nstr+1
c                  strks(nstr,1)=i+1
c                  if(sht(i+1).ne.' ') shtks(nstr,1)=sht(i+1)
c               else if(sstrks(i+1).eq.'E') then
c                  nstr=nstr+1
c                  strks(nstr,1)=i
c                  if(sht(i).ne.' ') shtks(nstr,1)=sht(i)
c               endif
c            endif
c         else if(sstrks(i).eq.'E') then
c            if(index(shchr,sstrks(i-1)).eq.0) then
c               nstr=nstr+1
c               strks(nstr,1)=i
c               if(sht(i).ne.' ') shtks(nstr,1)=sht(i)
c            else if(index(shchr,sstrks(i+1)).eq.0.or.
c     ~              i.eq.numres) then
c               strks(nstr,2)=i
c            else if(sht(i).ne.' ') then
c               if(shtks(nstr,1).eq.' ') then
c                  shtks(nstr,1)=sht(i)
c               else
c                  if(shtks(nstr,1).ne.sht(i)) shtks(nstr,2)=sht(i)
c               endif
c            endif
c         endif
c     MODIFY TO BE STRICT KABSCH AND SANDER 2.6.97
         if(sstrks(i).eq.'E') then
            if(sstrks(i-1).ne.'E') then
               nstr=nstr+1
               strks(nstr,1)=i
               if(sht(i).ne.0) shtks(nstr,1)=sht(i)
            else
               if(sht(i).ne.0) then
                  if(shtks(nstr,1).eq.0) then
                     shtks(nstr,1)=sht(i)
                  else
                     if(shtks(nstr,1).ne.sht(i)) shtks(nstr,2)=sht(i)
                  endif
               endif
            endif
            if(sstrks(i+1).ne.'E'.or.i.eq.numres) then
               strks(nstr,2)=i
            endif
         endif
 10   continue
      return
      end
*****************************************************************
      subroutine sheets(pnsht,pshtlt,strks,shtks,nstr,ibp,
     ~     sheet,abpart,hbe,hbp,hpin)
c     PURPOSE - DEFINE STRANDS AND PARTNERS IN SHEETS
      include 'p_hairpin.par'

      integer abpart(maxstr,6),hbp(maxres,4),i,ibp(maxres,2),
     ~     j,k,l,n,nres,nstr,numstr(maxsht),pnsht,posn(maxres),
     ~     sheet(ressht,6),il,shpart(nsmax,6),shtid,shtstr(maxsht),
     ~     ss,strks(maxstr,2),pshtlt(maxsht),shtks(maxstr,2),shtlt
      character hpin(maxstr)
      real hbe(maxres,4)
     
C     ZERO ARRAYS
      do 10 i=1,pnsht
         numstr(i)=0
 10   continue
      do 15 i=1,maxstr
         do 14 j=1,6
            abpart(i,j)=0
 14      continue
 15   continue
      do 1000 shtid=1,pnsht
         nres=0
         ss=0
         do 20 i=1,maxsht
            shtstr(i)=0
 20      continue
         do 30 i=1,ressht
            do 25 j=1,6
               sheet(i,j)=0
 25         continue
 30      continue
         do 32 i=1,maxres
            posn(i)=0
 32      continue
         shtlt=pshtlt(shtid)
         do 200 j=1,2
            do 190 k=1,nstr
               if(shtks(k,j).eq.shtlt) then
c     COUNT NUMBER OF STRANDS IN THE SHEET
                  numstr(shtid)=numstr(shtid)+1
c     RECORD THE STRANDS WITHIN THE SHEET
                  ss=numstr(shtid)
                  shtstr(ss)=k
C     RECORD THE RESIDUES IN THE SHEET AND THEIR BRIDGE PARTNERS
                  do 100 i=strks(k,1),strks(k,2)
                     nres=nres+1
                     sheet(nres,1)=i
                     posn(i)=nres
                     sheet(nres,4)=ss
                     do 90 l=1,2
                        if(ibp(i,l).gt.0) then
                           n=ibp(i,l)
                      sheet(nres,l+1)=n
                   endif
 90             continue
 100         continue
          endif
 190   continue
 200  continue
C     NOW ALLOCATE STRAND NUMBERS TO PARTNERS
      if(numstr(shtid).gt.1) then
         do 210 i=1,nres
            do 205 j=2,3
               if(sheet(i,j).ne.0) then
                  l=sheet(i,j)
                  il=posn(l)
                  if(il.ne.0) sheet(i,j+3)=sheet(il,4)
               endif
 205        continue
 210     continue
         call combin(sheet,numstr(shtid),nres,
     ~        abpart,shtstr,shpart)      
         call topcal(strks,nres,numstr(shtid),hbe,hbp,
     ~        sheet,shtstr,hpin,shpart)
      endif
 1000 continue
      return
      end        
c**********************************************************
      subroutine combin(sheet,numstr,nres,
     ~     abpart,shtstr,shpart)
c     CALCULATES ALL STRANDS AND PARTNERS

      INCLUDE 'p_hairpin.par'
      
      integer strp,sheet(ressht,6),nres
CPDBSUM RAL 29/8/19 -->
C      integer shtstr(maxstr),numstr,i,j,k,l,ip,np
      integer shtstr(maxsht),numstr,i,j,k,l,ip,np
CPDBSUM RAL 29/8/19 <--
      integer ibuff(2),mc,jmc,carr(nsmax,3)
      integer abpart(maxstr,6),shpart(nsmax,6)

      mc=0
      jmc=0
      do 10 i=1,nsmax
         do 5 j=1,3
            carr(i,j)=0
 5       continue
         do 6 j=1,6
           shpart(i,j)=0
 6      continue
 10   continue
      do 15 i=1,2
         ibuff(i)=0
 15   continue
            
c     COLLECT SERIES OF UNIQUE CONNECTED PAIRS
      do 300 i=1,nres
         do 295 j=5,6
            jmc=0
            ibuff(1)=sheet(i,4)
            ibuff(2)=sheet(i,j)
            if(ibuff(2).ne.0.and.ibuff(2).ne.ibuff(1)) then
               call compar(ibuff,mc,carr,jmc)
               if(jmc.eq.0) then
                  mc=mc+1
                  carr(mc,1)=min(ibuff(1),ibuff(2))
                  carr(mc,2)=max(ibuff(1),ibuff(2))
               else
                  carr(jmc,3)=carr(jmc,3)+1           
               endif
           endif
295         continue
300        continue

c        CALCULATE THE STRAND PARTNERS OF EACH STRAND
        do 350 i=1,numstr
           np=0
           k=shtstr(i)
           do 340 j=1,mc
              do 330 l=1,2
                 if(carr(j,l).eq.i) then
                    np=np+1
                    ip=strp(carr,j,l)
c     shpart contains relative strand numbers within sheet, abpart absolute
c     values
                    shpart(i,np)=ip
                    abpart(k,np)=shtstr(ip)
                 endif
 330          continue
 340       continue
 350    continue
        
        return
        end
C     

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine result(pname,nstr,strks,bseqno,hbp,phi,psi,AA,
     ~     map,sstrks,ibp,chnbrk,ftt,ttt,hpin,pdbnam,coption)
      include 'p_hairpin.par'
   
      integer nstr,strks(maxstr,2),hbp(maxres,4),map(36,36),
     ~     ibp(maxres,2),ftt,ttt,jnar(20,2),lhb(5,3),iloop1,
     ~     iloop2,i,ik1,ik2,jk1,ixlen,iylen,ii,ji,ntime,ils1,ils2,
     ~     ipin,jk2,n,k,j,jsng(3),option,i1,i2,j1,j2
      character pname*4,bseqno(maxres)*6,wknull*6,
     ~     idnull*6,seqs2*(MXSLEN),codnam,aa(maxres),sstrks(maxres),
CPDBSUM RAL 21/4/05 -->
C     ~     phpsn*300,sstn*300,seqn*300,seqs1*(MXSLEN),nulcl*3,type*3,
     ~     phpsn*(MXHLEN),sstn*(MXHLEN),seqn*(MXHLEN),seqs1*(MXSLEN),
     -     nulcl*3,type*3,
CPDBSUM RAL 21/4/05 <--
     ~     nar(20,2)*11,chnbrk(maxres),hpin(maxstr),pdbnam*78,
     ~     pdbcod*4,coption
CPDBSUM RAL 29/3/05 -->
      CHARACTER*80 TMP1
      INTEGER      LENSTR
CPDBSUM RAL 29/3/05 <--
      real phi(maxres),psi(maxres),aphi(5),apsi(5)
CPDBSUM RAL 24/3/05 -->
      CHARACTER*3 RSNAM1, RSNAM2, RSNAM3, RSNAM4
CPDBSUM RAL 24/3/05 <--
      data idnull/'-0000-'/,nulcl/'   '/
      wknull=idnull
      iloop1=0
      iloop2=0

      i=index(pdbnam,'.')
      if(i.gt.4) then
         pdbcod=pdbnam(i-4:i-1)
      else
         if(i.ne.0) then
            pdbcod=pdbnam(1:4)
         else
            k=index(pdbnam,' ')
            pdbcod=pdbnam(k-4:k-1)
         endif
      endif

      do 100 i=1,nstr
         do 1 k=1,3
            jsng(k)=0
 1       continue
c     DEFINE IK1,IK2:START AND END OF CURRENT STRAND; JK1 AND JK2 FOR NEXT STRAND
         ik1=strks(i,1)
         ik2=strks(i,2)
         jk1=strks(i+1,1)
         jk2=strks(i+1,2)
c     IDENTIFY HAIRPINS
         if(hpin(i).eq.'Y') then
         if(bseqno(ik2)(1:1).eq.bseqno(jk1)(1:1)) then
            ixlen=999
            iylen=999
C     DEFINE HAIRPIN CLASS - FIRST DEFINE Y -NUMBER OF RESIDUES BETWEEN 
C     LAST RESIDUE IN FIRST STRAND AND FIRST RESIDUE IN LAST STRAND TO
C     HAVE A BRIDGE PARTNER IN OTHER STRAND
            do 252 ii=ik2,ik1,-1
               do 251 ji=1,2
                  if(ibp(ii,ji).ge.jk1.and.ibp(ii,ji).le.jk2+1) then
                     iylen=ibp(ii,ji)-ii-1
                     iloop1=ii+1
                     iloop2=ibp(ii,ji)-1
                     go to 253
                  endif
 251           continue
 252         continue
C     NOW DEFINE X -NUMBER OF RESIDUES BETWEEN LAST RESIDUE IN FIRST 
c     STRAND AND FIRST RESIDUE IN LAST STRAND TO FORM H BOND TO OTHER STRAND
        NTIME=0
 253        do 255 ii=ik2+1,ik1,-1
c     253        do 255 ii=ik2+1,ik1,-1
c     do 254 ji=1,4
               do 254 ji=1,3,2
c                  if(hbp(ii,ji).ge.jk1.and.hbp(ii,ji).le.jk2) then
                  if(hbp(ii,ji).ge.jk1-1.and.hbp(ii,ji).le.jk2) then
                     if(ji.eq.1.or.ixlen.eq.999) then
                        ixlen=hbp(ii,ji)-ii-1
                     else
                        if((hbp(ii,ji)-ii-1).lt.ixlen)
     ~                      ixlen=hbp(ii,ji)-ii-1
                     endif
                     if(ji.eq.3) go to 256
                  endif
 254           continue
               if(ixlen.ne.999) go to 256
 255        continue
 256        continue
c     NO CLASSIFICATION FOR HAIRPINS WITH CHAIN BREAKS IN THE LOOP
            do 356 ii=iloop1,iloop2
               if(chnbrk(ii).eq.'!') then
                  ixlen=0
                  iylen=0
               endif
 356        continue
            IF(NTIME.EQ.1)THEN
               IK2=IK2-1
               jk1=jk1+1
               GO TO 357
            ENDIF
C*********************
C     22.2.93
            if(ixlen.ne.iylen.and.ixlen.ne.iylen-2.and.ixlen.gt.1) then
c     CLASSIFICATION IS WRONG IN THIS CASE - IYLEN OK BUT REDO IXLEN
               NTIME=1
               IK2=IK2+1
c     changed 13/5/96 to sort out 1asz
               jk1=jk1-1
c     end of change 13/5/96
               GO TO 253
            ENDIF
C*******************
c     CALCULATE NUMBER OF RESIDUES IN TWO STRANDS ACCORDING TO SIBANDA + THORNTON
c 357        ils1=iloop1-ik1
c            ils2=jk2-iloop2
c     change 11.1.95 to be number of residues in strands full stop
 357        ils1=ik2-ik1+1
            ils2=jk2-jk1+1
c     now record sequence and secondary structure for all residues in 'loop'
CPDBSUM RAL 21/4/05 -->
C            do 266 ii=1,300
            do 266 ii=1,MXHLEN
CPDBSUM RAL 21/4/05 <--
               seqn(ii:ii)=' '
               sstn(ii:ii)=' '
               phpsn(ii:ii)=' '
 266        continue
CPDBSUM RAL 21/4/05 -->
C            do 2661 ii=1,30
            do 2661 ii=1,MXSLEN
CPDBSUM RAL 21/4/05 <--
               seqs1(ii:ii)=' '
               seqs2(ii:ii)=' '
 2661       continue
            ipin=1
c3/6/97            do 267 ii=ik2+1,jk1-1
c            do 267 ii=iloop1,iloop2
            do 267 ii=strks(i,2)+1,strks(i+1,1)-1
CPDBSUM RAL 21/4/05 -->
               IF (ipin.le.MXHLEN) THEN
CPDBSUM RAL 21/4/05 <--
                   seqn(ipin:ipin)=aa(ii)
                   sstn(ipin:ipin)=sstrks(ii)
                   call efimov(codnam,phi(ii),psi(ii),map)
                   phpsn(ipin:ipin)=codnam
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
               ipin=ipin+1
 267        continue
            if(ipin.gt.80) write(*,*) 
     ~           "WARNING: LOOP SIZE EXCEEDS 80 - ONLY 1st 80 RESIDUES
     ~           WILL BE RECORDED IN OUTPUT FILE"
            ipin=23
c3/6/97            do 268 ii=ik2,ik1,-1
c            do 268 ii=iloop1-1,ik1,-1
            do 268 ii=strks(i,2),strks(i,1),-1
CPDBSUM RAL 21/4/05 -->
               IF (ipin.GT.0) THEN
CPDBSUM RAL 21/4/05 <--
                   seqs1(ipin:ipin)=aa(ii)
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
               ipin=ipin-1
 268        continue
            ipin=1
c3/6/97            do 269 ii=jk1,jk2
c            do 269 ii=iloop2+1,jk2
            do 269 ii=strks(i+1,1),strks(i+1,2)
CPDBSUM RAL 21/4/05 -->
               IF (ipin.LE.MXSLEN) THEN
CPDBSUM RAL 21/4/05 <--
                   seqs2(ipin:ipin)=aa(ii)
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
               ipin=ipin+1
 269        continue
c     CLASSIFY TYPE 2 TURNS INTO I,IP,II,IIP
            type=nulcl
c********************************************************
c     GET MORE DETAILED INFORMATION ABOUT 2:2,3:5,4:4
            n=iloop1
            if(ixlen.eq.2.and.iylen.eq.2) then
               do 305 n=iloop1,iloop2
                  aphi(n-iloop1+1)=phi(n)
                  apsi(n-iloop1+1)=psi(n)
 305           continue
               call twotwo(phpsn,type)
            else if(ixlen.eq.3.and.iylen.eq.5) then
               do 310 n=iloop1,iloop2
                  aphi(n-iloop1+1)=phi(n)
                  apsi(n-iloop1+1)=psi(n)
 310           continue
               call thrfiv(phpsn,type,ttt)
            else if(ixlen.eq.4.and.iylen.eq.4) then
               do 315 n=iloop1,iloop2
                  aphi(n-iloop1+1)=phi(n)
                  apsi(n-iloop1+1)=psi(n)
 315           continue
               call foufou(phpsn,type,ftt)
            endif
c********************************************************
c     DETERMINE ALL NARROW HYDROGEN BONDING PAIRS IN HAIRPIN
            do 2515 k=1,20
               do 2514 j=1,2
c                  nar(k,j)='----/-0000-'
                   nar(k,j)='           '
                   jnar(k,j)=0
 2514          continue
 2515       continue
c            call narrow(ik1,iloop1-1,iloop2+1,jk2,hbp,nar,bseqno,pname,
c     ~           jnar)
c            call narrow(ik1,ik2,jk1,jk2,hbp,nar,bseqno,pname,
c     ~           jnar,jsng)
            call narrow(strks(i,1),strks(i,2),strks(i+1,1),
     ~           strks(i+1,2),hbp,nar,bseqno,pname,jnar,jsng)
c     DETERMINE ANY LOOP HYDROGEN BONDS CLOSE TO THE STRANDS
c            call loophb(iloop1-1,iloop2+1,hbp,lhb)

c            call loophb(ik2,jk1,hbp,lhb)
            call loophb(strks(i,2),strks(i+1,1),hbp,lhb)
c     WRITE OUT DATA TO HAIRPINS FILE
c     THIS BIT IS NEW 3.6.97
            i1=strks(i,1)
            i2=strks(i,2)
            j1=strks(i+1,1)
            j2=strks(i+1,2)
cccccccccccEND OF NEW BIT
c            IF(OPTION.EQ.2) THEN
            if(COPTION.eq.'L'.or.COPTION.eq.'l') then
c               write(hpinfi,26) pdbcod,i,i+1,bseqno(ik1),
c     ~              bseqno(ik2),bseqno(jk1),bseqno(jk2),ixlen,
c               write(hpinfi,26) pdbcod,i,i+1,bseqno(strks(i,1)),
c     ~              bseqno(strks(i,2)),bseqno(strks(i+1,1)),
c     ~              bseqno(strks(i+1,2)),ixlen,
c     ~              iylen,type,ils1,jk2-jk1+1,seqs1,seqn,seqs2,
               write(hpinfi,26) pdbcod,i,i+1,bseqno(i1),
     ~              bseqno(i2),bseqno(j1),bseqno(j2),
     ~              ixlen,iylen,type,i2-i1+1,j2-j1+1
               write(hpinfi,27) seqs1
               write(hpinfi,28) seqn(1:80)
               write(hpinfi,29) seqs2
               write(hpinfi,30) phpsn(1:80)
               write(hpinfi,31) ((jnar(k,j),j=1,2),k=1,6),
     ~              ((lhb(k,j),j=1,3),k=1,4),(jsng(k),k=1,3)
 26            format(a4,1x,2(i3,1x),4(a6,1x),i3,':',i3,1x,a3,
     ~              1x,2(i2,1x))
 27            format(a30)
 28            format(a80)
 29            format(a30)
 30            format(a80)
 31            format(12(i2,1x),12(i1,1x),3(i2,1x))
            ELSE
c               write(hpinfi,25) i,i+1,bseqno(ik1),bseqno(ik2),
c     ~              bseqno(jk1),bseqno(jk2),ixlen,iylen,type,
c               write(hpinfi,25) i,i+1,bseqno(strks(i,1)),
c     ~              bseqno(strks(i,2)),bseqno(strks(i+1,1)),
c     ~              bseqno(strks(i+1,2)),ixlen,iylen,type,
c     ~              ils1,jk2-jk1+1,seqs1,seqn,seqs2,phpsn,
               write(hpinfi,25) i,i+1,bseqno(i1),bseqno(i2),
     ~              bseqno(j1),bseqno(j2),ixlen,iylen,type,
     ~              i2-i1+1,j2-j1+1
               write(hpinfi,27)seqs1
               write(hpinfi,28)seqn(1:80)
               write(hpinfi,29)seqs2
               write(hpinfi,30)phpsn(1:80)
               write(hpinfi,32)((jnar(k,j),j=1,2),k=1,6),
     ~              ((lhb(k,j),j=1,3),k=1,4),(jsng(k),k=1,3),
     ~              ik1,jk2
 25            format(2(i3,1x),4(a6,1x),i3,':',i3,1x,a3,
     ~              1x,2(i2,1x))
 32            format(12(i2,1x),12(i1,1x),3(i2,1x),i4,1x,i4)

CPDBSUM RAL 24/3/05 -->
C-----         Convert single-letter a.a. codes to 3-letter residue names
               CALL AA2NAM(aa(i1),RSNAM1)
               CALL AA2NAM(aa(i2),RSNAM2)
               CALL AA2NAM(aa(j1),RSNAM3)
               CALL AA2NAM(aa(j2),RSNAM4)

C----          Truncate the first strand sequence and phi-psi codes
               TMP1 = seqn(1:80)

C----          Write out data to .promotif data file
               WRITE(outpin,1120) bseqno(i1)(1:1),
     -             bseqno(i1)(2:6), RSNAM1,
     -             bseqno(i2)(2:6), RSNAM2,
     -             bseqno(j1)(2:6), RSNAM3,
     -             bseqno(j2)(2:6), RSNAM4, ixlen, iylen, type,
     ~             i2-i1+1, j2-j1+1, seqs1(1:LENSTR(seqs1)),
     -             TMP1(1:LENSTR(TMP1)), seqs2(1:LENSTR(seqs2)),
     -             ((jnar(k,j),j=1,2),k=1,6),
     ~             ((lhb(k,j),j=1,3),k=1,4),(jsng(k),k=1,3),
     ~             ik1,jk2
 1120          FORMAT(A1,'::',4(A5,'::',A3,'::'),2(I3,'::'),A3,'::',
     -             2(I2,'::'),3(A,'::'),
     -             12(I2,1X),12(I1,1X),3(I2,1X),I4,1X,I4,'::')
CPDBSUM RAL 24/3/05 <--
            ENDIF
         endif
         endif
 100  continue
      return
      end
c******************************************************************
C
        subroutine compar(ibuff,mc,carr,jmc)

        INCLUDE 'p_hairpin.par'

        integer ibuff(2),mc,jmc,carr(nsmax,3),i
        jmc=0
        do 10 i=1,mc
          if(carr(i,1).eq.ibuff(1).and.carr(i,2)
     1     .eq.ibuff(2)) then
           jmc=i
           go to 11
          else
           if(carr(i,1).eq.ibuff(2).and.carr(i,2)
     1            .eq.ibuff(1)) then
            jmc=i
            go to 11
           endif
          endif
10        continue
11        return
        end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        integer function strp(carr,m,n)

        INCLUDE 'p_hairpin.par'

        integer carr(nsmax,3),m,n
        if(n.eq.2) then
           strp=carr(m,1)
        else
           strp=carr(m,2)
        endif
        return
        end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE EFIMOV (codnam,F,S,map)
c     adapted from louise's program (/errors/progs/db.for)31.10.90

      REAL F, S
      INTEGER I, J, M, MAP(36,36),T
      CHARACTER*1 CODE(8), CODNAM

      DATA(CODE(M),M=1,8)/'O','B','P','A','a','L','G','E'/

C
C ----1.End of chain values removed and counted
C
      IF ((F.EQ.999.9).OR.(S.EQ.999.9)) THEN
        CODNAM='*'
      ELSE
C
C ----2. Phi becomes J
C
         if(f.eq.-180.0) j=1
         if(f.eq.180.0) j=36
         if(f.gt.-180.0.and.f.le.0.0) then
            j=int(f/10)+18
         else if(f.gt.0.0.and.f.lt.180.0) then
            j=int(f/10)+19
            if(mod(f,10.0).eq.0.0)j=j-1
         endif
C
C ----3. Psi becomes I
C
         if(s.eq.180.0) i=1
         if(s.eq.-180.0) i=36
         if(s.lt.180.0.and.s.ge.0.0) then
            i=18-int(s/10)
            if(mod(s,10.0).eq.0.0) i=i+1
         else if(s.lt.0.0.and.s.gt.-180.0) then
            i=abs(int(s/10))+19
         endif
C
C ----Assigns the map array contents to the name array
C
          T=MAP(I,J)
          CODNAM=CODE(T+1)
      ENDIF 
       
      RETURN
      END     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine topcal(strks,nres,numstr,hbe,hbp,
     ~     sheet,shtstr,hpin,shpart)
      
      include 'p_hairpin.par'
      
      integer brp(40,2),hb(40,2),hbp(maxres,4),ihb,m,n,
CPDBSUM RAL 29/8/19 -->
C     ~     nbr,nres,numstr,sheet(ressht,6),shtstr(maxstr),
     ~     nbr,nres,numstr,sheet(ressht,6),shtstr(maxsht),
CPDBSUM RAL 29/8/19 <--
     ~     strks(maxstr,2),shpart(nsmax,6),i,j,k1,k2,jj,jk,np
      real hbe(maxres,4)
      character*1 hpin(maxstr)
      
c     LOOP THROUGH THE STRANDS
      do 100 i=1,numstr
         do 90 j=1,6
            if(shpart(i,j).eq.i+1) then
c     ALSO NEED K1 AND K2 TO BE CONSECUTIVE I.E. CONSECUTIVE STRANDS IN PROTEIN
               k1=shtstr(i)
               k2=shtstr(i+1)
               if(k2.eq.k1+1) then
c     2 CONSECUTIVE STRANDS HYDROGEN BONDED
               ihb=0
               nbr=0
               do 8 jj=1,2
                  do 7 jk=1,40
                     brp(jk,jj)=0
                     hb(jk,jj)=0
 7                continue
 8             continue
c     COLLECT BRIDGES BETWEEN THESE TWO STRANDS
               do 15 m=1,nres
                  if(sheet(m,4).eq.i) then
                     do 14 n=5,6
                        if(sheet(m,n).eq.i+1) then
                           nbr=nbr+1
                           brp(nbr,1)=sheet(m,1)
                           brp(nbr,2)=sheet(m,n-3)
                        endif
 14                  continue
                  endif
 15            continue
c     FIND DIRECTION OF SECOND STRAND

               if(nbr.gt.1) then
                  if(brp(1,2).gt.brp(2,2)) hpin(k1)='Y'
               else
                  do 20 n=strks(k1,1)-1,strks(k1,2)+1
                     do 18 m=1,4
                        if(hbe(n,m).le.-0.5) then
                           np=hbp(n,m)
                           if(np.ge.strks(k2,1)-1.and.
     ~                          np.le.strks(k2,2)+1) then
                              ihb=ihb+1
                              hb(ihb,1)=n
                              hb(ihb,2)=np
                           endif
                        endif
 18                  continue
 20               continue
c                  write(*,*)'only one hydrogen bond'
c                  write(*,*) hb(1,1),hb(1,2)
c                  write(*,*) hb(2,1),hb(2,2)
c     if(hb(2,2).gt.hb(1,2).and.hb(2,1).gt.hb(1,1)) then
c     sgn=cdir(i,1)
c     else if (hb(2,2).lt.hb(1,2).and.hb(2,1).lt.hb(1,1))
c     ~                    then
c     sgn=cdir(i,1)
c     else if(hb(2,2).eq.hb(1,2).and.abs(hb(2,1)-hb(1,1)).
c     ~                    ge.1) then
c     sgn=cdir(i,1)
c     else
c     sgn=oppos(cdir(i,1))
c     endif
               endif
            endif
         endif
 90      continue
 100  continue

      return
      end
C
c*********************************************************************       
      subroutine rdsst(pname,numres,iseqno,bseqno,aa,sstrks,sht,
     ~     ibp,hbp,hbe,brgl,phi,psi,chnbrk)

      INCLUDE 'p_hairpin.par'
      
      integer hbp(maxres,4),ibp(maxres,2),iseqno(maxres),
     ~     numres,m,n,i,sht(maxres)
      real hbe(maxres,4),phi(maxres),psi(maxres)
      character aa(maxres),sstrks(maxres),pname*4,
     ~     bseqno(maxres)*6,brgl(maxres,2),
     ~     chnbrk(maxres)
      
      numres=0
      read(7,*)
      read(7,*)
      read(7,'(1X,A4)') pname
      read(7,'(19X,I4)') numres
      read(7,*)
      read(7,*)
      read(7,*)

c     ZERO RELEVANT PARTS OF ARRAYS      
      do 101 m=1,numres
         iseqno(m)=0
         bseqno(m)='      '
         aa(m)=' '
         sstrks(m)=' '
         sht(m)=0
         do 90 n=1,2
            ibp(m,n)=0
 90      continue
         do 95 n=1,4
            hbp(m,n)=0
            hbe(m,n)=0.0
 95      continue
 101  continue
      do 660 m=1,numres
c         read(7,650,end=1000)iseqno(m),chnbrk(m),bseqno(m),aa(m),
c     ~        sstrks(m),sht(m),brgl(m,1),brgl(m,2),ibp(m,1),ibp(m,2),
c     ~        phi(m),psi(m),(hbp(m,i),hbe(m,i),i=1,4)
c 650     format(1x,i4,a1,a6,7x,a1,5x,a1,2x,
c     1        3a1,5x,2(1x,i4),2(1x,f6.1),26x,4(1x,i4,1x,f4.1))
         read(7,650,end=1000) iseqno(m),chnbrk(m),bseqno(m),aa(m),
     ~        sstrks(m),sht(m),brgl(m,1),brgl(m,2),ibp(m,1),ibp(m,2),
     ~        phi(m),psi(m),(hbp(m,i),hbe(m,i),i=1,4)
 650     format(i5,a1,a6,7x,a1,5x,a1,2x,i3,2a1,5x,2(i5),2(1x,f6.1)
     ~        ,26x,4(i5,1x,f4.1))
 660  continue
 1000 continue
      
      return
      end
c************************************************************
      subroutine rdphps(map)
      INCLUDE 'p_hairpin.par'
      integer map(36,36),i,j
      open(mapfi,file='phipsi.mat',status='old')
C     
C     ----Reads in the digit map array from an external file
C     
      do 10 i=1,36
         read(mapfi,'(36i2)',end=20) (map(i,j),j=1,36)
c         do 5 j=1,36
c            if(map(i,j).eq.'  ') map(i,j)=0
c 5       continue
 10   continue
 20   close(mapfi)
      return
      end
c***************************************************************
      SUBROUTINE LOOPHB(IK2,JK1,HBP,LHB)
c     TO FIND ANY LOOP HYDROGEN BONDS CLOSE TO THE STRANDS

      INCLUDE 'p_hairpin.par'

      INTEGER HBP(MAXRES,4),IK2,JK1,LHB(5,3),IN,J,NHB,i

      DO 5 I=1,5
         DO 4 J=1,3
            LHB(I,J)=0
 4       CONTINUE
 5    CONTINUE
      NHB=0

      DO 20 IN=IK2,IK2+2
         if(in.lt.jk1) then
c      do 20 in=ik2-1,ik2+2
         DO 10 J=1,4
            IF(HBP(IN,J).GE.JK1-2.AND.HBP(IN,J).LT.JK1) THEN
c            if(hbp(in,j).ge.jk1-2.and.hbp(in,j).lt.jk1+1) then
               NHB=NHB+1
               LHB(NHB,1)=IN-IK2
               LHB(NHB,2)=JK1-HBP(IN,J)
               LHB(NHB,3)=MOD(J,2)
            ENDIF
 10      CONTINUE
         endif
 20   CONTINUE
      
      RETURN
      END
C***************************************************************
      subroutine narrow(ik1,ik2,jk1,jk2,hbp,nar,bseqno,pname,jnar,jsng)
c     SUBROUTINE TO DETERMINE THE NARROW HYDROGEN BONDING PARTNERS IN THE
C     HAIRPIN 14.12.92
      include 'p_hairpin.par'

      integer hbp(maxres,4),ik1,ik2,jk1,jk2,inar,narmax,jnar(20,2),
     ~        in,jsng(3),nhres,i
      character bseqno(maxres)*6,nar(20,2)*11,pname*4
      data narmax/0/
      inar=0
      do 10 in=ik1,ik2
         if(hbp(in,1).eq.hbp(in,2).and.hbp(in,1).ge.jk1.and.hbp(in,1).
     ~        le.jk2) then
            call narres(inar,nar,in,hbp(in,1),bseqno,pname,ik2,
     ~           jk1,jnar)
         else if(hbp(in,3).eq.hbp(in,4).and.hbp(in,3).ge.jk1.and.
     ~           hbp(in,3).le.jk2) then
            call narres(inar,nar,in,hbp(in,3),bseqno,pname,ik2,
     ~           jk1,jnar)
         else if(hbp(in,1).eq.hbp(in,4).and.hbp(in,1).ge.jk1.and.
     ~           hbp(in,1).le.jk2) then
            call narres(inar,nar,in,hbp(in,1),bseqno,pname,ik2,
     ~           jk1,jnar)
         else if(hbp(in,2).eq.hbp(in,3).and.hbp(in,2).ge.jk1.and.
     ~           hbp(in,2).le.jk2) then
            call narres(inar,nar,in,hbp(in,2),bseqno,pname,ik2,
     ~           jk1,jnar)
         endif
 10   continue
c     look for a single hydrogen bond at the end
      IF(INAR.eq.0) THEN
         jsng(1)=0
         jsng(2)=0
         jsng(3)=0
c         do 20 i=ik1,ik2-1
c            if(nar(inar,1)(6:11).eq.bseqno(i)) then
c               nhres=i
c               go to 21
c            endif
c 20      continue
         nhres=ik1
 21      do 30 i=nhres,ik2
            if(hbp(i,1).ge.jk1.and.hbp(i,1).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,1)-jk1+1
               jsng(3)=1
            else if(hbp(i,2).ge.jk1.and.hbp(i,2).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,2)-jk1+1
               jsng(3)=2
            else if(hbp(i,3).ge.jk1.and.hbp(i,3).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,3)-jk1+1
               jsng(3)=1
            else if(hbp(i,4).ge.jk1.and.hbp(i,4).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,4)-jk1+1
               jsng(3)=2
            endif
 30      continue
      else if(nar(inar,1)(6:11).ne.bseqno(ik2)) then
         jsng(1)=0
         jsng(2)=0
         jsng(3)=0
         do 40 i=ik1,ik2-1
            if(nar(inar,1)(6:11).eq.bseqno(i)) then
               nhres=i
               go to 41
            endif
 40      continue
 41      do 50 i=nhres,ik2
            if(hbp(i,1).ge.jk1.and.hbp(i,1).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,1)-jk1+1
               jsng(3)=1
            else if(hbp(i,2).ge.jk1.and.hbp(i,2).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,2)-jk1+1
               jsng(3)=2
            else if(hbp(i,3).ge.jk1.and.hbp(i,3).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,3)-jk1+1
               jsng(3)=1
            else if(hbp(i,4).ge.jk1.and.hbp(i,4).le.jk2) then
               jsng(1)=ik2-i+1
               jsng(2)=hbp(i,4)-jk1+1
               jsng(3)=2
            endif
 50      continue
      endif
      if(inar.gt.narmax) then
         narmax=inar
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine narres(inar,nar,in,k,bseqno,pname,ik2,jk1,jnar)
      include 'p_hairpin.par'

      integer inar,in,k,ik2,jk1,jnar(20,2)
      character bseqno(maxres)*6,nar(20,2)*11,pname*4
      
      inar=inar+1
      nar(inar,1)(6:11)=bseqno(in)
      nar(inar,2)(6:11)=bseqno(k)
      nar(inar,1)(1:4)=pname
      nar(inar,2)(1:4)=pname
      jnar(inar,1)=ik2-in+1
      jnar(inar,2)=k-jk1+1
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine twotwo(phpsn,type)
      integer n22,iclass,j
      character phpsn*300,twotyp(19)*2,twocl(19)*3,type*3
      data twotyp/'AA','Aa','aA','aa','LG','LL','GL','GG','LO',
     ~     'OG','GO','OL','PG','PL','PO','EA','Ea','OA','Oa'/
      data twocl/'  I','  I','  I','  I',' IP',' IP',' IP',' IP',
     ~     ' IP',' IP',' IP',' IP',' II',' II',' II','IIP','IIP',
     ~     'IIP','IIP'/

      iclass=0
      type='   '
      n22=0
c     DETERMINE CLASS OF TWO:TWO HAIRPIN
      do 5 j=1,19
         if(phpsn(1:2).eq.twotyp(j)) then
            n22=j
            go to 6
         endif
 5    continue
 6    if(n22.ne.0) type=twocl(n22)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine thrfiv(phpsn,type,ttt)
c     SUBROUTINE TO DEAL WITH 3:5 HAIRPINS
      integer ttt
      character phpsn*300,a*2,gl*2,type*3
      data a/'Aa'/,gl/'GL'/
c      character bp*2
c      data bp/'BP'/

      type='   '
c     CHECK IF IT CORRESPONDS TO TYPE 1 + G1 BULGE
c11.1.95
c      if(index(bp,phpsn(1:1)).ne.0.and.index(a,phpsn(2:2)).ne.0.
c     ~     .and.index(a,phpsn(3:3)).ne.0.and.index(gl,phpsn(4:4))
c     ~     .ne.0.and.index(bp,phpsn(5:5)).ne.0) then
      if(index(a,phpsn(1:1)).ne.0.and.index(a,phpsn(2:2)).ne.0.and.
     ~     index(gl,phpsn(3:3)).ne.0) then
         type=' IG'
         ttt=ttt+1
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine foufou(phpsn,type,ftt)

c     SUBROUTINE TO DEAL WITH 4:4 HAIRPINS
      integer ftt,nj,j
      character phpsn*300,type*3
      
c     DIVIDE INTO TYPE I AND OTHER - TYPE I HAVE EFIMOV REGIONS (AAAL etc)
      nj=0
      type='   '
      do 5 j=1,3
         if(phpsn(j:j).eq.'A'.or.phpsn(j:j).eq.'a') then
            nj=nj+2
         else if(phpsn(j:j).eq.'O') then
            nj=nj+1
         endif
 5    continue
      if(nj.ge.5) then
         if(phpsn(4:4).eq.'G'.or.phpsn(4:4).eq.'L'.or.
     ~        phpsn(4:4).eq.'O') then
            type='  I'
            ftt=ftt+1
         endif
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

c********************************************************************
CPDBSUM RAL 24/3/05 -->
C**************************************************************************
C
C  SUBROUTINE PINHED  -  Write header records to hairpins.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE PINHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 18)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN',
     -     'RESNUM1', 'RESNAM1',
     -     'RESNUM2', 'RESNAM2',
     -     'RESNUM3', 'RESNAM3',
     -     'RESNUM4', 'RESNAM4',
     -     'NHPIN_RES', 'NBRIDGE_RES',
     -     'TYPE',
     -     'NRES_STRAND1', 'NRES_STRAND2',
     -     'STRAND1_SEQ', 'HAIRPIN_SEQ', 'STRAND2_SEQ',
     -     'PLOT_DATA' /

C---- Define field id tag
      IDTAG = 'HP_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:BETA_HAIRPINS'
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
CPDBSUM RAL 24/3/05 <--
