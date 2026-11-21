      PROGRAM sheet_ensem
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: p_sheet_ensem.f                                          *
C*  FUNCTION: generates sheets and strand data for promotif            * 
C* COPYRIGHT: (1991-1996) University College London                    *
C*---------------------------------------------------------------------*
C*    AUTHOR: GAIL HUTCHINSON                                          *
C*      DATE: 1991                                                     *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
c  26/11/91  GH     add an extra column in strands table to
c                   allow for alternative connectivities within barrels
c  21/12/93  GH     from database to promotif format
c  24/10/94  GH     to increase strand array sizes to cope with up to
c                   30 strands per sheet
c  14/02/95  GH     to allow running for more than one protein
c  02/05/95  GH     corrected detection of psi-loops
C  09/05/95  GH     allow print out of single file for multiple protein sets
c                   or one set of files per protein
c  27/03/96  GH     reorganize multiple/single file handling
c  29/04/96  GH     increase size of ibuff arrays
c  29/04/96  GH     this version created to allow command line options
c  11/06/96  GH     increased size of ibuff array in subroutine sins
c  11/06/96  GH     increased maximum number of chains allowed to 25
c  25/09/96  GH     changed number of helical residues in beta alpha beta
c                   units to include those with 'h' and 'H' codes, for 
c                   consistency with helix data
c  28/09/99  GH     changed size of numstr array to maxsht(100) from 62
c c                 increase maxchn to 30 in .par file
c  28/09/99  GH     change array size of ksord to maxk (currently 800)
c  28/09/99  GH     change array size of sheet to ressht (currently 800)
c  28/09/99  GH     changed to generate more data on betaalphabeta and psi loop
c  28/09/99  GH     increased size of ibuff to (80,80)
c  01/11/00  GH     removed 'x' specifier from format statements
c                   changed ierror from integer to logical
c  21/11/00  GH     increased maxsht to 150
c  12/12/00  GH     change ibuff dimensions to bufx x bufy (80 x 120)
c                   modify method for calculating branched topologies
c  09/01/01  GH     modify output format to exclude pdbcode for runs
c                   where this is obvious from the output file name
c=====================================================================
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -o p_sheet3 p_sheet3.f
C
c***********************************************************************
      
      INCLUDE 'p_sheet.par'
      integer abpart(maxstr,6),endnfi,finish,fistat,hbond(4*maxres,2),
     ~     hbp(maxres,4),i,ibp(maxres,2),icnstr(maxstr),iocode,
     ~     irec(maxstr),iseqno(maxres),isht,lenstr(maxstr),namfi,
     ~     namlen,nblg(maxstr),nstr,numpart(maxstr),numres,ok,opnera,
     ~     opnerb,pnsht,posni(maxres,2),sheet(ressht,6),strks(maxstr,2),
     ~     icstr1(maxstr),nsrow,istart,iend,nps,OPTION,sht(maxres),
     ~     shtks(maxstr,2),pshtlt(150)
CPDBSUM RAL 30/3/05 -->
      INTEGER   LENST, OPNERG, OPNERH, OPNERI, OPNERJ,OUTBAB, OUTPSI,
     -          OUTSHT, OUTSTR, OUTSUM
CPDBSUM RAL 30/3/05 <--
      real hbe(maxres,4)
      character aa(maxres),absdir(maxstr,6),brgl(maxres,2),
     ~     bseqno(maxres)*6,constr(maxstr),indirn,
CPDBSUM RAL 21/4/05 -->
C     ~     space*1,sseque(maxstr)*25,sstrks(maxres),sstrau(maxres),
     ~     space*1,sseque(maxstr)*(MXSSEQ),sstrks(maxres),
     -     sstrau(maxres),
CPDBSUM RAL 21/4/05 <--
     ~     coption
      character*(fnamln) pdbcod,pdbnam,finam,sstnam,shtnam,strnam,
     ~     babnam,psinam
      character*(maxk) ksord
CPDBSUM RAL 30/3/05 -->
      CHARACTER*512 OUTNAM
CPDBSUM RAL 30/3/05 <--
      logical caonly,ENSEMB,FINUSE,IERROR
      parameter(endnfi=6,finish=2,indirn='@',namfi=1,ok=0,opnera=1,
CPDBSUM RAL 30/3/05 -->
C     ~     opnerb=3)
     ~     opnerb=3,OPNERG=4,OPNERH=5,OUTSHT=15,OUTSTR=16,OUTBAB=17,
     -         OUTPSI=18,OUTSUM=19)
CPDBSUM RAL 30/3/05 <--

C     ***************
      data space/' '/,FINUSE/.FALSE./,option/1/,fistat/0/
C     ***************
     
      READ(*,'(A)') PDBNAM
      IF(PDBNAM(1:1).EQ.'%') THEN
         READ(*,'(A1)') COPTION
         ENSEMB=.TRUE.
         IF(COPTION.EQ.'L'.OR.COPTION.EQ.'l') then
            SHTNAM='SHEETS'
            STRNAM='STRANDS'
            BABNAM='BETAALPHABETA'
            PSINAM='PSILOOPS'
            open(unit=9,file=SHTNAM,err=1500)
            OPEN(UNIT=10,FILE=STRNAM,ERR=1500)
            OPEN(UNIT=11,FILE=BABNAM,ERR=1500)
            OPEN(UNIT=12,FILE=PSINAM,ERR=1500)
         ENDIF
      ELSE
         ENSEMB=.FALSE.
         COPTION='p'
      ENDIF

      caonly=.false.

C     if ensemble remove % from filename
      IF(ENSEMB) THEN
         FINAM=PDBNAM
         FINUSE=.TRUE.
         NAMLEN=INDEX(PDBNAM,SPACE)-1
         DO 180 I=1,NAMLEN-1
            FINAM(I:I)=FINAM(I+1:I+1)
 180     CONTINUE
         FINAM(NAMLEN:NAMLEN)=SPACE
         OPEN(NAMFI,FILE=FINAM,STATUS='OLD',IOSTAT=IOCODE)
         IF(NAMLEN.LE.0) GO TO 1500
c     FISTAT=FINISH
      ENDIF

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
         PDBCOD=PDBNAM(istart:iend)
         NPS=index(pdbcod,' ')
         SSTNAM=pdbcod(1:nps-1)//'.sst'
         OPEN(unit=7,file=sstnam,status='old',iostat=iocode)
         IF(iocode.ne.0) fistat=opnerb
      ENDIF
      
c     read in data from .sst file
      IF(FISTAT.EQ.OK) THEN
         call rdsst(numres,iseqno,bseqno,aa,sstrks,sht,posni,ibp,
     ~        hbp,hbe,brgl,sstrau,caonly)   
C     if(.not.caonly) then
c     EXTRACT HYDROGEN BONDS
         call hbext(hbond,hbe,numres,hbp)
C     endif
C     
c     ZERO ARRAYS FOR WHOLE PROTEIN
         isht=0
         nsrow=0
         do 1540 i=1,maxstr
            icnstr(i)=0
            constr(i)=' '
 1540    continue
CPDBSUM RAL 15/4/05 -->
C----    Count residues in strands and helices and write out to summary
C        file
         CALL SUMMAR(OUTSUM,BSEQNO,SSTRKS,MAXRES,NUMRES)
CPDBSUM RAL 15/4/05 <--

         call strands(sstrks,sstrau,sht,bseqno,numres,
     ~        aa,brgl,caonly,pnsht,pshtlt,strks,shtks,
     ~        nstr,irec,lenstr,sseque,ksord,nblg)
C     if(.not.caonly.and.nstr.gt.0) then
         IF(NSTR.GT.0) THEN
c            IF(OPTION.ne.2) THEN
            IF(COPTION.NE.'L'.AND.COPTION.NE.'l') then
               SHTNAM=PDBCOD(1:NPS-1)//'.sht'
               STRNAM=PDBCOD(1:NPS-1)//'.str'
CPDBSUM RAL 14/4/05 -->
C               OPEN(unit=9,file=SHTNAM,err=1500)
C               OPEN(UNIT=10,FILE=STRNAM,ERR=1500)
               OPEN(unit=9,status='scratch',err=1500)
               OPEN(UNIT=10,status='scratch',ERR=1500)
CPDBSUM RAL 14/4/05 <--
CPDBSUM RAL 30/3/05 -->
C----       Open sheets data file
            outnam = "sheets.promotif"
            open(unit=outsht, file=outnam, status='unknown',
     -          iostat=iocode)
            if (iocode.ne.0)  fistat = opnerg
            call shthed(outsht)

C----       Open strands data file
            outnam = "strands.promotif"
            open(unit=outstr, file=outnam, status='unknown',
     -          iostat=iocode)
            if (iocode.ne.0)  fistat = opnerh
            call strhed(outstr)

C----       Open beta-alpha-beta data file
            outnam = "bab.promotif"
            open(unit=outbab, file=outnam, status='unknown',
     -          iostat=iocode)
            if (iocode.ne.0)  fistat = opneri
            call babhed(outbab)

C----       Open psi loops data file
            outnam = "psi.promotif"
            open(unit=outpsi, file=outnam, status='unknown',
     -          iostat=iocode)
            if (iocode.ne.0)  fistat = opnerj
            call psihed(outpsi)
CPDBSUM RAL 30/3/05 <--
            endif
c     SHTNAM=pdbcod(1:nps-1)//'.sht'
c     write(9,123) pdbcod(nps-4:nps-1)
c     123        format('Promotif output for protein :',a4)
c     write(9,*)
            call sheets(pnsht,strks,shtks,nstr,ibp,
     ~           sheet,abpart,numpart,absdir,hbe,hbp,icnstr,
CPDBSUM RAL 30/3/05 -->
C     ~           constr,posni,icstr1,pdbnam,coption)
     ~           constr,posni,icstr1,pdbnam,coption,bseqno,outsht)
CPDBSUM RAL 30/3/05 <--
            call result(nstr,strks,bseqno,irec,shtks,lenstr,
     ~           sseque,nblg,ksord,abpart,numpart,absdir,icnstr,
CPDBSUM RAL 30/3/05 -->
C     ~           constr,icstr1,pdbnam,coption)
     ~           constr,icstr1,pdbnam,coption,outstr,aa)
CPDBSUM RAL 30/3/05 <--
            call motifs(strks,bseqno,icnstr,constr,nstr,sstrks,
CPDBSUM RAL 11/4/05 -->
C     ~           lenstr,pdbnam,coption,istart,iend)
     ~           lenstr,pdbnam,coption,istart,iend,outbab,outpsi,aa)
CPDBSUM RAL 11/4/05 <--
         endif
      ELSE IF(FISTAT.EQ.OPNERB) THEN
         write(*,*)'Unable to open the 2D structure file ',
     +        sstnam(1:namlen)
      endif
 1500 continue
      
      IF(ENSEMB) THEN
         if(fistat.ne.finish) then
            fistat=ok
            if(finuse) go to 200
         endif
      ENDIF

      stop
      end
c     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      subroutine topcal(iconn,cdir,strks,nres,numstr,hbe,hbp,
     ~     jjoin,sheet,shtstr,shetyp,posni)
      
      include 'p_sheet.par'
      
      integer brp(40,2),col(nsmax),dst,hb(40,2),hbp(maxres,4),
     ~     ibuff(bufx,bufy),iconn(nsmax,3),ihb,iscen,jjoin(nsmax,3,4)
     ~     ,jx,jys,m,n,nanti,nbr,nparr,nres,numstr,
     ~     posni(maxres,2),rp(40),sheet(ressht,6),shtstr(maxstr),
     ~     strks(maxstr,2),vac,y1(40),y2(40),i,j,k,j1,j2,k1,k2,
     ~     np,len,mindst,mm,l,ik,ii,iscol,jk,jy,ilen,ji,iydel,
     ~     iposn,ll
      real hbe(maxres,4)
      character*1 cdir(nsmax,3),oppos,plus,sgn,sgni,sgnj,
     ~     shetyp
      
      data plus/'+'/
      
      nanti=0
      nparr=0
      do 10 j=1,3
         do 9 i=1,nsmax
            cdir(i,j)=' '
 9       continue
 10   continue
      do 12 j=1,bufy
         do 11 i=1,bufx
            ibuff(i,j)=0
 11      continue
 12   continue
      do 13 i=1,nsmax
         col(i)=0
 13   continue
C     ASSIGN A DIRECTION OF +1 TO THE FIRST STRAND IN ICONN ARRAY
      cdir(1,1)=plus
c     INSERT THIS STRAND IN IBUFF ARRAY
      jx=10
      j=iconn(1,1)
      k=shtstr(j)
      jys=iscen(strks(k,1),strks(k,2))
      call sins(strks(k,1),strks(k,2),jx,jys,ibuff,cdir(1,1),posni)
      col(j)=jx
c     LOOP THROUGH THE REMAINDER OF THE STRANDS
      do 100 i=1,numstr
         ihb=0
         nbr=0
         do 8 j=1,2
            do 7 k=1,40
               brp(k,j)=0
               hb(k,j)=0
 7          continue
 8       continue
         j1=iconn(i,1)
         if(iconn(i+1,1).ne.0) then
            j2=iconn(i+1,1)
         else
            go to 101
         endif
c     COLLECT BRIDGES BETWEEN THESE TWO STRANDS
         do 15 m=1,nres
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
c     FIND DIRECTION OF SECOND STRAND
         k1=shtstr(j1)
         k2=shtstr(j2)
         if(brp(1,2).gt.brp(2,2).and.nbr.gt.1) then
            sgn=oppos(cdir(i,1))
         else if(brp(1,2).le.brp(2,2).and.nbr.gt.1) then
            sgn=cdir(i,1)
         else
            do 20 n=strks(k1,1)-1,strks(k1,2)+1
               do 18 m=1,4
                  if(hbe(n,m).le.-0.5) then
                     np=hbp(n,m)
                     if(np.ge.strks(k2,1)-1.and.
     ~                    np.le.strks(k2,2)+1) then
                        ihb=ihb+1
                        hb(ihb,1)=n
                        hb(ihb,2)=np
                     endif
                  endif
 18            continue
 20         continue
            if(hb(2,2).gt.hb(1,2).and.hb(2,1).gt.hb(1,1)) then
               sgn=cdir(i,1)
            else if (hb(2,2).lt.hb(1,2).and.hb(2,1).lt.hb(1,1))
     ~              then
               sgn=cdir(i,1)
            else if(hb(2,2).eq.hb(1,2).and.abs(hb(2,1)-hb(1,1)).
     ~              eq.2) then
                 sgn=cdir(i,1)
              else
                 sgn=oppos(cdir(i,1))
              endif
           endif
           if(sgn.eq.'+') then
              do 30 m=1,nbr
                 rp(m)=brp(m,2)-strks(k2,1)+1
 30           continue
           else
              do 31 m=1,nbr
                 rp(m)=strks(k2,2)-brp(m,2)+1
 31           continue
           endif
           len=strks(k2,2)-strks(k2,1)+1
           mindst=100
           do 33 m=1,nbr
              mm=brp(m,1)
              y1(m)=posni(mm,2)
 33        continue
           do 40 l=1,60
              do 34 m=1,nbr
                 y2(m)=rp(m)+l
 34           continue
              dst=0
              do 35 m=1,nbr
                 dst=dst+(y2(m)-y1(m))**2
 35           continue
              if(dst.lt.mindst) then
                 mindst=dst
                 jys=l+1
              endif
 40        continue
           cdir(i+1,1)=sgn
           jx=jx+2
           if(jx.gt.bufx) write(*,*) 'ibuff array too small'
           call sins(strks(k2,1),strks(k2,2),jx,jys,ibuff,sgn,posni)
           col(j2)=jx
 100    continue
 101    continue
c     NOW INSERT EXTRA STRANDS FROM JJOIN ARRAY
        do 200 ik=1,3
           do 195 l=1,4
           do 190 i=1,numstr
              if(jjoin(i,ik,l).ne.0) then
                 j1=i
                 j2=jjoin(i,ik,l)
                 nbr=0
                 ihb=0
                 do 110 j=1,2
                    do 105 ii=1,40
                       brp(ii,j)=0
 105                continue
 110             continue
                 do 115 m=1,nres
                    if(sheet(m,4).eq.j1) then
                       do 114 n=5,6
                          if(sheet(m,n).eq.j2) then
                             nbr=nbr+1
                             brp(nbr,1)=sheet(m,1)
                             brp(nbr,2)=sheet(m,n-3)
                          endif
 114                   continue
                    endif
 115             continue
                 iposn=0
c     FIND DIRECTION OF EXISTING STRAND
                 do 117 k=1,numstr
                    do 116 j=1,3
                       if(iconn(k,j).eq.j1) then
                          iposn=k
                          sgni=cdir(k,j)
                          go to 118
                       endif
 116                continue
 117             continue
 118             continue
                 k1=shtstr(j1)
                 k2=shtstr(j2)
                 if(brp(1,2).gt.brp(2,2).and.brp(2,2).ne.0) then
                    sgnj=oppos(sgni)
                    nanti=nanti+1
                 else if(brp(1,2).le.brp(2,2).and.brp(2,2).ne.0) 
     ~                   then
                    sgnj=sgni
                    nparr=nparr+1
                 else
                    do 120 n=strks(k1,1)-1,strks(k1,2)+1
                       do 119 m=1,4
                          if(hbe(n,m).le.-0.5) then
                             np=hbp(n,m)
                             if(np.ge.strks(k2,1)-1.and.
     ~                            np.le.strks(k2,2)+1) then
                                ihb=ihb+1
                                hb(ihb,1)=n
                                hb(ihb,2)=np
                             endif
                          endif
 119                   continue
 120                continue
                    if(hb(2,2).gt.hb(1,2).and.
     ~                   hb(2,1).gt.hb(1,1)) then
                       sgnj=sgni
                       nparr=nparr+1
                    else
                       sgnj=oppos(sgni)
                       nanti=nanti+1
                    endif
                 endif
                 if(sgnj.eq.'+') then
                    do 1100 m=1,nbr
                       rp(m)=brp(m,2)-strks(k2,1)+1
 1100               continue
                 else
                    do 111 m=1,nbr
                       rp(m)=strks(k2,2)-brp(m,2)+1
 111                continue
                 endif
                 len=strks(k2,2)-strks(k2,1)
                 mindst=100
                 do 113 m=1,nbr
                    mm=brp(m,1)
                    y1(m)=posni(mm,2)
 113             continue
                 do 1200 ll=1,60
                    do 1140 m=1,nbr
                       y2(m)=rp(m)+ll
 1140               continue
                    dst=0
                    do 1150 m=1,nbr
                       dst=dst+(y2(m)-y1(m))**2
 1150               continue
                    if(dst.lt.mindst) then
                       mindst=dst
                       jys=ll+1
                    endif
 1200            continue
c     CALCULATE BEST COLUMN TO PUT STRAND IN AND BEST POSITION
c     following line : TEST definition of ilen 18.1.95
                 ilen=strks(k2,2)-strks(k2,1)
                 iscol=col(j1)
                 jk=0
                 jx=0
                 jy=0
                 do 130 k=0,ilen+10
                    do 129 ji=iscol-2,iscol+2,4
                       if(ji.gt.0.and.ji.le.60) then
                       vac=0
                       do 125 jy=jys-2+k,jys+len+k+1
                          if(jy.lt.bufy.and.jy.gt.0) then
                             if(ibuff(ji,jy).ne.0) vac=vac+1
                          endif
 125                   continue
                       if(vac.eq.0) then
                          iydel=k
                          jx=ji
                          go to 131
                       endif
                       vac=0
                       do 127 jy=jys-2-k,jys+len-k+1
                          if(jy.lt.bufy.and.jy.gt.0) then
                             if(ibuff(ji,jy).ne.0) then
                                vac=vac+1
                                endif
                          endif
 127                   continue
                       if(vac.eq.0) then
                          iydel=-k
                          jx=ji
                          go to 131
                       endif
                       endif
 129                continue
 130             continue
 131             continue
                 if(jx.gt.0) then
                 call sins(strks(k2,1),strks(k2,2),jx,jys+iydel,
     ~                ibuff,sgnj,posni)
                 col(j2)=jx
                 if(jx.eq.iscol-2) then
                    if(iposn.gt.1) then
                       if(iconn(iposn-1,2).eq.0) then
                          iconn(iposn-1,2)=j2
                          cdir(iposn-1,2)=sgnj
                       else
                          iconn(iposn-1,3)=j2
                          cdir(iposn-1,3)=sgnj
                       endif
                    else
                       do 150 k=2,1,-1
                          do 149 j=nsmax,2,-1
                             iconn(j,k)=iconn(j-1,k)
                             cdir(j,k)=cdir(j-1,k)
 149                      continue
                          iconn(1,k)=0
                          cdir(1,k)=' '
 150                   continue
                       if(iconn(1,2).eq.0) then
                          iconn(1,2)=j2
                          cdir(1,2)=sgnj
                       else
                          iconn(1,3)=j2
                          cdir(1,3)=sgnj
                       endif
                    endif
                 else
                    if(jx.eq.iscol+2) then
                       if(iconn(iposn+1,2).eq.0) then
                          iconn(iposn+1,2)=j2
                          cdir(iposn+1,2)=sgnj
                       else
                          iconn(iposn+1,3)=j2
                          cdir(iposn+1,3)=sgnj
                       endif
                    endif
                 endif
              endif
              endif
 190          continue
 195          continue
 200       continue
           
c     CALCULATE WHETHER THIS IS A PARALLEL, ANTIPARALLEL OR MIXED SHEET
      do 300 i=2,numstr
         if(cdir(i,1).ne.' ') then
            if(cdir(i,1).eq.cdir(i-1,1)) then
               nparr=nparr+1
            else
               if(cdir(i,1).eq.oppos(cdir(i-1,1))) then
                  nanti=nanti+1
               endif
            endif
         endif
 300  continue
      if(nanti.eq.0) then
         shetyp='P'
      else if(nparr.eq.0) then
         shetyp='A'
      else
         shetyp='M'
      endif
      return
      end
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c
      subroutine sins(m,n,jx,jys,ibuff,dir,posni)

      include 'p_sheet.par'
c     INSERTS STRANDS INTO IBUFF ARRAY
      integer i1,i2,ibuff(bufx,bufy),jx,jy,jys,k,m,n,posni(maxres,2),
     ~     sgn
      character dir*1


      jy=jys
      if(dir.eq.'+') then
         i1=min0(m,n)
         i2=max0(m,n)
         sgn=1
      else
         i1=max0(m,n)
         i2=min0(m,n)
         sgn=-1
      endif
      do 10 k=i1,i2,sgn
         if(jx.gt.0.and.jx.le.bufx.and.jy.gt.0.and.jy.le.bufy) then
            ibuff(jx,jy)=k
            posni(k,1)=jx
            posni(k,2)=jy
         else
            write(*,*) 'warning: trying to insert outside array'
            write(*,9) jx, jy, bufx, bufy, jys
 9          format('          jx, jy, bufx, bufy, jys: ',5I8)
         endif
         jy=jy+1
 10   continue

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer function iscen(m,n)

      include 'p_sheet.par'

      integer m,n
      
      iscen=int(bufy/2-abs((m-n)/2))
      return
      end
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      subroutine sheets(pnsht,strks,shtks,nstr,ibp,
     ~     sheet,abpart,numpart,absdir,hbe,hbp,icnstr,constr,posni,
CPDBSUM RAL 30/3/05 -->
C     ~     icstr1,pdbnam,cOPTION)
     ~     icstr1,pdbnam,cOPTION,bseqno,outsht)
CPDBSUM RAL 30/3/05 <--
      
      include 'p_sheet.par'

      integer abpart(maxstr,6),hbp(maxres,4),i,ibp(maxres,2),
     ~     icnstr(maxstr),iconn(nsmax,3),ijane(nsmax),j,
     ~     jjoin(nsmax,3,4),k,l,mc,n,nres,nstr,numpart(maxstr),
     ~     numstr(maxsht),pnsht,posn(maxres),posni(maxres,2),
     ~     sheet(ressht,6),
     ~     shpart(nsmax,6),shtid,shtstr(maxstr),ss,strks(maxstr,2),
     ~     icstr1(maxstr),il,kk,j1,m,jj,OPTION,ll,nweird,
     ~     adj(maxstr,maxstr),shtks(maxstr,2)
      character absdir(maxstr,6),barrel,cdir(nsmax,3),constr(maxstr),
     ~     coption,dir1,dir2,jane(nsmax),shetyp*1,
CPDBSUM RAL 21/4/05 -->
C     ~     shtlt,strdir*80,strord*120,tplgy*120,pdbcod*4
     ~     shtlt,strdir*(TOPLEN),strord*(TOPLEN),tplgy*(TOPLEN),pdbcod*4
CPDBSUM RAL 21/4/05 <--
      character*(fnamln) pdbnam
      real hbe(maxres,4)
CPDBSUM RAL 30/3/05 -->
      CHARACTER*1   CHAIN, CHAINS(100), LSTCHN
      CHARACTER*6   bseqno(maxres)
      INTEGER       ICHAIN, LEN, LENST, NCHAIN, outsht
      LOGICAL       NEW
CPDBSUM RAL 30/3/05 <--
      
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

      do 10 i=1,pnsht
         numstr(i)=0
 10   continue
      do 15 i=1,maxstr
         do 14 j=1,6
            abpart(i,j)=0
            absdir(i,j)=' '
 14      continue
         numpart(i)=0
         icstr1(i)=0
         icnstr(i)=0
 15   continue
c      write(9,*)
c      write(9,151) pnsht
c 151  format('The protein contains ',i3,' beta sheets')
c     RUN THROUGH ALL SHEETS
      do 1000 shtid=1,pnsht
         nres=0
         ss=0
         mc=0
         nweird=0
         barrel='N'
         do 20 i=1,maxstr
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
CPDBSUM RAL 21/4/05 -->
C         do 34 i=1,120
         do 34 i=1,TOPLEN
CPDBSUM RAL 21/4/05 <--
            strord(i:i)=' '
 34      continue
CPDBSUM RAL 21/4/05 -->
C         do 35 i=1,80
         do 35 i=1,TOPLEN
CPDBSUM RAL 21/4/05 <--
            strdir(i:i)=' '
 35      continue
         do 36 i=1,nsmax
            jane(i)=' '
            ijane(i)=0
 36      continue
c         shtlt=pshtlt(shtid:shtid)
         do 200 j=1,2
            do 190 k=1,nstr
               if(shtks(k,j).eq.shtid) then
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
      do 210 i=1,nres
         do 205 j=2,3
            if(sheet(i,j).ne.0) then
               l=sheet(i,j)
               il=posn(l)
               if(il.ne.0) sheet(i,j+3)=sheet(il,4)
            endif
 205     continue
 210  continue
cccccccccccccccccccccccccccccccccccccc14.11.00      
c     GENERATE ADJACENCY MATRIX FOR SHEET
      do 212 i=1,maxstr
         do 211 j=1,maxstr
            adj(i,j)=0
 211     continue
 212  continue
      do 215 i=1,nres
c         write(*,*) sheet(i,4),sheet(i,5),sheet(i,6)
         do 214 j=1,2
            if(sheet(i,4).ne.0.and.sheet(i,j+4).ne.0) then
               adj(sheet(i,4),sheet(i,j+4))=
     ~              adj(sheet(i,4),sheet(i,j+4))+1
            endif
 214     continue
 215  continue
c      do 217 i=1,24
c         write(*,216)(adj(i,j),j=1,24)
c 216     format(24(i2,1x))
c 217     continue
cccccccccccccccccccccccccccccccccccccc
      call combin(sheet,strks,iconn,numstr(shtid),nres,jjoin,barrel,
     ~     abpart,numpart,shtstr,shpart,nweird)      
CPDBSUM RAL 23/4/05 -->
C      if(numstr(shtid).gt.1.and.nweird.ne.numstr(shtid)) then
      if(numstr(shtid).gt.1.and.nweird.ne.numstr(shtid).and.
     -   iconn(1,1).ne.0) then
CPDBSUM RAL 23/4/05 <--
      call topcal(iconn,cdir,strks,nres,numstr(shtid),hbe,hbp,
     ~     jjoin,sheet,shtstr,shetyp,posni)
      call topology(iconn,cdir,numstr(shtid),barrel,jane
     ~     ,ijane,shtstr,icnstr,constr,icstr1)
      do 300 i=1,numstr(shtid)
         do 290 j=1,6
            if(shpart(i,j).ne.0) then
               do 280 k=1,numstr(shtid)
                  do 270 kk=1,3
                     if(iconn(k,kk).eq.i) then
                        dir1=cdir(k,kk)
                     else
                        if(iconn(k,kk).eq.shpart(i,j))then
                           dir2=cdir(k,kk)
                        endif
                     endif
 270              continue
 280           continue
               j1=shtstr(i)
               if(dir1.eq.dir2) then
                  absdir(j1,j)='+'
               else
                  absdir(j1,j)='-'
               endif
            endif
 290     continue
 300  continue
      n=1
      m=1
      do 400 i=1,nsmax
         j=iconn(i,1)
 301     if(j.ne.0) then
            if(iconn(i,2).eq.0.or.(iconn(i,2).ne.0.and.
     ~           iconn(i,1).eq.0)) then
CPDBSUM RAL 21/4/05 -->
               IF (n + 2.LE.TOPLEN) THEN
CPDBSUM RAL 21/4/05 <--
                   write(strord(n:n+2),'(i3)') shtstr(j)
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
               n=n+3
CPDBSUM RAL 21/4/05 -->
               IF (m + 1.LE.TOPLEN) THEN
CPDBSUM RAL 21/4/05 <--
                   if(iconn(i,2).eq.0) then
                      write(strdir(m:m+1),'(a1)') cdir(i,1)
                   else
                      write(strdir(m:m+1),'(a1)') cdir(i,2)
                   endif
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
               m=m+1
            else
               jj=iconn(i,2)
CPDBSUM RAL 21/4/05 -->
               IF (n + 8.LE.TOPLEN) THEN
CPDBSUM RAL 21/4/05 <--
                   write(strord(n:n+8),900) shtstr(j),shtstr(jj)
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
 900           format('(',i3,',',i3,')')
               n=n+9
CPDBSUM RAL 21/4/05 -->
               IF (m + 4.LE.TOPLEN) THEN
CPDBSUM RAL 21/4/05 <--
                   write(strdir(m:m+4),901) cdir(i,1),cdir(i,2)
 901               format('(',a1,',',a1,')')
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
               m=m+5
               if(iconn(i,3).ne.0) then
CPDBSUM RAL 21/4/05 -->
                   IF (n + 4.LE.TOPLEN) THEN
CPDBSUM RAL 21/4/05 <--
                       write(strord(n:n+4),902) shtstr(iconn(i,3))
CPDBSUM RAL 21/4/05 -->
                   ENDIF
CPDBSUM RAL 21/4/05 <--
 902              format('(',i3,')')
CPDBSUM RAL 21/4/05 -->
                  IF (m + 2.LE.TOPLEN) THEN
CPDBSUM RAL 21/4/05 <--
                      write(strdir(m:m+2),903) cdir(i,3)
 903                  format('(',a1,')')
CPDBSUM RAL 21/4/05 -->
                  ENDIF
CPDBSUM RAL 21/4/05 <--
                  n=n+5
                  m=m+3
               endif
            endif
         else
            if(i.eq.1) then
               j=iconn(i,2)
               go to 301
            endif
         endif
 400  continue
c     WRITE OUT DATA TO SHEETS TABLE
      endif
CPDBSUM RAL 21/4/05 -->
C      do 410 i=2,120,2
      do 410 i=2,TOPLEN,2
CPDBSUM RAL 21/4/05 <--
         if(strord(i:i).eq.'0'.and.strord(i-1:i-1).eq.' ')
     ~        strord(i:i)=' '
 410  continue
      write(tplgy(1:120),'(30(i2,a1,1x))') (ijane(i),jane(i),
     ~     i=1,30)
      do 420 i=2,120,4
         if(tplgy(i:i).eq.'0') then
            if(tplgy(i+1:i+1).eq.' ') then
               tplgy(i:i)=' '
               tplgy(i+1:i+1)=' '
            else if(tplgy(i+1:i+1).eq.'Y') then
               tplgy(i+1:i+1)=' '
            else
               if(tplgy(i+1:i+1).eq.'Z') then
                  tplgy(i+1:i+1)='X'
               endif
            endif
         endif
 420  continue
c      IF(OPTION.EQ.2) THEN
      if(coption.eq.'L'.or.coption.eq.'l') then
         write(9,422) pdbcod,shtid,numstr(shtid),shetyp,barrel
 422     format(a4,1x,i3,1x,i2,1x,a1,1x,a1)
         write(9,'(a4,1x,i3,1x,a120)') pdbcod,shtid,tplgy(1:120)
      ELSE
         write(9,4221) shtid,numstr(shtid),shetyp,barrel
 4221    format(i3,1x,i2,1x,a1,1x,a1,1x)
c         write(9,'(a4,1x,i3,1x,a120)') pdbcod,shtid,tplgy(1:120)
         write(9,'(a120)') tplgy(1:120)
CPDBSUM RAL 30/3/05 -->
C----    Initialise
         NCHAIN = 0
         LSTCHN = ' '

C----    Loop over the strands in this sheet to write out record for each
C        strand belonging to a unique chain
         DO 700, i = 1, numstr(shtid)

C----        Get the strand number
             j = shtstr(i)

C----        Get the strand's start residue
             k = strks(j,1)

C----        Get the chain that this strand belongs to
             CHAIN = bseqno(k)(1:1)
             NEW = .FALSE.

C----        If this is the first chain then treat as new
             IF (NCHAIN.EQ.0) THEN
                 NEW = .TRUE.

C----        If this chain differs from the previous one, check whether
C            it is a new chain
             ELSE IF (CHAIN.NE.LSTCHN) THEN

C----            Assume chain is new
                 NEW = .TRUE.

C----            Check any stored chains to see whether we already have
C                this one
                 DO 680, ICHAIN = 1, NCHAIN
                     IF (CHAINS(ICHAIN).EQ.CHAIN) NEW = .FALSE.
 680             CONTINUE
             ENDIF

C----        If have a new chain, then store in chains list
             IF (NEW) THEN

C----            Increment chain count
                 NCHAIN = NCHAIN + 1

C----            Store chain
                 CHAINS(NCHAIN) = CHAIN
             ENDIF

C----        Save current chain
             LSTCHN = CHAIN
 700     CONTINUE

C----    Loop over the unique chains in this sheet to write out sheet
C        record for each
         DO 800, ICHAIN = 1, NCHAIN

C----        Check that topology is not blank
             LEN = LENST(tplgy)
             IF (LEN.EQ.0) LEN = 1

C----        Write out data to sheets.promotif data file
             WRITE(outsht,720) CHAINS(ICHAIN), shtid, numstr(shtid),
     -           shetyp, barrel, tplgy(1:LEN)
 720         FORMAT(A1,'::',I3,'::',I2,'::',A1,'::',A1,'::',A,'::')
 800     CONTINUE
CPDBSUM RAL 30/3/05 <--
      ENDIF
c      write(9,'(a120)') strord(1:120)
c      write(9,'(a80)')strdir
c      write(9,*)
 1000 continue
      return
      end        
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        subroutine strands(sstrks,sstrau,sht,bseqno,numres,
     ~                     aa,brgl,caonly,pnsht,pshtlt,strks,
     ~                     shtks,nstr,irec,lenstr,sseque,ksord,nblg)

        INCLUDE 'p_sheet.par'
        integer irec(maxstr),lenstr(maxstr),n,nblg(maxstr),nchn,nhel,
     ~       nord,nsht(maxchn),nstr,ntrn,pnsht,resstr,strks(maxstr,2),
     ~       numres,ngel,i,j,sht(maxres),shtks(maxstr,2),
     ~       cursht(maxchn),shtlt(maxchn,maxsht),pshtlt(150)
CPDBSUM RAL 21/4/05 -->
        INTEGER NLOST
CPDBSUM RAL 21/4/05 <--
        character aa(maxres),brgl(maxres,2),bseqno(maxres)*6,
     ~       sp25*25,sp80*80,chnlt(maxchn),
CPDBSUM RAL 21/4/05 -->
C     ~       sscod*11,sseque(maxstr)*25,sstrau(maxres),sstrks(maxres)
     ~       sscod*11,sseque(maxstr)*(MXSSEQ),sstrau(maxres),
     -       sstrks(maxres)
CPDBSUM RAL 21/4/05 <--
c        character sp100*100
       character*(maxk) ksord
        logical caonly,cnew,pnew

        data sp25/'                         '/
        data sp80(1:40)/'                                        '/,
     ~       sp80(41:80)/'                                        '/
c        data sp100(1:40)/'                                         '/,
c     &       sp100(41:80)/'                                        '/,
c     &       sp100(81:100)/'                    '/
        data sscod/'EHJIKGLMNOP'/
        nstr=0
        nhel=0
        ngel=0
        do 6 j=1,2
           do 5 i=1,maxstr
              strks(i,j)=0
              shtks(i,j)=0
 5         continue
 6      continue
        do 8 i=1,maxstr
           sseque(i)=sp25
           lenstr(i)=0
           nblg(i)=0
           irec(i)=0
 8      continue
        do 81 i=1,maxk
           ksord(i:i)=' '
 81     continue
c        ksord(1:100)=sp100
c        ksord(101:200)=sp100
c        ksord(201:300)=sp100
c        ksord(301:400)=sp100
c        ksord(401:403)='   '
        resstr=0
        nord=4
        ntrn=0
        nchn=1
        pnsht=0
CPDBSUM RAL 21/4/05 -->
        NLOST = 0
CPDBSUM RAL 21/4/05 <--
        do 9 i=1,maxchn
           nsht(i)=0
           cursht(i)=0
           chnlt(i)='-'
 9      continue
        if(bseqno(1)(1:1).eq.' ') then
           chnlt(1)='-'
        else
           chnlt(1)=bseqno(1)(1:1)
        endif
        do 10 i=2,numres
c     RECORD WHETHER NEW CHAIN - FROM BROOKHAVEN CODE
           if(bseqno(i)(1:1).ne.bseqno(i-1)(1:1)) then
              nchn=nchn+1        
CPDBSUM RAL 21/4/05 -->
C              if(bseqno(i)(1:1).eq.' ') then
C                 chnlt(nchn)='-'
C              else
C                 chnlt(nchn)=bseqno(i)(1:1)
C              endif
CPDBSUM RAL 21/4/05 <--
           endif
           
c     CALCULATE NUMBER OF RESIDUES IN STRAND/HELIX/G-HELIX
           if(.not.caonly) then
              if(sstrks(i).eq.'E') then
                 resstr=resstr+1
                 cnew=.false.
                 pnew=.false.
                 if(sht(i).ne.0.and.sht(i).ne.cursht(nchn)) then
                    cnew=.true.
                    pnew=.true.
                    do 11 n=1,nsht(nchn)
                       if(sht(i).eq.shtlt(nchn,n)) cnew=.false.
 11                 continue
                    if(cnew) then
                       nsht(nchn)=nsht(nchn)+1
                       shtlt(nchn,nsht(nchn))=sht(i)
                       cursht(nchn)=sht(i)
                    endif
                    do 12 n=1,pnsht
                       if(sht(i).eq.pshtlt(n)) pnew=.false.
 12                 continue
                    if(pnew) then
                       pnsht=pnsht+1
                       pshtlt(pnsht)=sht(i)
                     endif
                 endif
         endif
c CALCULATE STRANDS
         if(sstrks(i).eq.'E') then
              if(sstrks(i-1).ne.'E') then
                  nstr=nstr+1
                  strks(nstr,1)=i
                  lenstr(nstr)=1
                  sseque(nstr)(1:1)=aa(i)
                  if(brgl(i,1).eq.'*'.or.brgl(i,2).eq.'*') 
     ~              nblg(nstr)=nblg(nstr)+1
                  if(sht(i).ne.0) shtks(nstr,1)=sht(i)
              else
                  lenstr(nstr)=lenstr(nstr)+1
                  j=lenstr(nstr)
CPDBSUM RAL 21/4/05 -->
                  if (j.LE.MXSSEQ - 3) then
CPDBSUM RAL 21/4/05 <--
                      sseque(nstr)(j:j)=aa(i)
CPDBSUM RAL 21/4/05 -->
                  else
                      sseque(nstr)(MXSSEQ - 2:MXSSEQ) = '...'
                  endif
CPDBSUM RAL 21/4/05 <--
                  if(brgl(i,1).eq.'*'.or.brgl(i,2).eq.'*') 
     ~             nblg(nstr)=nblg(nstr)+1
                  if(sht(i).ne.0) then
                    if(shtks(nstr,1).eq.0) then
                      shtks(nstr,1)=sht(i)
                     else
                      if(shtks(nstr,1).ne.sht(i)) shtks(nstr,2)=sht(i)
                    endif
                  endif
              endif
              if((sstrks(i+1).ne.'E').or.i.eq.numres) then
                 strks(nstr,2)=i
                 nord=nord+1
CPDBSUM RAL 21/4/05 -->
                 IF (nord.GT.maxk) THEN
                     NLOST = NLOST + 1
                 ELSE
CPDBSUM RAL 21/4/05 <--
                     ksord(nord:nord)='E'
CPDBSUM RAL 21/4/05 -->
                 ENDIF
CPDBSUM RAL 21/4/05 <--
                 irec(nstr)=nord
              endif
           endif
           
cc    TEST FOR KS+ HELIX RESIDUES
c-----------------------------------
         if(sstrks(i).eq.'H') then
            if((sstrks(i+1).ne.'H').or.i.eq.numres) then
               nhel=nhel+1
               nord=nord+1
CPDBSUM RAL 21/4/05 -->
               IF (nord.GT.maxk) THEN
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   ksord(nord:nord)='H'
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            endif
         endif
c     3,10 helices
c---------------
         
         if(sstrks(i).eq.'G') then
            if(sstrks(i+1).ne.'G'.or.i.eq.numres) then
               ngel=ngel+1
               nord=nord+1
CPDBSUM RAL 21/4/05 -->
               IF (nord.GT.maxk) THEN
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   ksord(nord:nord)='G'
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            endif
         endif
         
c     Turns
c---------------
         if(sstrks(i).eq.'T') then
            if(sstrks(i+1).ne.'T'.or.i.eq.numres) then
               ntrn=ntrn+1
               nord=nord+1
CPDBSUM RAL 21/4/05 -->
               IF (nord.GT.maxk) THEN
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   ksord(nord:nord)='T'
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            endif
         endif
c     THE FOLLOWING ENDIF APPLIES TO END OF INSTRUCTIONS FOR NON CA-ONLY FILES
      endif
 10   continue

CPDBSUM RAL 21/4/05 -->
C---- If maximum number of sec str elements exceeded, show warning
      IF (NLOST.GT.0) THEN
          PRINT 103, maxk, NLOST
 103      FORMAT('Warning. Maximum number of sec struc elements ',
     -         '(',I6,') exceeded. ',I6,' elements lost')
      ENDIF
CPDBSUM RAL 21/4/05 <--
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine result(nstr,strks,bseqno,irec,shtks,lenstr,
     ~      sseque,nblg,ksord,abpart,numpart,absdir,icnstr,constr,
CPDBSUM RAL 30/3/05 -->
C     ~      icstr1,pdbnam,cOPTION)
     ~      icstr1,pdbnam,cOPTION,outstr,aa)
CPDBSUM RAL 30/3/05 <--
       
       include 'p_sheet.par'
       
       integer abpart(maxstr,6),i,icnstr(maxstr),ik1,ik2,
     ~      irec(maxstr),lenstr(maxstr),n,nblg(maxstr),
     ~      nstr,numpart(maxstr),strks(maxstr,2),strtot,
     ~      icstr1(maxstr),j,m,OPTION,shtks(maxstr,2)
       character absdir(maxstr,6),bseqno(maxres)*6,constr(maxstr),
     ~      idnull*6,
CPDBSUM RAL 21/4/05 -->
C     ~      sseque(maxstr)*25,wknull*6,coption,pdbcod*4
     ~      sseque(maxstr)*(MXSSEQ),wknull*6,coption,pdbcod*4
CPDBSUM RAL 21/4/05 <--
       character*(fnamln) pdbnam
       character*(maxk) ksord
CPDBSUM RAL 30/3/05 -->
       CHARACTER*1   aa(maxres)
       CHARACTER*3   RSNAM1, RSNAM2
       CHARACTER*4   LEFT, RIGHT
       INTEGER       LENST, outstr
CPDBSUM RAL 30/3/05 <--
       
       data idnull/'-0000-'/
       data strtot/0/
       
       wknull=idnull
       
       i=index(pdbnam,'.')
       if(i.gt.4) then
          pdbcod=pdbnam(i-4:i-1)
       else
          if(i.ne.0) then
             pdbcod=pdbnam(1:4)
          else
             j=index(pdbnam,' ')
             pdbcod=pdbnam(j-4:j-1)
          endif
       endif
C     WRITE OUT STRAND DATA
c       write(9,*)
c       write(9,2) nstr
c 2     format('The protein contains ',i3,' beta strands:')
c       write(9,3)
 3     format
     ~(' NUM   BEG    END SHT LEN SEQUENCE',15x,'NBLG NPRT PARTNERS')
       do 100 i=1,nstr
          if(icnstr(i).eq.0) then
             if(constr(i).eq.'Y') then
                constr(i)=' '
             else if(constr(i).eq.'Z') then
                constr(i)='X'
             else
                if(constr(i).eq.' ') then
                   icnstr(i)=99
                endif
             endif
          endif
c     IK1 AND IK2 ARE BEGINNING AND END OF STRAND
          ik1=strks(i,1)
          ik2=strks(i,2)
          strtot=strtot+1
          n=irec(i)
c          if(OPTION.eq.2) then
CPDBSUM RAL 21/4/05 -->
          j = 0
          DO 108, m = n - 4, n - 1
              j = j + 1
              IF (m.LT.0 .OR. m.GT.maxk) THEN
                  LEFT(j:j) = ' '
              ELSE
                  LEFT(j:j) = ksord(m:m)
              ENDIF
 108      CONTINUE
          j = 0
          DO 109, m = n + 1, n + 4
              j = j + 1
              IF (m.LT.0 .OR. m.GT.maxk) THEN
                  RIGHT(j:j) = ' '
              ELSE
                  RIGHT(j:j) = ksord(m:m)
              ENDIF
 109          CONTINUE
CPDBSUM RAL 21/4/05 <--
          if(coption.eq.'L'.or.coption.eq.'l') then
             write(10,31) pdbcod,i,bseqno(ik1),bseqno(ik2),
     ~            shtks(i,1),shtks(i,2),lenstr(i),sseque(i),
     ~            numpart(i),(abpart(i,j),j=1,4),
     ~            (absdir(i,j), j=1,4), icnstr(i),constr(i),icstr1(i),
CPDBSUM RAL 21/4/05 -->
C     ~            constr(i),(ksord(m:m),m=n-4,n-1),
C     ~            (ksord(m:m),m=n+1,n+4)
     ~            constr(i), LEFT, RIGHT
CPDBSUM RAL 21/4/05 <--
          else
             write(10,30) i,bseqno(ik1),bseqno(ik2),
     ~            shtks(i,1),shtks(i,2),lenstr(i),sseque(i),
     ~            numpart(i),(abpart(i,j),j=1,4),
     ~            (absdir(i,j),j=1,4),icnstr(i),constr(i),icstr1(i),
CPDBSUM RAL 21/4/05 -->
C     ~            constr(i),(ksord(m:m),m=n-4,n-1),
C     ~            (ksord(m:m),m=n+1,n+4),ik1,ik2
     ~            constr(i), LEFT, RIGHT, ik1,ik2
CPDBSUM RAL 21/4/05 <--
CPDBSUM RAL 30/3/05 -->
C----        Convert single-letter a.a. codes to 3-letter residue names
             CALL AA2NAM(aa(ik1),RSNAM1)
             CALL AA2NAM(aa(ik2),RSNAM2)

C----        Write out data to sheets.promotif data file
             WRITE(outstr,640) bseqno(ik1)(1:1), i,
     -           bseqno(ik1)(2:6), RSNAM1, bseqno(ik2)(2:6), RSNAM2,
     -           shtks(i,1), shtks(i,2), lenstr(i),
     -           sseque(i)(1:LENST(sseque(i))), numpart(i),
     -           (abpart(i,j), j=1,4), (absdir(i,j), j=1,4), icnstr(i),
     -           constr(i), icstr1(i), LEFT, RIGHT
 640         FORMAT(A1,'::',I3,'::',2(A5,'::',A3,'::'),2(I3,'::'),I3,
     -           '::',A,'::',I3,'::',4(I3,'::'),4(A1,'::'),I2,'::',A1,
     -           '::',I2,'::',2(A4,'::'))
CPDBSUM RAL 30/3/05 <--
          endif
 100      continue
 31       format(a4,1x,i4,1x,2(a6,1x),2(i3,1x),1x,i2,2x,a25,
     ~      1x,i1,4(1x,i3),4a1,i2,a1,i2,a1,1x,8a1)
 30       format(i4,1x,2(a6,1x),2(i3,1x),1x,i2,2x,a25,
     ~      1x,i1,4(1x,i3),4a1,i2,a1,i2,a1,1x,8a1,1x,i4,1x,i4)
       return
       end
c******************************************************************
       subroutine combin(sheet,strks,iconn,numstr,nres,
     ~      jjoin,barrel,abpart,numpart,shtstr,shpart,nweird)
       
       INCLUDE 'p_sheet.par'
       
       integer abpart(maxstr,6),barr,carr(nsmax,3),freq(nsmax),i,
     ~      ibuff(2),icol(nsmax),iconn(nsmax,3),ie,iord,istr,
     ~      iuse(nsmax),j,jcol(nsmax),jconn(nsmax,3),jjoin(nsmax,3,4),
     ~      jmc,jp(20),jr(2),juse(nsmax),lplen,mc,minlen,n,nres,
     ~      numpart(maxstr),sheet(ressht,6),shpart(nsmax,6),
     ~      shtstr(maxstr),
     ~      strks(maxstr,2),strp,numstr,nend,iend,islen,np,k,l,ip,il,m,
     ~      ibnd,id1,id2,id,mk,nweird
       character barrel

       mc=0
       jmc=0
       ie=0
       barr=0
       minlen=100
       lplen=0
       do 10 i=1,nsmax
          freq(i)=0
          iuse(i)=0
          juse(i)=0
          icol(i)=0
          jcol(i)=0
          do 5 j=1,3
             iconn(i,j)=0
             jconn(i,j)=0
             carr(i,j)=0
 5        continue
          do 6 j=1,6
             shpart(i,j)=0
 6        continue
 10    continue
       do 15 i=1,2
          ibuff(i)=0
 15    continue
       do 20 j=1,nsmax
          do 19 l=1,4
             do 18 i=1,3
                jjoin(j,i,l)=0
 18          continue
 19       continue
 20    continue
       
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
 295      continue
 300   continue
       nend=0
       iend=0
       islen=0
       nweird=0
c     CALCULATE THE STRAND PARTNERS OF EACH STRAND
       do 350 i=1,numstr
          np=0
          k=shtstr(i)
          do 340 j=1,mc
             do 330 l=1,2
                if(carr(j,l).eq.i) then
                   freq(i)=freq(i)+1
                   np=np+1
                   ip=strp(carr,j,l)
c     shpart contains relative strand numbers within sheet, abpart absolute
c     values
                   shpart(i,np)=ip
                   abpart(k,np)=shtstr(ip)
                endif
 330         continue
 340      continue
          numpart(k)=np
c***  LOOK FOR END STRAND
          if(np.eq.1) then
             il=strks(k,2)-strks(k,1)
             nend=nend+1
             if(il.gt.islen.or.islen.eq.0) then
                islen=il
                iend=i
             endif
             else if(np.eq.0) then
                nweird=nweird+1
          endif
 350   continue
c     28.9.99 IF ALL STRANDS HAVE NO PARTNERS, THEN WEIRD SHEET, DON'T
c     CALCULATE (e.g. 1wat)
       if(nweird.eq.numstr) then
          write(*,*) 'Warning: Weird sheet, no strands have partners!!'
          go to 501
       else if (nweird.ne.0) then
          write(*,*) 'Warning: some strands have no partners'
       endif
c***  DEFINE BARRELS AS SHEETS WITH LESS THAN TWO END STRANDS
       if(nend.lt.2) barrel='Y'
       iord=1
       if(iend.ne.0) then
          iconn(1,1)=iend
          icol(iend)=iord
          istr=1
          iuse(iend)=1
          call conect(iord,carr,iconn,mc,iuse,istr,icol,numstr
     ~         ,jjoin,freq,strks)
       else
c     BARREL
          do 500 i=1,numstr
             if(freq(i).eq.2) then
                minlen=100
                n=0
                do 360 j=1,mc
                   do 355 k=1,2
                      if(carr(j,k).eq.i) then
                         n=n+1
                         jp(n)=strp(carr,j,k)
                         jr(n)=j
                      endif
                      if(n.eq.2) go to 365
 355               continue
 360            continue
 365            continue
                do 400 j=1,2
                   do 370 m=1,nsmax
                      iuse(m)=0
                      juse(m)=0
                      jcol(m)=0
                      icol(m)=0
 370               continue
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
                     call conect(iord,carr,jconn,mc,juse,istr,jcol,
     ~                    numstr,jjoin,freq,strks)
                     do 380 n=1,numstr-1
                        id1=jcol(n)
                        id2=jcol(n+1)
                        id=abs(id2-id1)
                        lplen=lplen+id
 380                 continue
                     lplen=lplen+ibnd
                     if(lplen.lt.minlen) then
                        minlen=lplen
                        do 390 k=1,nsmax
                           do 385 mk=1,3
                              iconn(k,mk)=jconn(k,mk)
 385                       continue
 390                    continue
                        iconn(istr+1,1)=iconn(1,1)
                     endif                
 400              continue
               endif
 500        continue
         endif        
 501     continue
         return
         end
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine conect(iord,carr,iconn,mc,iuse,istr
     1        ,icol,numstr,jjoin,freq,strks)
         
         INCLUDE 'p_sheet.par'
         
         integer carr(nsmax,3),freq(nsmax),i,icol(nsmax),iconn(nsmax,3),
     ~        idone(nsmax),iord,ip,istr,iuse(nsmax),j,jjoin(nsmax,3,4),
     ~        ki,kk,kb(nsmax,3),mc,mconn,n,nconn,noption,numstr,
     ~        strks(maxstr,2),strp,k,nr,m,ii,jlen,ilen,kj,l

 1       kk=0
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
 9          continue
 10      continue
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
C     CHECK WHETHER ANY OF THESE PARTNER STRANDS IS AT THE END OF A BRANCH
            noption=0
            ki=kb(1,1)
            do 11 i=1,kk
               if(kb(i,3).gt.1) then
                  noption=noption+1
                  ki=kb(i,1)
               endif
 11         continue
c     IF >1 OF THESE NOT AT END OF BRANCH SEE WHICH ONE TO USE
            if(noption.gt.1) then
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
 13                 continue
 14              continue
 116          continue
              jlen=strks(m,2)-strks(m,1)
                 if(nconn.gt.mconn.or.mconn.eq.0.or.
     ~                (nconn.eq.mconn.and.jlen.gt.ilen)) then
                    mconn=nconn
                    ki=kb(i,1)
                    ilen=strks(ki,2)-strks(ki,1)
                    kj=kb(i,2)
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
           DO 20 I=1,numstr
              IF(IUSE(I).EQ.0) THEN
c     find a strand not connected and find its partner
                 IP=0
                 DO 17 M=1,MC
                    DO 16 N=1,2
                       IF(CARR(M,N).EQ.I) ip=strp(carr,m,n)
                       IF(IP.NE.0) THEN
                          IF(IUSE(IP).EQ.1) THEN
                             DO 1116 L=1,4
                                IF(JJOIN(IP,1,L).EQ.0) THEN
                                   JJOIN(IP,1,L)=I
                                   GO TO 1117
                                ENDIF
 1116                        CONTINUE
 1117                        CONTINUE
                             IUSE(I)=2
                             ISTR=ISTR+1
                             IF(ISTR.EQ.numstr) THEN
                                GO TO 50
                             ELSE
                                GO TO 19
                             ENDIF
                          ENDIF
                       ENDIF
 16                 CONTINUE
 17              CONTINUE
 19              CONTINUE
              ENDIF
 20        CONTINUE
c     second level of branching
           DO 30 I=1,numstr
              IF(IUSE(I).EQ.0) THEN
c     find a strand not connected + find its partner
                 IP=0
                 DO 27 M=1,MC
                    DO 26 N=1,2
                       IF(CARR(M,N).EQ.I) ip=strp(carr,m,n)
                       IF(IP.NE.0) THEN
                          IF(IUSE(IP).EQ.2) THEN
                             DO 1119 L=1,4
                                IF(JJOIN(IP,2,L).EQ.0) THEN
                                   JJOIN(IP,2,L)=I
                                   GO TO 1120
                                ENDIF
 1119                        CONTINUE                            
 1120                        CONTINUE
                             IUSE(I)=3
                             ISTR=ISTR+1
                             GO TO 29
                          ELSE
                             IF(IUSE(IP).EQ.3) THEN
                                DO 1121 L=1,4
                                IF(JJOIN(IP,3,L).EQ.0) THEN
                                   JJOIN(IP,3,L)=I
                                   GO TO 1122
                                ENDIF
 1121                        CONTINUE                            
 1122                        CONTINUE
                                IUSE(I)=3
                                ISTR=ISTR+1
                                GO TO 29
                             ENDIF
                          ENDIF
                       ENDIF
 26                 CONTINUE
 27              CONTINUE
              ENDIF
 29           CONTINUE
 30        CONTINUE
           
C     THIRD LEVEL OF BRANCHING
           DO 40 I=1,numstr
              IF(IUSE(I).EQ.0) THEN
c     find a strand not connected and find its partner
                 IP=0
                 DO 37 M=1,MC
                    DO 36 N=1,2
                       IF(CARR(M,N).EQ.I) ip=strp(carr,m,n)
                       IF(IP.NE.0) THEN
                          IF(IUSE(IP).EQ.3) THEN
                                DO 1123 L=1,4
                                IF(JJOIN(IP,3,L).EQ.0) THEN
                                   JJOIN(IP,3,L)=I
                                   write(*,*) ip,l,'THREE','JJOIN'
                                   GO TO 1124
                                ENDIF
 1123                        CONTINUE                            
 1124                        CONTINUE
                             IUSE(I)=4
                             ISTR=ISTR+1
                             GO TO 39
                          ENDIF
                       ENDIF
 36                 CONTINUE
 37              CONTINUE
              ENDIF
 39           CONTINUE
 40        CONTINUE
        endif
 50     return
        end
        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        subroutine compar(ibuff,mc,carr,jmc)

        INCLUDE 'p_sheet.par'

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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        integer function strp(carr,m,n)
        include 'p_sheet.par'

        integer carr(nsmax,3),m,n

        if(n.eq.2) then
           strp=carr(m,1)
        else
           strp=carr(m,2)
        endif
        return
        end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine topology(iconn,cdir,numstr,barrel,jane
     ~                ,ijane,shtstr,icnstr,constr,icstr1)

        INCLUDE 'p_sheet.par'

        integer i,icnstr(maxstr),iconn(nsmax,3),iend,ijane(nsmax),ipsi,
     ~       j,join,posn(nsmax,2),shtstr(maxstr),toppos,
     ~       icstr1(maxstr),numstr,k,l,nalt,join1,join2,kk
        character barrel,cdir(nsmax,3),cdir1,cdir2,constr(maxstr),
     ~       cposn(nsmax),cross,jane(nsmax)
c     
        do 5 i=1,nsmax
           jane(i)=' '
           ijane(i)=0
           do 4 j=1,2
              posn(i,j)=0
 4         continue
          cposn(i)=' '
 5     continue
       toppos=1
       
c     CALCULATE THE LAST OCCUPIED STRAND COLUMN IN THIS SHEET
        do 10 j=1,nsmax
           if(iconn(j,1).eq.0) then
              iend=j-2
              go to 15
          endif
 10    continue
 15    continue
c        NOW CALCULATE THE DISTANCE OF EACH STRAND FROM EACH END OF ITS SHEET
       do 50 j=1,numstr
          do 40 k=1,nsmax
             do 30 l=1,3
                if(iconn(k,l).eq.j) then
                posn(j,1)=k
                posn(j,2)=iend+1-k
                cposn(j)=cdir(k,l)
c     COLUMN 1:DISTANCE FROM 1 END; COLUMN 2:DISTANCE FROM OTHER END; 
C     COLUMN 3 ORIENTATION OF STRAND
                go to 35
              endif
 30        continue
 35        continue
 40     continue
 50   continue
c     
      do 70 j=2,numstr
         ipsi=0
         nalt=0
         cdir1=cposn(j)
          cdir2=cposn(j-1)
          if(barrel.eq.'N') then
             join=posn(j,1)-posn(j-1,1) 
          else
c        FOR BARRELS TAKE SMALLER OF TWO POSSIBILITIES        
             if(posn(j,1).gt.posn(j-1,1)) then
                join1=posn(j,1)-(posn(j-1,1)+iend)
             else
c                join1=posn(j-1,1)-(posn(j,1)+iend)
                join1=posn(j,1)+iend-posn(j-1,1)
             endif
             join2=posn(j,1)-posn(j-1,1)
             if(abs(join1).lt.abs(join2)) then
                join=join1
             else
                join=join2
                if(abs(join1).eq.abs(join2)) nalt=1
             endif
          endif
            if(join.eq.0) then
               if(cdir1.eq.cdir2) then
                  cross='Y'
               else
                  cross='Z'
               endif
            else
               if(cdir1.eq.cdir2) then 
                  cross='X'
               else
                  cross=' '
               endif
            endif
                toppos=toppos+1
            ijane(toppos)=join
            jane(toppos)=cross
            if(shtstr(j).eq.shtstr(j-1)+1) then
               kk=shtstr(j-1)
               if(nalt.ne.1) then
                  icnstr(kk)=join
               else
                  icnstr(kk)=join1
                  icstr1(kk)=join2
               endif
               constr(kk)=cross
            endif
 70      continue
         return
         end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        character function oppos(c)
        character c

        if(c.eq.'+') then
           oppos='-'
        else
           if(c.eq.'-') then
              oppos='+'
           endif
        endif
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdsst(numres,iseqno,bseqno,aa,sstrks,sht,posni,ibp,
     ~     hbp,hbe,brgl,sstrau,caonly)

      INCLUDE 'p_sheet.par'
      
      integer hbp(maxres,4),ibp(maxres,2),iseqno(maxres),
     ~     posni(maxres,2),numres,m,n,i,sht(maxres)
      real hbe(maxres,4)
      character aa(maxres),sstrks(maxres),
     ~     bseqno(maxres)*6,brgl(maxres,2),sstrau(maxres)
      logical caonly

      numres=0
      read(7,*)
      read(7,*)
      read(7,*)
      read(7,'(19X,I4)') numres
      read(7,*)
      read(7,*)
      read(7,*)
c     ZERO RELEVANT PARTS OF ARRAYS
CPDBSUM RAL 12/4/05 -->
      IF (numres.gt.maxres - 1) THEN
          PRINT 110, numres, maxres
 110      FORMAT('*** ERROR. Max. no. of residues (maxres) exceeded: ',
     -        I10,' (',I6,')')
          STOP
      ENDIF
CPDBSUM RAL 12/4/05 <--

      do 101 m=1,numres
         iseqno(m)=0
         bseqno(m)='      '
         aa(m)=' '
         sstrks(m)=' '
         sstrau(m)=' '
         sht(m)=0
         do 90 n=1,2
            posni(m,n)=0
            ibp(m,n)=0
 90      continue
         do 95 n=1,4
            hbp(m,n)=0
            hbe(m,n)=0.0
 95      continue
 101  continue
      if(.not.caonly) then
         do 660 m=1,numres
c            read(7,650,end=1000)iseqno(m),bseqno(m),aa(m),sstrau(m),
c     ~           sstrks(m),sht(m),brgl(m,1),brgl(m,2),
c     ~           ibp(m,1),ibp(m,2),(hbp(m,i),hbe(m,i),i=1,4)
c 650        format(1x,i4,1x,a6,7x,a1,1x,a1,3x,a1,2x,
c     ~           3a1,5x,2(1x,i4),40x,4(1x,i4,1x,f4.1))
           read(7,650,end=1000)iseqno(m),bseqno(m),aa(m),sstrau(m),
     ~           sstrks(m),sht(m),brgl(m,1),brgl(m,2),
     ~           ibp(m,1),ibp(m,2),(hbp(m,i),hbe(m,i),i=1,4)
 650        format(i5,1x,a6,7x,a1,1x,a1,3x,a1,2x,
     ~           i3,2a1,5x,2(i5),40x,4(i5,1x,f4.1))
 660     continue
 1000    continue
      else
         do 680 m=1,numres
            read(7,670,end=690) iseqno(m),bseqno(m),aa(m),sstrau(m)
 670        format(2x,i3,1x,a6,7x,a1,3x,a1)
 680     continue
 690     continue
      endif
      return
      end

c*************************************************
c     subroutine motifs - test routine added 5.10.94 to search for
c     simple motifs
c**************************************************
      subroutine motifs(strks,bseqno,icnstr,constr,nstr,sstrks,
CPDBSUM RAL 11/4/05 -->
C     ~     lenstr,pdbnam,cOPTION,istart,iend)
     ~     lenstr,pdbnam,cOPTION,istart,iend,outbab,outpsi,aa)
CPDBSUM RAL 11/4/05 <--
      include 'p_sheet.par'

CPDBSUM RAL 12/4/05 -->
      INTEGER MAXHEL
      PARAMETER (MAXHEL=20)
CPDBSUM RAL 12/4/05 <--

      integer nbab,nps,npsi,nstr,icnstr(maxstr),i1,i2,j1,j2,
     ~     strks(maxstr,2),
     ~     lenstr(maxstr),lplen,i,nhel,k,OPTION,istart,iend,
CPDBSUM RAL 12/4/05 -->
C     ~     BABHEL(20,2),MHEL,l,j
     ~     BABHEL(MAXHEL,2),MHEL,l,j
CPDBSUM RAL 12/4/05 <--
      character bseqno(maxres)*6,constr(maxstr),sstrks(maxres),
     ~     coption,pdbcod*4
      character*(fnamln) pdbnam,babnam,psinam,pcod
CPDBSUM RAL 11/4/05 -->
       CHARACTER*1   aa(maxres)
       CHARACTER*3   RSNAM1, RSNAM2, RSNAM3, RSNAM4
       INTEGER       LENST, outbab, outpsi
CPDBSUM RAL 11/4/05 <--

      bseqno(maxres)='000000'
      nbab=0
      npsi=0
      i=index(pdbnam,'.')
      if(i.gt.4) then
         pdbcod=pdbnam(i-4:i-1)
      else
         if(i.ne.0) then
            pdbcod=pdbnam(1:4)
         else
            j=index(pdbnam,' ')
            pdbcod=pdbnam(j-4:j-1)
         endif
      endif
      pcod=pdbnam(istart:iend)
C     SEARCH FOR BETA-ALPHA-BETA UNITS
c      write(9,*)
c      write(9,*)'BETA ALPHA BETA UNITS:'
      do 30 i=1,nstr
         if((icnstr(i).eq.1.or.icnstr(i).eq.-1).
     ~        and.constr(i).eq.'X') then
            i1=strks(i,1)
            i2=strks(i,2)
            j1=strks(i+1,1)
            j2=strks(i+1,2)
            nhel=0
C     ZERO ARRAYS FOR STORING HELICES WITHIN BETA ALPHA BETA UNITS
            mhel=0
CPDBSUM RAL 12/4/05 -->
C            DO 22 k=1,20
            DO 22 k=1,MAXHEL
CPDBSUM RAL 12/4/05 <--
               DO 21 L=1,2
CPDBSUM RAL 12/4/05 -->
                  BABHEL(K,L)=maxres
CPDBSUM RAL 12/4/05 <--
 21            CONTINUE
 22         CONTINUE
            lplen=j1-i2-1
            do 20 k=i2+1,j1-1
               if(sstrks(k).eq.'H'.or.sstrks(k).eq.'h') nhel=nhel+1
c     DETECT BEGINNING AND END OF HELICES IN BETA ALPHA BETA UNITS
c     STORE IN ARRAY BABHEL; NUMBER OF HELICES IS MHEL
               if(sstrks(k).eq.'H') then
                  if(sstrks(k-1).ne.'H') then
                     mhel=mhel+1
CPDBSUM RAL 12/4/05 -->
                     if (mhel.le.MAXHEL) then
CPDBSUM RAL 12/4/05 <--
                     babhel(mhel,1)=k
CPDBSUM RAL 12/4/05 -->
                     endif
CPDBSUM RAL 12/4/05 <--
                  endif
                  if(sstrks(k+1).ne.'H') then
CPDBSUM RAL 12/4/05 -->
                     if (mhel.le.MAXHEL) then
CPDBSUM RAL 12/4/05 <--
                     babhel(mhel,2)=k
CPDBSUM RAL 12/4/05 -->
                     endif
CPDBSUM RAL 12/4/05 <--
                  endif
               endif
 20         continue
            if(nhel.gt.0) then
               nbab=nbab+1
c               if(OPTION.eq.2) then
               if(coption.eq.'L'.or.coption.eq.'l') then
                   write(11,15) pdbcod,bseqno(i1),bseqno(i2),
     ~                 bseqno(j1),bseqno(j2),(BSEQNO(BABHEL(K,1)),
     ~                 BSEQNO(BABHEL(K,2)),K=1,5),lenstr(i),lenstr(i+1),
     ~                 lplen,nhel,mhel
               else
                  if(nbab.eq.1) then
                     nps=index(pcod,' ')
                     BABNAM=PCOD(1:NPS-1)//'.bab'
CPDBSUM RAL 14/4/05 -->
C                     OPEN(UNIT=11,FILE=BABNAM,ERR=1500)
                     OPEN(UNIT=11,status='scratch',ERR=1500)
CPDBSUM RAL 14/4/05 <--
                  endif
                  write(11,151) bseqno(i1),bseqno(i2),bseqno(j1),
     ~                 bseqno(j2),lenstr(i),lenstr(i+1),lplen,
     ~                 nhel,i1,j2,
     ~                 (bseqno(babhel(k,1)),bseqno(babhel(k,2)),k=1,5)
CPDBSUM RAL 11/4/05 -->
C----             Convert single-letter a.a. codes to 3-letter residue names
                  CALL AA2NAM(aa(i1),RSNAM1)
                  CALL AA2NAM(aa(i2),RSNAM2)
                  CALL AA2NAM(aa(j1),RSNAM3)
                  CALL AA2NAM(aa(j2),RSNAM4)

C----             Write out data to sheets.promotif data file
                  WRITE(outbab,640) bseqno(i1)(1:1),
     -               bseqno(i1)(2:6), RSNAM1, bseqno(i2)(2:6), RSNAM2,
     -               bseqno(j1)(2:6), RSNAM3, bseqno(j2)(2:6), RSNAM4,
     -               lenstr(i), lenstr(i+1), lplen, mhel, nhel,
     -               (bseqno(babhel(k,1))(2:6),
     -                bseqno(babhel(k,2))(2:6), k=1,5)
 640              FORMAT(A1,'::',4(A5,'::',A3,'::'),5(I3,'::'),
     -               10(A5,'::'))
CPDBSUM RAL 11/4/05 <--
               endif
            endif
 15         format(a4,1x,a6,1x,a6,3x,a6,1x,a6,1x,10(1x,a6),1x,4(1x,i3),1x,
     ~           i2)
 151        format(a6,1x,a6,3x,a6,1x,a6,4(1x,i3),1x,i4,1x,i4,10(1x,a6))
         endif
 30   continue
c     SEARCH FOR PSI-LOOPS
c      write(9,*)
c      write(9,*)'PSI LOOPS: '
      do 40 i=1,nstr
         if((icnstr(i).eq.2.or.icnstr(i).eq.-2).
     ~        and.constr(i).eq.' ') then
            i1=strks(i,1)
            i2=strks(i,2)
            j1=strks(i+1,1)
            j2=strks(i+1,2)
            lplen=j1-i2-1
            npsi=npsi+1
c            if(OPTION.eq.2) then
            if(coption.eq.'L'.or.coption.eq.'l') then
               write(12,35) pdbcod,bseqno(i1),bseqno(i2),
     ~              bseqno(j1),bseqno(j2),lenstr(i),lenstr(i+1),
     ~              lplen
            else
               if(npsi.eq.1) then
                  nps=index(pcod,' ')
                  PSINAM=PCOD(1:NPS-1)//'.psi'
CPDBSUM RAL 14/4/05 -->
C                  OPEN(UNIT=12,FILE=PSINAM,ERR=1500)
                  OPEN(UNIT=12,status='scratch',ERR=1500)
CPDBSUM RAL 14/4/05 <--
               endif
               write(12,351) bseqno(i1),bseqno(i2),bseqno(j1),
     ~              bseqno(j2),i1,j2,lenstr(i),lenstr(i+1),lplen
CPDBSUM RAL 12/4/05 -->
C----          Convert single-letter a.a. codes to 3-letter residue names
               CALL AA2NAM(aa(i1),RSNAM1)
               CALL AA2NAM(aa(i2),RSNAM2)
               CALL AA2NAM(aa(j1),RSNAM3)
               CALL AA2NAM(aa(j2),RSNAM4)

C----          Write out data to sheets.promotif data file
               WRITE(outpsi,660) bseqno(i1)(1:1),
     -             bseqno(i1)(2:6), RSNAM1, bseqno(i2)(2:6), RSNAM2,
     -             bseqno(j1)(2:6), RSNAM3, bseqno(j2)(2:6), RSNAM4,
     -             lenstr(i), lenstr(i+1), lplen
 660           FORMAT(A1,'::',4(A5,'::',A3,'::'),3(I3,'::'))
CPDBSUM RAL 12/4/05 <--
            endif
 35         format(a4,1x,a6,1x,a6,3x,a6,1x,a6,3(1x,i3))
 351        format(a6,1x,a6,3x,a6,1x,a6,1x,i4,1x,i4,3(1x,i3))
         endif
 40   continue

 1500 return
      end
c*******************************************************
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
c***************************************************************************
      subroutine hbext(hbond,hbe,numres,hbp)
      include 'p_sheet.par'
      
      integer ihb,i,j,hbond(4*maxres,2),np,m,n,hbp(maxres,4),
     ~     numres
      real hbe(maxres,4),energy
      
      data energy/-0.5/
      ihb=0
      do 10 i=1,2
         do 5 j=1,4*maxres
            hbond(j,i)=0
 5       continue
 10   continue
      do 30 m=1,numres
         do 20 n=1,4
            if(hbe(m,n).le.energy) then
               np=hbp(m,n)
               if((np.lt.m.and.np.ge.1).or.np.eq.0) go to 20
               if(hbp(m,n).gt.0) then
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
 20    continue
 30   continue
      return
      end
c**************************************************************
CPDBSUM RAL 15/4/05 -->
C**************************************************************************
C
C  SUBROUTINE SUMMAR  -  Write summary information on numbers and
C                        of strands and helices in each chain
C
C----------------------------------------------------------------------+---

      SUBROUTINE SUMMAR(IUNIT,SEQNUM,SST,MAXRES,NRES)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 9)

      INTEGER       MAXRES

      CHARACTER*1   CHAIN, LSTCHN, SST(MAXRES)
      CHARACTER*3   IDTAG
      CHARACTER*6   SEQNUM(MAXRES)
      CHARACTER*20  FLDNAM(NFIELD)
      CHARACTER*512 OUTNAM

      INTEGER       IFIELD, IRES, IUNIT, LENST, NHEL, NGEL, NSTR, NRES,
     -              NREST, NTOTAL

      DATA FLDNAM / 'CHAIN', 'NRES_HELIX', 'NRES_3_10', 'NRES_STRAND',
     -     'NRES_REST', 'PERCEN_HELIX', 'PERCEN_3_10', 'PERCEN_STRAND',
     -     'PERCEN_REST' /

C---- Initialise variables
      LSTCHN = SEQNUM(1)(1:1)
      NHEL = 0
      NGEL = 0
      NSTR = 0
      NREST = 0
      NTOTAL = 0

C---- Open summary file
      outnam = "summary.promotif"
      OPEN(UNIT=IUNIT,FILE=OUTNAM,STATUS='UNKNOWN',ERR=900)

C---- Define field id tag
      IDTAG = 'SU_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:SUMMARY'
 10   FORMAT(500A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENST(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENST(FLDNAM(NFIELD)))

C---- Loop through all the residues counting the secondary structure
C     types
      DO 200, IRES = 1, NRES

C----     Get the chain identifier
          CHAIN = SEQNUM(IRES)(1:1)

C----     If this is a new chain, write out old chain's statistics
C         and re-initialise
          IF (CHAIN.NE.LSTCHN) THEN

C----         If have data for a previous chain, write out
              IF (NTOTAL.GT.0) THEN
                  CALL WRISUM(IUNIT,LSTCHN,NHEL,NGEL,NSTR,NREST,NTOTAL)
              ENDIF

C----         Save the chain identifier
              LSTCHN = CHAIN
          ENDIF

C----     Increment appropriate count
          IF (SST(IRES).EQ.'H') THEN
              NHEL = NHEL + 1
          ELSE IF (SST(IRES).EQ.'G') THEN
              NGEL = NGEL + 1
          ELSE IF (SST(IRES).EQ.'E') THEN
              NSTR = NSTR + 1
          ELSE
              NREST = NREST + 1
          ENDIF
          NTOTAL = NTOTAL + 1

 200  CONTINUE

C---- Write out stats for last chain
      IF (NTOTAL.GT.0) THEN
          CALL WRISUM(IUNIT,LSTCHN,NHEL,NGEL,NSTR,NREST,NTOTAL)
      ENDIF

      GO TO 999

 900  CONTINUE
      PRINT*, 'ERROR. Unable to open output file: ',
     -     OUTNAM(1:LENST(OUTNAM))

 999  CONTINUE
      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE WRIMAR  -  Calculate percentages and write out
C
C----------------------------------------------------------------------+---

      SUBROUTINE WRISUM(IUNIT,CHAIN,NHEL,NGEL,NSTR,NREST,NTOTAL)

      CHARACTER*1   CHAIN

      INTEGER       IUNIT, NHEL, NGEL, NREST, NSTR, NTOTAL
      REAL          PEREST, PERHEL, PERGEL, PERSTR

C---- Calculate the percentages
      IF (NTOTAL.GT.0) THEN
          PERHEL = 100.0 * REAL(NHEL) / REAL(NTOTAL)
          PERGEL = 100.0 * REAL(NGEL) / REAL(NTOTAL)
          PERSTR = 100.0 * REAL(NSTR) / REAL(NTOTAL)
          PEREST = 100.0 * REAL(NREST) / REAL(NTOTAL)
      ELSE
          PERHEL = 0.0
          PERGEL = 0.0
          PERSTR = 0.0
          PEREST = 0.0
      ENDIF

C---- Write the data out to the summary.promotif file
      WRITE(IUNIT,20) CHAIN, NHEL, NGEL, NSTR, NREST, PERHEL, PERGEL,
     -    PERSTR, PEREST
 20   FORMAT(A1,'::',4(I6,'::'),4(F5.1,'::'))

C---- Re-initialise counts
      NHEL = 0
      NGEL = 0
      NREST = 0
      NSTR = 0
      NTOTAL = 0

      RETURN
      END

C--------------------------------------------------------------------------
CPDBSUM RAL 15/4/05 <--
CPDBSUM RAL 30/3/05 -->
C**************************************************************************
C
C  SUBROUTINE SHTHED  -  Write header records to sheets.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE SHTHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 6)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENST, IUNIT

      DATA FLDNAM / 'CHAIN', 'SHEET_ID', 'NSTRANDS', 'TYPE',
     -     'BARREL_FLAG', 'TOPOLOGY' /

C---- Define field id tag
      IDTAG = 'SH_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:SHEETS'
 10   FORMAT(500A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENST(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENST(FLDNAM(NFIELD)))

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE STRHED  -  Write header records to strands.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE STRHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 25)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENST, IUNIT

      DATA FLDNAM / 'CHAIN', 'STRAND_NO',
     -              'RES_NUM1', 'RES_NAM1',
     -              'RES_NUM2', 'RES_NAM2',
     -              'SHEET_ID1', 'SHEET_ID2', 'NRESID',
     -              'SEQUENCE', 'NPARTNERS',
     -              'ABPART1', 'ABPART2', 'ABPART3', 'ABPART4',
     -              'ABSDIR1', 'ABSDIR2', 'ABSDIR3', 'ABSDIR4',
     -              'ICNSTR', 'CONSTR', 'ICSTR1', 'CONSTR1',
     -              'SECSTR_LEFT', 'SECSTR_RIGHT' /

C---- Define field id tag
      IDTAG = 'ST_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:STRANDS'
 10   FORMAT(500A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENST(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENST(FLDNAM(NFIELD)))

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE BABHED  -  Write header records to bab.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE BABHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 24)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENST, IUNIT

      DATA FLDNAM / 'CHAIN',
     -              'RES_NUM1', 'RES_NAM1', 'RES_NUM2', 'RES_NAM2',
     -              'RES_NUM3', 'RES_NAM3', 'RES_NUM4', 'RES_NAM4',
     -              'LEN_STRAND1', 'LEN_STRAND2',
     -              'LEN_LOOP', 'NHELICES', 'LEN_HELIX',
     -              'START_HELIX1', 'END_HELIX1',
     -              'START_HELIX2', 'END_HELIX2',
     -              'START_HELIX3', 'END_HELIX3',
     -              'START_HELIX4', 'END_HELIX4',
     -              'START_HELIX5', 'END_HELIX5' /

C---- Define field id tag
      IDTAG = 'BA_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:BAB_MOTIFS'
 10   FORMAT(500A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENST(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENST(FLDNAM(NFIELD)))

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE PSIHED  -  Write header records to psi.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE PSIHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 12)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENST, IUNIT

      DATA FLDNAM / 'CHAIN',
     -              'RES_NUM1', 'RES_NAM1', 'RES_NUM2', 'RES_NAM2',
     -              'RES_NUM3', 'RES_NAM3', 'RES_NUM4', 'RES_NAM4',
     -              'LEN_STRAND1', 'LEN_STRAND2',
     -              'LEN_LOOP' /

C---- Define field id tag
      IDTAG = 'PS_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:PSI_LOOPS'
 10   FORMAT(500A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENST(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENST(FLDNAM(NFIELD)))

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
C  FUNCTION LENST  -  Return length of string 'CHARS' excluding
C                      trailing blanks
C
C----------------------------------------------------------------------+---

      INTEGER FUNCTION LENST(CHARS)

      CHARACTER*(*) CHARS
      INTEGER J, MVAL

      MVAL = LEN(CHARS)
      DO 10, J = MVAL, 1, -1
          IF (CHARS(J:J).NE.' ') GO TO 20
   10 CONTINUE
      J = 0
   20 CONTINUE
      LENST = J

      RETURN
      END

C--------------------------------------------------------------------------
CPDBSUM RAL 30/3/05 <--
