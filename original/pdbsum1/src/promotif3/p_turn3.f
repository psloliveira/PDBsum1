      PROGRAM turn
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: p_turn_ensem.f                                           *
C*  FUNCTION: Program to identify and classify turns in proteins       * 
C* COPYRIGHT: (1989) David Keith Smith                                 *
C*---------------------------------------------------------------------*
C*    AUTHOR: D. K. SMITH and GAIL HUTCHINSON                          *
C*      DATE: 1989                                                         *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
c                    now multiple protein version for promotif
c  03/05/95  GH      separate out beta and gamma turn output
c  03/05/95  GH      include option to either print out data for each 
c                    protein separately or put all data in one file
c  27/03/96  GH      reorganise multiple/single file handling
c  29/04/96  GH      this version created to allow command line options
c  24/11/97  GH      altered to allow '-' to be a chain letter (eg in 1gav)
c  01/11/00  GH      removed 'x' specifier in format statements
c                    changed ierror from integer to logical
c  27/11/00  GH      modified to read in new sstruc format
c***********************************************************************
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -o p_turn3 p_turn3.f
C
c***********************************************************************
c
c     global variables
c
      include 'p_turn.par'
c
      integer   turnst(maxaas), hbondp(2,maxaas), nturns, seqlen,
     +          turnty(3,maxaas),ppset(nphis,nphis),carrie(2,maxaas),
     +          namlen,i,iseqno(maxaas),option
      real      phi(maxaas), psi(maxaas), cadp3(maxaas), cadp2(maxaas),
     +          hbonde(2,maxaas), turnen(maxaas), omega(maxaas),
     +          chi1(maxaas)
      character aaids(maxaas)*6, brcode*4, names1(maxaas)*1,
     +          trnpat(3,maxaas)*1, misflg(2,maxaas)*1, strk(maxaas)*1,
     +          hbpres(maxaas)*1,austrk(maxaas)*1,coption
      logical   ppok,ENSEMB
c
      integer   ok, finish, endfi, opnera, opnerb, opnerc, endnfi,
     ~     nps,istart,iend
c
      parameter (ok = 0, finish = 2, endfi = 4, 
     +           opnera = 1, opnerb = 3, opnerc = 5, endnfi = 6)
c
      character*(fnamln) sstnam, ppfil, pdbnam,pdbcod,finam,outnam,
     ~           outnmb,outnmg
      integer   fistat, iocode,nxpsi
      integer   strlen
CPDBSUM RAL 14/3/05 -->
      integer   lenstr
CPDBSUM RAL 14/3/05 <--
      logical   finuse,ierror

      DATA OPTION/1/,FISTAT/OK/,FINUSE/.FALSE./

      READ(*,'(A)')PDBNAM
      IF(PDBNAM(1:1).EQ.'%') THEN
         READ(*,'(A1)') COPTION
         ENSEMB=.TRUE.
         IF(COPTION.EQ.'L'.OR.COPTION.EQ.'l') THEN
            OUTNMB='BETATURNS'
            OUTNMG='GAMMATURNS'
            OPEN(unit=outfi, file=outnmb, status='unknown', err=902)
            OPEN(unit=outfig, file=outnmg, status='unknown', err=902)
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
            READ(NAMFI,'(A)',IOSTAT=IOCODE,END=999) PDBNAM
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
         NPS=INDEX(pdbcod,' ')
         SSTNAM=PDBCOD(1:nps-1)//'.sst'
         OPEN(UNIT=STRKFI,FILE=SSTNAM,STATUS='OLD',IOSTAT=IOCODE)
         IF(IOCODE.NE.0) FISTAT=OPNERB
      ENDIF

c     if survived this far have a file to process otherwise
c     provide error messages

      IF(FISTAT.EQ.OK) THEN
c     change name of o/p file 13.2.95
c         IF(OPTION.NE.2) THEN
         if(coption.ne.'L'.and.coption.ne.'l') then
            OUTNMB=PDBCOD(1:NPS-1)//'.bturns'
            OUTNMG=PDBCOD(1:NPS-1)//'.gturns'
         ENDIF
c     end of change
         
         ppfil = 'phipsi.mat'
         open(unit=ppfi, file=ppfil, status='old', iostat=iocode)
         if (iocode .ne. 0) then
            ppok = .false.
            write(*,7040) ppfil(1:strlen(ppfil))
 7040 format ('PPset file "',A,'" not found -',/,10X,
     +        'Wilmot classification will not be done')
         else
            ppok = .true.
            do 300, nxpsi = 1, nphis
            read(ppfi, 250, err=903, end=904) (ppset(nxpsi,i),i=1,nphis)
  250       format(36(1x,i1))
  300       continue
            close(unit=ppfi)
         end if
c
         call rdstrk(brcode, cadp3, cadp2, phi, psi, aaids, hbondp, 
     +               omega, trnpat, names1, seqlen, misflg, strk, 
     +               hbonde,austrk,chi1,iseqno)
         call gentrn(cadp3, phi, psi, hbondp, trnpat, nturns, 
     +               turnty, turnst, seqlen, misflg, strk, hbonde, 
     +               turnen, ppset, carrie, hbpres)
         if (nturns .gt. 0) then
c            if(option.ne.2) then
            if(coption.ne.'L'.and.coption.ne.'l') then
CPDBSUM RAL 14/4/05 -->
C               open(unit=outfi, file=outnmb, status='unknown', err=902)
C               open(unit=outfig, file=outnmg, status='unknown', err=902)
               open(unit=outfi, status='scratch', err=902)
               open(unit=outfig, status='scratch', err=902)
CPDBSUM RAL 14/4/05 <--
               i=index(pdbnam,'.')
c               if(i.gt.4) then
c                  write(outfi,123) pdbnam(i-4:i-1)
c                  write(outfig,123) pdbnam(i-4:i-1)
c               else
c                  write(outfi,123) pdbnam(1:4)
c                  write(outfig,123) pdbnam(1:4)
c               endif
                  write(outfi,123) pdbcod(1:nps-1)
                  write(outfig,123) pdbcod(1:nps-1)
 123           format('Promotif output for protein :',a)
            endif
            call result(nturns, turnty, turnst, phi, psi, cadp3,
     +                  cadp2, ppset, aaids, names1, carrie,
     +                  hbpres,strk,austrk,omega,chi1,iseqno,
     +                  coption,pdbnam)
c            if(option.ne.2) then
            if(coption.ne.'L'.and.coption.ne.'l') then
               close(unit=outfi)
               close(unit=outfig)
            endif
         endif
      ELSE IF(FISTAT.EQ.OPNERB) THEN
CPDBSUM RAL 15/3/05 -->
C         write(*,*)'Unable to open the 2D structure file ',
C     +              sstnam(1:namlen)
         write(*,*)'*** ERROR. Unable to open the 2D structure file ',
     +              sstnam(1:lenstr(sstnam))
CPDBSUM RAL 15/3/05 <--
      ENDIF
c
c
      IF(ENSEMB) THEN
         IF(FISTAT.NE.FINISH) THEN
            FISTAT=OK
            IF(FINUSE) THEN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF
         goto 999
 7000 format (' Error opening output file "',A,'"')
CPDBSUM RAL 14/3/05 -->
C  902 write (*, 7000) outnam(1:5)
  902 write (*, 7000) outnam(1:lenstr(outnam))
CPDBSUM RAL 14/3/05 <--
      goto 999
  903 write (*, 7010) ppfil(1:strlen(ppfil))
 7010 format ('Error reading PP input file "',A,'"')
      goto 999
  904 write (*,7020) ppfil(1:strlen(ppfil))
 7020 format ('Unexpected EOF reading PP input file "',A,'"')
  999 continue
      END
c
c*****************************************************************************
c
      Subroutine rdstrk(brcode, cadp3, cadp2, phi, psi, aaids, hbondp, 
     +                  omega, trnpat, names1, seqlen, misflg, strk, 
     +                  hbonde,austrk,chi1,iseqno)
c
c     routine to read the secondary structure files and place the 
c     turn details in the appropriate arrays
c
      include 'p_turn.par'
c
      integer   hbondp(2,maxaas), seqlen,iseqno(maxaas)
      character brcode*4, aaids(maxaas)*6, trnpat(3,maxaas)*1, 
     +          names1(maxaas)*1, misflg(2,maxaas)*1, strk(maxaas)*1,
     +          austrk(maxaas)*1
      real      cadp3(maxaas), phi(maxaas), psi(maxaas), cadp2(maxaas),
     +          hbonde(2,maxaas), omega(MAXAAS),chi1(maxaas)
c
c     local variables
c
      integer   nxaa, i
      character chnbrk(maxaas)*1, strkln*133
c
c     read structure file to get the brcode and seqlen
c
      do 10, i = 1, brcdpl
      read(strkfi, '(a)', err=901, end=901) strkln
   10 continue
      read(strkln, brform) brcode
      read(strkfi, slform, err=902, end=902) seqlen
      do 20 i = 1, hdsz
         read(strkfi, *, err=903, end=903) 
   20 continue
      nxaa=0
c 100  format(1x,i4,a1,a6,7x,a1,1x,a1,3x,a1,1x,a1,3x,2a1,13x,
c     +     3(1x,f6.1),19x,2(11x,i4,1x,f4.1))
 100  format(i5,a1,a6,7x,a1,1x,a1,3x,a1,1x,a1,5x,2a1,13x,
     +     3(1x,f6.1),19x,2(10x,i5,1x,f4.1))
      do 201 nxaa=1,seqlen
         read(strkfi, 100, err=904, end=904) iseqno(nxaa),chnbrk(nxaa), 
     +        aaids(nxaa), 
     +        names1(nxaa), austrk(nxaa),strk(nxaa), 
     +        (trnpat(i,nxaa), i = 1, 3), 
     +        phi(nxaa), psi(nxaa), omega(nxaa), (hbondp(i,nxaa), 
     +        hbonde(i,nxaa), i = 1, 2)
 201  continue
      do 320 i = 1, intset
         read(strkfi, *, err=905, end=905)
 320  continue
c
c     now at place to get cadp3
c
      nxaa=0
         do 1000, nxaa = 1, seqlen
            read(strkfi, 500, err=906, end=906) chi1(nxaa),cadp2(nxaa), 
     +           cadp3(nxaa)
 500        format(29x,f6.1,49x,f4.1,1x,f4.1)
 1000    continue
      close(unit=strkfi)
      do 2000, nxaa = 1, seqlen - 1
      misflg(beta,nxaa) = space
      misflg(gamma,nxaa) = space
      if (chnbrk(nxaa) .eq. brksym .and.
     +   chnbrk(nxaa+1) .eq. brksym) then
         do 1500, i = nxaa, nxaa - btrnsz + 1, -1
         if (i .gt. 0) misflg(beta,i) = misch
 1500    continue
         do 1600, i = nxaa, nxaa - gtrnsz + 1, -1
         if (i .gt. 0) misflg(gamma,i) = misch
 1600    continue
      end if
 2000 continue
      misflg(beta,seqlen) = space
      misflg(gamma,seqlen) = space
      return
c
  901 print*,'Error 1'
      stop
  902 print*,'Error 2'
      stop
  903 print*,'Error 3'
      stop
  904 print*,'Error 4'
      stop
  905 print*,'Error 5'
      stop
  906 print*,'Error 6'
      stop
      end
c
c****************************************************************************
c
      subroutine gentrn(cadp3, phi, psi, hbondp, trnpat, nturns, 
     +                 turnty, turnst, seqlen, misflg, strk, hbonde, 
     +                  turnen, ppset, carrie, hbpres)
c

c     extract the various turn types
c
      include 'p_turn.par'
c
      integer   hbondp(2,maxaas), nturns, turnst(maxaas), seqlen,
     +          turnty(3,maxaas),ppset(nphis,nphis),carrie(2,maxaas)
      real      cadp3(maxaas), phi(maxaas), psi(maxaas), 
     +          hbonde(2,maxaas), turnen(maxaas)
      character trnpat(3,maxaas)*1, misflg(2,maxaas)*1, strk(maxaas)*1,
     ~          hbpres(maxaas)*1

      integer   nxaa, nxturn, hlxsz(3), nxhlx, typenm, stntyp,
     +          libtyp, wiltyp, cartyp, trnsta, angchk, nlge, i 
      logical   helicl(maxaas), hlxsht(maxaas)
      real      phipsi(4,nbeta-1), bndeng, gampp(2,2), gameng, gamtol
      character bonded*1
      logical   xeqy
c
      data hlxsz / 4, 3, 5/,gameng/-1.1/,gamtol/40.0/
      data phipsi / -60.0,  -30.0,  -90.0,   0.0,
     +               60.0,   30.0,   90.0,   0.0,
     +              -60.0,  120.0,   80.0,   0.0,
     +               60.0, -120.0,  -80.0,   0.0,
     +              -60.0,  120.0,  -90.0,   0.0,
     +             -120.0,  120.0,  -60.0,   0.0,
     +              -60.0,  -30.0, -120.0,  120.0,
     +             -135.0,  135.0,   -75.0,   160.0/
      data (gampp(1,i),i=1,2) / 75.0, -64.0/,
     +     (gampp(2,i),i=1,2) /-79.0, 69.0/
      save hlxsz, phipsi, gampp

      nturns = 0
      do 100, nxaa = 1, seqlen
      helicl(nxaa) = .false.
  100 continue

      do 600, nxturn = 1, hlxtyp
      do 500, nxaa = 2, seqlen - hlxsz(nxturn)
      if ((trnpat(nxturn,nxaa) .eq. bndbeg .or.
     +   trnpat(nxturn,nxaa) .eq. bndbth) .and.
     +   (trnpat(nxturn,nxaa-1) .eq. bndbeg .or.
     +   trnpat(nxturn,nxaa-1) .eq. bndbth)) then
         do 400, nxhlx = 0, hlxsz(nxturn) - 1
         helicl(nxaa+nxhlx) = .true.
  400    continue
      end if
  500 continue
  600 continue
      do 800, nxaa = 1, seqlen
      if (strk(nxaa).eq.'E'.or.helicl(nxaa)) then
         hlxsht(nxaa) = .true.
      else
         hlxsht(nxaa) = .false.
      end if
  800 continue
      do 2000, nxaa = 1, seqlen - btrnsz
      if (misflg(beta,nxaa) .eq. space) then
         if (.not. helicl(nxaa+1) .and. .not. helicl(nxaa+2)) then
            stntyp = 0
            libtyp = 0
            if (cadp3(nxaa).gt.0.0.and.cadp3(nxaa).lt.btrnds) then
               typenm = 1
               trnsta = goodtn
 1000          continue
               nlge=0
               trnsta = angchk(phi(nxaa+1),phipsi(1,typenm))
               if (trnsta .ne. abandn) then
                  if(trnsta.eq.lgetn) nlge=nlge+1
                  trnsta = angchk(psi(nxaa+1),phipsi(2,typenm))
               endif
               if (trnsta .ne. abandn) then
                  if(trnsta.eq.lgetn) nlge=nlge+1
                  trnsta = angchk(phi(nxaa+2),phipsi(3,typenm))
               endif
               if (trnsta .ne. abandn) then
                  if(trnsta.eq.lgetn) nlge=nlge+1
                  trnsta = angchk(psi(nxaa+2),phipsi(4,typenm))
               endif
               if (trnsta.eq.abandn.and.typenm.lt.nbeta - 1) then
                  typenm = typenm + 1
                  go to 1000
               else
                  if(trnsta.ne.abandn) then
                     if(trnsta.eq.lgetn) nlge=nlge+1
                     if(nlge.le.1) then
                        libtyp=typenm
                     else
                        libtyp=nbeta
                     endif
                  endif
               endif
c              if (trnsta .eq. lgetn) then
c                 nlge=nlge+1
c                 stntyp = nbeta
c                 libtyp = typenm
c              end if
c              if (trnsta .eq. goodtn) then
c                 stntyp = typenm
c                 libtyp = typenm
c              end if
c           end if
               if (trnsta .eq. abandn) then
                  stntyp = nbeta
                  libtyp = nbeta
               end if
               cartyp = wiltyp(carrie(1,nxaa),carrie(2,nxaa),
     +                  ppset, PHI(NXAA+1),
     +                  phi(nxaa+2),PSI(NXAA+1),psi(nxaa+2))
            end if
            bonded = 'N'
            bndeng = 0.0
            if (trnpat(ks3,nxaa) .eq. bndbeg .or.
     +         trnpat(ks3,nxaa) .eq. bndbth) then
               bonded = 'Y'
               if (hbondp(1,nxaa) .EQ. nxaa + btrnsz) then
                  bndeng = hbonde(1,nxaa)
               else
                  bndeng = hbonde(2,nxaa)
               end if
            end if
c           if (bonded .eq. 'Y' .and. stntyp .eq. 0) 
c    +         stntyp = mxbgdf + ks3
            if (bonded.eq.'Y'.and.libtyp.eq.0) libtyp=mxbgdf+ks3
c           if (stntyp .ne. 0) then
            if (libtyp.ne.0) then
               nturns = nturns + 1
               turnst(nturns) = nxaa
               turnty(1,nturns) = stntyp
               turnty(2,nturns) = libtyp
               turnty(3,nturns) = cartyp
               turnen(nturns) = bndeng
c     RECORD WHETHER HYDROGEN BONDED TURN OR NOT 8.8.94
               hbpres(nturns)=bonded
C     END TEST
            end if
         end if
      end if
 2000 continue
      do 4000, nxturn = ks4, ks5, ks5 - ks4
      do 3900, nxaa = 1, seqlen - hlxsz(nxturn)
      if (trnpat(nxturn,nxaa) .eq. bndbeg .or.
     +   trnpat(nxturn,nxaa) .eq. bndbth) then
         if ((.not. helicl(nxaa+1) .and. .not. helicl(nxaa+2)) .or.
     +      (.not. helicl(nxaa+hlxsz(nxturn)-1) .and.
     +      .not. helicl(nxaa+hlxsz(nxturn)-2))) then
            nturns = nturns + 1
            turnst(nturns) = nxaa
c           turnty(1,nturns) = mxbgdf + nxturn
            turnty(2,nturns)=mxbgdf+nxturn
            if (hbondp(1,nxaa) .eq. nxaa + hlxsz(nxturn)) then
               turnen(nturns) = hbonde(1,nxaa)
            else
               turnen(nturns) = hbonde(2,nxaa)
            end if
         end if
      end if
 3900 continue
 4000 continue

      do 6000, nxaa = 1, seqlen - gtrnsz
      if (.not. hlxsht(nxaa+1)) then
         if ((hbondp(1,nxaa) .eq. nxaa + gtrnsz .or.
     +      hbondp(2,nxaa) .eq. nxaa + gtrnsz ).and.
     +      ((abs(abs(phi(nxaa+1))-abs(gampp(1,1))) .lt. gamtol .and.
     +      abs(abs(psi(nxaa+1))-abs(gampp(1,2))) .lt. gamtol) .or.
     +      (abs(abs(phi(nxaa+1))-abs(gampp(2,1))) .lt. gamtol .and.
     +      abs(abs(psi(nxaa+1))-abs(gampp(2,2)))  .lt.gamtol) .or.
     +      xeqy(phi(nxaa+1), nulang) .or. xeqy(psi(nxaa+1),nulang))) 
     +      then
            if (misflg(gamma,nxaa) .eq. space) then
               if (hbondp(1,nxaa) .eq. nxaa + gtrnsz) then
                  bndeng = hbonde(1,nxaa)
               else
                  bndeng = hbonde(2,nxaa)
               end if
               if (bndeng .lt. gameng) then
                  nturns = nturns + 1
                  turnst(nturns) = nxaa
                  turnen(nturns) = bndeng
                  if (xeqy(phi(nxaa+1), nulang) .or. 
     +               xeqy(psi(nxaa+1), nulang)) then
c                    turnty(1,nturns) = gamodd
                     turnty(2,nturns)=gamodd
                  else
                     if (phi(nxaa+1) .gt. 0.0) then
                        turnty(2,nturns)=gamcla
c                       turnty(1,nturns) = gamcla
                     else
                        turnty(2,nturns)=gaminv
c                       turnty(1,nturns) = gaminv
                     end if
                  end if
               end if
            end if
         end if
      end if
 6000 continue

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function angchk (chkang, basang)
c
c     check if the angle is within the tolerance
c
      include 'p_turn.par'
c
      real chkang, basang

      If (abs(chkang-basang) .lt. angtol) then
         angchk = goodtn
      else
         If (abs(chkang-basang) .lt. lgetol) then
            angchk = lgetn
         else
            angchk = abandn
         end if
      end if
      if(angchk.eq.abandn.and.chkang.lt.0) then
         chkang=chkang+360.0
         if(abs(chkang-basang).lt.angtol) then
            angchk=goodtn
         else
            if(abs(chkang-basang).lt.lgetol) then
               angchk=lgetn
            else
               angchk=abandn
            endif
         endif
         chkang=chkang-360.0
      endif
      END
c
c******************************************************************************
c

      integer function wiltyp (carrie1, carrie2, ppset, phi1, phi2,
     +                        psi1, psi2)
c
c     routine to classify turns according to the wilmot definitions
c     this function is here only for historycal reasons.
 
      include 'p_turn.par'
 
      integer   ppset(nphis,nphis), carrie1, carrie2
      real      phi1, phi2, psi1, psi2

      integer   phipl, psipl, nxcode, angpl, ccode(4,nbeta-1), i
c      character wilcod(nbeta-1)*2
c      data wilcod/'AA', 'LG', 'PG', 'EA', 'PA', 'BA', 'AB', 'BP'/
      data (ccode(1,i),i=1,nbeta-1)/3,5,2,7,2,1,3,1/
      data (ccode(2,i),i=1,nbeta-1)/4,5,2,7,2,1,4,1/
      data (ccode(3,i),i=1,nbeta-1)/3,6,6,3,3,3,1,2/
      data (ccode(4,i),i=1,nbeta-1)/4,6,6,4,4,4,1,2/

      phipl   = angpl(phi1)
      psipl   = nphis + 1 - angpl(psi1)
      carrie1 = ppset(psipl,phipl)
      phipl   = angpl(phi2)
      psipl   = nphis + 1 - angpl(psi2)
      carrie2 = ppset(psipl,phipl)
      nxcode  = 1
  500 continue
c      if (carrie1.ne. wilcod(nxcode)) then
c        nxcode = nxcode + 1
c        if (nxcode .lt. nbeta) goto 500
c      end if
      if ((carrie1.ne.ccode(1,nxcode).and.carrie1.ne.ccode(2,nxcode))
     +   .or.(carrie2.ne.ccode(3,nxcode).and.carrie2.ne.
     +   ccode(4,nxcode)))then
         nxcode=nxcode+1
         if(nxcode.lt.nbeta) go to 500
      endif
      wiltyp = nxcode
      END
c
c******************************************************************************
c

      integer function angpl (tstang)
c
      include 'p_turn.par'
c
c     get position in ppset for angle
c
      real     tstang
c
      angpl = int((tstang + 180.0) / 10.0) + 1
      if (angpl .gt. nphis) angpl = nphis
      END
c
c******************************************************************************

      subroutine result (nturns, turnty, turnst, phi, psi, cadp3, cadp2,
     +                   ppset, aaids, names1, carrie, hbpres,strk,
     +                   austrk,omega,chi1,iseqno,coption,pdbnam)
c
c     generate the fields for each turn
c
      include 'p_turn.par'
c
      integer   nbt, ngt
      integer   nturns, turnty(3,maxaas), turnst(maxaas), 
     +          carrie(2,maxaas), ppset(nphis,nphis),iseqno(maxaas),
     +          option,j
      real      phi(maxaas), psi(maxaas), cadp3(maxaas), cadp2(maxaas),
     +          omega(maxaas),chi1(maxaas)
      character aaids(maxaas)*6, names1(maxaas)*1, hbpres(maxaas)*1,
     +          austrk(maxaas)*1,om*4,strk(maxaas)*1,pdbcod*4,coption
      character*(fnamln) pdbnam
      integer   nxturn, nxaa, i,ist,k
      character wknull*6, type(mxbgdf)*8,
     +          efimov*2, eclass(7), chnid*1, efivom*1, getefi*1
CPDBSUM RAL 14/3/05 -->
      CHARACTER*3 RSNAM1, RSNAM2, RSNAM3, RSNAM4
      CHARACTER*5 HBFLAG
      CHARACTER*(fnamln) OUTNAM

      INTEGER   LENSTR

      LOGICAL   HAVBET, HAVGAM


      data type /'I', 'I''', 'II', 'II''', 'VIa1', 'VIa2', 'VIII', 
     +           'VIb',
     +           'IV',' ODD', ' CLASSIC', ' INVERSE'/
      data eclass/'B','P','A','a','L','G','E'/
      save type

CPDBSUM RAL 25/4/05 -->
      HAVBET = .FALSE.
      HAVGAM = .FALSE.
CPDBSUM RAL 25/4/05 <--
      wknull = idnull
      nbt    = 0
      ngt    = 0
      j=index(pdbnam,'.')
      if(j.gt.4) then
         pdbcod=pdbnam(j-4:j-1)
      else
         if(j.ne.0) then
            pdbcod=pdbnam(1:4)
         else
            k=index(pdbnam,' ')
            pdbcod=pdbnam(k-4:k-1)
         endif
      endif
c      if(option.ne.2) write(outfi,*)'Beta Turns'
      if(coption.ne.'L'.and.coption.ne.'l') write(outfi,*)'Beta Turns'
      do 2000, nxturn = 1, nturns
      nxaa = turnst(nxturn)
c     BETA TURN
c     if (turnty(1,nxturn) .le. nbeta) then
      if (turnty(2,nxturn).le.nbeta) then
         do 1005 i=1,2
         efimov(i:i)=' '
         if (carrie(i,nxaa).ne.0) efimov(i:i)=
     +      eclass(carrie(i,nxaa))
 1005    continue
         nbt = nbt + 1
         chnid = aaids(nxaa)(1:1)
c         if (chnid.eq.'-')  then
c            chnid = ' '
c         end if 
c     CHECK OMEGA ANGLES - ESPECIALLY FOR CIS PROLINES
         DO 1111 I=1,4
            OM(I:I)='T'
 1111    CONTINUE
         IF(NXAA.GT.1) THEN
            IST=NXAA-1
         ELSE
            IST=NXAA
         ENDIF
         DO 1112 I=IST,NXAA+2
            IF(OMEGA(I).GT.-90.0.AND.OMEGA(I).LT.90.0)
     ~           OM(I-NXAA+2:I-NXAA+2)='C'
 1112    CONTINUE
C     17.11.94 IMPOSE CIS-PROLINE IF TYPE VI TURNS
         IF(TURNTY(2,NXTURN).EQ.5.OR.TURNTY(2,NXTURN).EQ.6.
     ~        OR.TURNTY(2,NXTURN).EQ.8) THEN
            IF(OM(3:3).NE.'C'.OR.NAMES1(NXAA+2).NE.'P') THEN
               TURNTY(2,NXTURN)=9
            ENDIF
         ENDIF
         IF((AUSTRK(NXAA+1).EQ.'H'.AND.AUSTRK(NXAA+2).EQ.'H')
     +        .OR.(STRK(NXAA).EQ.'E'.AND.STRK(NXAA+1).EQ.'E'.
     +        AND.STRK(NXAA+2).EQ.'E'.AND.STRK(NXAA+3).EQ.'E')
     +        .OR.(AUSTRK(NXAA+1).EQ.'G'.AND.AUSTRK(NXAA+2).EQ.
     +        'G')) THEN
C     DISTORTED HELIX
         ELSE
c            if(option.eq.2) then
            if(coption.eq.'L'.or.coption.eq.'l') then
               write(outfi,1101) pdbcod,chnid,aaids(nxaa)(2:6),
     +              names1(nxaa),aaids(nxaa+1)(2:6),
     +              names1(nxaa+1),aaids(nxaa+2)(2:6),
     +              names1(nxaa+2),aaids(nxaa+3)(2:6),
     +              names1(nxaa+3),type(turnty(2,nxturn)),efimov,
     +              phi(nxaa+1),psi(nxaa+1),phi(nxaa+2),psi(nxaa+2),
     +              cadp3(nxaa),hbpres(nxturn),chi1(nxaa+1),
     +              chi1(nxaa+2),nxaa
 1101          format(a4,1x,a1,4(a5,a1,1x),a4,1x,a2,1x,4f6.1,1x,f4.1,
     +              1x,a1,2(f6.1),1x,i4)
            else
               write(outfi,1100)chnid,
c     +        ((aaids(j)(2:6),names1(j)),j=nxaa,nxaa+btrnsz)
     +        aaids(nxaa)(2:6),names1(nxaa),
     +        aaids(nxaa+1)(2:6),names1(nxaa+1),
     +        aaids(nxaa+2)(2:6),names1(nxaa+2),
     +        aaids(nxaa+3)(2:6),names1(nxaa+3)
     +        ,type(turnty(2,nxturn)), efimov, phi(nxaa+1), 
     +        psi(nxaa+1),phi(nxaa+2),psi(nxaa+2),cadp3(nxaa),
     +        hbpres(nxturn),chi1(nxaa+1),chi1(nxaa+2),nxaa
 1100         format(a1,4(a5,a1,1x),a4,1x,a2,1x,4f6.1,1x,f4.1,1x,a1,
     +        2(f6.1),1x,i4)
CPDBSUM RAL 14/3/05 -->
C-----        Convert single-letter a.a. codes to 3-letter residue names
              CALL AA2NAM(names1(nxaa),RSNAM1)
              CALL AA2NAM(names1(nxaa+1),RSNAM2)
              CALL AA2NAM(names1(nxaa+2),RSNAM3)
              CALL AA2NAM(names1(nxaa+3),RSNAM4)

C----         Set H-bond flag
              IF (hbpres(nxturn).EQ.'Y') THEN
                  HBFLAG = 'TRUE'
              ELSE
                  HBFLAG = 'FALSE'
              ENDIF

C----         If this is the very first beta turn, open the
C             bturns.promotif file and write headers
              IF (.NOT.HAVBET) THEN
                  OUTNAM = "bturns.promotif"
                  OPEN(UNIT=OUTBET, FILE=OUTNAM, STATUS='UNKNOWN',
     -                ERR=902)
                  CALL BETHED(outbet)

C----             Set flag
                  HAVBET = .TRUE.
              ENDIF

C----         Write out data to .promotif data file
              WRITE(outbet,1120) chnid, aaids(nxaa)(2:6), RSNAM1,
     -            aaids(nxaa+1)(2:6), RSNAM2,
     -            aaids(nxaa+2)(2:6), RSNAM3,
     -            aaids(nxaa+3)(2:6), RSNAM4, names1(nxaa),
     -            names1(nxaa+1),names1(nxaa+2),names1(nxaa+3),
     -            type(turnty(2,nxturn)), efimov, phi(nxaa+1), 
     -            psi(nxaa+1),phi(nxaa+2),psi(nxaa+2),cadp3(nxaa),
     -            HBFLAG,chi1(nxaa+1),chi1(nxaa+2),nxaa
 1120         FORMAT(A1,'::',4(A5,'::',A3,'::'),4A1,'::',A4,'::',A2,
     -            '::',4(F6.1,'::'),F4.1,'::',A5,'::',2(F6.1,'::'),I4)
CPDBSUM RAL 14/3/05 <--
           endif
         ENDIF
      endif
 2000 continue
c      write(outfi,*)
c      if(option.ne.2) write(outfig,*)'Gamma Turns'
      if(coption.ne.'L'.and.coption.ne.'l') write(outfig,*)'Gamma Turns'
      do 2005 nxturn = 1,nturns
         nxaa=turnst(nxturn)
         if(turnty(2,nxturn).gt.nbeta) then
c     else
C     GAMMA TURN
         chnid = aaids(nxaa)(1:1)
c         if (chnid.eq.'-')  then
c            chnid = ' '
c         end if 
         if (turnty(2,nxturn).le.nbeta+ngamma) then
            efivom = getefi(phi(nxaa+1),psi(nxaa+1),ppset,nphis)
            ngt = ngt + 1
c            if(option.eq.2) then
            if(coption.eq.'L'.or.coption.eq.'l') then 
               write(outfig,1201) pdbcod,chnid,aaids(nxaa)(2:6),
     +              names1(nxaa),aaids(nxaa+1)(2:6),
     +              names1(nxaa+1),aaids(nxaa+gtrnsz)(2:6),
     +              names1(nxaa+gtrnsz),type(turnty(2,nxturn))(2:8),
     +              efivom,phi(nxaa+1),psi(nxaa+1),cadp2(nxaa),
     +              chi1(nxaa+1),nxaa
 1201          format(a4,1x,a1,3(a5,a1,1x),a7,1x,a1,2f6.1,f4.1,1x,f6.1,
     +              1x,i4)
            else
               write(outfig,1200)chnid, aaids(nxaa)(2:6), 
     +           names1(nxaa), aaids(nxaa+1)(2:6),
     +           names1(nxaa+1),aaids(nxaa+gtrnsz)(2:6), 
     +           names1(nxaa+gtrnsz), 
     +           type(turnty(2,nxturn))(2:8), efivom,
     +           phi(nxaa+1), psi(nxaa+1),
     +           cadp2(nxaa),chi1(nxaa+1),nxaa
 1200            format(a1,3(a5,a1,1x),a7,1x,a1,2f6.1,f4.1,1x,f6.1,
     +              1x,i4)
CPDBSUM RAL 14/3/05 -->
C-----         Convert single-letter a.a. codes to 3-letter residue names
               CALL AA2NAM(names1(nxaa),RSNAM1)
               CALL AA2NAM(names1(nxaa+1),RSNAM2)
               CALL AA2NAM(names1(nxaa+gtrnsz),RSNAM3)

C----         If this is the very first gamma turn, open the
C             gturns.promotif file and write headers
              IF (.NOT.HAVGAM) THEN
                  OUTNAM = "gturns.promotif"
                  OPEN(UNIT=OUTGAM, FILE=OUTNAM, STATUS='UNKNOWN',
     -                ERR=902)
                  CALL GAMHED(outgam)

C----             Set flag
                  HAVGAM = .TRUE.
              ENDIF

C----         Write out data to .promotif data file
              WRITE(outgam,1220) chnid, aaids(nxaa)(2:6), RSNAM1,
     -            aaids(nxaa+1)(2:6), RSNAM2, aaids(nxaa+gtrnsz)(2:6),
     -            RSNAM3, names1(nxaa), names1(nxaa+1),
     -            names1(nxaa+gtrnsz), type(turnty(2,nxturn))(2:8),
     -            efivom, phi(nxaa+1), psi(nxaa+1),
     -            cadp2(nxaa),chi1(nxaa+1),nxaa
 1220         FORMAT(A1,'::',3(A5,'::',A3,'::'),3A1,'::',A7,'::',A1,
     -            '::',2(F6.1,'::'),F4.1,'::',F6.1,'::',I4)
CPDBSUM RAL 14/3/05 <--
              endif
c        else
         end if
      end if
c2000 continue
 2005 continue

CPDBSUM RAL 25/4/05 -->
C---- Close .promotif file, if open
      IF (HAVBET) CLOSE(UNIT=OUTBET)
      IF (HAVGAM) CLOSE(UNIT=OUTGAM)
      GO TO 999

C---- File errors
  902 CONTINUE
      WRITE (*, 7000) OUTNAM(1:LENSTR(OUTNAM))
 7000 FORMAT ('*** ERROR opening output file "',A,'"')

 999  CONTINUE
CPDBSUM RAL 25/4/05 <--

      end
c*******************************************************************************
c
      logical  function xeqy (x, y)
c
c     implicit none
c     implicit complex (a-z)
      real     x, y
c===============================================================================
c     This function returns TRUE is its two arguments (X and Y) are "equal",
c     else it returns "FALSE". How "equal" X and Y have to be is determined
c     by the parameter ACURCY, currently set to 0.0001.
c===============================================================================
      real       acurcy
      parameter (acurcy = 0.0001)
c
      xeqy = (abs(x-y) .lt. acurcy) 
      return
      end
c
c*******************************************************************************
c
c
      subroutine getstr (prompt, defans, string, ilen, ifail)
c
c     implicit complex(a-z)
c     implicit none
c===============================================================================
c     General purpose routine to prompt the user to type in a character
c     string. The prompt that will be written out to the user should be passed
c     to the routine in variable PROMPT, along with a default response in
c     DEFANS. The variable STRING will return the user's input to the calling
c     program, and ILEN will hold its length. The flag IFAIL will have one of
c     the following values:
c        IFAIL = 0 : everything OK;
c        IFAIL = 1 : variable PROMPT is empty;
c        IFAIL = 2 : variable DEFANS is empty;
c        IFAIL = 3 : user's response is too long to fit in STRING;
c        IFAIL = 4 : error occurred writing out prompt to terminal;
c        IFAIL = 5 : error occurred reading in reply from keyboard;
c        IFAIL = 6 : EOF detected reading in reply from keyboard.
c===============================================================================
      character*(*) prompt, defans, string
      integer       ilen, ifail
c
      character     line*255
      integer       lenp, lens, strlen, lend
c
 1000 format (A, ' <', A, '> - ')
 1010 format (A)
c
      ilen  = 0
      ifail = 0
c
      lenp  = strlen(prompt)
      if (lenp .le. 0)  then
         ifail = 1
         return
      end if
c
      lend  = strlen(defans)
      if (lend .le. 0)  then
         ifail = 2
         return
      end if
c
      lens  = len(string)
c
c     write (6, 1000, err = 904) prompt(1:lenp), defans(1:lend)
      read  (5, 1010, err = 905, end = 906) line
      ilen = strlen(line)
      if (ilen .eq. 0) then
         line = defans(1:lend)
         ilen = lend
      end if
      if (ilen .gt. lens)  then
         ifail = 3
         return
      end if
      string = line(1:ilen)
      ifail = 0
      return
c
  904 ifail = 4
      return
  905 ifail = 5
      return
  906 ifail = 6
      return
c
      end 
c
c*******************************************************************************
c
      integer  function strlen (string)
c
c     implicit complex (a-z)
c     implicit none
c
c     find out the real length of the string
c
      character*(*) string
c
      integer nxch
c
      do 200 nxch = len(string), 1, -1
      if (string(nxch:nxch) .ne. ' ') goto 300
  200 continue
      nxch = 0
  300 continue
      strlen = nxch
      return
      end
c
c*******************************************************************************
c
      character*1 function getefi (phi, psi, ppset, nphis)
c
c     implicit none
c     implicit complex (a-z)
c
      real      phi, psi
      integer   ppset, nphis
      dimension ppset(nphis,nphis)
c
      character*1 eclass(7)
      integer     phipl, psipl, angplt
      logical     xeqy
c
      data eclass /'B', 'P', 'A', 'a', 'L', 'G', 'E'/
c
      if (xeqy(phi,999.9) .or. xeqy(psi,999.9)) then
         getefi = ' '
         return
      end if
c
      phipl = angplt(phi, nphis)
      psipl = nphis + 1 - angplt(psi, nphis)
      if (phipl .lt. 1 .or. phipl .gt. nphis .or.
     +    psipl .lt. 1 .or. psipl .gt. nphis .or.
     +    ppset(psipl,phipl) .eq. 0) then
         getefi = ' '
      else
         getefi = eclass(ppset(psipl,phipl))
      end if
      return
      end
c
c*******************************************************************************
c
      integer function angplt (angle, nphis)
c
c     implicit none
c     implicit complex (a-z)
c
      real    angle
      integer nphis
c
      angplt = int((angle + 180.0) / 10.0) + 1
      if (angplt .gt. nphis) angplt = nphis
      return
      end
c
c*******************************************************************************
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
CPDBSUM RAL 14/3/05 -->
C**************************************************************************
C
C  SUBROUTINE BETHED  -  Write header records to bturns.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE BETHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 21)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN', 'RESNUM1', 'RESNAM1', 'RESNUM2', 'RESNAM2',
     -     'RESNUM3', 'RESNAM3', 'RESNUM4', 'RESNAM4', 'SEQUENCE',
     -     'TURNTYPE', 'EFIMOV_CLASS', 'PHI1', 'PSI1', 'PHI2', 'PSI2',
     -     'CA_DIST', 'HBOND_FLAG', 'CHI1_1', 'CHI1_2', 'SEQ_POS' /

C---- Define field id tag
      IDTAG = 'BT_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:BETA_TURNS'
 10   FORMAT(100A)

C---- Write the field names and identifiers
      WRITE(IUNIT,10) (IDTAG, FLDNAM(IFIELD)(1:LENSTR(FLDNAM(IFIELD))),
     -    '::', IFIELD = 1, NFIELD - 1), IDTAG,
     -    FLDNAM(NFIELD)(1:LENSTR(FLDNAM(NFIELD)))

      RETURN
      END

C--------------------------------------------------------------------------
C**************************************************************************
C
C  SUBROUTINE GAMHED  -  Write header records to gturns.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE GAMHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 15)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN', 'RESNUM1', 'RESNAM1', 'RESNUM2', 'RESNAM2',
     -     'RESNUM3', 'RESNAM3', 'SEQUENCE', 'TURNTYPE', 'EFIMOV_CLASS',
     -     'PHI', 'PSI', 'CA_DIST', 'CHI1', 'SEQ_POS' /

C---- Define field id tag
      IDTAG = 'GT_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:GAMMA_TURNS'
 10   FORMAT(100A)

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
CPDBSUM RAL 14/3/05 <--
