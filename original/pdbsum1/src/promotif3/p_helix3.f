      PROGRAM p_helix3
C ***********************************
C *** Modified version for PDBsum ***
C ===================================
C
C***********************************************************************
C*      NAME: p_helix3.f                                             *
C*  FUNCTION: Program to identify helices in proteins                  * 
C* COPYRIGHT: (1990-1996) University College London                    *
C*---------------------------------------------------------------------*
C*    AUTHORS: T. P. FLORES, D. NAYLOR, F. M. G. RICHARDSON and        *
C*             E. G. HUTCHINSON                                        *
C*      DATE:  31/1/90                                                 *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
c* 17/07/91  DN     a) included routine LOADHAT;
c                   b) included routine FNDAT2;
c                   c) replaced calls to FINDAT within helix interaction
c                      loop with calls to FNDAT2.
c  26/08/91  DN     d) dotprod format increased
c  11/12/91  DN     e) changes made to enable program to use list of input
c                      files
c  20/12/91  DN     f) changes (e) completed! 
c  04/02/92  DN     g) alterations for -ve residue numbers in READST
c  21/04/92  DN     h) major corrections to workings of program - fixing of
c                      bugs in/rewriting of READBK, READST, FINDAT, LOADHAT
c                      New routines GETFIRSTATOM, GETLASTATOM. 
c  01/05/92  DN     i) major corrections to workings of program - fixing of
c                      bugs in/rewriting READST, TRNANG. Subroutines SECSTR
c                      and CAPTL are deleted
c  15/05/92  DN     j) major corrections in order to recognize actual
c                      beginings and ends of helices. 
c  20/05/92  DN        continued . . . 
c  26/05/92  DN     k) completed (re)implementation of SND.
c  26/05/92  DN     l) started calculations of helix axes within subroutine
c                        FINDHELICES
c  27/05/92  DN          continued . . . 
c  28/05/92  DN          continued . . . 
c  29/05/92  DN          continued . . . 
c                   m) altered format 500 and other modifications to include
c                        the nos of interacting residues from BOTH helices
c  24/10/94  FR     correct error in the closest approach routine. A
c                   new routine newcls was added. (Ph.D thesis of Frances 
C                   Richardson)
c  25/10/94  GH     add in sequence information about helix (up to 55 residues)
c  14/02/95  GH     to allow use for more than one protein
c  09/05/95  GH     either write data all to one file or to write separate 
c                   files for each protein
c  13/03/96  GH     correct naming of file in multiple protein option
c  20/03/96  GH     write out helix sequence for every option
c  27/03/96  GH     reorganize multiple/single file handling
c  03/04/96  GH     only read in 1st model for NMR structures
c  29/04/96  GH     modify to allow command line options
c  31/07/96  GH     added in chain letter to interaction data
c  21/08/96  GH     changed meth to 3 i.e use modified Kabsch and Sander 
c                   assignments for consistency
c  04/09/96  GH     changed msplits in .par to 40
c                   changed meth back to 2 for testing 2gls output
c  06/11/96  GH     changed maxaph to 700 in p_helix.par -to cope with 1jun
c  06/11/96  GH     temporarily output hlxdebug.dat with helix vector info
c  06/11/96  GH     temporarily output p vector for Alex
c  18/12/96  GH     now includes PI helices
c  03/06/97  GH     changed meth back to 2 i.e. use strict Kabsch and Sander
c                   assignments for consistency within promotif_html
c  01/11/00  GH     removed 'x' specifier from format statements
c                   replaced atan2d with atan2
c                   changed ierror from integer to logical
c  21/11/00  GH     updated to read new sstruc3 format
c  20/03/01  GH     removed double definition of k in main program
c==============================================================================
C
C Compiling under g77:-
C
C f77 -Wimplicit -fbounds-check -o p_helix3 p_helix3.f
C
c***********************************************************************

      include 'p_helix.par'
c
      integer    ires(matm), iatm(matm), ihlx(2,mhlx), cres(mnch),
     +           gap, lc(matm), kc(matm), i, ii,j, nhlx, nlt, k,
     +           ip, nss,  m, n, lnlt, knlt, ih, meth, natm, nres,
     +           nch, ier, strlen, namlen,
     +           atomlist(mhlx,maxaph), iatom, jatom, debug,
     +           pairlist(mintp,2), fistat, ok, finish, opnera, opnerb,
     +           opnerc, opnerd, opnere, opnerf, iocode, namfi, endnfi,
     +           brkfi, strkfi, endfi, num1(mhlx), num2(mhlx), seqlen,
     +           gethlxpos, ijh, iih,option,istart,iend,nps
CPDBSUM RAL 29/3/05 -->
      INTEGER   HELIX1, HELIX2, HLXNO, LENSTR, NLOST, NRESID, OPNERG,
     -          OPNERH, OUTHEL, OUTHIN
CPDBSUM RAL 29/3/05 <--
c     integer    inhl(mhlx,mhlx),ifail
c
      real       x(3,matm), e(3), st(3), fn(3), p(3,matm), pc(3,matm),
     +           a(3), b(3), bg(3,mhlx), ed(3,mhlx), bj(3), bk(3),
     +           ej(3), ek(3), pck(3), pcj(3), omega, rise, rpt, dotprd,
     +           ddstnc, pi, rmsr, rmsp, dstnce, dist, sk, sj, rlen,
     +           pitch, r1, rc, tlen, dis, crv, dna(3), dnb(3),
     +           dndot, dnmaga, dnmagb, dnmagc, dndum1(3), dndum2(3),
     +           dnomg1, dnomg2, dntht1, dntht2, dnaxb(3), dntemp,
     +           haxis(3), val1, val2, val3, val4
      real       ekbk(3), ejbj(3), fj, fk, phk(3), phj(3)
c     real       val5, val6
c
      character  atmn*4, resn*3, chain*1, hch*1, snd*5000,lsnd*5008, 
     +           brkid*4, rpc*1, line*130, ins*1,lb1*1, lb2*1, 
     +           htype(mhlx)*1, chnid*1, chainlist(mintp,2)*1,
     +           indirn*1, space*1, inslist(mintp,2)*1, 
CPDBSUM RAL 21/4/05 -->
C     +           reslist*(12*mintp), hseq(mhlx)*300
     +           reslist*(24*mintp), hseq(mhlx)*300
CPDBSUM RAL 21/4/05 <--
      character  ks(maxres)*1, kspl(maxres)*1, spah(maxres)*1,
     +           sp310h(maxres)*1, seqbcd(maxres)*6, chnbrk(maxres)*1,
     +           seq(maxres)*1,pdbnam*4,coption, sppih(maxres)*1
      character*78 brknam,sstnam,pdbcod,finam,outfil,filen,outfl2
CPDBSUM RAL 29/3/05 -->
      CHARACTER*1   AA1, AA2
      CHARACTER*3   RSNAM1, RSNAM2
      CHARACTER*512 OUTNAM
CPDBSUM RAL 29/3/05 <--

      logical    plot, newpair, finuse, xeqy, allok,ENSEMB,IERROR
      real       phiang(maxres), psiang(maxres)
      logical    hlxok(mhlx)
c
      dimension  atmn(matm), resn(matm), chain(matm), hch(mhlx),
     +           rpc(2), ins(matm), lb1(mhlx), lb2(mhlx)
c
      parameter (pi = 3.1415926536, indirn = '@', space = ' ')
      parameter (ok = 0, finish = 2, endfi = 4, endnfi = 6,
     +           opnera = 1, opnerb = 3, opnerc = 5, opnerd = 7,
CPDBSUM RAL 29/3/05 -->
C     +           opnere = 9, opnerf = 11, namfi = 15, brkfi = 16,
C     +           strkfi = 17)
     -           opnere = 9, outhel = 10, outhin = 11, opnerf = 11,
     -           opnerg = 12, opnerh = 13, namfi = 15, brkfi = 16,
     -           strkfi = 17)
CPDBSUM RAL 29/3/05 <--

      DATA FISTAT/OK/,FINUSE/.FALSE./,OPTION/1/
c
  120 format ( a80 )
c 510 format (a4, i3, i3, a1, a1, f5.1, a1, a1, i2, f6.1 )
  510 format (i3, 1x,i3, 1x,a1, 1x,a1, 1x,f5.1, 1x,a1, 1x,a1, 1x,
     ~     f6.1, 2x,3(1x,i2),1x,a1,a1 )
 511  format (a4,1x,i3, 1x,i3, 1x,a1, 1x,a1, 1x,f5.1, 1x,a1, 1x,a1, 1x,
     ~     f6.1, 2x,3(1x,i2) )

c 515  format(50i5)
  530 format (a125)
  531  format (a130)
c 1010 format (/,'Main chain atom numbers for helix ',i3,
c     + ' determined by FINDAT')
c 1020 format (20(16I5,/))
c 1024 format ('Actual atom id numbers are:')
c 1021 format ('Structure too short, axis not calculated for hlx ',i4,/)
c 1030 format ('hx1 hx2 rs1 rs2 at1 at2 d:',2i4, 2(' (',a1,i4,a1,')'),
c     +         2i6, f6.2)
 1040 format ('Error opening file of file names "',a,'"')
 1050 format ('Error opening Brookhaven file "',a,'"')
 1060 format ('Error opening 2D structure file "',a,'"')
 1070 format ('Error opening output file "',a,'"')
 1234 format ('Using meth = ',i1, '     ',a)
c 7020 format ('Error in IPROMPT: Ifail = ', I2)
c 4000 format (/,'contents of SND: "',a,'"')
c 4005 format ('Residue ',i4,': Ks="',a1,'" Kspl="',a1,'"',
c     + ' Spah="',a1,'" Sp3h="',a1,'" phi,psi= ',2f10.1)
 4010 format ('Error looking through secondary structure list SND',
     +        'for helix ',i3)
 3334 format ('   interacting helices are ',i3,' and ',i3)
 3336 format ('   for helices ',i4,' and ',i4,', dnmaga = ',
     +        f10.4,' dnmagb = ',f10.4)
 3434 format ('Too many pairs of interacting residues for helices ',
     +        i3,' and ',i3,/,
     +        ' MAX allowed = ',I3)
 123  format('Promotif output for protein :',a)
c
c

c   - Use standard Kabsch and Sander structure assignments
      meth = 2
c     21/8/96 Use modified Kabsch and Sander assignments
c      meth=3
c      call iprompt ('Enter value for METH', meth, 1, 3, ifail)
c      if (ifail.ne.0) then
c         write (6, 7020) ifail
c      end if
c
c   - Distance for interacting helices (ANGSTROMS)
      dist = 4.5
c
c   - Use all atoms on main chain for calculation of axes
      gap  = 3
      plot = .false.
c
      debug = 0
c      open(7, file='hlxtorsn.dat', status='unknown', err=55)
c      open(8, file='hlxdebug.dat', status='unknown', err=55)
c      open(9, file='hlxinfo.dat', status='unknown', err=55)
   55 continue
c
      read (*, '(A)', iostat=iocode) brknam
      IF(BRKNAM(1:1).EQ.'%') THEN
         READ(*,'(A1)') COPTION
         ENSEMB=.TRUE.
         IF(COPTION.EQ.'L'.OR.COPTION.EQ.'l') THEN
            OUTFIL='HELICES'
            OPEN (UNIT=2, file=OUTFIL,status='unknown', iostat=iocode)
            IF(IOCODE.NE.0) FISTAT=OPNERB
            OUTFL2='HELIX_INT'
            OPEN(UNIT=3,file=outfl2,status='unknown',iostat=iocode)
            if(iocode.ne.0) fistat=opnerb
         ENDIF
      ELSE
         ENSEMB=.FALSE.
         coption='p'
      ENDIF

c     if ensemble remove % from file name
      IF(ENSEMB) THEN
         FINAM=BRKNAM
         FINUSE=.TRUE.
         NAMLEN=INDEX(BRKNAM,SPACE)-1
         DO 180 I=1,NAMLEN-1
            FINAM(I:I)=FINAM(I+1:I+1)
 180     CONTINUE
         FINAM(NAMLEN:NAMLEN)=SPACE
         OPEN(NAMFI,FILE=FINAM,STATUS='OLD',IOSTAT=IOCODE)
      ENDIF

      if (iocode.lt.0)  then
         close(namfi)
         finuse = .false.
         fistat = endnfi
      else
C         namlen = strlen(brknam)
      end if
      
c
c     open Brookhaven file
 200  CONTINUE
      IF(FISTAT.EQ.OK) THEN
         IF(ENSEMB) THEN
            READ(NAMFI,'(A)',IOSTAT=IOCODE) BRKNAM
            NAMLEN=INDEX(BRKNAM,SPACE)-1
            IF(IOCODE.LT.0) THEN
               CLOSE(NAMFI)
               FINUSE=.FALSE.
               FISTAT=ENDNFI
            ELSE
               WRITE(*,*) BRKNAM
            ENDIF
c         ELSE
c            BRKNAM=FINAM
         ENDIF
      ENDIF

      IF(FISTAT.EQ.OK)  then
         open (brkfi, file=brknam, status='old', iostat=iocode)
         if (iocode.ne.0)  fistat = opnerb
      end if
c
c     open SST file
      if (fistat.eq.ok)  then
         call getnam(brknam,istart,iend,ierror)
         pdbcod=brknam(istart:iend)
         nps=index(pdbcod,' ')
         sstnam=pdbcod(1:nps-1)//'.sst'
         i=index(brknam,'.')
         if(i.gt.4) then
            pdbnam=brknam(i-4:i-1)
c            pdbnam=brknam(i-5:i-2)
         else
            if(i.ne.0) then
               pdbnam=brknam(1:4)
            else
               k=index(brknam,' ')
               pdbnam=brknam(k-4:k-1)
            endif
         endif
         open (strkfi, file=sstnam, status='old', iostat=iocode)
         if (iocode.ne.0)  fistat = opnerc
      end if

c     Read coordinates etc. and find where the helices are
      IF(FISTAT.EQ.OK) THEN
         call readbk (x, atmn, resn, chain, iatm, ires, brkfi,
     +                natm, nres, cres, nch, ins, brkid)
         call readst (strkfi, nhlx, brkid, phiang, psiang,ks, kspl, 
     +                spah, sp310h, sppih, seqbcd, chnbrk, seqlen, seq)
         close (strkfi)
         if (debug.eq.1)  then
            i = 1
c            write (8, 4005) i, ks(i), kspl(i), spah(i),
c     +                      sp310h(i), phiang(i), psiang(i)
            i = 10
c            write (8, 4005) i, ks(i), kspl(i), spah(i),
c     +                      sp310h(i), phiang(i), psiang(i)
            i = seqlen - 10
c            write (8, 4005) i, ks(i), kspl(i), spah(i),
c     +                      sp310h(i), phiang(i), psiang(i)
            i = seqlen
c            write (8, 4005) i, ks(i), kspl(i), spah(i),
c     +                      sp310h(i), phiang(i), psiang(i)
         end if
         call findhelices (ks, kspl, spah, sp310h, sppih, seqbcd, 
     ~                     chnbrk, seqlen, meth, gap,
     +                     ihlx, nhlx, hch, lb1, lb2, num1, num2,
     +                     htype, chain, ires, ins, atmn, x,
     +                     natm, snd, nss, seq, hseq)
c
c         if (debug.eq.1) write (8, 4000) snd(1:strlen(snd))
c      end if
c
c      if (fistat.ne.ok)  goto 456
c
CPDBSUM RAL 21/4/05 -->
      if (nhlx .le. 0)  then
         print*, ' no helices in this structure!'
         goto 456
      endif
CPDBSUM RAL 21/4/05 <--
      lsnd(1:nss+8) = '    '//snd(1:nss)//'    '
c
      ip   = 0
      crv  = 0.0e0
      r1   = 0.0e0
      rc   = 0.0e0
c
c     Load ATOMLIST with the id numbers of all atoms in all helices
      call loadhat (atmn, ires, chain, natm, lb1, lb2, ins,
     +              num1, num2, hch, nhlx, atomlist, 2)
c
c     -------------------
c     part I: Helix data
c     -------------------

      if (nhlx .le. 0)  then
         print*, ' no helices in this structure!'
         goto 456
      else
c         IF(OPTION.ne.2) THEN
         IF(COPTION.NE.'L'.AND.COPTION.NE.'l') THEN
            outfil=pdbcod(1:nps-1)//'.hlx'
c     outfil='HELIX'
CPDBSUM RAL 14/4/05 -->
C            open (2, file=outfil, status='unknown', iostat=iocode)
            open (2, status='scratch', iostat=iocode)
CPDBSUM RAL 14/4/05 <--
            if (iocode.ne.0)  fistat = opnere
            write(2,123) pdbcod(1:nps-1)
            write(2,*)
            if(nhlx.ge.1) write(2,*)'Helix data:'
CPDBSUM RAL 29/3/05 -->
C----       Open helices data file
            outnam = "helices.promotif"
            open(unit=outhel, file=outnam, status='unknown',
     -          iostat=iocode)
            if (iocode.ne.0)  fistat = opnerg
            call helhed(outhel)

C----       Open helix interactions data file
            outnam = "helint.promotif"
            open(unit=outhin, file=outnam, status='unknown',
     -          iostat=iocode)
            if (iocode.ne.0)  fistat = opnerh
            call hinhed(outhin)
CPDBSUM RAL 29/3/05 <--
         ENDIF
      end if
c
      do 10 i = 1, nhlx
      rlen  = 0.0e0
      rise  = 0.0e0
      rpt   = 0.0e0
      pitch = 0.0e0
      rmsp  = 0.0e0
      crv   = 0.0e0
      r1    = 0.0e0
      rc    = 0.0e0
c
c     Check each helix for dodgy torsion angles
      call checktorsions (phiang, psiang, ihlx, i, htype(i), allok,
     +                    meth, brkid)
      hlxok(i) = allok
      if (allok)  then
         call findat (atmn, ires, chain, ins, natm, hch(i),
     +                num1(i), lb1(i), hch(i), num2(i), lb2(i), gap,
     +                lc, nlt)
         if (debug.eq.2) then
c            write (8, 1010) i
c            write (8, 1020) (lc(j), j=1,nlt)
c            write (8, 1024)
c            write (8, 1020) (iatm(lc(j)), j=1,nlt)
         end if
      else
         goto 7
      end if
c
c     Calculate helix axis
      if ((ihlx(2,i)-ihlx(1,i)+1).ge.minrph)  then
         call hlxaxs (x, lc, nlt, gap, i, e, haxis, st, fn, p, pc)
      else
c         write (8,1021) i
         goto 7
      end if
c
c     Calculate average rise
      rise = 0.0e0
      do 1, j = 1,nlt-gap
      rise = rise + dstnce(pc,j,j+gap)
    1 continue
      rise = rise / real(nlt-gap)
c
c     Calculate residues per turn
      rpt = 0.0e0
      do 3, j= 1, nlt-gap
      do 2, k= 1, 3
      a(k) = x(k,lc(j)) - pc(k,j)
      b(k) = x(k,lc(j+gap)) - pc(k,j+gap)
    2 continue
      if (debug.eq.3) then
         print*,'acos call 1: i,j=', i, j
         if (i.ge.25.and.j.ge.53) then
         print*,'    dotprd(a,b)=',dotprd(a,b)
         print*,'  ddstnc xj-pcj=',ddstnc( x(1,lc(j)), pc(1,j) ) 
         print*,'ddstnc xjg-pcjg=',ddstnc( x(1,lc(j+gap)), pc(1,j+gap)) 
         endif
      endif
      val1 = dotprd(a,b)
      val2 = ddstnc(x(1,lc(j)),pc(1,j))
      val3 = ddstnc(x(1,lc(j+gap)),pc(1,j+gap))
      val4 = val1/(val2*val3)
      if (val4.gt.1.0D0) val4 = 1.0D0
      if (val4.lt.-1.0D0) val4 = -1.0D0
      rpt = rpt +acos(val4)*180.0/pi
    3 continue
      rpt = real(nlt-gap)*360.0e0 / rpt
c
c     Calculate pitch
      pitch = rise*rpt
c
c     Calculate RMS deviation of vector for ideal
      rmsp = 0.0e0
      do 5, j= gap+1, nlt-gap
      do 4, k= 1, 3
      a(k) = p(k,j-gap)
      b(k) = pc(k,j) - x(k,lc(j))
    4 continue
c     if (debug.eq.3) print*,'acos call 2'
      rmsr = acos(dotprd(a,b) / (dstnce(p,j-gap,j-gap) *
     +             ddstnc(x(1,lc(j)), pc(1,j))))*180.0e0/pi
c     val5 = dstnce(p,j-gap,j-gap)
c     val6 = val1/(val5*val2)
c     rmsr = acos(val6)*180.0/pi
      rmsp = rmsp + rmsr*rmsr
    5 continue
      rmsp = sqrt( rmsp/real(nlt-2*gap) )
c
c     Store start and ends of helices for later use and calculate length
      rlen = 0.0e0
      do 6, j= 1, 3
      tlen = fn(j) - st(j)
      rlen = rlen + tlen*tlen
      bg(j,i) = st(j)
      ed(j,i) = fn(j)
    6 continue
      rlen = sqrt( rlen )
c
c     Calculate type of helix as Barlow.
c     call hltype ( x, lc, nlt, gap, crv, r1, rc, plot )
c
c     Find local secondary structure around helix
    7 continue
      ip = gethlxpos(i,snd)
      if (ip .eq. 0)  then
          write (6, 4010) i
          ip = nss + 5
      end if
c
c     Write output

c 160 format (a4, a1, i3, i3, a1, i3, i3, a9, 8f5.1 )
c     write (1,160) brkid, snd(ip:ip), num1(i), num2(i),
c    +              hch(i), ihlx(2,i)-ihlx(1,i)+1, i, lsnd(ip:ip+8),
c    +              rlen, rise, rpt, pitch, rmsp, crv, r1, rc
      chnid = hch(i)
      if (chnid.eq.'-')  then
         chnid = ' '
      end if
  500 format (i3,1x,a1, 2( 1x,a1, i4, a1), 1x,a1,1x,i3, 1x,a9,
     +        4(1x,f6.2),1x, f5.1,1x,a55)
 501  format (a4,1x,i3,1x,a1, 2( 1x,a1, i4, a1), 1x,a1,1x,i3, 1x,a9,
     +        4(1x,f6.2),1x, f5.1,1x,a55)
c      if(option.eq.2) then
      if(coption.eq.'L'.or.coption.eq.'l') then
         write (line,501) pdbnam,i,chnid, hch(i), num1(i), lb1(i),
     +        hch(i), num2(i), lb2(i), htype(i),
     +        ihlx(2,i)-ihlx(1,i)+1, lsnd(ip:ip+8), rlen, rise, 
     +        rpt, pitch, rmsp, hseq(i)(1:55)
      write (2,531) line
      else
         write(line,500) i,chnid, hch(i), num1(i), lb1(i),
     +        hch(i), num2(i), lb2(i), htype(i),
     +        ihlx(2,i)-ihlx(1,i)+1,lsnd(ip:ip+8), rlen, rise, 
     +        rpt, pitch, rmsp,hseq(i)(1:55)
      write (2,530) line(1:125)
CPDBSUM RAL 29/3/05 -->
C---- Get the start and end residues of this helix
      AA1 = hseq(i)(1:1)
      NRESID = ihlx(2,i)-ihlx(1,i)+1
      IF (NRESID.GT.0 .AND. NRESID.LE.300) THEN
          AA2 = hseq(i)(NRESID:NRESID)
          IF (AA2.EQ.' ') AA2 = 'X'
      ELSE
          AA2 = 'X'
      ENDIF

C---- Convert single-letter a.a. codes to 3-letter residue names
      CALL AA2NAM(AA1,RSNAM1)
      CALL AA2NAM(AA2,RSNAM2)

C---- Determine this helix's number in its chain
      CALL GETHNO(chnid,hch,mhlx,i,HLXNO)

C---- Write out data to .promotif data file
      WRITE(outhel,620) chnid, HLXNO, num1(i), lb1(i), RSNAM1, num2(i),
     -    lb2(i), RSNAM2, htype(i), NRESID, lsnd(ip:ip+8), rlen, rise, 
     -    rpt, pitch, rmsp, hseq(i)(1:LENSTR(hseq(i)))
 620  FORMAT(A1,'::',I3,'::',2(I4,A1,'::',A3,'::'),A1,'::',I3,'::',
     -    A9,'::',4(F6.2,'::'),F5.1,'::',A,'::')
CPDBSUM RAL 29/3/05 <--
      endif
c
c     next helix
   10 continue
c
c
C     if (plot) call endplt
c
c     --------------------------------
c     PART II: Helix interactions
c     --------------------------------
c
      call loadhat (atmn, ires, chain, natm, lb1, lb2, ins,
     +              num1, num2, hch, nhlx, atomlist, 1)
c
c     Calculate helix interactions
      if (nhlx .lt. 2) then
         print*, ' only one helix in this structure!'
         goto 456
      end if
      if(nhlx.gt.1) then
c         if(option.ne.2) then
         if(coption.ne.'L'.and.coption.ne.'l') then
            write(2,*)
            write(2,*)'Helix Interactions:'
         endif
      endif
c
c     begin loop over 1st helix
      do 25 i = 1, nhlx-1
         if ((ihlx(2,i)-ihlx(1,i)+1).lt.minrph.or.
     +       (.not.hlxok(i))) goto 25
         call fndat2 (atomlist, i, lc, lnlt)
c
c        begin loop over 2nd helix
         do 20 j = i+1, nhlx
            if ((ihlx(2,j)-ihlx(1,j)+1).lt.minrph.or.
     +          (.not.hlxok(j)))  goto 20
            call fndat2 (atomlist, j, kc, knlt)
c
c           do the helices interact? 
c           IH counts the no of interacting pairs of residues for this pair 
c              of helices: as each new pair of residues is identified, their
c              chain letters, residue nos, and insertion codes are stored in
c              the CHAINLIST(), PAIRLIST(), and INSLIST() arrays. Once all 
c              pairs are identified, these arrays are examined and the number
c              of unique residues from each helix are counted in IIH and IJH.
            ih = 0

CPDBSUM RAL 21/4/05 -->
C----       Initialise count of lost interactions
            NLOST = 0
CPDBSUM RAL 21/4/05 <--
c
c           begin loop over each atom in 1st helix
            do 111 m = 1, lnlt
               iatom = lc(m)
c
c              begin loop over each atom in 2nd helix
               do 211 n = 1, knlt
                  jatom = kc(n)
c
c                 determine the distance between the two atoms
                  dntemp = dstnce(x,iatom,jatom)
                  if (dntemp .le. dist) then
                     newpair = .true.
                     do 311 k = 1,ih
                     if (chainlist(k,1) .eq. chain(iatom)  .and.
     +                    pairlist(k,1) .eq. ires(iatom)   .and.
     +                     inslist(k,1) .eq. ins(iatom)    .and.
     +                   chainlist(k,2) .eq. chain(jatom)  .and.
     +                    pairlist(k,2) .eq. ires(jatom)   .and.
     +                     inslist(k,2) .eq. ins(jatom) )
     +                                       newpair = .false.
  311                continue
                     if (newpair) then
c                        if (debug.eq.2) write (8, 1030) i, j,
c     +                               chain(iatom), ires(iatom),
c     +                               ins(iatom),
c     +                               chain(jatom), ires(jatom),
c     +                               ins(jatom),
c     +                               iatm(iatom), iatm(jatom), dntemp
CPDBSUM RAL 21/4/05 -->
                        IF (ih.GE.mintp) THEN
                            NLOST = NLOST + 1
                        ELSE
CPDBSUM RAL 21/4/05 <--
                            ih = ih + 1
c                       print*,' IH = ',ih
CPDBSUM RAL 21/4/05 -->
C                            if (ih .ge. mintp) then
C                               write (6,3434) i,j, mintp
C                               stop 
C                            end if
CPDBSUM RAL 21/4/05 <--
                            chainlist(ih,1) = chain(iatom)
                            pairlist(ih,1)  = ires(iatom)
                            inslist(ih,1)   = ins(iatom)
                            chainlist(ih,2) = chain(jatom)
                            pairlist(ih,2)  = ires(jatom)
                            inslist(ih,2)   = ins(jatom)
CPDBSUM RAL 21/4/05 -->
                        ENDIF
CPDBSUM RAL 21/4/05 <--
                     end if
                  end if
  211          continue
  111       continue
c           If helices don't interact skip them
            if (ih .lt. 1) goto 20
CPDBSUM RAL 21/4/05 -->
C----       If maximum number of interactions exceeded, show warning
            IF (NLOST.GT.0) THEN
                PRINT 103, mintp, NLOST
 103            FORMAT('Warning. Maximum number of helix-helix interac',
     -              'tions (',I6,') exceeded. ',I6,' interactions lost')
            ENDIF
CPDBSUM RAL 21/4/05 <--
            Call Countuniqueresidues(chainlist, pairlist, inslist,
     +                               ih, iih, ijh, brkid, reslist)
c           Calculate closest approach
            do 13, k= 1, 3
               bk(k) = bg(k,i)
               ek(k) = ed(k,i)
               bj(k) = bg(k,j)
               ej(k) = ed(k,j)
   13       continue
            call clsaph (i, j, bk, ek, bj, ej, pck, pcj,
     +                   sk, sj, .false., ekbk, ejbj)
c           Calculate omega angle
            do 14, k= 1, 3
               a(k) = pck(k) + bk(k) - ek(k)
               b(k) = pcj(k) + bj(k) - ej(k)
   14       continue
            call trnang (a, pck, pcj, b, omega, ier)
            if (ier.ne.0)  then
                write (6,3334) i,j
c                write (8,3334) i,j
            end if
c           Calculate distance between points of closest approach, adjust if
c           outside helix
            if (sk .lt. 0.0e0) then
               rpc(1) = 'N'
               fk = 0.0
            else if (sk .gt. 1.e0) then
               rpc(1) = 'C'
               fk = 1.0
            else
               rpc(1) = 'I'
               fk = sk
            end if
            if (sj .lt. 0.0e0) then
               rpc(2) = 'N'
               fj = 0.0
            else if (sj .gt. 1.e0) then
               rpc(2) = 'C'
               fj = 1.0
            else
               rpc(2) = 'I'
               fj = sj
            end if
c     distance between helices is recalculated by fmgr 2.11.93
            dis = ddstnc(pck, pcj)
            if(rpc(1).EQ.'C'.OR.rpc(1).EQ.'N'.OR.
     -           rpc(2).EQ.'C'.OR.rpc(2).EQ.'N') then
               do 15 ii=1, 3
                  phk(ii) = bk(ii) + fk*ekbk(ii)
                  phj(ii) = bj(ii) + fj*ejbj(ii)
 15            continue
               call newcls(bk, ek, bj, ej, pck, pcj, ekbk, ejbj, phk, 
     +                     phj, rpc, dis)
            endif

c
c           calculate angle and dot product between this interacting pair
c           of helices.
c
            do 6001 k=1,3
            dna(k)    = ed(k,i) - bg(k,i)
            dnb(k)    = ed(k,j) - bg(k,j)
            dndum1(k) = 0.0
            dndum2(k) = 0.0
 6001       continue
            dndot   = dotprd (dna, dnb)
            dnmaga  = sqrt( dna(1)**2 + dna(2)**2 + dna(3)**2)
            dnmagb  = sqrt( dnb(1)**2 + dnb(2)**2 + dnb(3)**2)
            if (xeqy(dnmaga,0.0).or.xeqy(dnmagb,0.0))  then
                write(6,3336) i, j, dnmaga, dnmagb
c                write(8,3336) i, j, dnmaga, dnmagb
            end if
c           if (debug.eq.3) print*,'acos call 3'
            dntht1 = acos( dndot / (dnmaga*dnmagb) )
            dntht1 = 180.0 * dntht1 / pi
c
            call axb (dna, dnb, dnaxb)
            dnmagc = sqrt( dnaxb(1)**2 + dnaxb(2)**2 + dnaxb(3)**2 )
            dntht2 = asin( dnmagc / (dnmaga*dnmagb) )
            dntht2 = 180.0 * dntht2 / pi
            dnomg1 = 0.0
            dnomg2 = 0.0
c
c           call trnang (bg(1,i), ed(1,i), bg(1,j), ed(1,j), dnomg1,
c    +                   ier)
c           call trnang (ed(1,i), bg(1,i), ed(1,j), bg(1,j), dnomg2,
c    +                   ier)
c------------------------------------------------------------------------------
c           Write out data 
c
c            if(option.ne.2) then
            if(coption.ne.'L'.and.coption.ne.'l') then
               write (2,510) i, j, htype(i), htype(j), dis,
     +              rpc(1), rpc(2), omega, ih,
     +              iih, ijh,hch(i),hch(j) 
CPDBSUM RAL 29/3/05 -->
C---- Determine the helix numbers in their respective chains
      CALL GETHNO(hch(i),hch,mhlx,i,HELIX1)
      CALL GETHNO(hch(j),hch,mhlx,j,HELIX2)

C---- Write out data to .promotif data file
      WRITE(outhin,640) hch(i), hch(i), HELIX1, htype(i), hch(j),
     -    HELIX2, htype(j), dis, rpc(1), rpc(2), omega, ih, iih, ijh
 640  FORMAT(A1,'::',2(A1,'::',I3,'::',A1,'::'),F5.1,'::',A1,'::',A1,
     ~     '::',F6.1,'::',3(I3,'::'))

C---- If helices belong to different chains, write out again for the
C     second chain
      IF (hch(i).ne.hch(J)) THEN
          WRITE(outhin,640) hch(j), hch(i), HELIX1, htype(i), hch(j),
     -        HELIX2, htype(j), dis, rpc(1), rpc(2), omega, ih, iih, ijh
      ENDIF
CPDBSUM RAL 29/3/05 <--
            else
               write(3,511) pdbnam, i, j, htype(i), htype(j), dis,
     +              rpc(1), rpc(2), omega, ih,
     +              iih, ijh 
            endif

c     +                    reslist(1:strlen(reslist))
c     this prints out the pairs of interactions
c     write(3,515) (pairlist(k,1),pairlist(k,2),k=1,ih)
c        next j helix
   20    continue
c     next i helix
   25 continue
      
  456 continue
      close (1)
C      IF(OPTION.ne.2) close (2)
      IF(COPTION.NE.'L'.AND.COPTION.NE.'l') close(2)
c      close (3)
      else if (fistat.eq.opnera)  then
          write(6, 1040) filen(1:strlen(filen))
          finuse =.false.
          fistat = endfi
      else if (fistat.eq.opnerb)  then
          write(6, 1050) brknam(1:strlen(brknam))
      else if (fistat.eq.opnerc)  then
         write(*,*)'Unable to open the 2D structure file ',
     +        sstnam(1:namlen)
      else if (fistat.eq.opnerd.or.fistat.eq.opnere.or.
     +         fistat.eq.opnerf)  then
          write(6, 1070) outfil(1:strlen(outfil))
      end if
      continue
c      close (7)
       close (8)
c      close (9)
      close(brkfi)
      close(strkfi)
c     loop back for next structure
      IF(ENSEMB) THEN
         IF(FISTAT.NE.FINISH) THEN
            FISTAT=OK
            IF(FINUSE) THEN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF

      end
c
c*******************************************************************************
c
      subroutine axb (a, b, cross)
c     implicit   none
c     implicit   complex (a-z)
c
c
      real a(3), b(3), cross(3)
c
      cross(1) = a(2)*b(3) - a(3)*b(2)
      cross(2) = a(3)*b(1) - a(1)*b(3)
      cross(3) = a(1)*b(2) - a(2)*b(1)
      return
      end
c
c*******************************************************************************
c
      subroutine blank (string)
c     implicit none
      character*(*) string
      integer i
c
      do 10 i = 1,len(string)
      string(i:i) = ' '
   10 continue
      return
      end
c
c*******************************************************************************
c
      subroutine checkihlx (ihlx, nstart, alright)
      include 'p_helix.par'
      integer ihlx(2,mhlx), nstart
      logical alright
c
      integer i
c
 1000 format (' helix ',I3,':  start = ',i4,'   stop = ',i4,' !')
c
      alright = .true.
      do 10 i = 1,nstart
      if (ihlx(1,i) .gt. ihlx(2,i) ) then
         write (6, 1000) i, ihlx(1,i), ihlx(2,i)
         alright = .false.
      end if
   10 continue
      return
      end
c
c*******************************************************************************
c
      subroutine checktorsions (phiang, psiang, ihlx, i, htype, allok,
     +                          meth, brkid)
      include   'p_helix.par'
      real       phiang(maxres), psiang(maxres)
      integer    ihlx(2,mhlx), i, meth
      character  htype*1
      logical    allok
      character  brkid*4
c
      real       hphiav, hpsiav, gphiav, gpsiav, minmax
      integer    j, shrink
c
 1000 format ('WARNING: odd torsion in hlx ', i3, ':',
     + ' res,phi,psi,brkid = ',i4,2f7.1,1x,a4)
c
      hphiav = -62.0
      hpsiav = -41.0
      gphiav = -71.0
      gpsiav = -18.0
      minmax =  50.0
      allok  = .true.
      if (meth.eq.2) shrink = 1
      if (meth.eq.3) shrink = 2
c
      do 10  j = ihlx(1,i)+shrink, ihlx(2,i)-shrink
      if (htype.eq.'H')  then
          if (phiang(j).gt.(hphiav+minmax).or.
     +        phiang(j).lt.(hphiav-minmax).or.
     +        psiang(j).gt.(hpsiav+minmax).or.
     +        psiang(j).lt.(hpsiav-minmax))  then
              write (6,1000) i, j, phiang(j), psiang(j), brkid
c              write (7,1000) i, j, phiang(j), psiang(j), brkid
c              write (8,1000) i, j, phiang(j), psiang(j), brkid
c             allok = .false.
          end if
      else if (htype.eq.'G')  then
          if (phiang(j).gt.(gphiav+minmax).or.
     +        phiang(j).lt.(gphiav-minmax).or.
     +        psiang(j).gt.(gpsiav+minmax).or.
     +        psiang(j).lt.(gpsiav-minmax))  then
              write (6,1000) i, j, phiang(j), psiang(j), brkid
c              write (7,1000) i, j, phiang(j), psiang(j), brkid
c              write (8,1000) i, j, phiang(j), psiang(j), brkid
c             allok = .false.
          end if
      end if
   10 continue
      return
      end
c
c*******************************************************************************
c
      subroutine clsaph (hxk, hxj, bk, ek, bj, ej, pck, pcj,
     +                   sk, sj, constrain, ekbk, ejbj)
c     implicit none
c     implicit complex (a-z)
      integer  hxk, hxj
      real     bk(3), ek(3), bj(3), ej(3), pck(3), pcj(3), sk, sj
      logical  constrain
c===============================================================================
c     Subroutine to calculate points of closest approach of two vectors.
c     Input parameters:
c        HXK       : id number of helix K
c        HXJ       : id number of helix J
c        BK()      : coordinates of beginning of helix/vector K
c        EK()      : coordinates of end of helix/vector K
c        BJ()      : coordinates of beginning of helix/vector J
c        EJ()      : coordinates of end of helix/vector J
c        constrain : whether or not to constrain the points of closest contact
c                    to lie within the length of each vector
c     Output parameters:
c        PCK() : coordinates of point on vector K closest to vector J
c        PCJ() : coordinates of point on vector J closest to vector K
c        SK    : scale factor defining position of PCK along K from BK
c        SJ    : scale factor defining position of PCJ along J from BJ
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c     Based on Chothia et al J.Mol.Biol 145,215,1981.
c
c     First version : 31st January 1990
c     Last modified : 29th May 1992
c===============================================================================
      integer i, debug
      real    w1, w2, u11, u22, u12, dotprd, det,
     +        bjbk(3), ekbk(3), ejbj(3), ddstnc
      logical xeqy
c
 2000 format ('HLXK ',i3,' BG: ', 3F8.2,'  END: ', 3F8.2,/,
     +        '     ',3x,'PCK: ', 3F8.2,'   SK: ',  F8.2,/,
     +        'HLXJ ',i3,' BG: ', 3F8.2,'  END: ', 3F8.2,/,
     +        '     ',3x,'PCJ: ', 3F8.2,'   SJ: ',  F8.2,/,
     +        ' distance PCJ-PCK = ',F8.2,/)
 3333 format ('WARNING from subroutine CLSAPH: vectors are parallel')
c
      debug = 0
c
c     Calculate constants
      do 1, i= 1, 3
         bjbk(i) = bk(i) - bj(i)
         ekbk(i) = ek(i) - bk(i)
         ejbj(i) = ej(i) - bj(i)
    1 continue
c
c     Calculate coefficients of scale factors
      w1  = dotprd( bjbk, ekbk )
      w2  = dotprd( bjbk, ejbj )
      u11 = dotprd( ekbk, ekbk )
      u22 = dotprd( ejbj, ejbj )
      u12 = dotprd( ejbj, ekbk )
      det = u11*u22 - u12*u12
c
c     Check for parallel vectors
      if (xeqy(det,0.0e0)) then
          write(6,3333)
      else
c        Calculate scaling factors
         sk = (w2*u12 - w1*u22 ) /det
         sj = (w2*u11 - w1*u12 ) /det
c
c        limit contact points within the length of each vector?
         if (constrain) then
            if (sk .gt. 1.0e0) sk=1.0e0
            if (sj .gt. 1.0e0) sj=1.0e0
            if (sk .lt. 0.0e0) sk=0.0e0
            if (sj .lt. 0.0e0) sj=0.0e0
         end if
c        Calculate points of intersection
         do 2, i= 1, 3
            pck(i) = bk(i) + sk*ekbk(i)
            pcj(i) = bj(i) + sj*ejbj(i)
    2    continue
      end if
c
      if (debug.eq.1) then
         write (9, 2000) hxk, bk, ek, pck, sk,
     +                   hxj, bj, ej, pcj, sj,
     +                   ddstnc(pck,pcj)
      end if
c
      return
      end
C
c*******************************************************************************
c
      real function ddstnc (x1, x2)
c     implicit none
c     implicit complex(a-z)
c
c     Function to calculate the distance between two points
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 21st February 1990
c
      integer      k
      real         dis, x1(3), x2(3)
c
      ddstnc = 0.0e0
      do 1, k= 1, 3
         dis = x2(k) - x1(k)
         ddstnc = ddstnc + dis*dis
    1 continue
      ddstnc = sqrt( ddstnc )
      return
      end
c
c*******************************************************************************
c
      real function dotprd (a, b)
c     implicit none
c     implicit complex (a-z)
c
c     Function to calculate the dot product between two cartesian vectors
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 30th January 1990
c
      integer  i
      real     a(3), b(3)
c
      dotprd = 0.0e0
      do 1, i= 1, 3
         dotprd = dotprd + a(i)*b(i)
    1 continue
      end
c
c*******************************************************************************
c
      real function dstnce (X, I, J)
c
      include 'p_helix.par'
c
c     Function to calculate the distance between two points
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 21st February 1990
c
      integer  i, j, k
      real     dis, x(3,matm)
c
      dstnce = 0.0e0
      if (i .lt. 1) i = j
      if (j .lt. 1) j = i
      if (i .ne. j) then
         do 1 k = 1,3
            dis = x(k,j) - x(k,i)
            dstnce = dstnce + dis*dis
    1    continue
      else
         do 2 k = 1,3
            dstnce = dstnce + x(k,i)*x(k,i)
    2    continue
      end if
      dstnce = sqrt( dstnce )
      return
      end
c
C*******************************************************************************
c
      subroutine eigsrt (d, v, n, np)
c     implicit none
c     implicit complex (a-z)
c
      integer   i, n, j, k, np
      real      d, p, v
      dimension d(np),v(np,np)
c
      do 13 i = 1,n-1
         k = i
         p = d(i)
         do 11 j = i+1,n
            if (d(j).ge.p) then
               k = j
               p = d(j)
            end if
   11    continue
         if (k .ne. i) then
            d(k)=d(i)
            d(i)=p
            do 12 j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
   12       continue
         end if
   13 continue
      return
      end
c
c*******************************************************************************
c
      subroutine findat (atmn, ires, chain, ins, natm, chn1, res1, ins1,
     +                   chn2, res2, ins2, flag, list, nlist)
c
      include 'p_helix.par'
c
      character atmn(matm)*4, chain(matm)*1, ins(matm)*1, chn1*1,
     +          ins1*1, chn2*1, ins2*1
      integer   ires(matm), natm, res1, res2, flag, list(matm), nlist
c===============================================================================
c     Loads into array LIST() the atoms that lie between the two residues RES1
c     and RES2 inclusive. FLAG controls the types of atoms to be considered:
c     1 ==> CA atoms only;
c     2 ==> N,C atoms only;
c     3 ==> N,CA,C atoms only;
c     4 ==> N,CA,C,O atoms only.
c===============================================================================
c
c     subroutine to create list of C-alphas(1), main chain(2) or all (3) in
c     a given range
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 20th February 1990
c-----------------------------------------------------------------------------
      integer istart, istop, getfirstatom, getlastatom, i
c
      nlist = 0
      istart = getfirstatom(chn1, res1, ins1, chain, ires, ins, natm)
      istop  = getlastatom(chn2, res2, ins2, chain, ires, ins, natm)
      do 100 i = istart,istop
      if ( (atmn(i) .eq. ' CA ' .and. flag.ne.2 ) .OR.
     +     (atmn(i) .eq. ' N  ' .and. flag.gt.1 ) .OR.
     +     (atmn(i) .eq. ' C  ' .and. flag.gt.1 ) .OR.
     +     (atmn(i) .eq. ' O  ' .and. flag.gt.3 ) ) then
         nlist = nlist + 1
         list(nlist) = i
      end if
  100 continue
      return
      end
c
c*******************************************************************************
c
      subroutine findhelices (ks, kspl, spah, sp310h, sppih, seqbcd, 
     ~                        chnbrk, seqlen, meth, gap,
     +                        ihlx, nhlx, hch, lb1, lb2, num1, num2,
     +                        htype, chain, ires, ins, atmn, x,
     +                        natm, snd, nss, seq, hseq)
      include 'p_helix.par'
c==============================================================================
c     Last modified : 29th May 1992
c
c     Arguments passed to findhelices are:
c     ks, kspl, spah, sp310h, sppih, chnbrk (as in readst)
c     seqbcd(j) - chain letter//resnum//inscod for residue j ("crrrri")
c     seqlen    - no of residuies defined in SST file
c     meth      - which column of sec. sruc. assignment to use in
c                 defining helices
c     ihlx(2,i) - to be loaded with the residue id numbers of the
c                 residues at the start and end of helix i
c     nhlx      - no of helices found
c     hch(i)    - the chain letter associated with helix(i)
c     lb1(i)/   - the insertion codes of the residues at the
c     lb2(i)      start + end of helix i
c     num1(i)/  - the .pdb defined residue numbers of the
c     num2(i)     residues at ends of helix i
c     htype(i)  - the type of helix i (H or G or I)
c     resn(j)   - residue name of atom j
c     chain(j)  - chain letter of atom j
c     ires(j)   - residue number of atom j
c     ins(j)    - insertion code of atom j
c     atmn(j)   - atom name of atom j
c     x(3,j)    - the XYZ coords of atom j
c     natm      - no of atoms read in from .pdb file
c=============================================================================
      character ks(maxres)*1, kspl(maxres)*1, spah(maxres)*1,
     +          sp310h(maxres)*1, seq(maxres)*1, sppih(maxres)*1
      integer   seqlen, meth, gap, ihlx(2,mhlx), nhlx
      character hch(mhlx)*1, lb1(mhlx)*1, lb2(mhlx)*1, htype(mhlx)*1,
     +          hseq(mhlx)*300
      integer   num1(mhlx), num2(mhlx)
      character seqbcd(maxres)*6, chain(matm)*1,
     +          ins(matm)*1, chnbrk(maxres)*1, atmn(matm)*4, snd*(*)
      integer   ires(matm), natm, nss
      real      x(3,matm)
c
      integer   nhh, ngh, i, j, nstart, nstop, debug, itemp, jtemp,
     +          getfirstatom, khlx, ipos, kstart, kstop, start1, stop1,
     +          start2, stop2, jpos, nsplits, splits(msplits),
     +          gethlxpos, list1(matm), list2(matm), kount1, kount2,
     +          cstart, cstop, rn1, rn2, strlen, turnsz, lenhlx,
     +          nxstart, nxstop,n
CPDBSUM RAL 21/4/05 -->
      INTEGER   NLOST
CPDBSUM RAL 21/4/05 <--
      real      haxis1(3), haxis2(3), vtemp(3), sttemp(3), fntemp(3),
     +          ptemp(3,matm), pctemp(3,matm), angle
      character start*1, stop*1, either*1, sslist*10, curch*1, oldch*1,
     +          hlxpat*(maxrph), ch1*1, ch2*1, ins1*1, ins2*1
      logical   inhelix, alright, chnstrt, atastart
c
c 2020 format ('Hlx ',i3, ' ranges from residue',i5,' (',a3,1x,a1,i4,a1,
c     +        ') to ', i4,' (',a3,1x,a1,i4,a1,')' 2x,a1)
c 2021 format (/,'Helices after splitting!')
c 2100 format (/,'The following helices were split:')
c 2110 format (5x,15i5)
c 2200 format ('Angle between vectors ',3f6.2 ' and ', 3f6.2,
c     + ' is ',f7.2,/)
 2300 format ('Warning: helix ',i4,' begins with "',a,'" !')
 2320 format (i4)
 2330 format ('Warning: hlx ',i3,' knt1=',i3,' knt2=',i3,
     + ' cstrt,cstp=',2i3, ' nstrt,nstp=',2i3)
c 2340 format (' looking at helix ',i4,'  HLXPAT = "',a,'"')
c 2350 format (/,' calculating axes for residues ',i4,' to ',i4,' and ',
c     + i4, ' to ', i4)
c 3000 format ('NSTART = ',i3,'   NSTOP = ',I3)
c 3010 format ('Start residues:', 100i5)
c 3020 format ('Stop  residues:', 100i5)
 3030 format ('ERROR: nstart (',i3,') is not equal to nstop (',i3,')!')
 3100 format (i4)
c
c     Initialise all flags and counters
c
      start   = '>'
      stop    = '<'
      either  = 'X'
      nhh     =  0
      ngh     =  0
      nhlx    =  0
      nstart  =  0
      nstop   =  0
      debug   =  0
      nsplits =  0
      inhelix = .false.
      do 10 i=1,mhlx
         do 9 j=1,6
            hseq(i)((j-1)*50+1:j*50)=
     ~ '                                                  '
 9       continue
 10   continue
CPDBSUM RAL 21/4/05 -->
      NLOST = 0
CPDBSUM RAL 21/4/05 <--
c
c     Begin by finding the starting points of all (easy to find!) helices
c
      do 5000 i = 1,seqlen
c
      if (meth.eq.2) then
c
c        first residue is a special case
         if (i.eq.1) then
            if (ks(i).eq.'H' .or. ks(i).eq.'G'.or.ks(i).eq.'I') then
CPDBSUM RAL 21/4/05 -->
               if (nstart.ge.mhlx) then
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   nstart = nstart + 1
                   ihlx(1,nstart) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 5000
         endif
c        residues immediately after chain breaks are also special cases
         if ( (chnbrk(i).eq.'!' .and. chnbrk(i-1).eq.'!') .OR.
     +        (seqbcd(i)(1:1).ne.seqbcd(i-1)(1:1)) ) then
            if (ks(i).eq.'H' .or. ks(i).eq.'G'.or.ks(i).eq.'I') then
CPDBSUM RAL 21/4/05 -->
               if (nstart.ge.mhlx) then
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   nstart = nstart + 1
                   ihlx(1,nstart) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 5000
         end if
c        all other residues are normal!
         if ( (ks(i).eq.'H'.or.ks(i).eq.'G'.or.ks(i).eq.'I') .AND.
     +         ks(i).ne.ks(i-1) ) then
CPDBSUM RAL 21/4/05 -->
            if (nstart.ge.mhlx) then
               NLOST = NLOST + 1
            ELSE
CPDBSUM RAL 21/4/05 <--
                nstart = nstart + 1
                ihlx(1,nstart) = i
CPDBSUM RAL 21/4/05 -->
            ENDIF
CPDBSUM RAL 21/4/05 <--
         end if
         goto 5000
c        end of method 2 bits
c
      else if (meth.eq.3) then
c
c        first residue is a special case
         if (i.eq.1) then
            if (kspl(i).eq.'H' .or. kspl(i).eq.'G' .or.
     +          kspl(i).eq.'h' .or. kspl(i).eq.'g'.or.
     +          kspl(i).eq.'I'.or.kspl(i).eq.'i') then
CPDBSUM RAL 21/4/05 -->
               if (nstart.ge.mhlx) then
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   nstart = nstart + 1
                   ihlx(1,nstart) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 5000
         endif
c        residues immediately after chain breaks are also special cases
         if ( (chnbrk(i).eq.'!' .and. chnbrk(i-1).eq.'!') .OR.
     +        (seqbcd(i)(1:1).ne.seqbcd(i-1)(1:1)) ) then
            if (kspl(i).eq.'H' .or. kspl(i).eq.'G' .or.
     +          kspl(i).eq.'h' .or. kspl(i).eq.'g'.or.
     +          kspl(i).eq.'I'.or.kspl(i).eq.'i') then
CPDBSUM RAL 21/4/05 -->
               if (nstart.ge.mhlx) then
                   NLOST = NLOST + 1
               ELSE
CPDBSUM RAL 21/4/05 <--
                   nstart = nstart + 1
                   ihlx(1,nstart) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 5000
         end if
c        all other residues are normal!
         if ( (kspl(i).eq.'H' .and.
     +           kspl(i-1).ne.'h' .and. kspl(i-1).ne.'H') .OR.
     +        (kspl(i).eq.'h' .and. kspl(i-1).ne.'H')     .OR.
     +        (kspl(i).eq.'G' .and.
     +           kspl(i-1).ne.'g' .and. kspl(i-1).ne.'G') .OR.
     +        (kspl(i).eq.'g' .and. kspl(i-1).ne.'G') .OR.
     +        (kspl(i).eq.'I'.and.
     +           kspl(i-1).ne.'i'.and.kspl(i-1).ne.'I') .OR.
     +        (kspl(i).eq.'i'.and.kspl(i-1).ne.'I'))   then
CPDBSUM RAL 21/4/05 -->
            if (nstart.ge.mhlx) then
                NLOST = NLOST + 1
            ELSE
CPDBSUM RAL 21/4/05 <--
                nstart = nstart + 1
                ihlx(1,nstart) = i
CPDBSUM RAL 21/4/05 -->
            ENDIF
CPDBSUM RAL 21/4/05 <--
         end if
         goto 5000
c        end of method 3 bits
c
      else
         print*,'Method 1 not yet implemented!'
         return
      endif
 5000 continue
c
c     Next, find the end points of all (easy to find!) helices
c
      do 6000 i = 1,seqlen
c
      if (meth.eq.2) then
c
c        last residue is a special case
         if (i.eq.seqlen) then
            if (ks(i).eq.'H' .or. ks(i).eq.'G'.or.ks(i).eq.'I') then
CPDBSUM RAL 21/4/05 -->
               if (nstop.lt.mhlx) then
CPDBSUM RAL 21/4/05 <--
                   nstop = nstop + 1
                   ihlx(2,nstop) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 6000
         endif
c        residues immediately before chain breaks are also special cases
         if ( (chnbrk(i).eq.'!' .and. chnbrk(i+1).eq.'!') .OR.
     +        (seqbcd(i)(1:1).ne.seqbcd(i+1)(1:1)) ) then
            if (ks(i).eq.'H' .or. ks(i).eq.'G' .or.ks(i).eq.'I') then
CPDBSUM RAL 21/4/05 -->
               if (nstop.lt.mhlx) then
CPDBSUM RAL 21/4/05 <--
                   nstop = nstop + 1
                   ihlx(2,nstop) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 6000
         end if
c        all other residues are normal!
         if ( (ks(i).eq.'H'.or.ks(i).eq.'G'.or.ks(i).eq.'I') .AND.
     +         ks(i).ne.ks(i+1) ) then
CPDBSUM RAL 21/4/05 -->
            if (nstop.lt.mhlx) then
CPDBSUM RAL 21/4/05 <--
                nstop = nstop + 1
                ihlx(2,nstop) = i
CPDBSUM RAL 21/4/05 -->
            ENDIF
CPDBSUM RAL 21/4/05 <--
         end if
         goto 6000
c        end of method 2 bits
c
      else if (meth.eq.3) then
c
c        last residue is a special case
         if (i.eq.seqlen) then
            if (kspl(i).eq.'H' .or. kspl(i).eq.'G' .or.
     +          kspl(i).eq.'h' .or. kspl(i).eq.'g'.or.
     +          kspl(i).eq.'I'.or.kspl(i).eq.'i') then
CPDBSUM RAL 21/4/05 -->
               if (nstop.lt.mhlx) then
CPDBSUM RAL 21/4/05 <--
                   nstop = nstop + 1
                   ihlx(2,nstop) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 6000
         endif
c        residues immediately before chain breaks are also special cases
         if ( (chnbrk(i).eq.'!' .and. chnbrk(i+1).eq.'!') .OR.
     +        (seqbcd(i)(1:1).ne.seqbcd(i+1)(1:1)) ) then
            if (kspl(i).eq.'H' .or. kspl(i).eq.'G' .or.
     +          kspl(i).eq.'h' .or. kspl(i).eq.'g'.or.
     +          kspl(i).eq.'I'.or.kspl(i).eq.'i') then
CPDBSUM RAL 21/4/05 -->
               if (nstop.lt.mhlx) then
CPDBSUM RAL 21/4/05 <--
                   nstop = nstop + 1
                   ihlx(2,nstop) = i
CPDBSUM RAL 21/4/05 -->
               ENDIF
CPDBSUM RAL 21/4/05 <--
            end if
            goto 6000
         end if
c        all other residues are normal!
         if ( (kspl(i).eq.'H' .and.
     +           kspl(i+1).ne.'h' .and. kspl(i+1).ne.'H') .OR.
     +        (kspl(i).eq.'h' .and. kspl(i+1).ne.'H')     .OR.
     +        (kspl(i).eq.'G' .and.
     +           kspl(i+1).ne.'g' .and. kspl(i+1).ne.'G') .OR.
     +        (kspl(i).eq.'g' .and. kspl(i+1).ne.'G').OR.
     +        (kspl(i).eq.'I'.and.
     +           kspl(i+1).ne.'i'.and.kspl(i+1).ne.'I').OR.
     +        (kspl(i).eq.'i'.and.kspl(i+1).ne.'I') )   then
CPDBSUM RAL 21/4/05 -->
            if (nstop.lt.mhlx) then
CPDBSUM RAL 21/4/05 <--
                nstop = nstop + 1
                ihlx(2,nstop) = i
CPDBSUM RAL 21/4/05 -->
            ENDIF
CPDBSUM RAL 21/4/05 <--
         end if
         goto 6000
c        end of method 3 bits

      else
         print*,'Method 1 not yet implemented!'
         return
      endif
 6000 continue
CPDBSUM RAL 21/4/05 -->
C---- If maximum number of helices exceeded, show warning
      IF (NLOST.GT.0) THEN
          PRINT 103, mhlx, NLOST
 103      FORMAT('Warning. Maximum number of helices (',I6,') exceeded',
     -        '. ',I6,' helices lost')
            ENDIF
CPDBSUM RAL 21/4/05 <--
c
c      if (debug.eq.1) then
c         write (8, 3000) nstart, nstop
c         write (8, 3010) (ihlx(1,i), i=1,nstart)
c         write (8, 3020) (ihlx(2,i), i=1,nstop)
c      end if
c
c     Verify that the same no of starts and stops were found! and that each
c     pair are in the right numerical sequence (i.e. stop>start)
c
      if (nstart.eq.nstop) then
         nhlx = nstart
         call checkihlx (ihlx, nhlx, alright)
         if (alright) then
c           load the arrays HTYPE, HCH etc
            do 300 i = 1,nhlx
            hch(i)   = seqbcd(ihlx(1,i))(1:1)
            lb1(i)   = seqbcd(ihlx(1,i))(6:6)
            lb2(i)   = seqbcd(ihlx(2,i))(6:6)
            if (meth.eq.2) htype(i) = ks(ihlx(1,i))
            if (meth.eq.3) htype(i) = kspl(ihlx(1,i))
            if (htype(i).eq.'h') htype(i) = 'H'
            if (htype(i).eq.'g') htype(i) = 'G'
            if (htype(i).eq.'i') htype(i) = 'I'
            read(seqbcd(ihlx(1,i))(2:5), 3100, err=901) num1(i)
            read(seqbcd(ihlx(2,i))(2:5), 3100, err=902) num2(i)
  300       continue
         else
            return
         endif
      else
         write (6, 3030) nstart, nstop
         return
      endif
c
c     Ignore these obsolete bits!
c
c     do 100 i = 1,seqlen
c     if (i .eq. 1)  then
c         if (spah(i) .eq. start)  then
c             nhh = nhh + 1
c             hstart(nhh) = 1
c         end if
c         if (sp310h(i) .eq. start)  then
c             ngh = ngh + 1
c             gstart(ngh) = 1
c         end if
c         goto 100
c     end if
c     if ( (spah(i).eq.start.and.spah(i-1).ne.start.and.
c    +      spah(i-1).ne.either) .OR.
c    +     (spah(i).eq.either.and.spah(i-1).eq.stop) )  then
c         nhh = nhh+1
c         hstart(nhh) = i
c     end if
c     if ( (sp310h(i).eq.start.and.sp310h(i-1).ne.start.and.
c    +      sp310h(i-1).ne.either) .OR.
c    +     (sp310h(i).eq.either.and.sp310h(i-1).eq.stop) )  then
c         ngh = ngh+1
c         gstart(ngh) = i
c     end if
c     if ( (spah(i).eq.stop.and.spah(i+1).ne.stop.and.
c    +      spah(i+1).ne.either) .OR.
c    +     (spah(i).eq.either.and.spah(i+1).eq.stop.or.
c    +      spah(seqlen).eq.start) ) then
c          hstop(nhh) = i
c     end if
c     if ( (sp310h(i).eq.stop.and.sp310h(i+1).ne.stop.and.
c    +      sp310h(i+1).ne.either) .OR.
c    +     (sp310h(i).eq.either.and.sp310h(i+1).eq.stop.or.
c    +      sp310h(seqlen).eq.start) )  then
c          gstop(ngh) = i
c     end if
c 100 continue
c
c     Carry on with some debug printing
c
      if (debug.eq.1) then
         do 200 i = 1,nhlx
         itemp = getfirstatom (hch(i), num1(i), lb1(i), chain, ires,
     +                         ins, natm)
         jtemp = getfirstatom (hch(i), num2(i), lb2(i), chain, ires,
     +                         ins, natm)
c         write (8, 2020) i, ihlx(1,i), resn(itemp), hch(i), num1(i),
c     +                      lb1(i), ihlx(2,i), resn(jtemp), hch(i),
c     +                      num2(i), lb2(i), htype(i)
  200    continue
      end if
c
c     Now check if any helix is actually two or more overlapping helices - this
c     cannot easily be done as a loop (from 1,nhlx) since the value of NHLX
c     may change during execution!
c
      khlx = 0
  150 continue
      khlx = khlx + 1
      if (khlx.gt.nhlx) goto 250
  160 continue
      kstart = ihlx(1,khlx)
      kstop  = ihlx(2,khlx)
      call blank (hlxpat)
      call loadhp (hlxpat, kstart, kstop, spah, sp310h, sppih, 
     ~     htype(khlx))
      lenhlx = strlen(hlxpat)
c      if (debug.eq.1) write (8, 2340) khlx, hlxpat(1:lenhlx)
c
c     ensure first element of HLXPAT is at the origin of an H-bond - if not,
c     skip this helix.
      if (hlxpat(1:1).ne.'>' .and. hlxpat(1:1).ne.'X') then
         write (6, 2300) khlx, hlxpat(1:1)
c         write (8, 2300) khlx, hlxpat(1:1)
         goto 150
      end if
c
c     for meth=3 ensure second element of HLXPAT is also at the origin of
c     an H-bond - if not, skip this helix
      if (meth.eq.3) then
         if (hlxpat(2:2).ne.'>' .and. hlxpat(2:2).ne.'X') then
            write (6, 2300) khlx, hlxpat(1:2)
c            write (8, 2300) khlx, hlxpat(1:2)
            goto 150
         end if
      end if
c
      turnsz = 4
      if (htype(khlx).eq.'G') turnsz = 3
      if (htype(khlx).eq.'I') turnsz = 5
c
c     Setup variables CSTART and CSTOP which point to the first and last
c     residues in the current helix. CSTART points at the 1st residue of the
c     pair of consecutive residues that define the start of the helix (this
c     will initially be position 0 for helices beginning ">?" or "X?");
c     CSTOP points at the 2nd residue of the pair at the end of the
c     helix. I.E. in the helix ">>44<<" CSTART points to "i" and CSTOP points
c     to "n".                   ijklmn
c
      cstart = 0
      if (hlxpat(2:2).eq.'>' .or. hlxpat(2:2).eq.'X') cstart = 1
      cstop  = cstart + turnsz + 1
c
c     Now look for the position of the start of the next mini-helix (i.e. the
c     position of the first pair of consecutive H-bond donors). Start the search
c     at the character/residue immediately following CSTART; terminate the
c     search two characters/residues short of the end of HLXPAT.
      nxstart = cstart
  210 continue
      nxstart = nxstart + 1
      if (nxstart.gt.lenhlx-2) goto 150
      if (atastart(hlxpat,nxstart)) then
         nxstop = nxstart + turnsz + 1
c                                            !this helix is part of previous one
         if ( (nxstop-cstop) .lt. turnsz) then
            cstop = nxstop
            if (cstop.ge.lenhlx) then
               goto 150
            else
               goto 210
            end if
         end if
c                                                          !helices are distinct
         if ( (nxstop-cstop) .gt. turnsz) then
            start1 = kstart
            stop1  = kstart + cstop - 1
            if (meth.eq.2) stop1 = stop1 - 1
            start2 = kstart + nxstart - 1
            if (meth.eq.2) start2 = start2 + 1
            stop2  = kstop
            call splith (ihlx, nhlx, hch, lb1, lb2, num1, num2, htype,
     +                   seqbcd, khlx, stop1, start2, stop2)
            nsplits = nsplits + 1
            if (nsplits.gt.msplits) goto 903
            splits(nsplits) = khlx
            goto 160
         end if
c                                     !helices overlap at their termini and may
c                                     !or may not be distinct depending on the
c                                     !angle between their respective axes.
         if ( (nxstop-cstop) .eq. turnsz) then
            start1 = kstart
            stop1  = kstart + cstop - 1
            if (meth.eq.2) stop1 = stop1 - 1
            start2 = kstart + nxstart - 1
            if (meth.eq.2) start2 = start2 + 1
            stop2  = kstop
            ch1  = seqbcd(start1)(1:1)
            ch2  = seqbcd(stop1)(1:1)
            ins1 = seqbcd(start1)(6:6)
            ins2 = seqbcd(stop1)(6:6)
            read (seqbcd(start1)(2:5), 2320, err=907) rn1
            read (seqbcd(stop1)(2:5), 2320, err=908) rn2
            call findat (atmn, ires, chain, ins, natm, ch1, rn1, ins1,
     +                   ch2, rn2, ins2, gap, list1, kount1)
            ch1  = seqbcd(start2)(1:1)
            ch2  = seqbcd(stop2)(1:1)
            ins1 = seqbcd(start2)(6:6)
            ins2 = seqbcd(stop2)(6:6)
            read (seqbcd(start2)(2:5), 2320, err=907) rn1
            read (seqbcd(stop2)(2:5), 2320, err=908) rn2
            call findat (atmn, ires, chain, ins, natm, ch1, rn1, ins1,
     +                   ch2, rn2, ins2, gap, list2, kount2)
            if (kount1.lt.9.or.kount2.lt.9) then
               write (6, 2330) khlx, kount1, kount2, cstart, cstop,
     +                         nxstart, nxstop
c               write (8, 2330) khlx, kount1, kount2, cstart, cstop,
c     +                         nxstart, nxstop
               goto 150
            end if
c            if (debug.eq.1) write (8, 2350) start1, stop1, start2, stop2
            itemp = 500 + khlx
            call hlxaxs (x, list1, kount1, gap, itemp, vtemp, haxis1,
     +                   sttemp, fntemp, ptemp, pctemp)
            itemp = 600 + khlx
            call hlxaxs (x, list2, kount2, gap, itemp, vtemp, haxis2,
     +                   sttemp, fntemp, ptemp, pctemp)
            call getangle (haxis2(1), haxis2(2), haxis2(3),
     +                     haxis1(1), haxis1(2), haxis1(3),
     +                     angle)
c            write (8, 2200) haxis2, haxis1, angle
c
            if (abs(angle).gt.kinkangle) then
c                                                      !helices are NOT colinear
c                                                      !and need to be split
               call splith (ihlx, nhlx, hch, lb1, lb2, num1, num2,
     +                      htype, seqbcd, khlx, stop1,start2, stop2)
               nsplits = nsplits + 1
               if (nsplits.gt.msplits) goto 903
               splits(nsplits) = khlx
               goto 160
            else
c                                                          !helices ARE colinear
               cstop = nxstop
               if (cstop.ge.lenhlx) then
                  goto 150
               else
                  goto 210
               end if
            end if
         end if
         goto 150
      end if
      goto 210
c
  250 continue
c
      if (debug.eq.1) then
         if (nsplits.gt.0) then
c            write (8, 2100)
c            write (8, 2110) (splits(i), i=1,nsplits)
c            write (8, 2021)
            do 201 i = 1,nhlx
            itemp = getfirstatom (hch(i), num1(i), lb1(i), chain, ires,
     +                            ins, natm)
            jtemp = getfirstatom (hch(i), num2(i), lb2(i), chain, ires,
     +                            ins, natm)
c            write (8, 2020) i, ihlx(1,i), resn(itemp), hch(i), num1(i),
c     +                      lb1(i), ihlx(2,i), resn(jtemp), hch(i),
c     +                      num2(i), lb2(i), htype(i)
  201    continue
         end if
      end if
c
c     Now load the SND character string with the sequence of secondary
c     structure elements as defined in the KS() and KSPL() arrays.
c
      sslist = 'EGHTIeghti'
      nss    = 0
      call blank(snd)
c
c     first residue is a special case
c
      if (meth.eq.2) curch = ks(1)
      if (meth.eq.3) curch = kspl(1)
      ipos = index(sslist, curch)
      if (ipos.ge.1.and.ipos.le.5) then
         nss = nss + 1
         if (nss.gt.len(snd)) goto 904
         snd(nss:nss) = curch
      end if
      if (ipos.ge.6.and.ipos.le.10) then
         nss = nss + 1
         if (nss.gt.len(snd)) goto 904
         snd(nss:nss) = sslist(ipos-5:ipos-5)
      end if
c
c     Now loop over all remaining residues
c
c     oldchn = seqbcd(1)(1:1)
c
      do 7000 i = 2,seqlen
      chnstrt = .false.
c     if (seqbcd(i)(1:1).ne.oldchn) then
c        snd(nss+1:nss+4) = '    '
c        nss = nss + 4
c        oldchn = seqbcd(i)(1:1)
c        chnstrt = .true.
c     end if
c
c     Check if this residue is at the start of a chain. If it is
c     add a buffer of 4 spaces to SND irrespective of the SSSUM.
c
      if (chnbrk(i).eq.'!'.and.chnbrk(i-1).eq.'!') then
         nss = nss + 4
         if (nss.gt.len(snd)) goto 904
         chnstrt = .true.
      end if
c
c     Extract the appropriate SSSUMs for residues i and i-1
c
      if (meth.eq.2) then
         curch = ks(i)
         oldch = ks(i-1)
      end if
      if (meth.eq.3) then
         curch = kspl(i)
         oldch = kspl(i-1)
      end if
c
c     Is residue i of interest?
c
      ipos = index(sslist, curch)
      if (ipos.eq.0) goto 7000
c
c     It is - if it is at the start of a chain, add in the SSSUM to
c     SND irrespective of i-1.
c
      if (chnstrt) then
         if (ipos.ge.1.and.ipos.le.5) then
            nss = nss + 1
            if (nss.gt.len(snd)) goto 904
            snd(nss:nss) = curch
         end if
         if (ipos.ge.6.and.ipos.le.10) then
            nss = nss + 1
            if (nss.gt.len(snd)) goto 904
            snd(nss:nss) = sslist(ipos-5:ipos-5)
         end if
         goto 7000
      end if
c
c     Was residue i-1 of interest as well?
c
      jpos = index(sslist, oldch)
c
c     If i is lower case, then i defines a new piece of secondary
c     structure provided that (i-1) was not the corresponding upper
c     case letter.
c
      if (ipos.ge.6.and.jpos.ne.ipos-5) then
         nss = nss + 1
         if (nss.gt.len(snd)) goto 904
         snd(nss:nss) = sslist(ipos-5:ipos-5)
         goto 7000
      end if
c
c     If i was uppercase and j was different then i defines
c     a new bit again.
c
      if (ipos.le.5.and.jpos.ne.ipos.and.jpos.ne.ipos+5) then
         nss = nss + 1
         if (nss.gt.len(snd)) goto 904
         snd(nss:nss) = sslist(ipos:ipos)
         goto 7000
      endif
c
 7000 continue
c
c     Now duplicate elements of SND that correspond to helices that were
c     have been split
c
      do 7100 i = 1,nsplits
      ipos = gethlxpos(splits(i), snd)
      if (ipos.eq.0) goto 905
      do 7200 j = nss, ipos+1, -1
      snd(j+1:j+1) = snd(j:j)
 7200 continue
      snd(ipos+1:ipos+1) = snd(ipos:ipos)
      nss = nss + 1
      if (nss.gt.len(snd)) goto 906
 7100 continue
c

c      record sequence of helix
      do 8002 i=1,nhlx
         n=1
         do 8001 j=ihlx(1,i),ihlx(2,i)
            hseq(i)(n:n)=seq(j)
            n=n+1
 8001    continue
 8002 continue
      return
c
  901 stop 'ERROR 1 in subroutine FINDHELICES'
  902 stop 'ERROR 2 in subroutine FINDHELICES'
  903 stop 'ERROR 3 in subroutine FINDHELICES'
  904 stop 'too much secondary structure - increase size of SND'
  905 stop 'ERROR 5 in subroutine FINDHELICES'
  906 stop 'ERROR 6 in subroutine FINDHELICES'
  907 stop 'ERROR 7 in subroutine FINDHELICES'
  908 stop 'ERROR 8 in subroutine FINDHELICES'
      end
c
c*******************************************************************************
c
      logical function atastart (hlxpat, nxstart)
      character hlxpat*(*)
      integer   nxstart
c
      character ch1*1, ch2*2
c
      ch1 = hlxpat(nxstart:nxstart)
      ch2 = hlxpat(nxstart+1:nxstart+1)
      if ( (ch1.eq.'>'.or.ch1.eq.'X') .AND.
     +     (ch2.eq.'>'.or.ch2.eq.'X') ) then
         atastart = .true.
      else
         atastart = .false.
      end if
      return
      end
c
c*******************************************************************************
c
      subroutine fndat2 (atomlist, ihelix, lc, nlt)
c
      include 'p_helix.par'
c
      integer   atomlist, ihelix, lc, nlt
      dimension atomlist(mhlx,maxaph), lc(matm)
c
      integer j
c
      nlt = 0
      do 100 j = 1,maxaph
      if (atomlist(ihelix,j).gt.0) then
         nlt = nlt + 1
         lc(nlt) = atomlist(ihelix,j)
      else
         return
      end if
  100 continue
      return
      end
c
c*******************************************************************************
c
      Subroutine GetAngle (Xa, Ya, Za, Xb, Yb, Zb, Ang)
c     Implicit None
      Real       Xa, Ya, Za, Xb, Yb, Zb, Ang
C===============================================================================
C     Filename      : GETANGLE.FOR
C     Last modified : 26th May 1992
C
C     Returns the angle between two vectors (Xa,Ya,Za), and (Xb,Yb,Zb) in
C     degrees.
C
C     Eg: Call GetAngle( XA-XB, YA-YB, ZA-ZB, XC-XB, YC-YB, ZC-ZB, ANG)
C     where (XA,YA,ZA),(XB,YB,ZB),(XC,YC,ZC) are the coordinates
C     of points A,B,C respectively returns ANG=ANGLE A-B-C in degrees.
C===============================================================================
      Real    Cosa, Dx, Dy, Dz, Sina
      Logical Xeqy
C
      If (Xeqy(Xa,Xb) .and. Xeqy(Ya,Yb) .and. Xeqy(Za,Zb)) Then
         Ang = 0.0
         Return
      End If
      If (Xeqy(Xa,0.0) .and. Xeqy(Ya,0.0) .and. Xeqy(Za,0.0)) Then
         Ang = 0.0
         Return
      End If
      If (Xeqy(Xb,0.0) .and. Xeqy(Yb,0.0) .and. Xeqy(Zb,0.0)) Then
         Ang = 0.0
         Return
      End If
      Cosa = Xa*Xb + Ya*Yb + Za*Zb
      Dx   = Ya*Zb - Za*Yb
      Dy   = Za*Xb - Xa*Zb
      Dz   = Xa*Yb - Ya*Xb
      Sina = Sqrt(Dx**2 + Dy**2 + Dz**2)
      Ang  = Atan2(Sina,Cosa)
      Return
      End
c
c*******************************************************************************
c
      integer function getfirstatom (chnlet, resnum, inscode,
     +                               chain,  ires,   ins,      natoms)
c
      include 'p_helix.par'
c
      character chnlet*1, inscode*1, chain(matm)*1, ins(matm)*1
      integer   resnum, ires(matm), natoms
c===============================================================================
c     Returns the position in the arrays CHAIN(), IRES(), and INS() of the
c     first atom belonging to the residue specified by CHNLET, RESNUM, and
c     INSCODE.
c===============================================================================
      integer i
      character temp*6, cnull*6
c
      data cnull /'-0000-'/
c
 1000 format('ERROR in GETFIRSTATOM: cannot find residue "',A6,'"')
 1010 format (A1,I4,A1)
c
      do 100 i = 1,natoms
      if (chnlet .eq. chain(i) .AND.
     +    resnum .eq. ires(i)  .AND.
     +    inscode.eq. ins(i) ) then
         getfirstatom = i
         return
      end if
  100 continue
      write (temp,1010) chnlet, resnum, inscode
      do 120 i = 1,6
      if (temp(i:i).eq.' ') temp(i:i) = cnull(i:i)
  120 continue
      write (6, 1000) temp
      stop
      end
c
c*******************************************************************************
c
      integer function gethlxpos (hlxid, snd)
c     implicit none
      integer hlxid
      character snd*(*)
c
c     returns the position in SND of the HLXIDth occurence of the characters
c     "H" and/or "G" and/or "I"
c
      integer i, kount, strlen
c
      kount = 0
      do 10 i = 1,strlen(snd)
      if (snd(i:i).eq.'H' .or. snd(i:i).eq.'G'.or.
     ~        snd(i:i).eq.'I') then
         kount = kount + 1
         if (kount.eq.hlxid) then
            gethlxpos = i
            return
         end if
      end if
   10 continue
      gethlxpos = 0
      return
      end
c
c*******************************************************************************
c
      integer function getlastatom (chnlet, resnum, inscode,
     +                               chain,  ires,   ins,      natoms)
c
      include 'p_helix.par'
c
      character chnlet*1, inscode*1, chain(matm)*1, ins(matm)*1
      integer   resnum, ires(matm), natoms
c===============================================================================
c     Returns the position in the arrays CHAIN(), IRES(), and INS() of the
c     last atom belonging to the residue specified by CHNLET, RESNUM, and
c     INSCODE.
c===============================================================================
      integer i, getfirstatom, istart
c
      istart = getfirstatom(chnlet, resnum, inscode, chain, ires, ins,
     +                      natoms)
      do 100 i = istart,natoms
      if (chnlet .ne. chain(i) .OR.
     +    resnum .ne. ires(i)  .OR.
     +    inscode.ne. ins(i) ) then
         getlastatom = i-1
         return
      end if
  100 continue
      getlastatom = natoms
      return
      end
c
c*******************************************************************************
c
      subroutine hlxaxs (x, lc, nlt, gap, iii, v, a, st, fn, p, pc)
c
      include 'p_helix.par'
c===============================================================================
c     Subroutine to calculate the helix axis of helix III: it
c     returns the position of the start ST(3) and finish FN(3)
c     coordinates of the helix axis along with the unit vector along the
c     direction of the axis A(3). The calculation is based on the algorithms
c     described in Chothia et al J.Mol.Biol 145,215,1981. The contents of A()
c     are adjusted immediately before returning to ensure that they agree with
c     ST() and FN() (which is not always the case!).
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 30th January 1990
c     Last modified : 28th May 1992
c===============================================================================
      integer  i, j, k, icv, njr, lc(*), nlt, gap, iii, debug
c
      real     x(3,matm), st(3), fn(3), p(3,matm), m(3,3), ap, rcv, ms,
     +         e(3,3), v(3), c(3), dotprd, a(3), sb, se, te, pc(3,matm),
     +         b(3), dnlen, atemp(3),alex(3,matm)
c
c     real     s
c
  100 format( //, 1X, '** Subroutine HLXAXS called **', // )
  110 format( 1X, 'Statistics:-', //,
     +        1X, 'Percent variance along principle axis :- ', F6.2, /,
     +        1X, 'Percent variance along other axis :- ', 2F6.2 )
  120 format (/, ' For helix ',I4,' :',/,
     + ' Unit vector along direction of helix :- ', 3F10.4,/,
     + ' Starting point                       :- ', 3F10.4,/,
     + ' Centre of mass of helix              :- ', 3F10.4,/,
     + ' Finishing point                      :- ', 3F10.4,/,
     + ' Distance between end points          :- ',  F10.4,/,
     + ' Modified unit vector being returned  :- ', 3F10.4,/)
c
      debug = 0
      icv = nlt - 2*gap
      rcv = real( icv )
c
c     Calculate the set of vectors P() orthogonal to the helix axis
      do 22, j= 1, 3
         ap = 0.0e0
         do 1, i= 1, icv
            p(j,i) = x(j,lc(i)) + x(j,lc(i+2*gap)) -
     +               2.e0*x(j,lc(i+1*gap))
            ap = ap + p(j,i)
    1    continue
         ap = ap / rcv
         do 2, i= 1, icv
            p(j,i) = p(j,i) - ap
    2    continue
   22 continue
c     Calculate covariance matrix, M
      do 44, j= 1, 3
         do 4, k= 1, 3
            ms = 0.0e0
            do 3, i= 1, icv
               ms = ms + p(j,i) * p(k,i)
    3       continue
            m(j,k) = ms / rcv
    4    continue
   44 continue
c
c     Diagonalise matrix, store vectors, E, and eigen values, V
      call jacobi( m, 3, 3, v, e, njr )
c
c     Sort by value ie. Helix axes in V(3)
      call eigsrt( v, e, 3, 3 )
c
c     Store unit vector of helix axis
      a(1) = e(1,3)
      a(2) = e(2,3)
      a(3) = e(3,3)
c
c     Calculate total variance, TE
      te = (v(1) + v(2) + v(3)) / 100.0e0
c
c     Calculate centroid of helix, C
      icv = icv + 2*gap
      rcv = real( icv )
      do 6, j= 1, 3
         c(j) = 0.0e0
         do 5, i= 1, icv
            c(j) = c(j) + x(j,lc(i))
    5    continue
         c(j) = c(j) /  rcv
         st(j) = x(j,lc(1)) - c(j)
         fn(j) = x(j,lc(icv)) - c(j)
    6 continue
c
c     Calculate start, B, and end, E, of helix axis
      sb = dotprd( a, st )
      se = dotprd( a, fn )
      do 7, i= 1, 3
         st(i) = c(i) + sb*a(i)
         fn(i) = c(i) + se*a(i)
    7 continue
c
c     Calculate ideal vectors
      do 99, i= 1, icv
         do 8, j= 1, 3
            b(j) = x(j,lc(i)) - c(j)
    8    continue
         sb = dotprd( a, b )
         do 9, j= 1, 3
            pc(j,i) = c(j) + sb*a(j)
    9    continue
   99 continue
c
c     Output statistics
c      write (6,100)
c      write (6,110) v(3)/te, v(1)/te, v(2)/te
c
c     Adjust sign of helix axis if necessary
      dnlen = sqrt( (st(1)-fn(1))**2 + (st(2)-fn(2))**2 +
     +               (st(3)-fn(3))**2)
      do 50 i = 1,3
      atemp(i) = (fn(i)-st(i)) / dnlen
   50 continue
c
c      if (debug .eq.1)  then
c         if (iii.gt.500) write (8,120) iii, a, st, c, fn, dnlen, atemp
c         write (8,120) iii, a, st, c, fn, dnlen, atemp
c         write (9,120) iii, a, st, c, fn, dnlen, atemp
c      end if
c
      do 55 i = 1,3
      a(i) = atemp(i)
   55 continue
c
      return
      end
c
c*******************************************************************************
c
      subroutine iprompt (prompt, i, lowlim, highlim, ifail)
c     implicit   none
      character  prompt*(*)
      integer    i, lowlim, highlim, ifail
C===============================================================================
C     Filename      : iprompt.f                                (Unix version)
C     Last modified : 28th May 1992
C
C     Prompts the user to type in a value for the integer variable I by
C     writing out the string held in PROMPT. The user's response must lie
C     in the range LOWLIM-HIGHLIM.
C     IFAIL = 0   if everything is OK
C           = 1   error converting original value of I to a string
C           = 2   if the combined prompt string is >255 chars
C           = 3   error converting the user's response to an integer
C           = 99  if the user responds with a CONTROL_Z
C===============================================================================
      character fullprompt*255, answer*80, string*80, space*1
      integer   lenp, leni, strlen, ireply
C
CPDBSUM RAL 18/11/13 -->
C 1000 Format ($,A,' ')
 1000 Format (A,' ')
CPDBSUM RAL 18/11/13 <--
 1010 Format (A)
 1020 format ('Value must lie in range ',I5,' to ',i5)
 2001 Format (I1)
 2002 Format (I2)
 2003 Format (I3)
 2004 Format (I4)
 2005 Format (I5)
C
      ifail = 0
      space = ' '
      call blank (fullprompt)
      call blank (string)
c
C     LENP is the length of the original prompt string
      lenp = strlen(prompt)
      if (lenp.eq.0) then
         prompt = 'Enter an integer'
         lenp = strlen(prompt)
      end if
C
C     Convert current value of I to its character equivalent
      if (I.ge.0.and.I.le.9) then
         write (string, 2001) I
      else If (I.ge.10.and.I.le.99) then
         write (string, 2002) I
      else If (I.ge.100.and.I.le.999) then
         write (string, 2003) I
      else If (I.ge.1000.and.I.le.9999) then
         write (string, 2004) I
      else If (I.ge.10000.and.I.le.99999) then
         write (string, 2005) I
      else
         ifail = 1
         return
      end If
C
C     LENI is the length of integer I as a string
      leni = strlen(string)
C
      if (lenp+leni+5.gt.255) then
        ifail = 2
        return
      end if
C
C     Now generate FULLPROMPT
C
      fullprompt = prompt(1:lenp) // space // '<' // string(1:leni)
     +                                     // '>' // space // '-'
      lenp = strlen(fullprompt)
C
C     Now write out FULLPROMPT and read in reply - LENP is now the total
C     length of the combined prompt string.
C
   20 write (6, 1000) fullprompt(1:lenp)
      read  (5, 1010, end=901, err=902) answer
      if (strlen(answer).eq.0) return
      read (answer, *, err=902) ireply
      if (ireply.ge.lowlim .and. ireply.le.highlim) then
         i = ireply
         return
      else
         write (6, 1020) lowlim, highlim
         goto 20
      end if
C
  901 Ifail = 99
      Return
  902 Ifail = 3
      Return
c
      End
c
c*******************************************************************************
c
      subroutine jacobi (a, n, np, d, v, nrot)
c     implicit none
c     implicit complex (a-z)
c
      integer  nmax, np, ip, n, iq, i, j, nrot
c
      real     a, d, v, b, z, h, t, g, tresh, theta, sm, c, s, tau
c
      logical xeqy
c
      parameter (nmax = 100)
      dimension a(np,np), d(np), v(np,np), b(nmax), z(nmax)
c
      do 12 ip = 1,n
         do 11 iq = 1,n
            v(ip,iq) = 0.0
   11    continue
         v(ip,ip) = 1.0
   12 continue
      do 13 ip = 1,n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0
   13 continue
      nrot = 0
      do 24 i = 1,50
         sm = 0.0
         do 15 ip = 1,n-1
            do 14 iq = ip+1,n
               sm = sm + abs(a(ip,iq))
   14       continue
   15    continue
         if (xeqy(sm,0.0)) return
         if (i.lt.4) then
            tresh = 0.2*sm/n**2
         else
            tresh = 0.0
         end if
         do 22 ip = 1,n-1
            do 21 iq = ip+1,n
               g = 100.0*abs(a(ip,iq))
               if (i .gt. 4 .and. xeqy(abs(d(ip))+g,abs(d(ip))) .and.
     +                   xeqy(abs(d(iq))+g,abs(d(iq))) ) then
                  a(ip,iq) = 0.0
               else if (abs(a(ip,iq)).gt.tresh) then
                  h = d(iq)-d(ip)
                  if (xeqy(abs(h)+g,abs(h))) then
                     t = a(ip,iq)/h
                  else
                     theta = 0.5*h/a(ip,iq)
                     t = 1.0/(abs(theta)+sqrt(1.0+theta**2))
                     if (theta.lt.0.0) t = -t
                  end if
                  c = 1.0/sqrt(1.0+t**2)
                  s = t*c
                  tau = s/(1.0+c)
                  h = t*a(ip,iq)
                  z(ip) = z(ip)-h
                  z(iq) = z(iq)+h
                  d(ip) = d(ip)-h
                  d(iq) = d(iq)+h
                  a(ip,iq) = 0.0
                  do 16 j = 1,ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
   16             continue
                  do 17 j = ip+1,iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g-s*(h+g*tau)
                     a(j,iq) = h+s*(g-h*tau)
   17             continue
                  do 18 j = iq+1,n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g-s*(h+g*tau)
                     a(iq,j) = h+s*(g-h*tau)
   18             continue
                  do 19 j = 1,n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g-s*(h+g*tau)
                     v(j,iq) = h+s*(g-h*tau)
   19             continue
                  nrot = nrot+1
               end if
   21       continue
   22    continue
         do 23 ip = 1,n
            b(ip) = b(ip)+z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0
   23    continue
   24 continue
      print*, '50 iterations should never happen'
      return
      end
c
c*******************************************************************************
c
      subroutine loadhat (atmn, ires, chain, natm, lb1, lb2, 
     +     ins, num1, num2, hch, nhlx, atomlist, iflag)
c
      include 'p_helix.par'
c
c      integer       iatm(matm)
      character*4   atmn(matm)
      integer       ires(matm)
      character*1   chain(matm)
      integer       natm
CPDBSUM RAL 29/8/19 -->
C      character*1   lb1(matm), lb2(matm), ins(matm), hch(mhlx)
      character*1   lb1(mhlx), lb2(mhlx), ins(matm), hch(mhlx)
CPDBSUM RAL 29/8/19 <--
      integer       num1(mhlx), num2(mhlx)
      integer       nhlx, atomlist(mhlx,maxaph), iflag
c===============================================================================
c     According to the value of IFLAG, this routine generates lists of either
c     all atoms (IFLAG=1) or only mainchain atoms (IFLAG=2) that belong to the
c     residues of each of the helices. The atoms are loaded into the
c     ATOMLIST() array. Currently, main chain atoms are N, CA, and C.
c===============================================================================
      integer       i, j, kount, istart, istop, getfirstatom,
     +              getlastatom, debug
c
 1000 format ('ERROR in LOADHAT: too many atoms in helix ',I3,':',
     + ' max allowed = ',i4, ' flag = ',i2)
c 1010 format (/,'Main chain atom numbers for helix ',I3,
c     + ' determined by LOADHAT')
c 1020 format (20(16I5,/))
c 1025 format ('Actual atom id numbers are:')
c 1030 format (/,'Entire set of atom numbers for helix ',I3)
c
      debug = 0
c
c     Initialise the ATOMLIST() array
c
      do 50 i = 1,mhlx
      do 60 j = 1,maxaph
      atomlist(i,j) = 0
   60 continue
   50 continue
c
c     Begin loop over each helix
c
      do 100 i = 1,nhlx
      kount = 0
      istart = getfirstatom (hch(i), num1(i), lb1(i), chain, ires,
     +                       ins, natm)
      istop  = getlastatom (hch(i), num2(i), lb2(i), chain, ires,
     +                      ins, natm)
      do 200 j = istart,istop
      if ( iflag.eq.1 .OR.
     +    (iflag.eq.2.and.(atmn(j).eq.' N  '.or.
     +                     atmn(j).eq.' CA '.or.
     +                     atmn(j).eq.' C  ' ) ) ) then
         kount = kount + 1
         if (kount.gt.maxaph) then
            write (6, 1000) i, maxaph, iflag
            stop
         end if
         atomlist(i,kount) = j
      end if
  200 continue
c      if (debug.eq.1) then
c         if (iflag.eq.1) then
c            write (8, 1030) i
c         else
c            write (8, 1010) i
c         end if
c         write (8, 1020) (atomlist(i,k), k=1,kount)
c         write (8, 1025)
c         write (8, 1020) (iatm(atomlist(i,k)), k=1,kount)
c      end if
  100 continue
c
      return
      end
c
c*******************************************************************************
c
      subroutine loadhp (hlxpat, kstart, kstop, spah, sp310h, sppih, 
     ~     htype)
      include 'p_helix.par'
      character hlxpat*(maxrph)
      integer   kstart, kstop
      character spah(maxres)*1, sp310h(maxres)*1, htype*1, 
     ~     sppih(maxres)*1
c
c     Routine to load the HLXPAT string with the structure pattern characters
c     from SPAH or SP310H (depending on HTYPE) for residues KSTART to KSTOP.
c
      integer i, j
c
      j = 0
      do 10 i = kstart,kstop
      j = j + 1
      if (htype.eq.'H') then
         hlxpat(j:j) = spah(i)
      else if(htype.eq.'G') then
         hlxpat(j:j) = sp310h(i)
      else
         hlxpat(j:j)=sppih(i)
      end if
   10 continue
      return
      end
c
c*******************************************************************************
c
      subroutine readbk (x, atmn, resn, chain, iatm, ires, brkfi,
     +                   natm, nres, cres, nch, ins, brkid)
c
      include 'p_helix.par'
c
      integer   iatm(matm), ires(matm), cres(mnch), natm, nres, nch,
     +          brkfi
      real      x(3,matm)
      character atmn(matm)*4, resn(matm)*3, chain(matm)*1,
     +          ins(matm)*1, brkid*4
c===============================================================================
c     subroutine to read brookhaven files
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 20th February 1990
c===============================================================================
      integer   j, oldresnum, strlen, len, debug, kres
      character oldchain*1, oldinscode*1, line*80, space*1,
     +          chnid(mnch)*1
c
  105 format ( a80 )
  110 format ( 6x, 1x, i4, 1x, a4, 1x, a3, 1x, a1, i4, a1, 3x, 3f8.3 )
  120 format ( 62x, a4)
 1000 format ('More than ',I5,' atoms from brkid ',A4)
 1010 format ('More than ',I5,' chains from brkid ',A4)
c 1020 format (A4, ' contains ',I5,' atoms in ',I3,' chain(s)')
c 1030 format ('chain ',I2,':  letter "',a1,'":  nres=',I4)
c 1040 format ('Total no of residues = ',I4)
c
      natm = 0
      nres = 0
      nch  = 0
      debug = 0
      space = ' '
      brkid = '    '
      oldchain = '?'
      oldresnum = -999
      oldinscode = '?'
c
   10 continue
      read (brkfi, 105, err=901, end=100) line
      len = strlen(line)
      if (len.le.5) goto 10
      if (line(1:6).eq.'HEADER' .and. line(10:10).eq.space)
     +                   read (line, 120, err=902) brkid
      if(line(1:6).eq.'ENDMDL') go to 100
      if (line(1:6).eq.'ATOM  ' .or. line(1:6).eq.'HETATM') then
         natm = natm + 1
         if (natm.gt.matm) then
            write (6, 1000) matm, brkid
            stop
         endif
         read (line, 110, err=903) iatm(natm), atmn(natm), resn(natm),
     +                             chain(natm), ires(natm), ins(natm),
     +                             (x(j,natm), j=1,3)
         if (chain(natm).ne.oldchain) then
            nch = nch + 1
            if (nch.gt.mnch) then
               write (6, 1010) mnch, brkid
               stop
            end if
            if (nch.gt.1) cres(nch-1) = nres
            nres = 0
            chnid(nch) = chain(natm)
         end if
         if (chain(natm).ne.oldchain  .or.
     +        ires(natm).ne.oldresnum .or.
     +         ins(natm).ne.oldinscode ) then
            oldchain = chain(natm)
            oldresnum = ires(natm)
            oldinscode = ins(natm)
            nres = nres + 1
         end if
      end if
      goto 10
c
  100 continue
      close (brkfi)
      cres(nch) = nres
      kres = 0
c      if (debug.eq.1) then
c         write (8, 1020) brkid, natm, nch
c         do 150 j = 1,nch
c         write (8, 1030) j, chnid(j), cres(j)
c         kres = kres + cres(j)
c  150    continue
c         write (8, 1040) kres
c      end if
      return
c
  901 print*,'ERROR 1 reading Brookhaven file'
      stop
  902 print*,'ERROR 2 reading Brookhaven file'
      stop
  903 print*,'ERROR 3 reading Brookhaven file'
      stop
      end
c
c*******************************************************************************
c
      subroutine readst (strkfi, nhlx, brkid, phiang, psiang,ks, kspl,
     +                   spah, sp310h, sppih,seqbcd, chnbrk, seqlen, 
     ~                   seq)
c
      include 'p_helix.par'
c
      integer   seqlen
      integer   strkfi, nhlx
      character brkid*4
      character ks(maxres)*1, kspl(maxres)*1,
     +          spah(maxres)*1, sp310h(maxres)*1, seqbcd(maxres)*6,
     +          chnbrk(maxres)*1,seq(maxres)*1,sppih(maxres)*1
      real      phiang(maxres), psiang(maxres)
c===============================================================================
c     subroutine to read structure file
c
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
c
c     First version : 20th February 1990
c     Last modified : 15th May 1992
c
c     4 new character*1 arrays of dimension maxres:
c       ks()     - holds the standard K+S secondary structure assignments
c       kspl()   - holds the modified K+S secondary structure assignments
c       spah()   - holds the column of structure patterns for alpha-helices
c       sp310h() - holds the column of structure patterns for 3,10-helices
c       sppih()  - holds the column of structure patterns for pi helices
c       seqbcd() - holds the character*6 seqnum for each of SEQLEN residues
c       chnbrk() - holds the "!" character if residue i is at a chain break
c
c===============================================================================
      character line*132, brcode*4, chn*1, inscode*1, ss(3)*1 
      integer   i, j, resnum, debug 
      logical   inhelix
c
 1000 format (a)
 1010 format (1x, a4)
 1020 format (19x, i4)
c 1030 format (6x,a1,i4,a1,8x,3(1x,a1),1x,a1,3x,a1,14x,2(1x, f6.1))
c 1030 format (6x,a1,i4,a1,6x,4(1x,a1),1x,a1,3x,a1,a1,13x,2(1x, f6.1))
 1030 format (6x,a1,i4,a1,6x,4(1x,a1),1x,a1,5x,a1,a1,13x,2(1x, f6.1))
 1035 format (5x, a1, a6)
c 1180 format ('residue "',a6,'" has chnbrk "',a1,'"')
c 2000 format ('More than ',I3,' helices for brkid ',A4)
 2010 format ('Mismatched brcodes! ',A4,' and ', a4)
c
      nhlx = 0
      debug = 0
      inhelix = .false.
c
      do 10 i = 1,3
      read (strkfi, 1000, err=901, end=902) line
   10 continue
      read (line, 1010, err=903) brcode
      if (brcode.ne.brkid) then
         write (6, 2010) brkid, brcode
         stop
      endif
      read (strkfi, 1020, err=904, end=905) seqlen
      do 20 i = 1,3
      read (strkfi, 1000, err=906, end=907) line
   20 continue
c
c     now begin reading data proper
c
c     oldchain = '?'
c     lastsschn = '?'
c     oldss(1) = '?'
c     oldss(2) = '?'
c     oldss(3) = '?'
c     oldinscode = '?'
c     oldresnum = -999
c
      do 100 i = 1,seqlen
      read (strkfi, 1000, err=908, end=909) line
      read (line, 1030, err=908, end=909) chn, resnum, inscode,seq(i),
     +      (ss(j),j=1,3), spah(i), sp310h(i), sppih(i), phiang(i), 
     ~      psiang(i)
      read (line, 1035, err=908, end=909) chnbrk(i), seqbcd(i)
c     write (6, 1180) seqbcd(i), chnbrk(i)
      ks(i)   = ss(2)
      kspl(i) = ss(3)
  100 continue
c
c
c     check if the previous residue was at the end of a helix
c
c     if ( (oldss(meth).eq.'H'.or.oldss(meth).eq.'G') .AND.
c    +     (ss(meth).ne.oldss(meth).or.chn.ne.oldchain) ) then
c        ihlx(2,nhlx) = i-1
c        lb2(nhlx) = oldinscode
c        num2(nhlx) = oldresnum
c        inhelix = .false.
c     end if
c
c     now check if the current residue is at the start of a helix
c
c     if ( (ss(meth).eq.'H'.or.ss(meth).eq.'G') .AND.
c    +     (ss(meth).ne.oldss(meth).or.chn.ne.oldchain) ) then
c        nhlx = nhlx + 1
c        if (nhlx.gt.mhlx) then
c           write (6, 2000) mhlx, brkid
c           stop
c        end if
c        hch(nhlx) = chn
c        lb1(nhlx) = inscode
c        num1(nhlx) = resnum
c        htype(nhlx) = ss(meth)
c        ihlx(1,nhlx) = i
c        inhelix = .true.
c     end if
c
c     update flags and read next residue
c
c     oldchain = chn
c     oldresnum = resnum
c     oldinscode = inscode
c     do 120 j = 1,3
c     oldss(j) = ss(j)
c 120 continue
c
c 100 continue
c
c     check whether or not the final residue is in a helix
c
c     if (inhelix) then
c        ihlx(2,nhlx) = seqlen
c        lb2(nhlx) = inscode
c        num2(nhlx) = resnum
c     end if
c
      return
c
  901 print*,'ERROR 1 reading SST file'
      stop
  902 print*,'ERROR 2 reading SST file'
      stop
  903 print*,'ERROR 3 reading SST file'
      stop
  904 print*,'ERROR 4 reading SST file'
      stop
  905 print*,'ERROR 5 reading SST file'
      stop
  906 print*,'ERROR 6 reading SST file'
      stop
  907 print*,'ERROR 7 reading SST file'
      stop
  908 print*,'ERROR 8 reading SST file'
      stop
  909 print*,'ERROR 9 reading SST file'
      stop
      end
c
c*******************************************************************************
c
      subroutine splith (ihlx, nhlx, hch, lb1, lb2, num1, num2, htype,
     +                   seqbcd, khlx, stop1, start2, stop2)
      include 'p_helix.par'
      integer   ihlx(2,mhlx), nhlx, num1(mhlx), num2(mhlx), khlx,
     +          stop1, start2, stop2
      character hch(mhlx)*1, lb1(mhlx)*1, lb2(mhlx)*1, htype(mhlx)*1,
     +          seqbcd(maxres)*6
c
c     On arrival here, helix KHLX has been determined to consist of two (or
c     more!) helices - START1 and STOP1 define the end residues of the 1st
c     helix, START2 and STOP2 the ends of the second helix. (The original KHLX
c     has ends START1 and STOP2). Note that STOP1 and START2 may overlap.
c
      integer i, debug
c
 1000 format (' helix ',i4,' is to be split')
 3100 format (i4)
c
      if (nhlx.ge.mhlx) stop 'ERROR in SPLITH: too many helices!'
c
      debug = 0
      if (debug.eq.1) write (6,1000) khlx
c
c     Begin by shunting all helices of higher number than KHLX one notch
c     down the list.
c
      do 10 i = nhlx,khlx+1,-1
      hch(i+1)    = hch(i)
      lb1(i+1)    = lb1(i)
      lb2(i+1)    = lb2(i)
      num1(i+1)   = num1(i)
      num2(i+1)   = num2(i)
      htype(i+1)  = htype(i)
      ihlx(1,i+1) = ihlx(1,i)
      ihlx(2,i+1) = ihlx(2,i)
   10 continue
c
c     Increase the number of helices by 1
c
      nhlx = nhlx + 1
c
c     Adjust the parameters for the two "new" helices KHLX and KHLX+1
c
      hch(khlx+1)    = hch(khlx)
      htype(khlx+1)  = htype(khlx)
      lb1(khlx+1)    = seqbcd(start2)(6:6)
      lb2(khlx+1)    = lb2(khlx)
      ihlx(1,khlx+1) = start2
      ihlx(2,khlx+1) = stop2
      num2(khlx+1)   = num2(khlx)
      read (seqbcd(start2)(2:5), 3100, err=901) num1(khlx+1)
c
      ihlx(2,khlx)   = stop1
      lb2(khlx)      = seqbcd(stop1)(6:6)
      read (seqbcd(stop1)(2:5), 3100, err=902) num2(khlx)
      return
c
  901 stop 'ERROR 1 in SPLITH'
  902 stop 'ERROR 2 in SPLITH'
c
      end
c
c*******************************************************************************
c
      integer function strlen (string)
c     implicit none
c     implicit complex (a-z)
c
c     find out the real length of the string
c
      character*(*) string
c
      integer nxch
c
      do 200 nxch = len(string), 1, -1
      if (string(nxch:nxch) .ne. ' ') then
         strlen = nxch
         return
      end if
  200 continue
      strlen = 0
      return
      end
c
c*******************************************************************************
c
      subroutine trnang (a, b, c, d, tor, ier)
c     implicit none
c     implicit complex(a-z)
c
c     Subroutine to calculate torsion angle between four points:
C     The four points describing the helix are stored as A, B, C and D
C     The torsion angle is output as TOR
C
C     Written by Tom Flores, Laboratory of Molecular Biology, Department of
C     Crystallography, Birkbeck College, Malet Street, London. WC1E 7HX.
C
c     First version : 31th January 1990
c
      integer      i, ier
c
      real         a(3), b(3), c(3), d(3), tor, ab(3), bc(3), cd(3),
     +             conv, abl, bcl, cdl, rab, rac, rbc, sab, sbc,
     +             costor, sintor, costsq, sign, dotprd
c
      logical xeqy, xney
c
 3333 format ('WARNING from subroutine TRNANG: coinciding points',/,
     +        '    point A is,',3F10.4,/,
     +        '    point B is,',3F10.4,/,
     +        '    point C is,',3F10.4,/,
     +        '    point D is,',3F10.4)
      data         conv/.01745329/
c
c     Calculate vectors and lengths A-B, B-C, C-D
      do 1, i= 1, 3
         ab(i) = b(i) - a(i)
         bc(i) = c(i) - b(i)
         cd(i) = d(i) - c(i)
    1 continue
      ier = 0
      abl = sqrt( dotprd( ab, ab ) )
      bcl = sqrt( dotprd( bc, bc ) )
      cdl = sqrt( dotprd( cd, cd ) )
c     Check for zero length vectors
      if (xeqy(abl,0.0) .or. xeqy(bcl,0.0) .or. xeqy(cdl,0.0)) then
          write (6, 3333)
c          write (8, 3333)
c        if A-B and C-D >0 but B-C=0 calculate dot product and flag
         if (xeqy(bcl,0.0) .and. xney(abl,0.0) .and.
     +                           xney(cdl,0.0)) then
c           if (debug.eq.3) print*,'acos call 4'
            tor = acos( dotprd( ab, cd ) / (abl*cdl) ) /conv
            ier = 1
         else
            ier = -1
         end if
      else
c        OK to proceed
c        Calculate dot products to form cosine
         rab = dotprd( ab, bc ) / (abl*bcl)
         rac = dotprd( ab, cd ) / (abl*cdl)
         rbc = dotprd( bc, cd ) / (bcl*cdl)
c        Form sines
         sab = sqrt(1.e0 - rab*rab)
         sbc = sqrt(1.e0 - rbc*rbc)
         costor = (rab*rbc - rac) / (sab*sbc)
         costsq = costor*costor
         if (costsq .ge. 1.0e0 .and. costor .lt. 0.0e0) tor = 180.0e0
c        If angle is not 180 deg we need to calculate sign using the sine
         if (costsq .lt. 1.0e0) then
            sintor = sqrt(1.0e0 - costsq)
            tor = atan2( sintor, costor )
            tor = tor / conv
c           Produce unit vectors
            do 2, i= 1, 3
               ab(i) = ab(i) / abl
               bc(i) = bc(i) / bcl
               cd(i) = cd(i) / cdl
    2       continue
c           Find determinant
            sign = ab(1)*(bc(2)*cd(3) - bc(3)*cd(2))
            sign = sign + ab(2)*(bc(3)*cd(1) - bc(1)*cd(3))
            sign = sign + ab(3)*(bc(1)*cd(2) - bc(2)*cd(1))
c           Change sign if necessary
            if (sign .lt. 0.0e0) tor = -1.0e0*tor
         end if
      end if
      return
      end
c
c*******************************************************************************
c
      logical  function xeqy (x, y)
c     implicit none
c     implicit complex (a-z)

      real     x, y
c===============================================================================
c     This function returns TRUE if its two arguments (X and Y) are "equal",
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
      logical  function xney (x, y)
c     implicit none
c     implicit complex (a-z)
c
      real     x, y
c===============================================================================
c     This function returns TRUE if its two arguments (X and Y) are "unequal",
c     else it returns "FALSE". How "unequal" X and Y have to be is determined
c     by the parameter ACURCY, currently set to 0.0001.
c===============================================================================
      real       acurcy
      parameter (acurcy = 0.0001)
c
      xney = (abs(x-y) .gt. acurcy)
      return
      end
c
c*******************************************************************************
c
      subroutine countuniqueresidues(chainlist, pairlist, inslist,
     +                               ih, iih, ijh, brkid, reslist)
c     implicit none
      include 'p_helix.par'
      character*1 chainlist(mintp,2), inslist(mintp,2)
      integer     pairlist(mintp,2), ih, iih, ijh
CPDBSUM RAL 21/4/05 -->
C      character   brkid*4, reslist*(12*mintp)
      character   brkid*4, reslist*(24*mintp)
CPDBSUM RAL 21/4/05 <--
c
      integer   i, j, ipos, junk
      character space*1, temp*6, slash*1, catresinfo*6
c
      space = ' '
      slash = '/'
      call blank(reslist)
      ipos = 0
c
c     Begin by concatenating together the AMINOIDs for all pairs of
c     interacting residues - these are loaded into RESLIST a pair at a
c     time (i.e. Ipair1 Jpair1 Ipair2 Jpair2 ....)
c
      do 50 i = 1,ih
         temp = catresinfo(chainlist(i,1), pairlist(i,1), 
     +                       inslist(i,1))
         reslist(ipos+1:ipos+12) = space // brkid // slash // temp 
         ipos = ipos + 12
         temp = catresinfo(chainlist(i,2), pairlist(i,2), 
     +                       inslist(i,2))
         reslist(ipos+1:ipos+12) = space // brkid // slash // temp 
         ipos = ipos + 12
   50 continue
c
      if (ih.eq.1) then
         iih = 1
         ijh = 1
         return
      end if
c
      iih = 1
      do 100 i = 2,ih
      do 150 j = 1,i-1
c     print*,'ih,i,j = ',ih,i,j
      if (chainlist(j,1).eq.chainlist(i,1) .AND.
     +     pairlist(j,1).eq. pairlist(i,1) .AND.
     +      inslist(j,1).eq.  inslist(i,1)     ) goto 99
  150 continue
      iih = iih + 1
   99 continue
      junk = 0
  100 continue
c
      ijh = 1
      do 200 i = 2,ih
      do 250 j = 1,i-1
      if (chainlist(j,2).eq.chainlist(i,2) .AND.
     +     pairlist(j,2).eq. pairlist(i,2) .AND.
     +      inslist(j,2).eq.  inslist(i,2)     ) goto 199
  250 continue
      ijh = ijh + 1
  199 continue
      junk = 0
  200 continue
c
      return
      end
c
c*****************************************************************************
c
      character*6 function catresinfo(chain, resnum, ins)
      character chain*1, ins*1
      integer   resnum
c
      character null*6, space*1
      integer   i
c
 1000 format(i4)
c
      space = ' '
      null  = '-0000-'
c
      catresinfo(1:1) = chain
      catresinfo(6:6) = ins
      write(catresinfo(2:5), 1000) resnum
      do 100 i = 1,6
      if (catresinfo(i:i).eq.space) catresinfo(i:i) = null(i:i)
  100 continue
      return
      end
c
c*****************************************************************************
c
      subroutine newcls(bk, ek, bj, ej, pck, pcj, ekbk, ejbj, phk,
     +                 phj, rpc, dis)
      save
c     implicit none
c     implicit complex(a-z)
      real bk(3), ek(3), bj(3), ej(3), pck(3), pcj(3), phk(3),
     -     phj(3), dist, dis
c      integer hxj,hxk
c     real sk,sj
c===================================================================
c     Subroutine to calculate points of closest approach of two vectors
c     Input parameters:
c     HXK         : id number of helix K
c     HXJ         : id number of helix J
c     BK()        : coordinates of beginning of helix/vector K
c     EK()        : coordinates of end of helix/vector K
c     BJ()        : coordinates of beginning of helix/vector J
c     EJ()        : coordinates of end of helix/vector J
c   Output parameters:
c     PCK()   : coordinates of point on vector K closest to vector J
c     PCJ()   : coordinates of point on vector J closest to vector K
c     SK      : scale factor defining position of PCK along K from BK
c     SJ      : scale factor defining position of PCJ along J from BJ
c     Written by Tom Flores, Laboratory of Molecular Biology, Department of
c     Crystallography, Birkbeck College, Malet Street, London, WC1E 7HX
c     Based on Chothia et al J.Mol.Biol 145,215,1981.
c
c     First version : 31st January 1990
c     Last modified : 29th May 1992
c======================================================
      character*1 rpc(2)
      integer i, count
      real dotprd, mag, pr(3), r(3),
     +     ekbk(3), ejbj(3), ddstnc, distk, distj,
     -     othend(3), othpt(3), point(3), pq(3), vector(3), mag_pq,
     -     hlen, dist1, dist2, dist4, min_dist
c     real dist3
      data count / 0/
c
c     print*, 'In NEWCLS. Helices:', hxk,' [', rpc(1), ']', hxj,
c    - '  [', rpc(2), ']'
c     print*, 'Coords:',pck
c     print*, '      :',pcj
c     print*, 'Ends of helices:', bk
c     print*, '               :', bj
c     print*, 'End-ratios:     ', sk
c     print*, '                ', sj
c     print*, 'Helix vector:   ', ekbk
c     print*, '                ', ejbj
c     print*, 'Helix end-points', phk
c     print*, '                ', phj
c
c     calculate distances between closest approach point and nearest
c     end of appropriate helix
      distk = ddstnc (pck, phk)
      distj = ddstnc (pcj, phj)

c     print*, 'Closest distance:', distk, distj

c---- Check which is the further distance
      hlen = 0.0
      if (distk.gt.distj) then
         do 200, i = 1, 3
            point(i) = phk(i)
            vector(i) = ejbj(i)
            othpt(i) = bj(i)
            othend(i) = ej(i)
            pq(i) = othpt(i) - point(i)
            hlen = hlen + (bj(i)-ej(i)) * (bj(i)-ej(i))
 200     continue
      else
         do 300, i = 1, 3
            point(i) = phj(i)
            vector(i) = ekbk(i)
            othpt(i) = bk(i)
            othend(i) = ek(i)
            pq(i) = othpt(i) - point(i)
            hlen = hlen + (bk(i)-ek(i)) * (bk(i) - ek(i))
 300     continue
      endif

C---- Compute the helix length
      if (hlen.ge.0.0) hlen = sqrt(hlen)

C---- Normalise the vector
      mag = vector(1) * vector(1) + vector(2) * vector(2)
     -     + vector(3) * vector(3)
      if (mag.gt.0.0) then
         mag = sqrt (mag)
         vector(1) = vector(1) /mag
         vector(2) = vector(2) /mag
         vector(3) = vector(3) /mag
      endif

C---- Calculate the distance between the point and the vector
      call axb (pq, vector, pr)

C---- Calculate the required distance
      mag = pr(1) * pr(1) + pr(2) * pr(2) + pr(3) * pr(3)
      if (mag.gt.0.0) then
         dist=sqrt (mag)
      endif

c     print*, ' DISTANCE', dist

C---- Calculate point of intersection on the helix
      mag_pq = -dotprd (vector,pq)
      do 400, i=1, 3
         r(i) = othpt(i) + mag_pq * vector(i)
 400  continue

      count = count + 1
c      write(9,410) count, 1, (point(i), i =1, 3)
c      write(9,410) count, 2, (othpt(i), i = 1, 3)
c      write(9,410) count, 3 (othend(i), i = 1, 3)
c      write(9,410) count, 4, (r(i), i=1, 3)
c 410  format('ATOM    ',I3,'  C   ALA I',I4,'    ',3F8.3)

c     Find whether the intersection point falls within the helix ends
      dist1 = ddstnc ( r, othpt)
      dist2 = ddstnc (r, othend)
c     The 0.1 value in the next line could change to 0.01 or 0.2
      if (dist1 + dist2.gt.hlen + 0.10) then
c       print*, count, '. ** Intersection outside helix'
c     print*, 'Helices:', hxk,' [', rpc(1),']',hxj,
c    -     ' [', rpc(2),']'
c     print*, 'difference:', hlen -dist1 -dist2
c     print*, ' DISTANCE', dist
         min_dist = ddstnc (bk, bj)
         rpc(1) = 'n'
         rpc(2) = 'n'
         dist4 = ddstnc (ek, ej)
         
         if (dist4 .lt. min_dist) then
            min_dist = dist4
            rpc(1) = 'c'
            rpc(2) = 'c'
         end if

         dist4 = ddstnc (bk, ej)
         if (dist4 .lt. min_dist) then
            min_dist = dist4
            rpc(1) = 'n'
            rpc(2) = 'c'
         end if
         
         dist4 = ddstnc (ek, bj)
         if (dist4 .lt. min_dist) then
            min_dist = dist4
            rpc(1) = 'c'
            rpc(2) = 'n'
         end if

         dis = min_dist
c     print*, 'dis 1:', dis
c     dist3 = ddstnc (bk, bj)
c     dist4 = ddstnc (ek, ej)
c     print*, ' DISTANCE', dist3, dist4
c     dist3 = ddstnc (bk, ej)
c     dist4 = ddstnc (ek, bj)
c     print*, ' DISTANCE', dist3, dist4
c     do 110, i = 1, 3
c     pck(i) = phk(i)
c     pcj(i) = phj(i)
c110  continue
c     dis = ddstnc(pck, pcj)
      else
         dis = dist
      endif
      return
      end
C
c**************************************************************

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

c**************************************************************************
CPDBSUM RAL 29/3/05 -->
C**************************************************************************
C
C  SUBROUTINE HELHED  -  Write header records to helices.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE HELHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 15)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN', 'HELIX_NO',
     -     'RESNUM1', 'RESNAM1',
     -     'RESNUM2', 'RESNAM2',
     -     'TYPE', 'NRESID', 'SST_CONTEXT', 'RUN_LENGTH', 'UNIT_RISE',
     -     'RES_PER_TURN', 'PITCH', 'RMSD_FROM_IDEAL', 'SEQUENCE' /

C---- Define field id tag
      IDTAG = 'HL_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:HELICES'
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
C  SUBROUTINE HINHED  -  Write header records to helint.promotif file
C                        giving database.dat type field names
C
C----------------------------------------------------------------------+---

      SUBROUTINE HINHED(IUNIT)

      INTEGER       NFIELD
      PARAMETER    (NFIELD = 14)

      CHARACTER*3   IDTAG
      CHARACTER*20  FLDNAM(NFIELD)

      INTEGER       IFIELD, LENSTR, IUNIT

      DATA FLDNAM / 'CHAIN',
     -     'CHAIN1', 'HELIX_NUM1', 'HELIX_TYPE1',
     -     'CHAIN2', 'HELIX_NUM2', 'HELIX_TYPE2',
     -     'DISTANCE', 'INT_TYPE1', 'INT_TYPE2', 'OMEGA',
     -     'NINT_PAIRS', 'NRES_INT1', 'NRES_INT2' /

C---- Define field id tag
      IDTAG = 'HI_'

C---- Write record indicating start of these data
      WRITE(IUNIT,10) 'START:HELIX_INTERACTIONS'
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
C  SUBROUTINE GETHNO  -  Get the helix number according to its location in
C                        its chain
C
C----------------------------------------------------------------------+---

      SUBROUTINE GETHNO(CHAIN,HCH,MHLX,HPOS,HLXNO)

      INTEGER       MHLX

      CHARACTER*1   CHAIN, HCH(MHLX)

      INTEGER       HLXNO, HPOS, I

C---- Initialise helix number
      HLXNO = 0

C---- Determine this helix's number in its chain
      DO 600, I = 1, HPOS
          IF (CHAIN.EQ.HCH(I) .OR.
     -        (CHAIN.EQ.' ' .AND. HCH(I).EQ.'-')) THEN
              HLXNO = HLXNO + 1
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
CPDBSUM RAL 29/3/05 <--
