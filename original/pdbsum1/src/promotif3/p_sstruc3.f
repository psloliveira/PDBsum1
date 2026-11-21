      PROGRAM sstruc
C***********************************************************************
C*      NAME: p_sstruc_ensem.f                                         *
C*  FUNCTION: Program to identify protein secondary structure          * 
C* COPYRIGHT: (1989) David Keith Smith                                 *
C*---------------------------------------------------------------------*
C*    AUTHOR: D. K. Smith                                              *
C*      DATE:  1989                                                        *
C*---------------------------------------------------------------------*
C*    INPUTS: NONE                                                     *
C*   OUTPUTS: NONE                                                     *
C*    LOCALS: NONE                                                     *
C*     CALLS: NONE                                                     *
C*---------------------------------------------------------------------*
C* MODIFICATION RECORD                                                 *
C* DATE     INITS   COMMENTS                                           *
c  11/03/96  GH     to avoid divide by zero in Subroutine mktb for 1mco
c                      don't know why
c  27/03/96  GH     reorganize multiple/single file handling
c  27/07/96  GH     add HYP (hydroxyproline) to hetcod list
c  04/09/96  GH     increase maxchn to 1000 in order to handle files with 
c                   1 chain complete and 1 chain Calpha only
c  10/02/98  GH     add extra modified nucleic acids to dudcod (gpn etc)
c  21/11/00  GH     modified to change output format for sheet from
c                   single character to a 3 digit number 
c  28/11/00  GH     small modification to correct running for peptidenucleic
c                   acid files e.g. 1pdt and 1pup
c  20/03/01  GH     small modifications to correct for 'x' formatting
c                   ierror changed from integer to logical
c  23/09/02  GH     add DIV (disovaline) to hetcod list
c****************************************************************************
c global variables c
      INCLUDE 'p_sstruc.par'

      character aacod1*(abetsz),aacod3*(3*abetsz),atname*((4*mnchna)-4),
     +          ssbond(nstruc,maxres)*1, sssumm(maxres)*1,
     +          summpl(maxres)*1, seqbcd(maxres)*6, 
     +          autsum(maxres)*1, ssord*21, brksym(maxres)*1,
     +          headln*132, altcod(maxres)*1, aastd(maxres)*1,
     +          atpres(3,maxres)*1, thmprs(maxres)*1, hetcod*102,
     +          outsty*1, dsdef(maxres)*1, dudcod*117, modcod*(nmod*3),
     +          stdmod*(nmod+1), nstdmd*(nmod+1), modpl*(nmod),
     +          batnam*((4*mnchna)-4), tphnam*((4*mnchna)-4),
     +          pphnam*((4*mnchna)-4), olenam*((4*mnchna)-4),
     +          lolnam*((4*mnchna)-4)
      integer   seqcod(maxres), seqlen, hbond(maxbnd,maxres),
     +          chnsz(maxchn), chnct, bridpt(2,maxres), hedlen,
     +          dspart(maxres), chnend(maxchn), oois(2,maxres),i,
     +          shtcod(maxres)
      real      mcaaat(ncoord,mnchna,maxres), hbonde(maxbnd,maxres),
     +          scaaat(ncoord,sdchna,maxres), mcangs(nmnang,maxres),
     +          scangs(nsdang,maxres), cadsts(11,maxres), 
     +          dsdsts(3,maxres)
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres), caonly, 
     +          readok, astruc
c
      data aacod1(1:29)  /'ABCDEFGHIXKLMNXPQRSTXVWXYZX--'/
      data aacod3( 1:39) /'ALAASXCYSASPGLUPHEGLYHISILEXXXLYSLEUMET'/
      data aacod3(40:78) /'ASNXXXPROGLNARGSERTHRXXXVALTRPXXXTYRGLX'/
      data aacod3(79:87) /'UNKPCAINI'/
      data hetcod(1:39)  /'AIBPHLSECALMMPRFRDLYMGLMPPHPGLOLETPHABA'/
      data hetcod(40:78) /'NLEB2VB2IB1FBNOB2AB2FIVALOVSTAPVLCALPHA'/
      data hetcod(79:102) /'BNADCIAHSCHSMSELOLHYPDIV'/
      data modcod(1:39)  /'ACEANIBZOCBZCLTFORNH2PHORHATFATOS NHMYR'/
      data modcod(40:48) /'BOCMSUPYR'/
      data dudcod( 1:39) /'  A  C  G  T  U1MA5MCOMC1MG2MGM2G7MGOMG'/
      data dudcod(40:78) /' YG  I +UH2U5MUPSU +C +AGALAGLGCUASGMAN'/
      data dudcod(79:117) /'GLCCEGG4SAGSNAGGLSNGSNAMEXCGPNCPNTPNAPN'/
      data stdmod /'ABCDEFGHIJKLMNZY'/
      data nstdmd /'abcdefghijklmnzy'/
      data modpl  /'+0++++0++++0+++'/
      data ssord /'HEheBbGgIiJKLMNOPTtS '/
      data atname / ' N   CA  C   O      ' /
      data tphnam / ' N   CA  C   O1  O2 ' /
      data pphnam / ' N   CA  P   OP1 OP2' /
      data olenam / ' ON  CA  C   O      ' /
      data batnam / ' N   CA  B   O1  O2 ' /
      data lolnam / ' N   CA  CH  OH     ' /
      data atomin /allatm * .false./, mcaaat /allcrd * 0/
c
c     local variables
c
      integer   finish, opnera, opnerb, opnerc, fiend, ok, endnfi
      character indirn*1
      parameter (finish = 2, opnera = 1, opnerb = 3, fiend = 4, ok = 0, 
     +           opnerc = 5, endnfi = 6, indirn='@') 
      logical   finuse,ENSEMB,ierror
      integer   fistat, namlen, iocode, corend,istart,iend,nps
      data fistat/0/, finuse/.false./
      character*(fnamln) finam, brknam, outnam,pdbcod
c
c     input file control loop
c
      Outsty = 'F'
      fistat = ok
c
      read(*,'(A)')brknam
      IF(BRKNAM(1:1).EQ.'%') THEN
         ENSEMB=.TRUE.
      ELSE
         ENSEMB=.FALSE.
      ENDIF
 
c     if ensemble remove % from file name
      IF(ENSEMB) THEN
         FINAM=BRKNAM
         FINUSE=.TRUE.
         namlen = index(brknam,space) - 1
         DO 180 I=1,NAMLEN-1
            FINAM(I:I)=FINAM(I+1:I+1)
 180     CONTINUE
         FINAM(NAMLEN:NAMLEN)=SPACE
         OPEN(NAMFI,FILE=FINAM,STATUS='OLD',IOSTAT=IOCODE)
c      ENDIF
      if (namlen .le. 0) then
         fistat = finish
      end if
      endif
c
c     open the brookhaven file if everything ok
c
 200  CONTINUE
      IF(FISTAT.EQ.OK) THEN
         IF(ENSEMB) THEN
            READ(NAMFI,'(A)',IOSTAT=IOCODE) BRKNAM
            namlen=index(brknam,space)-1
            IF(IOCODE.LT.0) THEN
               CLOSE(NAMFI)
               FINUSE=.FALSE.
               FISTAT=ENDNFI
            ELSE
               WRITE(*,*)BRKNAM(1:namlen)
            ENDIF
c         ELSE
c            BRKNAM=FINAM
c            write(*,*) 'brknam is ',brknam
         ENDIF
      ENDIF
         
      If (fistat .eq. ok) then
         open (Unit=brkfi, file=brknam, status='old', iostat=iocode)
	 if (iocode .ne. 0) fistat = opnerb
      end if
      if(fistat.eq.ok) then
         call getnam(brknam,istart,iend,ierror)
         pdbcod=brknam(istart:iend)
         nps=index(pdbcod,' ')
         outnam=pdbcod(1:nps-1)//'.sst'
      end if
c
c     if survived this far have a file to process otherwise
c     provide error messages
c
      if (fistat .eq. ok) then
c        **Process file
         open(unit = errfi, status = 'SCRATCH')
         readok = .true.
         caonly = .true.
         astruc = .false.
         call rdbrk(mcaaat, atname, atomin, aacod3, hetcod, seqcod,
     +        seqbcd, scaaat, scatin, altcod, aastd,  atpres,
     +        thmprs, caonly, astruc, headln, hedlen, readok,
     +        dudcod, modcod, stdmod, nstdmd, modpl,  seqlen,
     +        batnam, tphnam, pphnam, olenam, lolnam)
         if (readok .and. .not. caonly) then
            open(unit=outfi,file=outnam,status='unknown',iostat=iocode)
            if (iocode .ne. 0) fistat = opnerc
c           ** look for breaks accross the peptide bonds
            call chnbrk(mcaaat, atomin, seqbcd, brksym, chnsz, chnct,
     +           chnend, caonly, seqlen)
C     6/3/96 TO STOP CRASHING FOR FILES WITH LOTS OF CALPHA-ONLY RESIDUES
            IF(CHNCT.GT.MAXCHN) then
               close(unit=outfi,status='DELETE')
               GO TO 1500
            endif
c           ** determine the calpha distances i+2 to i+11
            if (outsty .eq. 'F')
     +           call mkcadt(mcaaat, atomin, cadsts, seqbcd, seqlen)
c           ** create hydrogen atoms on the main chain nitrogens
            call creatH(mcaaat, atomin, chnsz, chnct)
c           ** determine all hydrogen bonds
            call mkhbnd(mcaaat, atomin, hbond, hbonde, seqcod, chnend, 
     +           seqlen)
c           ** calculate OOI numbers
            call ooi(atomin, mcaaat, seqbcd, oois, seqlen)
c           ** find main chain torsion angles
            call mkangl(mcaaat, mcangs, atomin, chnsz, chnct, caonly)
c           ** find sidechain chi angles
            if (outsty .eq. 'F')
     +           call mksang(mcaaat, atomin, scaaat, scangs, scatin, 
     +           seqcod, aacod3, hetcod, seqlen)
c           ** determine turns, bridges and sheets
            call mktb(hbond, ssbond, mcangs, bridpt, chnend, seqlen,
     ~           shtcod)
c           ** determine bend residues
            call mkbend(mcangs, ssbond, seqlen)
c           ** make Kabsch and Sander type 2ndry structure summary
            call mksumm(ssbond, sssumm, summpl, seqlen)
c           ** make the authors atructure
            call mkauth(seqbcd, autsum, ssord, astruc, seqlen)
c           ** find all disulphide bonds
            if (outsty .eq. 'F')
     +           call disulf(seqbcd, dsdsts, scangs, scaaat, scatin, 
     +           atomin, mcaaat, seqcod, dsdef, dspart, 
     +           astruc, seqlen)
c           **  produce the output file
            call result(hbonde, hbond,  aacod3, aacod1, seqcod, seqbcd, 
     +           brksym, ssbond, bridpt, sssumm, summpl, autsum, 
     +           mcangs, altcod, aastd,  atpres, thmprs, scangs,
     +           cadsts, brknam, headln, hedlen, outsty, dsdsts,
     +           dsdef,  dspart, hetcod, oois,   seqlen, shtcod)
            close(unit=outfi)
         else
c           close(unit=outfi)
            if (readok .and. caonly .and. outsty .eq. 'F') then
c           write(*,*)'Warning: this is a Calpha-only file'
               corend=index(brknam,'.')
               if(corend.ge.5) then
                  outnam=brknam(corend-4:corend)//'sca'
               else
                  outnam=brknam(1:corend)//'sca'
               endif
               open(unit=outfi,file=outnam,status='unknown',
     +              iostat=iocode)
               if (iocode. ne. 0) then
                  write(*,*)'Unable to open output file ',
     +                 outnam
		  close(unit=errfi)
               else
                  call chnbrk(mcaaat, atomin, seqbcd, brksym, chnsz, 
     +                 chnct,  chnend, caonly, seqlen)
                  call mkcadt(mcaaat, atomin, cadsts, seqbcd, seqlen)
                  call ooi(atomin, mcaaat, seqbcd, oois, seqlen)
                  call mkangl(mcaaat, mcangs, atomin, chnsz, chnct, 
     +                 caonly)
                  call mkauth(seqbcd, autsum, ssord, astruc, seqlen)
                  call carslt(aacod3, aacod1, seqcod, seqbcd, mcangs, 
     +                 cadsts, brksym, altcod, aastd,  brknam, 
     +                 headln, hedlen, autsum, hetcod, oois, 
     +                 seqlen)
               end if
            else
               close(unit=errfi)
            end if
            if (astruc) close(unit=authfi)
            if (caonly .and. readok) write(*,*)'C-alpha only file'
            if (.not. readok) write(*,*)
     +           'Error while reading file: No amino acids found'
         end if
 1500    CONTINUE
      else
	 if (fistat .eq. opnera) then
            write(*,*)'Unable to open the file of file names ',
     +           finam(1:namlen)
            finuse = .false.
         else          
            if (fistat .eq. opnerb) then
               write(*,*)'Unable to open the Brookhaven format file ',
     +              brknam(1:namlen)
            else
               if (fistat .eq. opnerc) then
                  write(*,*)'Unable to open output file ',
     +                 outnam
               end if
            end if
         end if
      end if
      IF(ENSEMB) THEN
         if(fistat.ne.finish) then
            fistat=ok
            if(finuse) then
               go to 200
            endif
         endif
      ENDIF
c     
      END
c****************************************************************************
      SUBROUTINE RDBRK( mcaaat, atname, atomin, aacod3, hetcod, seqcod,
     +                  seqbcd, scaaat, scatin, altcod, aastd, atpres,
     +                  thmprs, caonly, astruc, headln, hedlen, readok,
     +                  dudcod, modcod, stdmod, nstdmd, modpl, seqlen,
     +			batnam, tphnam, pphnam, olenam,lolnam)

      INCLUDE 'p_sstruc.par'

      Integer seqlen, seqcod(maxres), hedlen
      Character aacod3*(abetsz*3), atname*((mnchna*4)-4),
     +          seqbcd(maxres)*6,
     +          headln*132, altcod(maxres)*1, aastd(maxres)*1,
     +          atpres(3,maxres)*1, thmprs(maxres)*1, hetcod*102,
     +          modcod*(nmod*3), stdmod*(nmod+1), nstdmd*(nmod+1),
     +          modpl*(nmod), dudcod*117, batnam*((mnchna*4)-4),
     +          tphnam*((4*mnchna)-4), pphnam*((4*mnchna)-4),
     +          olenam*((4*mnchna)-4), lolnam*((4*mnchna)-4)
      Real mcaaat(ncoord,mnchna,maxres), scaaat(ncoord,sdchna,maxres)
      Logical Atomin(mnchna,maxres), scatin(sdchna,maxres), Caonly,
     +        astruc, readok
      Integer aaconv
c
c     local variables
c
      integer   aapln
      character allspc*56
      parameter (aapln = 13, allspc =
     +  '                                                        ')
      real      acoord(ncoord), oldocc, occup, therm, oldthm
      logical   altset, xney
      integer   atcnt,attype, I, aacnt, wkpl,
     +          clslen, cmplen, soulen, sctype, altpos(sideaa), respl,
     +          scatct, thmcnt, aacode, nscats(sideaa), hedend, sqrsct
      character aain(maxres)*3, atmnam*4, resnam*3, altloc*1,
     +          brkrec*80, oldchn*1, seqnum*6,
     +          innum*6, keywd*(keylen), oldatm*4, setloc*1,
     +          class*40, compnd*50, source*50, brcode*4, cbet*4,
     +          scname(sideaa)*56, altnam(3)*56, contch*1, incode*4,
     +          wkstr*50, inseqr*150, setatm*4
      integer   strlen, resfnd
      data nscats / 1, 4, 2, 4, 5, 7, 0, 6, 4,-1, 5, 4, 4, 4,-1, 
     +              3, 5, 7, 2, 3,-1, 3,10,-1, 8, 5,-1, 4,12, 3,
     +              7, 4, 2, 2, 7, 6, 2, 7, 4, 5, 9, 2, 4, 3, 4,
     +              7, 4, 1, 7, 3,11, 7, 1,15, 7,-1, 4,15,10, 4, 
     +              4, 3, 3  /
      data scname(1) /allspc/
      data scname(2)( 1:40) /' CG  AD1 AD2                            '/
      data scname(2)(41:56) /'                '/
      data scname(3)( 1:40) /' SG                                     '/
      data scname(3)(41:56) /'                '/
      data scname(4)( 1:40) /' CG  OD1 OD2                            '/
      data scname(4)(41:56) /'                '/
      data scname(5)( 1:40) /' CG  CD  OE1 OE2                        '/
      data scname(5)(41:56) /'                '/
      data scname(6)( 1:40) /' CG  CD1 CD2 CE1 CE2 CZ                 '/
      data scname(6)(41:56) /'                '/
      data scname(7) /allspc/
      data scname(8)( 1:40) /' CG  ND1 CD2 CE1 NE2                    '/
      data scname(8)(41:56) /'                '/
      data scname(9)( 1:40) /' CG1 CD1                                '/
      data scname(9)(41:56) /'                '/
      data scname(10)/allspc/
      data scname(11)( 1:40)/' CG  CD  CE  NZ                         '/
      data scname(11)(41:56)/'                '/
      data scname(12)( 1:40)/' CG  CD1 CD2                            '/
      data scname(12)(41:56)/'                '/
      data scname(13)( 1:40)/' CG  SD  CE                             '/
      data scname(13)(41:56)/'                '/
      data scname(14)( 1:40)/' CG  OD1 ND2                            '/
      data scname(14)(41:56)/'                '/
      data scname(15)/allspc/ 
      data scname(16)( 1:40)/' CG                                     '/
      data scname(16)(41:56)/'                '/
      data scname(17)( 1:40)/' CG  CD  OE1 NE2                        '/
      data scname(17)(41:56)/'                '/
      data scname(18)( 1:40)/' CG  CD  NE  CZ                         '/
      data scname(18)(41:56)/'                '/
      data scname(19)( 1:40)/' OG                                     '/
      data scname(19)(41:56)/'                '/
      data scname(20)( 1:40)/' OG1 CG2                                '/
      data scname(20)(41:56)/'                '/
      data scname(21)/allspc/
      data scname(22)( 1:40)/' CG1 CG2                                '/
      data scname(22)(41:56)/'                '/
      data scname(23)( 1:40)/' CG  CD1 CD2                            '/
      data scname(23)(41:56)/'                '/
      data scname(24)/allspc/
      data scname(25)( 1:40)/' CG  CD1 CD2                            '/
      data scname(25)(41:56)/'                '/
      data scname(26)/allspc/
      data scname(27)/allspc/
      data scname(28)( 1:40)/' CG  CD  OE                             '/
      data scname(28)(41:56)/'                '/
      data scname(29)( 1:40)/' CG  CD  CE  NZ                         '/
      data scname(29)(41:56)/'                '/
      data scname(30)( 1:40)/' CB1 CB2                                '/
      data scname(30)(41:56)/'                '/
      data scname(31)( 1:40)/' CG  CD1 CD2                            '/
      data scname(31)(41:56)/'                '/
      data scname(32)( 1:40)/'SEG  OD1 OD2                            '/
      data scname(32)(41:56)/'                '/
      data scname(33)( 1:40)/' CM                                     '/
      data scname(33)(41:56)/'                '/
      data scname(34)( 1:40)/' SG                                     '/
      data scname(34)(41:56)/'                '/
      data scname(35)( 1:40)/' CG  CD1 CD2                            '/
      data scname(35)(41:56)/'                '/
      data scname(36)/allspc/
      data scname(37)( 1:40)/' CG  CD  CE  NZ                         '/
      data scname(37)(41:56)/'                '/
      data scname(38)( 1:40)/' CM                                     '/
      data scname(38)(41:56)/'                '/
      data scname(39)( 1:40)/' CG  CD1 CD2                            '/
      data scname(39)(41:56)/'                '/
      data scname(40)/allspc/
      data scname(41)( 1:40)/' CG  CD1 CD2                            '/
      data scname(41)(41:56)/'                '/
      data scname(42)( 1:40)/' CG  CD1 CD2                            '/
      data scname(42)(41:56)/'                '/
      data scname(43)( 1:40)/' CG                                     '/
      data scname(43)(41:56)/'                '/
      data scname(44)( 1:40)/' CG  CD  CE                             '/
      data scname(44)(41:56)/'                '/
      data scname(45)( 1:40)/' CG1 CG2                                '/
      data scname(45)(41:56)/'                '/
      data scname(46)( 1:40)/' CG1 CD1                                '/
      data scname(46)(41:56)/'                '/
      data scname(47)( 1:40)/' CG  CD1 CD2                            '/
      data scname(47)(41:56)/'                '/
      data scname(48)( 1:40)/' CG  CD  CE                             '/
      data scname(48)(41:56)/'                '/
      data scname(49)/allspc/
      data scname(50)( 1:40)/' CG  CD1 CD2                            '/
      data scname(50)(41:56)/'                '/
      data scname(51)( 1:40)/' CG1 CG2                                '/
      data scname(51)(41:56)/'                '/
      data scname(52)( 1:40)/' CG  CD1 CD2 CH  OH  CM                 '/
      data scname(52)(41:56)/'                '/
      data scname(53)/allspc/
      data scname(54)( 1:40)/' CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  CA2'/
      data scname(54)(41:56)/' CB2 CG2 CD3 CD4'/
      data scname(55)( 1:40)/' CG  CD1 CD2 CE1 CE2 CZ                 '/
      data scname(55)(41:56)/'                '/
      data scname(56)/allspc/
      data scname(57)( 1:40)/' CG1 CG2 CD1                            '/
      data scname(57)(41:56)/'                '/
      data scname(58)( 1:40)/' CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM  N1 '/
      data scname(58)(41:56)/' CB2 CG2 CD3 CD4'/
      data scname(59)( 1:40)/' CG  CD1 CD2 CE1 CE2 CZ  CH  OH  CM     '/
      data scname(59)(41:56)/'                '/
      data scname(60)( 1:40)/' CG SED  CE                             '/
      data scname(60)(41:56)/'                '/
      data scname(61)( 1:40)/' CG  CD1 CD2                            '/
      data scname(61)(41:56)/'                '/
      DATA scname(62)(1:40)/' CG  CD1                                '/
      DATA scname(62)(41:56)/'                '/

      data altpos / 7*0, 1, 5*0, 2, 2*0, 3, 46*0 /
      data altnam(1) ( 1:40)/' CG  AD1 AD2 AE1 AE2                    '/
      data altnam(1) (41:56)/'                '/
      data altnam(2) ( 1:40)/' CG  AD1 AD2                            '/
      data altnam(2) (41:56)/'                '/
      data altnam(3) ( 1:40)/' CG  CD  AE1 AE2                        '/
      data altnam(3) (41:56)/'                '/
      data cbet /' CB '/
      save scname, altnam, altpos, nscats, cbet

c
c     Routine to read the brookhaven data file and extract the sequence
c     and the atom names and coordinates. Checking is also done for
c     residue compatibility between sequence and atom. C-alpha only
c     files are flagged. The sequence residues are expected together,
c     then the secondary structure records (sheets together), and then
c     the atom records which are expected in residue order.
c
      seqlen = 0
      hedlen = 0
      cmplen=0
      clslen=0
      soulen=0
      inseqr = space
      sqrsct = 0
      brcode = '    '
      astruc = .false.
  100 continue
      read(brkfi,'(A)',end=3000) brkrec
      keywd = brkrec(1:keylen)
      if (keywd .eq. seqkey) go to 200
      if (keywd .eq. hedkey) go to 500
      if (keywd .eq. cmpkey) go to 600
      if (keywd .eq. soukey) go to 700
      if (keywd .eq. atmkey) go to 900
      if (keywd .eq. htakey) go to 900
      if (keywd .eq. modkey) 
     +    write(*,*) 
     +   'This is an NMR file. Data generated for first structure'
      if (keywd .eq. hlxkey .or. keywd .eq. shtkey .or.
     +    keywd .eq. trnkey .or. keywd .eq. dsfkey) then
         if (.not. astruc) then
            open(unit = authfi, status = 'SCRATCH')
            astruc = .true.
         end if
         write(authfi,'(a)') brkrec
      end if
      go to 100

  200 continue
      read (brkrec,220) (aain(i),i = 1, aapln)
  220 format(18x,13(1x,a3))
      do 250, i = 1, aapln
      if (aain(i) .ne. '   ') then
         if (resfnd(inseqr,aain(i)) .eq. 0) then
            sqrsct = sqrsct + 1
            inseqr((sqrsct-1)*3+1:sqrsct*3) = aain(i)
         end if
      end if
  250 continue
      go to 100
c
c     set up header records
c
  500 continue
      read(brkrec,520) contch, wkstr, incode
  520 format(9x,a1,a40,12x,a4)
      if (contch .eq. space) then
         brcode = incode
         class = wkstr
         clslen = strlen(class)
      end if
      go to 100
  600 continue
      read(brkrec,620) contch, wkstr
  620 format(9x,a1,a50)
      if (contch .eq. space) then
         compnd = wkstr
         cmplen = strlen(compnd)
      end if
      go to 100
  700 continue
      read(brkrec,720) contch, wkstr
  720 format(9x,a1,a50)
      if (contch .eq. space) then
         source = wkstr
         soulen = strlen(source)
      end if
      go to 100
c
c     set up header if any exists
c
  900 continue
      headln(1:6) = brcode // '  '
      hedlen = 6
      if (clslen .gt. 0) then
         headln(7:8+clslen) = class(1:clslen) // '  '
         hedlen = hedlen + clslen + 2
      end if
      if (cmplen .gt. 0) then
         hedend = hedlen + cmplen + 2
         if (hedend .gt. 132) hedend = 132
         headln(hedlen+1:hedend) = compnd(1:cmplen) // '  '
         hedlen = hedend
      end if
      if (soulen .gt. 0) then
         hedend = hedlen + soulen+ 2
         if (hedend .gt. 132) hedend = 132
         headln(hedlen+1:hedend)  = source(1:soulen)
         hedlen = hedend
      end if
      if (hedlen .gt. 132) hedlen = 132
c
c     look for next atom records to process
c
 1000 continue
      if (brkrec(1:keylen) .ne. atmkey .and.
     +   brkrec(1:keylen) .ne. htakey) then
         read(brkfi,'(A)',end=3000) brkrec
         go to 1000
      end if
c
c     have first atom record for an aa
c     for each aa get the seqnum and resnam and check against the
c     read in sequence. If ok initialise to no atoms present and
c     get the coords
c
      aacnt = 0
      do 1234 i = 1, maxres
      aastd(i) = space
 1234 continue
      read(brkrec,1080)oldchn
 1080 format(21x,a1)
 
 1100 continue
      keywd = brkrec(1:keylen)
      aacode = 0
      altset = .false.
      setatm = '    '
      read(brkrec,1120) oldatm, setloc, resnam, seqnum, oldocc, oldthm
 1120 format(12x,a4,a1,a3,1x,a6,27x,2f6.2)
      if (sqrsct .ne. 0 .and. resfnd(inseqr,resnam) .eq. 0) go to 1400
      if (resfnd(dudcod,resnam) .ne. 0) go to 1400
      respl = resfnd(modcod,resnam)
      if (respl .ne. 0.and .aacnt .ne. 0) then
         respl = (respl - 1) / 3 + 1
         if (Modpl(respl:respl) .eq. '+') then
            aastd(aacnt+1) = stdmod(respl:respl)
         else
            if (aastd(aacnt) .eq. 'S') then
               aastd(aacnt) = stdmod(respl:respl)
            else
               if (aastd(aacnt) .eq. 'O') then
                  aastd(aacnt) = nstdmd(respl:respl)
               else
                  if (llt(aastd(aacnt),'a')) then
                     aastd(aacnt) = stdmod(nmod+1:nmod+1)
                  else
                     aastd(aacnt) = nstdmd(nmod+1:nmod+1)
                  end if
               end if
            end if
         end if
         go to 1400
      end if
      aacode = AAconv(resnam, aacod3, hetcod)

      if (aacode .le. 0) go to 1400
      aacnt = aacnt + 1
      if (aastd(aacnt) .eq. space) then
         aastd(aacnt) = 'S'
         if (aacode .gt. stdaa) aastd(aacnt) = 'O'
      else
         if (aacode .gt. stdaa) then
            wkpl = index(stdmod,aastd(aacnt))
            aastd(aacnt) = nstdmd(wkpl:wkpl)
         end if
      end if
      seqcod(aacnt) = aacode
      seqbcd(aacnt) = seqnum
      atcnt = 0
      scatct = 0
      thmcnt = 0
      thmprs(aacnt) = space
      do 1200 i = 1, mnchna
      atomin(i,aacnt) = .false.
 1200 continue
      do 1220 i = 1, sdchna
      scatin(i,aacnt) = .false.
 1220 continue
      do 1240 i = 1, 3
      atpres(i,aacnt) = space
 1240 continue
c
c     process cordinates find atom type and if main chain set in coords
c
 1300 continue
      read(brkrec,1320)atmnam,altloc,(acoord(i),i=1,3), occup, therm
 1320 format(12x,a4,a1,13x,3f8.3,2f6.2)
      if (atmnam .eq. ' OXT' .or. atmnam .eq. ' NXT' .or.
     +    atmnam(2:2) .eq. 'H' .or. atmnam(2:2) .eq. 'D')
     +    go to 1400
      attype = index(atname,atmnam)
      if (attype.eq.0)  then
         if (aacode.gt.43.and.aacode.le.49) then
            attype = index(batnam,atmnam)
         else if (aacode.eq.41)  then
            attype = index(tphnam,atmnam)
         else if (aacode.eq.40.or.aacode.eq.53)  then
            attype = index(olenam,atmnam)
         else if (aacode.eq.61)  then
            attype = index(lolnam,atmnam)
         else if (aacode.eq.38)  then
            attype = index(pphnam,atmnam)
         endif
      endif
      if (altloc .ne. space) then
         if (setatm .eq. '    ') then
            setatm = atmnam
            setloc = altloc
            oldocc = occup
            oldthm = therm
         end if
         if (.not. altset) then
            if (setatm .eq. atmnam) then
               if (occup .gt. oldocc) then
                  setloc = altloc
                  oldocc = occup
                  if (xney(oldthm,0.0)) thmcnt = thmcnt - 1
                  oldthm = therm
                  if (attype .ne. 0) then
                     atcnt = atcnt - 1
                  else
                     scatct = scatct - 1
                  end if
               end if
            else
               altset = .true.
            end if
         end if
      end if
      if (xney(therm,0.0).and. (altloc .eq. setloc .or.
     +          altloc .eq. space)) thmcnt = thmcnt + 1
      if (attype .ne. 0 .and. (altloc .eq. setloc .or.
     +   altloc .eq. space)) then
         attype = attype / 4 + 1
         if (attype.le.4) then
            do 1350 i = 1, ncoord
            mcaaat(i,attype,aacnt) = acoord(i)
 1350       continue
            atomin(attype,aacnt) = .true.
         endif
         if (attype .ne. calph) caonly = .false.
         atcnt = atcnt + 1
      end if
      if (attype .eq. 0 .and. (altloc .eq. setloc .or.
     +   altloc .eq. space)) then
         scatct = scatct + 1
         if (aacode .le. sideaa) then
            if (atmnam .eq. cbet) then
               sctype = 1
            else
               sctype = index(scname(aacode),atmnam)
               if (sctype .eq. 0 .and. altpos(aacode) .ne. 0)
     +            sctype = index(altnam(altpos(aacode)),atmnam)
               if (sctype .ne. 0) sctype = sctype / 4 + 2
            end if
            if (sctype .ne. 0) then     
               scatin(sctype,aacnt) = .true.
               do 1375 i = 1, ncoord
               scaaat(i,sctype,aacnt) = acoord(i)
 1375          continue
            end if
         end if
      end if
c
c     search for the next atom record and treat as appropriate
c     process coords only if same aa 
c     if new aa start aa process again after setting the count fields
c
 1400 continue
      read(brkfi,'(a)',end=2000) brkrec
      If (brkrec(1:keylen) .ne. endkey.and.brkrec(1:keylen).ne.
     +     'ENDMDL') then
         If (brkrec(1:keylen) .eq. atmkey .or.
     +      brkrec(1:keylen) .eq. htakey) then
            read(brkrec,'(21x,a6)') innum
            if (innum.eq.seqnum.and.keywd.eq.brkrec(1:keylen)) then
               if (aacode .le. 0) go to 1400
               go to 1300
            else
	       if (aacode .gt. 0)
     +            call cntset(atpres, thmprs, altcod, nscats, aacode, 
     +                        aacnt, atcnt, scatct, thmcnt, setloc)
               If (aacnt .ge. maxres) then
                  Write(errfi,*)'Too many amino acids',
     +                          '  - input terminated'
                  goto 2100
               else
                  if (aacode .gt. 0) aastd(aacnt + 1) = space
                  go to 1100
               end if
	    end if
         end if
         go to 1400
      end if
 2000 continue
      if (aacode .gt. 0)
     +   call cntset(atpres, thmprs, altcod, nscats, aacode, 
     +               aacnt, atcnt, scatct, thmcnt, setloc)
 2100 continue
      Seqlen = aacnt
 3000 continue
      if (seqlen .eq. 0) readok = .false.
      close(unit=brkfi)
      END
c*****************************************************************************
      Subroutine cntset(atpres, thmprs, altcod, nscats, aacode, 
     +                  aacnt,  atcnt,  scatct, thmcnt, setloc)
c
c     perform end of amino acid flag setting
c
      include 'p_sstruc.par'
      character atpres(3,maxres)*1, thmprs(maxres)*1, altcod(maxres)*1,
     +          setloc*1
      integer   nscats(sideaa), aacode, aacnt, scatct, thmcnt, atcnt
      integer   totats
      character atmchk*1

      altcod(aacnt) = setloc
      atpres(1,aacnt) = space
      atpres(2,aacnt) = space
      atpres(3,aacnt) = space
      thmprs(aacnt) = space
      totats = atcnt + scatct
      if (aacode .le. sideaa) then
         if (nscats(aacode) .ge. 0) then
            atpres(1,aacnt) = atmchk(mnchna+nscats(aacode)-2,totats)
            atpres(2,aacnt) = atmchk(mnchna-2,atcnt)
            atpres(3,aacnt) = atmchk(nscats(aacode),scatct)
            thmprs(aacnt) = atmchk(mnchna+nscats(aacode)-2,thmcnt)
         end if
      end if
      end
c*****************************************************************************
      character*1 function atmchk(Total,count)
c
c     set 'present' fields to all none or some as appropriate
c
      include 'p_sstruc.par'
      integer total, count
      if (total .eq. 0) then
         atmchk = space
      else
         if (count .eq. 0) then
            atmchk = 'N'
         else
            if (count .ge. total) then
               atmchk = 'A'
            else
               atmchk = 'S'
            end if
         end if
      end if
      END
c*****************************************************************************
      Integer function Strlen(string)
c
c     find out the real length of the string
c
      include 'p_sstruc.par'
      character*(*) string

      integer nxch
      do 200, nxch = len(string), 1, -1
        if (string(nxch:nxch) .ne. space) goto 300
  200 continue
      nxch = 0
  300 continue
      strlen = nxch
      end
c****************************************************************************
      Integer function AAconv(aacode, aacod3, hetcod)
c
c     find the numerical position of an amino acid code
c
      include 'p_sstruc.par'
      character aacode*3, aacod3*(3*abetsz), hetcod*102
      integer   resfnd
      integer   aaval
      aaval = resfnd(aacod3,aacode)
      if (aaval .eq. 0) then
         aaval = resfnd(hetcod,aacode)
         if (aaval .eq. 0) then
            aaval = unkcod
            write(errfi,530) aacode, ' unknown amino acid code'
  530       format(a3,a)
         else
            aaval = abetsz + aaval / 3 + 1
         end if
      else
         aaval = aaval / 3 + 1
      end if
      AAconv = aaval
      END
c*****************************************************************************
      Integer function resfnd(soustr,questr)
c
c     locate query in source on 3 character boundaries
c
      character soustr*(*), questr*3
      integer   respl, nxpl

      respl = Index(soustr,questr)
  100 continue
      If (mod(respl,3) .ne. 1 .and. respl .ne. 0) then
         nxpl = index(soustr(respl+1:),questr)
         if (nxpl .ne. 0) then
            respl = respl + nxpl
         else
            respl = nxpl
         end if
         go to 100
      end if
      resfnd = respl
      End
c*****************************************************************************
      Subroutine Chnbrk(mcaaat, atomin, seqbcd, brksym, chnsz, chnct,
     +                  chnend, caonly, seqlen)
c
c     search through looking for distant peptide bonds. If any over
c     the limit are found or if a c or n atom does not exist then a
c     chain break is deemed to occur
c
      include 'p_sstruc.par'
      real      mcaaat(ncoord,mnchna,maxres)
      logical   atomin(mnchna,maxres), caonly
      Integer   chnsz(maxchn), chnct, seqlen, chnend(maxchn)
      character seqbcd(maxres)*6, brksym(maxres)*1
      Real      atdist

      Integer   nxres, ninchn

      do 200 nxres = 1, seqlen
      brksym(nxres) = space
  200 continue
      ninchn = 1
      chnct = 1
      do 1000, nxres = 2, seqlen
      if (caonly) then
         if (atdist(mcaaat(1,calph,nxres-1),mcaaat(1,calph,nxres))
     +      .gt. cadbnd) then
            call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                  nxres, seqlen)
C     6/3/96 TO STOP CRASHING FOR FILES WITH LOTS OF CALPHA-ONLY RESIDUES
            IF(CHNCT.GT.MAXCHN) GO TO 1500
         else
            ninchn = ninchn + 1
            if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) 
     +         Write(errfi,500) nxres-1, seqbcd(nxres-1),
     +                          nxres, seqbcd(nxres)
  500       format('chain letter change but no chain break between ',
     +             i4,2x,a6,' and',i4,2x,a6)
         end if
      else
         If (atomin(carb,nxres-1) .and. atomin(nmain,nxres)) then
            if (atdist(mcaaat(1,carb,nxres-1),mcaaat(1,nmain,nxres))
     +          .gt. pepbnd) then
               call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                     nxres, seqlen)
C     6/3/96 TO STOP CRASHING FOR FILES WITH LOTS OF CALPHA-ONLY RESIDUES
               IF(CHNCT.GT.MAXCHN) GO TO 1500
            else
               ninchn = ninchn + 1
               if (seqbcd(nxres-1)(1:1) .ne. seqbcd(nxres)(1:1)) 
     +            Write(errfi,500) nxres-1, seqbcd(nxres-1),
     +                  nxres, seqbcd(nxres)
            end if
         else
            call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                 nxres, seqlen)
C     6/3/96 TO STOP CRASHING FOR FILES WITH LOTS OF CALPHA-ONLY RESIDUES
            IF(CHNCT.GT.MAXCHN) GO TO 1500
         end if
      end if
 1000 continue
      call chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd,
     +            seqlen+1, seqlen)
C     6/3/96 TO STOP CRASHING FOR FILES WITH LOTS OF CALPHA-ONLY RESIDUES
      IF(CHNCT.GT.MAXCHN) GO TO 1500
 1500 CONTINUE
      END
c*****************************************************************************
      SUBROUTINE Chnset(brksym, chnsz, chnct, chnend, ninchn, seqbcd, 
     +                  nxres, seqlen)
      include 'p_sstruc.par'
      integer   chnsz(maxchn), chnct, ninchn, seqlen, nxres,
     +          chnend(maxchn)
      character seqbcd(maxres)*6, brksym(maxres)*1

      if(chnct.gt.maxchn) then
         close(unit=outfi,status='DELETE')
         write(*,*) 'too many chains in this protein'
         go to 150
      endif
      chnsz(chnct) = ninchn
      if (chnct .eq. 1) then
         chnend(chnct) = ninchn
      else
         chnend(chnct) = chnend(chnct-1) + ninchn
      end if
      ninchn = 1
      if (nxres .le. seqlen) then
         brksym(nxres-1) = brkch
         brksym(nxres) = brkch
         write(errfi,100) nxres-1, seqbcd(nxres-1), nxres, seqbcd(nxres)
  100    format('chain break between ',i4,'(',a6,') and ',i4,'(',a6,')')
         chnct = chnct + 1
         if(chnct.gt.maxchn) then
            close(unit=outfi,status='DELETE')
            write(*,*) 'too many chains in this protein'
            go to 150
c            stop
         endif
      end if
 150  end
c*****************************************************************************
      Subroutine mkcadt(mcaaat, atomin, cadsts, seqbcd, seqlen)
c
c     generate the calpha distances for I+1 to I+11. take note of nulls
c     accross brookhaven chain breaks but not accross disconections at 
c     this point
c
      include 'p_sstruc.par'
      integer   seqlen
      real      mcaaat(ncoord,mnchna,maxres), cadsts(11,maxres)
      logical   atomin(mnchna,maxres)
      character seqbcd(maxres)*6

      integer   nxres, nxdist
      real      mxdist, adist
      parameter (mxdist = 99.9)
      real      atdist

      do 1000, nxres = 1, seqlen
      do 500, nxdist = 1, 11
      if (nxres+nxdist .gt. seqlen .or.
     +   .not. atomin(calph,nxres)) then
         cadsts(nxdist,nxres) = -1
      else
         if (seqbcd(nxres)(1:1) .ne. seqbcd(nxres+nxdist)(1:1)
     +      .or. .not. atomin(calph,nxres+nxdist)) then
            cadsts(nxdist,nxres) = -1
         else
            adist = atdist(mcaaat(1,calph,nxres),
     +                     mcaaat(1,calph,nxres+nxdist))
            if (adist .gt. mxdist) adist = mxdist
            cadsts(nxdist,nxres) = adist
         end if
      end if
  500 continue
 1000 continue
      END
c*****************************************************************************
      Subroutine CreatH(mcaaat, atomin, chnsz, chnct)
c
c     Procedure to create hydrogen atoms. These are deemed to be
c     one angstrom from the nitrogen atom and angled so that N-H
c     is parallel to C=O before. By taking the vector from O to C
c     the N-H bond is generated on the other side of C-N to the O.
c
      INCLUDE 'p_sstruc.par'

      integer chnsz(maxchn), chnct
      real    mcaaat(ncoord,mnchna,maxres)
      logical atomin(mnchna,maxres)
c
c     local variables
c
      Integer nxres, i, chnst, nxchn
      Real    atvec(ncoord), colen
c
      chnst = 0
      do 2000, nxchn = 1, chnct
      do 1900, nxres = chnst + 2, chnst + chnsz(nxchn)
      If (atomin(nmain,nxres) .and. atomin(carb,nxres-1) .and.
     +   atomin(oxyg,nxres-1)) then
         Colen = 0.0
         do 300 i = 1,ncoord
         atvec(i) = mcaaat(i,carb,nxres-1) - mcaaat(i,oxyg,nxres-1)
         colen = colen + atvec(i) ** 2
  300    continue
         colen = sqrt(colen)
         do 500 i = 1, ncoord
         atvec(i) = atvec(i) / colen
  500    continue
         do 700 i = 1, ncoord
         mcaaat(i,nhyd,nxres) = mcaaat(i,nmain,nxres) + atvec(i)
  700    continue
         atomin(nhyd,nxres) = .true.
      else
         write(errfi,1100) nxres
 1100    format('hydrogen atom not generated for residue ',i4)
         atomin(nhyd,nxres) = .false.
      end if
 1900 continue
      chnst = chnst + chnsz(nxchn)
 2000 continue
      END
c*****************************************************************************
      Subroutine MkHBnd(Mcaaat, atomin, Hbond, Hbonde, seqcod, chnend,
     +                  seqlen)
c
c     Procedure to discover the possible hydrogen bonds between two
c     amino acids. Each amino acid is allowed 4 bonds - 2 from the 
c     C=O and 2 from the N-H. Should it be that 3 or more bonds are
c     possible the most favourable are kept. A hydrogen
c     bond is deemed to exist if the Kabsch and Sander energy formula
c     gives a value less than the cutoff. The angle and distance
c     approach is implicitly built in to this method. The formula
c     corresponds to a maximum distance of 5.2 A for in line atoms
c     and a maximum angular deviation of 63 degrees. A check on the
c     maximum distance is used to abandon unnecessary calculations.
c     the formula is
c     E = q1*q2*(1/d(ON)+1/d(CH)-1/d(OH)-1/d(CN))*f
c     Residues must be separated by at least one residue to bond
c     Proline cannot be a donor of a hydrogen bond
c
      INCLUDE 'p_sstruc.par'

      Integer Seqlen
      Real    Mcaaat(ncoord,mnchna,maxres), Hbonde(maxbnd,maxres)
      Integer Hbond(maxbnd,maxres), seqcod(maxres), chnend(maxchn)
      logical Atomin(mnchna,maxres)
      Real    Atdist
c
c     local variables
c
      integer    ksf, mndist
      real       mxdist, maxeng, ksq1, ksq2, maxcad, englim
      parameter (mxdist = 5.2, maxeng = -0.5, mndist = 0.5,
     +           englim = -9.9, ksq1 = 0.42, ksq2 = 0.2, ksf = 332, 
     +           maxcad = 8.0)
      Integer    Nxres, Othres, i, nbonds, fstchn, othchn
      logical    nhpres, xeqy
      real       don, doh, dch, dcn, energy, cadist

      do 210 nxres = 1, seqlen
      do 200 i = 1,maxbnd
      hbond(i,nxres) = 0
      hbonde(i,nxres) = 0.0
  200 continue
  210 continue
      nbonds = 0
      fstchn = 1
      do 2000 nxres = 1, seqlen
      if (nxres .gt. chnend(fstchn)) fstchn = fstchn + 1
      nhpres = .false.
      if (atomin(nmain,nxres) .and. atomin(nhyd,nxres) .and.
     +   seqcod(nxres) .ne. procod) nhpres = .true.
      if (nhpres) then
         othchn = 1
         do 1900 othres = 1, seqlen
         if (othres .gt. chnend(othchn)) othchn = othchn + 1
         If ((iabs(nxres - othres) .eq. 1 .and. fstchn .ne. othchn)
     +      .or. iabs(nxres - othres) .ge. 2) then
            If (atomin(calph,nxres) .and. atomin(calph,othres)) then
               cadist = atdist(mcaaat(1,calph,othres),
     +                  mcaaat(1,calph,nxres))
               If (cadist .lt. maxcad) then
                  If (atomin(carb,othres) .and. 
     +               atomin(oxyg,othres))then
                     don = atdist(mcaaat(1,oxyg,othres),
     +                     mcaaat(1,nmain,nxres))
                     doh = atdist(mcaaat(1,oxyg,othres),
     +                     mcaaat(1,nhyd,nxres))
                     dch = atdist(mcaaat(1,carb,othres),
     +                     mcaaat(1,nhyd,nxres))
                     dcn = atdist(mcaaat(1,carb,othres),
     +                     mcaaat(1,nmain,nxres))
                     If (xeqy(don,0.0).or.xeqy(doh,0.0).or.
     +                   xeqy(dch,0.0).or.xeqy(dcn,0.0)) then
                        write(errfi,250)nxres, othres
  250                   format('coincident atoms in hydrogen bonding',
     +                         ' donor ',I4,' acceptor ',i4)
                     else
                        energy = ksq1 * ksq2 * ksf * 
     +                          (1.0/don + 1.0/dch - 1.0/doh - 1.0/dcn)
                        If (energy .lt. englim) then
                           write(errfi,300) nxres, othres, don, dch,
     +                                doh, dcn
  300                      format('atoms too close ',I4,2x,i4,
     +                           2x,4(2x,f8.3))
                           energy = englim
                        end if
                        If (energy .lt. maxeng) call engset(hbonde,
     +                     hbond, energy, nxres, othres, nbonds)
                     end if
                  end if
               end if
            end if
         end if
 1900    continue
      end if
 2000 continue
      write(errfi,2100) 'number of hydrogen bonds is ',nbonds
 2100 format(a,i5)
      END
c*****************************************************************************
      Subroutine Engset(hbonde, hbond, energy, donres, accres, bndct)
c
c     find a place to put this h-bond. discard third bonds if necessary
c
      Include 'p_sstruc.par'
      real    hbonde(maxbnd,maxres), energy
      integer hbond(maxbnd,maxres), donres, accres, bndct

      integer acpl, donpl, rejptr
      integer engpl

      donpl = engpl(energy,hbonde(nst,donres))
      acpl = engpl(energy,hbonde(cst,accres))
      If (acpl .eq. 0 .or. donpl .eq. 0) then
         write(errfi,400) donres, accres, energy
  400    format('third (+) Hbond (N-C) ',I4,2x,i4,
     +                          ' energy ',f6.2,' abandoned')
      else
         if (hbond(nst+1,donres) .ne. 0) then
            write(errfi,400) donres, hbond(nst+1,donres),
     +                    hbonde(nst+1,donres)
            rejptr = hbond(nst+1,donres) 
            call rejbnd(donres,hbond(cst,rejptr),hbonde(cst,rejptr))
            bndct = bndct - 1
         end if
         if (hbond(cst+1,accres) .ne. 0) then
            write(errfi,400) hbond(cst+1,accres), accres, 
     +                    hbonde(cst+1,accres)
            rejptr = hbond(cst+1,accres) 
            call rejbnd(accres,hbond(nst,rejptr),hbonde(nst,rejptr))
            bndct = bndct - 1
         end if
         If (acpl .eq. 1) then
            hbond(cst+1,accres) = hbond(cst,accres)
            hbonde(cst+1,accres) = hbonde(cst,accres)
         end if
         If (donpl .eq. 1) then
            hbond(nst+1,donres) = hbond(nst,donres)
            hbonde(nst+1,donres) = hbonde(nst,donres)
         end if
         hbond(cst+acpl-1,accres) = donres
         hbonde(cst+acpl-1,accres) = energy
         hbond(nst+donpl-1,donres) = accres
         hbonde(nst+donpl-1,donres) = energy
         bndct = bndct + 1
      end if
      END
c*****************************************************************************
      SUBROUTINE rejbnd(rejptr,bond,bonde)
      include 'p_sstruc.par'
      integer rejptr, bond(2)
      real    bonde(2)
      integer rejpl

      rejpl = 1
      if (bond(1).ne. rejptr) rejpl = rejpl + 1
      if (rejpl .eq. 1) then
         bond(1) = bond(2)
         bonde(1) = bonde(2)
      end if
      bond(2) = 0
      bonde(2) = 0.0
      end
c*****************************************************************************
      INTEGER FUNCTION Engpl(energy, bonde)
      include 'p_sstruc.par'
      real energy, bonde(2)
      engpl = 0
      If (energy .lt. bonde(1)) then
         engpl = 1
      else
         If (energy .lt. bonde(2)) engpl = 2
      end if
      end
c****************************************************************************   
      REAL Function Atdist(coorda,coordb)
c
c     routine to determine the distance in space between the two
c     atoms represented by the coordinate arrays.
c
      INCLUDE 'p_sstruc.par'

      Real    Coorda(ncoord),Coordb(ncoord)
c
      Integer i
      Real    Dist
      dist = 0.0
      do 200 i = 1,ncoord
      dist = dist + (coorda(i) - coordb(i)) ** 2
  200 continue
      atdist = sqrt(dist)
      END
c*****************************************************************************
      INTEGER Function Aplace(Bondpl)
c
c     find free place in bond table or return 0
c
      include 'p_sstruc.par'
      integer bondpl(2), i, wkpl
      do 200 i = 1,2
      If (Bondpl(i) .eq. 0) then
         wkpl = i
         go to 250
      end if
  200 continue
      wkpl = 0
  250 continue
      aplace = wkpl
      END
c*****************************************************************************
      Subroutine Ooi(atomin, mcaaat, seqbcd, oois, seqlen)
c
c     routine to calculate Ooi numbers. these are done for n8 and n14.
c     all atoms are included as there is noway to handle i+/-2 over
c     missing densities. Numbers are for within brookhaven chain only.
c
      include 'p_sstruc.par'
      logical   atomin(mnchna,maxres)
      real      mcaaat(ncoord,mnchna,maxres)
      character seqbcd(maxres)*6
      integer   oois(2,maxres), seqlen

      integer   nxaa, stchn, wkchn(maxchn), stpl, nxchn, nxpart,
     +          numchn
      character chnch*1
      real      wkdist
      real      atdist

      stchn = 1
      chnch = seqbcd(1)(1:1)
      numchn = 0
      do 200, nxaa = 1, seqlen
      if (atomin(calph,nxaa)) then
         oois(1,nxaa) = 0
         oois(2,nxaa) = 0
      else
         oois(1,nxaa) = -1
         oois(2,nxaa) = -1
      end if
      if (seqbcd(nxaa)(1:1) .ne. chnch) then
         numchn = numchn + 1
         wkchn(numchn) = nxaa - stchn
         chnch = seqbcd(nxaa)(1:1)
         stchn = nxaa
      end if
  200 continue
      numchn = numchn + 1
      wkchn(numchn) = seqlen + 1 - stchn
      stpl = 0
      do 3000, nxchn = 1, numchn
      do 2000, nxaa = stpl + 1, stpl + wkchn(nxchn) - 1
      do 1000, nxpart = nxaa + 1, stpl + wkchn(nxchn)
      if (atomin(calph,nxaa) .and. atomin(calph,nxpart)) then
         wkdist = atdist(mcaaat(1,calph,nxaa),
     +            mcaaat(1,calph,nxpart))
         if (wkdist .lt. ooi2rd) then
            oois(ooi2,nxaa) = oois(ooi2,nxaa) + 1
            oois(ooi2,nxpart) = oois(ooi2,nxpart) + 1
            if (wkdist .lt. ooi1rd) then
               oois(ooi1,nxaa) = oois(ooi1,nxaa) + 1
               oois(ooi1,nxpart) = oois(ooi1,nxpart) + 1
            end if
         end if
      end if
 1000 continue
 2000 continue
      stpl = stpl + wkchn(nxchn)
 3000 continue
      end      
c*****************************************************************************
      Subroutine MkAngl(McAAat, McAngs, Atomin, chnsz, chnct, caonly)
c
c     Procedure to calculate the main chain dihedral angles
c     if all the atoms are present.
c
      INCLUDE 'p_sstruc.par'

      Real    McAAat(ncoord,mnchna,maxres), McAngs(nmnang,maxres)
      Logical Atomin(mnchna,maxres), caonly
      Integer chnsz(maxchn), chnct
      Real    Dihed 
c
c     local variables
c
      Integer angrng(nmnang,3), angatm(nmnang,4), angoff(nmnang,4),
     +        angchn(nmnang), I, j, Nxres, Nxset, chnst, nxchn, stang,
     +        finang

      data ((angrng(I,j),j=1,3), i=1,nmnang)
     +      /2, 0, phi, 1, -1, psi, 1, -1, omega, 2, -2, chiral,
     +       3, -2, kappa, 2, 0, tco/
      data ((angatm(I,j),j=1,4), i=1,nmnang)
     +     /carb, nmain, calph, carb, nmain, calph, carb, nmain,
     +      calph, carb, nmain, calph, calph, calph, calph, calph,
     +      calph, calph, calph, calph, carb, oxyg, carb, oxyg/
      data ((angoff(I,j),j=1,4), i=1,nmnang)
     +     /-1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, -1, 0, 1, 2,
     +      -2, 0, 0, 2, 0, 0, -1, -1/
      data angchn / 2, 2, 2, 4, 5, 2/
      SAVE angrng, angatm, angoff, angchn
c
      stang = 1
      finang = nmnang
      if (caonly) then
         stang = chiral
         finang = kappa
      end if
      do 4000, Nxset = stang, finang
      chnst = 0
      do 3000, nxchn = 1, chnct
      if (chnsz(nxchn) .lt. angchn(nxset)) then
         do 1000, nxres = chnst + 1, chnst + chnsz(nxchn)
         McAngs(angrng(nxset,3),Nxres) = Null
 1000    continue
      else
         do 2000, Nxres = chnst + angrng(nxset,1),
     +                    chnst + chnsz(nxchn) + angrng(nxset,2)
         If (atomin(angatm(nxset,1), nxres+angoff(nxset,1)) .and.
     +      atomin(angatm(nxset,2), nxres+angoff(nxset,2)) .and.
     +      atomin(angatm(nxset,3), nxres+angoff(nxset,3)) .and.
     +      atomin(angatm(nxset,4), nxres+angoff(nxset,4))) then
            McAngs(angrng(nxset,3),Nxres) = Dihed(nxset,
     +          Mcaaat(1,angatm(nxset,1), nxres+angoff(nxset,1)),
     +          Mcaaat(1,angatm(nxset,2), nxres+angoff(nxset,2)),
     +          Mcaaat(1,angatm(nxset,3), nxres+angoff(nxset,3)),
     +          Mcaaat(1,angatm(nxset,4), nxres+angoff(nxset,4)))
         else
            McAngs(angrng(nxset,3),Nxres) = Null
         end if
 2000    continue
         do 2400, nxres = chnst + 1, chnst + angrng(nxset,1) - 1
         McAngs(angrng(nxset,3),nxres) = Null
 2400    continue
         do 2600, nxres = chnst + chnsz(nxchn) + angrng(nxset,2) + 1,
     +                    chnst + chnsz(nxchn)
         McAngs(angrng(nxset,3),nxres) = Null
 2600    continue
      end if
      chnst = chnst + chnsz(nxchn)
 3000 continue
 4000 continue
      END
c*****************************************************************************
      subroutine Mksang(mcaaat, atomin, scaaat, scangs, scatin, 
     +                  seqcod, aacod3, hetcod, seqlen)
c
c     calculate chi angles from atom data. check for that alternates
c     are closest to c alpha
c
      include 'p_sstruc.par'
      real      mcaaat(ncoord,mnchna,maxres), 
     +          scaaat(ncoord,sdchna,maxres),
     +          scangs(nsdang,maxres)
      character aacod3*(abetsz*3), hetcod*102
      logical   atomin(mnchna,maxres), scatin(sdchna,maxres)
      integer   seqcod(maxres), seqlen

      integer   natsdc(sideaa), nxang, i, j, nxchi, nxres, nchi, aacod,
     +          extprc(sideaa), swatm, extatm, chipl, nswap, wkswap
      logical   athere(maxchi+3), xney
      real      relats(ncoord,maxchi+3), distmn, distal, ang1, ang2
      character outerr*80, wkres*3
      real      atdist, dihed
      data natsdc / 0, 0, 1, 2, 3, 2, 0, 2, 2, 0, 4, 2, 3,
     +              2, 0, 1, 3, 5, 1, 1, 0, 1, 2, 0, 2, 0, 0, 3, 4,
     +              0, 2, 1, 1, 1, 2, 4, 0, 2, 0, 2, 2, 1, 3, 1,
     +              2, 2, 3, 0, 2, 1, 5, 2, 0, 6, 2, 0, 2, 6, 2, 3,
     +              2, 1, 1 /
      data extprc / 0, 0, 0, 4, 4, 4, 0, 0, 5, 0, 0, 2, 0,
     +              0, 0, 0, 0, 4, 0, 3, 0, 2, 0, 0, 4, 0, 0, 0, 0,
     +              0, 4, 0, 0, 0, 4, 0, 0, 4, 0, 2, 4, 0, 0, 2,
     +	            5, 4, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     +              0, 0, 0 /
      save natsdc, extprc

      nswap = 0
      outerr = space
      do 5000 nxres = 1, seqlen
      aacod = seqcod(nxres)
      nchi = 0
      nchi = natsdc(aacod)
      do 100, nxang = nchi + 1, nsdang
      scangs(nxang,nxres) = null
  100 continue
      if (nchi .gt. 0) then
         athere(1) = atomin(nmain,nxres)
         athere(2) = atomin(calph,nxres)
         do 200, i = 1, ncoord
         relats(i,1) = mcaaat(i,nmain,nxres)
         relats(i,2) = mcaaat(i,calph,nxres)
  200    continue
         extatm = 0
         if (extprc(aacod) .gt. 1) extatm = 1
         do 400 i = 1, nchi + 1 + extatm
         athere(2+i) = scatin(i,nxres)
         do 300, j = 1, ncoord
         relats(j,2+i) = scaaat(j,i,nxres)
  300    continue
  400    continue
         swatm = 0
         if (extprc(aacod) .eq. 1) then
            if (atomin(calph,nxres) .and. scatin(nchi+1,nxres)
     +         .and. scatin(nchi+2,nxres)) then
               distmn = atdist(mcaaat(1,calph,nxres),
     +                  scaaat(1,nchi+1,nxres))
               distal = atdist(mcaaat(1,calph,nxres),
     +                  scaaat(1,nchi+2,nxres))
               if (distal .lt. distmn) then 
                  swatm = 1
                  athere(nchi+3) = scatin(nchi+2,nxres)
                  do 500, i = 1, ncoord
                  relats(i,nchi+3) = scaaat(i,nchi+2,nxres)
  500             continue
               end if
            end if
         end if
         do 1000, nxchi = 1, nchi
         If (athere(nxchi) .and. athere(nxchi+1) .and.
     +      athere(nxchi+2) .and. athere(nxchi+3)) then
            scangs(nxchi,nxres) = dihed(nmnang+nxchi,
     +            relats(1,nxchi), relats(1,nxchi+1),
     +            relats(1,nxchi+2), relats(1,nxchi+3))
         else
            scangs(nxchi,nxres) = null
         end if
 1000    continue
         if (extprc(aacod) .gt. 1 .and. 
     +      xney(scangs(nchi,nxres), null)) then
            If (athere(nchi) .and. athere(nchi+1) .and.
     +         athere(nchi+2) .and. athere(nchi+3+extatm)) then
               scangs(nchi+1,nxres) = dihed(nmnang+nchi+extatm,
     +            relats(1,nchi), relats(1,nchi+1),
     +            relats(1,nchi+2), relats(1,nchi+3+extatm))
            else
               scangs(nchi+1,nxres) = null
            end if
            if (xney(scangs(nchi+1,nxres), null)) then
               If (extprc(aacod) .eq. 2 .or. extprc(aacod) .eq. 3) then
                  ang1 = scangs(nchi,nxres)
                  ang2 = scangs(nchi+1,nxres)
                  if (ang1 .lt. 0.0) ang1 = 360.0 + ang1
                  if (ang2 .lt. 0.0) ang2 = 360.0 + ang2
                  chipl = nchi
                  if (abs(ang1 - ang2) .gt. 180.0) then
                     if (ang2 .gt. ang1 .and. extprc(aacod) .eq. 2)
     +                                                chipl = nchi + 1
                     if (ang1 .gt. ang2 .and. extprc(aacod) .eq. 3)
     +                                                chipl = nchi + 1
                  else
                     if (ang2 .lt. ang1 .and. extprc(aacod) .eq. 2)
     +                                                chipl = nchi + 1
                     if (ang1 .lt. ang2 .and. extprc(aacod) .eq. 3)
     +                                                chipl = nchi + 1
                  end if
                  if (chipl .ne. nchi) then
                     swatm = 1
                     scangs(nchi,nxres) = scangs(nchi+1,nxres)
                  end if
               end if
               if (extprc(aacod) .eq. 4) then
                  if (abs(scangs(nchi+1,nxres)) .lt. 
     +                         abs(scangs(nchi,nxres))) then
                     swatm = 1
                     scangs(nchi,nxres) = scangs(nchi+1,nxres)
                  end if
               end if
               scangs(nchi+1,nxres) = null
            end if
         end if
         if (swatm .eq. 1) then
            if (nswap .eq. 0) 
     +         write(errfi,'(A)') 'side chain atoms swapped for '
            wkswap = mod(nswap,8) + 1
            if (seqcod(nxres) .le. abetsz) then
               wkres = aacod3(seqcod(nxres)*3-2:seqcod(nxres)*3)
            else
               wkres = hetcod((seqcod(nxres)-abetsz)*3-2:
     +                       (seqcod(nxres)-abetsz)*3)
            end if
            write(outerr((wkswap-1)*10+1:wkswap*10),'(a3,1x,i4,2x)')
     +        wkres, nxres
            if (wkswap .eq. 8) then
               write(errfi,'(A)') outerr
               outerr = space
            end if
            nswap = nswap + 1
         end if
      end if
 5000 continue
      if (mod(nswap,8) .ne. 0) write(errfi,'(A)') outerr
      END
c**************************************************************************
      Real Function Dihed(angnum, AtomA, AtomB, AtomC, AtomD)
c
c     routine to calculate the dihedral angle between the 4 atoms
c     given
c
      INCLUDE 'p_sstruc.par'

      Real    Atoma(ncoord), Atomb(ncoord), Atomc(ncoord), atomd(ncoord)
      integer angnum
      Real    Atdist
      Logical xeqy
      Real    dihatm(ncoord,dihdat), codist(ncoord,dihdat-1),
     +        atmdst(dihdat-1), dotprd(dihdat-1,dihdat-1),
     +        ang, cosang, detant
      Integer i, j, k

      do 100, I = 1, ncoord
      Dihatm(i,1) = atoma(i)
      Dihatm(i,2) = atomb(i)
      Dihatm(i,3) = atomc(i)
      Dihatm(i,4) = atomd(i)
  100 continue

      do 200, i = 1, dihdat-1
      atmdst(i) = atdist(Dihatm(1,i),Dihatm(1,I+1))
      do 150, j = 1,ncoord
      codist(j,i) = dihatm(j,I+1) - dihatm(j,i)
  150 continue
  200 continue

      do 310, i = 1, dihdat-1
      do 300, j = 1, dihdat-1
      dotprd(i,j) = 0.0
  300 continue
  310 continue
      do 420, i = 1, dihdat - 2
      do 410, j = i+1, dihdat - 1
      do 400, k = 1, ncoord
      dotprd(i,j) = dotprd(i,j) + codist(k,i) * codist(k,j)
  400 continue
  410 continue
  420 continue
      do 510, i = 1, dihdat-2
      do 500, j = i+1, dihdat - 1
      if (xeqy(atmdst(i), 0.0) .or. xeqy(atmdst(j), 0.0)) then
         dotprd(i,j) = 1.0
      else
         dotprd(i,j) = dotprd(i,j) / atmdst(i) / atmdst(j)
      end if
  500 continue
  510 continue

      if (angnum .gt. ndihed .and. angnum .le. nmnang) then
         cosang =  dotprd(1,3)
      else
         if (xeqy(dotprd(1,2), 1.0) .or. xeqy(dotprd(2,3), 1.0)) then
            cosang = 1.0
         else
            cosang = (dotprd(1,2) * dotprd(2,3) - dotprd(1,3)) /
     +               (sqrt(1.0-dotprd(1,2)**2) * 
     +                sqrt(1.0-dotprd(2,3)**2))
         end if
      end if
      if (abs(cosang) .gt. 1.0) cosang = anint(cosang/abs(cosang))
      ang = acos(cosang) * radian

      if (angnum .le. ndihed .or. angnum .gt. nmnang) then
         detant = codist(1,1) * codist(2,2) * codist(3,3) -
     +            codist(1,1) * codist(2,3) * codist(3,2) +
     +            codist(1,2) * codist(2,3) * codist(3,1) -
     +            codist(1,2) * codist(2,1) * codist(3,3) +
     +            codist(1,3) * codist(2,1) * codist(3,2) -
     +            codist(1,3) * codist(2,2) * codist(3,1)
         If (detant .lt. 0.0) ang = -ang
      end if

      dihed = ang
      END
c****************************************************************************
      Subroutine MkTB(Hbond, SSbond, mcangs, bridpt, chnend, seqlen,
     ~     shtcod)
c
c     routine to establish the turns and bridges
c
      INCLUDE 'p_sstruc.par'

      Integer   Hbond(maxbnd,maxres), seqlen, bridpt(2,maxres),
     +          chnend(maxchn)
      Real      Mcangs(nmnang,maxres)
      Character SSbond(nstruc,maxres)*1
      Integer   strdpl

      Integer   nrules, parlel
      Parameter (nrules = 3, parlel = 1)
      Integer   nxres, nxbond, nxturn, turnsz(nturns), turnpl(nturns),
     +          bsofst(nrules,4), i, j, nxrule, bndi, bndj, bndptr, 
     +          brpl, brdgpt, bsrng(nrules,2), pbrdg(nrules), brpla,
     +          brdg(4,maxres), brptr, brdir, newptr, nxptr, nxpl,
     +          exrule(2,nrules), brptpl, nwptpl, nxstrn, 
     +          lststd, stsht, finsht, nxstr, chcode, stchpl, thstpl,
     +          nxstpl, nxblg, nxbrdg, strpl, brdpt, nxbrpl,
     +          strcd(2,maxres), bstval, shtcod(maxres), fstchn

      Character turnch(nturns)*1, wkch*1, strdch*1, brdch*1,
     +          stchap*(nstchs), stchpa*(nstchs),
     +	        nstchp*(nstchs*2+10)
      logical   xney
c
c     for turn calculations
c
      data turnch / '3', '4', '5'/
      data turnsz / 3, 4, 5/
      data turnpl / thrtnh, alphah, pih/
c
c     for bridge calculations
c
      data ((bsofst(i,j),j = 1,4),i = 1,nrules)
     +     / -1, 0, 0, 1, 0, 0, 0, 0, -1, -1, -2, 1/
      data ((bsrng(i,j),j = 1,2), i = 1,nrules)
     +     / 2, -1, 1, 0, 2, -1/
      data pbrdg / 1, -1, -1/
      data exrule /-1, 1, 1, 1, -1, -1/
c
c     for strand and bridge labelling
c
      data stchap /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data stchpa /'abcdefghijklmnopqrstuvwxyz'/
      data nstchp(1:26) /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/,
     +     nstchp(27:52) /'abcdefghijklmnopqrstuvwxyz'/,
     +     nstchp(53:62) /'1234567890'/

      SAVE turnch, turnsz, turnpl, bsofst, bsrng, pbrdg, stchap, stchpa

      do 110, nxres = 1, seqlen
      do 100, i = 1, nstruc
      ssbond(i,nxres) = space
  100 continue
  110 continue

      fstchn = 1
      do 1000, nxres = 1, seqlen
      if (nxres .gt. chnend(fstchn)) fstchn = fstchn + 1
      do 800, nxbond = 1, 2
      do 600, nxturn = 1, nturns
      if (Hbond(nxbond,nxres) .ne. 0) then
         If (Hbond(nxbond,nxres) - nxres .eq. turnsz(nxturn) .and. 
     +      Hbond(nxbond,nxres) .le. chnend(fstchn)) then
            If (SSbond(turnpl(nxturn),nxres) .eq. bndend) then
               SSbond(turnpl(nxturn),nxres) = bndbth
            else
               SSbond(turnpl(nxturn),nxres) = bndbeg
            end if
            SSbond(turnpl(nxturn),Hbond(nxbond,nxres)) = bndend
            do 400, i = nxres + 1, nxres + turnsz(nxturn) - 1
            If (SSbond(turnpl(nxturn),i) .eq. space)
     +         SSbond(turnpl(nxturn),i) = turnch(nxturn)
  400       continue
            do 500 i = nxres, nxres + turnsz(nxturn)
            SSbond(turn,i) = trnch
  500       continue
         end if
      end if
  600 continue
  800 continue
 1000 continue

      do 2010, nxres = 1,seqlen
      do 2000, nxbond = 1, 4
      brdg(nxbond,nxres) = 0
 2000 continue
 2010 continue

      do 4000, nxrule = 1, nrules
      do 3000, Nxres = bsrng(nxrule,1), seqlen + bsrng(nxrule,2)
      do 2400, bndi = 1, 2
      bndptr = hbond(bndi,nxres+bsofst(nxrule,1))
      if (nxrule .le. parlel .or. bndptr .gt. nxres) then
         brdgpt = bndptr + bsofst(nxrule,2)
         If (bndptr .gt. - bsofst(nxrule,3) .and.
     +      iabs(nxres - brdgpt) .gt. bdgsep)  then
            do 2200, bndj = 1, 2
            If (hbond(bndj,bndptr+bsofst(nxrule,3)) .eq. 
     +         nxres+bsofst(nxrule,4)) then
               brpl = 1
               If (brdg(brpl,nxres) .ne. 0) brpl = 2
               Brdg(brpl,nxres) = pbrdg(nxrule) * brdgpt
               brdg(brpl+2,nxres) = exrule(1,nxrule)
               brpla = 1
               If (brdg(brpla,brdgpt) .ne. 0) brpla = 2
               Brdg(brpla,brdgpt) = pbrdg(nxrule) * nxres
               brdg(brpla+2,brdgpt) = exrule(2,nxrule)
            end if
 2200       continue
         end if
      end if
 2400 continue
 3000 continue
 4000 continue

      do 4500, nxres = 1, seqlen
      if (xney(mcangs(chiral,nxres), null)) then
         ssbond(chisgn,nxres) = '+'
         if (mcangs(chiral,nxres) .lt. 0.0) ssbond(chisgn,nxres) = '-'
      end if
 4500 continue

      do 5100, nxres = 1, seqlen
      shtcod(nxres) = 0
      strcd(1,nxres) = 0
      strcd(2,nxres) = 0
      bridpt(1,nxres) = 0
      bridpt(2,nxres) = 0
 5100 continue
      nxstrn = 0

      do 7000, nxres = 1, seqlen
      do 6900, brpl = 1, 2
      If (brdg(brpl,nxres) .ne. 0 .and. 
     +   iabs(brdg(brpl,nxres)) .gt. nxres) then
         brptr = iabs(brdg(brpl,nxres))
         brdir = brdg(brpl,nxres) / brptr
         brptpl = 2
         if (iabs(brdg(1,brptr)) .eq. nxres) brptpl = 1
         do 5600, nxptr = nxres + 1, min(seqlen,nxres+blgszl)
         do 5500, nxpl = 1,2
         If (brdg(nxpl,nxptr) * brdir .gt. 0) then
            newptr = iabs(brdg(nxpl,nxptr))
c     AVOID DIVIDE BY ZERO 11.3.96
            IF(IABS(NEWPTR-BRPTR).NE.0) THEN
            if ((newptr-brptr)/iabs(newptr-brptr) .eq. brdir) then
               if ((nxptr - nxres .gt. blgszs .and.
     +            iabs(newptr - brptr) .le. blgszs) .OR.
     +            (nxptr - nxres .le. blgszs .and.
     +            iabs(newptr - brptr) .le. blgszl)) then
                  nwptpl = 2
                  if (iabs(brdg(1,newptr)) .eq. nxptr) nwptpl = 1
                  do 5300, i = nxres, nxptr
                  if (ssbond(sheet,i) .eq. space)
     +               ssbond(sheet,i) = shtsym
 5300             continue
                  if (brdg(brpl+2,nxres) .lt. 0)
     +               ssbond(sheet,nxres) = shtsma
                  if (brdg(nxpl+2,nxptr) .lt. 0)
     +               ssbond(sheet,nxptr) = shtsma
                  do 5400, i = brptr, newptr, brdir
                  if (ssbond(sheet,i) .eq. space)
     +               ssbond(sheet,i) = shtsym
 5400             continue
                  if (brdg(brptpl+2,brptr) .lt. 0)
     +               ssbond(sheet,brptr) = shtsma
                  if (brdg(nwptpl+2,newptr) .lt. 0)
     +               ssbond(sheet,newptr) = shtsma
                  if (strcd(brpl,nxres) .eq. 0) then
                     nxstrn = nxstrn + 1
                     strcd(brpl,nxres) = nxstrn * brdir
                     strcd(brptpl,brptr) = nxstrn * brdir
                  end if
                  strcd(nxpl,nxptr) = strcd(brpl,nxres)
                  strcd(nwptpl,newptr) = strcd(brptpl,brptr)
                  go to 5700
               end if
               ENDIF
            end if
         end if
 5500    continue
 5600    continue
 5700    continue
         wkch = pbch
         if (brdir .lt. 0) wkch = apbch
         ssbond(bridge,nxres) = wkch
         ssbond(bridge,brptr) = wkch
      end if
 6900 continue
 7000 continue

      nxres = 1
 7500 continue
      if (ssbond(sheet,nxres) .eq. space) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) goto 10000
         goto 7500
      end if

      stsht = nxres
      nxres = stsht + 1
 7700 continue
      if (ssbond(sheet,nxres) .ne. space) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) goto 7800
         goto 7700
      end if
 7800 continue
      finsht = nxres - 1
      lststd = 0

 8000 continue
      call whstrd(strcd, stsht, finsht, nxstr, nxpl, lststd, bstval)
      if (bstval .ne. 0) then
         lststd = abs(strcd(nxpl,nxstr))
         chcode = mod(lststd,nstchs)
         if (chcode .eq. 0) chcode = nstchs
         if (strcd(nxpl,nxstr) .lt. 0) then
            strdch = stchap(chcode:chcode)
         else
            strdch = stchpa(chcode:chcode)
         end if
         stchpl = bridg1
         if (ssbond(bridg1,nxstr) .ne. space) stchpl = bridg2
         if (stchpl .ne. bridg2) then
            thstpl = nxstr
 8200       continue
            nxstpl = strdpl(thstpl,finsht,strcd(nxpl,nxstr),strcd)
            if (nxstpl .le. 0) goto 8300
            if (ssbond(bridg1,nxstpl) .ne. space) stchpl = bridg2
            thstpl = nxstpl
            if (stchpl .ne. bridg2) goto 8200
         end if
 8300    continue
         ssbond(stchpl,nxstr) = strdch
         bridpt(stchpl-bridg1+1,nxstr) = abs(brdg(nxpl,nxstr))
         thstpl = nxstr
 8500    continue
         nxstpl = strdpl(thstpl,finsht,strcd(nxpl,nxstr),strcd)
         if (nxstpl .le. 0) goto 8700
         do 8600, nxblg = thstpl + 1, nxstpl - 1
         ssbond(stchpl,nxblg) = bulgch
 8600    continue
         ssbond(stchpl,nxstpl) = strdch
         strpl = 1
         if (strcd(1,nxstpl) .ne. strcd(nxpl,nxstr)) strpl = 2
         bridpt(stchpl-bridg1+1,nxstpl) = abs(brdg(strpl,nxstpl))
         thstpl = nxstpl
         go to 8500
 8700    continue
         go to 8000
      end if
      nxres = finsht + 1
      if (nxres .le. seqlen) goto 7500
10000 continue
      If (lststd .gt. nstchs)
     +   Write(errfi,10100) lststd
10100 format('there are ',I3,' strands. Labels have restarted')

      Call shtlab(ssbond, bridpt, shtcod, seqlen)
      do 10500 nxres = 1, seqlen
      if (shtcod(nxres) .ne. 0) then
c        chcode = mod(shtcod(nxres),nstchs)
         chcode = shtcod(nxres)
c        if (chcode .eq. 0) chcode = nstchs
c         if(chcode.ge.0) then
c            ssbond(sheet1,nxres) = nstchp(chcode:chcode)
c         else
c            ssbond(sheet1,nxres)=' '
c         endif
      end if
10500 continue

      nxbrdg = 0
      do 11000, nxres = 1, seqlen
      do 10900, nxpl = 1, 2
      If (brdg(nxpl,nxres) .ne. 0 .and. strcd(nxpl,nxres) .eq. 0
     +   .and. abs(brdg(nxpl,nxres)) .gt. nxres) then
         nxbrdg = nxbrdg + 1
         brdpt = abs(brdg(nxpl,nxres))
         chcode = mod(nxbrdg,nstchs)
         if (chcode .eq. 0) chcode = nstchs
         if (brdg(nxpl,nxres) .lt. 0) then
            brdch = stchap(chcode:chcode)
         else
            brdch = stchpa(chcode:chcode)
         end if
         nxbrpl = bridg1
         if (ssbond(bridg1,nxres) .ne. space) nxbrpl = bridg2
         ssbond(nxbrpl,nxres) = brdch
         bridpt(nxbrpl-bridg1+1,nxres) = brdpt
         nxbrpl = bridg1
         if (ssbond(bridg1,brdpt) .ne. space) nxbrpl = bridg2
         ssbond(nxbrpl,brdpt) = brdch
         bridpt(nxbrpl-bridg1+1,brdpt) = nxres
      end if
10900 continue
11000 continue
      If (nxbrdg .gt. nstchs)
     +  Write(errfi,11100) nxbrdg
11100 format('there are ',I3,' bridges. Labels have restarted')

      END
c***************************************************************************
      Integer Function strdpl(stpl, finpl, thcode, strand)
c
c     find the next occurrence of the code if any in the array
c
      include 'p_sstruc.par'
      integer stpl, finpl, thcode, strand(2,maxres)
      integer nxpl

      do 1000, nxpl = stpl + 1, finpl
      If (strand(1,nxpl) .eq. thcode) go to 1100
      If (strand(2,nxpl) .eq. thcode) go to 1100
 1000 continue
      nxpl = 0
 1100 continue
      strdpl = nxpl
      END
c****************************************************************************
      Subroutine Whstrd(strcd, stsht, finsht, stres, stpl, lststd,
     +                  bstval)
c
c     find the starting position of the lowest numbered strand that
c     is greater than the last
c
      Include 'p_sstruc.par'
      integer strcd(2,maxres), stsht, finsht, stres, stpl, lststd,
     +        bstval

      integer  nxstr, nxpl

      bstval = 0
      do 1010, nxstr = stsht, finsht
      do 1000, nxpl = 1, 2
      If (abs(strcd(nxpl,nxstr)) .gt. lststd) then
         if (bstval .eq. 0) then
            bstval = abs(strcd(nxpl,nxstr))
            stpl = nxpl
            stres = nxstr
         else
            if (abs(strcd(nxpl,nxstr)) .lt. bstval) then
               bstval = abs(strcd(nxpl,nxstr))
               stpl = nxpl
               stres = nxstr
            end if
         end if
      end if
 1000 continue
 1010 continue
      END
c****************************************************************************
      Subroutine Shtlab(ssbond, bridpt, shtcod, seqlen)
c
c     label the sheets using numbers which are turned into
c     alphabetic symbols later
c
      include 'p_sstruc.par'
      character ssbond(nstruc,maxres)*1
      integer   bridpt(2,maxres), shtcod(maxres), seqlen

      integer   nxres, stsht, finsht, shtnum, nxisht, nxpl
      logical   found

      do 100, nxres = 1, seqlen
      shtcod(maxres) = 0
  100 continue

      shtnum = 0
      nxres = 1
 1000 continue
      call fdshus(nxres, ssbond, shtcod, stsht, finsht, seqlen)
      if (stsht .ne. 0) then
         shtnum = shtnum + 1
         do 1500, nxisht = stsht, finsht
         shtcod(nxisht) = shtnum
         do 1400, nxpl = 1,2
         If (bridpt(nxpl,nxisht) .ne. 0) then
            if (shtcod(bridpt(nxpl,nxisht)) .eq. 0)
     +         call setsht(bridpt(nxpl,nxisht), shtnum, shtcod, 
     +                     ssbond, seqlen)
         end if
 1400    continue
 1500    continue

 2000    continue
         call fdshps(finsht+1, shtnum, ssbond, shtcod, bridpt, found,
     +               seqlen)
         if (found) go to 2000

         nxres = finsht + 1
         if (nxres .le. seqlen) go to 1000
      end if
      If (shtnum .gt. nstchs)
     +   Write(errfi,2100) shtnum
 2100 format('there are ',I3,' sheets. Labels have restarted')
      END
c*****************************************************************************
      Subroutine fdshus(stpt, ssbond, shtcod, stsht, finsht, seqlen)
c
c     find the first unset sheet starting at stpt
c
      include 'p_sstruc.par'
      integer   stpt, stsht, finsht, shtcod(maxres), seqlen
      character ssbond(nstruc,maxres)*1

      integer   nxres

      stsht = 0
      nxres = stpt
 1000 continue
      if (ssbond(sheet,nxres) .eq. space .or. shtcod(nxres) .ne. 0) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) go to 2000
         go to 1000
      end if
      stsht = nxres
 1500 continue
      nxres = nxres + 1
      if (nxres .le. seqlen) then
         if (ssbond(sheet,nxres) .ne. space) goto 1500
      end if
      finsht = nxres - 1
 2000 continue
      END
c****************************************************************************
      subroutine fdshps(stpt, shtnum, ssbond, shtcod, bridpt, found,
     +                  seqlen)
c
c     find the next partially set sheet, if any (neg shtnum), set it
c     to fully set and then set its partners to partly set if not
c     already set
c
      include 'p_sstruc.par'
      integer   stpt, shtnum, shtcod(maxres), bridpt(2,maxres), seqlen
      character ssbond(nstruc,maxres)
      logical   found

      integer   nxres, stsht, finsht, nxpl

      found = .false.
      nxres = stpt
 1000 continue
      if (shtcod(nxres) .ne. - shtnum) then
         nxres = nxres + 1
         if (nxres .gt. seqlen) goto 5000
         goto 1000
      end if
      stsht = nxres
      found = .true.
 1500 continue
      if (shtcod(nxres) .eq. -shtnum) then
         nxres = nxres + 1
         if (nxres .le. seqlen) goto 1500
      end if
      finsht = nxres - 1

      do 2000, nxres = stsht, finsht
      shtcod(nxres) = shtnum
      do 1900, nxpl = 1,2
      if (bridpt(nxpl,nxres) .ne. 0) then
         if (shtcod(bridpt(nxpl,nxres)) .eq. 0)
     +      call setsht(bridpt(nxpl,nxres), shtnum, shtcod, ssbond, 
     +                  seqlen)
      end if
 1900 continue
 2000 continue
 5000 continue
      END
c*****************************************************************************
      Subroutine Setsht(refpl, shtnum, shtcod, ssbond, seqlen)
c
c     finds the limits of the super sheet containing refpl and
c     sets this to partially set
c
      include 'p_sstruc.par'
      integer   refpl, shtnum, shtcod(maxres), seqlen
      character ssbond(nstruc,maxres)*1

      integer   nxres, stsht, finsht

      nxres = refpl
  200 continue
      if (ssbond(sheet,nxres) .ne. space) then
	 nxres = nxres + 1
         if (nxres .le. seqlen) go to 200
      end if
      finsht = nxres - 1
      nxres = refpl
  400 continue
      if (ssbond(sheet,nxres) .ne. space) then
         nxres = nxres - 1
         if (nxres .ge. 1) go to 400
      end if
      stsht = nxres + 1

      do 800, nxres = stsht, finsht
      shtcod(nxres) = -shtnum
  800 continue
      END
c******************************************************************************
      Subroutine MkBend(mcangs, ssbond, seqlen)
c
c     routine to determine which main chain residues are bent at over
c     the requisite angle
c
      INCLUDE 'p_sstruc.par'

      Real      mcangs(nmnang,maxres)
      Character Ssbond(nstruc,maxres)*1
      Integer   seqlen
c
      real      bendsz
      parameter (bendsz = 70.0)
      Integer   nxres
      logical   xney

      do 1000, nxres = 1, seqlen
      if (mcangs(kappa,nxres) .gt. bendsz .and.
     +   xney(mcangs(kappa,nxres), null)) ssbond(bend,nxres) = bendch
 1000 continue
      END
c******************************************************************************
      Subroutine MkSumm(ssbond, sssumm, summpl, seqlen)
c
c     generate a structure summary using the priority rules
c
      include 'p_sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), summpl(maxres)
      integer     seqlen

      character   stsym(nstruc)*1, altsym(nstruc)*1
      data stsym /'H', ' ', 'E', 'B', ' ', ' ', 'G', 'I', 'T', 'S', '+'/
      data altsym/'h', ' ', 'e', 'b', ' ', ' ', 'g', 'i', 't', 's', '+'/
      SAVE stsym, altsym

      Integer     nxres

      do 200, nxres = 1, seqlen
      SSsumm(nxres) = space
      summpl(nxres) = space
  200 continue
      call hlxset(ssbond, sssumm, summpl, alphah, ahsz, stsym(alphah),
     +            altsym(alphah), seqlen)
      do 400, nxres = 1, seqlen
      If (ssbond(sheet,nxres) .ne. space) then
         if (sssumm(nxres) .eq. space) sssumm(nxres) = stsym(sheet)
         if (summpl(nxres) .eq. space .or. summpl(nxres) .eq.
     +      altsym(alphah) .or. summpl(nxres) .eq. altsym(sheet))
     +      summpl(nxres) = stsym(sheet)
         if (ssbond(sheet,nxres) .eq. shtsma) then
            if (summpl(nxres-1) .eq. space .or. summpl(nxres-1) .eq.
     +         stsym(bridge)) summpl(nxres-1) = altsym(sheet)
            if (summpl(nxres+1) .eq. space .or. summpl(nxres-1) .eq.
     +         stsym(bridge)) summpl(nxres+1) = altsym(sheet)
         end if
      end if
      If (ssbond(bridge,nxres) .ne. space) then
         if (sssumm(nxres) .eq. space) sssumm(nxres) = stsym(bridge)
         if (summpl(nxres) .eq. space) summpl(nxres) = stsym(bridge)
      end if
  400 continue
      call hlxset(ssbond, sssumm, summpl, thrtnh, thrtsz, stsym(thrtnh),
     +            altsym(thrtnh), seqlen)
      call hlxset(ssbond, sssumm, summpl, pih, pihsz, stsym(pih),
     +            altsym(pih), seqlen)
      call hlxres(sssumm, stsym(thrtnh), altsym(thrtnh), stsym(turn),
     +            altsym(turn), thrtsz, seqlen)
      call hlxres(summpl, stsym(thrtnh), altsym(thrtnh), stsym(turn),
     +            altsym(turn), thrtsz, seqlen)
      call hlxres(sssumm, stsym(pih), altsym(pih), stsym(turn),
     +            altsym(turn), pihsz, seqlen)
      call hlxres(summpl, stsym(pih), altsym(pih), stsym(turn),
     +            altsym(turn), pihsz, seqlen)
      call turnst(ssbond, sssumm, summpl, stsym(turn), altsym(turn),
     +            seqlen)
      call bendst(ssbond, sssumm, summpl, seqlen)

      END
c*****************************************************************************
      Subroutine Hlxset(ssbond, sssumm, summpl, hlxpos, hlxsz, hlxch,
     +                  altch,  seqlen)
c
c     go through the bond tables and note any helices into the
c     summary line
c
      include 'p_sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), hlxch, altch,
     +            summpl(maxres)
      integer     hlxpos, hlxsz, seqlen

      integer     nxres, nxsym

      do 1000, nxres = 2, seqlen - hlxsz
      if (ssbond(hlxpos,nxres) .eq. bndbeg .or.
     +   ssbond(hlxpos,nxres) .eq. bndbth) then
         if (ssbond(hlxpos,nxres-1) .eq. bndbeg .or.
     +      ssbond(hlxpos,nxres-1) .eq. bndbth) then
            do 500, nxsym = 0, hlxsz - 1
            if (sssumm(nxres+nxsym) .eq. space)
     +         sssumm(nxres+nxsym) = hlxch
            if (summpl(nxres+nxsym) .eq. space .OR. 
     +         summpl(nxres+nxsym) .eq. altch)
     +         summpl(nxres+nxsym) = hlxch
  500       continue
            if (summpl(nxres-1) .eq. space) summpl(nxres-1) = altch
            if (summpl(nxres+hlxsz) .eq. space)
     +             summpl(nxres+hlxsz) = altch
         end if
      end if
 1000 continue
      END
c*****************************************************************************
      Subroutine hlxres(summ, hlxch, altch, turnch, alturn, hlxsz,
     +                  seqlen)
c
c     to remove helix character runs of less than minimal helix size
c
      include 'p_sstruc.par'
      character*1 summ(maxres), hlxch, turnch, altch, alturn
      integer     hlxsz, seqlen

      integer     nxres, stpos, nxsym

      stpos = 0
      do 1000, nxres = 1, seqlen
      if (summ(nxres) .eq. hlxch .or. summ(nxres) .eq. altch) then
         If (stpos .eq. 0) stpos = nxres
      else
         if (stpos .ne. 0) then
            if (nxres - stpos .lt. hlxsz) then
               do 500, nxsym = stpos, nxres - 1
               if (summ(nxsym) .eq. hlxch) then
                  summ(nxsym) = turnch
               else
                  summ(nxsym) = alturn
               end if
  500          continue
            end if
            stpos = 0
         end if
      end if
 1000 continue
      if (stpos .ne. 0) then
         if (seqlen - stpos + 1 .lt. hlxsz) then
            do 1500, nxsym = stpos, seqlen
            if (summ(nxres) .eq. hlxch) then
               summ(nxsym) = turnch
            else
               summ(nxsym) = alturn
            end if
 1500       continue
         end if
      end if
      END
c**************************************************************************
      Subroutine turnst(ssbond, sssumm, summpl, turnch, alturn, seqlen)
c
c     search for solitary turn type bonds and set in the symbol
c
      include 'p_sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), turnch,
     +            summpl(maxres), alturn
      integer     seqlen

      integer     nxres, nxsym, nxhlx, hlxpri(nhelix), hlxsz(nhelix)

      hlxpri(1) = alphah
      hlxpri(2) = thrtnh
      hlxpri(3) = pih
      hlxsz(1) = ahsz
      hlxsz(2) = thrtsz
      hlxsz(3) = pihsz

      do 2000, nxhlx = 1, nhelix
      If (ssbond(hlxpri(nxhlx),1) .eq. bndbeg .and.
     +   ssbond(hlxpri(nxhlx),2) .ne. bndbeg) then
         do 200, nxsym = 2, hlxsz(nxhlx)
         If (sssumm(nxsym) .eq. space) sssumm(nxsym) = turnch
         If (summpl(nxsym) .eq. space .or. summpl(nxsym) .eq.
     +      alturn) summpl(nxsym) = turnch
  200    continue
         If (summpl(1) .eq. space) summpl(1) = alturn
         If (summpl(hlxsz(nxhlx)+1) .eq. space)
     +      summpl(hlxsz(nxhlx)+1) = alturn
      end if
      do 1000, nxres = 2, seqlen - hlxsz(nxhlx)
      If ((ssbond(hlxpri(nxhlx),nxres) .eq. bndbeg .or.
     +   ssbond(hlxpri(nxhlx),nxres) .eq. bndbth) .AND.
     +   (ssbond(hlxpri(nxhlx),nxres-1) .ne. bndbeg .and.
     +   ssbond(hlxpri(nxhlx),nxres-1) .ne. bndbth) .AND.
     +   (ssbond(hlxpri(nxhlx),nxres+1) .ne. bndbeg .and.
     +   ssbond(hlxpri(nxhlx),nxres+1) .ne. bndbth)) then
         do 500, nxsym = nxres + 1, Nxres + hlxsz(nxhlx) - 1
         If (sssumm(nxsym) .eq. space) sssumm(nxsym) = turnch
         If (summpl(nxsym) .eq. space .or. summpl(nxsym) .eq.
     +      alturn) summpl(nxsym) = turnch
  500    continue
         If (summpl(nxres) .eq. space) summpl(nxres) = alturn
         If (summpl(nxres+hlxsz(nxhlx)) .eq. space)
     +      summpl(nxres+hlxsz(nxhlx)) = alturn
      end if
 1000 continue
 2000 continue
      END
c*****************************************************************************
      Subroutine Bendst(ssbond, sssumm, summpl, seqlen)
c
c     enter bend symbols into summary line
c
      Include 'p_sstruc.par'
      character*1 ssbond(nstruc,maxres), sssumm(maxres), summpl(maxres)
      integer     seqlen
      Integer     nxres
c
      do 1000, nxres = 1, seqlen
      If (ssbond(bend,nxres) .ne. space) then
         if (sssumm(nxres) .eq. space) sssumm(nxres) =
     +      ssbond(bend,nxres)
         If (summpl(nxres) .eq. space) summpl(nxres) =
     +      ssbond(bend,nxres)
      end if
 1000 continue
      END
c****************************************************************************
      Subroutine mkauth(seqbcd, autsum, ssord, astruc, seqlen)
c
c     reread the author structure records and create the structure
c
      include 'p_sstruc.par'
      character seqbcd(maxres)*6, autsum(maxres)*1, ssord*21
      logical   astruc
      integer   seqlen
      integer   respl

      integer   nxsym, stpos, finpos, lstpos, hlxno, hlxcls, nxres,
     +          strno, sense, strpri, hlxcds
      parameter (hlxcds = 10)
      character stres*6, finres*6, autrec*80, strkch*1, shtid*3,
     +          hlxch(hlxcds)*1
      data hlxch/'H','J','I','K','G','L','M','N','O','P'/
      save hlxch

      do 100, nxres = 1, seqlen
      autsum(nxres) = space
  100 continue

      if (astruc) then
         Rewind(unit=authfi)
         lstpos = 0
  200    continue
         read(authfi,'(a)',end=2000) autrec
         if (autrec(1:keylen) .eq. hlxkey) then
            read(autrec,400) hlxno, stres(1:1), stres(2:6), finres(1:1),
     +          finres(2:6), hlxcls
  400       format(7x,i3,9x,a1,1x,a5,5x,a1,1x,a5,i2)
            strkch = hlxch(1)
            if (hlxcls .le. hlxcds .and. hlxcls .gt. 0) 
     +         strkch = hlxch(hlxcls)
         else
            if (autrec(1:keylen) .eq. shtkey) then
               read(autrec,600) strno, shtid, stres, finres, sense
  600          format(7x,i3,1x,a3,7x,a6,5x,a6,i2)
               strkch = shtsym
            else
               if (autrec(1:keylen) .eq. trnkey) then
                  read(autrec,800) stres, finres
  800             format(19x,a6,5x,a6)
                  strkch = trnch
               else
                  go to 200
               end if
            end if
         end if
         stpos = respl(seqbcd, stres, lstpos, seqlen)
         finpos = respl(seqbcd, finres, lstpos, seqlen)
         If (finpos .lt. 0 .or. stpos .lt. 0 .or.
     +      finpos .lt. stpos .or. (strkch .eq. hlxch(1) .and. 
     +      (hlxcls .le. 0 .or. hlxcls .gt. hlxcds))) then
            Write(errfi,*)'structure record in error'
            Write(errfi,*)autrec
         else
            strpri = index(ssord,strkch)
            do 1000, nxsym = stpos, finpos
            if (index(ssord,autsum(nxsym)) .gt. strpri) 
     +         autsum(nxsym) = strkch
 1000       continue
         end if
         go to 200
 2000    continue
      end if
      end
c*****************************************************************************
      Subroutine disulf(seqbcd, dsdsts, scangs, scaaat, scatin, atomin,
     +                  mcaaat, seqcod, dsdef, dspart,  astruc, seqlen)
c
c     read the author file and extract any author definitions of
c     disulphide bridges. check for possible bridges on a distance
c     criterion as well.
c
      include 'p_sstruc.par'
      character seqbcd(maxres)*6, dsdef(maxres)*1
      real      scangs(nsdang,maxres), dsdsts(3,maxres),
     +          scaaat(ncoord,sdchna,maxres), 
     +          mcaaat(ncoord,mnchna,maxres)
      integer   dspart(maxres), seqcod(maxres), seqlen
      logical   astruc, scatin(sdchna,maxres), atomin(mnchna,maxres)
      real      atdist
      integer   respl

      character autdef*1, disdef*1
      parameter (autdef = 'A', disdef = 'D')
      integer   nxdist, nxres, pos1, pos2, errset, partpl, nxcys
      real      partds, ssdist
      logical   cystyp, xney
      character autrec*80, code1*6, code2*6, typefl*1

      do 110, nxres = 1, seqlen
      dsdef(nxres) = space
      dspart(nxres) = 0
      do 100, nxdist = 1, 3
      dsdsts(nxdist,nxres) = 0.0
  100 continue
  110 continue

      If (astruc) then
         rewind(unit=authfi)
         pos1 = 0
  200    continue
         if (pos1 .lt. 0) pos1 = 0
         read (authfi,'(A)', end = 2000) autrec
         if (autrec(1:keylen) .ne. dsfkey) go to 200
         read(autrec,400) code1(1:1), code1(2:6), code2(1:1), code2(2:6)
  400    format(15x,a1,1x,a5,7x,a1,1x,a5)
         pos1 = respl(seqbcd, code1, pos1, seqlen)
         errset = 0
         if (pos1 .gt. 0) pos2 = respl(seqbcd, code2, pos1, seqlen)
         If (pos1 .gt. 0 .and. pos2 .gt. 0) then
            if (Cystyp(seqcod(pos1)) .and. cystyp(seqcod(pos2))) then
               typefl = autdef
               call dsset(pos1, pos2, dsdsts, scangs, scaaat, scatin,
     +                    atomin, mcaaat, dsdef, dspart, typefl)
            else
               errset = 1
            end if
         else
            errset = 1
         end if
         If (errset .eq. 1) then
            write(errfi,'(A)') 'error in SSBOND record'
            write(errfi,'(A)') autrec
         end if
         go to 200

 2000    continue
         close(unit = authfi)
      end if

      nxres = 1
 2500 continue
      if (cystyp(seqcod(nxres))) then
         if (scatin(sgamma,nxres)) then
            partpl = 0
            partds = 0.0
            do 3000, nxcys = nxres + 1, seqlen
            if (cystyp(seqcod(nxcys))) then
               if (scatin(sgamma,nxcys)) then
                  ssdist = atdist(scaaat(1,sgamma,nxres),
     +                     scaaat(1,sgamma,nxcys))
                  if (ssdist .lt. sssep) then
                     if (xney(partds, 0.0)) then
                        if (partds .lt. ssdist) then
                           Write(errfi,2700) nxres, nxcys, ssdist
                        else
                           Write(errfi,2700) nxres, partpl, partds
 2700                      format('multifurcation of disulphides ',
     +                            I4, 1x, i4, 1x, f3.1)
                           partpl = nxcys
                           partds = ssdist
                        end if
                     else
                        partpl = nxcys
                        partds = ssdist
                     end if
                  end if
               end if
            end if
 3000       continue
            if (partpl .ne. 0) then
               typefl = disdef
               call dsset(nxres, partpl, dsdsts, scangs, scaaat, scatin,
     +                    atomin, mcaaat, dsdef, dspart, typefl)
            end if
         end if
      end if
      nxres = nxres + 1
      if (nxres .le. seqlen) go to 2500
      END
c****************************************************************************
      logical function cystyp(seqcod)
c
c     check that we do have a potential disulfide former
c
      include 'p_sstruc.par'
      integer  seqcod

      If (seqcod .eq. cyscod .or. seqcod .eq. mprcod) then
         cystyp = .true.
      else
         cystyp = .false.
      end if
      end
c*****************************************************************************
      Subroutine dsset(part1,  part2,  dsdsts, scangs, scaaat, scatin,
     +                 atomin, mcaaat, dsdef,  dspart, typefl)
c
c     calculate all the features for a disulphide bridge. that is 3 chi
c     angles, calpha, cbeta, and ss distances. set in the partner and
c     definition fields. Check for conflict with author defined records
c     or previously defined bridges
c
      include 'p_sstruc.par'
      integer   part1, part2, dspart(maxres)
      real      dsdsts(3,maxres), scangs(nsdang,maxres), 
     +          scaaat(ncoord,sdchna,maxres), 
     +          mcaaat(ncoord,mnchna,maxres)
      logical   scatin(sdchna,maxres), atomin(mnchna,maxres)
      character dsdef(maxres)*1, typefl*1
      real      atdist, dihed
  
      character*1 botdef, autdef
      real        mxdis
      parameter  (botdef = 'B', autdef = 'A', mxdis = 9.9)
      integer     nxdis

      If ((dspart(part1) .ne. 0 .and. dspart(part1) .ne. part2) .or.
     +   (dspart(part1) .eq. 0 .and. dspart(part2) .ne. 0)) then
         if (dsdef(part1) .ne. space .or. dsdef(part2) .ne. space) then
            write(errfi,100) part1, part2
  100       format('conflict with author defined disulphide -',
     +             ' abandonning ',I4,1x,I4)
         else
            write(errfi,200) part1, part2
  200       format('conflict with previously defined disulphide -',
     +             ' abandoning ',I4,1x,I4)
         end if
      else
         if (dsdef(part1) .eq. autdef) then
            dsdef(part1) = botdef
            dsdef(part2) = botdef
         else
            dspart(part1) = part2
            dspart(part2) = part1
            dsdef(part1) = typefl
            dsdef(part2) = typefl
            if (atomin(calph,part1) .and. scatin(cbeta,part1) .and. 
     +         scatin(sgamma,part1) .and. scatin(sgamma,part2))
     +         scangs(2,part1) = dihed(nmnang+2,
     +         mcaaat(1,calph,part1), scaaat(1,cbeta,part1),
     +         scaaat(1,sgamma,part1), scaaat(1,sgamma,part2))
            if (scatin(cbeta,part1) .and. scatin(sgamma,part1) .and. 
     +         scatin(sgamma,part2) .and. scatin(cbeta,part2)) then
               scangs(3,part1) = dihed(nmnang+3,
     +         scaaat(1,cbeta,part1), scaaat(1,sgamma,part1), 
     +         scaaat(1,sgamma,part2), scaaat(1,cbeta,part2))
               scangs(3,part2) = scangs(3,part1)
            end if
            if (atomin(calph,part2) .and. scatin(cbeta,part2) .and. 
     +         scatin(sgamma,part2) .and. scatin(sgamma,part1))
     +         scangs(2,part2) = dihed(nmnang+2,
     +         mcaaat(1,calph,part2), scaaat(1,cbeta,part2),
     +         scaaat(1,sgamma,part2), scaaat(1,sgamma,part1))
            do 1000, nxdis = 1, 3
            dsdsts(nxdis,part1) = -1.0
 1000       continue
            If (atomin(calph,part1) .and. atomin(calph,part2))
     +         dsdsts(1,part1) = atdist(mcaaat(1,calph,part1),
     +         mcaaat(1,calph,part2))
            If (scatin(cbeta,part1) .and. atomin(cbeta,part2))
     +         dsdsts(2,part1) = atdist(scaaat(1,cbeta,part1),
     +         scaaat(1,cbeta,part2))
            If (scatin(sgamma,part1) .and. scatin(sgamma,part2))
     +         dsdsts(3,part1) = atdist(scaaat(1,sgamma,part1),
     +         scaaat(1,sgamma,part2))
            do 1200, nxdis = 1, 3
            if (dsdsts(nxdis,part1) .gt. mxdis) then
               write(errfi,1100) part1, part2, dsdsts(nxdis,part1)
 1100          format('disulphide distance too large for',
     +                I4,1x,I4,1x,f5.1)
               dsdsts(nxdis,part1) = mxdis
            end if
 1200       continue
            do 1500, nxdis = 1,3
            dsdsts(nxdis,part2) = dsdsts(nxdis,part1)
 1500       continue
         end if
      end if
      END
c*****************************************************************************
      Integer function respl(seqbcd, theres, lstpos, seqlen)
c
c     search through the residue identifiers to find the position
c     of the relevant code. Codes will normally be in order so the
c     search starts from the last position.
c
      include 'p_sstruc.par'
      character seqbcd(maxres)*6, theres*6
      integer   lstpos, seqlen

      integer   nxres

      do 100, nxres = lstpos + 1, seqlen
      if (seqbcd(nxres) .eq. theres) go to 300
  100 continue
      do 200, nxres = 1, lstpos
      if (seqbcd(nxres) .eq. theres) go to 300
  200 continue
      nxres = -1
  300 continue
      respl = nxres
      end
c*****************************************************************************
      Subroutine result(hbonde, hbond,  aacod3, aacod1, seqcod, seqbcd, 
     +                  brksym, ssbond, bridpt, sssumm, summpl, autsum, 
     +                  mcangs, altcod, aastd,  atpres, thmprs, scangs,
     +                  cadsts, brknam, headln, hedlen, outsty, dsdsts,
     +                  dsdef,  dspart, hetcod, oois,   seqlen, shtcod)
c
c     Results are in 4 parts. The standard output set of details
c     one line to an amino acid. The second set makes the full set
c     again one line per amino acid. The third set is the structure
c     patterns and summaries at 100 amino acids to the line
c     with finally error and information messages
c
      include 'p_sstruc.par'
      integer   hbond(maxbnd,maxres), seqcod(maxres), seqlen,
     +          bridpt(2,maxres), hedlen, dspart(maxres), 
     ~          oois(2,maxres),shtcod(maxres)
      real      hbonde(maxbnd,maxres), mcangs(nmnang,maxres),
     +          scangs(nsdang,maxres), cadsts(11,maxres), 
     +          dsdsts(3,maxres)
      character aacod3*(abetsz*3), ssbond(nstruc,maxres)*1,
     +          sssumm(maxres)*1, brknam*(fnamln), summpl(maxres),
     +          seqbcd(maxres)*6, aacod1*(abetsz), autsum(maxres)*1,
     +          brksym(maxres)*1, headln*132, altcod(maxres)*1,
     +          aastd(maxres)*1, atpres(3,maxres), thmprs(maxres),
     +          outsty*1, dsdef(maxres)*1, hetcod*102

      integer   nxres, i, bpl(maxbnd), staa, nxstr, nxang, nxst, nxrow,
     +          aathln, strpts(9), j, hednum
      character aa3*3, aa1*1, lefthd(13)*9, rthd(13)*9,
     +          wkln1*(aaonln), wkln2*122, sumln(3)*(aaonln)
      data bpl /3,1,4,2/
      data strpts/chisgn, bend, turn, pih, thrtnh, bridg2, bridg1,
     +            sheet1, alphah/
      data lefthd /'chirality', '    bends', '    turns', '  5-turns',
     +             '  3-turns', ' bridge-2', ' bridge-1', '   sheets',
     +             '  4-turns', '   author', 'Kabs/Sand', '  summary',
     +             ' sequence'/
      data rthd   /'chirality', 'bends    ', 'turns    ', '5-turns  ',
     +             '3-turns  ', 'bridge-2 ', 'bridge-1 ', 'sheets   ',
     +             '4-turns  ', 'author   ', 'Kabs/Sand', 'summary  ',
     +             'sequence '/
      save bpl, strpts, lefthd, rthd

      hednum = 1
      call hedset(brknam, headln, hedlen, hednum, seqlen)
      do 500, nxres = 1, seqlen
      if (seqcod(nxres) .le. abetsz) then
         staa = (seqcod(nxres) - 1) * 3 + 1
         aa3 = aacod3(staa:staa+2)
         aa1 = aacod1(seqcod(nxres):seqcod(nxres))
      else
         staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
         aa3 = hetcod(staa:staa+2)
         aa1 = '-'
      end if
c      write(outfi,300) nxres, brksym(nxres), seqbcd(nxres),
c     +     altcod(nxres), aa3, aa1, autsum(nxres), sssumm(nxres),
c     +     summpl(nxres), (ssbond(nxstr,nxres), nxstr = alphah, sheet1),
c     +     (ssbond(nxstr,nxres), nxstr = bridg1, chisgn),
c     +     (bridpt(nxstr,nxres), nxstr = 1, 2),
c     +     (mcangs(nxang,nxres), nxang = 1, nmnang),
c     +     (hbond(bpl(i),nxres), hbonde(bpl(i),nxres), i=1,4),
c     +     oois(ooi1,nxres), oois(ooi2,nxres)
c  300 format(1x,i4,a1,a6,x,a1,x,a3,x,a1,3(x,a1),x,9a1,2(x,i4),
c     +       4(x,f6.1),2(x,f5.1),4(x,i4,x,f4.1),2(x,i2))
      write(outfi,301) nxres,brksym(nxres),seqbcd(nxres),altcod(nxres),
     ~     aa3,aa1,autsum(nxres),sssumm(nxres),summpl(nxres),
     ~     ssbond(1,nxres),shtcod(nxres),
     ~     (ssbond(nxstr,nxres),nxstr=bridg1,chisgn),
     +     (bridpt(nxstr,nxres), nxstr = 1, 2),
     +     (mcangs(nxang,nxres), nxang = 1, nmnang),
     +     (hbond(bpl(i),nxres), hbonde(bpl(i),nxres), i=1,4),
     +     oois(ooi1,nxres), oois(ooi2,nxres)
 301  format(i5,a1,a6,1x,a1,1x,a3,1x,a1,3(1x,a1),1x,a1,i3,7a1,2(i5),
     ~     4(1x,f6.1),2(1x,f5.1),4(i5,1x,f4.1),2(1x,i2))
  500 continue
      write(outfi,*)char(12)

      if (outsty .eq. 'F') then
         hednum = 2
         call hedset(brknam, headln, hedlen, hednum, seqlen)
         do 1000, nxres = 1, seqlen
         if (seqcod(nxres) .le. abetsz) then
            staa = (seqcod(nxres) - 1) * 3 + 1
            aa3 = aacod3(staa:staa+2)
            aa1 = aacod1(seqcod(nxres):seqcod(nxres))
         else
            staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
            aa3 = hetcod(staa:staa+2)
            aa1 = '-'
         end if
         write(outfi,900)nxres, brksym(nxres), seqbcd(nxres), aa3, aa1,
     +      aastd(nxres), (atpres(i,nxres),i=1,3),
     +      thmprs(nxres), (scangs(i,nxres),i=1,4), dspart(nxres), 
     +      dsdef(nxres), (dsdsts(i,nxres),i=1,3), 
     +      (cadsts(i,nxres),i=1,11)
  900    format(1x,i4,a1,a6,1x,a3,1x,a1,5(1x,a1),4(1x,f6.1),1x,i4,
     +        1x,a1,3(1x,f4.1),11(1x,f4.1))
 1000    continue
         write(outfi,*) char(12)
      end if

      write(outfi,1001) brknam
 1001 format(1x,a)
      write(outfi,1001) headln(1:hedlen)
      do 2000, nxst = 1, seqlen, aaonln
      write(outfi,*)
      aathln = aaonln
      if (nxst + aathln - 1 .gt. seqlen) aathln = seqlen - nxst + 1
      do 1100, nxstr = 1,aathln
      sumln(1)(nxstr:nxstr) = autsum(nxst+nxstr-1)
      sumln(2)(nxstr:nxstr) = sssumm(nxst+nxstr-1)
 1100 continue
      do 1300, nxrow = 1, 2
      wkln2 = lefthd(9+nxrow) // '  ' // sumln(nxrow)(1:aathln) //
     +        '  ' // rthd(9+nxrow)
      write(outfi,1200)wkln2(1:aathln+22)
 1200 format(1x,A)
 1300 continue
      do 1500, nxrow = 1, 9
      do 1400, nxstr = 1, aathln
      wkln1(nxstr:nxstr) = ssbond(strpts(nxrow),nxst+nxstr-1)
 1400 continue
      wkln2 = lefthd(nxrow) // '  ' // wkln1(1:aathln) // '  '
     +        // rthd(nxrow)
      write(outfi,1200)wkln2(1:aathln+22)
 1500 continue
      do 1600, nxstr = 1,aathln
      sumln(3)(nxstr:nxstr) = summpl(nxst+nxstr-1)
      if (seqcod(nxst+nxstr-1) .le. abetsz) then
         wkln1(nxstr:nxstr) = 
     +   aacod1(seqcod(nxst+nxstr-1):seqcod(nxst+nxstr-1))
      else
         wkln1(nxstr:nxstr) = '-'
      end if
 1600 continue
      wkln2 = lefthd(12) // '  ' // sumln(3)(1:aathln) // '  '
     +        // rthd(12)
      write(outfi,1200) wkln2(1:aathln+22)
      wkln2 = lefthd(13) // '  ' // wkln1(1:aathln) // '  ' //rthd(13)
      write(outfi,1200) wkln2(1:aathln+22)
      write(outfi,1700) (nxst-1+j*10,j=1,aathln/10)
 1700 format(12x,10(i10))
 2000 continue

      call msgprt
      write(outfi,*) char(12)
      end
c*****************************************************************************
      subroutine carslt(aacod3, aacod1, seqcod, seqbcd, mcangs, cadsts,
     +                  brksym, altcod, aastd,  brknam, headln, hedlen,
     +                  autsum, hetcod, oois,   seqlen)
c
c     results for c-alpha only files. these are placed in a differently
c     named version of the output file
c
      include 'p_sstruc.par'
      character aacod3*(abetsz*3), aacod1*(abetsz), seqbcd(maxres)*6,
     +          altcod(maxres)*1, aastd(maxres)*1, headln*132,
     +          brknam*(fnamln), brksym(maxres)*1, autsum(maxres)*1,
     +          hetcod*102
      integer   seqcod(maxres), hedlen, seqlen, oois(2,maxres)
      real      mcangs(nmnang,maxres), cadsts(11,maxres)

      integer   hednum, nxres, staa, nxang, nxdis
      character aa3*3, aa1*1

      hednum = 3
      call hedset(brknam, headln, hedlen, hednum, seqlen)
      do 1000, nxres = 1, seqlen
      if (seqcod(nxres) .le. abetsz) then
         staa = (seqcod(nxres) - 1) * 3 + 1
         aa3 = aacod3(staa:staa+2)
         aa1 = aacod1(seqcod(nxres):seqcod(nxres))
      else
         staa = (seqcod(nxres) - abetsz - 1) * 3 + 1
         aa3 = hetcod(staa:staa+2)
         aa1 = '-'
      end if
      write(outfi,300) nxres, brksym(nxres), seqbcd(nxres),
     +      altcod(nxres), aa3, aa1, aastd(nxres), autsum(nxres),
     +      (mcangs(nxang,nxres), nxang = chiral,kappa),
     +      (cadsts(nxdis,nxres), nxdis = 1,11),
     +      oois(ooi1,nxres), oois(ooi2,nxres)
  300 format(1x,i4,a1,a6,1x,a1,1x,a3,3(1x,a1),1x,f6.1,1x,f5.1,
     +       11(1x,f4.1),2(1x,i2))
 1000 continue
      call msgprt
      write(outfi,*) char(12)
      end
c*****************************************************************************
      SUBROUTINE msgprt
c
c     write the error messages to the output file
c
      include 'p_sstruc.par'
      character wkrec*80
      integer   errflg

      rewind(unit=errfi)
      errflg = 0
  200 continue
      read(errfi,300,end=1000) wkrec
  300 format(A)
      if (errflg .eq. 0) then
         write(outfi,500)'Messages'
  500    format(////,1x,A)
      end if
      errflg = 1
      write(outfi,*) wkrec
c      write(*,*) wkrec
      go to 200
 1000 continue
      close(unit=errfi)
      END
c****************************************************************************
      Subroutine hedset(brknam, headln, hedlen, hednum, seqlen)
c
c     produce headings from results parts 1 and 2 and carslt
c
      include 'p_sstruc.par'
      character brknam*(fnamln), headln*132
      Integer   hedlen, hednum, seqlen

      write(outfi,101)
 101  format(' Secondary structure calculation program - ',
     +              'copyright by David Keith Smith, 1989')
      write(outfi,100) brknam
  100 format(1x,a)
      write(outfi,102) headln(1:hedlen)
 102  format(1x,a)
      write(outfi,120) seqlen
  120 format(1x,'Sequence length -',I5)
 
      if (hednum .eq. 1) then
         write(outfi,1)
 1       format
     +  ('             A       A K K                                  ',
     +   '                                      hydrogen bonding',
     +   '             Ooi''s')
         write(outfi,2)
 2       format
     + (' strk chain/ l amino u & S structure   bridge        dihedral',
     +  ' angles     ',
     +  '                donor    acceptor   donor    acceptor  N  N')
         write(outfi,3)
 3       format
     +  (' num  seq.no t acids t S + patterns   partners    phi    psi',
     +  '  omega  alpha ',
     +  'kappa   tco to/energy fr/energy to/energy fr/energy  8 14')
      else if (hednum .eq. 2) then
         Write(outfi,4)
 4       format('                   S Atoms T')
         Write(outfi,5)
 5       format
     +  (' strk chain/ amino t prsnt h  side chain dihedral angles ',
     +  '  disulphide bridges                  C alpha distances for')
         Write(outfi,6)
 6       format
     +  (' num  seq.no acids d W M S m  Chi1   Chi2   Chi3   Chi4  ',
     +  ' ptr T CACA CBCB  S/S  i+1  i+2  i+3  i+4  i+5  i+6  i+7',
     +  '  i+8  i+9 i+10 i+11')
      else if (hednum .eq. 3) then
         write (outfi,7)
 7       format
     +  ('             A       S A    dihedral    ',
     +  '                                                     Ooi''s')
         write(outfi,8)
 8       format
     +  (' strk chain/ l amino t u     angles                   ',
     +  'C alpha distances for                   N  N')
         write(outfi,9)
 9       format
     +  (' num  seq.no t acids d t  alpha kappa  i+1  i+2  i+3  i+4  ',
     +  'i+5  i+6  i+7  i+8  i+9 i+10 i+11  8 14')
      end if
      end
c
c*******************************************************************************
c
      logical  function xeqy(x, y)
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
      logical  function xney(x, y)
c     implicit none
c     implicit complex (a-z)
      real     x, y
c===============================================================================
c     This function returns TRUE is its two arguments (X and Y) are not "equal",
c     else it returns "FALSE". How "equal" X and Y have to be is determined
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
      integer function lastsl (string)
c     implicit complex(a-z)
c     implicit none
      character *(*) string
      integer   i, lens
c
      lens = len(string)
      do 10 i = lens, 1, -1
      if (string(i:i).eq.'/')  then
         lastsl = i
         return
      end if
   10 continue
      lastsl = 0
      return
      end
c
c*****************************************************************************

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

