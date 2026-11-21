c     *** subroutines concerning VDW radii ***
c     from S. Hubbards' NACCESS
c     broken down for customization
c     $Id: vdwaals.f,v 1.2 2005/04/01 11:40:46 hpo Exp $
c     $Log: vdwaals.f,v $
c     Revision 1.2  2005/04/01 11:40:46  hpo
c     Fixed a bug in vguess() where all guesses where overriden by a radius of 1.80.
c
c     Revision 1.1  2000/05/19 18:04:35  hpo
c     Initial revision
c
c     set default Van-der-Waals radii
c     ------------------------------
      subroutine vguess ( atom, vdw )
      real         vdw
      character*4  atom
      vdw=1.80
c
c -- Make a guess then !
c
      if(atom(2:2).eq.'C')vdw=1.80
      if(atom(2:2).eq.'N')vdw=1.60
      if(atom(2:2).eq.'S')vdw=1.85
      if(atom(2:2).eq.'O')vdw=1.40
      if(atom(2:2).eq.'P')vdw=1.90
      if(atom(1:2).eq.'CA')vdw=2.07
      if(atom(1:2).eq.'FE')vdw=1.47
      if(atom(1:2).eq.'CU')vdw=1.78
      if(atom(1:2).eq.'ZN')vdw=1.39
      if(atom(1:2).eq.'MG')then vdw=1.73
      return
      end

c     ****************************************
c     read VDW radii from external file
c     ---------------------------------
      subroutine vanin ( 
     -     vname,
     -     vlen,
     -     nacids, 
     -     aacids, 
     -     anames, 
     -     numats, 
     -     vradii,
     -     spolar,
     -     rtype
     -     )
c     
c -- Read in van der Waal radii from external file "vfile"
c -- nacids = number of residues read in
c -- aacids = *3 character array containing amino acid names
c -- anames = *4 character array containing atom names for each residue
c -- numats = array containing number of atoms for each residue
c -- radii  = vdw radii for each atom, indexed identically to anames
c
      include 'accall.pars'
      integer       readstring, parse, fopen, readint
      real          readfloat
      integer       nacids, 
     -              numats(maxr),
     -              l(256),
     -              nlabs,
     -              spolar(maxr,maxa),
     -              rtype(maxr),
     -              vlen
      real          vradii(maxr,maxa)
      character*3   aacids(maxr), aa3
      character*4   anames(maxr,maxa), aa4
      character*256 card, c(256), vname
      logical       ok
      if( fopen(10,vname,vlen,'old').eq.0 )then
         STOP'ERROR: unable to open "vdw.radii"'
      endif
      nacids = 0
      nlabs  = 0      
      do while ( readstring(10,card,ilen).ge.0 )
         n = parse(card,ilen,' ',c,l)
         if( c(1).eq.'RESIDUE' )then
            nacids = nacids + 1
            rtype(nacids)=1
            if(c(2)(1:4).eq.'NUCL')rtype(nacids)=2               
            if(c(2)(1:4).eq.'HETA')rtype(nacids)=3               
            if( nacids.gt.maxr )STOP'ERROR: increase maxr'
            aa3 = c(3)(1:3)
            do i = 1, 3
               if(aa3(i:i).eq.'_')aa3(i:i)=' '
            enddo
            aacids(nacids) = aa3
            numats(nacids) = 0
         endif
         if(c(1).eq.'ATOM')then
            numats(nacids) = numats(nacids) + 1
            if( numats(nacids).gt.maxa )STOP'ERROR: increase maxa'
c            aa4 = card(6:9)
c            do i = 1, 4
c               if(aa4(i:i).eq.'_')aa4(i:i)=' '
c            enddo
            anames(nacids,numats(nacids)) = card(6:9)
            vradii(nacids,numats(nacids)) = readfloat(card(11:14),4)
c            write(6,'(2a,1x,a,f8.3)')
c     -           'AT> ',anames(nacids,numats(nacids)),
c     -           aacids(nacids),vradii(nacids,numats(nacids))
            if(n.ge.4)then
               spolar(nacids,numats(nacids)) = readint(card(16:16),1)
            else
               spolar(nacids,numats(nacids)) = -1
            endif
         endif
      enddo
      close(10)
      return
      end
c     ****************************************
