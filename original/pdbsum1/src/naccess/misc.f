c     *** miscellaneous routines  
c     of S. Hubbards' NACCESS
c     broken down for customization purposes
c     $Id: misc.f,v 1.1 2000/05/19 18:02:57 hpo Exp $
c     $Log: misc.f,v $
c     Revision 1.1  2000/05/19 18:02:57  hpo
c     Initial revision
c
      integer function fopen( iochan, filename, flen, fstat )
c     open file
      integer       iochan, flen
      character*(*) filename, fstat
      open (
     -      unit   = iochan,
     -      file   = filename(1:flen),
     -      status = fstat,
     -      err    = 100
     -     )
      fopen = 1
      return
 100  fopen = 0
      return
      end
c     ******************************
      integer function readstring(file,card,flen)
c     read a string from unit FILE
      integer       file, flen
      character*256 card

      if(file.gt.200)STOP'ERROR: file number too large'
      read(file,'(a)',err=100,end=100)card
      flen=256
      do while(card(flen:flen).eq.' ') 
         flen=flen-1
      enddo
      readstring=flen
      return
 100  readstring=-1
      return
      end
c     ******************************
      integer function parse(card,length,separator,chars,clen)
c     partition a string into fields separated by the character SEPARATOR
c     return the number of fields as dunction value
c     return fields in CHARS[i] and field length in CLEN[i]
      integer       charpos(256,2), clen(256), length
      character*1   separator
      character*(*) card, chars(256)
      logical       search

      parse=0
      i=0
      search=.false.
      do while(i.lt.length)
         i = i + 1
         if(.not.search)then
            if(card(i:i).ne.separator)then
               parse = parse + 1
               charpos(parse,1) = i
               search=.true.
            endif
         else
            if(card(i:i).eq.separator)then
               charpos(parse,2) = i-1
               search=.false.
            endif
         endif
      enddo
      if(search)charpos(parse,2)=length
      do i = 1, parse
         chars(i)=card(charpos(i,1):charpos(i,2))
         clen(i)=charpos(i,2)-charpos(i,1)+1
      enddo
      return
      end
c     ******************************
      integer function readint(card,clen)
      integer       clen  
      real          value
      character*(*) card
      read(card(1:clen),*,err=100)value
      readint=int(value)
      return
 100  readint=-999
      return
      end
c     ******************************
      real function readfloat(card,clen)
      integer       clen
      real          value
      character*(*) card

      read(card(1:clen),*,err=100)value
      readfloat=value
      return
 100  readfloat=-999.9
      return
      end
c     ******************************
      subroutine tolow(text,clen)
      integer       clen
      character*(*) text
      character*1   alph(26,2)
      logical       ok

      data alph/'A','B','C','D','E','F','G','H','I','J',
     -          'K','L','M','N','O','P','Q','R','S','T',
     -          'U','V','W','X','Y','Z',
     -          'a','b','c','d','e','f','g','h','i','j',
     -          'k','l','m','n','o','p','q','r','s','t',
     -          'u','v','w','x','y','z'/

      do i = 1, clen
         j=0
         ok=.true.
         do while( j.lt.26.and.ok ) 
            j = j + 1
            if ( alph(j,1).eq.text(i:i) ) then
               ok=.false.
               text(i:i)=alph(j,2)
            endif
         enddo
      enddo

      return
      end
c     ******************************
      subroutine polguess(atom,ip)
c     is atom ATOM a backbone atom (N,O,CA,HA)? ip=1 TRUE, ip=0 false
      character*4 atom
      ip=0
      if(atom(2:2).eq.'O'.or.atom(2:2).eq.'N'.or.atom(2:2).eq.'A')ip=1
      return
      end
c     ******************************
      integer function what_atom(atom,ir,n)
c     test membershipf of atom in set depending on ir
c     return what_atmo=0 if part of set, =1 if not part of set
c     ir=1: set is main chains
c      ir=2: sugar atoms ?
      character*4 atom, mc(5), nc(11)
      data mc/' N  ',' C  ',' O  ',' OXT',' CA '/        
      data nc/' P  ',' O1P',' O2P',' O5*',' C5*',' C4*',
     -        ' O4*',' C3*',' O3*',' C2*',' C1*'/
      what_atom=0
      if(ir.eq.1)then
         do i = 1, n
            if(atom.eq.mc(i))return
         enddo
      elseif(ir.eq.2)then
         do i = 1, 11
            if(atom.eq.nc(i))return
         enddo
      else
         do i = 1, 4
            if(atom.eq.mc(i))return
         enddo
         do i = 1, 11
            if(atom.eq.nc(i))return
         enddo
      endif
      what_atom=1
      return
      end
 
