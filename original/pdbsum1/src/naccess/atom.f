c     *** subroutines handling atoms ***
c     from S. Hubbards' NACCESS
c     broken down for customisation purposes
c     $Id: atom.f,v 1.1 2000/05/19 17:56:35 hpo Exp $
c     $Log: atom.f,v $
c     Revision 1.1  2000/05/19 17:56:35  hpo
c     Initial revision
c
      integer function chain(nc, names, c)
c --  chain-id c encountered already?
c --  yes: return number in list of encountered chain-ids
c --  no: append to list & return number
      include 'accall.pars'
      integer     nc
      character*1 names(maxc), c
      do i = 1, nc
         if(c.eq.names(i))then
            chain = i 
            return
         endif
      enddo      
      nc = nc + 1
      names(nc) = c
      chain = nc
      return
      end
c     **********************************
      subroutine gatom(atom,anames,nres,nats,ok,ir,ia)
c --  find atom name ATOM among the standard residues
c --  return ok=.true. if found, ir=standard residue index, ia=atom index
      include 'accall.pars'
      integer      nres, nats(maxr), ir, ia
      character*4  atom, anames(maxr,maxa)
      logical      ok
      ok=.false.
      ir=0
      do while ( ir.lt.nres .and. .not.ok )
         ir=ir+1
         call ratom(atom,anames,ir,nats,ok,ia)
      enddo
      return
      end
c     **********************************
      subroutine ratom(atom,anames,ires,nats,ok,find)
c     test if atom ATOM of residue with index IRES is standard
c
c --  Checks to see if a "standard" atom name has been read in.
c --  Standard atom names are read from the file "vdw.radii", for
c --  the defined residue types therein. OK=.true. if found, and
c --  residue type = ires, and find = integer identifier of atom.
c
      include 'accall.pars'
      integer      ires, nats(maxr), find
      character*4  atom, anames(maxr,maxa)
      logical      ok
      ok=.false.
      find=0
      if(ires.eq.0.or.ires.gt.maxr)return
      do i = 1, nats(ires)
         if(atom.eq.anames(ires,i))then
            find=i
            ok=.true.
            return
         endif
      enddo
      return
      end
c     **********************************
      integer function restype(card)
c     return the CARD type as integer ATOM(1),HETATM(2),waterHOH(3)
      character*20 card
      restype=0
      if(card(1:4).eq.'ATOM')restype=1
      if(card(1:6).eq.'HETATM')restype=2
      if(card(18:20).eq.'HOH')restype=3      
      return
      end
c     **********************************
      subroutine which3(res,acids,ires,nacids,ok)
c
c -- Search array "acids" for existence of residue "res".
c -- OK = .true. if found, and index of res returned in "ires".
c
      include 'accall.pars'
      integer      ires, nacids
      character*3  acids(maxr), res
      logical      ok
      ires=nacids+1
      do i = 1, nacids
         if(res.eq.acids(i))then
            ires=i
            ok=.true.
            return
         endif
      enddo
      ok=.false.
      return
      end
c     **********************************
