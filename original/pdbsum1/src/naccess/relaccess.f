c     from S. Hubbards NACCESS program
c     broken sown for customization
c     $Id: relaccess.f,v 1.1 2000/05/19 18:14:57 hpo Exp $
c     $Log: relaccess.f,v $
c     Revision 1.1  2000/05/19 18:14:57  hpo
c     Initial revision
c
      subroutine summer(
     -     sname,
     -     slen,
     -     nats,
     -     accs,
     -     backbone,
     -     polstats,
     -     rtype,
     -     resindex,
     -     resnam,
     -     ressums,
     -     tsums
     -     )

c     -- 	program to sum atomic accessibilities by residue.
c     --	copes with atom and hetatom records
c     --	produces relative accessibilities for the 20 common aminos
c     --	ouput written to .rsa file (channel 4)      

      include 'accall.pars'
      integer       fopen, readstring
      integer       nats, slen,
     -              backbone(maxs), 
     -              polstats(maxs), 
     -              resindex(maxs),
     -              rindex(maxx),
     -              rtype(maxr),
     -              nacids
      real	    ressums(maxx,7,2),
     -              tsums(7), 
     -              accs(maxs),
     -		    standarea(maxr,7)
      character*256 line, sname
      character*10  resnam(maxx)
      character*3   res, acids(maxr)
      logical       stand, ok
c
c     if "standard.data" exists in, read them in.
c
      stand=.false.
      if( fopen(1,sname,slen,'old').ne.0 )then
         write(3,'(4a)')'REM  Relative accessibilites read from',
     -                  ' external file "',sname(1:slen),'"'
         stand=.true.
         i=0
         do while( readstring(1,line,ilen).ge.0.and.i.lt.maxr )
            if(line(1:4).eq.'ATOM')then
               i=i+1
               acids(i)=line(13:15)
               read(line(17:23),'(f7.2)')standarea(i,1)
               read(line(30:36),'(f7.2)')standarea(i,2)
               read(line(43:49),'(f7.2)')standarea(i,3)
               read(line(56:62),'(f7.2)')standarea(i,4)
               read(line(69:75),'(f7.2)')standarea(i,5)
               read(line(82:88),'(f7.2)')standarea(i,6)
               read(line(95:108),'(f7.2)')standarea(i,7)
            endif
         enddo
         write(4,'(a,i3,a)')
     -        ' RELATIVE (STANDARD) ACCESSIBILITIES READFOR ',
     -        i,' AMINO ACIDS'
         close(1)
      else
         write(4,'(a)')' NO STANDARD VALUES INPUT'
      endif

      nacids=i
      do i = 1, resindex(nats)
         rindex(i)=0
         if(stand)then
            res=resnam(i)(1:3)
            call which3(res,acids,ires,nacids,ok)
            if(ok)rindex(i)=ires
         endif
      enddo
c
c     -- sum the values
c     
      do i = 1, nats
         ir=resindex(i)
         tsums(1)=tsums(1)+accs(i)
         ressums(ir,1,1)=ressums(ir,1,1)+accs(i)
         if(backbone(i).eq.0)then
            ressums(ir,5,1)=ressums(ir,5,1)+accs(i)
            tsums(5)=tsums(5)+accs(i)
         else
            ressums(ir,4,1)=ressums(ir,4,1)+accs(i)
            tsums(4)=tsums(4)+accs(i)
            if(polstats(i).eq.0)then
               ressums(ir,2,1)=ressums(ir,2,1)+accs(i)
               tsums(2)=tsums(2)+accs(i)
            elseif(polstats(i).eq.1)then
               ressums(ir,3,1)=ressums(ir,3,1)+accs(i)
               tsums(3)=tsums(3)+accs(i)
            endif
         endif
         if(polstats(i).eq.0)then
            ressums(ir,6,1)=ressums(ir,6,1)+accs(i)
            tsums(6)=tsums(6)+accs(i)
         else
            ressums(ir,7,1)=ressums(ir,7,1)+accs(i)
            tsums(7)=tsums(7)+accs(i)
         endif
      enddo
c
c     -- calculate relative accessibilities
c
      do i = 1, resindex(nats)
         ires=rindex(i)
         if(stand.and.ires.ne.0)then
            do j = 1, 7
               if(standarea(ires,j).gt.0.0)then
                  ressums(i,j,2)=100.0*ressums(i,j,1)/standarea(ires,j)
               else
                  ressums(i,j,2)=0.0
               endif
            enddo
         else
            do j = 1, 7
               ressums(i,j,2)=-99.9
            enddo               
         endif
      enddo
      return
      end
