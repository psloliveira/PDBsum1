      subroutine delta_solva(nats,xyz,rads,achain,accs,probe,zslice)
c** $Id: delta_solva.f,v 1.4 2002/01/17 14:53:14 hpo Exp $
c** $Log: delta_solva.f,v $
c** Revision 1.4  2002/01/17 14:53:14  hpo
c** save before making code g77 compliant
c**
c** Revision 1.3  2000/10/24 18:59:52  hpo
c** *** empty log message ***
c**
c** Revision 1.2  2000/10/24 18:57:37  hpo
c** minor bug fix
c**
c** Revision 1.1  2000/05/19 17:53:59  hpo
c** Initial revision
c**
c*************************************************************
c**  SOLVA - LEE & RICHARDS TYPE ACCESSIBLITY CALCULATIONS  **
c**  MODIFIED to give accessibility hidden between chains   **            
c*************************************************************
c**
c** Calculate accessible surface area for a group of atoms.
c** The accessible area for a given atom is calculated by the 
c** formula:
c**     (arcsum) x (atom radius+probe radius) x (deltaz)
c** Numerical integration is carried out over z. in each z-
c** section, the arcsum for a given atom is the arclength of 
c** the circle (intersection of the atom sphere with the z-
c** section) that is not interior to any other atom circles
c** in the same z-section.
c**
c*************************************************************
c**
c**  error parameter  - this gives accuracy of calculation
c**                   - suitable values are 0.01 (high 
c**                     accuracy) to 0.1 (low accuracy)
c**                   - in detail the z sections are spaced
c**                     at about 1/zslice*diameter of atom
c**
c**  probe size       - radius of probe in angstroms
c**                   - suitable value for water = 1.4
c**
c*************************************************************
      implicit none
      include 'accall.pars'
      real PRECIS
      parameter(PRECIS=1.0e-6)
c
c     the following are dimensioned to the max no of atoms (maxs)
c
      integer achain(maxs)
      real    xyz(maxs,3),
     -        rads(maxs),
     -        rad(maxs),
     -        radsq(maxs) 
c     
c     the following are dimensioned to the max no of intersections
c     of neighbouring spheres (nint)
c     
      dimension
     -        inov(nint),
     -        tag(nint),
     -        arc(nint),
     -        dx(nint),
     -        dy(nint),
     -        d(nint),
     -        dsq(nint)

      integer cubidx(3),cub2at(nac,ncube), natincub(ncube), 
     -     at2cub(maxs)
      integer nats, tag , inov, karc,
     -        i,j,k,l,n,m,nm
      integer ir, in, io, nzp
      real    accs(maxs),
     -        arc,dx,dy,d,dsq, b, alpha, beta,
     -        probe, zslice, trig_test, ti, tf
      real    pi, pix2, coomin(3), coomax(3)
      real    zres,zgrid,area,rsecn,rsec2n,rsecr,rsec2r
      real    delta_arc, darc, delarc
c     flag intra-chain intersections
      logical dfl(nint), chfl,find_near_atoms
c
c     initialise variables, constants
      call set_solva_par(nats,xyz,rads,probe,coomin,coomax,
     -     cubidx, cub2at, at2cub, natincub,radsq,rad)

      pi=acos(-1.0)
      pix2=2.0*pi
      karc=nint
c
c     -- Process each atom in turn
c     number of z-slices per angstroem 
c     (=! number of z-slices/atom)
      nzp=1./zslice+0.5
      do  ir = 1, nats
         darc=0.
c     find neighbouring atoms, 
c     return true if multiple chain id
         chfl = find_near_atoms(ir,xyz,achain,rad,
     -        cubidx, cub2at, at2cub, natincub,
     -        io, inov,dx,dy,dsq,d)
c        print*,"io: ",io
c        print*,(inov(l),l=1,io)
c        print*
         
         if ((.not.chfl).or.(io.lt.1)) goto 18
c     
c     z resolution determined
c     the following procedure keeps the number of Z-slices through 
c     the atom constant, is this an attempt to ensure constant 
c     relative error? doesn't work - but keep it that way
            zres=2.*rad(ir)/nzp
            zgrid=xyz(ir,3)-rad(ir)-zres/2. 

c     
c     section atom spheres perpendicular to the z axis
         do i = 1, nzp
            zgrid=zgrid+zres
c            print*,'**** new zslice:',zgrid
c
c     find the radius of the circle of intersection of 
c     the ir sphere on the current z-plane
            rsec2r=radsq(ir)-(zgrid-xyz(ir,3))**2           
            rsecr=sqrt(rsec2r)
c            write(*,'(a,f7.4,a,f7.4,a,f7.4)')
c     +           'rad:',rad(ir),'; z:',xyz(ir,3),'rsecr:',rsecr

c     -- reset the arcs and arc flags
            do k = 1, karc
               arc(k)=0.0
c              reset flag for intr-chain intersections
               dfl(k)=.false.
            enddo
            karc=0

c     -- loop over all neighbouring atoms
            do j = 1, io
               in=inov(j)
    
c     find radius of the circle of intersection of 
c     of the in sphere on the current z-plane
               rsec2n=radsq(in)-(zgrid-xyz(in,3))**2
               if (rsec2n.le.0.0) goto 10
               rsecn=sqrt(rsec2n)
c
c     find intersections of n.circles with ir circles in section
c     too far apart -> goto 10, next atom
               if (d(j).gt.rsecr+rsecn) goto 10
               chfl=achain(ir).ne.achain(in)
c
c     do the circles intersect, or is one circle completely inside the other?
               b=rsecr-rsecn
               if (d(j).gt.abs(b)) goto 20
c     here, the circles are inside each other
c     ir encloses in -> nothing to add, next atom
               if (b.gt.0.0) goto 10
c              ir competeley covered by in 
c              if in same chain -> nothing to add
               if(.not.chfl) goto 9
               karc=karc+1
               if(karc.ge.nint)
     -              STOP'SOLVA_ERROR: max intersections exceeded2'
c     begin of a full arc
               arc(karc)=0
               dfl(karc)=chfl
               karc=karc+1
               dfl(karc)=chfl
               arc(karc)=pix2
               goto 10
c     arc(1),arc(2) -> begin, end of arc 1
c     dfl(1) flag for arc(1), arc(2) indicating intersection with different chain
c     dfl(2) flag for arc(3), arc(4)
c     int(karc/2+.5) == 1 => karc points to end of arc
c     int(karc/2+.5) == 0 => karc points to start of arc
c
c     if the circles intersect, find the points of intersection
c
 20            karc=karc+1
               if(karc.ge.nint-1)then
                  STOP'SOLVA_ERROR: max intersections exceeded2'
               endif
c
c     Initial and final arc endpoints are found for the ir circle intersected
c     by a neighboring circle contained in the same plane.     
               trig_test=(dsq(j)+rsec2r-rsec2n)/(2.*d(j)*rsecr)
               if((trig_test-1.0).gt.PRECIS)trig_test=1-PRECIS
               if((-1.0-trig_test).gt.PRECIS)trig_test=-1+PRECIS
               alpha=acos(trig_test)
c    
c     alpha is the angle between a line containing a point of intersection and
c     the reference circle center and the line containing both circle centers
c     tan(y/x)=atna2(y,x)
c     origin of circle on neg. x-axis, anti-clockwise
               beta=atan2(dy(j),dx(j))+pi
c               print*,'ir:',ir,'; in:',in,
c     +              '; beta=',beta,'; alpha=',alpha,
c     +              '; dx=',dx(j),'; dy=',dy(j)
c     
c     beta is the angle between the line containing both circle centers and the x-axis 
c     cut the hidden segment! (ti<->tf) otherwise)
               ti=beta-alpha
               tf=beta+alpha
               if(tf.gt.pix2)tf=tf-pix2
               if(ti.lt.0.0) ti=ti+pix2
               arc(karc)=ti
               dfl(karc)=chfl
               if(tf.lt.ti) then 
c     
c     if the arc crosses zero, then it is broken into two segments.
c     the first ends at pix2 and the second begins at zero     
                  karc=karc+1
                  if(karc.ge.nint-2)then
                     STOP'SOLVA_ERROR: max intersections exceeded2'
                  endif
                  arc(karc)=pix2
                  dfl(karc)=chfl
c     new beginning at 0.0
                  karc=karc+1
                  arc(karc)=0.0
                  dfl(karc)=chfl
               endif
               karc=karc+1
               arc(karc)=tf
               dfl(karc)=chfl
c     next atom
 10         enddo
cxx   testing
cxx            do k=1, karc
cxx               write(*,'(f6.3,1x,i1)')  arc(karc), dfl(karc)
cxx            enddo
cxx   end testing

            if (karc.gt.0) then
               delarc = delta_arc(karc,arc,dfl)      
               darc = darc + delarc
c
c     find the in-accessible surface area for the sphere ir on this section
c               write(*,'(i2,a,i2,a,f6.3,a,f7.4,a,f7.4,a,f7.4)') 
c     -              i,') karc', karc,
c     -              '; zgrid: ', zgrid, '; delarc: ', delarc,
c     -              '; absarc: ', (delarc*rsecr), 'rad(ir)',rad(ir)
            endif
c            print*,'+++ arcsum:', delarc, 
c     +           'totsum:',darc
c     next zslice
 9       enddo
c     
c     scale area to vdw shell
c     
 18      accs(ir)=darc*zres*rad(ir)
c        print*, 'atom nr. ',ir,': dasa =', accs(ir)
c------------------------------------------------------------------
c The following line converts from accessible to contact surface
c         c=(b*(rad(ir)-probe)**2)/(rad(ir)**2)
c------------------------------------------------------------------
c     next ir (next central atom)
      enddo

      write(4,'(a)')' SOLVA: PROGRAM ENDS CORRECTLY'
      return
      end
c     ***********************************************************
      subroutine set_solva_par(nats,xyz,rads,probe,coomin,coomax,
     -     cubidx, cub2at, at2cub, natincub,radsq,rad)
c     set xyz coordinate minima & maxima 
c          -> coomin(3), coomax(3)
c     set the jumps used to convert from 3D to single cube index 
c          -> cubidx(3)
c     set pointers from cube i to atoms j in cube
c           -> natincub(i) number of atoms in cube i
c           -> cub2at(j,i) atom index of atom j in cube
c     set pointer from atom j to cube
c           -> at2cub(j)
      implicit none
      include 'accall.pars'
c     input:
      integer nats
      real probe, xyz(maxs,3), rads(maxs)
c     return:
      integer cubidx(3), cub2at(nac,ncube), at2cub(maxs),
     -     natincub(nac)
      real radsq(maxs), rad(maxs), coomin(3), coomax(3)

c     other variables
      integer i, l, n, kji, j, k
      real rmax

c
c     -- Radius of an atom sphere = atom radius + probe radius
c     -- Find maxima and minima
c
      rmax=0.0

c     set initial values
      do l=1, 3
         coomin(l)=xyz(1,l)
         coomax(l)=xyz(1,l)
      enddo
      do i = 1, nats
         rad(i)   = rads(i) + probe
         radsq(i) = rad(i)**2
         if (rad(i).gt.rmax)rmax = rad(i)
         do l=1, 3
            if(coomin(l).gt.xyz(i,l)) coomin(l)=xyz(i,l)
            if(coomax(l).lt.xyz(i,l)) coomax(l)=xyz(i,l)
         enddo
      enddo
c
c     rmax = max diameter
c
      rmax = rmax*2.
c
c     -- 1D array of Cubicals containing the atoms are setup. 
c     -- length of an edge = largest atom sphere diamete
c     -- Minimum of 3 x 3 cubic grid, EXIT if max cubes exceeded
c
      do l=1, 3
         cubidx(l)=1+(coomax(l)-coomin(l))/rmax
         if(cubidx(l).lt.3) cubidx(l)=3
         if (l.gt.1) cubidx(l)=cubidx(l)*cubidx(l-1)
      enddo
      if(cubidx(3).ge.ncube)then
         write(*,*) 'cubidx = ',cubidx(3)
         STOP'SOLVA_ERROR: max cubes exceeded'
      endif
c
c     -- Prepare upto ncube cubes each containing upto nac atoms. The cube index
c     -- is kji. The atom index for each cube is cub2at(ncube)
c
      do l = 1, cubidx(3)
        natincub(l)=0
      enddo

      do l=1, nats
      
c     calculate index of cube containing atom l
         i = 1+(xyz(l,1)-coomin(1))/rmax
         j =   (xyz(l,2)-coomin(2))/rmax
         k =   (xyz(l,3)-coomin(3))/rmax
         kji = k*cubidx(2) + j*cubidx(1) + i
         if (kji.lt.0) STOP'SOLVA_ERROR: wrong limits'

c     increment natincub(kji)
         n = natincub(kji)+1
         if(n.gt.nac)STOP'SOLVA_ERROR: max atoms per cube exceeded'
         natincub(kji) = n
c     setting pointers
         cub2at(n,kji) = l
         at2cub(l) = kji
      enddo
      return
      end
c **************************************************************
      logical function find_near_atoms(ir,xyz,chid,rad,
c          variables concerning cubicles
     -     cubidx, cub2at, at2cub, natincub,
c          return variables
     -     io, inov,dx,dy,dsq,d)
c     find neighbouring atoms to central atom index ir
c     return true if multiple chain ids encountered
c     return number of neighbours io
c     atom indices in array inov(nac)
c     also return distances to central atom
      implicit none
      include 'accall.pars'
c     input variables
      integer ir, cubidx(3), cub2at(nac,ncube), at2cub(maxs),
     -     natincub(nac), chid(maxs)
      real xyz(maxs,3), rad(maxs)
c     return variables
      integer io, inov(nint)
      real dx(nint), dy(nint), dsq(nint), d(nint)
c     internal variables
      integer i,j,k,m,mkji,nm,in
      real dzsq, drsq

      find_near_atoms=.false.
      io=0
      

c
c     -- Find the 'mkji' cubes neighboring the kji cube
c
      do k = -1, 1, 1
         do j = -1, 1, 1
            do i = -1, 1, 1
               mkji=k*cubidx(2)+j*cubidx(1)+i+at2cub(ir)
               if(mkji.ge.1) then
                  if(mkji.gt.cubidx(3)) goto 14
                  nm=natincub(mkji)
                  if(nm.ge.1)then
c     
c     -- record the atoms in inov that neighbor atom ir
c     -- flag different chain id
c
                     do m = 1, nm
                        in=cub2at(m,mkji)
                        if (in.ne.ir)then
                           io=io+1  
                           if (io.gt.nint)then
                              STOP'SOLVA_ERROR: intrsctns > max'
                           endif
                           dx(io)=xyz(in,1)-xyz(ir,1)
                           dy(io)=xyz(in,2)-xyz(ir,2)
                           dzsq = (xyz(in,3)-xyz(ir,3))**2
                           dsq(io)=dx(io)**2+dy(io)**2
                           drsq =(rad(ir)+rad(in))**2 
                           d(io)=sqrt(dsq(io))
                           if (drsq.lt.dsq(io)+dzsq) then
                              io=io-1
                           else
                              inov(io)=in
c     -- flag different chain
                              if(chid(ir).ne.chid(in)) 
     -                             find_near_atoms=.true.
                           endif
                        endif
                     enddo
c
                  endif
               endif
            enddo
         enddo
      enddo

 14   return
      end
c **************************************************************
      real function delta_arc(karc,arc,dfl)
c     calculate sum of arcs
c     cut by inter- but not by intra-chain neighbouring atoms
c     arc(k) has two types:
c       inter (dfl(k)=.true.) and intra-chain (dfl(k)=.false.)
c     k even: beginnig of an arc, k odd: end of an arc
      implicit none
      include 'accall.pars'

      integer k,ev,karc, tag(karc), hetctr, homctr, hfk
      real arc(nac), alpha, od
      logical dfl(nac), ifl

      delta_arc=0.0
      if (karc.eq.0) goto 21
      if ((karc-2*int(karc/2.)).gt.0) then 
         STOP'SOLVA:DELTA_ARC ERROR: ends of intersections do not match'
      endif
       
c     sort by size
      call sortag(arc(1),karc,tag)

      alpha =0.0
      hetctr=0
      homctr=0

      do k=1, karc
c         write(*,'(i2,a,f7.4,a,i2,a,i2)') k,') ',arc(k),'::',tag(k),
c     +   '; type:', dfl(tag(k)) 
c         print*,dfl(tag(k))
         if((tag(k)-2*int(tag(k)/2.)).gt.0) then
c     beginning of arc (odd k)
            if (dfl(tag(k))) then
c     ... belongs to hetero arc, increment hetero counter
               hetctr=hetctr+1
               if((hetctr.eq.1).and.(homctr.eq.0)) alpha = arc(k)
            else
c     ... belongs to homo arc, increment homo counter
               homctr = homctr+1
               if((hetctr.gt.0).and.(homctr.eq.1)) then
                  delta_arc = delta_arc + arc(k) - alpha
c                  write (*,'(a,f7.4)') 'delta_arc=',arc(k) - alpha
               endif
            endif            
         else
c     end of an arc (even k)
            if (dfl(tag(k))) then
c     ... belongs to a hetero arc, decrement hetero counter
               hetctr=hetctr-1
               if((hetctr.eq.0).and.(homctr.eq.0)) then
                  delta_arc = delta_arc + arc(k) - alpha
c                  write (*,'(a,f7.4)') 'delta_arc=',arc(k) - alpha
               endif
            else
c     ... belongs to homo arc, decrement homo counter
               homctr = homctr-1
               if((hetctr.gt.0).and.(homctr.eq.0)) alpha = arc(k)
            endif         
         endif
      enddo
 21   return
      end
            
         
      
      
