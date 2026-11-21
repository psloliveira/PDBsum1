c     *** sorting routine from S. Hubbards' NACCESS
c     broken down for customization purposes
c     $Id: sorting.f,v 1.1 2000/05/19 17:58:19 hpo Exp $
c     $Log: sorting.f,v $
c     Revision 1.1  2000/05/19 17:58:19  hpo
c     Initial revision
c
      subroutine sortag(a,n,tag)
c     sort array a(n), return the permutation in array tag(n)
      integer tag,tg
      dimension a(n),iu(16),il(16),tag(n)
      do i = 1, n
         tag(i)=i
      enddo
      m=1
      i=1
      j=n
 5    if(i.ge.j) go to 70
 10   k=i
      ij=(j+i)/2
      t=a(ij)
      if(a(i).le.t) go to 20
      a(ij)= a(i)
      a(i)=t
      t=a(ij)
      tg=tag(ij)
      tag(ij)=tag(i)
      tag(i)=tg
 20   l=j
      if(a(j).ge.t) go to 40
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      tg=tag(ij)
      tag(ij)=tag(j)
      tag(j)=tg
      if(a(i).le.t) go to 40
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      tg=tag(ij)
      tag(ij)=tag(i)
      tag(i)=tg
      go to 40
 30   a(l)=a(k)
      a(k)=tt
      tg=tag(l)
      tag(l)=tag(k)
      tag(k)=tg
 40   l=l-1
      if(a(l).gt.t) go to 40
      tt=a(l)
 50   k=k+1
      if(a(k).lt.t) go to 50
      if(k.le.l) go to 30
      if(l-i.le.j-k) go to 60
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 80
 60   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 80
 70   m=m-1
      if(m.eq.0) return
      i=il(m)
      j=iu(m)
 80   if(j-i.ge.1) go to 10
      if(i.eq.1) go to 5
      i=i-1
 90   i=i+1
      if(i.eq.j) go to 70
      t=a(i+1)
      if(a(i).le.t) go to 90
      tg=tag(i+1)
      k=i
 100  a(k+1)=a(k)
      tag(k+1)=tag(k)
      k=k-1
      if(t.lt.a(k)) go to 100
      a(k+1)=t
      tag(k+1)=tg
      go to 90
      end
