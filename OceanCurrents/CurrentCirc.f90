! =====================================================================================
! ASCIDIAN
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License,or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful,but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not,write to the Free Software Foundation,Inc.,59 Temple
! Place,Suite 330,Boston,MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  CurrentCirc.f90
!
!    Description:  Computes winds induced circulation.
!
!        Version:  1.0
!        Created:  21/09/15 14:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module currents_circ

  use currents_io
  use currents_data
  use parallel_circ

  implicit none

  real::g1,g2,g3,g4,g5,g6,g7

  public

contains

  ! =====================================================================================
  subroutine currents_circulation(fst)

    integer::kt,fst
    real::akate,afilt,bfilt,tpp

    call build_computationalgrid(fst)
    if(iam/=0) return

    if(.not.circon) return
    tpp=20.

    ! First compute the ocean circulation
    do kt=1,kend
       akate=(kt/real(kend))*100.
       alin=1.
       alin=real(kt)/niner
       if(kt>=niner) alin=1.
       if(akate>=tpp)then
          if(iam==0) write(*,13) int(akate)
          tpp=tpp+20.
       endif
       if(kt==kend-1) call build_output
       call update_bounds(kt)
       afilt=real(kt)/real(ifilt)
       bfilt=real(kt/ifilt)
       if(afilt==bfilt) call filter
       call sweepx_fct
       call update_bounds(kt)
       call sweepy_fct
    enddo

13  format('+',' currents circulation completion: ',i5,' %')

    return

  end subroutine currents_circulation
  ! =====================================================================================
  subroutine update_bounds(kt)

    integer::kt,km,kk,k1,k2,kn,k,halo

    real::amu1,zcrit,amu,amumu

    if(kt>=niner) goto 7
    halo=2

    do km=1,ne
       if(ie(km)<=halo+1)then
          if(kt==1) a(ie(km),je(km))=(float(kt)/niner)*a(ie(km),je(km))
          aa(ie(km),je(km))=(float(kt)/niner)*aa(ie(km),je(km))
       elseif(ie(km)>=mp-halo+1)then
          if(kt==1) a(ie(km),je(km))=(float(kt)/niner)*a(ie(km),je(km))
          aa(ie(km),je(km))=(float(kt)/niner)*aa(ie(km),je(km))
       endif
    enddo

    ! Solves the partially clamped explicit radiation condition
    if(kt==1) goto 7

    amu1=dt/atf
    zcrit=4427.0
    do km=1,ne
       ! Bottom margin
       if(ije(km)/=-1) goto 13
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)<mp.and.je(km)>1.and.je(km)<np) &
            amu=dt*sqrt(grav*0.5*(a(ie(km)+1,je(km)-1)+a(ie(km)+1,je(km)+1)))/dx
       if(ie(km)==ip(ie(km)).and.je(km)==ip(je(km)).and.je(km)>1.and.je(km)<np.and.ie(km)>0) &
            amu=dt*sqrt(grav*0.5*(a(ie(km),je(km)-1)+a(ie(km),je(km)+1)))/dx
       if(ie(km)>0.and.je(km)>0.and.ie(km)+4<mp)then
          aa(ie(km),je(km))=(a(ie(km),je(km))*(1.0-0.25*amu1)+ &
               amu*0.125*(4.0*a(ie(km)+2,je(km))-a(ie(km)+4,je(km))- &
               3.0*a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*dx/dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) aa(ie(km),je(km))=0.0

       ! Upper margin
13     if(ije(km)/=1) goto 14
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>1.and.je(km)>1.and.je(km)<np) &
            amu=dt*sqrt(grav*0.5*(a(ie(km)-1,je(km)-1)+a(ie(km)-1,je(km)+1)))/dx
       if(ie(km)==ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>0.and.je(km)>1.and.je(km)<np) &
            amu=dt*sqrt(grav*0.5*(a(ie(km),je(km)-1)+a(ie(km),je(km)+1)))/dx
       if(ie(km)<np.and.ie(km)>0.and.je(km)>0.and.ie(km)-4>0)then
          aa(ie(km),je(km))=(a(ie(km),je(km))*(1.0-0.25*amu1)+ &
               amu*0.125*(4.0*a(ie(km)-2,je(km))-a(ie(km)-4,je(km))- &
               3.0*a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*dx/dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) aa(ie(km),je(km))=0.0

       ! Left margin
14     if(ije(km)/=-2) goto 15
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>1.and.ie(km)<mp.and.je(km)<np) &
            amu=dt*sqrt(grav*0.5*(a(ie(km)+1,je(km)+1)+a(ie(km)-1,je(km)+1)))/dx
       if(ie(km)/=ip(ie(km)).and.je(km)/=ip(je(km)).and.ie(km)>1.and.ie(km)<mp.and.je(km)>0) &
            amu=dt*sqrt(grav*0.5*(a(ie(km)+1,je(km))+a(ie(km)-1,je(km))))/dx
       if(ie(km)>0.and.je(km)>0.and.je(km)+4<np)then
          aa(ie(km),je(km))=(a(ie(km),je(km))*(1.0-0.25*amu1)+ &
               amu*0.125*(4.0*a(ie(km),je(km)+2)-a(ie(km),je(km)+4)- &
               3.0*a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*dx/dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) aa(ie(km),je(km))=0.0

       ! Right margin
15     if(ije(km)/=2) goto 11
       if(ie(km)/=ip(ie(km)).and.je(km)==ip(je(km)).and.ie(km)>1.and.ie(km)<mp.and.je(km)>1) &
            amu=dt*sqrt(grav*0.5*(a(ie(km)+1,je(km)-1)+a(ie(km)-1,je(km)-1)))/dx
       if(ie(km)/=ip(ie(km)).and.je(km)/=ip(je(km)).and.ie(km)>1.and.ie(km)<mp.and.je(km)>0)&
            amu=dt*sqrt(grav*0.5*(a(ie(km)+1,je(km))+a(ie(km)-1,je(km))))/dx
       if(ie(km)>0.and.je(km)>0.and.je(km)-4>0)then
          aa(ie(km),je(km))=(a(ie(km),je(km))*(1.0-0.25*amu1)+&
               amu*0.125*(4.0*a(ie(km),je(km)-2)-a(ie(km),je(km)-4)- &
               3.0*a(ie(km),je(km))))/(1.0+0.25*amu1)
       endif
       amumu=amu*dx/dt
       if(amumu>zcrit.and.ie(km)>0.and.je(km)>0) aa(ie(km),je(km))=0.0

11     continue
    enddo

7   do kk=1,nnn
       if(ina(kk)>0.and.ina(kk)+2<=mp.and.jnab(kk)>0) &
            a(ina(kk),jnab(kk))=a(ina(kk)+2,jnab(kk))
       if(inb(kk)>2.and.inb(kk)<mp.and.jnab(kk)>0) &
            a(inb(kk),jnab(kk))=a(inb(kk)-2,jnab(kk))
    enddo

    do kk=1,mmm
       if(imab(kk)>0.and.jma(kk)>0.and.jma(kk)+2<=np)&
            a(imab(kk),jma(kk))=a(imab(kk),jma(kk)+2)
       if(imab(kk)>0.and.jmb(kk)>2)&
            a(imab(kk),jmb(kk))=a(imab(kk),jmb(kk)-2)
    enddo

    do k1=1,mmm
       do k2=1,nnn
          if(ina(k2)==imab(k1).and.jnab(k2)==jma(k1).and.ina(k2)>0.and.jnab(k2)>0.and.ina(k2)+2<=mp.and.jnab(k2)+2<=np) &
               a(ina(k2),jnab(k2))=(a(ina(k2)+2,jnab(k2))+a(ina(k2),jnab(k2)+2))/2.
          if(inb(k2)==imab(k1).and.jnab(k2)==jma(k1).and.inb(k2)>2.and.jnab(k2)>0.and.jnab(k2)+2<=np) &
               a(inb(k2),jnab(k2))=(a(inb(k2)-2,jnab(k2))+a(inb(k2),jnab(k2)+2))/ 2.
          if(ina(k2)==imab(k1).and.jnab(k2)==jmb(k1).and.ina(k2)>0.and.ina(k2)+2<=mp.and.jnab(k2)>2) &
               a(ina(k2),jnab(k2))=(a(ina(k2),jnab(k2)-2)+a(ina(k2)+2,jnab(k2)))/2.
          if(inb(k2)==imab(k1).and.jnab(k2)==jmb(k1).and.inb(k2)>2.and.jnab(k2)>2) &
               a(inb(k2),jnab(k2))=(a(inb(k2)-2,jnab(k2))+a(inb(k2),jnab(k2)-2))/2.
       enddo
    enddo

    do kn=1,nn
       do k1=1,ne
          if(in1(kn)==ie(k1).and.jn(kn)==je(k1)) goto 60
       enddo
       if(in1(kn)>0.and.jn(kn)>0.and.in1(kn)<=mp.and.jn(kn)<=np)then
          a(in1(kn),jn(kn))=0.0
          aa(in1(kn),jn(kn))=0.0
       endif
60     do k1=1,ne
          if(in2(kn)==ie(k1).and.jn(kn)==je(k1)) goto 30
       enddo
       if(in2(kn)>0.and.jn(kn)>0.and.in2(kn)<=mp.and.jn(kn)<=np)then
          a(in2(kn),jn(kn))=0.0
          aa(in2(kn),jn(kn))=0.0
       endif
30     continue
    enddo

    do km=1,mm
       do k1=1,ne
          if(im(km)==ie(k1).and.jm1(km)==je(k1)) goto 90
       enddo
       if(im(km)>0.and.jm1(km)>0.and.im(km)<=mp.and.jm1(km)<=np)then
          a(im(km),jm1(km))=0.0
          aa(im(km),jm1(km))=0.0
       endif
90     do k1=1,ne
          if(im(km)==ie(k1).and.jm2(km)==je(k1)) goto 40
       enddo
       if(im(km)>0.and.jm2(km)>0.and.im(km)<=mp.and.jm2(km)<=np)then
          a(im(km),jm2(km))=0.0
          aa(im(km),jm2(km))=0.0
       endif
40     continue
    enddo

    do k=1,nnnn
       if(inn1(k)/=ip(inn1(k)).and.inn1(k)>0.and.jnn(k)>0.and.inn1(k)+2<=mp.and.jnn(k)<=np) &
            a(inn1(k),jnn(k))=a(inn1(k)+2,jnn(k))
       if(inn2(k)/=ip(inn2(k)).and.inn2(k)>2.and.jnn(k)>0.and.inn2(k)<=mp.and.jnn(k)<=np) &
            a(inn2(k),jnn(k))=a(inn2(k)-2,jnn(k))
    enddo

    do k=1,mmmm
       if(jmm1(k)==ip(jmm1(k)).and.imm(k)>0.and.jmm1(k)>0.and.imm(k)<=mp.and.jmm1(k)+2<=np) &
            a(imm(k),jmm1(k))=a(imm(k),jmm1(k)+2)
       if(jmm2(k)==ip(jmm2(k)).and.imm(k)>0.and.jmm2(k)>2.and.imm(k)<=mp.and.jmm2(k)<=np) &
            a(imm(k),jmm2(k))=a(imm(k),jmm2(k)-2)
    enddo

    return

  end subroutine update_bounds
  ! =====================================================================================
  subroutine sweepx_fct

    integer::kn,i,j,i1,i2,ide,idid

    real::X

    real::U1(0:mp),U2(0:mp),W1(0:mp),W2(0:mp)

    U1=0.
    U2=0.
    W1=0.
    W2=0.

    do kn=1,nn
       j=jn(kn)
       i1=in1(kn)
       i2=in2(kn)-2
       if(ip(i1)/=i1.and.ip(i2)==i2) ide=1
       if(ip(i1)==i1.and.ip(i2)==i2) ide=2
       if(ip(i1)==i1.and.ip(i2)/=i2) ide=3
       if(ip(i1)/=i1.and.ip(i2)/=i2) ide=4

       if(ide==1.or.ide==3) i1=i1+1

       lp: do  i=i1,i2,2

          if(i > i2) exit lp

          idid=0
          if(ide==2.and.i==i1) idid=1
          if(i>2.and.i<mp-2)then
             call gegeX(i,j,ide)
          else
             goto 40
          endif
          if(ide==2.and.i==i1) goto 50
          if(ide==4.and.i==i1) goto 60
          if(ide==3.or.ide==4) goto 45

          if(i==i1) U2(i-2)=aa(i-1,j)
          if(i==i1) W2(i-2)=0.0
          X=g1+g2*W2(i-2)-g3*g5
          U1(i)=(g4+g2*U2(i-2)-g3*g7)/X
          W1(i)=g3*g6/X
          U2(i)=(g7*(g1+g2*W2(i-2))-g5*(g4+g2*U2(i-2)))/X
          W2(i)=g6*(g1+g2*W2(i-2))/X
          goto 40

45        if(i==i1) U1(i-2)=aa(i-1,j)
          if(i==i1) W1(i-2)=0.0
          X=g1*(1-g5*W1(i-2))+g2*g6
          U1(i)=(g2*(g7-g5*U1(i-2))+g4*(1-g5*W1(i-2)))/X
          W1(i)=g3*(1-g5*W1(i-2))/X
          U2(i)=(g1*(g7-g5*U1(i-2))-g6*g4)/X
          W2(i)=g3*g6/X
          goto 40

50        U2(i)=g7-g5*aa(i,j)
          W2(i)=g6
          U1(i)=0.0
          W1(i)=0.0
          goto 40

60        X= g1
          U1(i)=(g4+g2*aa(i,j))/X
          W1(i)=g3/X
40        continue
       enddo lp

       !  Computes UU and EE values
       lp2: do i=i2,i1,-2
          if(i < i1) exit lp2
          if(i2>mp-2) goto 70
          if(ide==4.and.i==i1) goto 80
          if(ide==2.and.i==i1) goto 90
          if(ide==3.or.ide==4) goto 85

          aa(i,j)=U1(i)+W1(i)*aa(i+2,j)
          aa(i+1,j)=U2(i)-W2(i)*aa(i+2,j)
          goto 70
85        aa(i+1,j)=U1(i)-W1(i)*aa(i+2,j)
          aa(i,j)=U2(i)+W2(i)*aa(i+2,j)
          goto 70
80        aa(i+1,j)=U1(i)-W1(i)*aa(i+2,j)
          goto 70
90        aa(i+1,j)=U2(i)-W2(i)*aa(i+2,j)
70        continue

       enddo lp2
    enddo
    call calculateV

    do i=1,mp
       do j=1,np
          if(i==ip(i).and.j/=ip(j)) goto 150
          a(i,j)=aa(i,j)
150       continue
       enddo
    enddo

    return

  end subroutine sweepx_fct
  ! =====================================================================================
  subroutine sweepy_fct

    integer::kn,i,j,j1,j2,ide,idid

    real::X

    real::U1(0:np),U2(0:np),W1(0:np),W2(0:np)

    U1=0.0
    U2=0.0
    W1=0.0
    W2=0.0

    do kn=1,mm

       i=im(kn)
       j1=jm1(kn)
       j2=jm2(kn)-2

       if(ip(j1)==j1.and.ip(j2)/=j2) ide=1
       if(ip(j1)/=j1.and.ip(j2)/=j2) ide=2
       if(ip(j1)/=j1.and.ip(j2)==j2) ide=3
       if(ip(j1)==j1.and.ip(j2)==j2) ide=4

       if(ide==1.or.ide==3) j1=j1+1

       do j=j1,j2,2

          idid=0
          if(ide==2.and.j==j1) idid=1

          if(i>2.and.i<mp-2)then
             call gegeY(i,j,ide)
          else
             goto 40
          endif

          if(ide==2.and.j==j1) goto 50
          if(ide==4.and.j==j1) goto 60
          if(ide==3.or.ide==4) goto 45

          if(j==j1) U2(j-2)=aa(i,j-1)
          if(j==j1) W2(j-2)=0.

          X=g1+g2*W2(j-2)-g3*g5
          U1(j)=(g4+g2*U2(j-2)-g3*g7)/X
          W1(j)=g3*g6/X
          U2(j)=(g7*(g1+g2*W2(j-2))-g5*(g4+g2*U2(j-2)))/X
          W2(j)=g6*(g1+g2*W2(j-2))/X
          goto 40

45        if(j==j1) U1(j-2)=aa(i,j-1)
          if(j==j1) W1(j-2)=0.
          X=g1*(1-g5*W1(j-2))+g2*g6
          U1(j)=(g2*(g7-g5*U1(j-2))+g4*(1-g5*W1(j-2)))/X
          W1(j)=g3*(1.0-g5*W1(j-2))/X
          U2(j)=(g1*(g7-g5*U1(j-2))-g6*g4)/X
          W2(j)=g3*g6/X
          goto 40

50        U2(j)=g7-g5*aa(i,j)
          W2(j)=g6
          U1(j)=0.
          W1(j)=0.
          goto 40

60        X=g1
          U1(j)=(g4+g2*aa(i,j))/X
          W1(j)=g3/X
40        continue

       enddo

       ! Computes the VV and ee values
       do j=j2,j1,-2

          if((ide==4).and.(j==j1)) goto 80
          if((ide==2).and.(j==j1)) goto 90
          if((ide==3).or. (ide==4)) goto 85

          aa(i,j)=U1(j)+W1(j)*aa(i,j+2)
          aa(i,j+1)=U2(j)-W2(j)*aa(i,j+2)
          goto 70
85        aa(i,j+1)=U1(j)-W1(j)*aa(i,j+2)
          aa(i,j)=U2(j)+W2(j)*aa(i,j+2)
          goto 70
80        aa(i,j+1)=U1(j)-W1(j)*aa(i,j+2)
          goto 70
90        aa(i,j+1)=U2(j)-W2(j)*aa(i,j+2)
70        continue
       enddo
    enddo

    call calculateU

    do i=1,mp
       do j=1,np
          if(i==ip(i).and.j/=ip(j)) goto 150
          a(i,j)=aa(i,j)
150       continue
       enddo
    enddo

    return

  end subroutine sweepy_fct
  ! =====================================================================================
  subroutine filter

    integer::i,j,k,i1,i2,j1,j2
    real::S

    S=0.1

    do k=1,nn
       i1=in1(k)
       i2=in2(k)
       if(ip(i1)/=i1) i1=i1+1
       if(ip(i2)/=i2) i2=i2-1
       j=jn(k)
       do i=i1,i2,2
          if(i>2.and.i<=mp-2)then
             aa(i,j)=a(i,j)+(S/2.0)*(1.0-S)*(a(i+2,j)+a(i-2,j)+a(i,j-2)+a(i,j+2) &
                -4.0*a(i,j))+(S*S/4.0)*(a(i-2,j+2)+a(i+2,j+2)+a(i+2,j-2) &
                 +a(i-2,j-2)-4.0*a(i,j))
          endif
       enddo
    enddo

    do k=1,mm
       j1=jm1(k)
       j2=jm2(k)
       if(ip(j1)==j1) j1=j1+1
       if(ip(j2)==j2) j2=j2-1
       i=im(k)
       do j=j1,j2,2
          if(i>2.and.i<=mp-2)then
             aa(i,j)=a(i,j)+(S/2.0)*(1.0-S)*(a(i+2,j)+a(i-2,j)+a(i,j-2)+a(i,j+2) &
                -4.0*a(i,j))+(S*S/4.0)*(a(i-2,j+2)+a(i+2,j+2)+a(i+2,j-2) &
                 +a(i-2,j-2)-4.0*a(i,j))
          endif
       enddo
    enddo

    do k=1,nn
       i1=in1(k)+2
       i2=in2(k)-2
       if(ip(i1)==i1) i1=i1+1
       if(ip(i2)==i2) i2=i2-1
       j=jn(k)
       do i=i1,i2,2
          if(i>2.and.i<=mp-2)then
             aa(i,j)=a(i,j)+(S/2.0)*(1.0-S)*(a(i+2,j)+a(i-2,j)+a(i,j-2)+a(i,j+2) &
                -4.0*a(i,j))+(S*S/4.0)*(a(i-2,j+2)+a(i+2,j+2)+a(i+2,j-2) &
                 +a(i-2,j-2)-4.0*a(i,j))
          endif
       enddo
    enddo
    do i=1,mp
       do j=1,np
          if(i==ip(i).and.j/=ip(j)) goto 70
          a(i,j)=aa(i,j)
70        continue
       enddo
    enddo

    return

  end subroutine filter
  ! =====================================================================================
  subroutine gegeY(i,j,ide)

    integer::i,j,ide,jaUXi

    real::S,b,dtx,ViXm,VimOD,Ca,uhm,vhm
    real::TSX,TSY,D,DDD,Um,f1,f2,f3,f4,f5,f6,f7,f8
    real::f9,f10,f11,f12,f13,f14,f15,f16,f17

    S=alfa
    b=beta
    dtx=dt/dx

    if(ide==3.or.ide==4) then
       jaUXi=j
       j=j+1
    endif

    ViXm=0.25*(p(i+1,j-1)+p(i+1,j+1)+p(i-1,j+1)+p(i-1,j-1))

    VimOD=sqrt(p(i,j)*p(i,j)+ViXm*ViXm)
    VimOD=VimOD*51.44

    Ca=0.0012875
    if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD

    TSX=0.001*Ca*51.44*ViXm*VimOD
    TSY=0.001*Ca*51.44*p(i,j)*VimOD
    TSX=TSX*alin
    TSY=TSY*alin

    D=0.5*(a(i-1,j)+a(i+1,j)+a(i,j+1)+a(i,j-1))
    if(D <= 0.) return

    DDD=0.5*(a(i-1,j)+a(i+1,j))

    Um=0.25*(a(i-1,j-1)+a(i+1,j-1)+a(i-1,j+1)+a(i+1,j+1))

    f1=0.5*S* a(i,j)/D

    f2=0.125*f1*dtx*(a(i+1,j+2)-a(i+1,j-2)+a(i-1,j+2) &
       -a(i-1,j-2))

    f3=0.25*S*dtx*Um*(a(i+1,j)-a(i-1,j))/D

    f4=0.0625*S*dtx*Um*(a(i+2,j+1)-a(i-2,j+1)+&
         a(i+2,j-1)-a(i-2,j-1))/D

    f5=0.125*S*dtx*(a(i+1,j+1)-a(i-1,j+1)+&
         a(i+1,j-1)-a(i-1,j-1))

    f6=0.125*S*dtx*Um*(a(i+2,j)-a(i-2,j))
    f7=0.5*dt*f*Um
    f8=0.25*dtx*(grav-S*(a(i,j)**2)/D)
    f9=0.5*dt*TSY/(rho*D)

    uhm=Um
    vhm=a(i,j)

    f10=0.5*dt*cfric*sqrt(uhm*uhm+vhm*vhm)/D
    f11=0.25*dtx*gamma*(p(i,j+1)-p(i,j-1))/(rho/1000.)
    f11=f11*alin
    f12=0.125*dtx*ah*(a(i+2,j)-2*a(i,j)+a(i-2,j))/dx
    f13=0.125*dtx*ah*(a(i,j+2)-2*a(i,j)+a(i,j-2))/dx

    if(ide==3.or.ide==4)then
       j=jaUXi
       j=j-1
    endif

    f14=0.0625*dtx*(a(i+1,j+1)+a(i-1,j+1))*(a(i+1,j+2) &
       -a(i-1,j+2)+a(i+1,j)-a(i-1,j)+a(i+2,j+1)-&
         a(i-2,j+1))

    f15=0.0625*dtx*(a(i+1,j+2)-a(i+1,j)+a(i-1,j+2) &
       -a(i-1,j)+a(i,j+3)-a(i,j-1))

    f16=0.0625 *dtx*(4*a(i,j+1)+a(i+1,j+2)+a(i+1,j)+&
         a(i-1,j)+a(i-1,j+2))*(a(i+1,j+1)-a(i-1,j+1))

    f17=0.0625*dtx*(4*a(i,j+1)+a(i+1,j+2)+a(i+1,j) &
        +a(i-1,j)+a(i-1,j+2))

    g1=1.0-b*(f2+f3+f4+f5-f10)
    g2=f1+b*f8
    g3=b*f8-f1

    if(ide==3.or.ide==4)then
       j=jaUXi
       j=j+1
    endif

    g4=a(i,j)*(1.0+(1.0-b)*(f2+f3+f4+f5-f10))+&
         a(i,j-1)*((1.0-b)*f8-f1)-a(i,j+1)*((1.0-b)*f8+&
         f1)-f6-f7+f9-f11+f12+f13
    g5=b*(f15-f17)
    g6=b*(f15+f17)

    if(ide==3.or.ide==4)then
       j=jaUXi
       j=j-1
    endif

    g7=a(i,j+1)+a(i,j)*(1.0-b)*(f17-f15)-a(i,j+2)*&
         (1.0-b)*(f15+f17)-f14-f16

    if(ide==3.or.ide==4) j=jaUXi

    return

  end subroutine gegeY
  ! =====================================================================================
  subroutine gegeX(i,j,ide)

    integer::i,j,ide,iaUXi

    real::S,b,dtx,ViYm,VimOD,Ca,uhm,vhm
    real::TSX,TSY,D,DDD,Vm,f1,f2,f3,f4,f5,f6,f7,f8
    real::f9,f10,f11,f112,f113,f12,f13,f14,f15

    S=alfa
    b=beta
    dtx=dt/dx

    if(ide==3.or.ide==4)then
       iaUXi=i
       i=i+1
    endif

    ViYm=0.25*(p(i+1,j-1)+p(i+1,j+1)+p(i-1,j+1)+p(i-1,j-1))
    VimOD=sqrt(p(i,j)*p(i,j)+ViYm*ViYm)
    VimOD=VimOD*51.44

    Ca=0.0012875
    if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD

    TSX=0.001*Ca*51.44*p(i,j)*VimOD
    TSY=0.001*Ca*51.44*ViYm*VimOD
    TSX=TSX*alin
    TSY=TSY*alin

    D=0.5*(a(i,j-1)+a(i,j+1)+a(i+1,j)+a(i-1,j))
    if(D <= 0.) return

    DDD=0.5*(a(i,j-1)+a(i,j+1))
    Vm=0.25*(a(i-1,j-1)+a(i+1,j-1)+a(i-1,j+1)+a(i+1,j+1))

    f1=0.5*S*a(i,j)/D

    f2=0.125*f1*dtx*(a(i+2,j-1)-a(i-2,j-1)+a(i+2,j+1) &
       -a(i-2,j+1))
    f3=0.25*S*dtx*Vm*(a(i,j+1)-a(i,j-1))/D

    f4=0.0625*S*dtx*Vm*(a(i+1,j+2)-a(i+1,j-2)+&
         a(i-1,j+2)-a(i-1,j-2))/D
    f5=0.125*S*dtx*Vm*(a(i,j+2)-a(i,j-2))

    f6=0.125*S*dtx *(a(i+1,j+1)-a(i+1,j-1)+a(i-1,j+1) &
       -a(i-1,j-1))

    f7=0.5*dt*f*Vm
    f8=0.25*dtx*(grav-S*((a(i,j))**2)/D)
    f9=0.5*dt*TSX/(rho*D)

    uhm=a(i,j)
    vhm=Vm

    f10=0.5*dt*cfric*sqrt(uhm*uhm+vhm*vhm)/D
    f11=0.25*dtx*gamma*(p(i+1,j)-p(i-1,j))/(rho/1000.)
    f11=f11*alin
    f112= 0.125*dtx*ah*(a(i+2,j)-2.0*a(i,j)+a(i-2,j))/dx
    f113= 0.125*dtx*ah*(a(i,j+2)-2.0*a(i,j)+a(i,j-2))/dx

    if(ide==3.or.ide==4)then
       i=iaUXi
       i=i-1
    endif

    f12=0.0625*dtx*(a(i+2,j-1)-a(i,j-1)+a(i+2,j+1) &
       -a(i,j+1)+a(i+3,j)-a(i-1,j))

    f13=0.0625*dtx*(4*a(i+1,j)+a(i,j+1)+a(i,j-1) &
        +a(i+2,j+1)+a(i+2,j-1))

    f14=0.0625*dtx*(a(i+1,j+1)+a(i+1,j-1))*(a(i+2,j+1) &
       -a(i+2,j-1)+a(i,j+1)-a(i,j-1)+a(i+1,j+2)-&
         a(i+1,j-2))

    f15=0.0625*dtx*(4*a(i+1,j)+a(i+2,j+1)+a(i+2,j-1) &
        +a(i,j+1)+a(i,j-1))*(a(i+1,j+1)-a(i+1,j-1))

    g1=1-b*(f2+f3+f4+f6-f10)
    g2=f1+b*f8
    g3=b*f8-f1

    if(ide==3.or.ide==4) then
       i=iaUXi
       i=i+1
    endif

    g4=a(i,j)*(1.0+(1.0-b)*(f2+f3+f4+f6-f10))+&
         a(i-1,j)*((1.0-b)*f8-f1)-a(i+1,j)*((1.0-b)*f8+&
         f1)-f5+f7+f9-f11+f112+f113
    g5=b*(f12-f13)
    g6=b*(f12+f13)

    if(ide==3.or.ide==4) then
       i=iaUXi
       i=i-1
    endif

    g7=a(i+1,j)+a(i,j)*(1.0-b)*(f13-f12)-a(i+2,j)*&
         (1.0-b)*(f12+f13)-f14-f15

    if(ide==3.or.ide==4) i=iaUXi

    return

  end subroutine gegeX
  ! =====================================================================================
  subroutine calculateV

    integer::j,i,kn,iSen1,iSen2,j1,j2,ideS1,ini,ifi

    real::S,b,dtx,Um,D,DDD,f1,f2,f3,f4,f5,f6,f7,f8,f9

    real::ViXm,VimOD,Ca,X,TSX,TSY,uhm,vhm

    real::U1(0:np),W1(0:np)

    U1=0.
    W1=0.

    dtx=dt/dx
    S=alfa
    b=beta

    do kn=1,mm
       iSen1=0
       iSen2=0
       i=im(kn)
       if(i<=2.or.i>mp-2) goto 10
       j1= jm1(kn)+2
       j2= jm2(kn)-2
       ideS1=jm2(kn)-jm1(kn)
       if(ideS1==2) goto 10

       ini=j1
       ifi=j2
       if(ip(j1)/=j1) goto 20

       j1=j1-1
       ini=j1
       iSen1=1

20     if(ip(j2)/=j2) goto 30
       j2=j2+1
       ifi=j2
       iSen2=1

30     j=ini
35     if(iSen1==0.and.iSen2==1.and.j==ini) j=ifi

       Um=0.25*(b*(aa(i+1,j+1)+aa(i+1,j-1)+aa(i-1,j-1)+&
            aa(i-1,j+1))+(1.0-b)*(a(i+1,j+1)+a(i+1,j-1)+&
            a(i-1,j-1)+a(i-1,j+1)))
       D=0.5*(a(i+1,j)+a(i-1,j)+(1.0-b)*(a(i,j+1)+a(i,j-1)) &
           +b*(aa(i,j+1)+aa(i,j-1)))
       if(D <= 0.) goto 10
       DDD=0.5*(a(i+1,j)+a(i-1,j))

       f1=0.125*dtx*S*Um*(a(i+2,j)-a(i-2,j))
       f2=0.125*S*dtx*a(i,j)
       f3=0.5*dt*f*Um
       f4=0.25*dtx*gamma*(p(i,j+1)-p(i,j-1))/(rho/1000.0)
       f4=f4*alin
       f5=0.25*dtx*grav*(b*(aa(i,j+1)-aa(i,j-1))+(1.0-b)*&
            (a(i,j+1)-a(i,j-1)))

       ViXm=0.25*(p(i+1,j-1)+p(i+1,j+1)+p(i-1,j+1) &
           +p(i-1,j-1))
       VimOD=sqrt(p(i,j)*p(i,j)+ViXm*ViXm)
       VimOD=VimOD*51.44

       Ca=0.0012875
       if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD

       TSX=0.001*Ca*51.44*ViXm*VimOD
       TSY=0.001*Ca*51.44*p(i,j)*VimOD
       TSX=TSX*alin
       TSY=TSY*alin

       f6=0.5*dt*TSY/(rho*D)

       uhm=Um
       vhm=b*aa(i,j)+(1.0-b)*a(i,j)

       f7=0.5*dt*cfric*sqrt(uhm*uhm+vhm*vhm)/D
       f8=0.125*ah*dtx*(a(i+2,j)-2*a(i,j)+a(i-2,j)) /dx
       f9=0.125*ah*dtx*(1.0-b)*(a(i,j+2)-2*a(i,j)+&
            a(i,j-2))/dx

       if(iSen1==1.and.j==ini) goto 50
       if(iSen2==1.and.j==ifi) goto 50

       goto 60

50     aa(i,j)=(a(i,j)*(1.0-(1.0-b)*f7)-f3-f4-f5+f6)/(1.0+b*f7)

       if(iSen1==1.and.j==ini)then
          ini=ini+2
          iSen1=0
          goto 30
       else
          ifi=ifi-2
          iSen2=0
          goto 30
       endif

60     g1=- b*(f2+0.125*dtx*ah/dx)
       g2=1.0+b*(f7+0.25*dtx*ah/dx)
       g3=b*(f2-0.125*dtx*ah/dx)
       g4=a(i,j)-f1-(1.0-b)*f2*(a(i,j+2)-a(i,j-2))-&
            f3-f4-f5+f6-f7*(1.0-b)*a(i,j)+f8+f9

       if(j/=ini) goto 55
       U1(j-2)=0.0
       W1(j-2)=aa(i,j-2)

55     X=g2-U1(j-2)*g1
       U1(j)=g3/X
       W1(j)=(g4-W1(j-2)*g1)/X
       if(j >= ifi) goto 15
       j=j+2
       goto 35

       ! Computes the VV values
15     do j=ifi,ini,-2
          aa(i,j)=W1(j)-U1(j) *aa(i,j+2)
       enddo

10     continue

    enddo

    return

  end subroutine calculateV
  ! =====================================================================================
  subroutine calculateU

    integer::j,i,kn,iSen1,iSen2,i1,i2,ideS1,ini,ifi

    real::S,b,dtx,Vm,D,DDD,f1,f2,f3,f4,f5,f6,f7,f8,f9
    real::ViYm,VimOD,Ca,X,TSX,TSY,uhm,vhm

    real::U1(0:mp),W1(0:mp)

    U1=0.
    W1=0.

    dtx=dt/dx
    S=alfa
    b=beta

    do kn=1,nn
       iSen1=0
       iSen2=0
       i1=in1(kn)+2
       i2=in2(kn)-2
       j=jn(kn)
       ideS1=in2(kn)-in1(kn)
       if(ideS1==2) goto 10
       ini=i1
       ifi=i2
       if(ip(i1)==i1) goto 20
       i1=i1-1
       ini=i1
       iSen1=1
20     if(ip(i2)==i2) goto 30
       i2=i2+1
       ifi=i2
       iSen2=1
30     i=ini

35     if(iSen1==0.and.iSen2==1.and.i==ini) i=ifi

       if(i<=2.or.i>mp-2) goto 10

       Vm=0.25*(b*(aa(i+1,j+1)+aa(i+1,j-1)+aa(i-1,j-1)+&
            aa(i-1,j+1))+(1-b)*(a(i+1,j+1)+a(i+1,j-1)+&
            a(i-1,j-1)+a(i-1,j+1)))

       D=0.5*(a(i,j+1)+a(i,j-1)+(1-b)*(a(i+1,j)+a(i-1,j)) &
           +b*(aa(i+1,j)+aa(i-1,j)))

       if(D <= 0.) goto 10
       DDD=0.5*(a(i,j+1)+a(i,j-1))

       f1=0.125*S*dtx*a(i,j)
       f2=0.125*S*dtx*Vm*(a(i,j+2)-a(i,j-2))
       f3=0.5*dt*f*Vm

       ViYm=0.25*(p(i+1,j-1)+p(i+1,j+1)+p(i-1,j+1) &
           +p(i-1,j-1))
       VimOD=sqrt(p(i,j)*p(i,j)+ViYm*ViYm)
       VimOD=VimOD*51.44
       Ca=0.0012875

       if(VimOD > 750.) Ca=0.0008+6.5e-07*VimOD
       TSX=0.001*Ca*51.44*p(i,j)*VimOD
       TSY=0.001*Ca*51.44*ViYm*VimOD
       TSX=TSX*alin
       TSY=TSY*alin

       f4=0.5*dt*TSX/(rho*D)

       uhm=b*aa(i,j)+(1.0-b)*a(i,j)
       vhm=Vm

       f5=0.5*dt*cfric*sqrt(uhm*uhm+vhm*vhm)/D
       f6=0.25*dtx*grav*(b*(aa(i+1,j)-aa(i-1,j))+(1.0-b)*&
            (a(i+1,j)-a(i-1,j)))
       f7=0.25*dtx*gamma*(p(i+1,j)-p(i-1,j))/(rho/1000.0)
       f7=f7*alin
       f8=0.125*ah*dtx*(1.0-b)*(a(i+2,j)-2.0*a(i,j)+a(i-2,j))/dx
       f9=0.125*ah*dtx *(a(i,j+2)-2.0*a(i,j)+a(i,j-2))/dx

       if(iSen1==1.and.i==ini) goto 50
       if(iSen2==1.and.i==ifi) goto 50
       goto 60

50     aa(i,j)=(a(i,j)*(1.0-(1.0-b)*f5)+f3+f4-f6-f7)/(1.0+b*f5)
       if(iSen1==1.and.i==ini) then
          ini=ini+2
          iSen1=0
          goto 30
       else
          ifi=ifi-2
          iSen2=0
          goto 30
       endif

60     g1=-b*(f1+0.125*dtx*ah/dx)
       g2=1.0+b*(f5+0.25*dtx*ah/dx)
       g3=b*(f1-0.125*dtx*ah/dx)
       g4=a(i,j)-f2-(1.0-b)*f1*(a(i+2,j)-a(i-2,j))+&
            f3+f4-f6-f7-f5*(1.0-b)*a(i,j)+f8+f9
       if(i/=ini) goto 55

       U1(i-2)=0.0
       W1(i-2)=aa(i-2,j)

55     X=g2-U1(i-2)*g1
       U1(i)=g3/X
       W1(i)=(g4-W1(i-2)*g1)/X
       if(i >= ifi) goto 15

       ! Computes the UU values
       i=i+2
       goto 35
15     do i=ifi,ini,-2
          aa(i,j)=W1(i)-U1(i)*aa(i+2,j)
       enddo

10     continue
    enddo

    return

  end subroutine calculateU
  ! =====================================================================================

end module currents_circ
