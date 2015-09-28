! =====================================================================================
! ASCIDIAN
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple
! Place, Suite 330, Boston, MA 02111-1307 USA
! =====================================================================================
! =====================================================================================
!
!       Filename:  CurrentIO.f90
!
!    Description:  Encapsulates the parameter module used for tides and winds induced
!
!        Version:  1.0
!        Created:  19/09/15 15:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module currents_io

  use mpidata
  use xml_reader
  use interpolation
  use parallel_circ
  use currents_data
  use precision_data

  implicit none

  ! Computing matrix
  real,dimension(:,:),allocatable::Umean,Vmean
  real,dimension(:,:),allocatable::Uave,Vave
  real,dimension(:,:),allocatable::pUmean,pVmean
  real,dimension(:,:),allocatable::lUmean,lVmean

  ! Forces parameters (wind and wave)
  real::sea_level,tnow,toutput

  real,dimension(:),allocatable::x1_array,y1_array
  real,dimension(:),allocatable::xc_array,yc_array
  real,dimension(:,:),allocatable::o1,o2

  ! Processes flags
  logical::circon,waveon,spmon

  public

contains

  ! =====================================================================================
  subroutine currents_init

    logical::found
    integer::i,j,l1,scn
    character(len=128)::command,stg,fildir
    real::ome,tiner,x,y,z,hmax

    pi=4.*atan(1.)
    maxdepth=50000.0

    ! Open input file
    spmon=.false.
    circon=.false.
    waveon=.false.
    outputflag=.false.
    outdir=' '
    outdir='outOcean'
    initDep=0.0
    circflag=0
    waveflag=0
    spmflag=0
    call xml_Data
    ! Conversion in cm
    maxdepth=maxdepth*100.
    if(circflag==1) circon=.true.
    if(waveflag==1) waveon=.true.
    if(spmflag==1) spmon=.true.
    if(spmon) initDep=0.0

    ! The coefficient beta takes any value between 0 and 1.
    ! A value of beta > 0.5 permits theoretically
    ! the use of higher Courant numbers.
    beta=0.6666666
    ! Value of 0 or 1. If 0 the atmospheric pressure gradients are
    ! not considered
    gamma=0.0
    ! SPM grid extent
    stratal_x=sp_n
    stratal_y=sp_m
    st_dx=stratal_dx
    ! Conversion in cm
    dx=stratal_dx*100.
    stratal_dx=dx
    nbPt=sp_n*sp_m

    ! Circulation computational grid
    n=int((sp_n*stratal_dx)/dx)+gap*2
    m=int((sp_m*stratal_dx)/dx)+gap*2

    ! Get maximum bathymetry value
    open(unit=17,file=trim(xyzfile))
    rewind(17)
    hmax=0.0_8
    do i=1,sp_m*sp_n
      read(17,*)j,x,y,z
      if(z<0.0) hmax=max(hmax,-z)
    enddo
    close(17)

    if(.not.spmon) hmax=hmax-initDep
    grav=981.0

    if(hmax==0.0)then
      dt=60.0
    else
      dt=Courant/(sqrt(grav/100*hmax)*sqrt(2/(st_dx**2)))
    endif
    ! Horizontal eddy viscosity coefficient (cm2/sec)
    ! and whose value should be about U L, where
    !   - U is a typical velocity and
    !   - L is horizontal grid spacing
    ah=5.*dx
    atf=100.*dt

    ! Time steps for circulation computation
    kend=int(kend*3600/dt)

    rho=1.0
    alat=-alat
    ome=2.*pi/(24.*3600.)
    f=2.*ome*sin(alat*pi/180.)
    tiner=2.*pi/abs(f)
    niner=tiner/dt

    ! Allocate forecasts
    forecast_param(1)=hindcast(1)%subgroup(1)%hs
    forecast_param(2)=hindcast(1)%subgroup(1)%per
    forecast_param(3)=hindcast(1)%subgroup(1)%dir
    forecast_param(4)=hindcast(1)%subgroup(1)%dd
    forecast_param(5)=hindcast(1)%subgroup(1)%wvel
    forecast_param(6)=hindcast(1)%subgroup(1)%wdir
    forecast_param(7)=forecast_param(5)*sin(forecast_param(6)*180./pi)
    forecast_param(8)=forecast_param(5)*cos(forecast_param(6)* 180./pi)

    ! Build partition
    call define_band_partitioning
    call mpi_barrier(ocean_comm_world,ierr)

    ! Allocate variables
    if(iam==0)then
      if(allocated(a)) deallocate(a)
      if(allocated(p)) deallocate(p)
      if(allocated(aa)) deallocate(aa)
      allocate(a(mp,np),aa(mp,np),p(mp,np))
    endif

    if(allocated(seafloor)) deallocate(seafloor)
    if(allocated(sp_topo)) deallocate(sp_topo)
    allocate(sp_topo(sp_m,sp_n),seafloor(m,n))

    if(allocated(Uave)) deallocate(Uave)
    if(allocated(Vave)) deallocate(Vave)
    allocate(Uave(sp_m,sp_n),Vave(sp_m,sp_n))
    Uave=0.
    Vave=0.

    if(iam==0)then
       ! if(allocated(pUmean)) deallocate(pUmean)
       ! if(allocated(pVmean)) deallocate(pVmean)
       ! allocate(pUmean(mp,np),pVmean(mp,np))
       ! if( allocated(lUmean)) deallocate(lUmean)
       ! if( allocated(lVmean)) deallocate(lVmean)
       ! allocate(lUmean(m,n),lVmean(m,n))
    endif
    if(allocated(Umean)) deallocate(Umean)
    if(allocated(Vmean)) deallocate(Vmean)
    allocate(Umean(m,n),Vmean(m,n))

    ! Arrays for interpolation
    if(allocated(x1_array)) deallocate(x1_array)
    if(allocated(y1_array)) deallocate(y1_array)
    allocate(x1_array(n),y1_array(m))

    if(allocated(xc_array)) deallocate(xc_array)
    if(allocated(yc_array)) deallocate(yc_array)
    allocate(xc_array(n),yc_array(m))

    ylen=0
    yclen=0
    x1_array=0.
    y1_array=0.
    xc_array=0.
    yc_array=0.
    do i=1,m
      if(i==ip(i))then
        ylen=ylen+1
        y1_array(ylen)=(i-1)*dx/100.
      endif
      yclen=yclen+1
      yc_array(yclen)=(i-1)*dx/100.
    enddo

    xlen=0
    xclen=0
    do j=1,n
      if(j==ip(j))then
        xlen=xlen+1
        x1_array(xlen)=(j-1)*dx/100.
      endif
      xclen=xclen+1
      xc_array(xclen)=(j-1)*dx/100.
    enddo

    if(allocated(o1)) deallocate(o1)
    if(allocated(o2)) deallocate(o2)
    allocate(o1(xlen,ylen),o2(xlen,ylen))
    if(iam==0)then
      if(allocated(p1)) deallocate(p1,p2,p3,p4,p5,p6,p7)
      allocate(p1(np*np,3),p2( mp*mp,3),p3(np*np,3),p4(mp*mp,3))
      allocate(p5(2*mp+2*np,3),p6(np*np,3),p7(np*np,3))
    endif
    if(allocated(sea_floor)) deallocate(sea_floor)
    allocate(sea_floor(nbPt,3))

    ! Create the output directory
    outdir1=''
    outdir2=''
    if(iam==0)then
      call noblnk(outdir)
      fildir=outdir
      stg='/outputs/checkout-fc1-fg1-id1.csv'
      call noblnk(stg)
      call append_str(fildir,stg)
      call noblnk(fildir)
      i=len_trim(fildir)
      inquire(file=fildir(1:i),exist=found)
      if(.not.found)then
        command=' '
        command(1:6)='mkdir '
        l1=len_trim(outdir)
        command(7:l1+7)=outdir
        call term_command(command)
        command(l1+7:l1+7)='/'
        stg=''
        stg(1:l1+7)=command(1:l1+7)
        stg(l1+8:l1+13)='inputs'
        call term_command( stg )
        stg=''
        stg(1:l1+7)=command(1:l1+7)
        stg(l1+8:l1+14)='outputs'
        call term_command(stg)
        outdir1(1:l1)=outdir(1:l1)
        outdir1(l1+1:l1+7)='/inputs'
        outdir2(1:l1)=outdir(1:l1)
        outdir2(l1+1:l1+8)='/outputs'
        call noblnk(outdir1)
        call noblnk(outdir2)
       else
        command=' '
        command(1:6)='rm -r '
        l1=len_trim(outdir)
        command(7:l1+7)=outdir
        call term_command(command)
        command=' '
        command(1:6)='mkdir '
        l1=len_trim(outdir)
        command(7:l1+7)=outdir
        call term_command(command)
        command(l1+7:l1+7)='/'
        stg=''
        stg(1:l1+7)=command(1:l1+7)
        stg(l1+8:l1+13)='inputs'
        call term_command(stg)
        stg=''
        stg(1:l1+7)=command(1:l1+7)
        stg(l1+8:l1+14)='outputs'
        call term_command(stg)
        outdir1(1:l1)=outdir(1:l1)
        outdir1(l1+1:l1+7)='/inputs'
        outdir2(1:l1)=outdir(1:l1)
        outdir2(l1+1:l1+8)='/outputs'
        call noblnk(outdir1)
        call noblnk(outdir2)
      endif
    endif
    call mpi_bcast(outdir,128,mpi_character,0,ocean_comm_world,ierr)
    call mpi_bcast(outdir1,128,mpi_character,0,ocean_comm_world,ierr)
    call mpi_bcast(outdir2,128,mpi_character,0,ocean_comm_world,ierr)

    if(iam==0.and.spmon)then
      maestro='poseidon'
      call noblnk(maestro)
      call addpath3(maestro)
      xyzfile='topSPM_surface.xyz'
      call noblnk(xyzfile)
      call addpath3(xyzfile)
    endif
    call mpi_bcast(maestro,128,mpi_character,0,ocean_comm_world,ierr)
    call mpi_bcast(xyzfile,128,mpi_character,0,ocean_comm_world,ierr)

    ! Allocate output arrays
    maxsubgroup=0
    do scn=1,forecast_nb
      maxsubgroup=max(maxsubgroup,hindcast(scn)%cnb)
    enddo
    totnodes=sp_m*sp_n
    allocate(nodes(3*totnodes))
    allocate(htopo(totnodes))
    allocate(curU(maxsubgroup,totnodes))
    allocate(curV(maxsubgroup,totnodes))
    allocate(wavU(maxsubgroup,totnodes))
    allocate(wavV(maxsubgroup,totnodes))

    return

  end subroutine currents_init
  ! =====================================================================================
  subroutine build_computationalgrid( fst )

    integer::i,j,li,lj,ig,jg,k,step,fst,ipp,procID

    real::x,y

    ! Initialise computational arrays
    if(iam==0)then
      a=0.0
      p=0.0
      aa=0.0
    endif

    ! Substract the sea level elevation to start with a realtive sea level
    ! elevation value
    open(unit=17,file=trim(xyzfile))
    rewind(17)
    sp_topo=0.
    step=int(dx/stratal_dx)
    li=1
    lj=1
    ig=0
    jg=0
    k=0
    seafloor=0.0
    do i=1,sp_m
       do j=1,sp_n
          if(spmon)then
             read(17,*)x,y,sp_topo(i,j)
             k=k+1
          else
             read(17,*)k,x,y,sp_topo(i,j)
          endif
          if(i==1.and.j==1)then
             stratal_xo=x
             stratal_yo=y
          elseif(i==sp_n.and.j==sp_m)then
             stratal_xm=x
             stratal_ym=y
          endif
          sea_floor(k,1)=x
          sea_floor(k,2)=y

          sp_topo(i,j)=sp_topo(i,j)+initDep-sea_level
          sea_floor(k,3)=sp_topo(i,j)
          if(ig==0.and.jg==0)then
             seafloor(li+gap,lj+gap)=sp_topo(i,j)
             if(sp_topo(i,j)<=-5900.0) seafloor(li+gap,lj+gap)=-5900.0
             lj=lj+1
          endif
          jg=jg+1
          if(jg==step) jg=0
       enddo
       ig=ig+1
       jg=0
       if(ig==step)then
          ig=0
          li=li+1
          lj=1
       endif
    enddo
    close(17)

    ! Update boundaries elevation
    do k=1,2
       do i=1,gap
          seafloor(i,1:n)=seafloor(gap+1,1:n)
          seafloor(m-i+1,1:n)=seafloor(m-gap,1:n)
       enddo
       do j=1,gap
          seafloor(1:m,j)=seafloor(1:m,gap+1)
          seafloor(1:m,n-j+1) = seafloor(1:m,n-gap)
       enddo
    enddo

    ! Specify the initial bathymetry.
    if(iam==0)then
       procID=iam
       do i=1,m
          ipp=i-cextent(procID+1,1)+1
          do j=1,n
             if(sea_level>seafloor(i,j))then
                ! Declare mean depth as positive if marine and
                ! convert it to centimetres.
                seafloor(i,j)=(sea_level-seafloor(i,j))*100.
             else
                seafloor(i,j)=0.
             endif
             if(i>=cextent(procID+1,1).and.i<=cextent(procID+1,2))then
                psea(ipp,j)=seafloor(i,j)
             endif
             if(i==ip(i).and.j/=ip(j))then
                if(i>=cextent(procID+1,1).and.i<=cextent(procID+1,2))then
                   a(ipp,j)=seafloor(i,j)
                endif
             endif
          enddo
       enddo
    endif

    ! Specify the meteorological forces.
    if(iam==0)then
       procID=iam
       do i=1,m
          ipp=i-cextent(procID+1,1)+1
          do j=1,n
             ! Get wind velocity along x which needs to be converted in knots
             if(i==ip(i).and.j==ip(j))then
                if(i>=cextent(procID+1,1).and.i<=cextent(procID+1,2))then
                   p(ipp,j)=forecast_param(7)*knots
                endif
                ! Get wind velocity along y which needs to be converted in knots
             elseif(i/=ip(i).and.j/=ip(j))then
                if(i>=cextent(procID+1,1).and.i<=cextent(procID+1,2))then
                   p(ipp,j)=forecast_param(8)*knots
                endif
             endif
          enddo
       enddo
    endif

    ! Build boundaries
    if(iam==0) call build_bounds

    return

  end subroutine build_computationalgrid
  ! =====================================================================================
  subroutine build_bounds

    logical::marine
    integer::i,j,st,minp1,minp2,halo

    halo=2

    ! Specify the positions i,j of the computational domain boundaries
    ! First compute boundaries for jn,in1,in2
    pt1=0
    p1=0
    minp1=n
    do j=gap+1,np-gap
       marine=.false.
       if(j==ip(j))then
          do i=gap+1,mp-gap
             if(seafloor(i,j)>0.and.seafloor(i,j)<maxdepth.and..not.marine.and.i/=mp-gap)then
                pt1=pt1+1
                p1(pt1,1)=j
                p1(pt1,2)=i
                marine=.true.
                minp1=min(minp1,j+1)
             elseif(marine.and.seafloor(i,j)<=0.)then
                marine=.false.
                p1(pt1,3)=i-1
             elseif(marine.and.i==mp-gap)then
                p1(pt1,3)=i
             endif
          enddo
       endif
    enddo

    ! Then compute boundaries for im,jm1,jm2
    pt2=0
    p2=0
    minp2=m
    do i=gap+1,mp-gap
       marine=.false.
       if(i/=ip(i))then
          do j=gap+1,np-gap
             if(seafloor(i,j)>0.and.seafloor(i,j)<maxdepth.and..not.marine.and.j/=np-gap)then
                pt2=pt2+1
                p2(pt2,1)=i
                p2(pt2,2)=j
                minp2=min(minp2,i+1)
                marine=.true.
             elseif(marine.and.seafloor(i,j)<=0.)then
                marine=.false.
                p2(pt2,3)=j-1
             elseif(marine.and.j==np-gap)then
                p2(pt2,3)=j
             endif
          enddo
       endif
    enddo

    ! Then compute boundaries for ina,inb,jnab
    pt3=0
    p3=0
    do j=gap+1,np-gap
       marine=.false.
       do i=gap+1,mp-gap
          if(seafloor(i,j)>0.and.seafloor(i,j)<maxdepth.and..not.marine.and.i/=mp-gap)then
             pt3=pt3+1
             p3(pt3,1)=j
             p3(pt3,2)=i-1
             marine=.true.
          elseif(marine.and.seafloor(i,j)<=0.)then
             marine=.false.
             p3(pt3,3)=i-1
          elseif(marine.and.i==mp-gap)then
             p3(pt3,3)=i-1
          endif
       enddo
    enddo
    ! Then compute boundaries for imab,jma,jmb
    pt4=0
    p4=0
    do i=gap+1,mp-gap
       marine=.false.
       do j=gap+1,np-gap
          if(seafloor(i,j)>0.and.seafloor(i,j)<maxdepth.and..not.marine.and.j/=np-gap)then
             pt4=pt4+1
             p4(pt4,1)=i
             p4(pt4,2)=j-1
             marine=.true.
          elseif(marine.and.seafloor(i,j)<=0.)then
             marine=.false.
             p4(pt4,3)=j
          elseif(marine.and.j==np-gap)then
             p4(pt4,3)=j
          endif
       enddo
    enddo
    ! Then compute boundaries for ie,je,ije
    ! Lower boundary
    pt5=0
    p5=0
    st=gap+1
    do j=gap+1,np-gap
       if(seafloor(st,j)>0.and.seafloor(st,j)<maxdepth.and.j==ip(j))then
          pt5=pt5+1
          p5(pt5,1)=st
          p5(pt5,2)=j
          p5(pt5,3)=-1
       endif
    enddo
    ! Upper bondary
    st=mp-gap
    do j=gap+1,np-gap
       if(seafloor(st,j)>0.and.seafloor(st,j)<maxdepth.and.j==ip(j))then
          pt5=pt5+1
          p5(pt5,1)=st
          p5(pt5,2)=j
          p5(pt5,3)=1
       endif
    enddo
    ! Left bondary
    st=gap+1
    do i=gap+1,mp-gap
       if(seafloor(i,st)>0.and.seafloor(i,st)<maxdepth.and.i==ip(i))then
          pt5=pt5+1
          p5(pt5,1)=i
          p5(pt5,2)=st
          p5(pt5,3)=-2
       endif
    enddo
    ! Right bondary
    st=np-gap
    do i=gap+1,mp-gap
       if(seafloor(i,st)>0.and.seafloor(i,st)<maxdepth.and.i==ip(i))then
          pt5=pt5+1
          p5(pt5,1)=i
          p5(pt5,2)=st
          p5(pt5,3)=2
       endif
    enddo
    ! In case there is no open boundary
    if(pt5==0)then
       pt5=1
       p5(1,1)=1
       p5(1,2)=1
       p5(1,3)=0
    endif
    ! Then compute boundaries for inn1,inn2,jnn
    pt6=0
    p6=0
    if(minp1==0)print*,'Problem defining jnn'
    minp1=gap+1
    minp2=gap+1
    do j=minp1,np-gap
       marine=.false.
       if(j/=ip(j))then
          do i=gap+1,mp-gap
             if(seafloor(i,j)>0.and.seafloor(i,j)<maxdepth.and. &
                  .not.marine.and.i/=mp-gap)then
                pt6=pt6+1
                p6(pt6,1)=j
                p6(pt6,2)=i
                marine=.true.
             elseif(marine.and.seafloor(i,j)<=0.)then
                marine=.false.
                p6(pt6,3)=i-1
             elseif(marine.and.i==mp-gap)then
                p6(pt6,3)=i
             endif
          enddo
       endif
    enddo
    ! Then compute boundaries for imm,jmm1,jmm2
    pt7=0
    p7=0
    if(minp2==0)print*,'Problem defining imm'
    do i=minp2,mp-gap
      marine=.false.
      if(i==ip(i))then
        do j=gap+1,np-gap
          if(seafloor(i,j)>0.and.seafloor(i,j)<maxdepth.and..not.marine.and.j/=np-gap)then
            pt7=pt7+1
            p7(pt7,1)=i
            p7(pt7,2)=j
            marine=.true.
          elseif(marine.and.seafloor(i,j)<=0.)then
            marine=.false.
            p7(pt7,3)=j-1
          elseif(marine.and.j==np-gap)then
            p7(pt7,3)=j
          endif
        enddo
      endif
    enddo

    ! Allocate computational grid boundaries
    nn=pt1
    if(allocated(jn)) deallocate(jn)
    if(allocated(in1)) deallocate(in1)
    if(allocated(in2)) deallocate(in2)
    allocate(in1(nn),in2(nn),jn(nn))
    in1=0
    in2=0
    jn=0
    jn(1:nn)=p1(1:nn,1)
    in1(1:nn)=p1(1:nn,2)
    in2(1:nn)=p1(1:nn,3)

    mm=pt2
    if(allocated(im)) deallocate(im)
    if(allocated(jm1)) deallocate(jm1)
    if(allocated(jm2)) deallocate(jm2)
    allocate(im(mm),jm1(mm),jm2(mm))
    im=0
    jm1=0
    jm2=0
    im(1:mm)=p2(1:mm,1)
    jm1(1:mm)=p2(1:mm,2)
    jm2(1:mm)=p2(1:mm,3)

    nnn=pt3
    if(allocated(ina)) deallocate(ina)
    if(allocated(inb)) deallocate(inb)
    if(allocated(jnab)) deallocate(jnab)
    allocate(ina(nnn),inb(nnn),jnab(nnn))
    ina=0
    inb=0
    jnab=0
    jnab(1:nnn)=p3(1:nnn,1)
    ina(1:nnn)=p3(1:nnn,2)
    inb(1:nnn)=p3(1:nnn,3)

    mmm=pt4
    if(allocated(jmb)) deallocate(jmb)
    if(allocated(jma)) deallocate(jma)
    if(allocated(imab)) deallocate(imab)
    allocate(imab(mmm),jma(mmm),jmb(mmm))
    imab=0
    jma=0
    jmb=0
    imab(1:mmm)=p4(1:mmm,1)
    jma(1:mmm)=p4(1:mmm,2)
    jmb(1:mmm)=p4(1:mmm,3)

    ne=pt5
    if(allocated(ie)) deallocate(ie)
    if(allocated(je)) deallocate(je)
    if(allocated(ije)) deallocate(ije)
    allocate(ie(ne),je(ne),ije(ne))
    ie=0
    je=0
    ije=0
    ie(1:ne)=p5(1:ne,1)
    je(1:ne)=p5(1:ne,2)
    ije(1:ne)=p5(1:ne,3)

    nnnn=pt6
    if(allocated(jnn)) deallocate(jnn)
    if(allocated(inn1)) deallocate(inn1)
    if(allocated(inn2)) deallocate(inn2)
    allocate(inn1(nnnn),inn2(nnnn),jnn(nnnn))
    inn1=0
    inn2=0
    jnn=0
    jnn(1:nnnn)=p6(1:nnnn,1)
    inn1(1:nnnn)=p6(1:nnnn,2)
    inn2(1:nnnn)=p6(1:nnnn,3)

    mmmm=pt7
    if(allocated(imm)) deallocate(imm)
    if(allocated(jmm1)) deallocate(jmm1)
    if(allocated(jmm2)) deallocate(jmm2)
    allocate(imm(mmmm),jmm1(mmmm), jmm2(mmmm))
    imm=0
    jmm1=0
    jmm2=0
    imm(1:mmmm)=p7(1:mmmm,1)
    jmm1(1:mmmm)=p7(1:mmmm,2)
    jmm2(1:mmmm)=p7(1:mmmm,3)

    return

  end subroutine build_bounds
  ! =====================================================================================
  subroutine build_output

    integer::k,i,i1,i2,j1,j

    real::Vma,Uma,x1,y1

    Umean=0.
    Vmean=0.
    if(iam==0)then
      do i=1,mp
        do j=1,np
           aa(i,j)= 0.0
        enddo
      enddo
      do k=1,nn
        i1=in1(k)
        i2=in2(k)
        j=jn(k)
        if(i1==ip(i1)) i1=i1+1
        if(i2==ip(i2)) i2=i2-1
        do i=i1,i2,2
           if(i>=1.and.i<=mp) aa(i,j)=a(i,j)
        enddo
      enddo
      ! Get velocities and directions
      do k=1,nn
        i1=in1(k)
        i2=in2(k)
        j=jn(k)
        if(i1/=ip(i1)) i1=i1-1
        if(i2/=ip(i2)) i2=i2+1
        do i=i1,i2,2
          Uma=0
          Vma=0
          if(i>=1.and.i<=mp)then
            Vmean(i,j)=a(i,j)
          endif
          if(i>1.and.i<mp)then
            Umean(i,j)=0.25*(a(i+1,j+1)+a(i+1,j-1)+a(i-1,j-1)+a(i-1,j+1))
          endif
        enddo
      enddo
      ! call mpi_allreduce(mpi_in_place,Umean,m*n,real_type,max_type,ocean_comm_world,ierr)
      ! call mpi_allreduce(mpi_in_place,Vmean,m*n,real_type,max_type,ocean_comm_world,ierr)
      ! Sort computed values points
      i1=0
      j1=0
      do i=1,m
        if(i==ip(i)) i1=i1+1
        do j=1,n
          if(i==ip(i).and.j==ip(j))then
            j1=j1+1
            o1(j1,i1)=Umean(i,j)
            o2(j1,i1)=Vmean(i,j)
          endif
        enddo
        j1=0
      enddo
      ! Interpolate velocities over the computational grid
      do i=1,m
        do j=1,n
          if((j-1)*dx/100>x1_array(1).and.(j-1)*dx/100<x1_array(xlen))then
            if((i-1)*dx/100>y1_array(1).and.(i-1)*dx/100<y1_array(ylen))then
              y1=(i-1)*dx/100
              x1=(j-1)*dx/100
              Umean(i,j)=interpolate(xlen,x1_array(1:xlen),ylen,y1_array(1:ylen),o1,x1,y1)
              Vmean(i,j)=interpolate(xlen,x1_array(1:xlen),ylen,y1_array(1:ylen),o2,x1,y1)
            endif
          endif
        enddo
      enddo
      ! Interpolate velocities over stratigraphic mesh
      Uave=0.
      Vave=0.
      do i=1,sp_m
        do j=1,sp_n
          if((j-1)*stratal_dx/100>=xc_array(1).and.(j-1)*stratal_dx/100<=xc_array(xclen))then
            if((i-1)*stratal_dx/100>=yc_array(1).and.(i-1)*stratal_dx/100<=yc_array(yclen))then
              if(sp_topo(i,j)-sea_level<0.0)then
                Uave(i,j)=Umean(i+gap,j+gap)
                Vave(i,j)=Vmean(i+gap,j+gap)
              elseif(sp_topo(i,j)<-5900.0)then
                Uave(i,j)=0.0
                Vave(i,j)=0.0
              endif
            endif
          endif
        enddo
      enddo
      do i=1,mp
        do j=1,np
          if(i==ip(i).and.j/=ip(j)) goto 500
            aa(i,j)=a(i,j)
500         continue
        enddo
      enddo
    endif

    return

  end subroutine build_output
  ! =====================================================================================
  subroutine current_final

    ! Circulation computational grid arrays
    if(allocated(inn1))then
       deallocate(inn1,inn2,jnn,imm,jmm1,jmm2)
       deallocate(ina,inb,jnab,imab,jma,jmb)
       deallocate(in1,in2,jn,im,jm1,jm2)
       deallocate(ie,je,ije)
       deallocate(a,aa,p)
    endif

    ! Circulation boundaries
    if(allocated(p1))then
       deallocate(p1,p2,p3,p4,p5,p6,p7)
    endif

    ! Deallocate dataset
    if(allocated(o1)) deallocate(o1,o2)
    if(allocated(Uave)) deallocate(Uave,Vave)
    if(allocated(Umean)) deallocate(Umean,Vmean)
    if(allocated(lUmean)) deallocate(lUmean,lVmean)
    if(allocated(pUmean)) deallocate(pUmean,pVmean)
    if(allocated(seafloor)) deallocate(seafloor,sp_topo)

    return

  end subroutine current_final
  ! =====================================================================================

end module currents_io
