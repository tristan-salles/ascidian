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
!       Filename:  CurrentClass.f90
!
!    Description:  Encapsulates main ocean classes.
!
!        Version:  1.0
!        Created:  21/09/15 12:17:19
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module currents_data

  ! SP grid model values
  integer::sp_n,sp_m,nbPt,faces_nb
  integer::xlen,ylen,xclen,yclen,stratal_x,stratal_y
  integer,parameter::gap=2

  logical::outputflag

  real::stratal_xo,stratal_yo,dx,st_dx,stratal_dx,stratal_xm,stratal_ym,initDep
  real,dimension(:,:),allocatable::grd_bathy,sea_floor,f_XYcoord

  ! Extension values for circulation grid computational values
  integer::nn,mm,n,m,nnn,mmm,ne,nnnn,mmmm

  ! Circulation grid computational arrays
  integer,dimension(:),allocatable::ina,inb,jnab,imab,jma,jmb
  integer,dimension(:),allocatable::in1,in2,jn,im,jm1,jm2
  integer,dimension(:),allocatable::inn1,inn2,jnn,imm,jmm1,jmm2
  integer,dimension(:),allocatable::ie,je,ije

  ! Circulation computation parameters
  integer::kend,ifilt,waveflag,circflag,contflag,spmflag

  real::dt,grav,rho,ah,cfric,atf
  real::alfa,beta,gamma,Courant
  real::alat,f,pi,niner,alin,maxdepth

  ! Conversion from m/s to knots
  real,parameter::knots=1.9438445

  ! Circulation computational arrays
  real,dimension(:,:),allocatable::a,aa,p

  ! Time parameters
  real::time_start,time_end,time_step,out_step

  ! Contour parameters
  real::vmean_cont,contour1,contour2
  integer::inside,outside

  ! Forecast parameters
  integer::forecast_nb
  real::wave_base
  real,dimension(8)::forecast_param

  ! Number of hindcast scenarios
  type hindcast_param
    ! Percentage of wind/wave subgroup.
    real::oc
    ! Significant wave height (in metres).
    real::hs
    ! Wave period of the energy spectrum
    real::per
    ! Peak wave direction.
    real::dir
    ! Coefficient of directional spreading.
    real::dd
    ! Wind velocity at 10 m elevation (m/s).
    real::wvel
    ! Wind direction.
    real::wdir
  end type hindcast_param

  type hindcast_def
    ! Class number
    integer::cnb
    ! Time start
    real::tstart
    ! Time end
    real::tend
    ! Subgroud parameters
    type(hindcast_param),dimension(:),allocatable::subgroup
  end type hindcast_def
  type(hindcast_def),dimension(:),allocatable::hindcast

  ! Hindcast group number
  integer,dimension(:),allocatable::hcast_gp

  ! Parallel parameters
  integer::phalol,phalor,mp,np
  integer,dimension(:,:),allocatable::cextent

  real,dimension(:,:),allocatable::seafloor,sp_topo,psea

  ! Output parameter arrays
  integer::maxsubgroup,totnodes
  real,dimension(:),allocatable::nodes,htopo
  real,dimension(:,:),allocatable::curU,curV,wavU,wavV

  public

contains

  ! ============================================================================
  integer function ip(L)
    integer::L
    ip=int(L/2)*2
    return
  end function ip
  ! =====================================================================================
  subroutine write_initialparam

    integer::i
    print*,nn,mm,n,m,ne
    print*,nnn,mmm,nnnn,mmmm
    open(44,file='meteorol.csv')
    rewind(44)
    do i=m,1,-1
      write(44,*)p(i,1:n)
    enddo
    close(44)
    open(44,file='boundary.csv')
    rewind(44)
    write(44,*)in1(1:nn)
    write(44,*)in2(1:nn)
    write(44,*)jn(1:nn)
    write(44,*)im(1:mm)
    write(44,*)jm1(1:mm)
    write(44,*)jm2(1:mm)
    write(44,*)ina(1:nnn)
    write(44,*)inb(1:nnn)
    write(44,*)jnab(1:nnn)
    write(44,*)imab(1:mmm)
    write(44,*)jma(1:mmm)
    write(44,*)jmb(1:mmm)
    write(44,*)ie(1:ne)
    write(44,*)je(1:ne)
    write(44,*)ije(1:ne)
    write(44,*)inn1(1:nnnn)
    write(44,*)inn2(1:nnnn)
    write(44,*)jnn(1:nnnn)
    write(44,*)imm(1:mmmm)
    write(44,*)jmm1(1:mmmm)
    write(44,*)jmm2(1:mmmm)
    close(44)
    return

  end subroutine write_initialparam
  ! ============================================================================

end module currents_data
! ============================================================================
module mpidata

#include "mpif.h"

  ! Error code
  integer::ierr
  ! Processor ID
  integer::iam
  ! Number of processors
  integer::nprocs

end module mpidata
! ============================================================================
module precision_data

  use mpidata

  ! Input / output files
  character(len=128)::finput
  character(len=128)::xyzfile
  character(len=128)::swaninput
  character(len=128)::swaninfo
  character(len=128)::swanbot
  character(len=128)::swanout
  character(len=128)::outdir
  character(len=128)::outdir1
  character(len=128)::outdir2
  character(len=128)::syncfolder
  character(len=128)::maestro
  character(len=128)::velofile

  ! INTEL COMPILER
  ! use ifport

  ! Single or double precision
  integer,parameter::singlep=kind(0.0)
  integer,parameter::doublep=kind(0.0d0)

  integer,parameter::tkind=doublep
  integer,parameter::sizesp=4
  integer,parameter::sizedp=8

  ! MPI integer type communicator
  integer::int_type
  ! MPI double type communicator
  integer::dbl_type
  ! MPI double type communicator
  integer::real_type
  ! MPI logical type communicator
  integer::lgc_type
  ! MPI max type communicator
  integer::max_type
  ! MPI min type communicator
  integer::min_type
  ! MPI sum type communicator
  integer::sum_type
  ! SPModel communicator
  integer::ocean_comm_world

contains

  ! ============================================================================
  subroutine term_command(cmds)

    logical(4)::result
    character(len=128)::cmds
    result=.false.
    ! INTEL FORTRAN COMPILER
    !result=systemqq(cmds)
    ! GNU FORTRAN COMPILER
    call system(cmds)

    return

  end subroutine term_command
  ! ============================================================================
  subroutine append_str2(stg1,stg2)

    integer::l1,l2
    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2

    return

  end subroutine append_str2
  ! ============================================================================
  subroutine append_str(stg1,stg2)

    integer::l1,l2
    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2
    call noblnk(stg1)

    return

  end subroutine append_str
  ! ============================================================================
  subroutine append_nb(stg1,i)

    integer::l1,l2,i
    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2
    call noblnk(stg1)

    return

  end subroutine append_nb
  ! ============================================================================
  subroutine append_zero(stg1,i)

    integer::l1,l2,i
    character(len=128)::stg1,stg2,str

    l2=len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l1=len_trim(stg2)
    str=''
    if(l1==1)then
      str(1:3)='000'
      call append_str(str,stg2)
    elseif(l1==2)then
      str(1:2)='00'
      call append_str(str,stg2)
    elseif(l1==3)then
      str(1:1)='0'
      call append_str(str,stg2)
    endif
    l1=len_trim(str)
    stg1(l2+1:l2+l1)=str
    call noblnk(stg1)

    return

  end subroutine append_zero
  ! ============================================================================
  subroutine append_nb2(stg1,i)

    integer::l1,l2,i
    character(len=128)::stg1,stg2
    l1=len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l2=len_trim(stg2)
    stg1(l1+2:l1+l2+1)=stg2

    return

  end subroutine append_nb2
  ! ============================================================================
  subroutine noblnk(string)

    integer::i,jlg
    character(len=128)::string

    lg=len(string)
    do
      if(lg<=0.or.string(lg:lg)/=' ') exit
      lg=lg-1
    enddo
    if(lg>0)then
      ! find first non-blank character
      i=1
      do
        if(i>lg.or.(string(i:i)/=' '.and.string /= ' ')) exit
        i=i+1
      enddo
      ! determine end of continuous (non-blank) string
      j=i
      do
        if(j>lg)then
          exit
        elseif(string(j:j)==' ')then
          exit
        elseif(string=='  ')then
          exit
        elseif(j==128)then
          exit
        endif
        j=j+1
      enddo
      ! j points to first blank position or past end of string; adjust to last
      ! non-blank position in string
      j=min(j-1,lg)
      string=string(i:j)
      if(j<len(string)) string(j+1:len(string))=' '
    else
       ! there were only blanks in string
       string=' '
    endif

    return

  end subroutine noblnk
  ! ============================================================================
  subroutine addpath(fname)

    integer:: pathlen,flen
    character(len=128)::fname,dname,dummy

    ! for files to be read,they'll be in the session path
    dname=' '
    call noblnk(outdir)
    pathlen=len_trim(outdir)
    dname(1:pathlen)=outdir
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath
  ! ============================================================================
  subroutine addpath1(fname)

    integer:: pathlen,flen
    character(len=128)::fname,dname,dummy

    ! for files to be read,they'll be in the session path
    dname=' '
    call noblnk(outdir1)
    pathlen=len_trim(outdir1)
    dname(1:pathlen)=outdir1
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath1
  ! ============================================================================
  subroutine addpath2(fname)

    integer::pathlen,flen
    character(len=128)::fname,dname,dummy

    ! for files to be read,they'll be in the session path
    dname=' '
    call noblnk(outdir2)
    pathlen=len_trim(outdir2)
    dname(1:pathlen)=outdir2
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath2
  ! ============================================================================
  subroutine addpath3(fname)

    integer:: pathlen,flen
    character(len=128)::fname,dname,dummy

    ! for files to be read,they'll be in the session path
    dname=' '
    call noblnk(syncfolder)
    pathlen=len_trim(syncfolder)
    dname(1:pathlen)=syncfolder
    dname(pathlen+1:pathlen+1)='/'
    pathlen=pathlen+1
    call noblnk(fname)
    flen=len_trim(fname)
    dummy=' '
    dummy=fname
    fname=' '
    fname(1:pathlen)=dname(1:pathlen)
    fname(pathlen+1:pathlen+flen)=dummy(1:flen)

    return

  end subroutine addpath3
  ! ============================================================================
  subroutine wait_function

    logical::found

    integer::iunit,ios,err

    character(len=1)::charac
    character(len=128)::stg,fildir

    if(iam == 0)then
      found=.false.
      do while(.not.found)
        call noblnk(syncfolder)
        fildir=syncfolder
        stg='/topSPM_surface.xyz'
        call noblnk(stg)
        call append_str(fildir,stg)
        call noblnk(fildir)
        inquire(file=fildir,exist=found)
        if(.not.found) call Sleep(1)
      enddo
    endif
    call mpi_barrier(ocean_comm_world,ierr)

    charac='L'

    if(iam==0)then
112   continue
      do while(charac/='O')
        iunit=121
110     continue
        inquire(file=maestro,exist=found)
        if(.not.found)then
          call Sleep(1)
          goto 110
        endif
        ! Read the maestro file
        open(iunit,file=maestro,status="old",action="read",iostat=ios)
        rewind(iunit)
        read(iunit,'(a1)',err=112,end=112,iostat=err) charac
        if(err/=0)  charac='L'
        close(iunit)
        call Sleep(1)
      enddo
    endif
    call mpi_barrier(ocean_comm_world,ierr)

    return

  end subroutine wait_function
  ! ============================================================================
  subroutine update_maestro

    integer::iunit,ios

    iunit=20

    ! Update/Create the maestro file
    open(iunit,file=maestro,status="replace",action="write",iostat=ios)
    rewind(iunit)

    write(iunit,'(a1)')'L'
    write(iunit,'(a1)')' '

    close(iunit)

    return

  end subroutine update_maestro
  ! ============================================================================

end module precision_data
