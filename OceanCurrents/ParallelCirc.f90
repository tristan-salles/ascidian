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
!    Description:  Encapsulates the parallelisation of the circulation module.
!
!        Version:  1.0
!        Created:  25/09/15 10:24:45
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module parallel_circ

  use mpidata
  use currents_data
  use precision_data

  implicit none

  ! Circulation grid computational parameters
  integer::pt1,pt2,pt3,pt4,pt5,pt6,pt7
  integer,dimension(:,:),allocatable::p1,p2,p3,p4,p5,p6,p7

  public

contains

  ! ============================================================================
  !> Subroutine define_band_partitioning
  !! Subroutine define_band_partitioning used to send the partitioned grid values.
  !<
  ! ============================================================================
  subroutine define_band_partitioning

    integer::avecol,extras,cols,k,halo,istart
    integer,dimension(nprocs)::circ_col

    if(allocated(cextent)) deallocate(cextent)
    allocate(cextent(nprocs,2))

    if(iam==0)then
      ! Do not use any partitioning for the circulation model
      ! allocate(cextent(1,2))
      mp=m
      np=n
      cextent(1,1)=1
      cextent(1,2)=m
      if(allocated(psea)) deallocate(psea)
      allocate(psea(mp,np))
      return
    else
      mp=0
      np=0
      cextent(iam+1,1:2)=0
      return
    endif

    ! Partitioning technique for parallel run based on row partitioning
    avecol=int(m/nprocs)
    extras=mod(m,nprocs)

    ! if(allocated(cextent)) deallocate(cextent)
    ! allocate(cextent(nprocs,2))
    phalor=2
    phalol=2
    halo=2

    ! Get the number of rows for each processor
    if(iam==0)then
      do k=0,nprocs-1
        cols=avecol
        if(k<extras) cols=avecol+1
        if(k==0.or.k==nprocs-1)then
          if(k==0)then
            cextent(k+1,1)=1
            istart=(cols+phalor)-phalol-phalor
            ! In case of even number
            if(ip(istart)==istart)then
              circ_col(k+1)=cols+phalor+1
            else
              circ_col(k+1)=cols+phalor
            endif
            cextent(k+1,2)=circ_col(k+1)
          else
            cextent(k+1,1)=cextent(k,2)-phalol-phalor
            cextent(k+1,2)=m
            circ_col(k+1)=cextent(k+1,2)-cextent(k+1,1)+1
          endif
        else
          cextent(k+1,1)=cextent(k,2)-phalol-phalor
          istart=cextent(k+1,1)+cols+phalol+phalor-1-phalol-phalor
          ! In case of even number
          if(ip(istart)==istart)then
            circ_col(k+1)=cols+phalol+phalor+1
          else
            circ_col(k+1)=cols+phalol+phalor
          endif
          cextent(k+1,2)=cextent(k+1,1)+circ_col(k+1)-1
        endif
      enddo
    endif

    call mpi_bcast(cextent(1:nprocs,1),nprocs,int_type,0,ocean_comm_world,ierr)
    call mpi_bcast(cextent(1:nprocs,2),nprocs,int_type,0,ocean_comm_world,ierr)

    mp=cextent(iam+1,2)-cextent(iam+1,1)+1
    np=n
    phalol=halo

    if(allocated(psea)) deallocate(psea)
    allocate(psea(mp,np))

    return

  end subroutine define_band_partitioning
  ! ============================================================================
  !> Subroutine parallel_halo
  !! Subroutine parallel_halo send/rcv halo borders to boundaries.
  !<
  ! ============================================================================
  subroutine parallel_halo(tp)

    integer::k,tag,req,req2,k1,k2,halo,tp
    integer,dimension(mpi_status_size)::stat1,stat2
    real,dimension(np)::s,r

    halo=2
    ! Pass a matrix halo values
    if(tp==1.or.tp==3)then
      if(iam+1<nprocs)then
        tag=0
        req=0
        k1=mp-phalor-halo
        k2=mp-phalor-1
        do k=k1,k2
          s(1:np)=a(k,1:np)
          call mpi_isend(s,np,mpi_real,iam+1,121+tag,ocean_comm_world,req,ierr)
          call mpi_request_free(req,ierr)
          tag=tag+1
          req=req+1
        enddo
        tag=0
        req2=300
        k1=mp-halo+1
        k2=mp
        do k=k1,k2
          call mpi_irecv(r,np,real_type,iam+1,131+tag,ocean_comm_world,req2,ierr)
          call mpi_wait(req2,stat2,ierr)
          a(k,1:np)=r(1:np)
          tag=tag+1
          req2=req2+1
        enddo
       endif
       if(iam>0)then
         tag=0
         req2=300
         k1=halo+2
         k2=halo+1+phalor
         do k=k1,k2
           s(1:np)=a(k,1:np)
           call mpi_isend(s,np,real_type,iam-1,131+tag,ocean_comm_world,req2,ierr)
           call mpi_request_free(req2,ierr)
           tag=tag+1
           req2=req2+1
         enddo
         tag=0
         req=0
         k1=1
         k2=halo
         do k=k1,k2
           call mpi_irecv(r,np,mpi_real,iam-1,121+tag,ocean_comm_world,req,ierr)
           call mpi_wait(req,stat1,ierr)
           a(k,1:np)=r(1:np)
           tag=tag+1
           req=req+1
         enddo
      endif
    endif

    ! Pass the aa matrix halo values
    if(tp==2.or.tp==3)then
      if(iam+1<nprocs)then
        tag=0
        req=0
        k1=mp-phalor-halo
        k2=mp-phalor-1
        do k=k1,k2
          s(1:np)=aa(k,1:np)
          call mpi_isend(s,np,mpi_real,iam+1,141+tag,ocean_comm_world,req,ierr)
          call mpi_request_free(req,ierr)
          tag=tag+1
          req=req+1
        enddo
        tag=0
        req2=300
        k1=mp-halo+1
        k2=mp
        do k=k1,k2
          call mpi_irecv(r,np,real_type,iam+1,151+tag,ocean_comm_world,req2,ierr)
          call mpi_wait(req2,stat2,ierr)
          aa(k,1:np)=r(1:np)
          tag=tag+1
          req2=req2+1
        enddo
       endif
       if(iam>0)then
        tag=0
        req2=300
        k1=halo+2
        k2=halo+1+phalor
        do k=k1,k2
          s(1:np)=aa(k,1:np)
          call mpi_isend(s,np,real_type,iam-1,151+tag,ocean_comm_world,req2,ierr)
          call mpi_request_free(req2,ierr)
          tag=tag+1
          req2=req2+1
        enddo
        tag=0
        req=0
        k1=1
        k2=halo
        do k=k1,k2
          call mpi_irecv(r,np,mpi_real,iam-1,141+tag,ocean_comm_world,req,ierr)
          aa(k,1:np)=r(1:np)
          call mpi_wait(req,stat1,ierr)
          tag=tag+1
          req=req+1
        enddo
      endif
    endif

    return

  end subroutine parallel_halo
  ! ============================================================================

end module parallel_circ
