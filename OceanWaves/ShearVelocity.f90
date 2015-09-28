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
!       Filename:  WaveCirculation.f90
!
!    Description:  Computes wave-current interaction based on Grant & Madsen,1986
!                   and Li & Amos,1995.
!
!        Version:  1.0
!        Created:  17/09/15 13:49:32
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module WaveCirculation

  use hdf5
  use mpidata
  use currents_io
  use currents_data
  use precision_data
  use swan_coupler_export
  use swan_coupler_parallel
  use swan_coupler_functions

  implicit none

  real::hs,period
  real,parameter::kappa=0.4
  real,dimension(:,:),allocatable::uCur,dCur
  real,dimension(:,:),allocatable::l_Uw,l_Dw
  real,dimension(:),allocatable::g_Uw,g_Dw
  real,dimension(:),allocatable::glob_Uw,glob_Dw
  !  Decomposition specification for SWAN Model
  type(SGFM_WADecomp),dimension(:),allocatable::wdcp

  public

contains

  ! ============================================================================
  subroutine combined_wavecurrent(sgp)

    integer::i,j,p,n,xpart,ypart,ii,jj,xo,yo,x1,y1,sgp
    real::depth

    xo=wdcp(iam+1)%minX_id
    yo=wdcp(iam+1)%minY_id
    x1=wdcp(iam+1)%maxX_id
    y1=wdcp(iam+1)%maxY_id
    xpart=wdcp(iam+1)%X_nb
    ypart=wdcp(iam+1)%Y_nb
    ! Combined contour currents and ocean circulation
    if(circon)then
      do i=1,sp_n
        do j=1,sp_m
          p=p+1
          if(sp_topo(j,i)<-5900.0)then
            Uave(j,i)=0.0
            Vave(j,i)=0.0
           endif
        enddo
      enddo
    endif
    do i=1,xpart
      do j=1,ypart
        depth=bathyfield(i,j)
        ! Wave orbital velocity
        if(depth<0) exportSWANfields(i,j,1)=0.0
        if(depth<0) exportSWANfields(i,j,2)=0.0
       enddo
    enddo
    ! Gather values to all processors
    if(waveon)then
      l_Uw=-1.e8
      l_Dw=-1000.0
    endif
    ii=0
    do i=xo,x1
      ii=ii+1
      jj=0
      do j=yo,y1
        jj=jj+1
        if(waveon)then
          l_Uw(i,j)=exportSWANfields(ii,jj,1)
          l_Dw(i,j)=exportSWANfields(ii,jj,2)
        endif
      enddo
    enddo
    p=0
    do i=1,sp_n
      do j=1,sp_m
        p=p+1
        if(waveon)then
          g_Uw(p)=l_Uw(i,j)
          g_Dw(p)=l_Dw(i,j)
        endif
      enddo
    enddo
    if(waveon)then
      call mpi_allreduce(g_Uw,glob_Uw,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)
      call mpi_allreduce(g_Dw,glob_Dw,sp_n*sp_m,real_type,max_type,ocean_comm_world,ierr)
    else
      glob_Uw=0.0
      glob_Dw=0.0
    endif

    ! Write the output
    if(iam==0.and.outputflag)then
      n=1
      do i=1,sp_m
        p=i
        do j=1,sp_n
          curU(sgp,n)=real(Uave(i,j)/100.)
          !curU(sgp,n)=min(5.0,curU(sgp,n))
          curV(sgp,n)=real(Vave(i,j)/100.)
          !curV(sgp,n)=min(5.0,curV(sgp,n))
          if(glob_Uw(p)>0.0)then
            wavU(sgp,n)=real(glob_Uw(p)*cos(glob_Dw(p)*pi/180.))
            wavV(sgp,n)=real(glob_Uw(p)*sin(glob_Dw(p)*pi/180.))
          else
            wavU(sgp,n)=0.0
            wavV(sgp,n)=0.0
          endif
          p=p+sp_m
          n=n+1
        enddo
      enddo
    endif

    if(iam==0.and.spmon)then
      ! Write H5 output for lecode
      call write_h5_spm(sgp)
    endif

    return

  end subroutine combined_wavecurrent
  ! ============================================================================
  subroutine write_h5_spm(sgp)

    logical::compression
    integer::sgp,hdferr,i,j,p,n,rank
    character(len=128)::stg,text
    real(tkind),dimension(2*sp_n*sp_m)::Vparamc
    real(tkind),dimension(2*sp_n*sp_m)::Vparamw
    real(tkind),dimension(4*sp_n*sp_m)::Vparam
    integer(hid_t)::file_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    velofile='ocean_scn'
    call noblnk(velofile)
    call noblnk(velofile)
    call append_nb(velofile,sgp)
    stg='.h5'
    call append_str(velofile,stg)
    call noblnk(velofile)
    call addpath3(velofile)

    ! Initialize predefined datatypes
    call h5open_f(hdferr)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,hdferr)

    ! Create the file collectively.
    call h5fcreate_f(velofile,h5f_acc_trunc_f,file_id,hdferr)

    dims(1)=4
    if(.not.waveon.and.circon) dims(1)=2
    if(waveon.and..not.circon) dims(1)=2
    dims(2)=sp_n*sp_m
    rank=2

    ! Create dataspace and opens it for access
    call h5screate_simple_f(rank,dims,filespace,hdferr)
    p=0
    do i=1,sp_m
      n=i
      do j=1,sp_n
        if(.not.waveon.and.circon)then
          if(Uave(i,j)==0..and.Vave(i,j)==0.)then
            Vparamc(p+1)=0.
            Vparamc(p+2)=0.
          else
            Vparamc(p+1)=sqrt((Uave(i,j)/100.)**2.0+(Vave(i,j)/100.)**2.0)
            Vparamc(p+1)=min(2.5,Vparamc(p+1))
            Vparamc(p+2)=atan2(Vave(i,j),Uave(i,j))
          endif
          p=p+2
        elseif(waveon.and..not.circon)then
          Vparamw(p+1)=glob_Uw(n)
          Vparamw(p+2)=glob_Dw(n)
          p=p+2
        else
          if(Uave(i,j)==0..and.Vave(i,j)==0.)then
            Vparam(p+1)=0.
            Vparam(p+2)=0.
          else
            Vparam(p+1)=sqrt((Uave(i,j)/100.)**2.0+(Vave(i,j)/100.)**2.0)
            Vparam(p+1)=min(2.5,Vparam(p+1))
            Vparam(p+2)=atan2(Vave(i,j),Uave(i,j))
          endif
          Vparam(p+3)=glob_Uw(n)
          Vparam(p+4)=glob_Dw(n)
          p=p+4
        endif
        n=n+sp_m
       enddo
    enddo
    text=''
    text="/WavesCirculation"

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,hdferr)
    ! Write the dataset collectively
    if(.not.waveon.and.circon) then
      call h5dwrite_f(dset_id,h5t_native_double,Vparamc,dims,hdferr)
    elseif(waveon.and..not.circon)then
      call h5dwrite_f(dset_id,h5t_native_double,Vparamw,dims,hdferr)
    else
      call h5dwrite_f(dset_id,h5t_native_double,Vparam,dims,hdferr)
    endif

    ! Close the dataset
    call h5dclose_f(dset_id,hdferr)
    call h5sclose_f(filespace,hdferr)

    ! Close the file.
    call h5fclose_f(file_id,hdferr)
    ! Close interface
    call h5close_f(hdferr)

    return

  end subroutine write_h5_spm
  ! ============================================================================

end module WaveCirculation
