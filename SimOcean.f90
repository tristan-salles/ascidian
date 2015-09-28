! =====================================================================================
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
!       Filename:  SimOcean.f90
!
!    Description:  Top level Ocean Application Driver
!
!        Version:  1.0
!        Created:  21/09/15 08:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================

program SimOcean

  use mpidata
  use SimWaves
  use currents_circ
  use currents_data
  use precision_data
  use OceanOut

  implicit none

  integer::arg,scn,sgp,fst,scn1,scn2,id,ido
  real(tkind)::time_st,time_ed

  ! start up MPI
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  call mpi_comm_rank(mpi_comm_world,iam,ierr)
  ocean_comm_world=mpi_comm_world

  int_type=mpi_integer
  real_type=mpi_real
  dbl_type=mpi_double_precision
  lgc_type=mpi_logical
  max_type=mpi_max
  min_type=mpi_min
  sum_type=mpi_sum

  time_st=mpi_wtime( )
  ! Get experiment file name
  if(iam==0)then
    arg=iargc()
    finput=' '
    if(arg<1)then
      write(6,*)'use: mpirun -np X ./ocean <input-file-name>'
      write(6,*)'where X refers to the number of processors for the run.'
      stop
    elseif(arg==1)then
      call getarg(1,finput)
    endif
  endif
  call mpi_bcast(finput,128,mpi_character,0,ocean_comm_world,ierr)

  ! Initialise currents and waves
  call currents_init
  call waves_init
  fst=0
  id=0
  ido=1
  ! Time loop
  tnow=time_start
  toutput=time_start
  do while(tnow<time_end)
    ! Wait for SPM model
    if(spmon) call wait_function
    ! Find the appropriate climatic forces
    do scn=1,forecast_nb
      if(tnow>=hindcast(scn)%tstart.and.tnow<hindcast(scn)%tend)then
        scn1=scn
        scn2=scn
        if(tnow+time_step>=hindcast(scn)%tend)then
          scn2=scn+1
          if(scn2 >forecast_nb) scn2=forecast_nb
        endif
      endif
    enddo
    id=id+1
    do sgp=1,hindcast(scn1)%cnb
      ! Allocate forecasts
      ! For the current circulation use current step
      forecast_param(7)=hindcast(scn1)%subgroup(sgp)%wvel*sin(hindcast(scn1)%subgroup(sgp)%wdir*pi/180.)
      forecast_param(8)=hindcast(scn1)%subgroup(sgp)%wvel*cos(hindcast(scn1)%subgroup(sgp)%wdir*pi/180.)
      ! For the wave simulation initialise the next step
      if(sgp<hindcast(scn1)%cnb)then
        forecast_param(1)=hindcast(scn1)%subgroup(sgp+1)%hs
        forecast_param(2)=hindcast(scn1)%subgroup(sgp+1)%per
        forecast_param(3)=hindcast(scn1)%subgroup(sgp+1)%dir
        forecast_param(4)=hindcast(scn1)%subgroup(sgp+1)%dd
        forecast_param(5)=hindcast(scn1)%subgroup(sgp+1)%wvel
        forecast_param(6)=hindcast(scn1)%subgroup(sgp+1)%wdir
      elseif(scn1==scn2)then
        forecast_param(1)=hindcast(scn1)%subgroup(1)%hs
        forecast_param(2)=hindcast(scn1)%subgroup(1)%per
        forecast_param(3)=hindcast(scn1)%subgroup(1)%dir
        forecast_param(4)=hindcast(scn1)%subgroup(1)%dd
        forecast_param(5)=hindcast(scn1)%subgroup(1)%wvel
        forecast_param(6)=hindcast(scn1)%subgroup(1)%wdir
      else
        forecast_param(1)=hindcast(scn2)%subgroup(1)%hs
        forecast_param(2)=hindcast(scn2)%subgroup(1)%per
        forecast_param(3)=hindcast(scn2)%subgroup(1)%dir
        forecast_param(4)=hindcast(scn2)%subgroup(1)%dd
        forecast_param(5)=hindcast(scn2)%subgroup(1)%wvel
        forecast_param(6)=hindcast(scn2)%subgroup(1)%wdir
      endif
      ! First compute current circulation
      call currents_circulation(fst)
      fst=1
      call mpi_barrier(ocean_comm_world,ierr)
      if(iam==0.and.circon)print*,' '
      ! Then compute the wave fields
      call waves_circulation(sgp)
      if(iam==0.and.waveon) print*,' '
    enddo

    ! Output ocean
    if( outputflag.and.abs(tnow-toutput)<=0.001)then
      call xdmf_output(ido,hindcast(scn1)%cnb)
      toutput=toutput+out_step
      ido=ido+1
    endif
    if(iam==0)then
      print*,'-----------------------------------------------'
      write(*,101)'Current time : ',tnow
      print*,'-----------------------------------------------'
      print*,' '
    endif
101  format(a15,f16.3)
    tnow=tnow+time_step
    ! Exchange with SPM model
    if(spmon.and.iam==0) call update_maestro
  enddo

  call mpi_barrier(ocean_comm_world,ierr)
  time_ed=mpi_wtime( )
  if(iam==0) print*,'Time elapse:',time_ed-time_st

  ! Finalisation
  call current_final
  call wave_final

  return

end program SimOcean
