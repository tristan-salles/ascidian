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
!       Filename:  SimWaves.f90
!
!    Description:  Computes wave field.
!
!        Version:  1.0
!        Created:  24/09/15 07:42:49
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module SimWaves

    use mpidata
    use currents_io
    use currents_data
    use precision_data
    use WaveCirculation
    use swan_coupler_export
    use swan_coupler_parallel
    use swan_coupler_functions

    implicit none

    real,dimension(6)::hcast

    public

contains

    ! ============================================================================
    subroutine waves_init

      integer::xpart,ypart

      if(.not.waveon)then
        if(.not.allocated(wdcp)) allocate(wdcp(nprocs))
        wdcp(iam+1)%minX_id=1
        wdcp(iam+1)%maxX_id=sp_n
        wdcp(iam+1)%minY_id=1
        wdcp(iam+1)%maxY_id=sp_m
        wdcp(iam+1)%X_nb=sp_n
        wdcp(iam+1)%Y_nb=sp_m
        if(.not.allocated(exportSWANfields)) allocate(exportSWANfields(sp_n,sp_m,2))
        exportSWANfields=0.0
        if(.not.allocated(bathyfield)) allocate(bathyfield(sp_n,sp_m))
        if(.not.allocated(glob_Uw)) allocate(glob_Uw(sp_n*sp_m))
        if(.not.allocated(glob_Dw)) allocate(glob_Dw(sp_n*sp_m))
        glob_Uw=0.0
        glob_Dw=0.0
        bathyfield=0.0
        if(allocated(l_Uw)) deallocate(l_Uw)
        allocate(l_Uw(sp_n,sp_m))
        if(allocated(g_Uw)) deallocate(g_Uw)
        allocate(g_Uw(sp_n*sp_m))
        if(allocated(uCur)) deallocate(uCur)
        if(allocated(dCur)) deallocate(dCur)
        allocate(uCur(sp_n,sp_m),dCur(sp_n,sp_m))
        return
      endif
      call swan_createinput
      call swan_initialize(ocean_comm_world,swaninput(1:80),swaninfo(1:80))
      if(allocated(wdcp)) deallocate(wdcp)
      allocate(wdcp(nprocs))
      call swan_decomposition_grid(nprocs,wdcp)
      xpart=wdcp(iam+1)%X_nb
      ypart=wdcp(iam+1)%Y_nb
      if(allocated(uCur)) deallocate(uCur)
      if(allocated(dCur)) deallocate(dCur)
      allocate(uCur(xpart,ypart),dCur(xpart,ypart))
      if(allocated(l_Uw)) deallocate(l_Uw)
      if(allocated(l_Dw)) deallocate(l_Dw)
      allocate(l_Uw(sp_n,sp_m),l_Dw(sp_n,sp_m))
      if(allocated(g_Uw)) deallocate(g_Uw)
      if(allocated(g_Dw)) deallocate(g_Dw)
      allocate(g_Uw(sp_n*sp_m),g_Dw(sp_n*sp_m))
      if(allocated(glob_Dw)) deallocate(glob_Dw)
      if(allocated(glob_Uw)) deallocate(glob_Uw)
      allocate(glob_Uw(sp_n*sp_m),glob_Dw(sp_n*sp_m))

      return

    end subroutine waves_init
    ! ============================================================================
    subroutine swan_createinput

      integer::ios,iu
      real::mx_x,mx_y
      character(len=128)::stg1,stg2

      ! Create the input for SWAN simulation
      swaninput='swan.swn'
      call noblnk(swaninput)
      call addpath1(swaninput)
      swaninfo='swanInfo'
      call noblnk(swaninfo)
      call addpath1(swaninfo)
      swanbot='swan.bot'
      call noblnk(swanbot)
      call addpath1(swanbot)
      swanout='swan.csv'
      call noblnk(swanout)
      call addpath2(swanout)
      if(iam==0)then
        iu=342
        mx_x=stratal_dx*(stratal_x-1)/100
        mx_y=stratal_dx*(stratal_y-1)/100
        open(iu,file=swaninput,status="replace",action="write",iostat=ios)
        write(iu,'(a17)') "PROJECT ' ' 'S01' "
        write(iu,102) "CGRID REG",stratal_xo,stratal_yo,0,mx_x,mx_y,stratal_x-1,stratal_y-1,"CIRCLE 36 0.05 1.0"
102     format(a9,1x,f12.3,1x,f12.3,1x,i1,1x,f12.3,1x,f12.3,1x,i4,1x,i4,1x,a19)
        write(iu,103) "INPGRID BOTTOM REG",stratal_xo,stratal_yo,0,1,1,mx_x,mx_y,'EXC -999999.000'
103     format(a18,1x,f12.3,1x,f12.3,1x,i1,1x,i1,1x,i1,1x,f12.3,1x,f12.3,1x,a15)
        stg1="READINP BOTTOM 1 '"
        call append_str2(stg1,swanbot)
        stg2="' 3 0 FREE"
        call append_str2(stg1,stg2)
        write(iu,*)trim(stg1)
        write(iu,'(a42)') "BOUNd SHAPESPEC JONSWAP 3.3 PEAK DSPR DEGR"
        write(iu,'(a7)') 'DIFFRAC'
        write(iu,'(a8)') 'FRICTION'
        !write(iu,'(a26)') 'BREAKING CONSTANT 1.1 0.73'
        !write(iu,'(a8)') 'WCAPPING'
        !write(iu,'(a8)') 'OFF QUAD'
        write(iu,'(a4)') 'GEN1'
        !write(iu,'(a4)') 'QUAD'
        !write(iu,'(a10)') 'GEN3 AGROW'
        !write(iu,'(a10)') 'OFF BNDCHK'
        write(iu,'(a15,1x,i1,1x,i4,1x,i1,1x,i4)')"GROUP 'gf' SUBG",0,stratal_x-1,0,stratal_y-1
        stg1="TABLE 'gf' IND '"
        call append_str2(stg1,swanout)
        stg2="' XP YP DIR UBOT" !HS PER WLEN UBOT"
        call append_str2(stg1,stg2)
        write(iu,*)trim(stg1)
        write(iu,'(a5,2f12.3)') 'WIND ',forecast_param(5:6)
!       write(iu,'(a9,4f12.3)')'INIT PAR ',forecast_param(1:4)
        write(iu,'(a7)') 'COMPUTE'
        close(iu)
        open(iu,file=swanbot,status="replace",action="write",iostat=ios)
        write(iu,'(a6)')'-1 -1'
        write(iu,'(a6)')'-1 -1'
        close(iu)
      endif
      call mpi_barrier(ocean_comm_world,ierr)
      return

    end subroutine swan_createinput
    ! ============================================================================
    subroutine waves_circulation(sgp)

      integer::i,j,ks,ii,jj,sgp

      ! Get the bathymetry from sp model
      ii=0
      ks=iam+1
      bathyfield=0.0
      do i=wdcp(ks)%minX_id,wdcp(ks)%maxX_id
        ii=ii+1
        jj=0
        do j=wdcp(ks)%minY_id,wdcp(ks)%maxY_id
          jj=jj+1
          bathyfield(ii,jj)=sea_level-sp_topo(j,i)
          uCur(ii,jj)=Uave(j,i) !/100.
          dCur(ii,jj)=Vave(j,i)
          if(bathyfield(ii,jj)<0.) bathyfield(ii,jj)=-999999.0
          if(bathyfield(ii,jj)>wave_base) bathyfield(ii,jj)=-999999.0
        enddo
      enddo
      if(.not.waveon) goto 15
      ! Define forcing waves parameters
      call import_bathymetry
      hcast(1:6)=forecast_param(1:6)
      call swan_run(hcast)
15    continue
      call combined_wavecurrent(sgp)

      return

    end subroutine waves_circulation
    ! ============================================================================
    subroutine wave_final

      if(waveon) call swan_finalize
      if(allocated(l_Uw)) deallocate(l_Uw)
      if(allocated(g_Uw)) deallocate(g_Uw)
      if(allocated(glob_Uw)) deallocate(glob_Uw)
      if(allocated(uCur)) deallocate(uCur)
      if(allocated(dCur)) deallocate(dCur)
      if(allocated(l_Dw)) deallocate(l_Dw)
      if(allocated(g_Dw)) deallocate(g_Dw)
      if(allocated(glob_Dw)) deallocate(glob_Dw)
      if(allocated(exportSWANfields)) deallocate(exportSWANfields)

      return

    end subroutine wave_final
    ! ============================================================================

end module SimWaves
