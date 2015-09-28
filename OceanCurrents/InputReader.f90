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
!       Filename:  InputReader.f90
!
!    Description:  Gather the information within the Ocean XmL input file.
!
!        Version:  1.0
!        Created:  22/09/15 09:33:51
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module xml_reader

  use mpidata
  use FoX_sax
  use mpidata
  use FoX_common
  use currents_data
  use precision_data

  implicit none


  integer::hclass,sclass

  ! Output Directory
  logical,save::in_OutDir=.false.
  ! Grid Parameters
  logical,save::toposection=.false.
  logical,save::in_Topo=.false.
  logical,save::in_Grid=.false.
  logical,save::in_GridXs=.false.
  logical,save::in_GridYs=.false.
  logical,save::in_GridSpaces=.false.
  logical,save::in_SyncFold=.false.
  logical,save::in_latGrid=.false.
  logical,save::in_spmForcing=.false.
  logical,save::in_initDep=.false.
  logical,save::in_circForcing=.false.
  logical,save::in_waveForcing=.false.
  logical,save::in_contourForcing=.false.
  character(len=128),save::grid_dxs,grid_dys
  character(len=128),save::grid_spaces,grid_initdep
  character(len=128),save::grid_contourForcing,grid_waveForcing
  character(len=128),save::grid_lat,grid_circForcing,grid_spmforcing
  ! Time parameters
  logical,save::timesection=.false.
  logical,save::in_Time=.false.
  logical,save::in_startTime=.false.
  logical,save::in_endTime=.false.
  logical,save::in_syncTime=.false.
  logical,save::in_outTime=.false.
  character(len=128),save::startTime,endTime,syncTime,outTime
  ! Circulation parameters
  logical,save::circsection=.false.
  logical,save::in_Circ=.false.
  logical,save::in_CircStep=.false.
  logical,save::in_CircTime=.false.
  logical,save::in_CircFric=.false.
  logical,save::in_CircFilter=.false.
  logical,save::in_max_depth=.false.
  character(len=128),save::circstep,circtime,circfric,circfilter,max_depth
  ! Contourite parameters
  logical,save::contsection=.false.
  logical,save::in_Cont=.false.
  logical,save::in_ContTop=.false.
  logical,save::in_ContBot=.false.
  logical,save::in_ContEnt=.false.
  logical,save::in_ContExit=.false.
  logical,save::in_ContVel=.false.
  character(len=128),save::conttop,contbot,contEnt,contExit,contVel
  ! Wave Parameters Region
  logical,save::in_wave=.false.
  logical,save::in_wbase=.false.
  logical,save::in_hindcast=.false.
  logical,save::in_wwcast=.false.
  logical,save::in_wtstart=.false.
  logical,save::in_wtend=.false.
  logical,save::in_groupnb=.false.
  logical,save::in_groupparam=.false.
  character(len=128),save::wbase,wwcast,wtstart,wtend,groupnb,groupparam

contains

  ! ============================================================================
  subroutine startDocument_handler

  end subroutine startDocument_handler
  ! ============================================================================
  subroutine startElement_handler(namespaceURI,localname,name,atts)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name
    type(dictionary_t),intent(in)::atts

    ! Output Element
    if(name=='OutputDirectory')then
      outputflag=.true.
      in_OutDir=.true.
    endif

    ! Topo Element
    if(name=='TopoGrid') in_Topo=.true.
    if(name=='TopoGrid') toposection=.true.
    if(in_Topo) call StopoElement_handler(name)

    ! Time Element
    if(name=='Time') in_Time=.true.
    if(in_Time) timesection=.true.
    if(in_Time) call StimeElement_handler(name)

    ! Circulation Element
    if(name=='CirculationParam') in_Circ=.true.
    if(in_Circ) circsection=.true.
    if(in_Circ) call ScircElement_handler(name)

    ! Contour Element
    if(name=='ContouriteParam') in_Cont=.true.
    if(in_Cont) contsection=.true.
    if(in_Cont) call ScontElement_handler(name)

    ! Wave Element
    if(name=='OceanForecast') in_wave=.true.
    if(in_wave) call SwaveElement_handler(name)
    if(name=='forecastClass') in_hindcast=.true.
    if(in_hindcast) call ShindcastElement_handler(name)

  end subroutine startElement_handler
  !============================================================================
  subroutine endElement_handler(namespaceURI,localname,name)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name

    ! Output Element
    if(name=='OutputDirectory') in_OutDir=.false.
    ! Topo Element
    call EtopoElement_handler(name)
    ! Time Element
    call EtimeElement_handler(name)
    ! Circulation Element
    call EcircElement_handler(name)
    ! Contour Element
    call EcontElement_handler(name)
    ! Wave Element
    call EwaveElement_handler(name)
    call EhindcastElement_handler(name)

  end subroutine endElement_handler
  ! ============================================================================
  subroutine characters_handler(chars)

    character(len=*),intent(in)::chars

    ! Get Output Dircetory Name
    if(in_OutDir)then
      outdir=''
      outdir=chars
    endif
    ! Topo Element
    if(in_Topo) call topo_characters_handler(chars)
    ! Get Time Parameters
    if(in_Time) call time_characters_handler(chars)
    ! Get Circulation Parameters
    if(in_Circ) call circ_characters_handler(chars)
    ! Get Contour Parameters
    if(in_Cont) call cont_characters_handler(chars)
    ! Get Wave Parameters
    if(in_wave) call wave_characters_handler(chars)
    if(in_hindcast) call hindcast_characters_handler(chars)

  end subroutine characters_handler
  ! ============================================================================
  subroutine StopoElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='Grid') in_Grid=.true.
    if(name=='GridX') in_GridXs=.true.
    if(name=='GridY') in_GridYs=.true.
    if(name=='GridSpace') in_GridSpaces=.true.
    if(name=='latGrid') in_latGrid=.true.
    if(name=='syncFolder') in_SyncFold=.true.
    if(name=='spmForcing') in_spmForcing=.true.
    if(name=='initDep') in_initDep=.true.
    if(name=='circForcing') in_circForcing=.true.
    if(name=='waveForcing') in_waveForcing=.true.
    if(name=='contourForcing') in_contourForcing=.true.

  end subroutine StopoElement_handler
  ! ============================================================================
  subroutine StimeElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='startTime') in_startTime=.true.
    if(name=='endTime') in_endTime=.true.
    if(name=='syncTime') in_syncTime=.true.
    if(name=='outputTime') in_outTime=.true.

  end subroutine StimeElement_handler
  ! ============================================================================
  subroutine ScircElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='Courant') in_CircStep=.true.
    if(name=='stormTime') in_CircTime=.true.
    if(name=='fricCoef') in_CircFric=.true.
    if(name=='filterStep') in_CircFilter=.true.
    if(name=='actionMaxDepth') in_max_depth=.true.

  end subroutine ScircElement_handler
  ! ============================================================================
  subroutine ScontElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='topContour') in_ContTop=.true.
    if(name=='botContour') in_ContBot=.true.
    if(name=='enteringBorder') in_ContEnt=.true.
    if(name=='exitingBorder') in_ContExit=.true.
    if(name=='meanVel') in_ContVel=.true.

  end subroutine ScontElement_handler
  ! ============================================================================
  subroutine SwaveElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='waveBase') in_wbase=.true.
    if(name=='nbForecast') in_wwcast=.true.
    if(name=='forecastClass') in_hindcast=.true.

  end subroutine SwaveElement_handler
  ! ============================================================================
  subroutine ShindcastElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='start') in_wtstart=.true.
    if(name=='end') in_wtend=.true.
    if(name=='subgroupNb') in_groupnb=.true.
    if(name=='subgroupParam') in_groupparam=.true.
    if(in_groupnb)then
      hclass=hclass+1
      sclass=0
    endif

  end subroutine ShindcastElement_handler
  ! ============================================================================
  subroutine EtopoElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='Grid') in_Grid=.false.
    if(name=='GridX') in_GridXs=.false.
    if(name=='GridY') in_GridYs=.false.
    if(name=='GridSpace') in_GridSpaces=.false.
    if(name=='latGrid') in_latGrid=.false.
    if(name=='syncFolder') in_SyncFold=.false.
    if(name=='spmForcing') in_spmForcing=.false.
    if(name=='initDep') in_initDep=.false.
    if(name=='circForcing') in_circForcing=.false.
    if(name=='waveForcing') in_waveForcing=.false.
    if(name=='contourForcing') in_contourForcing=.false.

  end subroutine EtopoElement_handler
  ! ============================================================================
  subroutine EtimeElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='startTime') in_startTime=.false.
    if(name=='endTime') in_endTime=.false.
    if(name=='syncTime') in_syncTime=.false.
    if(name=='outputTime') in_outTime=.false.

  end subroutine EtimeElement_handler
  ! ============================================================================
  subroutine EcircElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='Courant') in_CircStep=.false.
    if(name=='stormTime') in_CircTime=.false.
    if(name=='fricCoef') in_CircFric=.false.
    if(name=='filterStep') in_CircFilter=.false.
    if(name=='actionMaxDepth') in_max_depth=.false.

  end subroutine EcircElement_handler
  ! ============================================================================
  subroutine EcontElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='topContour') in_ContTop=.false.
    if(name=='botContour') in_ContBot=.false.
    if(name=='enteringBorder') in_ContEnt=.false.
    if(name=='exitingBorder') in_ContExit=.false.
    if(name=='meanVel') in_ContVel=.false.

  end subroutine EcontElement_handler
  ! ============================================================================
  subroutine EwaveElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='waveBase') in_wbase=.false.
    if(name=='nbForecast') in_wwcast=.false.
    if(name=='forecastClass') in_hindcast=.false.

  end subroutine EwaveElement_handler
  ! ============================================================================
  subroutine EhindcastElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='start') in_wtstart=.false.
    if(name=='end') in_wtend=.false.
    if(name=='subgroupNb') in_groupnb=.false.
    if(name=='subgroupParam') in_groupparam=.false.

  end subroutine EhindcastElement_handler
  ! ============================================================================
  subroutine topo_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_Grid)then
      xyzfile=chars
    elseif(in_GridXs) then
      grid_dxs=chars
      call rts(grid_dxs,sp_n)
    elseif(in_GridYs) then
      grid_dys=chars
      call rts(grid_dys,sp_m)
    elseif(in_GridSpaces) then
      grid_spaces=chars
      call rts(grid_spaces,stratal_dx)
    elseif(in_latGrid) then
      grid_lat=chars
      call rts(grid_lat,alat)
    elseif(in_SyncFold) then
      syncfolder=chars
    elseif(in_spmForcing) then
      grid_spmForcing=chars
      call rts(grid_spmForcing,spmflag)
    elseif(in_initDep) then
      grid_initDep=chars
      call rts(grid_initDep,initDep)
    elseif(in_circForcing) then
      grid_circForcing=chars
      call rts(grid_circForcing,circflag)
    elseif(in_waveForcing) then
      grid_waveForcing=chars
      call rts(grid_waveForcing,waveflag)
    elseif(in_contourForcing) then
      grid_contourForcing=chars
      call rts(grid_contourForcing,contflag)
    endif

  end subroutine topo_characters_handler
  ! ============================================================================
  subroutine time_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_startTime) then
      startTime=chars
      call rts(startTime,time_start)
    elseif(in_endTime)then
      endTime=chars
      call rts(endTime,time_end)
    elseif(in_syncTime)then
      syncTime=chars
      call rts(syncTime,time_step)
    elseif(in_outTime)then
      outTime=chars
      call rts(outTime,out_step)
    endif

  end subroutine time_characters_handler
  ! ============================================================================
  subroutine circ_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_CircStep) then
      circstep=chars
      call rts(circstep,Courant)
    elseif(in_CircTime)then
      circtime=chars
      call rts(circtime,kend)
    elseif(in_CircFric)then
      circfric=chars
      call rts(circfric,cfric)
    elseif(in_CircFilter)then
      circfilter=chars
      call rts(circfilter,ifilt)
    elseif(in_max_depth)then
      max_depth=chars
      call rts(max_depth,maxdepth)
    endif

  end subroutine circ_characters_handler
  ! ============================================================================
  subroutine cont_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_ContTop) then
      conttop=chars
      call rts(conttop,contour1)
    elseif(in_ContBot)then
      contbot=chars
      call rts(contbot,contour2)
    elseif(in_ContEnt)then
      contEnt=chars
      call rts(contEnt,inside)
    elseif(in_ContExit)then
      contExit=chars
      call rts(contExit,outside)
    elseif(in_ContVel)then
      contVel=chars
      call rts(contVel,vmean_cont)
    endif

  end subroutine cont_characters_handler
  ! ============================================================================
  subroutine wave_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_wbase)then
      wbase=chars
      call rts(wbase,wave_base)
    elseif(in_wwcast)then
      wwcast=chars
      call rts(wwcast,forecast_nb)
      allocate(hindcast(forecast_nb))
    endif

  end subroutine wave_characters_handler
  ! ============================================================================
  subroutine hindcast_characters_handler(chars)

    character(len=*),intent(in)::chars
    real::sbwi( 3 )

    if(in_wtstart)then
      wtstart=chars
      call rts(wtstart,hindcast(hclass)%tstart)
    elseif(in_wtend)then
      wtend=chars
      call rts(wtend,hindcast(hclass)%tend)
    elseif(in_groupnb)then
      groupnb=chars
      call rts(groupnb,hindcast(hclass)%cnb)
      allocate(hindcast(hclass)%subgroup(hindcast(hclass)%cnb))
    elseif(in_groupparam)then
      sclass=sclass+1
      groupparam=chars
      if(waveflag==1)then
        call rts(groupparam,sbwi)
        hindcast(hclass)%subgroup(sclass)%oc=sbwi(1)
        hindcast(hclass)%subgroup(sclass)%wvel=sbwi(2)
        hindcast(hclass)%subgroup(sclass)%wdir=sbwi(3)
        hindcast(hclass)%subgroup(sclass)%hs=0.0_8
        hindcast(hclass)%subgroup(sclass)%per=0.0_8
        hindcast(hclass)%subgroup(sclass)%dir=0.0_8
        hindcast(hclass)%subgroup(sclass)%dd=0.0_8
      elseif(waveflag==0.and.circflag==1)then
        call rts(groupparam,sbwi)
        hindcast(hclass)%subgroup(sclass)%oc=sbwi(1)
        hindcast(hclass)%subgroup(sclass)%wvel=sbwi(2)
        hindcast(hclass)%subgroup(sclass)%wdir=sbwi(3)
        hindcast(hclass)%subgroup(sclass)%hs=0.0_8
        hindcast(hclass)%subgroup(sclass)%per=0.0_8
        hindcast(hclass)%subgroup(sclass)%dir=0.0_8
        hindcast(hclass)%subgroup(sclass)%dd=0.0_8
      endif
    endif

  end subroutine hindcast_characters_handler
  ! ============================================================================
  subroutine endDocument_handler

  end subroutine endDocument_handler
  ! ============================================================================
  subroutine xml_Data

    type(xml_t)::xf
    integer::ios

    hclass=0
    sclass=0
    ! Open file
    call open_xml_file(xf,finput ,ios)
    if(ios/=0) then
      print*,'---------------------'
      print*,'Error opening input file for parsing XmL'
      print*,'---------------------'
    endif
    ! Parse file at first and find data
    call parse(xf,&
      startDocument_handler=startDocument_handler,&
      startElement_handler=startElement_handler,&
      endElement_handler=endElement_handler,&
      characters_handler=characters_handler,&
      endDocument_handler=endDocument_handler)
    ! Close file
    call close_xml_t(xf)

  end subroutine xml_Data
  ! ============================================================================

end module xml_reader
