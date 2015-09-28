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
!      filename:  OceanOut.f90
!
!    Description:  Implements the XdmF and HdF5 SPM output of surface velocity evolution
!
!        Version:  1.0
!        Created:  21/09/15 09:45:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================
module OceanOut

    use hdf5
    use currents_io
    use currents_data
    use precision_data
    use FoX_wxml
    use mpidata

    implicit none

    public

contains

    ! =====================================================================================
    subroutine xdmf_output( iter, sgp )

      integer::iter,sgp

      ! Output is generated in serial
      if(iam==0)then
        call Ocean_hdf5(iter,sgp)
        call Ocean_xmf(iter)
        call Ocean_series(iter)
      endif
      call mpi_barrier(ocean_comm_world,ierr)

      return

    end subroutine xdmf_output
    ! =====================================================================================
    subroutine Ocean_hdf5(iter,sgp)

      logical::compression
      integer::i,j,p,n,id,k,iter,sgp,rank
      character(len=128)::text,stg,file
      integer(hid_t)::file_id,plist_id
      integer(hid_t)::filespace,dset_id
      integer(hsize_t),dimension(2)::dims

     file=''
     file='OceanGrid.'
      call noblnk(file)
      call append_zero(file,iter-1)
      stg='.h5'
      call append_str(file,stg)
      call addpath2(file)
      ! Create nodes arrays
      id=1
      n=1
      do i=1,sp_m
        p=i
        do j=1,sp_n
          nodes(id)=real((j-1)*stratal_dx/100+stratal_xo)
          nodes(id+1)=real((i-1)*stratal_dx/100+stratal_yo)
          nodes(id+2)=0.0
          htopo(n)=real(sp_topo(i,j))
          do k=1,maxsubgroup
            if(k>sgp)then
              curU(k,n)=0.0
              curV(k,n)=0.0
              wavU(k,n)=0.0
              wavV(k,n)=0.0
            endif
          enddo
          n=n+1
          p=p+sp_m
          id=id+3
        enddo
      enddo
      ! Initialize predefined datatypes
      call h5open_f(ierr)
      call h5zfilter_avail_f(h5z_filter_deflate_f,compression,ierr)
      ! Setupfileaccess property list for MPI-IO access.
      call h5pcreate_f(h5p_file_access_f,plist_id,ierr)
      ! Create thefilecollectively.
      call h5fcreate_f(file,h5f_acc_trunc_f,file_id,ierr,access_prp=plist_id)
      ! ========================
      ! The Coordinates - vertices
      ! ========================
      dims(1)=3
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/vertices"
      if(.not.compression)then
        dims(1)=1
        dims(2)=totnodes*3
        ! Create the dataset with default properties
        call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr)
        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f,plist_id,ierr)
        ! Write the dataset collectively
        call h5dwrite_f(dset_id,h5t_native_real,nodes,dims,ierr,file_space_id=filespace,xfer_prp=plist_id)
      else
        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
        call h5pset_deflate_f(plist_id,9,ierr)
        call h5pset_chunk_f(plist_id,rank,dims,ierr)
        dims(1)=1
        dims(2)=totnodes*3
        ! Create the dataset with default properties
        call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr,plist_id)
        ! Write the dataset collectively
        call h5dwrite_f(dset_id,h5t_native_real,nodes,dims,ierr)
        call h5pclose_f(plist_id,ierr)
      endif
      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)
      dims(1)=1
      dims(2)=totnodes
      rank=2
      call h5screate_simple_f(rank,dims,filespace,ierr)
      text=''
      text="/Bathy"
      if(.not.compression)then
        ! Create the dataset with default properties
        call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr)
        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f,plist_id,ierr)
        ! Write the dataset
        call h5dwrite_f(dset_id,h5t_native_real,htopo(1:totnodes),dims,ierr,file_space_id=filespace,xfer_prp=plist_id)
      else
        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
        call h5pset_deflate_f(plist_id,9,ierr)
        call h5pset_chunk_f(plist_id,rank,dims,ierr)
        ! Create the dataset with default properties
        call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr,plist_id)
        ! Write the dataset collectively
        call h5dwrite_f(dset_id,h5t_native_real,htopo(1:totnodes),dims,ierr)
        call h5pclose_f(plist_id,ierr)
      endif
      ! Close the dataset
      call h5dclose_f(dset_id,ierr)
      call h5sclose_f(filespace,ierr)

      ! Loop over the number of subgroup and attribute ocean circulation adequately
      do p=1,maxsubgroup
        ! ========================
        ! Current velocity attribute
        ! ========================
        dims(1)=1
        dims(2)=totnodes
        rank=2
        call h5screate_simple_f(rank,dims,filespace,ierr)
        text=''
        text="/CurU"
        call append_nb(text,p)
        if(.not.compression)then
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr)
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_xfer_f,plist_id,ierr)
          ! Write the dataset
          call h5dwrite_f(dset_id,h5t_native_real,curU(p,1:totnodes),dims,ierr,file_space_id=filespace,xfer_prp=plist_id)
        else
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
          call h5pset_deflate_f(plist_id,9,ierr)
          call h5pset_chunk_f(plist_id,rank,dims,ierr)
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr,plist_id)
          ! Write the dataset collectively
          call h5dwrite_f(dset_id,h5t_native_real,curU(p,1:totnodes),dims,ierr)
          call h5pclose_f(plist_id,ierr)
        endif
        ! Close the dataset
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(filespace,ierr)
        ! ========================
        ! Current direction attribute
        ! ========================
        dims(1)=1
        dims(2)=totnodes
        rank=2
        call h5screate_simple_f(rank,dims,filespace,ierr)
        text=''
        text="/CurV"
        call append_nb(text,p)
        if(.not.compression)then
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr)
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_xfer_f,plist_id,ierr)
          ! Write the dataset
          call h5dwrite_f(dset_id,h5t_native_real,curV(p,1:totnodes),dims,ierr,file_space_id=filespace,xfer_prp=plist_id)
        else
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
          call h5pset_deflate_f(plist_id,9,ierr)
          call h5pset_chunk_f(plist_id,rank,dims,ierr)
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr,plist_id)
          ! Write the dataset collectively
          call h5dwrite_f(dset_id,h5t_native_real,curV(p,1:totnodes),dims,ierr)
          call h5pclose_f(plist_id,ierr)
        endif
        ! Close the dataset
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(filespace,ierr)

        ! ========================
        ! Wave velocity attribute
        ! ========================
        dims(1)=1
        dims(2)=totnodes
        rank=2
        call h5screate_simple_f(rank,dims,filespace,ierr)
        text=''
        text="/WavU"
        call append_nb(text,p)
        if(.not.compression)then
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr)
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_xfer_f,plist_id,ierr)
          ! Write the dataset
          call h5dwrite_f(dset_id,h5t_native_real,wavU(p,1:totnodes),dims,ierr,file_space_id=filespace,xfer_prp=plist_id)
        else
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
          call h5pset_deflate_f(plist_id,9,ierr)
          call h5pset_chunk_f(plist_id,rank,dims,ierr)
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr,plist_id)
          ! Write the dataset collectively
          call h5dwrite_f(dset_id,h5t_native_real,wavU(p,1:totnodes),dims,ierr)
          call h5pclose_f(plist_id,ierr)
        endif
        ! Close the dataset
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(filespace,ierr)
        ! ========================
        ! Wave direction attribute
        ! ========================
        dims(1)=1
        dims(2)=totnodes
        rank=2
        call h5screate_simple_f(rank,dims,filespace,ierr)
        text=''
        text="/WavV"
        call append_nb(text,p)
        if(.not.compression)then
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr)
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_xfer_f,plist_id,ierr)
          ! Write the dataset
          call h5dwrite_f(dset_id,h5t_native_real,wavV(p,1:totnodes),dims,ierr,file_space_id=filespace,xfer_prp=plist_id)
        else
          ! Create property list for collective dataset write
          call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
          call h5pset_deflate_f(plist_id,9,ierr)
          call h5pset_chunk_f(plist_id,rank, dims,ierr)
          ! Create the dataset with default properties
          call h5dcreate_f(file_id,trim(text),h5t_native_real,filespace,dset_id,ierr,plist_id)
          ! Write the dataset collectively
          call h5dwrite_f(dset_id,h5t_native_real,wavV(p,1:totnodes),dims,ierr)
          call h5pclose_f(plist_id,ierr)
        endif
        ! Close the dataset
        call h5dclose_f(dset_id,ierr)
        call h5sclose_f(filespace,ierr)
      enddo

      ! Close thefile.
      call h5fclose_f(file_id,ierr)
      ! Close interface
      call h5close_f(ierr)


      return

    end subroutine Ocean_hdf5
    ! =====================================================================================
    subroutine Ocean_xmf(iter)

        ! Parameters Declaration
        type(xmlf_t)::xf

        integer::iter,p
        character(len=128)::str,txt,filename,filename0,filename1,filename2,file
        character(len=128)::filename3,filename4,filename5,filename6

        file=''
        file='Ocean.'
        call noblnk(file)
        call append_zero(file,iter-1)
        str='.xmf'
        call append_str(file,str)
        call addpath2(file)
        call xml_OpenFile(file,xf)
        ! Header
        call xml_AddDOCTYPE(xf,"Xdmf","Xdmf.dtd")
        call xml_DeclareNamespace(xf,"http://www.w3.org/2001/XInclude","xi")
        call xml_NewElement(xf,"Xdmf")
        call xml_AddAttribute(xf,"Version","2.0")
        call xml_NewElement(xf,"Domain")
        call xml_NewElement(xf,"Grid")
        call xml_AddAttribute(xf,"GridType","Collection")
        call xml_AddAttribute(xf,"CollectionType","Spatial")
        call xml_NewElement(xf,"Time")
        call xml_AddAttribute(xf,"Type","Single")
        call xml_AddAttribute(xf,"Value",tnow)
        call xml_EndElement(xf,"Time")
        filename=''
        filename='OceanGrid.'
        call noblnk(filename)
        call append_zero(filename,iter-1)
        str='.h5'
        call append_str(filename,str)
        filename0=filename
        filename1=filename
        filename2=filename
        filename3=filename
        filename4=filename
        filename5=filename
        str=':/Bathy'
        call append_str(filename0,str)
        str=':/vertices'
        call append_str(filename1,str)
        str=':/CurU'
        call append_str(filename2,str)
        str=':/CurV'
        call append_str(filename3,str)
        str=':/WavU'
        call append_str(filename4,str)
        str=':/WavV'
        call append_str(filename5,str)

        ! Block begin
        call xml_NewElement(xf,"Grid")
        str='OBlock.'
        call append_zero(str,iter-1)
        call xml_AddAttribute(xf,"Name",trim(str))
        call xml_NewElement(xf,"Topology")
        call xml_AddAttribute(xf,"Type","3DSMesh")
        str=' '
        call append_nb2(str,1)
        call append_nb2(str,sp_m)
        call append_nb2(str,sp_n)
        call xml_AddAttribute(xf,"Dimensions",trim(str))
        call xml_EndElement(xf,"Topology")

        ! Geometry
        call xml_NewElement(xf,"Geometry")
        call xml_AddAttribute(xf,"Type","XYZ")
        call xml_NewElement(xf,"DataItem")
        call xml_AddAttribute(xf,"Format","HDF")
        call xml_AddAttribute(xf,"NumberType","Float")
        call xml_AddAttribute(xf,"Precision","4")
        str=' '
        call append_nb2(str,sp_m)
        call append_nb2(str,sp_n)
        call append_nb2(str,1)
        call append_nb2(str,3)
        call xml_AddAttribute(xf,"Dimensions",trim(str))
        call xml_AddCharacters(xf,trim(filename1))
        call xml_EndElement(xf,"DataItem")
        call xml_EndElement(xf,"Geometry")
        str=' '
        call append_nb2(str,sp_m)
        call append_nb2(str,sp_n)
        call append_nb2(str,1)

        ! Bathymetry
        txt='Bathy'
        call xml_NewElement(xf,"Attribute")
        call xml_AddAttribute(xf,"Type","Scalar")
        call xml_AddAttribute(xf,"Center","Node")
        call xml_AddAttribute(xf,"Name",trim(txt))
        call xml_NewElement(xf,"DataItem")
        call xml_AddAttribute(xf,"Format","HDF")
        call xml_AddAttribute(xf,"NumberType","Float")
        call xml_AddAttribute(xf,"Precision","4")
        call xml_AddAttribute(xf,"Dimensions",trim(str))
        call xml_AddCharacters(xf,trim(filename0))
        call xml_EndElement(xf,"DataItem")
        call xml_EndElement(xf,"Attribute")

        do p=1,maxsubgroup
          ! Current velocity
          filename6=filename2
          call append_nb(filename6,p)
          txt='CurU_grp'
          call append_nb(txt,p)
          call xml_NewElement(xf,"Attribute")
          call xml_AddAttribute(xf,"Type","Scalar")
          call xml_AddAttribute(xf,"Center","Node")
          call xml_AddAttribute(xf,"Name",trim(txt))
          call xml_NewElement(xf,"DataItem")
          call xml_AddAttribute(xf,"Format","HDF")
          call xml_AddAttribute(xf,"NumberType","Float")
          call xml_AddAttribute(xf,"Precision","4")
          call xml_AddAttribute(xf,"Dimensions",trim(str))
          call xml_AddCharacters(xf,trim(filename6))
          call xml_EndElement(xf,"DataItem")
          call xml_EndElement(xf,"Attribute")

          ! Current direction
          filename6=filename3
          call append_nb(filename6,p)
          txt='CurV_grp'
          call append_nb(txt,p)
          call xml_NewElement(xf,"Attribute")
          call xml_AddAttribute(xf,"Type","Scalar")
          call xml_AddAttribute(xf,"Center","Node")
          call xml_AddAttribute(xf,"Name",trim(txt))
          call xml_NewElement(xf,"DataItem")
          call xml_AddAttribute(xf,"Format","HDF")
          call xml_AddAttribute(xf,"NumberType","Float")
          call xml_AddAttribute(xf,"Precision","4")
          call xml_AddAttribute(xf,"Dimensions",trim(str))
          call xml_AddCharacters(xf,trim(filename6))
          call xml_EndElement(xf,"DataItem")
          call xml_EndElement(xf,"Attribute")

          ! Current velocity
          filename6=filename4
          call append_nb(filename6,p)
          txt='WavU_grp'
          call append_nb(txt,p)
          call xml_NewElement(xf,"Attribute")
          call xml_AddAttribute(xf,"Type","Scalar")
          call xml_AddAttribute(xf,"Center","Node")
          call xml_AddAttribute(xf,"Name",trim(txt))
          call xml_NewElement(xf,"DataItem")
          call xml_AddAttribute(xf,"Format","HDF")
          call xml_AddAttribute(xf,"NumberType","Float")
          call xml_AddAttribute(xf,"Precision","4")
          call xml_AddAttribute(xf,"Dimensions",trim(str))
          call xml_AddCharacters(xf,trim(filename6))
          call xml_EndElement(xf,"DataItem")
          call xml_EndElement(xf,"Attribute")

          ! Current direction
          filename6=filename5
          call append_nb(filename6,p)
          txt='WavV_grp'
          call append_nb(txt,p)
          call xml_NewElement(xf,"Attribute")
          call xml_AddAttribute(xf,"Type","Scalar")
          call xml_AddAttribute(xf,"Center","Node")
          call xml_AddAttribute(xf,"Name",trim(txt))
          call xml_NewElement(xf,"DataItem")
          call xml_AddAttribute(xf,"Format","HDF")
          call xml_AddAttribute(xf,"NumberType","Float")
          call xml_AddAttribute(xf,"Precision","4")
          call xml_AddAttribute(xf,"Dimensions",trim(str))
          call xml_AddCharacters(xf,trim(filename6))
          call xml_EndElement(xf,"DataItem")
          call xml_EndElement(xf,"Attribute")
        enddo
        call xml_EndElement(xf,"Grid")

        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)

        return

    end subroutine Ocean_xmf
    ! =====================================================================================
    subroutine Ocean_series(iter)

      ! Parameters Declaration
      type(xmlf_t)::xf

      integer::i,iter,it0
      character(len=128)::filename,str,fname

      filename='Ocean_series.xdmf'
      call addpath2(filename)
      ! Header
      call xml_OpenFile(filename,xf)
      call xml_AddDOCTYPE(xf,"Xdmf","Xdmf.dtd")
      call xml_DeclareNamespace(xf,"http://www.w3.org/2001/XInclude","xi")
      call xml_NewElement(xf,"Xdmf")
      call xml_AddAttribute(xf,"Version","2.0")
      call xml_NewElement(xf,"Domain")
      call xml_NewElement(xf,"Grid")
      call xml_AddAttribute(xf,"GridType","Collection")
      call xml_AddAttribute(xf,"CollectionType","Temporal")
      it0=1
      ! Loop over time step
      do i=it0,iter
        ! Grid name
        fname=''
        fname='Ocean.'
        call noblnk(fname)
        call append_zero(fname,i-1)
        str='.xmf'
        call append_str(fname,str)
        call xml_NewElement(xf,"xi:include")
        call xml_AddAttribute(xf,"href",trim(fname))
        call xml_AddAttribute(xf,"xpointer","xpointer(//Xdmf/Domain/Grid)")
        call xml_EndElement(xf,"xi:include")
      enddo
      ! Footer
      call xml_EndElement(xf,"Grid")
      call xml_EndElement(xf,"Domain")
      call xml_EndElement(xf,"Xdmf")
      call xml_Close(xf)

      return

    end subroutine Ocean_series
    ! =====================================================================================

end module OceanOut
! =====================================================================================
