!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP3_decompDiagMod
  implicit none
  private
  public :: HYMAP3_decomp_count,HYMAP3_decomp_write_nc
contains

  !ag(23Nov2026) Count hydrologically connected interprocessor river links.
  ! Each upstream -> downstream pair is ONE link, independent of direction
  ! or number of MPI messages exchanged for the routing physics.
  ! Flag per river element: 0 = local only, 1 = remote downstream,
  ! 2 = remote upstream, 3 = both. Outlets are not counted as links.
  subroutine HYMAP3_decomp_count(next,outlet,owner,np,load,nout,nin,ncell,flag,ncut,ierr)
    implicit none
    integer, intent(in) :: next(:),outlet(:),owner(:),np
    integer, intent(out) :: load(:),nout(:),nin(:),ncell(:),flag(:),ncut,ierr
    integer :: ng,i,j,src,dst,p

    ierr=0
    ncut=0
    load=0
    nout=0
    nin=0
    ncell=0
    flag=0
    ng=size(owner)
    if(np<1.or.ng<1.or.size(next)/=ng.or.size(outlet)/=ng.or. &
         size(flag)/=ng.or.size(load)/=np.or.size(nout)/=np.or. &
         size(nin)/=np.or.size(ncell)/=np)then
       ierr=1
       return
    endif
    if(any(owner<0).or.any(owner>=np))then
       ierr=2
       return
    endif
    do i=1,ng
       p=owner(i)+1
       load(p)=load(p)+1
    enddo
    do i=1,ng
       if(outlet(i)/=0)cycle
       j=next(i)
       if(j<1.or.j>ng)then
          ierr=3
          return
       endif
       src=owner(i)+1
       dst=owner(j)+1
       if(src==dst)cycle
       ncut=ncut+1
       nout(src)=nout(src)+1
       nin(dst)=nin(dst)+1
       if(iand(flag(i),1)==0)then
          if(flag(i)==0)ncell(src)=ncell(src)+1
          flag(i)=ior(flag(i),1)
       endif
       if(iand(flag(j),2)==0)then
          ! Count a grid cell only once, even when it is both a sender
          ! and a receiver or receives from multiple remote tributaries.
          if(flag(j)==0)ncell(dst)=ncell(dst)+1
          flag(j)=ior(flag(j),2)
       endif
    enddo
    if(sum(nout)/=ncut.or.sum(nin)/=ncut.or.sum(load)/=ng)then
       ierr=4
    endif
  end subroutine HYMAP3_decomp_count

  !ag(23Nov2026) One initialization-only, optional spatial NetCDF diagnostic.
  ! Outputs two integer grid variables, NOT routing prognostic variables.
  ! Neither processor attribution nor model physics is altered.
  subroutine HYMAP3_decomp_write_nc(file,sindex,owner,flag,ierr,lon,lat,errmsg)
    !ag(23Nov2026) NetCDF is required by LIS output on Discover.
    ! Do not hide this writer behind USE_NETCDF3/4: in the reported
    ! LIS build neither symbol was defined, causing ierr=3 at runtime.
    use netcdf
    implicit none
    character(*), intent(in) :: file
    integer, intent(in) :: sindex(:,:),owner(:),flag(:)
    integer, intent(out) :: ierr
    real, optional, intent(in) :: lon(:),lat(:)
    !ag(23Nov2026) Return the NetCDF operation and its actual error string.
    ! Diagnostic output must not abort otherwise valid routing simulations.
    character(*), optional, intent(out) :: errmsg
    integer, allocatable :: processor_grid(:,:),link_grid(:,:)
    integer :: nx,ny,ix,iy,k
    integer :: ncid,dx,dy,pid,linkid,latid,lonid,st,st2
    logical :: hasgeo
    character(96) :: stage
    ierr=0
    if(present(errmsg))errmsg=''
    nx=size(sindex,1)
    ny=size(sindex,2)
    if(nx<1.or.ny<1.or.size(owner)/=size(flag))then
       ierr=1
       if(present(errmsg))errmsg='Invalid grid or ownership dimensions'
       return
    endif
    allocate(processor_grid(nx,ny),link_grid(nx,ny))
    processor_grid=-1
    link_grid=-1
    do iy=1,ny
       do ix=1,nx
          k=sindex(ix,iy)
          if(k<=0)cycle
          if(k>size(owner))then
             ierr=2
             if(present(errmsg))errmsg='Invalid global river ID in sindex'
             return
          endif
          processor_grid(ix,iy)=owner(k)
          link_grid(ix,iy)=flag(k)
       enddo
    enddo
    hasgeo=.false.
    if(present(lon).and.present(lat))then
       hasgeo=size(lon)==nx*ny.and.size(lat)==nx*ny
    endif
    ncid=-1
    stage='nf90_create (check destination directory permissions)'
    st=nf90_create(trim(file),NF90_CLOBBER,ncid)
    if(st/=nf90_noerr)goto 900
    stage='nf90_def_dim x'
    st=nf90_def_dim(ncid,'x',nx,dx)
    if(st/=nf90_noerr)goto 900
    stage='nf90_def_dim y'
    st=nf90_def_dim(ncid,'y',ny,dy)
    if(st/=nf90_noerr)goto 900
    stage='nf90_def_var processor_id'
    st=nf90_def_var(ncid,'processor_id',NF90_INT,(/dx,dy/),pid)
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att processor long_name'
    st=nf90_put_att(ncid,pid,'long_name','MPI rank owning HyMAP routing grid cell')
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att processor _FillValue'
    st=nf90_put_att(ncid,pid,'_FillValue',-1)
    if(st/=nf90_noerr)goto 900
    stage='nf90_def_var cross_processor_connection'
    st=nf90_def_var(ncid,'cross_processor_connection',NF90_INT,(/dx,dy/),linkid)
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att connection long_name'
    st=nf90_put_att(ncid,linkid,'long_name','Hydrologic link to another MPI processor')
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att connection flag_values'
    st=nf90_put_att(ncid,linkid,'flag_values',(/0,1,2,3/))
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att connection flag_meanings'
    st=nf90_put_att(ncid,linkid,'flag_meanings', &
         'none remote_downstream remote_upstream both_directions')
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att connection _FillValue'
    st=nf90_put_att(ncid,linkid,'_FillValue',-1)
    if(st/=nf90_noerr)goto 900
    if(hasgeo)then
       stage='nf90_def_var longitude'
       st=nf90_def_var(ncid,'longitude',NF90_FLOAT,(/dx,dy/),lonid)
       if(st/=nf90_noerr)goto 900
       stage='nf90_put_att longitude units'
       st=nf90_put_att(ncid,lonid,'units','degrees_east')
       if(st/=nf90_noerr)goto 900
       stage='nf90_def_var latitude'
       st=nf90_def_var(ncid,'latitude',NF90_FLOAT,(/dx,dy/),latid)
       if(st/=nf90_noerr)goto 900
       stage='nf90_put_att latitude units'
       st=nf90_put_att(ncid,latid,'units','degrees_north')
       if(st/=nf90_noerr)goto 900
       stage='nf90_put_att processor coordinates'
       st=nf90_put_att(ncid,pid,'coordinates','longitude latitude')
       if(st/=nf90_noerr)goto 900
       stage='nf90_put_att connection coordinates'
       st=nf90_put_att(ncid,linkid,'coordinates','longitude latitude')
       if(st/=nf90_noerr)goto 900
    endif
    stage='nf90_put_att global description'
    st=nf90_put_att(ncid,NF90_GLOBAL,'description', &
         'HyMAP gridded-domain MPI ownership; computed once at initialization')
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_att global link_definition'
    st=nf90_put_att(ncid,NF90_GLOBAL,'link_definition', &
         'Cross-processor directed river edges, excluding outlet cells')
    if(st/=nf90_noerr)goto 900
    stage='nf90_enddef'
    st=nf90_enddef(ncid)
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_var processor_id'
    st=nf90_put_var(ncid,pid,processor_grid)
    if(st/=nf90_noerr)goto 900
    stage='nf90_put_var cross_processor_connection'
    st=nf90_put_var(ncid,linkid,link_grid)
    if(st/=nf90_noerr)goto 900
    if(hasgeo)then
       stage='nf90_put_var longitude'
       st=nf90_put_var(ncid,lonid,reshape(lon,(/nx,ny/)))
       if(st/=nf90_noerr)goto 900
       stage='nf90_put_var latitude'
       st=nf90_put_var(ncid,latid,reshape(lat,(/nx,ny/)))
       if(st/=nf90_noerr)goto 900
    endif
    stage='nf90_close'
    st=nf90_close(ncid)
    if(st/=nf90_noerr)then
       ierr=100+abs(st)
       if(present(errmsg))errmsg=trim(stage)//': '//trim(nf90_strerror(st))
    endif
    return
900 continue
    ierr=100+abs(st)
    if(present(errmsg))errmsg=trim(stage)//': '//trim(nf90_strerror(st))
    if(ncid>=0)then
       st2=nf90_close(ncid)
    endif
  end subroutine HYMAP3_decomp_write_nc
end module HYMAP3_decompDiagMod
