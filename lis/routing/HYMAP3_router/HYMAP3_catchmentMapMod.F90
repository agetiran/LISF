!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP3_catchmentMapMod
  implicit none
  private
  public :: HYMAP3_catchment_transfer_dec
  public :: HYMAP3_catchment_map_sources, HYMAP3_catchment_build_transfer
  public :: HYMAP3_catchment_pack_grid, HYMAP3_catchment_unpack_vector
  public :: HYMAP3_catchment_pack_vector, HYMAP3_catchment_unpack_grid

  !ag(22Nov2026)
  ! A deterministic LIS-rectangle <-> catchment transfer plan.
  ! Peer lists share increasing global-ID order. The communication
  ! buffers are allocated once by HYMAP3_catchmentMPIMod and reused.
  type :: HYMAP3_catchment_transfer_dec
     integer, allocatable :: send_count(:),recv_count(:)
     integer, allocatable :: send_offset(:),recv_offset(:)
     integer, allocatable :: send_gid(:),recv_gid(:)
     integer, allocatable :: recv_local(:)
     real, allocatable :: buffer_out(:,:),buffer_in(:,:)
     integer :: nbuffer_fields = 0
  end type HYMAP3_catchment_transfer_dec
contains
  ! Use HyMAP's existing global sindex grid to avoid allocating a
  ! replicated nx*ny temporary processor grid on every MPI rank.
  ! This is O(Ngrid + Nriver) work and O(Nriver) additional memory.
  subroutine HYMAP3_catchment_map_sources(sindex,np,xstart,xend,ystart,yend,source,ierr,message)
    integer, intent(in) :: sindex(:,:),np
    integer, intent(in) :: xstart(:),xend(:),ystart(:),yend(:)
    integer, intent(out) :: source(:),ierr
    character(*), intent(out) :: message
    integer :: p,ix,iy,gid,nx,ny
    ierr=0
    message='OK'
    nx=size(sindex,1)
    ny=size(sindex,2)
    source=-1
    if(nx<1.or.ny<1.or.np<1.or. &
         size(xstart)/=np.or.size(xend)/=np.or. &
         size(ystart)/=np.or.size(yend)/=np)then
       ierr=1
       message='Invalid LIS source-grid dimensions'
       return
    endif
    do p=1,np
       if(xstart(p)<1.or.xend(p)>nx.or.ystart(p)<1.or.yend(p)>ny.or. &
            xstart(p)>xend(p).or.ystart(p)>yend(p))then
          ierr=2
          message='Invalid LIS processor grid bounds'
          return
       endif
       do iy=ystart(p),yend(p)
          do ix=xstart(p),xend(p)
             gid=sindex(ix,iy)
             if(gid<=0)cycle
             if(gid>size(source))then
                ierr=4
                message='HyMAP global sequence index outside network'
                return
             endif
             if(source(gid)/=-1)then
                ierr=3
                message='Overlapping LIS interior grid ownership'
                return
             endif
             source(gid)=p-1
          enddo
       enddo
    enddo
    if(any(source<0))then
       ierr=5
       message='Active HyMAP cell has no LIS grid owner'
    endif
  end subroutine HYMAP3_catchment_map_sources

  subroutine HYMAP3_catchment_build_transfer(source,owner,np,rank,plan,ierr,message)
    integer, intent(in) :: source(:),owner(:),np,rank
    type(HYMAP3_catchment_transfer_dec), intent(out) :: plan
    integer, intent(out) :: ierr
    character(*), intent(out) :: message
    integer, allocatable :: scursor(:),rcursor(:),local(:)
    integer :: i,p,j,n,k
    ierr=0
    message='OK'
    n=size(source)
    if(np<1.or.rank<0.or.rank>=np.or.size(owner)/=n)then
       ierr=1
       message='Invalid processor ranks or map lengths'
       return
    endif
    if(any(source<0).or.any(source>=np).or.any(owner<0).or.any(owner>=np))then
       ierr=2
       message='Invalid source or catchment processor ID'
       return
    endif
    allocate(plan%send_count(0:np-1),plan%recv_count(0:np-1))
    allocate(plan%send_offset(0:np),plan%recv_offset(0:np))
    plan%send_count=0
    plan%recv_count=0
    do i=1,n
       if(source(i)==rank)plan%send_count(owner(i))=plan%send_count(owner(i))+1
       if(owner(i)==rank)plan%recv_count(source(i))=plan%recv_count(source(i))+1
    enddo
    plan%send_offset(0)=0
    plan%recv_offset(0)=0
    do p=0,np-1
       plan%send_offset(p+1)=plan%send_offset(p)+plan%send_count(p)
       plan%recv_offset(p+1)=plan%recv_offset(p)+plan%recv_count(p)
    enddo
    allocate(plan%send_gid(plan%send_offset(np)))
    allocate(plan%recv_gid(plan%recv_offset(np)))
    allocate(plan%recv_local(plan%recv_offset(np)))
    allocate(scursor(0:np-1),rcursor(0:np-1),local(n))
    scursor=plan%send_offset(0:np-1)
    rcursor=plan%recv_offset(0:np-1)
    local=0
    k=0
    do i=1,n
       if(owner(i)==rank)then
          k=k+1
          local(i)=k
       endif
    enddo
    do i=1,n
       if(source(i)==rank)then
          p=owner(i)
          scursor(p)=scursor(p)+1
          plan%send_gid(scursor(p))=i
       endif
       if(owner(i)==rank)then
          p=source(i)
          rcursor(p)=rcursor(p)+1
          j=rcursor(p)
          plan%recv_gid(j)=i
          plan%recv_local(j)=local(i)
       endif
    enddo
  end subroutine HYMAP3_catchment_build_transfer

  !ag(22Nov2026) Packing the LIS->HyMAP exchange in stable global-ID order.
  ! Actual MPI communication is in HYMAP3_catchmentMPIMod.
  subroutine HYMAP3_catchment_pack_grid(plan,gx,gy,west,south,grid,buffer,ierr)
    type(HYMAP3_catchment_transfer_dec), intent(in) :: plan
    integer, intent(in) :: gx(:),gy(:),west,south
    real, intent(in) :: grid(:,:,:)
    real, intent(out) :: buffer(:,:)
    integer, intent(out) :: ierr
    integer :: i,gid,ix,iy
    ierr=0
    if(size(buffer,1)/=size(grid,3).or.size(buffer,2)/=size(plan%send_gid))then
       ierr=1
       return
    endif
    do i=1,size(plan%send_gid)
       gid=plan%send_gid(i)
       if(gid<1.or.gid>size(gx).or.gid>size(gy))then
          ierr=2
          return
       endif
       ix=gx(gid)-west+1
       iy=gy(gid)-south+1
       if(ix<1.or.ix>size(grid,1).or.iy<1.or.iy>size(grid,2))then
          ierr=3
          return
       endif
       buffer(:,i)=grid(ix,iy,:)
    enddo
  end subroutine HYMAP3_catchment_pack_grid

  subroutine HYMAP3_catchment_unpack_vector(plan,buffer,vector,ierr)
    type(HYMAP3_catchment_transfer_dec), intent(in) :: plan
    real, intent(in) :: buffer(:,:)
    real, intent(out) :: vector(:,:)
    integer, intent(out) :: ierr
    integer :: i,k
    integer, allocatable :: seen(:)
    ierr=0
    if(size(buffer,1)/=size(vector,1).or.size(buffer,2)/=size(plan%recv_gid).or. &
         size(vector,2)/=size(plan%recv_gid))then
       ierr=1
       return
    endif
    allocate(seen(size(vector,2)))
    seen=0
    do i=1,size(plan%recv_gid)
       k=plan%recv_local(i)
       if(k<1.or.k>size(vector,2))then
          ierr=2
          return
       endif
       seen(k)=seen(k)+1
       vector(:,k)=buffer(:,i)
    enddo
    if(any(seen/=1))ierr=3
  end subroutine HYMAP3_catchment_unpack_vector

  ! Reverse transfer for LIS-grid history or two-way LSM coupling.
  subroutine HYMAP3_catchment_pack_vector(plan,vector,buffer,ierr)
    type(HYMAP3_catchment_transfer_dec), intent(in) :: plan
    real, intent(in) :: vector(:,:)
    real, intent(out) :: buffer(:,:)
    integer, intent(out) :: ierr
    integer :: i,k
    ierr=0
    if(size(buffer,1)/=size(vector,1).or.size(buffer,2)/=size(plan%recv_gid).or. &
         size(vector,2)/=size(plan%recv_gid))then
       ierr=1
       return
    endif
    do i=1,size(plan%recv_gid)
       k=plan%recv_local(i)
       if(k<1.or.k>size(vector,2))then
          ierr=2
          return
       endif
       buffer(:,i)=vector(:,k)
    enddo
  end subroutine HYMAP3_catchment_pack_vector

  subroutine HYMAP3_catchment_unpack_grid(plan,gx,gy,west,south,buffer,grid,imis,ierr)
    type(HYMAP3_catchment_transfer_dec), intent(in) :: plan
    integer, intent(in) :: gx(:),gy(:),west,south
    real, intent(in) :: buffer(:,:),imis
    real, intent(out) :: grid(:,:,:)
    integer, intent(out) :: ierr
    integer :: i,gid,ix,iy
    ierr=0
    grid=imis
    if(size(buffer,1)/=size(grid,3).or.size(buffer,2)/=size(plan%send_gid))then
       ierr=1
       return
    endif
    do i=1,size(plan%send_gid)
       gid=plan%send_gid(i)
       if(gid<1.or.gid>size(gx).or.gid>size(gy))then
          ierr=2
          return
       endif
       ix=gx(gid)-west+1
       iy=gy(gid)-south+1
       if(ix<1.or.ix>size(grid,1).or.iy<1.or.iy>size(grid,2))then
          ierr=3
          return
       endif
       grid(ix,iy,:)=buffer(:,i)
    enddo
  end subroutine HYMAP3_catchment_unpack_grid
end module HYMAP3_catchmentMapMod
