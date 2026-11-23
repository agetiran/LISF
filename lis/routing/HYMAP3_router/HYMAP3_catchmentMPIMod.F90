!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module HYMAP3_catchmentMPIMod
  !ag(22Nov2026) Reuse the existing network-MPI approach for LIS/HyMAP
  ! forcing and history exchanges. All peer lists are calculated ONCE
  ! by HYMAP3_catchment_build_transfer, in increasing global river ID.
  ! Data are moved without interpolating or modifying runoff rates.
  use HYMAP3_catchmentMapMod
  use LIS_mpiMod
  implicit none
  private
  public :: HYMAP3_catchment_forward, HYMAP3_catchment_reverse

contains

  subroutine HYMAP3_catchment_forward(plan,gx,gy,west,south,grid,vector,ierr)
    type(HYMAP3_catchment_transfer_dec), intent(inout) :: plan
    integer, intent(in) :: gx(:),gy(:),west,south
    real, intent(in) :: grid(:,:,:)
    real, intent(out) :: vector(:,:)
    integer, intent(out) :: ierr
    integer :: nfield

    ierr=0
    nfield=size(grid,3)
    if(nfield<1.or.size(vector,1)/=nfield.or. &
         size(vector,2)/=size(plan%recv_gid))then
       ierr=101
       return
    endif
    call HYMAP3_catchment_workspace(plan,nfield)
    call HYMAP3_catchment_pack_grid(plan,gx,gy,west,south,grid, &
         plan%buffer_out(:,:size(plan%send_gid)),ierr)
    if(ierr/=0)then
       ierr=110+ierr
       return
    endif
    call HYMAP3_catchment_exchange(plan%send_count,plan%send_offset, &
         plan%recv_count,plan%recv_offset, &
         plan%buffer_out(:,:size(plan%send_gid)), &
         plan%buffer_in(:,:size(plan%recv_gid)),ierr)
    if(ierr/=0)return
    call HYMAP3_catchment_unpack_vector(plan, &
         plan%buffer_in(:,:size(plan%recv_gid)),vector,ierr)
    if(ierr/=0)ierr=120+ierr
  end subroutine HYMAP3_catchment_forward

  subroutine HYMAP3_catchment_reverse(plan,gx,gy,west,south,vector,grid,imis,ierr)
    type(HYMAP3_catchment_transfer_dec), intent(inout) :: plan
    integer, intent(in) :: gx(:),gy(:),west,south
    real, intent(in) :: vector(:,:),imis
    real, intent(out) :: grid(:,:,:)
    integer, intent(out) :: ierr
    integer :: nfield

    ierr=0
    nfield=size(grid,3)
    if(nfield<1.or.size(vector,1)/=nfield.or. &
         size(vector,2)/=size(plan%recv_gid))then
       ierr=201
       return
    endif
    call HYMAP3_catchment_workspace(plan,nfield)
    call HYMAP3_catchment_pack_vector(plan,vector, &
         plan%buffer_out(:,:size(plan%recv_gid)),ierr)
    if(ierr/=0)then
       ierr=210+ierr
       return
    endif
    call HYMAP3_catchment_exchange(plan%recv_count,plan%recv_offset, &
         plan%send_count,plan%send_offset, &
         plan%buffer_out(:,:size(plan%recv_gid)), &
         plan%buffer_in(:,:size(plan%send_gid)),ierr)
    if(ierr/=0)return
    call HYMAP3_catchment_unpack_grid(plan,gx,gy,west,south, &
         plan%buffer_in(:,:size(plan%send_gid)),grid,imis,ierr)
    if(ierr/=0)ierr=220+ierr
  end subroutine HYMAP3_catchment_reverse

  subroutine HYMAP3_catchment_workspace(plan,nfield)
    type(HYMAP3_catchment_transfer_dec), intent(inout) :: plan
    integer, intent(in) :: nfield
    integer :: nbuf
    if(plan%nbuffer_fields==nfield)return
    if(allocated(plan%buffer_out))deallocate(plan%buffer_out)
    if(allocated(plan%buffer_in))deallocate(plan%buffer_in)
    nbuf=max(size(plan%send_gid),size(plan%recv_gid))
    allocate(plan%buffer_out(nfield,nbuf),plan%buffer_in(nfield,nbuf))
    plan%nbuffer_fields=nfield
  end subroutine HYMAP3_catchment_workspace

  subroutine HYMAP3_catchment_exchange(send_count,send_offset, &
       recv_count,recv_offset,sendbuf,recvbuf,ierr)
    ! Direct peer-to-peer communication, following HYMAP3_network_*.
    ! Peer segments are contiguous; MPI moves nfield consecutive REALs
    ! per cell. Self sends are legal, and empty peers are skipped.
    integer, intent(in) :: send_count(0:),send_offset(0:)
    integer, intent(in) :: recv_count(0:),recv_offset(0:)
    real, intent(in) :: sendbuf(:,:)
    real, intent(out) :: recvbuf(:,:)
    integer, intent(out) :: ierr
    integer :: nf,np,p
#if (defined SPMD)
    integer :: nreq,status
    integer, allocatable :: request(:),statuses(:,:)
    external :: MPI_IRECV, MPI_ISEND, MPI_WAITALL
#endif
    ierr=0
    nf=size(sendbuf,1)
    np=size(send_count)
    if(nf<1.or.size(recvbuf,1)/=nf.or.size(recv_count)/=np.or. &
         size(send_offset)/=np+1.or.size(recv_offset)/=np+1)then
       ierr=301
       return
    endif
    if(send_offset(0)/=0.or.recv_offset(0)/=0.or. &
         any(send_count<0).or.any(recv_count<0))then
       ierr=302
       return
    endif
    do p=0,np-1
       if(send_offset(p+1)-send_offset(p)/=send_count(p).or. &
            recv_offset(p+1)-recv_offset(p)/=recv_count(p))then
          ierr=303
          return
       endif
    enddo
    if(send_offset(np)/=size(sendbuf,2).or. &
         recv_offset(np)/=size(recvbuf,2))then
       ierr=304
       return
    endif
    if(max(size(sendbuf,2),size(recvbuf,2))>huge(nf)/nf)then
       ierr=305
       return
    endif
#if (defined SPMD)
    allocate(request(2*np),statuses(MPI_STATUS_SIZE,2*np))
    nreq=0
    do p=0,np-1
       if(recv_count(p)==0)cycle
       nreq=nreq+1
       call MPI_IRECV(recvbuf(1,recv_offset(p)+1),nf*recv_count(p), &
            MPI_REAL,p,22640,LIS_mpi_comm,request(nreq),status)
       if(status/=MPI_SUCCESS)then
          ierr=310
          return
       endif
    enddo
    do p=0,np-1
       if(send_count(p)==0)cycle
       nreq=nreq+1
       call MPI_ISEND(sendbuf(1,send_offset(p)+1),nf*send_count(p), &
            MPI_REAL,p,22640,LIS_mpi_comm,request(nreq),status)
       if(status/=MPI_SUCCESS)then
          ierr=311
          return
       endif
    enddo
    if(nreq>0)then
       call MPI_WAITALL(nreq,request,statuses,status)
       if(status/=MPI_SUCCESS)ierr=312
    endif
#else
    ! Valid non-SPMD execution has one rank and identical source/dest.
    if(np/=1.or.send_count(0)/=recv_count(0))then
       ierr=313
       return
    endif
    recvbuf=sendbuf
#endif
  end subroutine HYMAP3_catchment_exchange
end module HYMAP3_catchmentMPIMod
