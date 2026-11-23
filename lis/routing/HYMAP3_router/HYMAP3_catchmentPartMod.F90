!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP3_catchmentPartMod
  implicit none
  private
  public :: HYMAP3_partition_catchments

contains

  !ag(22Nov2026)
  ! Partition an existing 1D HyMAP river network, independently of its
  ! original 2D grid. A link next(i)>0 flows from i to next(i).
  ! next(i)=0, -9 or -9999 denotes an outlet; there may be many outlets.
  !
  ! The connected components are built bottom-up and never exceed the
  ! target number of river cells. Whole components are assigned to the
  ! least-loaded processor (largest first). No geometry is modified.
  ! Graph construction O(N), processor assignment O(N log N + N log P).
  ! IMPORTANT: this routine computes an ownership PLAN only; it does not
  ! redistribute LIS fields, routing state, or ESMF routing tiles.
  subroutine HYMAP3_partition_catchments(next, np, owner, load, ncut, ierr, message)
    integer, intent(in) :: next(:), np
    integer, intent(out) :: owner(:), load(:), ncut, ierr
    character(*), intent(out) :: message
    integer, allocatable :: indeg(:), queue(:), parent(:), csize(:)
    integer, allocatable :: component(:), proc_heap(:)
    integer :: n, i, j, k, first, last, processed, target, ncomp
    integer :: root, rank

    ierr=0
    message='OK'
    ncut=0
    owner=-1
    load=0
    n=size(next)
    if(n<1.or.np<1.or.size(owner)/=n.or.size(load)/=np)then
       ierr=1
       message='Invalid array sizes or processor count'
       return
    endif
    if(np>n)then
       ierr=2
       message='More processors than river elements: empty MPI routing ranks'
       return
    endif

    allocate(indeg(n),queue(n),parent(n),csize(n),component(n),proc_heap(np))
    indeg=0
    do i=1,n
       j=next(i)
       if(j==0.or.j==-9.or.j==-9999)cycle
       if(j<1.or.j>n)then
          ierr=3
          message='Downstream index outside the river network'
          return
       endif
       if(j==i)then
          ierr=4
          message='Self-connected river element'
          return
       endif
       indeg(j)=indeg(j)+1
    enddo

    first=1
    last=0
    do i=1,n
       if(indeg(i)==0)then
          last=last+1
          queue(last)=i
       endif
    enddo
    processed=0
    do while(first<=last)
       i=queue(first)
       first=first+1
       processed=processed+1
       queue(processed)=i
       j=next(i)
       if(j>0)then
          indeg(j)=indeg(j)-1
          if(indeg(j)==0)then
             last=last+1
             queue(last)=j
          endif
       endif
    enddo
    if(processed/=n)then
       ierr=5
       message='River network contains a directed cycle'
       return
    endif

    ! A topological pass always visits tributaries before their recipient.
    ! Accepted links merge two connected pieces; rejected links become
    ! inter-component boundaries, never splitting the river continuity.
    target=(n+np-1)/np
    do i=1,n
       parent(i)=i
       csize(i)=1
    enddo
    do k=1,n
       i=queue(k)
       j=next(i)
       if(j<=0)cycle
       root=HYMAP3_find_root(i,parent)
       if(csize(j)+csize(root)<=target)then
          parent(root)=j
          csize(j)=csize(j)+csize(root)
       endif
    enddo

    ncomp=0
    do i=1,n
       if(parent(i)==i)then
          ncomp=ncomp+1
          component(ncomp)=i
       endif
    enddo
    ! Sort connected pieces largest first, with a deterministic ID tie-break.
    call HYMAP3_sort_components(component,ncomp,csize)
    do i=1,np
       proc_heap(i)=i-1
    enddo
    do k=ncomp,1,-1
       root=component(k)
       rank=proc_heap(1)
       owner(root)=rank
       load(rank+1)=load(rank+1)+csize(root)
       call HYMAP3_shift_processor(proc_heap,np,load)
    enddo
    do i=1,n
       root=HYMAP3_find_root(i,parent)
       owner(i)=owner(root)
    enddo
    if(any(owner<0).or.sum(load)/=n.or.any(load==0))then
       ierr=6
       message='Internal inconsistency in ownership or processor loads'
       return
    endif
    do i=1,n
       j=next(i)
       if(j>0)then
          if(owner(i)/=owner(j))ncut=ncut+1
       endif
    enddo
  end subroutine HYMAP3_partition_catchments

  integer function HYMAP3_find_root(i,parent) result(root)
    integer, intent(in) :: i
    integer, intent(inout) :: parent(:)
    integer :: j, k
    root=i
    do while(parent(root)/=root)
       root=parent(root)
    enddo
    j=i
    do while(parent(j)/=root)
       k=parent(j)
       parent(j)=root
       j=k
    enddo
  end function HYMAP3_find_root

  logical function HYMAP3_larger(a,b,csize) result(yes)
    integer, intent(in) :: a,b,csize(:)
    yes=csize(a)>csize(b)
    if(csize(a)==csize(b))yes=a<b
  end function HYMAP3_larger

  subroutine HYMAP3_sort_components(component,ncomp,csize)
    integer, intent(in) :: ncomp,csize(:)
    integer, intent(inout) :: component(:)
    integer :: i, last, tmp
    do i=ncomp/2,1,-1
       call sift(i,ncomp)
    enddo
    do last=ncomp,2,-1
       tmp=component(1)
       component(1)=component(last)
       component(last)=tmp
       call sift(1,last-1)
    enddo
  contains
    subroutine sift(start,finish)
      integer, intent(in) :: start,finish
      integer :: root, child, t
      root=start
      do while(root<=finish/2)
         child=2*root
         if(child<finish)then
            if(HYMAP3_larger(component(child+1),component(child),csize)) &
                 child=child+1
         endif
         if(.not.HYMAP3_larger(component(child),component(root),csize))exit
         t=component(root)
         component(root)=component(child)
         component(child)=t
         root=child
      enddo
    end subroutine sift
  end subroutine HYMAP3_sort_components

  subroutine HYMAP3_shift_processor(heap,np,load)
    integer, intent(in) :: np,load(:)
    integer, intent(inout) :: heap(:)
    integer :: root, child, t
    root=1
    do while(root<=np/2)
       child=2*root
       if(child<np)then
          if(load(heap(child+1)+1)<load(heap(child)+1).or. &
               (load(heap(child+1)+1)==load(heap(child)+1).and. &
               heap(child+1)<heap(child)))child=child+1
       endif
       if(load(heap(root)+1)<load(heap(child)+1))exit
       if(load(heap(root)+1)==load(heap(child)+1))then
          if(heap(root)<heap(child))exit
       endif
       t=heap(root)
       heap(root)=heap(child)
       heap(child)=t
       root=child
    enddo
  end subroutine HYMAP3_shift_processor
end module HYMAP3_catchmentPartMod
