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

module HYMAP3_netcdf_utilsMod

  use LIS_coreMod
  use LIS_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  private

  public :: HYMAP3_nc_get
  public :: HYMAP3_nc_get_dimlen
  public :: HYMAP3_nc_var_exists
  public :: HYMAP3_nc_dim_exists
  public :: HYMAP3_nc_get_gatt
  public :: HYMAP3_nc_get_vatt

  !----------------------------
  ! Generic variable reader API
  !----------------------------
  interface HYMAP3_nc_get
     ! integer
     module procedure HYMAP3_nc_get_var_int_0d
     module procedure HYMAP3_nc_get_var_int_1d
     module procedure HYMAP3_nc_get_var_int_2d
     module procedure HYMAP3_nc_get_var_int_3d

     module procedure HYMAP3_nc_get_var_real8_0d
     module procedure HYMAP3_nc_get_var_real8_1d
     module procedure HYMAP3_nc_get_var_real8_2d
     module procedure HYMAP3_nc_get_var_real8_3d
     module procedure HYMAP3_nc_get_var_real8_1d_from_2d
  end interface HYMAP3_nc_get

  !----------------------------
  ! Generic attribute reader API
  !   HYMAP3_nc_get_gatt(file, attname, val)
  !   HYMAP3_nc_get_vatt(file, varname, attname, val)
  !----------------------------
  interface HYMAP3_nc_get_gatt
     module procedure HYMAP3_nc_get_gatt_int
     module procedure HYMAP3_nc_get_gatt_real8
     module procedure HYMAP3_nc_get_gatt_char
  end interface HYMAP3_nc_get_gatt

  interface HYMAP3_nc_get_vatt
     module procedure HYMAP3_nc_get_vatt_int
     module procedure HYMAP3_nc_get_vatt_real8
     module procedure HYMAP3_nc_get_vatt_char
  end interface HYMAP3_nc_get_vatt

contains

  !===========================================================
  ! Low-level open/close (LIS-style fatal on missing file)
  !===========================================================
  subroutine HYMAP3_nc_open_read(ncfile, ncid)
    implicit none
    character(*), intent(in)  :: ncfile
    integer,      intent(out) :: ncid
    logical :: ex

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    inquire(file=trim(ncfile), exist=ex)
    if (.not. ex) then
       write(LIS_logunit,*) '[ERR] NetCDF input file not found: ', trim(ncfile)
       write(LIS_logunit,*) '[ERR] failed in HYMAP3_nc_open_read in HYMAP3_netcdf_utilsMod'
       call LIS_endrun()
    endif

    call LIS_verify(nf90_open(path=trim(ncfile), mode=NF90_NOWRITE, ncid=ncid), &
         'nf90_open failed in HYMAP3_nc_open_read in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled in this build.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_open_read


  subroutine HYMAP3_nc_close(ncid, where)
    implicit none
    integer,      intent(in) :: ncid
    character(*), intent(in) :: where
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call LIS_verify(nf90_close(ncid), 'nf90_close failed in '//trim(where))
#else
    ! unreachable because open_read endruns
#endif
  end subroutine HYMAP3_nc_close


  !===========================================================
  ! Dimension length (NOT part of HYMAP3_nc_get generic)
  !===========================================================
  subroutine HYMAP3_nc_get_dimlen(ncfile, dimname, dimlen)
    implicit none
    character(*), intent(in)  :: ncfile
    character(*), intent(in)  :: dimname
    integer,      intent(out) :: dimlen
    integer :: ncid, dimid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    dimlen = -999
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_dimid(ncid, trim(dimname), dimid), &
         'nf90_inq_dimid failed for dim '//trim(dimname)//' in HYMAP3_nc_get_dimlen in HYMAP3_netcdf_utilsMod')

    call LIS_verify(nf90_inquire_dimension(ncid, dimid, len=dimlen), &
         'nf90_inquire_dimension failed for dim '//trim(dimname)//' in HYMAP3_nc_get_dimlen in HYMAP3_netcdf_utilsMod')

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_dimlen in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_dimlen.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' dim=', trim(dimname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_dimlen


  !===========================================================
  ! Variable existence check
  !===========================================================
  subroutine HYMAP3_nc_var_exists(ncfile, varname, exists)

    implicit none

    character(*), intent(in)  :: ncfile
    character(*), intent(in)  :: varname
    logical,      intent(out) :: exists

    integer :: ncid
    integer :: varid
    integer :: status

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    exists = .false.

    call HYMAP3_nc_open_read(ncfile, ncid)

    status = nf90_inq_varid(ncid, trim(varname), varid)

    if (status == nf90_noerr) then
       exists = .true.
    endif

    call HYMAP3_nc_close(ncid, &
         'HYMAP3_nc_var_exists in HYMAP3_netcdf_utilsMod')

#else

    exists = .false.

    write(LIS_logunit,*) &
         '[ERR] NetCDF support not enabled for HYMAP3_nc_var_exists.'
    call LIS_endrun()

#endif

  end subroutine HYMAP3_nc_var_exists


  !===========================================================
  ! Dimension existence check
  !===========================================================
  subroutine HYMAP3_nc_dim_exists(ncfile, dimname, exists)

    implicit none

    character(*), intent(in)  :: ncfile
    character(*), intent(in)  :: dimname
    logical,      intent(out) :: exists

    integer :: ftn
    integer :: dimid
    integer :: status

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    exists = .false.

    call HYMAP3_nc_open_read(ncfile, ftn)

    status = nf90_inq_dimid(ftn, trim(dimname), dimid)

    if (status == nf90_noerr) then
      exists = .true.
    endif

    call HYMAP3_nc_close(ftn, &
         'HYMAP3_nc_dim_exists in HYMAP3_netcdf_utilsMod')

#else

    exists = .false.

#endif

  end subroutine HYMAP3_nc_dim_exists


  !===========================================================
  ! INTEGER variable readers
  !===========================================================
  subroutine HYMAP3_nc_get_var_int_0d(ncfile, varname, val, start, count)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    integer,      intent(out)          :: val
    integer,      intent(in), optional :: start, count
    integer :: ncid, varid
    integer :: s(1), c(1)
    integer :: tmp(1)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = 0
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_0d in HYMAP3_netcdf_utilsMod')

    if (present(start) .and. present(count)) then
       s(1) = start
       c(1) = count
       call LIS_verify(nf90_get_var(ncid, varid, tmp, start=s, count=c), &
            'nf90_get_var failed for var '//trim(varname)//' (int scalar slice) in HYMAP3_nc_get_var_int_0d in HYMAP3_netcdf_utilsMod')
       val = tmp(1)
    else
       call LIS_verify(nf90_get_var(ncid, varid, val), &
            'nf90_get_var failed for var '//trim(varname)//' (int scalar) in HYMAP3_nc_get_var_int_0d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_int_0d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_int_0d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_int_0d


  subroutine HYMAP3_nc_get_var_int_1d(ncfile, varname, arr, start1, count1)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    integer,      intent(inout)        :: arr(:)
    integer,      intent(in), optional :: start1, count1
    integer :: ncid, varid
    integer :: s(1), c(1)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_1d in HYMAP3_netcdf_utilsMod')

    if (present(start1) .and. present(count1)) then
       s(1) = start1
       c(1) = count1
       call LIS_verify(nf90_get_var(ncid, varid, arr, start=s, count=c), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_1d in HYMAP3_netcdf_utilsMod')
    else
       call LIS_verify(nf90_get_var(ncid, varid, arr), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_1d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_int_1d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_int_1d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_int_1d


  subroutine HYMAP3_nc_get_var_int_2d(ncfile, varname, arr, start, count)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    integer,      intent(inout)        :: arr(:,:)
    integer,      intent(in), optional :: start(2), count(2)
    integer :: ncid, varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_2d in HYMAP3_netcdf_utilsMod')

    if (present(start) .and. present(count)) then
       call LIS_verify(nf90_get_var(ncid, varid, arr, start=start, count=count), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_2d in HYMAP3_netcdf_utilsMod')
    else
       call LIS_verify(nf90_get_var(ncid, varid, arr), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_2d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_int_2d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_int_2d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_int_2d


  subroutine HYMAP3_nc_get_var_int_3d(ncfile, varname, arr, start, count)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    integer,      intent(inout)        :: arr(:,:,:)
    integer,      intent(in), optional :: start(3), count(3)
    integer :: ncid, varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_3d in HYMAP3_netcdf_utilsMod')

    if (present(start) .and. present(count)) then
       call LIS_verify(nf90_get_var(ncid, varid, arr, start=start, count=count), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_3d in HYMAP3_netcdf_utilsMod')
    else
       call LIS_verify(nf90_get_var(ncid, varid, arr), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_int_3d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_int_3d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_int_3d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_int_3d


  !===========================================================
  !===========================================================
  subroutine HYMAP3_nc_get_var_real8_0d(ncfile, varname, val, start, count)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    real*8,       intent(out)          :: val
    integer,      intent(in), optional :: start, count
    integer :: ncid, varid
    integer :: s(1), c(1)
    real*8  :: tmp(1)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = 0.0d0
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_0d in HYMAP3_netcdf_utilsMod')

    if (present(start) .and. present(count)) then
       ! IMPORTANT: scalar-slice read via 1-element array (portable)
       s(1) = start
       c(1) = count
       call LIS_verify(nf90_get_var(ncid, varid, tmp, start=s, count=c), &
            'nf90_get_var failed for var '//trim(varname)//' (real*8 scalar slice) in HYMAP3_nc_get_var_real8_0d in HYMAP3_netcdf_utilsMod')
       val = tmp(1)
    else
       call LIS_verify(nf90_get_var(ncid, varid, val), &
            'nf90_get_var failed for var '//trim(varname)//' (real*8 scalar) in HYMAP3_nc_get_var_real8_0d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_real8_0d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_real8_0d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_real8_0d


  subroutine HYMAP3_nc_get_var_real8_1d(ncfile, varname, arr, start1, count1)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    real*8,       intent(inout)        :: arr(:)
    integer,      intent(in), optional :: start1, count1
    integer :: ncid, varid
    integer :: s(1), c(1)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_1d in HYMAP3_netcdf_utilsMod')

    if (present(start1) .and. present(count1)) then
       s(1) = start1
       c(1) = count1
       call LIS_verify(nf90_get_var(ncid, varid, arr, start=s, count=c), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_1d in HYMAP3_netcdf_utilsMod')
    else
       call LIS_verify(nf90_get_var(ncid, varid, arr), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_1d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_real8_1d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_real8_1d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_real8_1d


  subroutine HYMAP3_nc_get_var_real8_2d(ncfile, varname, arr, start, count)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    real*8,       intent(inout)        :: arr(:,:)
    integer,      intent(in), optional :: start(2), count(2)
    integer :: ncid, varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_2d in HYMAP3_netcdf_utilsMod')

    if (present(start) .and. present(count)) then
       call LIS_verify(nf90_get_var(ncid, varid, arr, start=start, count=count), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_2d in HYMAP3_netcdf_utilsMod')
    else
       call LIS_verify(nf90_get_var(ncid, varid, arr), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_2d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_real8_2d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_real8_2d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_real8_2d


  subroutine HYMAP3_nc_get_var_real8_3d(ncfile, varname, arr, start, count)
    implicit none
    character(*), intent(in)           :: ncfile, varname
    real*8,       intent(inout)        :: arr(:,:,:)
    integer,      intent(in), optional :: start(3), count(3)
    integer :: ncid, varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_3d in HYMAP3_netcdf_utilsMod')

    if (present(start) .and. present(count)) then
       call LIS_verify(nf90_get_var(ncid, varid, arr, start=start, count=count), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_3d in HYMAP3_netcdf_utilsMod')
    else
       call LIS_verify(nf90_get_var(ncid, varid, arr), &
            'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_3d in HYMAP3_netcdf_utilsMod')
    endif

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_real8_3d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_real8_3d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_real8_3d


  !-----------------------------------------------------------
  ! 2-D -> 1-D row-slice (mandatory start/count rank=2)
  !   Used by calls like:
  !     call HYMAP3_nc_get(file,'mon_min_sto', arr(:), start=(/row,1/), count=(/1,12/))
  !-----------------------------------------------------------
  subroutine HYMAP3_nc_get_var_real8_1d_from_2d(ncfile, varname, arr, start, count)
    implicit none
    character(*), intent(in)    :: ncfile, varname
    real*8,       intent(inout) :: arr(:)
    integer,      intent(in)    :: start(2), count(2)
    integer :: ncid, varid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    call HYMAP3_nc_open_read(ncfile, ncid)

    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_1d_from_2d in HYMAP3_netcdf_utilsMod')

    call LIS_verify(nf90_get_var(ncid, varid, arr, start=start, count=count), &
         'nf90_get_var failed for var '//trim(varname)//' in HYMAP3_nc_get_var_real8_1d_from_2d in HYMAP3_netcdf_utilsMod')

    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_var_real8_1d_from_2d in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_var_real8_1d_from_2d.'
    write(LIS_logunit,*) '[ERR] file=', trim(ncfile), ' var=', trim(varname)
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_var_real8_1d_from_2d


  !===========================================================
  ! Attribute readers (global + variable)
  !===========================================================
  subroutine HYMAP3_nc_get_gatt_int(ncfile, attname, val)
    implicit none
    character(*), intent(in) :: ncfile, attname
    integer,      intent(out):: val
    integer :: ncid
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = 0
    call HYMAP3_nc_open_read(ncfile, ncid)
    call LIS_verify(nf90_get_att(ncid, NF90_GLOBAL, trim(attname), val), &
         'nf90_get_att failed for global att '//trim(attname)//' in HYMAP3_nc_get_gatt_int in HYMAP3_netcdf_utilsMod')
    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_gatt_int in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_gatt_int.'
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_gatt_int


  subroutine HYMAP3_nc_get_gatt_real8(ncfile, attname, val)
    implicit none
    character(*), intent(in) :: ncfile, attname
    real*8,       intent(out):: val
    integer :: ncid
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = 0.0d0
    call HYMAP3_nc_open_read(ncfile, ncid)
    call LIS_verify(nf90_get_att(ncid, NF90_GLOBAL, trim(attname), val), &
         'nf90_get_att failed for global att '//trim(attname)//' in HYMAP3_nc_get_gatt_real8 in HYMAP3_netcdf_utilsMod')
    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_gatt_real8 in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_gatt_real8.'
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_gatt_real8


  subroutine HYMAP3_nc_get_gatt_char(ncfile, attname, val)
    implicit none
    character(*), intent(in)  :: ncfile, attname
    character(*), intent(out) :: val
    integer :: ncid
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = ''
    call HYMAP3_nc_open_read(ncfile, ncid)
    call LIS_verify(nf90_get_att(ncid, NF90_GLOBAL, trim(attname), val), &
         'nf90_get_att failed for global att '//trim(attname)//' in HYMAP3_nc_get_gatt_char in HYMAP3_netcdf_utilsMod')
    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_gatt_char in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_gatt_char.'
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_gatt_char


  subroutine HYMAP3_nc_get_vatt_int(ncfile, varname, attname, val)
    implicit none
    character(*), intent(in) :: ncfile, varname, attname
    integer,      intent(out):: val
    integer :: ncid, varid
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = 0
    call HYMAP3_nc_open_read(ncfile, ncid)
    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_vatt_int in HYMAP3_netcdf_utilsMod')
    call LIS_verify(nf90_get_att(ncid, varid, trim(attname), val), &
         'nf90_get_att failed for var att '//trim(attname)//' on '//trim(varname)//' in HYMAP3_nc_get_vatt_int in HYMAP3_netcdf_utilsMod')
    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_vatt_int in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_vatt_int.'
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_vatt_int


  subroutine HYMAP3_nc_get_vatt_real8(ncfile, varname, attname, val)
    implicit none
    character(*), intent(in) :: ncfile, varname, attname
    real*8,       intent(out):: val
    integer :: ncid, varid
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = 0.0d0
    call HYMAP3_nc_open_read(ncfile, ncid)
    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_vatt_real8 in HYMAP3_netcdf_utilsMod')
    call LIS_verify(nf90_get_att(ncid, varid, trim(attname), val), &
         'nf90_get_att failed for var att '//trim(attname)//' on '//trim(varname)//' in HYMAP3_nc_get_vatt_real8 in HYMAP3_netcdf_utilsMod')
    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_vatt_real8 in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_vatt_real8.'
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_vatt_real8


  subroutine HYMAP3_nc_get_vatt_char(ncfile, varname, attname, val)
    implicit none
    character(*), intent(in)  :: ncfile, varname, attname
    character(*), intent(out) :: val
    integer :: ncid, varid
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    val = ''
    call HYMAP3_nc_open_read(ncfile, ncid)
    call LIS_verify(nf90_inq_varid(ncid, trim(varname), varid), &
         'nf90_inq_varid failed for var '//trim(varname)//' in HYMAP3_nc_get_vatt_char in HYMAP3_netcdf_utilsMod')
    call LIS_verify(nf90_get_att(ncid, varid, trim(attname), val), &
         'nf90_get_att failed for var att '//trim(attname)//' on '//trim(varname)//' in HYMAP3_nc_get_vatt_char in HYMAP3_netcdf_utilsMod')
    call HYMAP3_nc_close(ncid, 'HYMAP3_nc_get_vatt_char in HYMAP3_netcdf_utilsMod')
#else
    write(LIS_logunit,*) '[ERR] NetCDF support not enabled for HYMAP3_nc_get_vatt_char.'
    call LIS_endrun()
#endif
  end subroutine HYMAP3_nc_get_vatt_char

end module HYMAP3_netcdf_utilsMod
