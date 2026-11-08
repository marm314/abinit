!!****m* ABINIT/m_rdmft_driver
!! NAME
!!  m_rdmft_driver
!!
!! FUNCTION
!!  Driver for Reduced Density Matrix Functional Theory (RDMFT) calculations (optdriver=11).
!!
!! COPYRIGHT
!!  Copyright (C) 2026 ABINIT group (MRM)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rdmft_driver

 use defs_basis
 use m_errors
 use m_xmpi
 use m_dtfil,            only : datafiles_type
 use m_dtset,            only : dataset_type
 use m_fstrings,         only : sjoin, ftoa
 use m_nctk,             only : nctk_try_fort_or_ncfile
 use m_rdmft,            only : rdmft_hello
 use m_rdmft_orbitals,   only : rdmft_orbs_t
 use m_specialmsg,       only : wrtout

 implicit none

 private
!!***

 public :: rdmft_driver
!!***

contains
!!***

!!****f* m_rdmft_driver/rdmft_driver
!! NAME
!!  rdmft_driver
!!
!! FUNCTION
!!  Primary routine for RDMFT calculations.
!!  1) Read the ground-state (DFT or HF) orbitals from the WFK file (getwfk/irdwfk).
!!  2) (to come) optimize the occupation numbers for the fixed orbitals.
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type(datafiles_type)> = infos about file names
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!
!! OUTPUT
!!  (only printing)
!!
!! SOURCE
subroutine rdmft_driver(codvsn,dtfil,dtset)

 !Arguments ------------------------------------
 !scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset

 !Local variables-------------------------------
 !scalars
 integer,parameter :: master=0
 integer :: ierr
 real(dp) :: maxdev
 character(len=500) :: msg
 character(len=fnlen) :: wfk_path
 type(rdmft_orbs_t) :: orbs

! *************************************************************************

 ABI_UNUSED(dtset%natom)

 write(msg,'(3a)') ch10,' RDMFT driver (code version '//trim(codvsn)//')',ch10
 call wrtout([ab_out,std_out],msg)

 call rdmft_hello()

 !Find the WFK file of the ground-state calculation (Fortran or netcdf format)
 wfk_path = dtfil%fnamewffk
 if (xmpi_comm_rank(xmpi_world) == master) then
   if (nctk_try_fort_or_ncfile(wfk_path,msg) /= 0) then
     ABI_ERROR(sjoin('RDMFT: cannot find the ground-state WFK file:',ch10,msg))
   end if
 end if
 call xmpi_bcast(wfk_path,master,xmpi_world,ierr)

 !Read the orbitals
 call orbs%read_wfk(wfk_path)

 maxdev = orbs%check_orthonormality()
 write(msg,'(2a,es12.4,a)') ch10,' RDMFT: maximum deviation from orthonormality of the orbitals = ',maxdev,ch10
 call wrtout([ab_out,std_out],msg)

 call orbs%free()

end subroutine rdmft_driver
!!***

end module m_rdmft_driver
!!***
