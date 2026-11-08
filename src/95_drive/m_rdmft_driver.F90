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
 use m_rdmft,      only : rdmft_hello
 use m_specialmsg, only : wrtout

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
!!  For the moment it only calls a placeholder that prints a message.
!!
!! INPUTS
!!  codvsn = code version
!!
!! OUTPUT
!!  (only printing)
!!
!! SOURCE
subroutine rdmft_driver(codvsn)

 !Arguments ------------------------------------
 !scalars
 character(len=8), intent(in) :: codvsn

 !Local variables-------------------------------
 !scalars
 character(len=500) :: msg

! *************************************************************************

 write(msg,'(3a)') ch10,' RDMFT driver (code version ',trim(codvsn)//')'
 call wrtout(ab_out,msg)
 call wrtout(std_out,msg)

 call rdmft_hello()

end subroutine rdmft_driver
!!***

end module m_rdmft_driver
!!***
