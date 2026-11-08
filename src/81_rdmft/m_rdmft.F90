!!****m* ABINIT/m_rdmft
!! NAME
!!  m_rdmft
!!
!! FUNCTION
!!  Contains the routines of the Reduced Density Matrix Functional Theory (RDMFT) library.
!!  For the moment it only contains a placeholder that prints a message.
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

module m_rdmft

 use defs_basis
 use m_specialmsg, only : wrtout

 implicit none

 private
!!***

 public :: rdmft_hello
!!***

contains
!!***

!!****f* m_rdmft/rdmft_hello
!! NAME
!!  rdmft_hello
!!
!! FUNCTION
!!  Placeholder routine: print "Hello world!" in the main output and log files.
!!
!! INPUTS
!!  (none)
!!
!! OUTPUT
!!  (only printing)
!!
!! SOURCE
subroutine rdmft_hello()

 !Arguments ------------------------------------

 !Local variables-------------------------------
 !scalars
 character(len=500) :: msg

! *************************************************************************

 write(msg,'(3a)') ch10,' Hello world!',ch10
 call wrtout(ab_out,msg)
 call wrtout(std_out,msg)

end subroutine rdmft_hello
!!***

end module m_rdmft
!!***
