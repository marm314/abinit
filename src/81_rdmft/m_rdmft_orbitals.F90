!!****m* ABINIT/m_rdmft_orbitals
!! NAME
!!  m_rdmft_orbitals
!!
!! FUNCTION
!!  Container for the Bloch orbitals (plane-wave coefficients, eigenvalues and occupations)
!!  of a ground-state (DFT or Hartree-Fock) calculation, read from a WFK file.
!!  The orbitals are kept fixed during the RDMFT optimization of the occupation numbers.
!!  Present limitations: norm-conserving pseudopotentials, nspinor=1, nsppol=1, istwfk=1 at all k-points.
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

module m_rdmft_orbitals

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_hdr
 use m_wfk

 use m_io_tools,   only : iomode_from_fname, get_unit
 use m_specialmsg, only : wrtout

 implicit none

 private
!!***

!!****t* m_rdmft_orbitals/rdmft_orbs_t
!! NAME
!! rdmft_orbs_t
!!
!! FUNCTION
!!  Ground-state Bloch orbitals read from a WFK file.
!!
!! SOURCE

 type,public :: rdmft_orbs_t

   integer :: nkpt = 0
   ! Number of k-points.

   integer :: nsppol = 0
   ! Number of spin polarizations.

   integer :: nspinor = 0
   ! Number of spinor components.

   integer :: mband = 0
   ! Maximum number of bands over k-points.

   integer :: mpw = 0
   ! Maximum number of plane waves over k-points.

   real(dp) :: nelect = zero
   ! Number of electrons.

   integer,allocatable :: nband(:,:)
   ! nband(nkpt,nsppol)
   ! Number of bands at each k-point and spin.

   integer,allocatable :: npw(:)
   ! npw(nkpt)
   ! Number of plane waves at each k-point.

   integer,allocatable :: kg(:,:,:)
   ! kg(3,mpw,nkpt)
   ! Reduced coordinates of the plane waves at each k-point.

   real(dp),allocatable :: kpt(:,:)
   ! kpt(3,nkpt)
   ! Reduced coordinates of the k-points.

   real(dp),allocatable :: wtk(:)
   ! wtk(nkpt)
   ! Weights of the k-points (normalized to one).

   real(dp),allocatable :: eig(:,:,:)
   ! eig(mband,nkpt,nsppol)
   ! Eigenvalues (Hartree).

   real(dp),allocatable :: occ(:,:,:)
   ! occ(mband,nkpt,nsppol)
   ! Occupation numbers of the ground-state calculation (between 0 and 2/nsppol/nspinor).

   real(dp),allocatable :: cg(:,:,:,:,:)
   ! cg(2,mpw*nspinor,mband,nkpt,nsppol)
   ! Plane-wave coefficients, normalized to one for each band.

 contains

   procedure :: free => rdmft_orbs_free
   ! Free memory.

   procedure :: read_wfk => rdmft_orbs_read_wfk
   ! Read the orbitals from a WFK file.

   procedure :: check_orthonormality => rdmft_orbs_check_orthonormality
   ! Compute the maximum deviation from orthonormality.

 end type rdmft_orbs_t
!!***

contains
!!***

!!****f* m_rdmft_orbitals/rdmft_orbs_free
!! NAME
!!  rdmft_orbs_free
!!
!! FUNCTION
!!  Free the memory of the orbitals.
!!
!! SOURCE
subroutine rdmft_orbs_free(orbs)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_orbs_t),intent(inout) :: orbs

! *************************************************************************

 ABI_SFREE(orbs%nband)
 ABI_SFREE(orbs%npw)
 ABI_SFREE(orbs%kg)
 ABI_SFREE(orbs%kpt)
 ABI_SFREE(orbs%wtk)
 ABI_SFREE(orbs%eig)
 ABI_SFREE(orbs%occ)
 ABI_SFREE(orbs%cg)

 orbs%nkpt = 0
 orbs%nsppol = 0
 orbs%nspinor = 0
 orbs%mband = 0
 orbs%mpw = 0
 orbs%nelect = zero

end subroutine rdmft_orbs_free
!!***

!!****f* m_rdmft_orbitals/rdmft_orbs_read_wfk
!! NAME
!!  rdmft_orbs_read_wfk
!!
!! FUNCTION
!!  Read the orbitals, eigenvalues and occupations of all k-points from a WFK file.
!!
!! INPUTS
!!  wfk_path = path of the WFK file (Fortran or netcdf format)
!!
!! OUTPUT
!!  orbs <type(rdmft_orbs_t)> = the orbitals (allocated here)
!!
!! SOURCE
subroutine rdmft_orbs_read_wfk(orbs,wfk_path)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_orbs_t),intent(inout) :: orbs
 character(len=*),intent(in) :: wfk_path

 !Local variables-------------------------------
 !scalars
 integer,parameter :: formeig0=0
 integer :: ik,isppol,nb,npw_k,sc_mode
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wfk_t) :: wfk
 !arrays
 integer :: band_block(2)
 integer,allocatable :: kg_k(:,:)
 real(dp),allocatable :: cg_k(:,:),eig_k(:),occ_k(:)

! *************************************************************************

 call orbs%free()

 call wfk%open_read(wfk_path,formeig0,iomode_from_fname(wfk_path),get_unit(),xmpi_comm_self,hdr_out=hdr)

 if (hdr%usepaw /= 0) then
   ABI_ERROR('RDMFT: PAW is not yet supported.')
 end if
 if (hdr%nspinor /= 1) then
   ABI_ERROR('RDMFT: nspinor=1 is required.')
 end if
 if (hdr%nsppol /= 1) then
   ABI_ERROR('RDMFT: nsppol=1 is required (closed-shell).')
 end if
 if (any(hdr%istwfk(:) /= 1)) then
   ABI_ERROR('RDMFT: istwfk=1 is required at all k-points (use istwfk *1 in the ground-state calculation).')
 end if

 orbs%nkpt = hdr%nkpt
 orbs%nsppol = hdr%nsppol
 orbs%nspinor = hdr%nspinor
 orbs%mband = maxval(hdr%nband)
 orbs%mpw = maxval(hdr%npwarr)
 orbs%nelect = hdr%nelect

 ABI_MALLOC(orbs%nband,(orbs%nkpt,orbs%nsppol))
 ABI_MALLOC(orbs%npw,(orbs%nkpt))
 ABI_MALLOC(orbs%kg,(3,orbs%mpw,orbs%nkpt))
 ABI_MALLOC(orbs%kpt,(3,orbs%nkpt))
 ABI_MALLOC(orbs%wtk,(orbs%nkpt))
 ABI_MALLOC(orbs%eig,(orbs%mband,orbs%nkpt,orbs%nsppol))
 ABI_MALLOC(orbs%occ,(orbs%mband,orbs%nkpt,orbs%nsppol))
 ABI_MALLOC(orbs%cg,(2,orbs%mpw*orbs%nspinor,orbs%mband,orbs%nkpt,orbs%nsppol))

 orbs%nband = 0
 orbs%kg = 0
 orbs%eig = zero
 orbs%occ = zero
 orbs%cg = zero
 orbs%npw(:) = hdr%npwarr(:)
 orbs%kpt(:,:) = hdr%kptns(:,:)
 orbs%wtk(:) = hdr%wtk(:)/sum(hdr%wtk)

 sc_mode = xmpio_single
 do isppol=1,orbs%nsppol
   do ik=1,orbs%nkpt
     nb = hdr%nband(ik+(isppol-1)*orbs%nkpt)
     orbs%nband(ik,isppol) = nb
     npw_k = orbs%npw(ik)
     band_block = [1,nb]

     ABI_MALLOC(kg_k,(3,npw_k))
     ABI_MALLOC(cg_k,(2,npw_k*orbs%nspinor*nb))
     ABI_MALLOC(eig_k,(wfk%mband))
     ABI_MALLOC(occ_k,(wfk%mband))

     call wfk%read_band_block(band_block,ik,isppol,sc_mode,kg_k=kg_k,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)

     orbs%kg(:,1:npw_k,ik) = kg_k(:,:)
     orbs%cg(:,1:npw_k*orbs%nspinor,1:nb,ik,isppol) = reshape(cg_k,[2,npw_k*orbs%nspinor,nb])
     orbs%eig(1:nb,ik,isppol) = eig_k(1:nb)
     orbs%occ(1:nb,ik,isppol) = occ_k(1:nb)

     ABI_FREE(kg_k)
     ABI_FREE(cg_k)
     ABI_FREE(eig_k)
     ABI_FREE(occ_k)
   end do
 end do

 call wfk%close()
 call hdr%free()

 write(msg,'(3a,i0,a,i0,a,i0,a,i0,a)') ' RDMFT: orbitals read from ',trim(wfk_path),ch10,&
   orbs%nkpt,' k-points, ',orbs%mband,' bands, up to ',orbs%mpw,' plane waves, nelect = ',nint(orbs%nelect),ch10
 call wrtout([ab_out,std_out],msg)

end subroutine rdmft_orbs_read_wfk
!!***

!!****f* m_rdmft_orbitals/rdmft_orbs_check_orthonormality
!! NAME
!!  rdmft_orbs_check_orthonormality
!!
!! FUNCTION
!!  Compute the maximum over k-points of |<psi_nk|psi_mk> - delta_nm|.
!!
!! OUTPUT
!!  maxdev = maximum deviation from orthonormality
!!
!! SOURCE
function rdmft_orbs_check_orthonormality(orbs) result(maxdev)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_orbs_t),intent(in) :: orbs
 real(dp) :: maxdev

 !Local variables-------------------------------
 !scalars
 integer :: ib,ik,isppol,jb,nb,npwso
 complex(dpc) :: ovlp

! *************************************************************************

 maxdev = zero
 do isppol=1,orbs%nsppol
   do ik=1,orbs%nkpt
     nb = orbs%nband(ik,isppol)
     npwso = orbs%npw(ik)*orbs%nspinor
     do ib=1,nb
       do jb=1,ib
         ovlp = sum(cmplx(orbs%cg(1,1:npwso,ib,ik,isppol),-orbs%cg(2,1:npwso,ib,ik,isppol),kind=dpc)* &
                    cmplx(orbs%cg(1,1:npwso,jb,ik,isppol), orbs%cg(2,1:npwso,jb,ik,isppol),kind=dpc))
         if (ib == jb) ovlp = ovlp - cone
         maxdev = max(maxdev,abs(ovlp))
       end do
     end do
   end do
 end do

end function rdmft_orbs_check_orthonormality
!!***

end module m_rdmft_orbitals
!!***
