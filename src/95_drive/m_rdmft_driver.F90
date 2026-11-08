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
 use defs_abitypes,        only : MPI_type
 use defs_datatypes,       only : pseudopotential_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_errors
 use m_xmpi
 use m_dtfil,            only : datafiles_type
 use m_dtset,            only : dataset_type
 use m_fstrings,         only : sjoin, ftoa
 use m_nctk,             only : nctk_try_fort_or_ncfile
 use m_rdmft,            only : rdmft_hello
 use m_pawfgr,           only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_rdmft_functionals, only : RDMFT_SD, RDMFT_MULLER_AS, RDMFT_PNOF5, RDMFT_GNOF, rdmft_active_space, &
                                 rdmft_energy_jk, rdmft_energy_pnof, rdmft_functional_name, rdmft_pnof_table_t
 use m_rdmft_hamiltonian, only : rdmft_one_electron
 use m_rdmft_integrals,  only : rdmft_ints_t
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
subroutine rdmft_driver(codvsn,dtfil,dtset,mpi_enreg,pawrad,pawtab,psps)

 !Arguments ------------------------------------
 !scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(inout) :: psps
 !arrays
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: master=0
 integer :: ierr,ik,iq,ip,jb,jk,ib,mgfftf,nfftf
 real(dp),allocatable :: h1diag(:,:),ekin(:,:)
 real(dp) :: e_corepsp,e_ewald,eone,etot,ekin_tot
 real(dp) :: ecut_eff,ecutdg_eff,ehartree,eexch,gsqcutc_eff,gsqcutf_eff,maxdev,np,nq
 character(len=500) :: msg
 character(len=fnlen) :: wfk_path
 type(pawfgr_type) :: pawfgr
 type(rdmft_ints_t) :: ints
 type(rdmft_pnof_table_t) :: table
 integer :: ierr2,ifunc,nbk,nocc_k
 integer,allocatable :: bandc(:),noccb(:)
 real(dp),allocatable :: eigc(:),h1c(:),nc(:)
 logical,allocatable :: actc(:),occd(:)
 type(rdmft_orbs_t) :: orbs
 !arrays
 integer :: ngfftc(18),ngfftf(18)
 real(dp),parameter :: k0(3)=zero

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

 !FFT mesh and G^2 cut-off containing the products of two orbitals (same as for the density)
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf, &
                  gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=orbs%gmet,k0=k0)

 !J and K integrals between all the orbitals of the full BZ
 call ints%compute(orbs,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut)

 !Hartree and exchange energies per unit cell with the occupations of the ground-state calculation
 ehartree = zero
 eexch = zero
 do ik=1,orbs%nkpt
   do ib=1,orbs%nband(ik,1)
     ip = ints%offset(ik)+ib
     np = orbs%occ(ib,ik,1)/two   ! occupation per spin
     do jk=1,orbs%nkpt
       do jb=1,orbs%nband(jk,1)
         iq = ints%offset(jk)+jb
         nq = orbs%occ(jb,jk,1)/two
         ehartree = ehartree+two*np*nq*ints%jint(ip,iq)
         eexch = eexch-np*nq*ints%kint(ip,iq)
       end do
     end do
   end do
 end do
 ehartree = ehartree/orbs%nkpt**2
 eexch = eexch/orbs%nkpt**2
 write(msg,'(2a,es22.14,a,es22.14,2a)') ch10,' RDMFT: Hartree energy (ground-state occupations) = ',ehartree,&
   ' Ha, exchange energy = ',eexch,' Ha',ch10
 call wrtout([ab_out,std_out],msg)


 !One-electron part: h_pp with getghc, Ewald and core energies
 call rdmft_one_electron(orbs,dtset,dtfil,gsqcutf_eff,mpi_enreg,pawfgr,pawrad,pawtab,psps,h1diag,ekin,e_ewald,e_corepsp)
 eone = zero
 ekin_tot = zero
 do ik=1,orbs%nkpt
   do ib=1,orbs%nband(ik,1)
     eone = eone+orbs%wtk(ik)*orbs%occ(ib,ik,1)*h1diag(ib,ik)
     ekin_tot = ekin_tot+orbs%wtk(ik)*orbs%occ(ib,ik,1)*ekin(ib,ik)
   end do
 end do
 etot = eone+ehartree+eexch+e_ewald+e_corepsp
 write(msg,'(2a,es22.14,2a,es22.14,2a,es22.14,2a,es22.14,2a,es22.14,2a,es22.14,a)') &
   ch10,' RDMFT: kinetic energy                         = ',ekin_tot,&
   ch10,' RDMFT: one-electron energy (T + Vloc + Vnl)    = ',eone,&
   ch10,' RDMFT: Ewald energy                           = ',e_ewald,&
   ch10,' RDMFT: core pseudopotential energy             = ',e_corepsp,&
   ch10,' RDMFT: Hartree + exchange                      = ',ehartree+eexch,&
   ch10,' RDMFT: total energy (SD, ground-state occ.)    = ',etot,ch10
 call wrtout([ab_out,std_out],msg)

 !Energies of all the functionals with the ground-state occupations. Except MULLER_AS, the functionals reduce
 !to the single determinant (Hartree-Fock) energy for integer occupations (up to the fit of ML and MLSIC).
 nbk = ints%nbk
 ABI_MALLOC(nc,(nbk))
 ABI_MALLOC(h1c,(nbk))
 ABI_MALLOC(eigc,(nbk))
 ABI_MALLOC(occd,(nbk))
 ABI_MALLOC(actc,(nbk))
 ABI_MALLOC(bandc,(nbk))
 ABI_MALLOC(noccb,(nbk))
 do ik=1,orbs%nkpt
   nocc_k = count(orbs%occ(1:orbs%nband(ik,1),ik,1) > one)
   do ib=1,orbs%nband(ik,1)
     ip = ints%offset(ik)+ib
     bandc(ip) = ib
     noccb(ip) = nocc_k
     nc(ip) = orbs%occ(ib,ik,1)/two
     h1c(ip) = h1diag(ib,ik)
     eigc(ip) = orbs%eig(ib,ik,1)
     occd(ip) = (nc(ip) > half)
   end do
 end do
 eexch = zero
 do ip=1,nbk
   if (occd(ip)) eexch = eexch+(ints%kint(ip,ip)-ints%jint(ip,ip))
 end do
 write(msg,'(2a,es14.6,a)') ch10,' RDMFT: sum over occupied orbitals of (K_pp - J_pp) / N_k^2 = ',eexch/orbs%nkpt**2,' Ha'
 call wrtout([ab_out,std_out],msg)
 write(msg,'(2a,i0,a)') ch10,' RDMFT: functionals evaluated with the ground-state occupations (rdmft_selfint = ',&
   dtset%rdmft_selfint,'; total energy, difference with SD):'
 call wrtout([ab_out,std_out],msg)
 do ifunc=RDMFT_SD,RDMFT_MULLER_AS
   eone = rdmft_energy_jk(ifunc,orbs%nkpt,nbk,nc,h1c,ints%jint,ints%kint,.not.occd,0.6_dp)+e_ewald+e_corepsp
   write(msg,'(2a,es22.14,a,es10.2)') '   ',rdmft_functional_name(ifunc),eone,'  diff =',eone-etot
   call wrtout([ab_out,std_out],msg)
 end do
 ! Active space (rdmft_nocc_act, rdmft_nvir_act, rdmft_ecut_act) and PNOF subspaces (rdmft_ncoupl)
 call rdmft_active_space(nbk,eigc,bandc,noccb,occd,dtset%rdmft_nocc_act,dtset%rdmft_nvir_act,dtset%rdmft_ecut_act,actc)
 write(msg,'(2a,4(a,i0),a)') ch10,' RDMFT: active space over all k-points: ',' frozen occupied = ',count(occd.and..not.actc),&
   ', active occupied = ',count(occd.and.actc),', active virtual = ',count(.not.occd.and.actc),&
   ', frozen virtual = ',count(.not.occd.and..not.actc),ch10
 call wrtout([ab_out,std_out],msg)
 call table%build(nbk,eigc,occd,actc,dtset%rdmft_ncoupl,ierr2)
 if (ierr2 /= 0) then
   write(msg,'(3a)') ' RDMFT: PNOF functionals skipped: no active occupied orbital or not enough virtual orbitals for',&
     ' rdmft_ncoupl (increase nband of the ground-state calculation or decrease rdmft_nocc_act/rdmft_ncoupl).',ch10
   call wrtout([ab_out,std_out],msg)
 else
   do ifunc=RDMFT_PNOF5,RDMFT_GNOF
     eone = rdmft_energy_pnof(ifunc,table,orbs%nkpt,nbk,nc,h1c,ints%jint,ints%kint,dtset%rdmft_selfint)+e_ewald+e_corepsp
     write(msg,'(2a,a,i0,a,es22.14,a,es10.2)') '   ',rdmft_functional_name(ifunc),' (rdmft_ncoupl=',dtset%rdmft_ncoupl,&
       ')',eone,'  diff =',eone-etot
     call wrtout([ab_out,std_out],msg)
   end do
 end if
 call table%free()
 ABI_FREE(nc)
 ABI_FREE(h1c)
 ABI_FREE(eigc)
 ABI_FREE(occd)
 ABI_FREE(actc)
 ABI_FREE(bandc)
 ABI_FREE(noccb)

 ABI_FREE(h1diag)
 ABI_FREE(ekin)

 call ints%free()
 call pawfgr_destroy(pawfgr)
 call orbs%free()

end subroutine rdmft_driver
!!***

end module m_rdmft_driver
!!***
