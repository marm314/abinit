!!****m* ABINIT/m_rdmft_driver
!! NAME
!!  m_rdmft_driver
!!
!! FUNCTION
!!  Driver for Reduced Density Matrix Functional Theory (RDMFT) calculations (optdriver=11).
!!  The occupation numbers of the active orbitals are optimized for the fixed Bloch orbitals of a ground-state
!!  (DFT or Hartree-Fock) calculation.
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
 use defs_abitypes,       only : MPI_type
 use defs_datatypes,      only : pseudopotential_type
 use m_errors
 use m_xmpi

 use m_dtfil,             only : datafiles_type
 use m_dtset,             only : dataset_type
 use m_fstrings,          only : sjoin
 use m_io_tools,          only : open_file
 use m_nctk,              only : nctk_try_fort_or_ncfile
 use m_pawfgr,            only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_pawrad,            only : pawrad_type
 use m_pawtab,            only : pawtab_type
 use m_rdmft_functionals, only : RDMFT_SD, RDMFT_MULLER_AS, RDMFT_PNOF5, RDMFT_GNOF, rdmft_active_space, rdmft_bg_t, &
                                 rdmft_energy_jk, rdmft_energy_pnof, rdmft_fermi_dirac, rdmft_functional_name, &
                                 rdmft_pnof_table_remap, rdmft_pnof_table_t
 use m_rdmft_hamiltonian, only : rdmft_one_electron
 use m_rdmft_integrals,   only : rdmft_ints_t
 use m_rdmft_optimizer,   only : rdmft_optimize_jk, rdmft_optimize_pnof
 use m_rdmft_orbitals,    only : rdmft_orbs_t, rdmft_rot_matrix
 use m_specialmsg,        only : wrtout

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
!!  1) Read the ground-state (DFT or HF) orbitals from the WFK file (getwfk/irdwfk). For the PNOF functionals the
!!     k-point mesh is completed with the time-reversed orbitals if it is not closed under k -> -k.
!!  2) Active space: frozen occupied orbitals (background, occupation 1), active orbitals (free) and frozen virtual
!!     orbitals (dropped). The J and K integrals are computed in the active space only.
!!  3) Energy of all the functionals with the ground-state occupations.
!!  4) Optimization of the occupation numbers (rdmft_functional): SQP for the JK-only functionals, L-BFGS over the
!!     angles of the subspaces for the PNOF functionals.
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type(datafiles_type)> = infos about file names
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <type(MPI_type)> = MPI information
!!  pawrad(ntypat*usepaw), pawtab(ntypat*usepaw) = PAW data (empty for norm-conserving pseudopotentials)
!!  psps <type(pseudopotential_type)> = pseudopotentials
!!
!! OUTPUT
!!  (only printing and the file *_RDMFT_OCC with the optimized occupations)
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
 integer :: ib,ierr,ierr_t,ifunc,ik,ik_w,ib_w,ip,iq,jp,mgfftf,nadded,nbk,nf,nfftf,niter,nocc_k
 logical :: conv,is_jk,is_pnof,need_table,ok_fd,table_ok
 real(dp),parameter :: eps_box=1.0e-6_dp
 real(dp) :: alpha,e_corepsp,e_ewald,eexch,ehartree,ekin_tot,eone,eopt,ecut_eff,ecutdg_eff,etot,gsqcutc_eff,gsqcutf_eff
 real(dp) :: devtr,maxdev,ntarget,nk2,sexch,shart
 character(len=500) :: msg
 character(len=fnlen) :: wfk_path
 type(pawfgr_type) :: pawfgr
 type(rdmft_bg_t) :: bgt
 type(rdmft_ints_t) :: ints
 type(rdmft_orbs_t) :: orbs
 type(rdmft_pnof_table_t) :: table
 !arrays
 integer :: ngfftc(18),ngfftf(18)
 integer,allocatable :: bandc(:),noccb(:),offset(:),partner(:),partf(:)
 real(dp),parameter :: k0(3)=zero
 real(dp),allocatable :: eigc(:),eigf(:),ekin(:,:),h1c(:),h1diag(:,:),h1f(:),n0(:),nc(:),nfull(:),nf_gs(:),nopt(:)
 complex(dpc),allocatable :: hmat(:,:,:)
 logical,allocatable :: actc(:),bg(:),free(:),occd(:),occf(:)

! *************************************************************************

 ABI_UNUSED(dtset%natom)
 alpha = dtset%rdmft_power
 is_jk = (dtset%rdmft_functional >= RDMFT_SD .and. dtset%rdmft_functional <= RDMFT_MULLER_AS)
 is_pnof = (dtset%rdmft_functional >= RDMFT_PNOF5 .and. dtset%rdmft_functional <= RDMFT_GNOF)

 write(msg,'(3a)') ch10,' RDMFT driver (code version '//trim(codvsn)//')',ch10
 call wrtout([ab_out,std_out],msg)

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

 !Alignment of the orbitals at -k with the time-reversed orbitals of k (assumed by the ties of the occupations of k and -k)
 call orbs%tr_alignment(devtr,ik_w,ib_w)
 if (ik_w > 0) then
   write(msg,'(2a,es10.2,2(a,i0),a)') ' RDMFT: time-reversal alignment of the orbitals at -k with those at k:',&
     ' largest 1-|<psi_-k|Theta psi_k>| = ',devtr,' (k-point ',ik_w,', band ',ib_w,')'
   call wrtout([ab_out,std_out],msg)
   if (devtr > 1.0e-4_dp) then
     write(msg,'(7a)') ' RDMFT: WARNING: the orbitals at -k are not the conjugates of the ones at k (up to a phase), for',&
       ' example in a degenerate shell.',ch10,&
       '        The energy functionals only use K_pq and do not depend on this, but the occupations of k and -k are tied:',&
       ch10,'        the pairing of degenerate states with different occupations is then ambiguous.',ch10
     call wrtout([ab_out,std_out],msg)
   end if
 end if

 !PNOF functionals need subspaces closed under k -> -k: complete the mesh with the time-reversed orbitals if needed
 if (is_pnof) then
   call orbs%complete_time_reversal(nadded)
   if (nadded > 0) then
     write(msg,'(2a,i0,a,i0,a,i0,2a)') ch10,' RDMFT: the k-point mesh is not closed under k -> -k; ',nadded,&
       ' k-points were added with the time-reversed orbitals (total ',orbs%nkpt,' k-points, weight 1/',orbs%nkpt,')',ch10
     call wrtout([ab_out,std_out],msg)
   end if
 end if
 nk2 = real(orbs%nkpt,dp)**2

 !FFT mesh and G^2 cut-off containing the products of two orbitals (same as for the density)
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf, &
                  gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=orbs%gmet,k0=k0)

 !One-electron part: h_pp with getghc, Ewald and core energies
 if (dtset%rdmft_check_rot > 0) then
   if (.not.(is_jk .or. is_pnof)) then
     ABI_ERROR('RDMFT: rdmft_check_rot needs rdmft_functional /= 0 (the occupations at which the gradient is checked).')
   end if
   if (dtset%rdmft_rot_k > orbs%nkpt) then
     ABI_ERROR('RDMFT: rdmft_rot_k is larger than the number of k-points.')
   end if
   call rdmft_one_electron(orbs,dtset,dtfil,gsqcutf_eff,mpi_enreg,pawfgr,pawrad,pawtab,psps,h1diag,ekin,e_ewald,e_corepsp,hmat)
 else
   call rdmft_one_electron(orbs,dtset,dtfil,gsqcutf_eff,mpi_enreg,pawfgr,pawrad,pawtab,psps,h1diag,ekin,e_ewald,e_corepsp)
 end if

 ! ------------------------------------------------------------------
 ! Orbitals with a composite index (band, k-point), active space and free/background sets
 ! ------------------------------------------------------------------
 nbk = sum(orbs%nband(:,1))
 ABI_MALLOC(offset,(orbs%nkpt))
 ABI_MALLOC(nc,(nbk))
 ABI_MALLOC(h1c,(nbk))
 ABI_MALLOC(eigc,(nbk))
 ABI_MALLOC(occd,(nbk))
 ABI_MALLOC(actc,(nbk))
 ABI_MALLOC(bandc,(nbk))
 ABI_MALLOC(noccb,(nbk))
 ABI_MALLOC(partner,(nbk))
 offset(1) = 0
 do ik=2,orbs%nkpt
   offset(ik) = offset(ik-1)+orbs%nband(ik-1,1)
 end do
 do ik=1,orbs%nkpt
   nocc_k = count(orbs%occ(1:orbs%nband(ik,1),ik,1) > one)
   do ib=1,orbs%nband(ik,1)
     ip = offset(ik)+ib
     bandc(ip) = ib
     noccb(ip) = nocc_k
     nc(ip) = orbs%occ(ib,ik,1)/two
     h1c(ip) = h1diag(ib,ik)
     eigc(ip) = orbs%eig(ib,ik,1)
     occd(ip) = (nc(ip) > half)
     ! Time-reversal partner (n,-k): the orbital itself for the time-reversal-invariant k-points and if -k is not in the mesh
     partner(ip) = ip
     if (orbs%tr_k(ik) > 0) partner(ip) = offset(orbs%tr_k(ik))+ib
   end do
 end do

 ! Active space (rdmft_nocc_act, rdmft_nvir_act, rdmft_ecut_act)
 call rdmft_active_space(nbk,eigc,bandc,noccb,occd,dtset%rdmft_nocc_act,dtset%rdmft_nvir_act,dtset%rdmft_ecut_act,actc)

 ! PNOF subspaces (rdmft_ncoupl), needed to optimize a PNOF and to evaluate them with the ground-state occupations
 need_table = (is_pnof .or. dtset%rdmft_functional == 0)
 ierr_t = 1
 if (need_table) call table%build(nbk,eigc,occd,actc,partner,dtset%rdmft_ncoupl,ierr_t)
 table_ok = (need_table .and. ierr_t == 0)
 if (is_pnof .and. .not.table_ok) then
   ABI_ERROR('RDMFT: no active occupied orbital or not enough virtual orbitals for rdmft_ncoupl (PNOF optimization).')
 end if

 ! Free orbitals: the ones whose occupation can change. Background orbitals: the other occupied orbitals (occupation 1).
 ! The frozen virtual orbitals have occupation 0 and are dropped.
 ABI_MALLOC(free,(nbk))
 ABI_MALLOC(bg,(nbk))
 if (is_pnof) then
   free = (table%orb_var > 0)
 else
   free = actc
   if (table_ok) free = free .or. (table%orb_var > 0)
 end if
 bg = occd .and. .not.free
 write(msg,'(2a,4(a,i0),a)') ch10,' RDMFT: active space over all k-points: ',' frozen occupied = ',count(occd.and..not.actc),&
   ', active occupied = ',count(occd.and.actc),', active virtual = ',count(.not.occd.and.actc),&
   ', frozen virtual = ',count(.not.occd.and..not.actc),ch10
 call wrtout([ab_out,std_out],msg)

 ! ------------------------------------------------------------------
 ! J and K integrals in the active space
 ! ------------------------------------------------------------------
 call ints%compute(orbs,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,free,bg,&
                   dtset%rdmft_core_ints == 1)
 nf = ints%nfree
 ABI_MALLOC(nf_gs,(nf))
 ABI_MALLOC(h1f,(nf))
 ABI_MALLOC(eigf,(nf))
 ABI_MALLOC(occf,(nf))
 ABI_MALLOC(partf,(nf))
 do ip=1,nf
   jp = ints%comp_of(ip)
   nf_gs(ip) = nc(jp)
   h1f(ip) = h1c(jp)
   eigf(ip) = eigc(jp)
   occf(ip) = occd(jp)
   partf(ip) = ints%free_of(partner(jp))
   if (partf(ip) == 0) then
     if (is_jk) then
       ABI_ERROR('RDMFT: the active window is not closed under k -> -k.')
     end if
     partf(ip) = ip
   end if
 end do
 bgt%sh = sum(h1c,mask=bg)
 bgt%sjo = ints%sjo; bgt%sko = ints%sko; bgt%sjd = ints%sjd; bgt%skd = ints%skd
 if (table_ok) call rdmft_pnof_table_remap(table,nbk,ints%free_of,nf)

 ! ------------------------------------------------------------------
 ! Energy with the occupations of the ground-state calculation
 ! ------------------------------------------------------------------
 eone = zero
 ekin_tot = zero
 do ik=1,orbs%nkpt
   do ib=1,orbs%nband(ik,1)
     eone = eone+orbs%wtk(ik)*orbs%occ(ib,ik,1)*h1diag(ib,ik)
     ekin_tot = ekin_tot+orbs%wtk(ik)*orbs%occ(ib,ik,1)*ekin(ib,ik)
   end do
 end do
 ! Hartree and exchange energies: free-free, free-background and background-background pairs
 shart = zero
 sexch = zero
 do ip=1,nf
   do iq=1,nf
     shart = shart+two*nf_gs(ip)*nf_gs(iq)*ints%jint(ip,iq)
     sexch = sexch-nf_gs(ip)*nf_gs(iq)*ints%kint(ip,iq)
   end do
   shart = shart+four*nf_gs(ip)*ints%bj(ip)
   sexch = sexch-two*nf_gs(ip)*ints%bk(ip)
 end do
 shart = shart+two*(ints%sjo+ints%sjd)
 sexch = sexch-(ints%sko+ints%skd)
 ehartree = shart/nk2
 eexch = sexch/nk2
 etot = eone+ehartree+eexch+e_ewald+e_corepsp
 write(msg,'(2a,es22.14,2a,es22.14,2a,es22.14,2a,es22.14,2a,es22.14,2a,es22.14,a)') &
   ch10,' RDMFT: kinetic energy                         = ',ekin_tot,&
   ch10,' RDMFT: one-electron energy (T + Vloc + Vnl)    = ',eone,&
   ch10,' RDMFT: Ewald energy                           = ',e_ewald,&
   ch10,' RDMFT: core pseudopotential energy             = ',e_corepsp,&
   ch10,' RDMFT: Hartree + exchange                      = ',ehartree+eexch,&
   ch10,' RDMFT: total energy (SD, ground-state occ.)    = ',etot,ch10
 call wrtout([ab_out,std_out],msg)
 if (dtset%rdmft_core_ints == 0) then
   write(msg,'(3a)') ' RDMFT: rdmft_core_ints=0, the exchange among the frozen occupied orbitals is not included in the ',&
     'energies above and below: only their differences are meaningful.',ch10
   call wrtout([ab_out,std_out],msg)
 else
   write(msg,'(a,es14.6,a)') ' RDMFT: sum over occupied orbitals of (K_pp - J_pp) / N_k^2 = ',diag_sum()/nk2,' Ha'
   call wrtout([ab_out,std_out],msg)
 end if

 write(msg,'(2a,i0,a)') ch10,' RDMFT: functionals evaluated with the ground-state occupations (rdmft_selfint = ',&
   dtset%rdmft_selfint,'; total energy, difference with SD):'
 call wrtout([ab_out,std_out],msg)
 do ifunc=RDMFT_SD,RDMFT_MULLER_AS
   eopt = rdmft_energy_jk(ifunc,orbs%nkpt,nf,nf_gs,h1f,ints%jint,ints%kint,alpha,ints%bj,ints%bk,bgt)+e_ewald+e_corepsp
   write(msg,'(2a,es22.14,a,es10.2)') '   ',rdmft_functional_name(ifunc),eopt,'  diff =',eopt-etot
   call wrtout([ab_out,std_out],msg)
 end do
 if (table_ok) then
   do ifunc=RDMFT_PNOF5,RDMFT_GNOF
     eopt = rdmft_energy_pnof(ifunc,table,orbs%nkpt,nf,nf_gs,h1f,ints%jint,ints%kint,dtset%rdmft_selfint,&
                              ints%bj,ints%bk,bgt)+e_ewald+e_corepsp
     write(msg,'(2a,a,i0,a,es22.14,a,es10.2)') '   ',rdmft_functional_name(ifunc),' (rdmft_ncoupl=',dtset%rdmft_ncoupl,&
       ')',eopt,'  diff =',eopt-etot
     call wrtout([ab_out,std_out],msg)
   end do
 else if (dtset%rdmft_functional == 0) then
   write(msg,'(3a)') ' RDMFT: PNOF functionals skipped: no active occupied orbital or not enough virtual orbitals for',&
     ' rdmft_ncoupl (increase nband of the ground-state calculation or decrease rdmft_nocc_act/rdmft_ncoupl).',ch10
   call wrtout([ab_out,std_out],msg)
 end if

 ! ------------------------------------------------------------------
 ! Optimization of the occupation numbers
 ! ------------------------------------------------------------------
 if (is_jk .or. is_pnof) then
   ntarget = orbs%nelect/two*orbs%nkpt-real(ints%nbg,dp)
   ABI_MALLOC(n0,(nf))
   ABI_MALLOC(nopt,(nf))
   n0 = nf_gs
   if (dtset%rdmft_start == 1) then
     if (is_pnof) then
       ABI_ERROR('RDMFT: rdmft_start=1 (Fermi-Dirac) is not available for the PNOF functionals.')
     end if
     call rdmft_fermi_dirac(nf,eigf,ntarget,dtset%tsmear,n0,ok_fd)
     if (.not.ok_fd) then
       write(msg,'(3a)') ' RDMFT: Fermi-Dirac start not possible (tsmear <= 0 or no free electrons), aufbau start used.',ch10
       call wrtout([ab_out,std_out],msg)
       n0 = nf_gs
     end if
   end if
   if (dtset%rdmft_restart == 1) call read_occ_file(nf,ints%comp_of,offset,orbs%nband(:,1),n0)

   if (is_jk) then
     call rdmft_optimize_jk(dtset%rdmft_functional,orbs%nkpt,nf,h1f,ints%jint,ints%kint,alpha,partf,n0,ntarget,eps_box,&
                            ints%bj,ints%bk,bgt,nopt,eopt,conv,niter)
   else
     call rdmft_optimize_pnof(dtset%rdmft_functional,table,orbs%nkpt,nf,h1f,ints%jint,ints%kint,dtset%rdmft_selfint,&
                              ints%bj,ints%bk,bgt,dtset%rdmft_restart == 1,n0,1.0e-8_dp,500,dtset%prtvol >= 10,&
                              nopt,eopt,conv,niter)
   end if
   write(msg,'(2a,a,i0,a,l1,3a,es22.14,a,es22.14,a,es10.2,a)') ch10,' RDMFT: occupation numbers optimized for the ',&
     'functional ',dtset%rdmft_functional,' (converged = ',conv,')',ch10,&
     ' RDMFT: optimized total energy         = ',eopt+e_ewald+e_corepsp,&
     '; energy with the ground-state occupations = ',etot,'  (optimized - ground state =',eopt+e_ewald+e_corepsp-etot,')'
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,i0,a,es14.6)') ' RDMFT: number of iterations = ',niter,'; sum of the occupations of the free orbitals = ',&
     sum(nopt)
   call wrtout([ab_out,std_out],msg)

   ! Occupations of all the orbitals: background 1, free from the optimization, frozen virtual 0
   ABI_MALLOC(nfull,(nbk))
   nfull = zero
   where (bg) nfull = one
   do ip=1,nf
     nfull(ints%comp_of(ip)) = nopt(ip)
   end do
   write(msg,'(a)') ' RDMFT: occupation numbers of the first k-point (band, energy, n):'
   call wrtout([ab_out,std_out],msg)
   do ib=1,orbs%nband(1,1)
     write(msg,'(i6,f14.6,f16.10)') ib,eigc(offset(1)+ib),nfull(offset(1)+ib)
     call wrtout([ab_out,std_out],msg)
   end do
   if (xmpi_comm_rank(xmpi_world) == master) call write_occ_file(trim(dtfil%filnam_ds(4))//'_RDMFT_OCC')
   if (dtset%rdmft_check_rot > 0) call check_rotations(eopt+e_ewald+e_corepsp)
   ABI_FREE(nfull)
   ABI_FREE(n0)
   ABI_FREE(nopt)
 end if

 ABI_FREE(offset)
 ABI_FREE(nc)
 ABI_FREE(h1c)
 ABI_FREE(eigc)
 ABI_FREE(occd)
 ABI_FREE(actc)
 ABI_FREE(bandc)
 ABI_FREE(noccb)
 ABI_FREE(partner)
 ABI_FREE(free)
 ABI_FREE(bg)
 ABI_FREE(nf_gs)
 ABI_FREE(h1f)
 ABI_FREE(eigf)
 ABI_FREE(occf)
 ABI_FREE(partf)
 ABI_FREE(h1diag)
 ABI_FREE(ekin)
 ABI_SFREE(hmat)
 call table%free()
 call ints%free()
 call pawfgr_destroy(pawfgr)
 call orbs%free()

contains

 function diag_sum() result(dsum)
   ! Sum over the occupied orbitals of K_pp - J_pp (free and background)
   real(dp) :: dsum
   integer :: jf
   dsum = ints%skd-ints%sjd
   do jf=1,nf
     if (occf(jf)) dsum = dsum+ints%kint(jf,jf)-ints%jint(jf,jf)
   end do
 end function diag_sum


 function energy_rot(orbs_x,h1c_x) result(e)
   ! Total energy of the functional at the occupations nfull for the orbitals orbs_x (core Hamiltonian diagonal h1c_x)
   type(rdmft_orbs_t),intent(in) :: orbs_x
   real(dp),intent(in) :: h1c_x(nbk)
   real(dp) :: e
   integer :: jf
   type(rdmft_ints_t) :: ix
   type(rdmft_bg_t) :: bx
   real(dp),allocatable :: h1x(:),nx(:)
   call ix%compute(orbs_x,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,free,bg,&
                   dtset%rdmft_core_ints == 1,quiet=.true.)
   ABI_MALLOC(h1x,(nf))
   ABI_MALLOC(nx,(nf))
   do jf=1,nf
     h1x(jf) = h1c_x(ix%comp_of(jf))
     nx(jf) = nfull(ix%comp_of(jf))
   end do
   bx%sh = sum(h1c_x,mask=bg)
   bx%sjo = ix%sjo; bx%sko = ix%sko; bx%sjd = ix%sjd; bx%skd = ix%skd
   if (is_jk) then
     e = rdmft_energy_jk(dtset%rdmft_functional,orbs%nkpt,nf,nx,h1x,ix%jint,ix%kint,alpha,ix%bj,ix%bk,bx)+e_ewald+e_corepsp
   else
     e = rdmft_energy_pnof(dtset%rdmft_functional,table,orbs%nkpt,nf,nx,h1x,ix%jint,ix%kint,dtset%rdmft_selfint,&
                           ix%bj,ix%bk,bx)+e_ewald+e_corepsp
   end if
   ABI_FREE(h1x)
   ABI_FREE(nx)
   call ix%free()
 end function energy_rot

 subroutine update_h(jk,ib,jb,umat,h1x)
   ! Diagonal elements of the core Hamiltonian of the rotated bands: h'_pp = sum_ij conjg(U_ip) h_ij U_jp
   integer,intent(in) :: jk,ib,jb
   complex(dpc),intent(in) :: umat(2,2)
   real(dp),intent(inout) :: h1x(nbk)
   integer :: ic,ii,jj,bnd(2)
   complex(dpc) :: hl(2,2),acc
   bnd = [ib,jb]
   do ii=1,2
     do jj=1,2
       hl(ii,jj) = hmat(bnd(ii),bnd(jj),jk)
     end do
   end do
   do ic=1,2
     acc = czero
     do ii=1,2
       do jj=1,2
         acc = acc+conjg(umat(ii,ic))*hl(ii,jj)*umat(jj,ic)
       end do
     end do
     h1x(offset(jk)+bnd(ic)) = real(acc,dp)
   end do
 end subroutine update_h

 function energy_rotated(ik,ib,jb,theta,imag) result(e)
   ! Energy after rotating the bands ib,jb of the k-point ik by theta (and the time-reversed ones at -k)
   integer,intent(in) :: ik,ib,jb
   real(dp),intent(in) :: theta
   logical,intent(in) :: imag
   real(dp) :: e
   integer :: kp
   complex(dpc) :: umat(2,2)
   real(dp),allocatable :: h1x(:)
   type(rdmft_orbs_t) :: ow
   ow = orbs
   ABI_MALLOC(h1x,(nbk))
   h1x = h1c
   umat = rdmft_rot_matrix(theta,imag,.false.)
   call ow%rotate_pair(ik,ib,jb,umat)
   call update_h(ik,ib,jb,umat,h1x)
   kp = orbs%tr_k(ik)
   if (kp > 0 .and. kp /= ik) then
     umat = rdmft_rot_matrix(theta,imag,.true.)
     call ow%rotate_pair(kp,ib,jb,umat)
     call update_h(kp,ib,jb,umat,h1x)
   end if
   e = energy_rot(ow,h1x)
   ABI_FREE(h1x)
   call ow%free()
 end function energy_rotated

 subroutine check_rotations(e_ref)
   ! Finite-difference gradient dE/dtheta and curvature of the energy with respect to the rotation of two states of the
   ! same k-point, at fixed occupations
   real(dp),intent(in) :: e_ref
   integer :: ik0,ip_b,iq_b,nb0,ndir
   real(dp) :: e0,ep,em,step
   logical :: imag
   character(len=500) :: cmsg
   ik0 = dtset%rdmft_rot_k
   nb0 = orbs%nband(ik0,1)
   step = dtset%rdmft_rot_step
   e0 = energy_rot(orbs,h1c)
   write(cmsg,'(2a,es20.12,a,es10.2,3a,i0,a,es9.2,a)') ch10,' RDMFT: rotation check, E(0) = ',e0,&
     '  (E(0) - optimized energy =',e0-e_ref,')',ch10,' k-point ',ik0,', step ',step,' radians (central differences)'
   call wrtout([ab_out,std_out],cmsg)
   write(cmsg,'(a)') '   band p  band q    n_p       n_q     direction      dE/dtheta        d2E/dtheta2'
   call wrtout([ab_out,std_out],cmsg)
   do ip_b=1,nb0
     do iq_b=ip_b+1,nb0
       if (dtset%rdmft_check_rot == 1) then
         if (ip_b /= min(dtset%rdmft_rot_p,dtset%rdmft_rot_q) .or. iq_b /= max(dtset%rdmft_rot_p,dtset%rdmft_rot_q)) cycle
       end if
       do ndir=1,2
         imag = (ndir == 2)
         ep = energy_rotated(ik0,ip_b,iq_b,step,imag)
         em = energy_rotated(ik0,ip_b,iq_b,-step,imag)
         write(cmsg,'(2i8,2f10.6,5x,a,2es18.8)') ip_b,iq_b,nfull(offset(ik0)+ip_b),nfull(offset(ik0)+iq_b),&
           merge('imaginary','real     ',imag),(ep-em)/(two*step),(ep-two*e0+em)/step**2
         call wrtout([ab_out,std_out],cmsg)
       end do
     end do
   end do
 end subroutine check_rotations

 subroutine write_occ_file(fname)
   ! File of the occupation numbers per spin of all the orbitals: first line nkpt and nband, then "ik ib n"
   character(len=*),intent(in) :: fname
   integer :: unt,jk,jb
   character(len=500) :: emsg
   if (open_file(fname,emsg,newunit=unt,status='replace',form='formatted') /= 0) then
     ABI_ERROR(emsg)
   end if
   write(unt,'(2i8)') orbs%nkpt,maxval(orbs%nband(:,1))
   do jk=1,orbs%nkpt
     do jb=1,orbs%nband(jk,1)
       write(unt,'(2i8,f22.14)') jk,jb,nfull(offset(jk)+jb)
     end do
   end do
   close(unt)
   write(emsg,'(3a)') ch10,' RDMFT: optimized occupations written in ',trim(fname)
   call wrtout([ab_out,std_out],emsg)
 end subroutine write_occ_file

 subroutine read_occ_file(nfree,comp_of,off,nband,nin)
   ! Read RDMFT_OCC.in: occupations of the free orbitals (the other ones keep their values)
   integer,intent(in) :: nfree,comp_of(nfree),off(:),nband(:)
   real(dp),intent(inout) :: nin(nfree)
   integer :: unt,jk,jb,nk_file,nb_file,ios,jf,nread
   real(dp) :: val
   real(dp),allocatable :: occ_file(:)
   logical,allocatable :: found(:)
   character(len=500) :: emsg
   if (open_file('RDMFT_OCC.in',emsg,newunit=unt,status='old',form='formatted') /= 0) then
     ABI_ERROR(sjoin('RDMFT: rdmft_restart=1 but the file RDMFT_OCC.in cannot be read:',ch10,emsg))
   end if
   read(unt,*) nk_file,nb_file
   if (nk_file /= size(off)) then
     ABI_ERROR('RDMFT: the number of k-points of RDMFT_OCC.in is different from the one of the calculation.')
   end if
   ABI_MALLOC(occ_file,(sum(nband)))
   ABI_MALLOC(found,(sum(nband)))
   found = .false.
   nread = 0
   do
     read(unt,*,iostat=ios) jk,jb,val
     if (ios /= 0) exit
     if (jk < 1 .or. jk > size(off) .or. jb < 1 .or. jb > nband(jk)) cycle
     occ_file(off(jk)+jb) = val
     found(off(jk)+jb) = .true.
     nread = nread+1
   end do
   close(unt)
   do jf=1,nfree
     if (found(comp_of(jf))) nin(jf) = min(one,max(zero,occ_file(comp_of(jf))))
   end do
   write(emsg,'(2a,i0,a)') ch10,' RDMFT: occupations read from RDMFT_OCC.in for ',nread,' orbitals'
   call wrtout([ab_out,std_out],emsg)
   ABI_FREE(occ_file)
   ABI_FREE(found)
 end subroutine read_occ_file

end subroutine rdmft_driver
!!***

end module m_rdmft_driver
!!***
