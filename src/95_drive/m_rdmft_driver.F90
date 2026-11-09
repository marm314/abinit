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
 use m_crystal,          only : crystal_t
 use m_hdr,               only : hdr_type
 use m_ioarr,             only : fftdatar_write
 use m_io_tools,          only : get_unit, open_file
 use m_nctk,              only : nctk_ncify, nctk_try_fort_or_ncfile
 use m_pawfgr,            only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_pawrad,            only : pawrad_type
 use m_pawtab,            only : pawtab_type
 use m_rdmft_adam,        only : rdmft_adam_t, rdmft_expm
 use m_rdmft_fock,        only : rdmft_gen_fock_k
 use m_rdmft_functionals, only : RDMFT_SD, RDMFT_MULLER_AS, RDMFT_PNOF5, RDMFT_GNOF, rdmft_active_space, rdmft_bg_t, &
                                 rdmft_coef_jk, rdmft_coef_pnof, rdmft_energy_jk, rdmft_energy_pnof, rdmft_fermi_dirac, &
                                 rdmft_functional_name, rdmft_pnof_table_remap, rdmft_pnof_table_t
 use m_rdmft_hamiltonian, only : rdmft_one_electron
 use m_rdmft_fock,        only : rdmft_fock_hartree_k
 use m_mpinfo,            only : destroy_mpi_enreg, initmpi_seq
 use m_rdmft_integrals,   only : orb_to_real_one, rdmft_fockreq_t, rdmft_ints_t, rdmft_mem_on, rdmft_tint, rdmft_tprof
 use m_rdmft_optimizer,   only : rdmft_grad_jk, rdmft_optimize_jk, rdmft_optimize_pnof, rdmft_tpnof
 use m_rdmft_orbitals,    only : rdmft_orbs_t, rdmft_rot_matrix
 use m_specialmsg,        only : wrtout
 use m_time,              only : abi_wtime

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
 integer :: ib,ierr,ierr_t,ifunc,ik,ik_w,ib_w,ip,iq,jk,jp,mgfftf,nadded,nbk,nf,nfftf,niter,nocc_k
 logical :: conv,is_jk,is_pnof,need_table,ok_fd,table_ok
 real(dp) :: alpha,e_corepsp,e_ewald,eexch,ehartree,ekin_tot,eone,eopt,ecut_eff,ecutdg_eff,etot,gsqcutc_eff,gsqcutf_eff
 real(dp) :: devtr,maxdev,ntarget,nk2,sexch,shart
 character(len=500) :: msg
 character(len=100) :: eline(6)
 character(len=fnlen) :: wfk_path
 type(pawfgr_type) :: pawfgr
 type(MPI_type) :: mpi_seq
 type(rdmft_bg_t) :: bgt
 type(rdmft_ints_t) :: ints
 type(rdmft_orbs_t) :: orbs
 type(rdmft_pnof_table_t) :: table
 !arrays
 integer :: ngfftc(18),ngfftf(18)
 integer,allocatable :: bandc(:),noccb(:),offset(:),partner(:),partf(:)
 real(dp),parameter :: k0(3)=zero
 real(dp),allocatable :: rho_ref(:)
 real(dp),allocatable :: eigc(:),eigf(:),ekin(:,:),h1c(:),h1diag(:,:),h1f(:),n0(:),nc(:),nfull(:),nf_gs(:),nopt(:)
 complex(dpc),allocatable :: hmat(:,:,:)
 logical,allocatable :: actc(:),bg(:),free(:),occd(:),occf(:)
 real(dp) :: tm(8),t0,tint_max(6),tint_min(6),t_start,tst(7)
 integer :: nrel_last,nmem_eval

! *************************************************************************

 ABI_UNUSED(dtset%natom)
 alpha = dtset%rdmft_power
 tm = zero
 rdmft_tint = zero
 rdmft_tprof = zero
 rdmft_tpnof = zero
 t_start = abi_wtime()
 tst = zero
 nmem_eval = 0
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
 tst(1) = abi_wtime()-t_start
 !The plane-wave basis of the orbitals is the one of the ground-state calculation: ecut must be the same
 if (abs(orbs%ecut-dtset%ecut) > tol8) then
   write(msg,'(2(a,es14.6),3a)') 'RDMFT: ecut = ',dtset%ecut,' in the RDMFT dataset but ',orbs%ecut,&
     ' in the ground-state WFK file. The cutoff (as well as the cell, the k-points and the pseudopotentials) must be those',&
     ' of the ground-state calculation: the orbitals are expanded in its plane-wave basis, and the FFT mesh and the',&
     ' one-electron operator depend on it.'
   ABI_ERROR(msg)
 end if
 if (dtset%rdmft_complete == 1) then
   call orbs%complete_basis(nadded)
   write(msg,'(2a,i0,a,i0,a)') ch10,' RDMFT: rdmft_complete = 1, bands added to complete the plane-wave basis of each k-point: ',&
     nadded,' (total number of orbitals ',sum(orbs%nband(:,1)),')'
   call wrtout([ab_out,std_out],msg)
 end if

 !Size of the orbital space with respect to the plane-wave basis of ecut
 write(msg,'(2a,i0,a,i0,a,i0,a,i0,a,f6.1,a)') ch10,' RDMFT: bands per k-point: ',minval(orbs%nband(:,1)),'-',&
   maxval(orbs%nband(:,1)),'; plane waves per k-point (ecut): ',minval(orbs%npw),'-',maxval(orbs%npw),' (the bands are ',&
   100.0_dp*real(sum(orbs%nband(:,1)),dp)/real(sum(orbs%npw),dp),' % of the plane-wave basis)'
 call wrtout([ab_out,std_out],msg)
 if (all(orbs%nband(:,1) == orbs%npw)) then
   write(msg,'(2a)') ' RDMFT: the plane-wave basis is complete: the occupations and the orbital rotations span all the',&
     ' states that ecut can produce.'
 else
   write(msg,'(3a)') ' RDMFT: the occupations and the orbital rotations are restricted to these bands (a variational ',&
     'restriction of the orbital space, the energies are upper bounds of the full-basis ones);',&
     ' use rdmft_complete = 1 to span the whole plane-wave basis.'
 end if
 call wrtout([ab_out,std_out],msg)

 maxdev = orbs%check_orthonormality()
 write(msg,'(2a,es12.4,a)') ch10,' RDMFT: maximum deviation from orthonormality of the orbitals = ',maxdev,ch10
 call wrtout([ab_out,std_out],msg)

 !Alignment of the orbitals at -k with the time-reversed orbitals of k (assumed by the ties of the occupations of k and -k)
 call orbs%tr_alignment(devtr,ik_w,ib_w)
 if (ik_w > 0) then
   write(msg,'(2a,es10.2,2(a,i0),a)') ' RDMFT: time-reversal alignment of the orbitals at -k with those at k:',&
     ' largest 1-Re<psi_-k|Theta psi_k> = ',devtr,' (k-point ',ik_w,', band ',ib_w,')'
   call wrtout([ab_out,std_out],msg)
   if (devtr > 1.0e-4_dp) then
     write(msg,'(7a)') ' RDMFT: WARNING: the orbitals at -k are not the conjugates of the ones at k (up to a phase), for',&
       ' example in a degenerate shell.',ch10,&
       '        The energy functionals only use K_pq and do not depend on this, but the occupations of k and -k are tied:',&
       ch10,'        the pairing of degenerate states with different occupations is then ambiguous.',ch10
     call wrtout([ab_out,std_out],msg)
   end if
 end if

 !Enforce the band-by-band time-reversal correspondence of the orbitals at k and -k (degenerate shells)
 if (ik_w > 0 .and. devtr > 1.0e-8_dp) then
   call orbs%align_time_reversal(nadded)
   call orbs%tr_alignment(devtr,ik_w,ib_w)
   write(msg,'(a,i0,a,es10.2,a)') ' RDMFT: the orbitals of ',nadded,' k-points were replaced by the time-reversed ones of -k;'//&
     ' alignment after the replacement: ',devtr,ch10
   call wrtout([ab_out,std_out],msg)
 end if

 !Complete the mesh with the time-reversed orbitals (c_{n,-k}(G) = conjg(c_{nk}(-G))) if the mesh is not closed under
 !k -> -k, for all the functionals (for example a WFK file computed with the k-points reduced by time reversal)
 call orbs%complete_time_reversal(nadded)
 if (nadded > 0) then
   write(msg,'(2a,i0,a,i0,a,i0,2a)') ch10,' RDMFT: the k-point mesh is not closed under k -> -k; ',nadded,&
     ' k-points were added with the time-reversed orbitals (total ',orbs%nkpt,' k-points, weight 1/',orbs%nkpt,')',ch10
   call wrtout([ab_out,std_out],msg)
 end if
 nk2 = real(orbs%nkpt,dp)**2
 tst(2) = abi_wtime()-t_start

 !FFT mesh and G^2 cut-off containing the products of two orbitals (same as for the density)
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf, &
                  gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=orbs%gmet,k0=k0)

 !Sequential MPI_type for the FFTs of the orbitals (density)
 call initmpi_seq(mpi_seq)
 call mpi_seq%distribfft%init_seq('c',ngfftf(2),ngfftf(3),'all')
 call mpi_seq%distribfft%init_seq('f',ngfftf(2),ngfftf(3),'all')

 !One-electron part: h_pp with getghc, Ewald and core energies
 if (dtset%rdmft_check_rot > 0 .or. dtset%rdmft_orb_macro > 0) then
   if (.not.(is_jk .or. is_pnof)) then
     ABI_ERROR('RDMFT: rdmft_check_rot and rdmft_orb_macro need rdmft_functional /= 0 (the occupations to use).')
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


 !With rdmft_tr_reduce = 1 the matrices of the integrals are filled by time reversal: the active space (free and
 !background orbitals) must be the same at k and -k
 if (dtset%rdmft_tr_reduce == 1) then
   do ik=1,orbs%nkpt
     jk = orbs%tr_k(ik)
     if (jk <= 0 .or. jk == ik) cycle
     if (orbs%nband(jk,1) /= orbs%nband(ik,1)) then
       write(msg,'(2(a,i0),2a)') 'RDMFT: the k-points ',ik,' and ',jk,' (time-reversed partners) have different ',&
         'numbers of bands. Use rdmft_tr_reduce = 0.'
       ABI_ERROR(msg)
     end if
     do ib=1,orbs%nband(ik,1)
       if ((free(offset(ik)+ib) .neqv. free(offset(jk)+ib)) .or. (bg(offset(ik)+ib) .neqv. bg(offset(jk)+ib))) then
         write(msg,'(3(a,i0),3a)') 'RDMFT: the active space of band ',ib,' differs between the k-point ',ik,&
           ' and its time-reversed partner ',jk,'. The active window (rdmft_ecut_act, rdmft_nocc_act, rdmft_nvir_act)',&
           ' must be symmetric under k -> -k for rdmft_tr_reduce = 1; use rdmft_tr_reduce = 0.'
         ABI_ERROR(msg)
       end if
     end do
   end do
 end if

 ! ------------------------------------------------------------------
 ! J and K integrals in the active space
 ! ------------------------------------------------------------------
 tst(3) = abi_wtime()-t_start
 call ints%compute(orbs,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,free,bg,&
                   dtset%rdmft_core_ints == 1,reduce=(dtset%rdmft_tr_reduce == 1))
 tst(4) = abi_wtime()-t_start
 call report_mem('after the integrals')
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
 ! Energies in a column: label, then the value at a fixed position
 write(eline(1),'(a,t52,a,es22.14)') ' RDMFT: kinetic energy',' = ',ekin_tot
 write(eline(2),'(a,t52,a,es22.14)') ' RDMFT: one-electron energy (T + Vloc + Vnl)',' = ',eone
 write(eline(3),'(a,t52,a,es22.14)') ' RDMFT: Ewald energy',' = ',e_ewald
 write(eline(4),'(a,t52,a,es22.14)') ' RDMFT: core pseudopotential energy',' = ',e_corepsp
 write(eline(5),'(a,t52,a,es22.14)') ' RDMFT: Hartree + exchange',' = ',ehartree+eexch
 write(eline(6),'(a,t52,a,es22.14)') ' RDMFT: total energy (SD, ground-state occ.)',' = ',etot
 write(msg,'(7a)') ch10,(trim(eline(ip))//ch10,ip=1,6)
 call wrtout([ab_out,std_out],msg)
 if (dtset%rdmft_core_ints == 0) then
   write(msg,'(3a)') ' RDMFT: rdmft_core_ints=0, the exchange among the frozen occupied orbitals is not included in the ',&
     'energies above and below: only their differences are meaningful.',ch10
   call wrtout([ab_out,std_out],msg)
 else
   write(msg,'(a,es14.6,a)') ' RDMFT: sum over occupied orbitals of (K_pp - J_pp) / N_k^2 = ',diag_sum()/nk2,' Ha'
   call wrtout([ab_out,std_out],msg)
 end if

 if (dtset%rdmft_functional == 0) then
   write(msg,'(2a,i0,a)') ch10,' RDMFT: all the functionals evaluated with the ground-state occupations (rdmft_selfint = ',&
     dtset%rdmft_selfint,'; total energy, difference with SD):'
 else
   write(msg,'(2a,i0,a)') ch10,' RDMFT: functional evaluated with the ground-state occupations (rdmft_selfint = ',&
     dtset%rdmft_selfint,'; total energy, difference with SD):'
 end if
 call wrtout([ab_out,std_out],msg)
 do ifunc=RDMFT_SD,RDMFT_MULLER_AS
   if (dtset%rdmft_functional /= 0 .and. ifunc /= dtset%rdmft_functional) cycle
   eopt = rdmft_energy_jk(ifunc,orbs%nkpt,nf,nf_gs,h1f,ints%jint,ints%kint,alpha,ints%bj,ints%bk,bgt)+e_ewald+e_corepsp
   write(msg,'(2a,es22.14,a,es10.2)') '   ',rdmft_functional_name(ifunc),eopt,'  diff =',eopt-etot
   call wrtout([ab_out,std_out],msg)
 end do
 if (table_ok) then
   do ifunc=RDMFT_PNOF5,RDMFT_GNOF
     if (dtset%rdmft_functional /= 0 .and. ifunc /= dtset%rdmft_functional) cycle
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
   ABI_MALLOC(nfull,(nbk))
   ABI_MALLOC(rho_ref,(nfftf))
   call compute_density(nc,rho_ref)
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

   call optimize_occupations(n0,eopt,conv,niter,.true.)
   tst(5) = abi_wtime()-t_start
   call report_mem('after the first occupation optimization')
   write(msg,'(2a,i0,a,l1,a)') ch10,' RDMFT: occupation numbers optimized for the functional ',dtset%rdmft_functional,&
     ' (converged = ',conv,')'
   call wrtout([ab_out,std_out],msg)
   write(eline(1),'(a,t52,a,es22.14)') ' RDMFT: optimized total energy',' = ',eopt+e_ewald+e_corepsp
   write(eline(2),'(a,t52,a,es22.14)') ' RDMFT: energy with the ground-state occupations',' = ',etot
   write(eline(3),'(a,t52,a,es22.14)') ' RDMFT: optimized - ground state',' = ',eopt+e_ewald+e_corepsp-etot
   write(msg,'(3a)') trim(eline(1))//ch10,trim(eline(2))//ch10,trim(eline(3))
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,i0,a,es14.6)') ' RDMFT: number of iterations = ',niter,'; sum of the occupations of the free orbitals = ',&
     sum(nopt)
   call wrtout([ab_out,std_out],msg)
   call assemble_nfull()

   ! Optimization of the orbitals (ADAM) alternated with the one of the occupations
   if (dtset%rdmft_orb_macro > 0) then
     call orbital_macro(eopt+e_ewald+e_corepsp)
     tst(6) = abi_wtime()-t_start
     eopt = energy_rot(orbs,h1c)-e_ewald-e_corepsp
   end if

   call print_occupations()
   call print_final_energies()
   if (dtset%prtwf > 0 .or. dtset%prtden > 0) call write_natural_orbitals()
   tst(7) = abi_wtime()-t_start
   call report_mem('end of the run')
   write(msg,'(a,7f9.2)') ' RDMFT: elapsed (s) after read WFK / completion+TR / orbital setup / integrals / first occupations'//&
     ' / orbital macro / output =',tst
   call wrtout([ab_out,std_out],msg)
   if (xmpi_comm_rank(xmpi_world) == master) call write_occ_file(trim(dtfil%filnam_ds(4))//'_RDMFT_OCC')
   if (dtset%rdmft_check_rot > 0) call check_rotations(eopt+e_ewald+e_corepsp)
   ABI_FREE(nfull)
   ABI_FREE(rho_ref)
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
 call destroy_mpi_enreg(mpi_seq)
 call orbs%free()

contains

 subroutine report_mem(label)
   ! Peak and current resident memory of this process (Linux /proc/self/status), printed by the master
   character(len=*),intent(in) :: label
   character(len=200) :: line
   character(len=60) :: hwm,rss
   integer :: ios,unt
   hwm = '?'; rss = '?'
   unt = get_unit()
   open(unit=unt,file='/proc/self/status',status='old',action='read',iostat=ios)
   if (ios == 0) then
     do
       read(unt,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (line(1:6) == 'VmHWM:') hwm = adjustl(line(7:))
       if (line(1:6) == 'VmRSS:') rss = adjustl(line(7:))
     end do
     close(unt)
   end if
   if (xmpi_comm_rank(xmpi_world) == master) then
     write(msg,'(4a)') ' RDMFT: memory of the master process ',trim(label),': peak '//trim(hwm)//', current ',trim(rss)
     call wrtout([ab_out,std_out],msg)
   end if
 end subroutine report_mem

 subroutine optimize_occupations(nstart,e_out,conv_out,nit_out,talkative)
   ! Optimization of the occupations (nopt) for the current integrals and one-electron elements
   real(dp),intent(in) :: nstart(:)
   real(dp),intent(out) :: e_out
   logical,intent(out) :: conv_out
   integer,intent(out) :: nit_out
   logical,intent(in) :: talkative
   logical :: done
   if (is_jk .and. dtset%rdmft_occ_skip > zero) then
     call optimize_occ_guarded(nstart,e_out,conv_out,nit_out,talkative,done)
     if (done) return
   end if
   if (is_jk) then
     call rdmft_optimize_jk(dtset%rdmft_functional,orbs%nkpt,nf,h1f,ints%jint,ints%kint,alpha,partf,nstart,ntarget,&
                            dtset%rdmft_occ_min,&
                            ints%bj,ints%bk,bgt,nopt,e_out,conv_out,nit_out,verbose=talkative,e_shift=e_ewald+e_corepsp)
   else
     call rdmft_optimize_pnof(dtset%rdmft_functional,table,orbs%nkpt,nf,h1f,ints%jint,ints%kint,dtset%rdmft_selfint,&
                              ints%bj,ints%bk,bgt,dtset%rdmft_restart == 1,nstart,1.0e-8_dp,500,dtset%prtvol >= 10,&
                              nopt,e_out,conv_out,nit_out,verbose=talkative,e_shift=e_ewald+e_corepsp)
   end if
 end subroutine optimize_occupations

 subroutine optimize_occ_guarded(nstart,e_out,conv_out,nit_out,talkative,done)
   ! Optimization of the occupations of the JK-only functionals with the (almost) empty orbitals held at the lower bound:
   ! an orbital is fixed if its couplings in the energy are below rdmft_occ_skip and its occupation is at the lower bound.
   ! The solution is accepted only if the exact gradient of the full energy satisfies the Karush-Kuhn-Tucker conditions
   ! for all the fixed orbitals (dE/dn_p - mu > 0, mu = multiplier of the electron count, from the orbitals that are in
   ! the interior of the box); the orbitals that violate them are freed and the problem is solved again. If it cannot
   ! be decided, done = .false. and the full problem is solved.
   real(dp),intent(in) :: nstart(:)
   real(dp),intent(out) :: e_out
   logical,intent(out) :: conv_out,done
   integer,intent(out) :: nit_out
   logical,intent(in) :: talkative
   integer :: ib2,jf,jk2,nb2,nfix,nr,nviol,r,round,it2
   real(dp) :: lb,mu,ntr,spread,gmin,er
   logical :: convr
   logical,allocatable :: atub(:),fixed(:),interior(:),newfree(:),relc(:)
   integer,allocatable :: idxr(:),posr(:),partr(:)
   real(dp),allocatable :: acf(:,:),ccf(:,:),nux(:),nsave(:),gorb(:),h1r(:),jr(:,:),kr(:,:),bjr(:),bkr(:),n0r(:),nrout(:)
   done = .false.
   conv_out = .false.
   nit_out = 0
   e_out = zero
   lb = dtset%rdmft_occ_min
   ! Relevance of each orbital at the starting occupations: its largest coupling coefficient
   ABI_MALLOC(nsave,(nbk))
   ABI_MALLOC(relc,(nbk))
   nsave = nfull
   nfull = zero
   where (bg) nfull = one
   do jf=1,nf
     nfull(ints%comp_of(jf)) = min(one-lb,max(lb,nstart(jf)))
   end do
   do jk2=1,orbs%nkpt
     nb2 = orbs%nband(jk2,1)
     ABI_MALLOC(acf,(nb2,nbk))
     ABI_MALLOC(ccf,(nb2,nbk))
     ABI_MALLOC(nux,(nb2))
     call build_coefs(jk2,acf,ccf,nux)
     do ib2=1,nb2
       relc(offset(jk2)+ib2) = (nk2*max(maxval(abs(acf(ib2,:))),maxval(abs(ccf(ib2,:)))) >= dtset%rdmft_occ_skip)
     end do
     ABI_FREE(acf)
     ABI_FREE(ccf)
     ABI_FREE(nux)
   end do
   nfull = nsave
   ABI_FREE(nsave)
   ABI_MALLOC(fixed,(nf))
   do jf=1,nf
     fixed(jf) = (.not.relc(ints%comp_of(jf))) .and. (nstart(jf) <= max(ten*lb,tol14))
   end do
   ABI_FREE(relc)
   ! the time-reversed partners are fixed together
   do jf=1,nf
     if (fixed(jf) .and. .not.fixed(partf(jf))) fixed(jf) = .false.
   end do
   do jf=1,nf
     if (fixed(jf) .and. .not.fixed(partf(jf))) fixed(jf) = .false.
   end do
   if (count(fixed) < max(10,nf/10)) then
     ABI_FREE(fixed)
     return
   end if
   ABI_MALLOC(gorb,(nf))
   ABI_MALLOC(posr,(nf))
   ABI_MALLOC(newfree,(nf))
   gmin = zero
   spread = zero
   mu = zero
   do round=1,10
     nfix = count(fixed)
     nr = nf-nfix
     ABI_MALLOC(idxr,(nr))
     r = 0
     posr = 0
     do jf=1,nf
       if (.not.fixed(jf)) then
         r = r+1
         idxr(r) = jf
         posr(jf) = r
       end if
     end do
     ABI_MALLOC(partr,(nr))
     ABI_MALLOC(h1r,(nr))
     ABI_MALLOC(jr,(nr,nr))
     ABI_MALLOC(kr,(nr,nr))
     ABI_MALLOC(bjr,(nr))
     ABI_MALLOC(bkr,(nr))
     ABI_MALLOC(n0r,(nr))
     ABI_MALLOC(nrout,(nr))
     do r=1,nr
       jf = idxr(r)
       partr(r) = posr(partf(jf))
       h1r(r) = h1f(jf)
       bjr(r) = ints%bj(jf)
       bkr(r) = ints%bk(jf)
       n0r(r) = nstart(jf)
     end do
     do r=1,nr
       jr(:,r) = ints%jint(idxr,idxr(r))
       kr(:,r) = ints%kint(idxr,idxr(r))
     end do
     ntr = ntarget-lb*real(nfix,dp)
     call rdmft_optimize_jk(dtset%rdmft_functional,orbs%nkpt,nr,h1r,jr,kr,alpha,partr,n0r,ntr,lb,bjr,bkr,bgt,nrout,er,&
                            convr,it2,verbose=talkative,e_shift=e_ewald+e_corepsp)
     nit_out = nit_out+it2
     nopt = lb
     do r=1,nr
       nopt(idxr(r)) = nrout(r)
     end do
     ! KKT conditions of the full problem for the fixed orbitals, with the exact gradient
     call rdmft_grad_jk(dtset%rdmft_functional,orbs%nkpt,nf,h1f,ints%jint,ints%kint,alpha,nopt,ints%bj,ints%bk,gorb)
     ABI_MALLOC(interior,(nf))
     ABI_MALLOC(atub,(nf))
     interior = .false.
     atub = .false.
     do r=1,nr
       interior(idxr(r)) = (nrout(r) > ten*lb+tol14 .and. nrout(r) < one-ten*lb-tol14)
       atub(idxr(r)) = (nrout(r) >= one-ten*lb-tol14)
     end do
     ABI_FREE(partr)
     ABI_FREE(h1r)
     ABI_FREE(jr)
     ABI_FREE(kr)
     ABI_FREE(bjr)
     ABI_FREE(bkr)
     ABI_FREE(n0r)
     ABI_FREE(nrout)
     ABI_FREE(idxr)
     ! Lower end of the multiplier of the electron count: the largest dE/dn of the orbitals that are in the interior of the
     ! box (they all must have dE/dn = mu) or at the upper bound (dE/dn <= mu)
     if (count(interior .or. atub) == 0) then
       ABI_FREE(interior)
       ABI_FREE(atub)
       exit
     end if
     if (count(interior) > 0) then
       mu = sum(gorb,mask=interior)/real(count(interior),dp)
       spread = maxval(abs(gorb-mu),mask=interior)
     else
       mu = maxval(gorb,mask=atub)
       spread = zero
     end if
     if (spread > 1.0e-6_dp .or. maxval(gorb,mask=(interior .or. atub)) > mu+1.0e-6_dp) then
       ABI_FREE(interior)
       ABI_FREE(atub)
       exit
     end if
     newfree = fixed .and. (gorb-mu < 1.0e-9_dp)
     nviol = count(newfree)
     gmin = minval(gorb-mu,mask=fixed)
     ABI_FREE(interior)
     ABI_FREE(atub)
     if (nviol == 0) then
       done = .true.
       exit
     end if
     fixed = fixed .and. .not.newfree
     do jf=1,nf
       if (fixed(jf) .and. .not.fixed(partf(jf))) fixed(jf) = .false.
     end do
     if (count(fixed) < max(10,nf/10)) exit
   end do
   if (done) then
     conv_out = convr
     e_out = rdmft_energy_jk(dtset%rdmft_functional,orbs%nkpt,nf,nopt,h1f,ints%jint,ints%kint,alpha,ints%bj,ints%bk,bgt)
     if (talkative) then
       write(msg,'(a,i0,a,i0,2a,es10.2,a,es10.2,a)') ' RDMFT: occupations optimized for ',nf-count(fixed),' of ',nf,&
         ' free orbitals, the others held at the lower bound (rdmft_occ_skip); KKT check of the fixed orbitals:',&
         ' min(dE/dn - mu) = ',gmin,', spread of dE/dn of the interior orbitals = ',spread,ch10
       call wrtout([ab_out,std_out],msg)
     end if
   else
     if (talkative) then
       write(msg,'(2a)') ' RDMFT: the reduced occupation problem could not be verified; the full problem is solved.',ch10
       call wrtout([ab_out,std_out],msg)
     end if
   end if
   ABI_FREE(fixed)
   ABI_FREE(gorb)
   ABI_FREE(posr)
   ABI_FREE(newfree)
 end subroutine optimize_occ_guarded

 subroutine assemble_nfull()
   ! Occupations of all the orbitals: background 1, free from the optimization, frozen virtual 0
   integer :: jf
   nfull = zero
   where (bg) nfull = one
   do jf=1,nf
     nfull(ints%comp_of(jf)) = nopt(jf)
   end do
 end subroutine assemble_nfull

 subroutine apply_u(jk,umat)
   ! Rotate the bands of the k-point jk with U, and the core Hamiltonian matrix: h' = U^dagger h U
   integer,intent(in) :: jk
   complex(dpc),intent(in) :: umat(:,:)
   integer :: jb2,nb2
   nb2 = orbs%nband(jk,1)
   call orbs%rotate_k(jk,umat)
   hmat(1:nb2,1:nb2,jk) = matmul(conjg(transpose(umat)),matmul(hmat(1:nb2,1:nb2,jk),umat))
   do jb2=1,nb2
     h1c(offset(jk)+jb2) = real(hmat(jb2,jb2,jk),dp)
   end do
 end subroutine apply_u

 subroutine rotate_all(kpar,dx)
   ! Rotation of all the k-points with the ADAM step dx: kappa_pq = dx (p<q), kappa_qp = -conjg(kappa_pq), U = exp(kappa);
   ! the time-reversed k-point -k gets the conjugate matrix
   integer,intent(in) :: kpar(:)
   complex(dpc),intent(in) :: dx(:)
   integer :: jk,kk,ipar,pp,qq,nb2
   complex(dpc),allocatable :: kap(:,:),um(:,:)
   do jk=1,orbs%nkpt
     if (kpar(jk) < 0) cycle
     nb2 = orbs%nband(jk,1)
     ABI_MALLOC(kap,(nb2,nb2))
     ABI_MALLOC(um,(nb2,nb2))
     kap = czero
     ipar = kpar(jk)
     do pp=1,nb2
       do qq=pp+1,nb2
         ipar = ipar+1
         kap(pp,qq) = dx(ipar)
         if (orbs%tr_k(jk) == jk) kap(pp,qq) = cmplx(real(dx(ipar),dp),zero,kind=dpc)
         kap(qq,pp) = -conjg(kap(pp,qq))
       end do
     end do
     um = rdmft_expm(kap)
     call apply_u(jk,um)
     kk = orbs%tr_k(jk)
     if (kk > 0 .and. kk /= jk) call apply_u(kk,conjg(um))
     ABI_FREE(kap)
     ABI_FREE(um)
   end do
 end subroutine rotate_all

 subroutine check_h()
   ! Compare the matrix h updated by the rotations with the one recomputed from the rotated orbitals
   real(dp) :: e2,ew2,cp2,dmax,dmaxd,e1
   integer :: jk,nb2
   real(dp),allocatable :: h1d2(:,:),ek2(:,:),h1c2(:)
   complex(dpc),allocatable :: hm2(:,:,:)
   call rdmft_one_electron(orbs,dtset,dtfil,gsqcutf_eff,mpi_enreg,pawfgr,pawrad,pawtab,psps,h1d2,ek2,ew2,cp2,hm2)
   dmax = zero
   do jk=1,orbs%nkpt
     nb2 = orbs%nband(jk,1)
     dmax = max(dmax,maxval(abs(hm2(1:nb2,1:nb2,jk)-hmat(1:nb2,1:nb2,jk))))
   end do
   ABI_MALLOC(h1c2,(nbk))
   do jk=1,orbs%nkpt
     h1c2(offset(jk)+1:offset(jk)+orbs%nband(jk,1)) = h1d2(1:orbs%nband(jk,1),jk)
   end do
   dmaxd = maxval(abs(h1c2-h1c))
   e1 = energy_rot(orbs,h1c)
   e2 = energy_rot(orbs,h1c2)
   write(msg,'(a,es10.2,a,es10.2,2(a,es22.14),a)') ' RDMFT: check of the rotated h: max|h_rot - h_getghc| = ',dmax,&
     ' (diagonal ',dmaxd,'); E with the rotated h = ',e1,', with the recomputed h = ',e2,ch10
   call wrtout([ab_out,std_out],msg)
   ABI_FREE(h1c2)
   ABI_FREE(h1d2)
   ABI_FREE(ek2)
   ABI_FREE(hm2)
 end subroutine check_h

 function kinetic_cur() result(t)
   ! Kinetic energy per cell of the current orbitals and occupations nfull: sum_kn 2 w_k n_kn <kn|T|kn>, T diagonal in G
   real(dp) :: t
   integer :: jk,jb,ig
   real(dp) :: kg(3),tk
   t = zero
   do jk=1,orbs%nkpt
     do jb=1,orbs%nband(jk,1)
       if (nfull(offset(jk)+jb) < tol14) cycle
       tk = zero
       do ig=1,orbs%npw(jk)
         kg(:) = orbs%kpt(:,jk)+real(orbs%kg(:,ig,jk),dp)
         tk = tk+(orbs%cg(1,ig,jb,jk,1)**2+orbs%cg(2,ig,jb,jk,1)**2)*dot_product(kg,matmul(orbs%gmet,kg))
       end do
       t = t+two*orbs%wtk(jk)*nfull(offset(jk)+jb)*two*pi**2*tk
     end do
   end do
 end function kinetic_cur

 function energy_with(jm,km,bjv,bkv,bgx) result(e)
   ! Electronic energy of the functional with the given integrals (zero them to isolate a term): the energy is linear in J and K
   real(dp),intent(in) :: jm(:,:),km(:,:),bjv(:),bkv(:)
   type(rdmft_bg_t),intent(in) :: bgx
   real(dp) :: e
   integer :: jf
   real(dp),allocatable :: nx(:)
   ABI_MALLOC(nx,(nf))
   do jf=1,nf
     nx(jf) = nfull(ints%comp_of(jf))
   end do
   if (is_jk) then
     e = rdmft_energy_jk(dtset%rdmft_functional,orbs%nkpt,nf,nx,h1f,jm,km,alpha,bjv,bkv,bgx)
   else
     e = rdmft_energy_pnof(dtset%rdmft_functional,table,orbs%nkpt,nf,nx,h1f,jm,km,dtset%rdmft_selfint,bjv,bkv,bgx)
   end if
   ABI_FREE(nx)
 end function energy_with

 subroutine print_final_energies()
   ! Energy components of the optimized state next to those of the starting state (SD with the ground-state occupations)
   real(dp) :: ef(7),es(7),base,eh,ex,e1,esum
   real(dp),allocatable :: zm(:,:),zv(:)
   type(rdmft_bg_t) :: b0,bjt,bkt
   integer :: ic
   character(len=100) :: lbl(7)
   ABI_MALLOC(zm,(nf,nf))
   ABI_MALLOC(zv,(nf))
   zm = zero
   zv = zero
   b0 = bgt
   b0%sjo = zero; b0%sko = zero; b0%sjd = zero; b0%skd = zero
   bjt = bgt
   bjt%sko = zero; bjt%skd = zero
   bkt = bgt
   bkt%sjo = zero; bkt%sjd = zero
   base = energy_with(zm,zm,zv,zv,b0)
   eh = energy_with(ints%jint,zm,ints%bj,zv,bjt)-base
   ex = energy_with(zm,ints%kint,zv,ints%bk,bkt)-base
   e1 = base
   ef(1) = kinetic_cur()
   ef(2) = e1
   ef(3) = eh
   ef(4) = ex
   ef(5) = e_ewald
   ef(6) = e_corepsp
   ef(7) = e1+eh+ex+e_ewald+e_corepsp
   esum = ef(7)-energy_cur()
   es = [ekin_tot,eone,ehartree,eexch,e_ewald,e_corepsp,etot]
   lbl(1) = ' RDMFT: kinetic energy'
   lbl(2) = ' RDMFT: one-electron energy (T + Vloc + Vnl)'
   lbl(3) = ' RDMFT: Hartree energy'
   lbl(4) = ' RDMFT: exchange (functional) energy'
   lbl(5) = ' RDMFT: Ewald energy'
   lbl(6) = ' RDMFT: core pseudopotential energy'
   lbl(7) = ' RDMFT: total energy'
   write(msg,'(3a)') ch10,' RDMFT: energy components (Ha) of the optimized state, and change with respect to the starting',&
     ' state (SD, ground-state occupations)'
   call wrtout([ab_out,std_out],msg)
   do ic=1,7
     write(msg,'(a,t47,2es22.14)') trim(lbl(ic)),ef(ic),ef(ic)-es(ic)
     call wrtout([ab_out,std_out],msg)
   end do
   write(msg,'(a,t47,es22.14,a)') ' RDMFT: sum of the components - total energy',esum,ch10
   call wrtout([ab_out,std_out],msg)
   ABI_FREE(zm)
   ABI_FREE(zv)
 end subroutine print_final_energies

 function energy_cur() result(e)
   ! Total energy at the occupations nfull for the current integrals, one-electron elements and background sums
   real(dp) :: e
   integer :: jf
   real(dp),allocatable :: nx(:)
   ABI_MALLOC(nx,(nf))
   do jf=1,nf
     nx(jf) = nfull(ints%comp_of(jf))
   end do
   if (is_jk) then
     e = rdmft_energy_jk(dtset%rdmft_functional,orbs%nkpt,nf,nx,h1f,ints%jint,ints%kint,alpha,ints%bj,ints%bk,bgt)&
         +e_ewald+e_corepsp
   else
     e = rdmft_energy_pnof(dtset%rdmft_functional,table,orbs%nkpt,nf,nx,h1f,ints%jint,ints%kint,dtset%rdmft_selfint,&
                           ints%bj,ints%bk,bgt)+e_ewald+e_corepsp
   end if
   ABI_FREE(nx)
 end function energy_cur

 subroutine sync_h1f()
   ! One-electron elements and background sums from h1c and the integrals
   integer :: jf
   do jf=1,nf
     h1f(jf) = h1c(ints%comp_of(jf))
   end do
   bgt%sh = sum(h1c,mask=bg)
   bgt%sjo = ints%sjo; bgt%sko = ints%sko; bgt%sjd = ints%sjd; bgt%skd = ints%skd
 end subroutine sync_h1f

 subroutine eval_point(kpar,e_out,grad,gmax)
   ! Energy and complex orbital gradient (dE/dRe kappa_pq + i dE/dIm kappa_pq, p<q) of the current orbitals at the occupations
   ! nfull, from ONE sweep over the pair densities: the exchange integrals and the exchange part of the generalized Fock
   ! matrices are built in the same loop. The time-reversed k-point gives the same contribution as k, and the rotations at
   ! the time-reversal-invariant k-points are real.
   integer,intent(in) :: kpar(:)
   real(dp),intent(out) :: e_out,gmax
   complex(dpc),intent(out) :: grad(:)
   integer :: ib2,iq2,jk,ipar,pp,qq,nb2,nbmax,ncol,nrow
   real(dp) :: fac
   logical,allocatable :: relev(:)
   complex(dpc) :: dg
   complex(dpc),allocatable :: fk(:,:),fh(:,:),fchk(:,:)
   real(dp),allocatable :: acf(:,:),ccf(:,:),nux(:)
   type(rdmft_fockreq_t) :: req
   nbmax = maxval(orbs%nband(:,1))
   if (nmem_eval < 2) call report_mem('ORB: eval_point start')
   ABI_MALLOC(req%want,(orbs%nkpt))
   ABI_MALLOC(req%fx,(nbmax,nbmax,orbs%nkpt))
   req%want = (kpar >= 0)
   ABI_MALLOC(relev,(nbk))
   ! An orbital is relevant if its largest coupling (coefficient of the energy, in electrons) reaches the threshold; the
   ! couplings are symmetric, so the row of each band gives its relevance. The dropped couplings are set to zero.
   do jk=1,orbs%nkpt
     nb2 = orbs%nband(jk,1)
     ABI_MALLOC(acf,(nb2,nbk))
     ABI_MALLOC(ccf,(nb2,nbk))
     ABI_MALLOC(nux,(nb2))
     call build_coefs(jk,acf,ccf,nux)
     do ib2=1,nb2
       relev(offset(jk)+ib2) = (nk2*max(maxval(abs(acf(ib2,:))),maxval(abs(ccf(ib2,:)))) >= dtset%rdmft_occ_skip)
     end do
     ABI_FREE(acf)
     ABI_FREE(ccf)
     ABI_FREE(nux)
   end do
   if (nmem_eval < 2) call report_mem('ORB: eval_point after the relevance pass')
   ! The couplings are stored only between relevant bands (rows) and relevant orbitals (columns): the memory does not grow
   ! with the square of the number of orbitals of the complete basis
   ABI_MALLOC(req%rowmap,(nbmax,orbs%nkpt))
   ABI_MALLOC(req%colmap,(nbk))
   req%rowmap = 0
   req%colmap = 0
   ncol = 0
   do ib2=1,nbk
     if (relev(ib2)) then
       ncol = ncol+1
       req%colmap(ib2) = ncol
     end if
   end do
   nrow = 0
   do jk=1,orbs%nkpt
     do ib2=1,orbs%nband(jk,1)
       if (relev(offset(jk)+ib2)) then
         nrow = nrow+1
         req%rowmap(ib2,jk) = nrow
       end if
     end do
   end do
   ABI_MALLOC(req%ccoef,(max(nrow,1),max(ncol,1)))
   req%ccoef = zero
   do jk=1,orbs%nkpt
     if (kpar(jk) < 0) cycle
     nb2 = orbs%nband(jk,1)
     ABI_MALLOC(acf,(nb2,nbk))
     ABI_MALLOC(ccf,(nb2,nbk))
     ABI_MALLOC(nux,(nb2))
     call build_coefs(jk,acf,ccf,nux)
     do ib2=1,nb2
       if (req%rowmap(ib2,jk) == 0) cycle
       do iq2=1,nbk
         if (req%colmap(iq2) > 0) req%ccoef(req%rowmap(ib2,jk),req%colmap(iq2)) = ccf(ib2,iq2)
       end do
     end do
     ABI_FREE(acf)
     ABI_FREE(ccf)
     ABI_FREE(nux)
   end do
   nrel_last = count(relev)
   if (nmem_eval < 2) call report_mem('ORB: eval_point after the coupling arrays')
   t0 = abi_wtime()
   rdmft_mem_on = (nmem_eval < 2 .and. dtset%prtvol >= 5)
   call ints%compute(orbs,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,free,bg,&
                     dtset%rdmft_core_ints == 1,quiet=.true.,comm=xmpi_world,req=req,rel=relev,&
                     reduce=(dtset%rdmft_tr_reduce == 1))
   tm(1) = tm(1)+abi_wtime()-t0
   rdmft_mem_on = .false.
   nmem_eval = nmem_eval+1
   if (nmem_eval <= 2) call report_mem('ORB: eval_point after the integrals and Fock exchange')
   t0 = abi_wtime()
   call sync_h1f()
   e_out = energy_cur()
   tm(2) = tm(2)+abi_wtime()-t0
   grad = czero
   do jk=1,orbs%nkpt
     if (kpar(jk) < 0) cycle
     nb2 = orbs%nband(jk,1)
     ABI_MALLOC(acf,(nb2,nbk))
     ABI_MALLOC(ccf,(nb2,nbk))
     ABI_MALLOC(nux,(nb2))
     ABI_MALLOC(fk,(nb2,nb2))
     ABI_MALLOC(fh,(nb2,nb2))
     call build_coefs(jk,acf,ccf,nux)
     do ib2=1,nbk
       if (.not.relev(ib2)) acf(:,ib2) = zero
     end do
     do ib2=1,nb2
       if (.not.relev(offset(jk)+ib2)) acf(ib2,:) = zero
     end do
     t0 = abi_wtime()
     call rdmft_fock_hartree_k(orbs,jk,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,offset,nbk,acf,xmpi_world,fh)
     tm(3) = tm(3)+abi_wtime()-t0
     do qq=1,nb2
       do pp=1,nb2
         fk(pp,qq) = fh(pp,qq)+req%fx(pp,qq,jk)+nux(qq)*hmat(qq,pp,jk)
       end do
     end do
     if (dtset%prtvol >= 10 .and. jk == 1) then
       ! Validation of the fused (one sweep) Fock matrix against the stand-alone one
       ABI_MALLOC(fchk,(nb2,nb2))
       call build_fock(jk,fchk)
       write(msg,'(a,es10.2,a,es10.2)') '   RDMFT: Fock matrix of the first k-point, fused vs stand-alone: max difference = ',&
         maxval(abs(fk-fchk)),', max |F| = ',maxval(abs(fchk))
       call wrtout(std_out,msg)
       ABI_FREE(fchk)
     end if
     fac = one
     if (orbs%tr_k(jk) > 0 .and. orbs%tr_k(jk) /= jk) fac = two
     ipar = kpar(jk)
     do pp=1,nb2
       do qq=pp+1,nb2
         ipar = ipar+1
         dg = fk(pp,qq)-conjg(fk(qq,pp))
         grad(ipar) = fac*cmplx(two*real(dg,dp),-two*aimag(dg),kind=dpc)
         if (orbs%tr_k(jk) == jk) grad(ipar) = cmplx(real(grad(ipar),dp),zero,kind=dpc)
       end do
     end do
     ABI_FREE(acf)
     ABI_FREE(ccf)
     ABI_FREE(nux)
     ABI_FREE(fk)
     ABI_FREE(fh)
   end do
   gmax = zero
   if (size(grad) > 0) gmax = maxval(abs(grad))
   if (nmem_eval <= 2) call report_mem('ORB: eval_point after the Hartree part')
   ABI_FREE(relev)
   ABI_FREE(req%want)
   ABI_FREE(req%ccoef)
   ABI_FREE(req%rowmap)
   ABI_FREE(req%colmap)
   ABI_FREE(req%fx)
 end subroutine eval_point

 subroutine spot_check(kpar)
   ! Check of rdmft_occ_skip: the orbital gradient of the first independent k-point computed with all the partners and with the
   ! skipping (cost of the order of one sweep over the partners of one k-point)
   integer,intent(in) :: kpar(:)
   integer :: ib2,jk,k0,nb2,pp,qq
   real(dp) :: fac,gall,gskip,gdiff
   complex(dpc) :: dg,ga,gs
   logical,allocatable :: relv(:)
   real(dp),allocatable :: acf(:,:),ccf(:,:),acs(:,:),ccs(:,:),nux(:)
   complex(dpc),allocatable :: fa(:,:),fs(:,:)
   k0 = 0
   do jk=1,orbs%nkpt
     if (kpar(jk) >= 0) then
       k0 = jk
       exit
     end if
   end do
   if (k0 == 0) return
   ABI_MALLOC(relv,(nbk))
   do jk=1,orbs%nkpt
     nb2 = orbs%nband(jk,1)
     ABI_MALLOC(acf,(nb2,nbk))
     ABI_MALLOC(ccf,(nb2,nbk))
     ABI_MALLOC(nux,(nb2))
     call build_coefs(jk,acf,ccf,nux)
     do ib2=1,nb2
       relv(offset(jk)+ib2) = (nk2*max(maxval(abs(acf(ib2,:))),maxval(abs(ccf(ib2,:)))) >= dtset%rdmft_occ_skip)
     end do
     ABI_FREE(acf)
     ABI_FREE(ccf)
     ABI_FREE(nux)
   end do
   nb2 = orbs%nband(k0,1)
   ABI_MALLOC(acf,(nb2,nbk))
   ABI_MALLOC(ccf,(nb2,nbk))
   ABI_MALLOC(acs,(nb2,nbk))
   ABI_MALLOC(ccs,(nb2,nbk))
   ABI_MALLOC(nux,(nb2))
   ABI_MALLOC(fa,(nb2,nb2))
   ABI_MALLOC(fs,(nb2,nb2))
   call build_coefs(k0,acf,ccf,nux)
   acs = acf
   ccs = ccf
   do ib2=1,nbk
     if (.not.relv(ib2)) then
       acs(:,ib2) = zero
       ccs(:,ib2) = zero
     end if
   end do
   do ib2=1,nb2
     if (.not.relv(offset(k0)+ib2)) then
       acs(ib2,:) = zero
       ccs(ib2,:) = zero
     end if
   end do
   call rdmft_gen_fock_k(orbs,k0,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,&
                         offset,nbk,hmat(1:nb2,1:nb2,k0),nux,acf,ccf,xmpi_world,fa)
   call rdmft_gen_fock_k(orbs,k0,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,&
                         offset,nbk,hmat(1:nb2,1:nb2,k0),nux,acs,ccs,xmpi_world,fs)
   fac = one
   if (orbs%tr_k(k0) > 0 .and. orbs%tr_k(k0) /= k0) fac = two
   gall = zero
   gskip = zero
   gdiff = zero
   do pp=1,nb2
     do qq=pp+1,nb2
       dg = fa(pp,qq)-conjg(fa(qq,pp))
       ga = fac*cmplx(two*real(dg,dp),-two*aimag(dg),kind=dpc)
       dg = fs(pp,qq)-conjg(fs(qq,pp))
       gs = fac*cmplx(two*real(dg,dp),-two*aimag(dg),kind=dpc)
       if (orbs%tr_k(k0) == k0) then
         ga = cmplx(real(ga,dp),zero,kind=dpc)
         gs = cmplx(real(gs,dp),zero,kind=dpc)
       end if
       gall = max(gall,abs(ga))
       gskip = max(gskip,abs(gs))
       gdiff = max(gdiff,abs(ga-gs))
     end do
   end do
   write(msg,'(a,i0,a,es10.2,a,es10.2,a,es10.2,a,i0,a,i0,a)') ' RDMFT: check of rdmft_occ_skip at the k-point ',k0,&
     ': max|g| with all the pairs = ',gall,', with the skipping = ',gskip,', largest difference = ',gdiff,&
     ' (orbitals kept in the sums: ',count(relv),' of ',nbk,')'
   call wrtout([ab_out,std_out],msg)
   ABI_FREE(relv)
   ABI_FREE(acf)
   ABI_FREE(ccf)
   ABI_FREE(acs)
   ABI_FREE(ccs)
   ABI_FREE(nux)
   ABI_FREE(fa)
   ABI_FREE(fs)
 end subroutine spot_check

 subroutine orbital_macro(e_start)
   ! Macro-iterations: ADAM rotation of the orbitals at fixed occupations, then optimization of the occupations
   real(dp),intent(in) :: e_start
   integer :: imac,icall,icall_max,extra,jk,nadam,npar,nstall,nit2,mi
   integer,allocatable :: kpar(:)
   real(dp) :: e_old,e_new,e_prev,e_best,e_run0,lr,gmax,gmax_b,de
   logical :: gconv,improved,restart,conv2
   character(len=120) :: stat
   complex(dpc),allocatable :: grad(:),dx(:),grad_b(:)
   complex(dpc),allocatable :: hmat_b(:,:,:)
   real(dp),allocatable :: h1c_b(:)
   type(rdmft_adam_t) :: adam
   type(rdmft_orbs_t) :: orbs_b

   ABI_MALLOC(kpar,(orbs%nkpt))
   npar = 0
   do jk=1,orbs%nkpt
     ! independent k-points: -k is represented by k (time reversal)
     if (orbs%tr_k(jk) > 0 .and. orbs%tr_k(jk) < jk) then
       kpar(jk) = -1
     else
       kpar(jk) = npar
       npar = npar+orbs%nband(jk,1)*(orbs%nband(jk,1)-1)/2
     end if
   end do
   ABI_MALLOC(grad,(npar))
   ABI_MALLOC(grad_b,(npar))
   ABI_MALLOC(dx,(npar))
   ABI_MALLOC(h1c_b,(nbk))
   ABI_MALLOC(hmat_b,(size(hmat,1),size(hmat,2),size(hmat,3)))
   call adam%init(npar,dtset%rdmft_adam_betam,dtset%rdmft_adam_betav)
   call report_mem('ORB: ADAM arrays allocated')
   write(msg,'(2a,i0,a,es9.2,a,es9.2,a,es9.2,a)') ch10,' RDMFT: ADAM optimization of the orbitals (',npar,&
     ' complex rotation parameters); energy tolerance ',dtset%rdmft_orb_etol,', gradient tolerance ',dtset%rdmft_orb_gtol,&
     ', learning rate ',dtset%rdmft_adam_lr,ch10
   call wrtout([ab_out,std_out],msg)
   write(msg,'(2a)') ' Orbital optimization, macro-iterations (ADAM rotations of the orbitals, then occupation numbers):',ch10//&
     '  MACRO  macro        Etotal(Ha)      deltaE(Ha)    max|g|  ADAM-its  occ-its  status'
   call wrtout([ab_out,std_out],msg)
   lr = dtset%rdmft_adam_lr
   extra = 0
   e_old = e_start
   nstall = 0
   mi = 0
   do imac=1,dtset%rdmft_orb_macro
     icall_max = dtset%rdmft_adam_iter+extra
     if (icall_max > dtset%rdmft_adam_maxit) then
       icall_max = dtset%rdmft_adam_iter
       extra = 0
     end if
     call adam%clean()
     ! energy and gradient at the start of the run (new occupations)
     call eval_point(kpar,e_best,grad,gmax)
     e_run0 = e_best
     e_prev = e_best
     write(msg,'(a,i0,a,i0,a,es9.2,a,i0,a,i0,a)') ch10//' ADAM orbital rotations, macro-iteration ',imac,' (at most ',&
       icall_max,' iterations, learning rate ',lr,'); orbitals kept in the pair sums: ',nrel_last,' of ',nbk,' (rdmft_occ_skip)'
     call wrtout([ab_out,std_out],msg)
     write(msg,'(a)') '  ORB  macro  iter        Etotal(Ha)      deltaE(Ha)    max|g|'
     call wrtout([ab_out,std_out],msg)
     write(msg,'(a,2i6,f20.12,a,es10.2)') '  ORB',imac,0,e_best,'                ',gmax
     call wrtout([ab_out,std_out],msg)
     nadam = 0
     orbs_b = orbs
     h1c_b = h1c
     hmat_b = hmat
     grad_b = grad
     gmax_b = gmax
     gconv = .false.
     improved = .false.
     de = zero
     do icall=1,icall_max
       if (gmax < dtset%rdmft_orb_gtol) then
         gconv = .true.
         exit
       end if
       t0 = abi_wtime()
       call adam%step(grad,lr,dx)
       tm(6) = tm(6)+abi_wtime()-t0
       t0 = abi_wtime()
       call rotate_all(kpar,dx)
       tm(4) = tm(4)+abi_wtime()-t0
       call eval_point(kpar,e_new,grad,gmax)
       nadam = icall
       write(msg,'(a,2i6,f20.12,es16.3,es10.2)') '  ORB',imac,icall,e_new,e_new-e_prev,gmax
       call wrtout([ab_out,std_out],msg)
       if (e_new < e_best-tol14) then
         e_best = e_new
         improved = .true.
         orbs_b = orbs
         h1c_b = h1c
         hmat_b = hmat
         grad_b = grad
         gmax_b = gmax
       end if
       de = e_new-e_prev
       e_prev = e_new
       if (icall > 1 .and. abs(de) < tol14) exit
     end do
     ! A run that ends on the gradient criterion keeps its last point (stationary), even if its energy equals the best one
     if (gconv .and. icall > 1) then
       e_best = e_prev
       orbs_b = orbs
       h1c_b = h1c
       hmat_b = hmat
       grad_b = grad
       gmax_b = gmax
     end if
     ! No improvement although the energy still moves: restart from the best point with a smaller learning rate
     ! A run that ends without converging the gradient and without lowering the energy by more than the tolerance is limited by
     ! the noise of ADAM at this learning rate: the learning rate is reduced for the next run
     restart = (.not.gconv) .and. ((.not.improved .and. abs(de) > dtset%rdmft_orb_etol) .or. &
                                   (e_run0-e_best < dtset%rdmft_orb_etol))
     if (restart) then
       lr = lr*dtset%rdmft_adam_rate
       extra = extra+dtset%rdmft_adam_incr
     else if (gconv) then
       lr = dtset%rdmft_adam_lr
       extra = 0
     end if
     if (imac == 1) call report_mem('ORB: macro-iteration 1, after the ADAM run')
     orbs = orbs_b
     h1c = h1c_b
     hmat = hmat_b
     gmax = gmax_b
     ! all the pairs of orbitals for the optimization of the occupations
     t0 = abi_wtime()
     call ints%compute(orbs,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,free,bg,&
                       dtset%rdmft_core_ints == 1,quiet=.true.,comm=xmpi_world,reduce=(dtset%rdmft_tr_reduce == 1))
     tm(7) = tm(7)+abi_wtime()-t0
     if (imac == 1) call report_mem('ORB: macro-iteration 1, full integrals for the occupations')
     call sync_h1f()
     ! Occupations for the new orbitals
     t0 = abi_wtime()
     call optimize_occupations(nopt,e_new,conv2,nit2,.false.)
     tm(5) = tm(5)+abi_wtime()-t0
     call assemble_nfull()
     e_new = e_new+e_ewald+e_corepsp
     stat = 'ADAM ran'
     if (gconv) stat = 'gradient converged'
     if (restart) stat = trim(stat)//', restart with a smaller learning rate'
     if (.not.conv2) stat = trim(stat)//'; occupations not converged'
     write(msg,'(a,i6,f20.12,es16.3,es10.2,2i10,2a)') '  MACRO',imac,e_new,e_new-e_old,gmax,nadam,nit2,'  ',trim(stat)
     call wrtout([ab_out,std_out],msg)
     nstall = merge(nstall+1,0,abs(e_new-e_old) < dtset%rdmft_orb_etol)
     mi = imac
     if (abs(e_new-e_old) < dtset%rdmft_orb_etol .and. gmax < dtset%rdmft_orb_gtol .and. (.not.restart .or. gconv)) then
       e_old = e_new
       exit
     end if
     if (lr < 1.0e-6_dp*dtset%rdmft_adam_lr) then
       write(msg,'(a,es9.2,a)') ' RDMFT: the learning rate of ADAM reached its floor (',lr,'); the orbital optimization is stopped'
       call wrtout([ab_out,std_out],msg)
       e_old = e_new
       exit
     end if
     if (nstall >= 10) then
       e_old = e_new
       exit
     end if
     e_old = e_new
   end do
   write(msg,'(a,i0,a)') ' RDMFT: orbital optimization finished after ',mi,' macro-iterations'
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,t52,a,es22.14,a)') ' RDMFT: total energy after the orbital optimization',' = ',e_old,ch10
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,5(a,f9.2))') ' RDMFT: timings (s):',' ints+Fock exchange =',tm(1),', energy =',tm(2),', Hartree Fock =',tm(3),&
     ', rotation =',tm(4),', occupations =',tm(5)
   call wrtout([ab_out,std_out],msg)
   ! Parallel efficiency report: wall times of the phases (max and min over the processes)
   tm(8) = abi_wtime()-t_start
   tint_max = rdmft_tint
   tint_min = -rdmft_tint
   call xmpi_max(tint_max,xmpi_world,ierr)
   call xmpi_max(tint_min,xmpi_world,ierr)
   tint_min = -tint_min
   write(msg,'(a,i0,a)') ' RDMFT: parallel timings (s) on ',xmpi_comm_size(xmpi_world),' processes (max over processes)'
   call wrtout([ab_out,std_out],msg)
   write(msg,'(2(a,f10.2),a,f10.2)') '   ADAM step =',tm(6),'   full ints for occupations =',tm(7),'   total RDMFT =',tm(8)
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,6f10.2)') '   ints: dens, sumdens, J, Ksweep, Kwait, Ksum (max) =',tint_max
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,6f10.2)') '   ints: dens, sumdens, J, Ksweep, Kwait, Ksum (min) =',tint_min
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,5f10.2,f10.2,es11.3)') '   K sweep profile, rank 0: toreal(k), toreal(k2), product, FFT, K value, Fock ZGEMM, n_FFT = ',&
     rdmft_tprof(1:6),rdmft_tprof(7)
   call wrtout([ab_out,std_out],msg)
   write(msg,'(a,4f10.2,f10.0)') '   PNOF evaluation profile: occupations+Jacobian, energy, gradient, contraction, n_eval =',rdmft_tpnof
   call wrtout([ab_out,std_out],msg)
   if (dtset%rdmft_occ_skip > zero) call spot_check(kpar)
   if (dtset%prtvol >= 10) call check_h()
   call adam%free()
   call orbs_b%free()
   ABI_FREE(kpar)
   ABI_FREE(grad)
   ABI_FREE(grad_b)
   ABI_FREE(dx)
   ABI_FREE(h1c_b)
   ABI_FREE(hmat_b)
 end subroutine orbital_macro

 subroutine compute_density(occ_all,rho)
   ! Electron density rho(r) = 2 sum_k w_k sum_n n_kn |u_kn(r)|^2 / Omega (electrons/bohr^3) on the dense FFT mesh,
   ! with n_kn the occupations per spin (between 0 and 1) of the composite index offset(k)+n
   real(dp),intent(in) :: occ_all(:)
   real(dp),intent(out) :: rho(:)
   integer :: jk,jb
   real(dp),allocatable :: ur(:,:)
   ABI_MALLOC(ur,(2,nfftf))
   rho = zero
   do jk=1,orbs%nkpt
     do jb=1,orbs%nband(jk,1)
       if (occ_all(offset(jk)+jb) < tol14) cycle
       call orb_to_real_one(orbs,jk,jb,ngfftf,mpi_seq,ur)
       rho(:) = rho(:)+two*orbs%wtk(jk)*occ_all(offset(jk)+jb)*(ur(1,:)**2+ur(2,:)**2)/orbs%ucvol
     end do
   end do
   ABI_FREE(ur)
 end subroutine compute_density

 subroutine print_occupations()
   ! Optimized occupation numbers k-point by k-point (spin summed, as in the GW natural orbitals), the electron count and the
   ! integrated density
   integer :: jk,jb,nb2,m,ncomp(3)
   real(dp) :: tk,tot,dv,drho,qint,rmax,rmin,idem
   real(dp),allocatable :: rho_end(:)
   character(len=100) :: ln(8)
   write(msg,'(2a)') ch10,' RDMFT: optimized occupation numbers (spin summed, between 0 and 2), k-point by k-point:'
   call wrtout([ab_out,std_out],msg)
   tot = zero
   idem = zero
   do jk=1,orbs%nkpt
     nb2 = orbs%nband(jk,1)
     write(msg,'(a,3f10.5)') ' Occupations at the k-point (reduced coordinates): ',orbs%kpt(:,jk)
     call wrtout([ab_out,std_out],msg)
     do jb=1,nb2,10
       m = min(10,nb2-jb+1)
       write(msg,'(f11.5,*(f10.5))') two*nfull(offset(jk)+jb:offset(jk)+jb+m-1)
       call wrtout([ab_out,std_out],msg)
     end do
     tk = two*sum(nfull(offset(jk)+1:offset(jk)+nb2))
     tot = tot+orbs%wtk(jk)*tk
     idem = idem+orbs%wtk(jk)*sum(nfull(offset(jk)+1:offset(jk)+nb2)*(one-nfull(offset(jk)+1:offset(jk)+nb2)))
     write(msg,'(a22,i5,a3,i5,a21,f10.5,a)') ' Total occ. from band ',1,' to',nb2,' at current k-point: ',tk,ch10
     call wrtout([ab_out,std_out],msg)
   end do
   ! Electron density of the natural orbitals: integral, change with respect to the ground state, extrema
   ABI_MALLOC(rho_end,(nfftf))
   call compute_density(nfull,rho_end)
   dv = orbs%ucvol/real(nfftf,dp)
   qint = sum(rho_end)*dv
   drho = sum(abs(rho_end-rho_ref))*dv
   rmax = maxval(rho_end)
   rmin = minval(rho_end)
   ABI_FREE(rho_end)
   ncomp(1) = count(nfull > one-1.0e-6_dp)
   ncomp(3) = count(nfull < 1.0e-6_dp)
   ncomp(2) = nbk-ncomp(1)-ncomp(3)
   write(ln(1),'(a,t52,a,f22.10)') ' RDMFT: total averaged occupation, all k-points',' = ',tot
   write(ln(2),'(a,t52,a,f22.10)') ' RDMFT: number of electrons (nelect)',' = ',orbs%nelect
   write(ln(3),'(a,t52,a,f22.10)') ' RDMFT: integrated electron density',' = ',qint
   write(ln(4),'(a,t52,a,f22.10)') ' RDMFT: integral of |rho - rho(ground state)|',' = ',drho
   write(ln(5),'(a,t52,a,es22.10)') ' RDMFT: maximum of the density (e/bohr^3)',' = ',rmax
   write(ln(6),'(a,t52,a,es22.10)') ' RDMFT: minimum of the density (e/bohr^3)',' = ',rmin
   write(ln(7),'(a,t52,a,es22.10)') ' RDMFT: sum_p n_p (1 - n_p), 0 for one determinant',' = ',idem
   write(ln(8),'(a,t52,a,3(i0,1x))') ' RDMFT: orbitals with n=1 / fractional / n=0',' = ',ncomp
   write(msg,'(4a)') trim(ln(1))//ch10,trim(ln(2))//ch10,trim(ln(3))//ch10,trim(ln(4))
   call wrtout([ab_out,std_out],msg)
   write(msg,'(4a)') trim(ln(5))//ch10,trim(ln(6))//ch10,trim(ln(7))//ch10,trim(ln(8))
   call wrtout([ab_out,std_out],msg)
 end subroutine print_occupations

 subroutine write_natural_orbitals()
   ! WFK file (natural orbitals, occupations 2 n between 0 and 2, eigenvalues = diagonal of the one-electron Hamiltonian)
   ! and DEN file (density of the natural orbitals) of the final state, written by the master as in the GW 1-RDM driver
   integer :: jk,nb2
   real(dp),allocatable :: eig_w(:,:),occ_w(:,:),rho_w(:,:)
   type(hdr_type) :: hdr_w
   type(crystal_t) :: cryst_w
   character(len=fnlen) :: wfk_out,den_out
   ABI_MALLOC(eig_w,(orbs%mband,orbs%nkpt))
   ABI_MALLOC(occ_w,(orbs%mband,orbs%nkpt))
   eig_w = zero
   occ_w = zero
   do jk=1,orbs%nkpt
     nb2 = orbs%nband(jk,1)
     eig_w(1:nb2,jk) = h1c(offset(jk)+1:offset(jk)+nb2)
     occ_w(1:nb2,jk) = two*nfull(offset(jk)+1:offset(jk)+nb2)
   end do
   ABI_MALLOC(rho_w,(nfftf,1))
   if (dtset%prtden > 0) call compute_density(nfull,rho_w(:,1))
   wfk_out = dtfil%fnameabo_wfk
   den_out = dtfil%fnameabo_den
   if (dtset%iomode == IO_MODE_ETSF) then
     wfk_out = nctk_ncify(wfk_out)
     den_out = nctk_ncify(den_out)
   end if
   if (xmpi_comm_rank(xmpi_world) == master) then
     call orbs%make_hdr(hdr_w,occ_w)
     if (dtset%prtwf > 0) then
       call orbs%write_wfk(dtfil%fnameabo_wfk,dtset%iomode,hdr_w,eig_w,occ_w)
       write(msg,'(2a)') ch10,' RDMFT: natural orbitals written to '//trim(wfk_out)//&
         ' (occupations 2n; the eigenvalues are the diagonal elements of the one-electron Hamiltonian)'
       call wrtout([ab_out,std_out],msg)
     end if
     if (dtset%prtden > 0) then
       hdr_w%ngfft(1:3) = ngfftf(1:3)
       cryst_w = hdr_w%get_crystal()
       call fftdatar_write("density",dtfil%fnameabo_den,dtset%iomode,hdr_w,cryst_w,ngfftf,1,nfftf,1,rho_w,mpi_seq)
       call cryst_w%free()
       write(msg,'(2a)') ch10,' RDMFT: density of the natural orbitals written to '//trim(den_out)
       call wrtout([ab_out,std_out],msg)
     end if
     call hdr_w%free()
   end if
   ABI_FREE(eig_w)
   ABI_FREE(occ_w)
   ABI_FREE(rho_w)
 end subroutine write_natural_orbitals

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
                   dtset%rdmft_core_ints == 1,quiet=.true.,reduce=(dtset%rdmft_tr_reduce == 1))
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


 subroutine build_coefs(ikx,acf,ccf,nux)
   ! Coefficients of the energy for the bands of the k-point ikx at the occupations nfull:
   ! acf = a_jq, ccf = c_jq (divided by N_k^2), nux = 2 n_j / N_k
   integer,intent(in) :: ikx
   real(dp),intent(out) :: acf(:,:),ccf(:,:),nux(:)
   integer :: jb2,oj,qq,nbx,jsub,jf
   real(dp) :: ac,cc
   real(dp),allocatable :: hp_c(:)
   integer,allocatable :: sub_c(:)
   logical,allocatable :: prin_c(:)
   nbx = orbs%nband(ikx,1)
   ABI_MALLOC(sub_c,(nbk))
   ABI_MALLOC(prin_c,(nbk))
   ABI_MALLOC(hp_c,(nbk))
   sub_c = 0; prin_c = .false.; hp_c = one
   if (is_pnof) then
     ! subspaces of the table; the background orbitals are trivial subspaces (n=1)
     do jf=1,table%ng
       oj = ints%comp_of(table%gem_orb(jf))
       sub_c(oj) = table%gem_sub(jf)
       prin_c(oj) = table%gem_principal(jf)
     end do
     do jf=1,table%ng
       oj = ints%comp_of(table%gem_orb(jf))
       if (prin_c(oj)) then
         do jsub=1,table%ng
           if (sub_c(ints%comp_of(table%gem_orb(jsub))) == sub_c(oj)) hp_c(ints%comp_of(table%gem_orb(jsub))) = one-nfull(oj)
         end do
       end if
     end do
     jsub = table%nsub
     do qq=1,nbk
       if (bg(qq)) then
         jsub = jsub+1
         sub_c(qq) = jsub; prin_c(qq) = .true.; hp_c(qq) = zero
       end if
     end do
   end if
   do jb2=1,nbx
     oj = offset(ikx)+jb2
     nux(jb2) = two*nfull(oj)/orbs%nkpt
     do qq=1,nbk
       if (is_jk) then
         call rdmft_coef_jk(dtset%rdmft_functional,nfull(oj),nfull(qq),qq==oj,alpha,ac,cc)
       else
         call rdmft_coef_pnof(dtset%rdmft_functional,nfull(oj),nfull(qq),qq==oj,sub_c(oj),sub_c(qq),prin_c(oj),prin_c(qq),&
                              hp_c(oj),hp_c(qq),dtset%rdmft_selfint,ac,cc)
       end if
       acf(jb2,qq) = ac/nk2
       ccf(jb2,qq) = cc/nk2
     end do
   end do
   ABI_FREE(sub_c)
   ABI_FREE(prin_c)
   ABI_FREE(hp_c)
 end subroutine build_coefs

 subroutine build_fock(ikx,fx)
   ! Generalized Fock matrix (see m_rdmft_fock) of the k-point ikx at the occupations nfull
   integer,intent(in) :: ikx
   complex(dpc),intent(out) :: fx(:,:)
   integer :: nbx
   real(dp),allocatable :: acf(:,:),ccf(:,:),nux(:)
   nbx = orbs%nband(ikx,1)
   ABI_MALLOC(acf,(nbx,nbk))
   ABI_MALLOC(ccf,(nbx,nbk))
   ABI_MALLOC(nux,(nbx))
   call build_coefs(ikx,acf,ccf,nux)
   call rdmft_gen_fock_k(orbs,ikx,ngfftf,gsqcutf_eff,dtset%icutcoul,dtset%rcut,dtset%fock_icutcoul,dtset%fock_rcut,&
                         offset,nbk,hmat(1:nbx,1:nbx,ikx),nux,acf,ccf,xmpi_world,fx)
   ABI_FREE(acf)
   ABI_FREE(ccf)
   ABI_FREE(nux)
 end subroutine build_fock

 subroutine check_rotations(e_ref)
   ! Finite-difference gradient dE/dtheta and curvature of the energy with respect to the rotation of two states of the
   ! same k-point, at fixed occupations
   real(dp),intent(in) :: e_ref
   integer :: ik0,ip_b,iq_b,kp,nb0,ndir
   integer :: ib_w,ik_w
   real(dp) :: devtr,e0,ep,em,gall,gan,step
   complex(dpc) :: dg
   complex(dpc),allocatable :: f0(:,:),fp(:,:)
   logical :: has_p
   logical :: imag
   character(len=500) :: cmsg
   ik0 = dtset%rdmft_rot_k
   nb0 = orbs%nband(ik0,1)
   step = dtset%rdmft_rot_step
   e0 = energy_rot(orbs,h1c)
   ! Analytic gradient from the generalized Fock matrices of the k-point (and of -k, rotated with the conjugate matrix)
   ABI_MALLOC(f0,(nb0,nb0))
   call build_fock(ik0,f0)
   kp = orbs%tr_k(ik0)
   has_p = (kp > 0 .and. kp /= ik0)
   ABI_MALLOC(fp,(nb0,nb0))
   fp = czero
   if (has_p) call build_fock(kp,fp)
   gall = zero
   do ip_b=1,nb0
     do iq_b=ip_b+1,nb0
       gall = max(gall,two*abs(f0(ip_b,iq_b)-conjg(f0(iq_b,ip_b))))
     end do
   end do
   write(cmsg,'(a,i0,a,es10.2)') ' RDMFT: rotation check, k-point ',ik0,&
     ': largest |gradient| over all the pairs (real direction) = ',gall
   if (has_p) then
     call orbs%tr_alignment(devtr,ik_w,ib_w)
     write(cmsg,'(a,es10.2,a,i0,a,i0,a,es10.2,a,es10.2)') ' RDMFT: rotation check: TR alignment of the orbitals ',devtr,&
       ' (k-point ',ik_w,', band ',ib_w,'); max|F(-k) - conjg(F(k))| = ',maxval(abs(fp-conjg(f0))),', max|F| = ',maxval(abs(f0))
     call wrtout([ab_out,std_out],cmsg)
     ! occupations and one-electron matrices of k and -k
     write(cmsg,'(a,es10.2,a,es10.2)') '   max|n(-k) - n(k)| = ',maxval(abs(nfull(offset(kp)+1:offset(kp)+nb0)-&
       nfull(offset(ik0)+1:offset(ik0)+nb0))),'; max|h(-k) - conjg(h(k))| = ',maxval(abs(hmat(1:nb0,1:nb0,kp)-&
       conjg(hmat(1:nb0,1:nb0,ik0))))
     call wrtout([ab_out,std_out],cmsg)
   end if
   call wrtout([ab_out,std_out],cmsg)
   write(cmsg,'(2a,es20.12,a,es10.2,3a,i0,a,es9.2,a)') ch10,' RDMFT: rotation check, E(0) = ',e0,&
     '  (E(0) - optimized energy =',e0-e_ref,')',ch10,' k-point ',ik0,', step ',step,' radians (central differences)'
   call wrtout([ab_out,std_out],cmsg)
   write(cmsg,'(2a)') '   band p  band q    n_p       n_q     direction      dE/dtheta (num.)   dE/dtheta (analytic)  ',&
     '  d2E/dtheta2'
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
         ! analytic: real direction 2 Re (F_pq - conjg F_qp); imaginary direction -2 Im (F_pq - conjg F_qp);
         ! at -k the rotation matrix is conjugated, which changes the sign of the imaginary direction
         dg = f0(ip_b,iq_b)-conjg(f0(iq_b,ip_b))
         if (.not.imag) then
           gan = two*real(dg,dp)
           if (has_p) gan = gan+two*real(fp(ip_b,iq_b)-conjg(fp(iq_b,ip_b)),dp)
         else
           gan = -two*aimag(dg)
           if (has_p) gan = gan+two*aimag(fp(ip_b,iq_b)-conjg(fp(iq_b,ip_b)))
         end if
         write(cmsg,'(2i8,2f10.6,5x,a,3es18.8)') ip_b,iq_b,nfull(offset(ik0)+ip_b),nfull(offset(ik0)+iq_b),&
           merge('imaginary','real     ',imag),(ep-em)/(two*step),gan,(ep-two*e0+em)/step**2
         call wrtout([ab_out,std_out],cmsg)
       end do
     end do
   end do
   ABI_FREE(f0)
   ABI_FREE(fp)
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
