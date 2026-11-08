!!****m* ABINIT/m_rdmft_hamiltonian
!! NAME
!!  m_rdmft_hamiltonian
!!
!! FUNCTION
!!  One-electron part of the RDMFT energy for fixed ground-state orbitals.
!!  The core Hamiltonian h = T + V_loc,psp + V_nl is applied to the orbitals with ABINIT's getghc,
!!  and the diagonal elements h_pp = <p|h|p> are returned, together with the Ewald energy and the
!!  pseudopotential core energy (the constant terms of the total energy).
!!  Present limitations: norm-conserving pseudopotentials, nspinor=1, nsppol=1, useylm=0, no band-FFT parallelism.
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

module m_rdmft_hamiltonian

 use defs_basis
 use defs_abitypes,      only : MPI_type
 use defs_datatypes,     only : pseudopotential_type
 use defs_wvltypes,      only : wvl_data, nullify_wvl_data
 use m_abicore
 use m_errors
 use m_xmpi

 use m_dtfil,            only : datafiles_type
 use m_dtset,            only : dataset_type
 use m_ewald,            only : ewald
 use m_geometry,         only : metric
 use m_getghc,           only : multithreaded_getghc
 use m_hamiltonian,      only : gs_hamiltonian_type, gspot_transgrid_and_pack
 use m_kg,               only : getph, mkkin
 use m_mkffnl,           only : mkffnl
 use m_mklocl,           only : mklocl
 use m_pawang,           only : pawang_type
 use m_pawcprj,          only : pawcprj_type
 use m_pawfgr,           only : pawfgr_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pspini,           only : pspini
 use m_rdmft_orbitals,   only : rdmft_orbs_t
 use m_specialmsg,       only : wrtout

 implicit none

 private
!!***

 public :: rdmft_one_electron
!!***

contains
!!***

!!****f* m_rdmft_hamiltonian/rdmft_one_electron
!! NAME
!!  rdmft_one_electron
!!
!! FUNCTION
!!  Compute h_pp = <p|T + V_loc,psp + V_nl|p> for all the orbitals, the Ewald energy and the core energy.
!!
!! INPUTS
!!  orbs <type(rdmft_orbs_t)> = ground-state orbitals
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  dtfil <type(datafiles_type)> = infos about file names
!!  gsqcut = cutoff on G^2 for the FFT mesh of the density
!!  mpi_enreg <type(MPI_type)> = MPI information
!!  pawfgr <type(pawfgr_type)> = FFT meshes
!!  pawrad(ntypat*usepaw), pawtab(ntypat*usepaw) = PAW data (empty for norm-conserving pseudopotentials)
!!  psps <type(pseudopotential_type)> = pseudopotentials (read here by pspini)
!!
!! OUTPUT
!!  h1diag(mband,nkpt) = diagonal elements of the core Hamiltonian (Hartree)
!!  ekin(mband,nkpt) = kinetic part of h1diag
!!  e_ewald = Ewald energy per unit cell (Hartree)
!!  e_corepsp = pseudopotential core energy per unit cell (Hartree)
!!  hmat(mband,mband,nkpt) = (optional) matrix elements <i|h|j> of the core Hamiltonian between the bands of each k-point
!!
!! SOURCE
subroutine rdmft_one_electron(orbs,dtset,dtfil,gsqcut,mpi_enreg,pawfgr,pawrad,pawtab,psps,&
                              h1diag,ekin,e_ewald,e_corepsp,hmat)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 type(dataset_type),intent(inout) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(inout) :: psps
 real(dp),intent(out) :: e_ewald,e_corepsp
 !arrays
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 real(dp),allocatable,intent(out) :: h1diag(:,:),ekin(:,:)
 complex(dpc),allocatable,optional,intent(out) :: hmat(:,:,:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: cpopt=0,sij_opt=0,tim_getghc=5,type_calc=0,ider=0,idir=0,dimffnl=1,nkpg=0,option_locl=1
 integer :: ib,ierr,ig,ik,iatom,indx,itypat,jb,n4,n5,n6,nb,nfft,npw_k,psp_gencond
 real(dp) :: ecore,eei,ucvol,zion
 character(len=500) :: msg
 type(gs_hamiltonian_type) :: ham_k
 type(wvl_data) :: wvl
 !arrays
 integer :: ngfft(18)
 integer,allocatable :: atindx(:),atindx1(:),nattyp(:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),dyfr_dum(3,3,dtset%natom),gr_dum(3,dtset%natom)
 real(dp) :: grewtn(3,dtset%natom),dummy6(6)
 real(dp),allocatable :: cwavef(:,:),ffnl(:,:,:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),rhog_dum(:,:),rhor_dum(:,:),vlocal(:,:,:,:),vpsp(:),vtrial(:,:,:)
 real(dp),allocatable :: xred(:,:),ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *************************************************************************

 if (psps%usepaw /= 0 .or. dtset%nspinor /= 1 .or. dtset%nsppol /= 1 .or. psps%useylm /= 0 .or. dtset%paral_kgb /= 0) then
   ABI_ERROR('RDMFT: norm-conserving pseudopotentials, nspinor=1, nsppol=1, useylm=0 and paral_kgb=0 are required.')
 end if

 ngfft = dtset%ngfft
 nfft = product(ngfft(1:3))
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)
 rprimd = orbs%rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 ABI_MALLOC(xred,(3,dtset%natom))
 xred(:,:) = dtset%xred_orig(:,:,1)

 !Read the pseudopotential files
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcut,gsqcut,pawrad,pawtab,psps,rprimd,comm_mpi=xmpi_world)
 if (dtset%icoulomb == 1) ecore = zero

 !Index table of the atoms, type after type
 ABI_MALLOC(atindx,(dtset%natom))
 ABI_MALLOC(atindx1,(dtset%natom))
 ABI_MALLOC(nattyp,(psps%ntypat))
 indx = 1
 zion = zero
 do itypat=1,psps%ntypat
   nattyp(itypat) = 0
   do iatom=1,dtset%natom
     if (dtset%typat(iatom) == itypat) then
       atindx(iatom) = indx
       atindx1(indx) = iatom
       indx = indx+1
       nattyp(itypat) = nattyp(itypat)+1
       zion = zion+psps%ziontypat(itypat)
     end if
   end do
 end do

 !Structure factor phases, Ewald energy and local pseudopotential
 ABI_MALLOC(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)

 call ewald(e_ewald,gmet,grewtn,gsqcut,dtset%icutcoul,dtset%natom,ngfft,dtset%nkpt,psps%ntypat,dtset%rcut,&
            rmet,rprimd,dtset%typat,ucvol,dtset%vcutgeo,xred,psps%ziontypat)

 ABI_MALLOC(vpsp,(nfft))
 ABI_MALLOC(rhog_dum,(2,nfft))
 ABI_MALLOC(rhor_dum,(nfft,dtset%nspden))
 rhog_dum = zero; rhor_dum = zero; eei = zero; dummy6 = zero
 call nullify_wvl_data(wvl)
 call mklocl(dtset,dyfr_dum,eei,gmet,gprimd,gr_dum,gsqcut,dummy6,dtset%mgfft,mpi_enreg,dtset%natom,nattyp,nfft,ngfft,&
             dtset%nspden,psps%ntypat,option_locl,pawtab,ph1d,psps,dtset%qprtrb,rhog_dum,rhor_dum,rprimd,ucvol,&
             dtset%vprtrb,vpsp,wvl%descr,wvl%den,xred)
 ABI_FREE(rhog_dum)
 ABI_FREE(rhor_dum)

 !Constant terms of the energy (as in ABINIT: usepotzero)
 select case (dtset%usepotzero)
 case (0,1)
   e_corepsp = ecore/ucvol
 case (2)
   e_corepsp = zero
   vpsp(:) = vpsp(:)+ecore/(zion*ucvol)
 end select

 !Hamiltonian: only the local pseudopotential (no Hartree, no xc) and the nonlocal pseudopotential
 call ham_k%init(psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,dtset%typat,xred,nfft,dtset%mgfft,&
                 ngfft,rprimd,dtset%nloalg,ph1d=ph1d,usecprj=0)
 ABI_MALLOC(vtrial,(nfft,dtset%nspden,1))
 vtrial = zero
 vtrial(:,1,1) = vpsp(:)
 ABI_MALLOC(vlocal,(n4,n5,n6,ham_k%nvloc))
 call gspot_transgrid_and_pack(1,psps%usepaw,dtset%paral_kgb,nfft,ngfft,nfft,dtset%nspden,ham_k%nvloc,1,pawfgr,&
                               mpi_enreg,vtrial,vlocal)
 call ham_k%load_spin(1,vlocal=vlocal,with_nonlocal=.true.)

 ABI_MALLOC(h1diag,(orbs%mband,orbs%nkpt))
 ABI_MALLOC(ekin,(orbs%mband,orbs%nkpt))
 h1diag = zero
 ekin = zero
 if (present(hmat)) then
   ABI_MALLOC(hmat,(orbs%mband,orbs%mband,orbs%nkpt))
   hmat = czero
 end if

 do ik=1,orbs%nkpt
   nb = orbs%nband(ik,1)
   npw_k = orbs%npw(ik)

   ABI_MALLOC(kg_k,(3,npw_k))
   kg_k(:,:) = orbs%kg(:,1:npw_k,ik)

   !Kinetic energy, nonlocal form factors and phases at this k-point
   ABI_MALLOC(kinpw,(npw_k))
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,orbs%kpt(:,ik),npw_k,0,0)
   ABI_MALLOC(kpg_k,(npw_k,nkpg))
   ABI_MALLOC(ylm_k,(npw_k,0))
   ABI_MALLOC(ylmgr_k,(npw_k,3,0))
   ABI_MALLOC(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,&
               orbs%kpt(:,ik),psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,&
               psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
   ABI_MALLOC(ph3d,(2,npw_k,ham_k%matblk))
   call ham_k%load_k(kpt_k=orbs%kpt(:,ik),istwf_k=1,npw_k=npw_k,kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,&
                     ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)

   !Apply h to all the bands of this k-point
   ABI_MALLOC(cwavef,(2,npw_k*nb))
   ABI_MALLOC(ghc,(2,npw_k*nb))
   ABI_MALLOC(gsc,(2,npw_k*nb))
   ABI_MALLOC(gvnlxc,(0,0))
   ABI_MALLOC(cwaveprj,(0,0))
   cwavef(:,:) = reshape(orbs%cg(:,1:npw_k,1:nb,ik,1),[2,npw_k*nb])
   call multithreaded_getghc(cpopt,cwavef,cwaveprj,ghc,gsc,ham_k,gvnlxc,one,mpi_enreg,nb,dtset%prtvol,sij_opt,&
                             tim_getghc,type_calc)
   if (present(hmat)) then
     do ib=1,nb
       do jb=1,nb
         ! <ib|h|jb>
         hmat(ib,jb,ik) = cmplx(sum(cwavef(1,(ib-1)*npw_k+1:ib*npw_k)*ghc(1,(jb-1)*npw_k+1:jb*npw_k)+ &
                                    cwavef(2,(ib-1)*npw_k+1:ib*npw_k)*ghc(2,(jb-1)*npw_k+1:jb*npw_k)), &
                                sum(cwavef(1,(ib-1)*npw_k+1:ib*npw_k)*ghc(2,(jb-1)*npw_k+1:jb*npw_k)- &
                                    cwavef(2,(ib-1)*npw_k+1:ib*npw_k)*ghc(1,(jb-1)*npw_k+1:jb*npw_k)),kind=dpc)
       end do
     end do
   end if
   do ib=1,nb
     h1diag(ib,ik) = sum(cwavef(:,(ib-1)*npw_k+1:ib*npw_k)*ghc(:,(ib-1)*npw_k+1:ib*npw_k))
     do ig=1,npw_k
       ekin(ib,ik) = ekin(ib,ik)+kinpw(ig)*sum(cwavef(:,(ib-1)*npw_k+ig)**2)
     end do
   end do

   ABI_FREE(cwavef)
   ABI_FREE(ghc)
   ABI_FREE(gsc)
   ABI_FREE(gvnlxc)
   ABI_FREE(cwaveprj)
   ABI_FREE(ph3d)
   ABI_FREE(ffnl)
   ABI_FREE(ylmgr_k)
   ABI_FREE(ylm_k)
   ABI_FREE(kpg_k)
   ABI_FREE(kinpw)
   ABI_FREE(kg_k)
 end do

 call ham_k%free()
 ABI_FREE(vlocal)
 ABI_FREE(vtrial)
 ABI_FREE(vpsp)
 ABI_FREE(ph1d)
 ABI_FREE(atindx)
 ABI_FREE(atindx1)
 ABI_FREE(nattyp)
 ABI_FREE(xred)

 write(msg,'(3a)') ch10,' RDMFT: one-electron part (kinetic + local and nonlocal pseudopotential) computed with getghc.',ch10
 call wrtout([ab_out,std_out],msg)
 ABI_UNUSED(ierr)

end subroutine rdmft_one_electron
!!***

end module m_rdmft_hamiltonian
!!***
