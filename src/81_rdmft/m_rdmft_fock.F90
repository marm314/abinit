!!****m* ABINIT/m_rdmft_fock
!! NAME
!!  m_rdmft_fock
!!
!! FUNCTION
!!  Generalized Fock matrix of an RDMFT energy for the rotation of the orbitals of one k-point, and the orbital
!!  gradient. The energy is written as
!!    E = sum_p nu_p h_pp + sum_{pq} [ a_pq J_pq - c_pq K_pq ] ,   a_pq=A_pq/N_k^2, c_pq=C_pq/N_k^2, nu_p=2 n_p/N_k
!!  (J_pq=<pq|pq>, K_pq=<pq|qp>, both with the Coulomb kernels of ABINIT, see m_rdmft_integrals), which covers the
!!  JK-only and the PNOF functionals (including the self-interaction terms, that are diagonal entries of a and c).
!!  Under phi_j -> phi_j + sum_p phi_p kappa_pj one has dE = 2 Re sum_{pj} kappa_pj F_pj with
!!    F_pj = nu_j h_jp + (2/Omega) sum_G v_J(G) conjg(M_pj(G)) W_j(G)
!!                     - 2 sum_q c_jq (1/Omega) sum_G v_K(k_q-k+G) conjg(M_pq(G)) M_jq(G) ,
!!    W_j(G) = sum_q a_jq M_qq(G) ,  M_ab(G) = Fourier components of conjg(u_a) u_b .
!!  This is the plane-wave form of the "RI" contraction: the auxiliary index is the reciprocal vector G, the
!!  Hartree part needs only the weighted density W_j and the pair densities of the orbitals of the same k-point, and the
!!  exchange part the pair densities of all the bands p of k with each orbital q of any k-point. F is block diagonal in k:
!!  only states of the same k-point are rotated.
!!  The gradient for the rotation of the pair (p,q) is dE/dtheta = 2 Re (F_pq - conjg(F_qp)) in the real direction
!!  (kappa_pq=-kappa_qp=theta) and -2 Im (F_pq - conjg(F_qp)) in the imaginary one (kappa_pq=kappa_qp=i theta).
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

module m_rdmft_fock

 use defs_basis
 use defs_abitypes,      only : MPI_type
 use m_abicore
 use m_errors

 use m_barevcoul,        only : barevcoul
 use m_fft,              only : fourdp
 use m_mpinfo,           only : destroy_mpi_enreg, initmpi_seq
 use m_xmpi
 use m_rdmft_integrals,  only : orbs_to_real, orb_to_real_one
 use m_rdmft_orbitals,   only : rdmft_orbs_t

 implicit none

 private
!!***

 public :: rdmft_gen_fock_k
 public :: rdmft_fock_hartree_k
 public :: rdmft_fock_exchange_k
!!***

contains
!!***

!!****f* m_rdmft_fock/rdmft_fock_hartree_k
!! NAME
!!  rdmft_fock_hartree_k
!!
!! FUNCTION
!!  Hartree part of the generalized Fock matrix of the k-point ik:
!!  (2/Omega) sum_G v_J(G) conjg(M_pj(G)) W_j(G), with W_j(G) = sum_q a_jq rho_qq(G).
!!  The work is shared among the processes of comm (the result is the same in all of them).
!!
!! INPUTS
!!  orbs <type(rdmft_orbs_t)> = orbitals
!!  ik = k-point
!!  ngfft(18), gsqcut = FFT mesh and G^2 cut-off of the density
!!  icutcoul_j, rcut_j = Coulomb kernel of J
!!  offset(nkpt) = composite index offsets, nbk = number of orbitals
!!  acoef(nb,nbk) = a_jq (already divided by N_k^2) for the bands j of ik and all the orbitals q
!!  comm = MPI communicator
!!
!! OUTPUT
!!  fock(nb,nb) = Hartree part of F_pj
!!
!! SOURCE
subroutine rdmft_fock_hartree_k(orbs,ik,ngfft,gsqcut,icutcoul_j,rcut_j,offset,nbk,acoef,comm,fock)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: ik,icutcoul_j,nbk,comm
 real(dp),intent(in) :: gsqcut,rcut_j
 !arrays
 integer,intent(in) :: ngfft(18),offset(orbs%nkpt)
 real(dp),intent(in) :: acoef(:,:)
 complex(dpc),intent(out) :: fock(:,:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: izero0=0
 integer :: ib,ierr,iq,j,jk,nb,nbj,nfft,nproc,me,cnt
 type(MPI_type) :: mpi_enreg
 !arrays
 real(dp) :: gmet(3,3),qpt(3)
 real(dp),allocatable :: mall(:,:,:),prod(:,:),rho(:,:),ujl(:,:),uk(:,:,:),vc0(:),wj(:,:,:),yy(:,:,:)

! *************************************************************************

 nb = orbs%nband(ik,1)
 nfft = product(ngfft(1:3))
 gmet = orbs%gmet
 nproc = xmpi_comm_size(comm)
 me = xmpi_comm_rank(comm)
 fock = czero
 if (size(acoef,2) /= nbk .or. size(acoef,1) /= nb) then
   ABI_BUG('RDMFT: inconsistent size of the coefficients of the generalized Fock matrix.')
 end if

 call initmpi_seq(mpi_enreg)
 call mpi_enreg%distribfft%init_seq('c',ngfft(2),ngfft(3),'all')
 call mpi_enreg%distribfft%init_seq('f',ngfft(2),ngfft(3),'all')

 ABI_MALLOC(prod,(2,nfft))
 ABI_MALLOC(ujl,(2,nfft))
 ABI_MALLOC(rho,(2,nfft))
 ABI_MALLOC(vc0,(nfft))
 ABI_MALLOC(uk,(2,nfft,nb))
 ABI_MALLOC(mall,(2,nfft,nb))
 ABI_MALLOC(wj,(2,nfft,nb))
 ABI_MALLOC(yy,(2,nfft,nb))
 call orbs_to_real(orbs,ik,ngfft,mpi_enreg,uk)

 qpt = zero
 call barevcoul(rcut_j,icutcoul_j,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc0)
 vc0(1) = zero

 ! W_j(G): the orbitals q are shared among the processes
 wj = zero
 cnt = 0
 do jk=1,orbs%nkpt
   nbj = orbs%nband(jk,1)
   if (.not.any(abs(acoef(:,offset(jk)+1:offset(jk)+nbj)) > zero)) cycle
   do ib=1,nbj
     iq = offset(jk)+ib
     if (.not.any(abs(acoef(:,iq)) > zero)) cycle
     cnt = cnt+1
     if (mod(cnt,nproc) /= me) cycle
     call orb_to_real_one(orbs,jk,ib,ngfft,mpi_enreg,ujl)
     prod(1,:) = ujl(1,:)**2+ujl(2,:)**2
     prod(2,:) = zero
     call fourdp(2,rho,prod,-1,mpi_enreg,nfft,1,ngfft,0)
     do j=1,nb
       if (acoef(j,iq) /= zero) wj(:,:,j) = wj(:,:,j)+acoef(j,iq)*rho(:,:)
     end do
   end do
 end do
 if (nproc > 1) call xmpi_sum(wj,comm,ierr)

 ! F_pj = (2/Omega) sum_G vc0 conjg(M_pj) W_j: the columns j without couplings (W_j = 0) are skipped and the others
 ! are shared among the processes
 cnt = 0
 do j=1,nb
   if (.not.any(abs(acoef(j,:)) > zero)) cycle
   cnt = cnt+1
   if (mod(cnt,nproc) /= me) cycle
   do ib=1,nb
     prod(1,:) = uk(1,:,ib)*uk(1,:,j)+uk(2,:,ib)*uk(2,:,j)
     prod(2,:) = uk(1,:,ib)*uk(2,:,j)-uk(2,:,ib)*uk(1,:,j)
     call fourdp(2,mall(:,:,ib),prod,-1,mpi_enreg,nfft,1,ngfft,0)
   end do
   do ib=1,nfft
     yy(1,ib,1) = vc0(ib)*wj(1,ib,j)
     yy(2,ib,1) = vc0(ib)*wj(2,ib,j)
   end do
   call zgemm('C','N',nb,1,nfft,dcmplx(two/orbs%ucvol,zero),mall,nfft,yy(1,1,1),nfft,czero,fock(1,j),nb)
 end do
 if (nproc > 1) call xmpi_sum(fock,comm,ierr)

 ABI_FREE(prod)
 ABI_FREE(ujl)
 ABI_FREE(rho)
 ABI_FREE(vc0)
 ABI_FREE(uk)
 ABI_FREE(mall)
 ABI_FREE(wj)
 ABI_FREE(yy)
 call destroy_mpi_enreg(mpi_enreg)

end subroutine rdmft_fock_hartree_k
!!***

!!****f* m_rdmft_fock/rdmft_fock_exchange_k
!! NAME
!!  rdmft_fock_exchange_k
!!
!! FUNCTION
!!  Exchange part of the generalized Fock matrix of the k-point ik:
!!  -2 sum_q c_jq (1/Omega) sum_G v_K(k_q-k+G) conjg(M_pq(G)) M_jq(G).
!!  The orbitals q are shared among the processes of comm. (The integrals routine rdmft_ints_t%compute can build the
!!  same matrix in the same loop as the K integrals, from the same pair densities.)
!!
!! INPUTS
!!  orbs, ik, ngfft, gsqcut, offset, nbk, comm = as in rdmft_fock_hartree_k
!!  icutcoul_k, rcut_k = Coulomb kernel of K
!!  ccoef(nb,nbk) = c_jq (already divided by N_k^2)
!!
!! OUTPUT
!!  fock(nb,nb) = exchange part of F_pj
!!
!! SOURCE
subroutine rdmft_fock_exchange_k(orbs,ik,ngfft,gsqcut,icutcoul_k,rcut_k,offset,nbk,ccoef,comm,fock)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: ik,icutcoul_k,nbk,comm
 real(dp),intent(in) :: gsqcut,rcut_k
 !arrays
 integer,intent(in) :: ngfft(18),offset(orbs%nkpt)
 real(dp),intent(in) :: ccoef(:,:)
 complex(dpc),intent(out) :: fock(:,:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: izero0=0
 integer :: ib,ierr,iq,j,jk,nb,nbj,nfft,nproc,me,cnt
 type(MPI_type) :: mpi_enreg
 !arrays
 real(dp) :: gmet(3,3),qpt(3)
 real(dp),allocatable :: mall(:,:,:),prod(:,:),ujl(:,:),uk(:,:,:),vc(:),yy(:,:,:)

! *************************************************************************

 nb = orbs%nband(ik,1)
 nfft = product(ngfft(1:3))
 gmet = orbs%gmet
 nproc = xmpi_comm_size(comm)
 me = xmpi_comm_rank(comm)
 fock = czero
 if (size(ccoef,2) /= nbk .or. size(ccoef,1) /= nb) then
   ABI_BUG('RDMFT: inconsistent size of the coefficients of the generalized Fock matrix.')
 end if

 call initmpi_seq(mpi_enreg)
 call mpi_enreg%distribfft%init_seq('c',ngfft(2),ngfft(3),'all')
 call mpi_enreg%distribfft%init_seq('f',ngfft(2),ngfft(3),'all')

 ABI_MALLOC(prod,(2,nfft))
 ABI_MALLOC(ujl,(2,nfft))
 ABI_MALLOC(vc,(nfft))
 ABI_MALLOC(uk,(2,nfft,nb))
 ABI_MALLOC(mall,(2,nfft,nb))
 ABI_MALLOC(yy,(2,nfft,nb))
 call orbs_to_real(orbs,ik,ngfft,mpi_enreg,uk)

 cnt = 0
 do jk=1,orbs%nkpt
   nbj = orbs%nband(jk,1)
   if (.not.any(abs(ccoef(:,offset(jk)+1:offset(jk)+nbj)) > zero)) cycle
   qpt(:) = orbs%kpt(:,jk)-orbs%kpt(:,ik)
   call barevcoul(rcut_k,icutcoul_k,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)
   do ib=1,nbj
     iq = offset(jk)+ib
     if (.not.any(abs(ccoef(:,iq)) > zero)) cycle
     cnt = cnt+1
     if (mod(cnt,nproc) /= me) cycle
     call orb_to_real_one(orbs,jk,ib,ngfft,mpi_enreg,ujl)
     ! M_pq(G) for all the bands p of k with the orbital q
     do j=1,nb
       prod(1,:) = uk(1,:,j)*ujl(1,:)+uk(2,:,j)*ujl(2,:)
       prod(2,:) = uk(1,:,j)*ujl(2,:)-uk(2,:,j)*ujl(1,:)
       call fourdp(2,mall(:,:,j),prod,-1,mpi_enreg,nfft,1,ngfft,0)
     end do
     ! Y_j = -2 c_jq / Omega v M_jq ;  F_pj += sum_G conjg(M_p) Y_j
     do j=1,nb
       yy(1,:,j) = (-two*ccoef(j,iq)/orbs%ucvol)*vc(:)*mall(1,:,j)
       yy(2,:,j) = (-two*ccoef(j,iq)/orbs%ucvol)*vc(:)*mall(2,:,j)
     end do
     call zgemm('C','N',nb,nb,nfft,cone,mall,nfft,yy,nfft,cone,fock,nb)
   end do
 end do
 if (nproc > 1) call xmpi_sum(fock,comm,ierr)

 ABI_FREE(prod)
 ABI_FREE(ujl)
 ABI_FREE(vc)
 ABI_FREE(uk)
 ABI_FREE(mall)
 ABI_FREE(yy)
 call destroy_mpi_enreg(mpi_enreg)

end subroutine rdmft_fock_exchange_k
!!***

!!****f* m_rdmft_fock/rdmft_gen_fock_k
!! NAME
!!  rdmft_gen_fock_k
!!
!! FUNCTION
!!  Generalized Fock matrix F_pj (p,j = bands of the k-point ik) of the energy above: one-electron, Hartree and
!!  exchange parts.
!!
!! INPUTS
!!  orbs <type(rdmft_orbs_t)> = orbitals
!!  ik = k-point
!!  ngfft(18), gsqcut = FFT mesh and G^2 cut-off of the density
!!  icutcoul_j, rcut_j, icutcoul_k, rcut_k = Coulomb kernels of J and K
!!  offset(nkpt) = composite index offsets (index of the band 1 of each k-point minus one), nbk = number of orbitals
!!  hk(nb,nb) = <i|h|j> of the bands of ik
!!  nu(nb) = 2 n_j / N_k
!!  acoef(nb,nbk), ccoef(nb,nbk) = a_jq and c_jq (already divided by N_k^2) for the bands j of ik and all the orbitals q
!!  comm = MPI communicator
!!
!! OUTPUT
!!  fock(nb,nb) = F_pj
!!
!! SOURCE
subroutine rdmft_gen_fock_k(orbs,ik,ngfft,gsqcut,icutcoul_j,rcut_j,icutcoul_k,rcut_k,offset,nbk,hk,nu,acoef,ccoef,comm,fock)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: ik,icutcoul_j,icutcoul_k,nbk,comm
 real(dp),intent(in) :: gsqcut,rcut_j,rcut_k
 !arrays
 integer,intent(in) :: ngfft(18),offset(orbs%nkpt)
 complex(dpc),intent(in) :: hk(:,:)
 real(dp),intent(in) :: nu(:),acoef(:,:),ccoef(:,:)
 complex(dpc),intent(out) :: fock(:,:)

 !Local variables-------------------------------
 !scalars
 integer :: j,nb,p
 !arrays
 complex(dpc),allocatable :: fx(:,:)

! *************************************************************************

 nb = orbs%nband(ik,1)
 ABI_MALLOC(fx,(nb,nb))
 call rdmft_fock_hartree_k(orbs,ik,ngfft,gsqcut,icutcoul_j,rcut_j,offset,nbk,acoef,comm,fock)
 call rdmft_fock_exchange_k(orbs,ik,ngfft,gsqcut,icutcoul_k,rcut_k,offset,nbk,ccoef,comm,fx)
 do j=1,nb
   do p=1,nb
     fock(p,j) = fock(p,j)+fx(p,j)+nu(j)*hk(j,p)
   end do
 end do
 ABI_FREE(fx)

end subroutine rdmft_gen_fock_k
!!***

end module m_rdmft_fock
!!***
