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
 use m_rdmft_integrals,  only : orbs_to_real
 use m_rdmft_orbitals,   only : rdmft_orbs_t

 implicit none

 private
!!***

 public :: rdmft_gen_fock_k
!!***

contains
!!***

!!****f* m_rdmft_fock/rdmft_gen_fock_k
!! NAME
!!  rdmft_gen_fock_k
!!
!! FUNCTION
!!  Generalized Fock matrix F_pj (p,j = bands of the k-point ik) of the energy above.
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
!!
!! OUTPUT
!!  fock(nb,nb) = F_pj
!!
!! SOURCE
subroutine rdmft_gen_fock_k(orbs,ik,ngfft,gsqcut,icutcoul_j,rcut_j,icutcoul_k,rcut_k,offset,nbk,hk,nu,acoef,ccoef,fock)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: ik,icutcoul_j,icutcoul_k,nbk
 real(dp),intent(in) :: gsqcut,rcut_j,rcut_k
 !arrays
 integer,intent(in) :: ngfft(18),offset(orbs%nkpt)
 complex(dpc),intent(in) :: hk(:,:)
 real(dp),intent(in) :: nu(:),acoef(:,:),ccoef(:,:)
 complex(dpc),intent(out) :: fock(:,:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: izero0=0
 integer :: ib,ig,iq,j,jk,nb,nbj,nfft,p
 real(dp) :: cj,fac
 complex(dpc) :: acc
 type(MPI_type) :: mpi_enreg
 !arrays
 real(dp) :: gmet(3,3),qpt(3)
 real(dp),allocatable :: mall(:,:,:),prod(:,:),rho(:,:),uj(:,:,:),uk(:,:,:),vc(:),vc0(:),wj(:,:,:),yy(:,:)

! *************************************************************************

 nb = orbs%nband(ik,1)
 nfft = product(ngfft(1:3))
 gmet = orbs%gmet
 fock = czero
 if (size(acoef,2) /= nbk .or. size(ccoef,2) /= nbk .or. size(acoef,1) /= nb) then
   ABI_BUG('RDMFT: inconsistent size of the coefficients of the generalized Fock matrix.')
 end if

 call initmpi_seq(mpi_enreg)
 call mpi_enreg%distribfft%init_seq('c',ngfft(2),ngfft(3),'all')
 call mpi_enreg%distribfft%init_seq('f',ngfft(2),ngfft(3),'all')

 ! One-electron part: nu_j h_jp
 do j=1,nb
   do p=1,nb
     fock(p,j) = nu(j)*hk(j,p)
   end do
 end do

 ABI_MALLOC(prod,(2,nfft))
 ABI_MALLOC(rho,(2,nfft))
 ABI_MALLOC(vc,(nfft))
 ABI_MALLOC(vc0,(nfft))
 ABI_MALLOC(uk,(2,nfft,nb))
 ABI_MALLOC(mall,(2,nfft,nb))
 ABI_MALLOC(wj,(2,nfft,nb))
 ABI_MALLOC(yy,(2,nfft))
 call orbs_to_real(orbs,ik,ngfft,mpi_enreg,uk)

 qpt = zero
 call barevcoul(rcut_j,icutcoul_j,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc0)
 vc0(1) = zero

 ! ------------------------------------------------------------------
 ! W_j(G) = sum_q a_jq rho_qq(G) and the Hartree part of F
 ! ------------------------------------------------------------------
 wj = zero
 do jk=1,orbs%nkpt
   nbj = orbs%nband(jk,1)
   if (.not.any(abs(acoef(:,offset(jk)+1:offset(jk)+nbj)) > zero)) cycle
   ABI_MALLOC(uj,(2,nfft,nbj))
   call orbs_to_real(orbs,jk,ngfft,mpi_enreg,uj)
   do ib=1,nbj
     iq = offset(jk)+ib
     if (.not.any(abs(acoef(:,iq)) > zero)) cycle
     prod(1,:) = uj(1,:,ib)**2+uj(2,:,ib)**2
     prod(2,:) = zero
     call fourdp(2,rho,prod,-1,mpi_enreg,nfft,1,ngfft,0)
     do j=1,nb
       if (acoef(j,iq) /= zero) wj(:,:,j) = wj(:,:,j)+acoef(j,iq)*rho(:,:)
     end do
   end do
   ABI_FREE(uj)
 end do
 do j=1,nb
   ! M_pj(G) of all the bands p of k with the band j
   do p=1,nb
     prod(1,:) = uk(1,:,p)*uk(1,:,j)+uk(2,:,p)*uk(2,:,j)
     prod(2,:) = uk(1,:,p)*uk(2,:,j)-uk(2,:,p)*uk(1,:,j)
     call fourdp(2,mall(:,:,p),prod,-1,mpi_enreg,nfft,1,ngfft,0)
   end do
   do p=1,nb
     ! (2/Omega) sum_G vJ(G) conjg(M_pj(G)) W_j(G)
     acc = czero
     do ig=1,nfft
       if (vc0(ig) == zero) cycle
       acc = acc+vc0(ig)*cmplx(mall(1,ig,p)*wj(1,ig,j)+mall(2,ig,p)*wj(2,ig,j), &
                               mall(1,ig,p)*wj(2,ig,j)-mall(2,ig,p)*wj(1,ig,j),kind=dpc)
     end do
     fock(p,j) = fock(p,j)+two*acc/orbs%ucvol
   end do
 end do

 ! ------------------------------------------------------------------
 ! Exchange part: -2 sum_q c_jq (1/Omega) sum_G vK conjg(M_pq(G)) M_jq(G)
 ! ------------------------------------------------------------------
 do jk=1,orbs%nkpt
   nbj = orbs%nband(jk,1)
   if (.not.any(abs(ccoef(:,offset(jk)+1:offset(jk)+nbj)) > zero)) cycle
   ABI_MALLOC(uj,(2,nfft,nbj))
   call orbs_to_real(orbs,jk,ngfft,mpi_enreg,uj)
   qpt(:) = orbs%kpt(:,jk)-orbs%kpt(:,ik)
   call barevcoul(rcut_k,icutcoul_k,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)
   do ib=1,nbj
     iq = offset(jk)+ib
     if (.not.any(abs(ccoef(:,iq)) > zero)) cycle
     ! M_pq(G) for all the bands p of k with the orbital q
     do p=1,nb
       prod(1,:) = uk(1,:,p)*uj(1,:,ib)+uk(2,:,p)*uj(2,:,ib)
       prod(2,:) = uk(1,:,p)*uj(2,:,ib)-uk(2,:,p)*uj(1,:,ib)
       call fourdp(2,mall(:,:,p),prod,-1,mpi_enreg,nfft,1,ngfft,0)
     end do
     do j=1,nb
       cj = ccoef(j,iq)
       if (cj == zero) cycle
       yy(1,:) = vc(:)*mall(1,:,j)
       yy(2,:) = vc(:)*mall(2,:,j)
       fac = -two*cj/orbs%ucvol
       do p=1,nb
         acc = czero
         do ig=1,nfft
           acc = acc+cmplx(mall(1,ig,p)*yy(1,ig)+mall(2,ig,p)*yy(2,ig),mall(1,ig,p)*yy(2,ig)-mall(2,ig,p)*yy(1,ig),kind=dpc)
         end do
         fock(p,j) = fock(p,j)+fac*acc
       end do
     end do
   end do
   ABI_FREE(uj)
 end do

 ABI_FREE(prod)
 ABI_FREE(rho)
 ABI_FREE(vc)
 ABI_FREE(vc0)
 ABI_FREE(uk)
 ABI_FREE(mall)
 ABI_FREE(wj)
 ABI_FREE(yy)
 call destroy_mpi_enreg(mpi_enreg)

end subroutine rdmft_gen_fock_k
!!***

end module m_rdmft_fock
!!***
