!!****m* ABINIT/m_rdmft_integrals
!! NAME
!!  m_rdmft_integrals
!!
!! FUNCTION
!!  Coulomb (J) and exchange (K) integrals between the Bloch orbitals of the full Brillouin zone.
!!  With p=(n,k) and q=(m,k'), cell-normalized integrals (N_k times the integrals of the whole crystal):
!!    J_pq = <pq|pq> = (1/ucvol) sum_{G/=0} 4pi/|G|^2 conjg(rho_pp(G)) rho_qq(G),
!!    K_pq = <pq|qp> = (1/ucvol) sum_G v(k'-k+G) |M_pq(G)|^2,
!!  where M_pq(G) is the Fourier transform of conjg(u_p(r)) u_q(r) (u: periodic part of the orbital)
!!  and v is the Coulomb kernel of ABINIT (barevcoul), which treats the q+G=0 singularity: icutcoul, rcut for J (as the
!!  Hartree term of the ground state) and fock_icutcoul, fock_rcut for K (as exact exchange). The G=0 term of J is dropped.
!!  Time-reversal symmetry (no spin-orbit, collinear spin) gives L_pq=K_pq (see RDMFT/L_vs_K.tex).
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

module m_rdmft_integrals

 use defs_basis
 use defs_abitypes,      only : MPI_type
 use m_abicore
 use m_errors
 use m_xmpi

 use m_barevcoul,        only : barevcoul
 use m_fft,              only : fourdp
 use m_mpinfo,           only : destroy_mpi_enreg, initmpi_seq
 use m_rdmft_orbitals,   only : rdmft_orbs_t
 use m_specialmsg,       only : wrtout

 implicit none

 private
!!***

!!****t* m_rdmft_integrals/rdmft_ints_t
!! NAME
!! rdmft_ints_t
!!
!! FUNCTION
!!  J and K integrals between all the orbitals (n,k) of the full BZ.
!!  The composite orbital index is ip = offset(ik) + n.
!!
!! SOURCE

 type,public :: rdmft_ints_t

   integer :: nbk = 0
   ! Total number of orbitals (sum over k-points of the number of bands).

   integer,allocatable :: offset(:)
   ! offset(nkpt)
   ! Index of the orbital (1,k) minus one.

   real(dp),allocatable :: jint(:,:)
   ! jint(nbk,nbk)
   ! Coulomb integrals J_pq (Hartree).

   real(dp),allocatable :: kint(:,:)
   ! kint(nbk,nbk)
   ! Exchange integrals K_pq (Hartree), equal to the time-inversion integrals L_pq.

 contains

   procedure :: free => rdmft_ints_free
   ! Free memory.

   procedure :: compute => rdmft_ints_compute
   ! Compute the J and K integrals.

 end type rdmft_ints_t
!!***

contains
!!***

!!****f* m_rdmft_integrals/rdmft_ints_free
!! NAME
!!  rdmft_ints_free
!!
!! FUNCTION
!!  Free the memory of the integrals.
!!
!! SOURCE
subroutine rdmft_ints_free(ints)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_ints_t),intent(inout) :: ints

! *************************************************************************

 ABI_SFREE(ints%offset)
 ABI_SFREE(ints%jint)
 ABI_SFREE(ints%kint)
 ints%nbk = 0

end subroutine rdmft_ints_free
!!***

!!****f* m_rdmft_integrals/orbs_to_real
!! NAME
!!  orbs_to_real
!!
!! FUNCTION
!!  Periodic part u(r) = sum_G c(G) exp(iG.r) of all the bands of one k-point on the FFT mesh
!!  (plane waves placed in the FFT box and transformed with fourdp).
!!
!! INPUTS
!!  orbs <type(rdmft_orbs_t)> = ground-state orbitals
!!  ik = index of the k-point
!!  ngfft(18) = FFT mesh
!!  mpi_enreg <type(MPI_type)> = sequential MPI information needed by fourdp
!!
!! OUTPUT
!!  ur(2,nfft,nband) = real and imaginary parts of u(r), with sum_r |u(r)|^2/nfft = 1
!!
!! SOURCE
subroutine orbs_to_real(orbs,ik,ngfft,mpi_enreg,ur)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: ik
 type(MPI_type),intent(in) :: mpi_enreg
 !arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: ur(:,:,:)

 !Local variables-------------------------------
 !scalars
 integer :: ib,ig,i1,i2,i3,nb,nfft
 !arrays
 real(dp),allocatable :: ug(:,:,:)

! *************************************************************************

 nfft = product(ngfft(1:3))
 nb = orbs%nband(ik,1)

 ABI_MALLOC(ug,(2,nfft,nb))
 ug = zero
 do ig=1,orbs%npw(ik)
   i1 = modulo(orbs%kg(1,ig,ik),ngfft(1))
   i2 = modulo(orbs%kg(2,ig,ik),ngfft(2))
   i3 = modulo(orbs%kg(3,ig,ik),ngfft(3))
   do ib=1,nb
     ug(:,1+i1+ngfft(1)*(i2+ngfft(2)*i3),ib) = orbs%cg(:,ig,ib,ik,1)
   end do
 end do

 call fourdp(2,ug,ur,+1,mpi_enreg,nfft,nb,ngfft,0)
 ABI_FREE(ug)

end subroutine orbs_to_real
!!***

!!****f* m_rdmft_integrals/rdmft_ints_compute
!! NAME
!!  rdmft_ints_compute
!!
!! FUNCTION
!!  Compute the J and K integrals for all pairs of orbitals of the full BZ.
!!
!! INPUTS
!!  orbs <type(rdmft_orbs_t)> = ground-state orbitals
!!  ngfft(18) = FFT mesh (must contain the products of two orbitals, boxcut=2)
!!  gsqcut = cutoff on G^2 for the FFT mesh (as in the density calculation)
!!  icutcoul_j, rcut_j = Coulomb cut-off option and radius for J (icutcoul, rcut: as the Hartree term of the ground state)
!!  icutcoul_k, rcut_k = Coulomb cut-off option and radius for K (fock_icutcoul, fock_rcut: as exact exchange)
!!
!! OUTPUT
!!  ints%jint, ints%kint
!!
!! SOURCE
subroutine rdmft_ints_compute(ints,orbs,ngfft,gsqcut,icutcoul_j,rcut_j,icutcoul_k,rcut_k)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_ints_t),intent(inout) :: ints
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: icutcoul_j,icutcoul_k
 real(dp),intent(in) :: gsqcut,rcut_j,rcut_k
 !arrays
 integer,intent(in) :: ngfft(18)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: izero0=0
 integer :: ib,ik,iq,ip,jb,jk,nb,nbj,nfft
 real(dp) :: kval
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg
 !arrays
 real(dp) :: qpt(3),gmet(3,3)
 real(dp),allocatable :: prod(:,:),pgrid(:,:),rhog(:,:,:),vc(:),vq(:,:),ui(:,:,:),uj(:,:,:)

! *************************************************************************

 call ints%free()
 nfft = product(ngfft(1:3))
 gmet = orbs%gmet

 !Sequential MPI_type, needed by the FFT routines
 call initmpi_seq(mpi_enreg)
 call mpi_enreg%distribfft%init_seq('c',ngfft(2),ngfft(3),'all')
 call mpi_enreg%distribfft%init_seq('f',ngfft(2),ngfft(3),'all')

 !Composite orbital index
 ABI_MALLOC(ints%offset,(orbs%nkpt))
 ints%nbk = 0
 do ik=1,orbs%nkpt
   ints%offset(ik) = ints%nbk
   ints%nbk = ints%nbk+orbs%nband(ik,1)
 end do
 ABI_MALLOC(ints%jint,(ints%nbk,ints%nbk))
 ABI_MALLOC(ints%kint,(ints%nbk,ints%nbk))
 ints%jint = zero
 ints%kint = zero

 ABI_MALLOC(prod,(2,nfft))
 ABI_MALLOC(pgrid,(2,nfft))
 ABI_MALLOC(vc,(nfft))

 ! ------------------------------------------------------------------
 ! Coulomb integrals: rho_pp(G) of every orbital, v(G/=0) at q=0
 ! ------------------------------------------------------------------
 ABI_MALLOC(rhog,(2,nfft,ints%nbk))
 ABI_MALLOC(vq,(2,nfft))
 qpt = zero
 call barevcoul(rcut_j,icutcoul_j,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)
 vc(1) = zero   ! G=0 is the first element of the FFT mesh: neutralizing background

 do ik=1,orbs%nkpt
   nb = orbs%nband(ik,1)
   ABI_MALLOC(ui,(2,nfft,nb))
   call orbs_to_real(orbs,ik,ngfft,mpi_enreg,ui)
   do ib=1,nb
     prod(1,:) = ui(1,:,ib)**2+ui(2,:,ib)**2
     prod(2,:) = zero
     call fourdp(2,rhog(:,:,ints%offset(ik)+ib),prod,-1,mpi_enreg,nfft,1,ngfft,0)
   end do
   ABI_FREE(ui)
 end do

 do iq=1,ints%nbk
   vq(1,:) = vc(:)*rhog(1,:,iq)
   vq(2,:) = vc(:)*rhog(2,:,iq)
   do ip=1,iq
     ! J_pq = (1/ucvol) sum_G v(G) conjg(rho_p(G)) rho_q(G) (real part)
     kval = sum(rhog(1,:,ip)*vq(1,:)+rhog(2,:,ip)*vq(2,:))/orbs%ucvol
     ints%jint(ip,iq) = kval
     ints%jint(iq,ip) = kval
   end do
 end do
 ABI_FREE(rhog)
 ABI_FREE(vq)

 ! ------------------------------------------------------------------
 ! Exchange integrals: pair densities of all the pairs (n,k), (m,k'), k <= k'
 ! ------------------------------------------------------------------
 do ik=1,orbs%nkpt
   nb = orbs%nband(ik,1)
   ABI_MALLOC(ui,(2,nfft,nb))
   call orbs_to_real(orbs,ik,ngfft,mpi_enreg,ui)
   do jk=ik,orbs%nkpt
     nbj = orbs%nband(jk,1)
     ABI_MALLOC(uj,(2,nfft,nbj))
     call orbs_to_real(orbs,jk,ngfft,mpi_enreg,uj)

     !Coulomb kernel for q = k' - k
     qpt(:) = orbs%kpt(:,jk)-orbs%kpt(:,ik)
     call barevcoul(rcut_k,icutcoul_k,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)

     do ib=1,nb
       do jb=1,nbj
         if (ik == jk .and. jb > ib) cycle
         ! conjg(u_i) u_j
         prod(1,:) = ui(1,:,ib)*uj(1,:,jb)+ui(2,:,ib)*uj(2,:,jb)
         prod(2,:) = ui(1,:,ib)*uj(2,:,jb)-ui(2,:,ib)*uj(1,:,jb)
         call fourdp(2,pgrid,prod,-1,mpi_enreg,nfft,1,ngfft,0)
         kval = sum(vc(:)*(pgrid(1,:)**2+pgrid(2,:)**2))/orbs%ucvol
         ip = ints%offset(ik)+ib
         iq = ints%offset(jk)+jb
         ints%kint(ip,iq) = kval
         ints%kint(iq,ip) = kval
       end do
     end do
     ABI_FREE(uj)
   end do
   ABI_FREE(ui)
 end do

 ABI_FREE(prod)
 ABI_FREE(pgrid)
 ABI_FREE(vc)
 call destroy_mpi_enreg(mpi_enreg)

 write(msg,'(3a,i0,a)') ch10,' RDMFT: J and K integrals computed for all pairs of orbitals of the full BZ,',&
   ' number of orbitals = ',ints%nbk,ch10
 call wrtout([ab_out,std_out],msg)

end subroutine rdmft_ints_compute
!!***

end module m_rdmft_integrals
!!***
