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

 public :: orbs_to_real

 type,public :: rdmft_ints_t

   integer :: nbk = 0
   ! Total number of orbitals (sum over k-points of the number of bands).

   integer :: nfree = 0
   ! Number of free orbitals (the active ones, whose occupation is optimized).

   integer :: nbg = 0
   ! Number of background orbitals (frozen occupied orbitals, occupation 1).

   integer,allocatable :: offset(:)
   ! offset(nkpt)
   ! Index of the orbital (1,k) minus one.

   integer,allocatable :: free_of(:)
   ! free_of(nbk)
   ! Index of the orbital among the free orbitals (0 if it is not free).

   integer,allocatable :: comp_of(:)
   ! comp_of(nfree)
   ! Composite index (offset+band) of each free orbital.

   real(dp),allocatable :: jint(:,:)
   ! jint(nfree,nfree)
   ! Coulomb integrals J_pq (Hartree) between the free orbitals.

   real(dp),allocatable :: kint(:,:)
   ! kint(nfree,nfree)
   ! Exchange integrals K_pq (Hartree) between the free orbitals, equal to the time-inversion integrals L_pq.

   real(dp),allocatable :: bj(:)
   ! bj(nfree)
   ! Sum over the background orbitals q of J_pq, for each free orbital p.

   real(dp),allocatable :: bk(:)
   ! bk(nfree)
   ! Sum over the background orbitals q of K_pq, for each free orbital p.

   real(dp) :: sjo = zero
   ! Sum of J_pq over the pairs p/=q of background orbitals (ordered pairs).

   real(dp) :: sko = zero
   ! Sum of K_pq over the pairs p/=q of background orbitals (ordered pairs); needs core_ints.

   real(dp) :: sjd = zero
   ! Sum of J_pp over the background orbitals.

   real(dp) :: skd = zero
   ! Sum of K_pp over the background orbitals; needs core_ints.

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
 ABI_SFREE(ints%free_of)
 ABI_SFREE(ints%comp_of)
 ABI_SFREE(ints%jint)
 ABI_SFREE(ints%kint)
 ABI_SFREE(ints%bj)
 ABI_SFREE(ints%bk)
 ints%nbk = 0
 ints%nfree = 0
 ints%nbg = 0
 ints%sjo = zero
 ints%sko = zero
 ints%sjd = zero
 ints%skd = zero

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
!!  Compute the J and K integrals restricted to the active space: J_pq and K_pq between the free orbitals, and
!!  for each free orbital p the sums B_J(p)=sum_q J_pq and B_K(p)=sum_q K_pq over the background (frozen
!!  occupied) orbitals q, which have occupation 1 and therefore only enter through these sums. The background-
!!  background pairs only give constants (sjo, sko, sjd, skd); the exchange ones need all the pairs of background
!!  orbitals and are computed only if core_ints is true.
!!  The other orbitals (frozen virtual) have occupation 0 and do not contribute.
!!
!! INPUTS
!!  orbs <type(rdmft_orbs_t)> = ground-state orbitals
!!  ngfft(18) = FFT mesh (must contain the products of two orbitals, boxcut=2)
!!  gsqcut = cutoff on G^2 for the FFT mesh (as in the density calculation)
!!  icutcoul_j, rcut_j = Coulomb cut-off option and radius for J (icutcoul, rcut: as the Hartree term of the ground state)
!!  icutcoul_k, rcut_k = Coulomb cut-off option and radius for K (fock_icutcoul, fock_rcut: as exact exchange)
!!  free(nbk), bg(nbk) = free and background orbitals (composite index), disjoint
!!  core_ints = if true compute the exchange constants of the background-background pairs
!!  quiet = (optional) if true, do not print the summary line
!!
!! OUTPUT
!!  ints
!!
!! SOURCE
subroutine rdmft_ints_compute(ints,orbs,ngfft,gsqcut,icutcoul_j,rcut_j,icutcoul_k,rcut_k,free,bg,core_ints,quiet)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_ints_t),intent(inout) :: ints
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: icutcoul_j,icutcoul_k
 real(dp),intent(in) :: gsqcut,rcut_j,rcut_k
 logical,intent(in) :: core_ints
 logical,optional,intent(in) :: quiet
 !arrays
 integer,intent(in) :: ngfft(18)
 logical,intent(in) :: free(:),bg(:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: izero0=0
 integer :: ib,ik,iq,ip,jb,jk,jp,nb,nbj,nfft
 real(dp) :: kval
 logical :: fp,fq,bp,bq
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg
 !arrays
 real(dp) :: qpt(3),gmet(3,3)
 real(dp),allocatable :: prod(:,:),pgrid(:,:),rhog(:,:,:),rbg(:,:),rtmp(:,:),vc(:),vq(:,:),ui(:,:,:),uj(:,:,:)

! *************************************************************************

 call ints%free()
 nfft = product(ngfft(1:3))
 gmet = orbs%gmet

 !Sequential MPI_type, needed by the FFT routines
 call initmpi_seq(mpi_enreg)
 call mpi_enreg%distribfft%init_seq('c',ngfft(2),ngfft(3),'all')
 call mpi_enreg%distribfft%init_seq('f',ngfft(2),ngfft(3),'all')

 !Composite orbital index and free orbitals
 ABI_MALLOC(ints%offset,(orbs%nkpt))
 ints%nbk = 0
 do ik=1,orbs%nkpt
   ints%offset(ik) = ints%nbk
   ints%nbk = ints%nbk+orbs%nband(ik,1)
 end do
 if (size(free) /= ints%nbk .or. size(bg) /= ints%nbk) then
   ABI_BUG('RDMFT: inconsistent size of the free/background masks.')
 end if
 if (any(free .and. bg)) then
   ABI_BUG('RDMFT: an orbital cannot be both free and background.')
 end if
 ints%nfree = count(free)
 ints%nbg = count(bg)
 ABI_MALLOC(ints%free_of,(ints%nbk))
 ABI_MALLOC(ints%comp_of,(ints%nfree))
 ints%free_of = 0
 ip = 0
 do jp=1,ints%nbk
   if (free(jp)) then
     ip = ip+1
     ints%free_of(jp) = ip
     ints%comp_of(ip) = jp
   end if
 end do
 ABI_MALLOC(ints%jint,(ints%nfree,ints%nfree))
 ABI_MALLOC(ints%kint,(ints%nfree,ints%nfree))
 ABI_MALLOC(ints%bj,(ints%nfree))
 ABI_MALLOC(ints%bk,(ints%nfree))
 ints%jint = zero
 ints%kint = zero
 ints%bj = zero
 ints%bk = zero

 ABI_MALLOC(prod,(2,nfft))
 ABI_MALLOC(pgrid,(2,nfft))
 ABI_MALLOC(vc,(nfft))

 ! ------------------------------------------------------------------
 ! Coulomb integrals: rho_pp(G) of the free orbitals, sum of rho_qq(G) over the background orbitals
 ! ------------------------------------------------------------------
 ABI_MALLOC(rhog,(2,nfft,max(ints%nfree,1)))
 ABI_MALLOC(rbg,(2,nfft))
 ABI_MALLOC(rtmp,(2,nfft))
 ABI_MALLOC(vq,(2,nfft))
 rbg = zero
 qpt = zero
 call barevcoul(rcut_j,icutcoul_j,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)
 vc(1) = zero   ! G=0 is the first element of the FFT mesh: neutralizing background

 do ik=1,orbs%nkpt
   nb = orbs%nband(ik,1)
   if (.not.any(free(ints%offset(ik)+1:ints%offset(ik)+nb) .or. bg(ints%offset(ik)+1:ints%offset(ik)+nb))) cycle
   ABI_MALLOC(ui,(2,nfft,nb))
   call orbs_to_real(orbs,ik,ngfft,mpi_enreg,ui)
   do ib=1,nb
     ip = ints%offset(ik)+ib
     if (.not.(free(ip) .or. bg(ip))) cycle
     prod(1,:) = ui(1,:,ib)**2+ui(2,:,ib)**2
     prod(2,:) = zero
     if (free(ip)) then
       call fourdp(2,rhog(:,:,ints%free_of(ip)),prod,-1,mpi_enreg,nfft,1,ngfft,0)
     else
       call fourdp(2,rtmp,prod,-1,mpi_enreg,nfft,1,ngfft,0)
       rbg = rbg+rtmp
       ints%sjd = ints%sjd+sum(vc(:)*(rtmp(1,:)**2+rtmp(2,:)**2))/orbs%ucvol
     end if
   end do
   ABI_FREE(ui)
 end do

 do iq=1,ints%nfree
   vq(1,:) = vc(:)*rhog(1,:,iq)
   vq(2,:) = vc(:)*rhog(2,:,iq)
   do ip=1,iq
     kval = sum(rhog(1,:,ip)*vq(1,:)+rhog(2,:,ip)*vq(2,:))/orbs%ucvol
     ints%jint(ip,iq) = kval
     ints%jint(iq,ip) = kval
   end do
   ints%bj(iq) = sum(vq(1,:)*rbg(1,:)+vq(2,:)*rbg(2,:))/orbs%ucvol
 end do
 ! Background-background pairs: |sum_q rho_q|^2 contains the diagonal terms
 ints%sjo = sum(vc(:)*(rbg(1,:)**2+rbg(2,:)**2))/orbs%ucvol-ints%sjd
 ABI_FREE(rhog)
 ABI_FREE(rbg)
 ABI_FREE(rtmp)
 ABI_FREE(vq)

 ! ------------------------------------------------------------------
 ! Exchange integrals: pair densities of the needed pairs (n,k), (m,k'), k <= k'
 ! ------------------------------------------------------------------
 do ik=1,orbs%nkpt
   nb = orbs%nband(ik,1)
   if (.not.any(free(ints%offset(ik)+1:ints%offset(ik)+nb) .or. bg(ints%offset(ik)+1:ints%offset(ik)+nb))) cycle
   ABI_MALLOC(ui,(2,nfft,nb))
   call orbs_to_real(orbs,ik,ngfft,mpi_enreg,ui)
   do jk=ik,orbs%nkpt
     nbj = orbs%nband(jk,1)
     if (.not.any(free(ints%offset(jk)+1:ints%offset(jk)+nbj) .or. bg(ints%offset(jk)+1:ints%offset(jk)+nbj))) cycle
     ABI_MALLOC(uj,(2,nfft,nbj))
     call orbs_to_real(orbs,jk,ngfft,mpi_enreg,uj)

     !Coulomb kernel for q = k' - k
     qpt(:) = orbs%kpt(:,jk)-orbs%kpt(:,ik)
     call barevcoul(rcut_k,icutcoul_k,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)

     do ib=1,nb
       ip = ints%offset(ik)+ib
       fp = free(ip); bp = bg(ip)
       if (.not.(fp .or. bp)) cycle
       do jb=1,nbj
         if (ik == jk .and. jb > ib) cycle
         iq = ints%offset(jk)+jb
         fq = free(iq); bq = bg(iq)
         if (.not.(fq .or. bq)) cycle
         if (bp .and. bq .and. .not.core_ints) cycle
         ! conjg(u_i) u_j
         prod(1,:) = ui(1,:,ib)*uj(1,:,jb)+ui(2,:,ib)*uj(2,:,jb)
         prod(2,:) = ui(1,:,ib)*uj(2,:,jb)-ui(2,:,ib)*uj(1,:,jb)
         call fourdp(2,pgrid,prod,-1,mpi_enreg,nfft,1,ngfft,0)
         kval = sum(vc(:)*(pgrid(1,:)**2+pgrid(2,:)**2))/orbs%ucvol
         if (fp .and. fq) then
           ints%kint(ints%free_of(ip),ints%free_of(iq)) = kval
           ints%kint(ints%free_of(iq),ints%free_of(ip)) = kval
         else if (fp .and. bq) then
           ints%bk(ints%free_of(ip)) = ints%bk(ints%free_of(ip))+kval
         else if (bp .and. fq) then
           ints%bk(ints%free_of(iq)) = ints%bk(ints%free_of(iq))+kval
         else
           ! background-background pair
           if (ip == iq) then
             ints%skd = ints%skd+kval
           else
             ints%sko = ints%sko+two*kval
           end if
         end if
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

 if (present(quiet)) then
   if (quiet) return
 end if
 write(msg,'(3a,i0,a,i0,a,i0,2a)') ch10,' RDMFT: J and K integrals computed in the active space,',&
   ' number of free orbitals = ',ints%nfree,', background orbitals = ',ints%nbg,', total orbitals = ',ints%nbk,ch10
 call wrtout([ab_out,std_out],msg)

end subroutine rdmft_ints_compute
!!***

end module m_rdmft_integrals
!!***
