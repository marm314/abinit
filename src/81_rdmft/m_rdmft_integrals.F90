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
 use m_rdmft_orbitals,   only : rdmft_orbs_t, rdmft_symtab_t
 use m_io_tools,         only : get_unit
 use m_specialmsg,       only : wrtout
 use m_time,              only : abi_wtime

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
 public :: rdmft_tint
 public :: rdmft_tprof
 public :: rdmft_mem_on
 public :: rdmft_mem_report

 !> Debug: print the memory of the process inside rdmft_ints_compute (set by the driver)
 logical,save :: rdmft_mem_on = .false.
 public :: orb_to_real_one

 !> Accumulated wall times (s) of the phases of rdmft_ints_compute on this process: 1 densities, 2 sum of the densities,
 !> 3 J products and sum, 4 K pair sweep, 5 wait at the end of the K sweep (load imbalance), 6 sum of the K results
 real(dp),save :: rdmft_tint(6) = 0.0_dp
 !> Profile of the K sweep: 1 orbitals to real space (k), 2 orbital to real space (k'), 3 product of the orbitals,
 !> 4 FFT of the pair density, 5 K value and storage, 6 Fock part (inverse FFT, accumulation, ZGEMM),
 !> 7 number of pair densities (FFTs)
 real(dp),save :: rdmft_tprof(7) = 0.0_dp

 !> Request of the exchange part of the generalized Fock matrices, built in the loop of the K integrals
 type,public :: rdmft_fockreq_t
   logical,allocatable :: want(:)
   ! want(nkpt): the exchange part of F is built for this k-point
   real(dp),allocatable :: ccoef(:,:)
   ! ccoef(nrow,ncol): c_jq (divided by N_k^2) between the relevant bands j (all the k-points) and the relevant orbitals q
   integer,allocatable :: rowmap(:,:)
   ! rowmap(nbmax,nkpt): row of ccoef of each band of each k-point (0 = no coupling); the rows of all the k-points are stacked
   integer,allocatable :: colmap(:)
   ! colmap(nbk): column of ccoef of each orbital (0 = no coupling)
   complex(dpc),allocatable :: fx(:,:,:)
   ! fx(nbmax,nbmax,nkpt): exchange part of F_pj (output)
 end type rdmft_fockreq_t

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

!!****f* m_rdmft_integrals/orb_to_real_one
!! NAME
!!  orb_to_real_one
!!
!! FUNCTION
!!  Periodic part of the band ib of the k-point ik on the real-space FFT mesh (one band).
!!
!! SOURCE
subroutine orb_to_real_one(orbs,ik,ib,ngfft,mpi_enreg,ur)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: ik,ib
 type(MPI_type),intent(in) :: mpi_enreg
 !arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: ur(:,:)

 !Local variables-------------------------------
 !scalars
 integer :: ig,i1,i2,i3,nfft
 !arrays
 real(dp),allocatable :: ug(:,:)

! *************************************************************************

 nfft = product(ngfft(1:3))
 ABI_MALLOC(ug,(2,nfft))
 ug = zero
 do ig=1,orbs%npw(ik)
   i1 = modulo(orbs%kg(1,ig,ik),ngfft(1))
   i2 = modulo(orbs%kg(2,ig,ik),ngfft(2))
   i3 = modulo(orbs%kg(3,ig,ik),ngfft(3))
   ug(:,1+i1+ngfft(1)*(i2+ngfft(2)*i3)) = orbs%cg(:,ig,ib,ik,1)
 end do
 call fourdp(2,ug,ur,+1,mpi_enreg,nfft,1,ngfft,0)
 ABI_FREE(ug)

end subroutine orb_to_real_one
!!***

!!****f* m_rdmft_integrals/rdmft_mem_report
!! NAME
!!  rdmft_mem_report
!!
!! FUNCTION
!!  Print the peak and the current resident memory of this process (Linux /proc/self/status), if rdmft_mem_on is set.
!!
!! SOURCE
subroutine rdmft_mem_report(label)

 character(len=*),intent(in) :: label
 character(len=200) :: line
 character(len=60) :: hwm,rss
 character(len=500) :: msg
 integer :: ios,unt

! *************************************************************************

 if (.not.rdmft_mem_on) return
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
 if (xmpi_comm_rank(xmpi_world) == 0) then
   write(msg,'(4a)') '     memory (master) ',trim(label),': peak '//trim(hwm)//', current ',trim(rss)
   call wrtout(std_out,msg)
 end if

end subroutine rdmft_mem_report
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
!!  comm = (optional) MPI communicator over which the pairs of orbitals are shared
!!  req = (optional) request of the exchange part of the generalized Fock matrices, built from the same pair densities
!!  rel(nbk) = (optional) relevant orbitals: the pairs with an orbital that is not relevant are not computed
!!  reduce = (optional) if true, use time reversal: only the independent blocks of k-points (a representative of each pair
!!    k, -k) are computed, the rest of the matrices J_pq, K_pq is filled by symmetry (J_pq = J_p,-q = J_-p,-q and
!!    K_-p,-q = K_pq). The free and background orbitals must be the same at k and -k.
!!  sym = (optional) classes of pairs of k-points under the symmetry operations of the ground state (rdmft_symtab_t): only
!!    the representative pair of each class is computed and the integrals of the pairs of the class are filled with the same
!!    value (replaces the time-reversal reduction). The free and background orbitals must be the same at the k-points
!!    related by the operations.
!!
!! OUTPUT
!!  ints
!!
!! SOURCE
subroutine rdmft_ints_compute(ints,orbs,ngfft,gsqcut,icutcoul_j,rcut_j,icutcoul_k,rcut_k,free,bg,core_ints,quiet,comm,req,&
                               rel,reduce,sym)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_ints_t),intent(inout) :: ints
 type(rdmft_orbs_t),intent(in) :: orbs
 integer,intent(in) :: icutcoul_j,icutcoul_k
 real(dp),intent(in) :: gsqcut,rcut_j,rcut_k
 logical,intent(in) :: core_ints
 logical,optional,intent(in) :: quiet
 integer,optional,intent(in) :: comm
 type(rdmft_fockreq_t),optional,intent(inout) :: req
 logical,optional,intent(in) :: rel(:)
 logical,optional,intent(in) :: reduce
 type(rdmft_symtab_t),optional,intent(in) :: sym
 !arrays
 integer,intent(in) :: ngfft(18)
 logical,intent(in) :: free(:),bg(:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: izero0=0
 real(dp) :: tph,tpk
 integer :: nblk,nloc,nrow
 integer,allocatable :: fk(:),fb(:)
 integer :: jlo,jhi,cnt,comm_,ib,ierr,ik,iq,ip,jb,jk,jp,me,nb,nbj,nfft,nproc
 logical :: canon,fused,needk,red,wantf,isff,anybg
 integer :: ipb,iqb,pk,cc,imem,ma,mb,ipm,iqm
 logical :: usesym,cselft
 integer,allocatable :: trk(:)
 real(dp) :: wgt
 real(dp),allocatable :: wacc(:,:,:),yk(:,:),wrs(:,:)
 logical :: wf
 real(dp) :: fcy
 real(dp) :: kval
 logical :: fp,fq,bp,bq
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg
 !arrays
 real(dp) :: qpt(3),gmet(3,3)
 real(dp),allocatable :: ui1(:,:),uj1(:,:),jtmp(:,:,:),vrho(:,:,:),prod(:,:),pgrid(:,:),rhog(:,:,:),rbg(:,:)
 real(dp),allocatable :: rtmp(:,:),vc(:),vq(:,:),ui(:,:,:),uj(:,:,:)

! *************************************************************************

 call rdmft_mem_report('ints: entry')
 call ints%free()
 comm_ = xmpi_world
 if (present(comm)) comm_ = comm
 nproc = xmpi_comm_size(comm_)
 me = xmpi_comm_rank(comm_)
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
 ! Use of time reversal: the orbital (n,-k) is the complex conjugate of (n,k)
 red = .false.
 if (present(reduce)) red = reduce
 usesym = present(sym)
 if (usesym) usesym = (sym%ncls > 0)
 if (usesym) red = .false.
 cselft = .false.
 ABI_MALLOC(trk,(orbs%nkpt))
 trk = orbs%tr_k
 if (red .and. any(trk == 0)) red = .false.
 if (red) then
   do ik=1,orbs%nkpt
     pk = trk(ik)
     if (pk == ik) cycle
     if (orbs%nband(pk,1) /= orbs%nband(ik,1)) then
       ABI_ERROR('RDMFT: different numbers of bands at k and -k (rdmft_tr_reduce = 1).')
     end if
     do ib=1,orbs%nband(ik,1)
       if ((free(ints%offset(ik)+ib) .neqv. free(ints%offset(pk)+ib)) .or. &
           (bg(ints%offset(ik)+ib) .neqv. bg(ints%offset(pk)+ib))) then
         ABI_ERROR('RDMFT: the active space is not the same at k and -k (use rdmft_tr_reduce = 0).')
       end if
     end do
   end do
 end if
 ABI_MALLOC(ints%jint,(ints%nfree,ints%nfree))
 ABI_MALLOC(ints%kint,(ints%nfree,ints%nfree))
 ABI_MALLOC(ints%bj,(ints%nfree))
 ABI_MALLOC(ints%bk,(ints%nfree))
 ints%jint = zero
 ints%kint = zero
 ints%bj = zero
 ints%bk = zero
 call rdmft_mem_report('ints: after jint, kint')

 call rdmft_mem_report('ints: before the J part')
 ABI_MALLOC(prod,(2,nfft))
 ABI_MALLOC(pgrid,(2,nfft))
 ABI_MALLOC(vc,(nfft))

 ! ------------------------------------------------------------------
 ! Coulomb integrals: rho_pp(G) of the free orbitals, sum of rho_qq(G) over the background orbitals
 ! ------------------------------------------------------------------
 ABI_MALLOC(rbg,(2,nfft))
 ABI_MALLOC(rtmp,(2,nfft))
 ABI_MALLOC(vq,(2,nfft))
 rbg = zero
 qpt = zero
 call barevcoul(rcut_j,icutcoul_j,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)
 vc(1) = zero   ! G=0 is the first element of the FFT mesh: neutralizing background

 ABI_MALLOC(ui1,(2,nfft))
 cnt = 0
 tph = abi_wtime()

 ! Background orbitals: sum of the densities (shared among the processes). The density of the orbital (n,-k) is the one
 ! of (n,k): the representative k-points are transformed (background orbitals count twice)
 do ik=1,orbs%nkpt
   if (red .and. trk(ik) < ik) cycle
   wgt = one
   if (red .and. trk(ik) > ik) wgt = two
   nb = orbs%nband(ik,1)
   do ib=1,nb
     ip = ints%offset(ik)+ib
     if (.not.bg(ip)) cycle
     cnt = cnt+1
     if (mod(cnt,nproc) /= me) cycle
     call orb_to_real_one(orbs,ik,ib,ngfft,mpi_enreg,ui1)
     prod(1,:) = ui1(1,:)**2+ui1(2,:)**2
     prod(2,:) = zero
     call fourdp(2,rtmp,prod,-1,mpi_enreg,nfft,1,ngfft,0)
     rbg = rbg+wgt*rtmp
     ints%sjd = ints%sjd+wgt*sum(vc(:)*(rtmp(1,:)**2+rtmp(2,:)**2))/orbs%ucvol
   end do
 end do
 rdmft_tint(1) = rdmft_tint(1)+abi_wtime()-tph
 tph = abi_wtime()
 if (nproc > 1) then
   call xmpi_sum(rbg,comm_,ierr)
   call xmpi_sum(ints%sjd,comm_,ierr)
 end if
 rdmft_tint(2) = rdmft_tint(2)+abi_wtime()-tph
 tph = abi_wtime()

 ! J_pq = Re sum_G conjg(rho_p) v rho_q / Omega. The densities of the free orbitals are never stored all at once: the columns
 ! q of this process (a block of the free orbitals) are kept with the Coulomb kernel applied (vrho), the densities of the
 ! rows p are computed block by block (by all the processes: a density costs one FFT, much less than the exchange
 ! integrals, and no communication of the densities is needed).
 if (ints%nfree > 0) then
   ! k-point and band of each free orbital; the representative of its time-reversal pair (the same density)
   ABI_MALLOC(fk,(ints%nfree))
   ABI_MALLOC(fb,(ints%nfree))
   do iq=1,ints%nfree
     jp = ints%comp_of(iq)
     do ik=orbs%nkpt,1,-1
       if (jp > ints%offset(ik)) exit
     end do
     fk(iq) = ik
     fb(iq) = jp-ints%offset(ik)
     if (red .and. trk(ik) < ik) fk(iq) = trk(ik)
   end do
   jlo = (me*ints%nfree)/nproc+1
   jhi = ((me+1)*ints%nfree)/nproc
   nloc = max(jhi-jlo+1,0)
   ABI_MALLOC(vrho,(2,nfft,max(nloc,1)))
   do iq=jlo,jhi
     call orb_to_real_one(orbs,fk(iq),fb(iq),ngfft,mpi_enreg,ui1)
     prod(1,:) = ui1(1,:)**2+ui1(2,:)**2
     prod(2,:) = zero
     call fourdp(2,rtmp,prod,-1,mpi_enreg,nfft,1,ngfft,0)
     vrho(1,:,iq-jlo+1) = vc(:)*rtmp(1,:)
     vrho(2,:,iq-jlo+1) = vc(:)*rtmp(2,:)
     ints%bj(iq) = sum(vrho(1,:,iq-jlo+1)*rbg(1,:)+vrho(2,:,iq-jlo+1)*rbg(2,:))/orbs%ucvol
   end do
   nblk = min(ints%nfree,128)
   ABI_MALLOC(rhog,(2,nfft,nblk))
   ABI_MALLOC(jtmp,(2,nblk,max(nloc,1)))
   do jb=1,ints%nfree,nblk
     nrow = min(nblk,ints%nfree-jb+1)
     do ib=1,nrow
       call orb_to_real_one(orbs,fk(jb+ib-1),fb(jb+ib-1),ngfft,mpi_enreg,ui1)
       prod(1,:) = ui1(1,:)**2+ui1(2,:)**2
       prod(2,:) = zero
       call fourdp(2,rhog(:,:,ib),prod,-1,mpi_enreg,nfft,1,ngfft,0)
     end do
     if (nloc > 0) then
       call zgemm('C','N',nrow,nloc,nfft,cone,rhog,nfft,vrho,nfft,czero,jtmp,nblk)
       ints%jint(jb:jb+nrow-1,jlo:jhi) = jtmp(1,1:nrow,1:nloc)/orbs%ucvol
     end if
   end do
   ABI_FREE(rhog)
   ABI_FREE(jtmp)
   ABI_FREE(vrho)
   ABI_FREE(fk)
   ABI_FREE(fb)
   if (nproc > 1) then
     call xmpi_sum(ints%bj,comm_,ierr)
     call xmpi_sum(ints%jint,comm_,ierr)
   end if
   ints%jint = half*(ints%jint+transpose(ints%jint))
 end if
 ABI_FREE(ui1)
 rdmft_tint(3) = rdmft_tint(3)+abi_wtime()-tph
 call rdmft_mem_report('ints: after the J part')
 ! Background-background pairs: |sum_q rho_q|^2 contains the diagonal terms
 ints%sjo = sum(vc(:)*(rbg(1,:)**2+rbg(2,:)**2))/orbs%ucvol-ints%sjd
 ABI_FREE(rbg)
 ABI_FREE(rtmp)
 ABI_FREE(vq)

 ! ------------------------------------------------------------------
 ! Exchange integrals: pair densities of the needed pairs (n,k), (m,k'), k <= k'. If a request of the exchange part of
 ! the generalized Fock matrix is present, the same pair densities are used for the k-points asked (all the k').
 ! The pairs (k,k',q) are shared among the processes of comm.
 ! ------------------------------------------------------------------
 fused = present(req)
 if (fused) then
   ABI_MALLOC(yk,(2,nfft))
   ABI_MALLOC(wrs,(2,nfft))
   if (allocated(req%fx)) then
     req%fx = czero
   end if
 end if
 ABI_MALLOC(uj1,(2,nfft))
 cnt = 0
 tph = abi_wtime()
 do ik=1,orbs%nkpt
   nb = orbs%nband(ik,1)
   if (.not.any(free(ints%offset(ik)+1:ints%offset(ik)+nb) .or. bg(ints%offset(ik)+1:ints%offset(ik)+nb))) cycle
   wantf = .false.
   if (fused) wantf = req%want(ik)
   ABI_MALLOC(ui,(2,nfft,nb))
   if (wantf) then
     ABI_MALLOC(wacc,(2,nfft,nb))
     wacc = zero
   end if
   tpk = abi_wtime()
   call orbs_to_real(orbs,ik,ngfft,mpi_enreg,ui)
   rdmft_tprof(1) = rdmft_tprof(1)+abi_wtime()-tpk
   do jk=1,orbs%nkpt
     nbj = orbs%nband(jk,1)
     if (.not.any(free(ints%offset(jk)+1:ints%offset(jk)+nbj) .or. bg(ints%offset(jk)+1:ints%offset(jk)+nbj))) cycle
     ! Independent blocks (k,k') of the exchange integrals between free orbitals. Without time reversal: k' >= k. With it
     ! (the orbit of a block under the transposition and k,k' -> -k,-k'): k is a representative; k' is a representative
     ! with k' >= k, or the partner of a representative k'' >= k (and then k is not time-reversal invariant).
     canon = (jk >= ik)
     if (red) then
       canon = .false.
       if (trk(ik) >= ik) then
         if (trk(jk) >= jk) then
           canon = (jk >= ik)
         else
           canon = (trk(jk) >= ik .and. trk(ik) /= ik)
         end if
       end if
     end if
     ! Pairs with a background orbital: k' >= k
     anybg = any(bg(ints%offset(ik)+1:ints%offset(ik)+nb)) .or. any(bg(ints%offset(jk)+1:ints%offset(jk)+nbj))
     if (usesym) then
       ! Symmetry classes: the representative pair of a class that is computed (not the transposed of a smaller class)
       cc = sym%cls(ik,jk)
       canon = (sym%rep(1,cc) == ik .and. sym%rep(2,cc) == jk .and. sym%use(cc))
       cselft = sym%selft(cc)
       if (.not.(canon .or. wantf)) cycle
     else
       if (.not.(canon .or. wantf .or. (anybg .and. jk >= ik))) cycle
     end if
     !Coulomb kernel for q = k' - k
     qpt(:) = orbs%kpt(:,jk)-orbs%kpt(:,ik)
     call barevcoul(rcut_k,icutcoul_k,qpt,gsqcut,gmet,nfft,orbs%nkpt,ngfft,orbs%ucvol,izero0,vc)

     do jb=1,nbj
       iq = ints%offset(jk)+jb
       fq = free(iq); bq = bg(iq)
       if (.not.(fq .or. bq)) cycle
       if (present(rel)) then
         if (.not.rel(iq)) cycle
       end if
       cnt = cnt+1
       if (mod(cnt,nproc) /= me) cycle
       tpk = abi_wtime()
       call orb_to_real_one(orbs,jk,jb,ngfft,mpi_enreg,uj1)
       rdmft_tprof(2) = rdmft_tprof(2)+abi_wtime()-tpk
       do ib=1,nb
         ip = ints%offset(ik)+ib
         fp = free(ip); bp = bg(ip)
         isff = (fp .and. fq)
         if (usesym) then
           ! The band pairs (n,m) and (m,n) of a representative pair that is its own transposed are equivalent
           needk = canon .and. (isff .or. fp .or. bp)
           if (cselft .and. jb < ib) needk = .false.
           if (.not.isff .and. bp .and. bq .and. .not.core_ints) needk = .false.
         else if (isff) then
           needk = canon
           if (jb > ib .and. (jk == ik .or. (red .and. jk == trk(ik)))) needk = .false.
         else
           needk = (jk >= ik) .and. (fp .or. bp)
           if (jk == ik .and. jb > ib) needk = .false.
           if (bp .and. bq .and. .not.core_ints) needk = .false.
         end if
         if (present(rel) .and. needk) needk = rel(ip)
         ! Exchange part of the Fock matrix of this pair: only if the coupling is not zero
         wf = wantf
         if (wf) then
           wf = (req%rowmap(ib,ik) > 0 .and. req%colmap(iq) > 0)
           if (wf) wf = (req%ccoef(req%rowmap(ib,ik),req%colmap(iq)) /= zero)
         end if
         if (.not.(needk .or. wf)) cycle
         ! conjg(u_i) u_j
         tpk = abi_wtime()
         prod(1,:) = ui(1,:,ib)*uj1(1,:)+ui(2,:,ib)*uj1(2,:)
         prod(2,:) = ui(1,:,ib)*uj1(2,:)-ui(2,:,ib)*uj1(1,:)
         rdmft_tprof(3) = rdmft_tprof(3)+abi_wtime()-tpk
         tpk = abi_wtime()
         call fourdp(2,pgrid,prod,-1,mpi_enreg,nfft,1,ngfft,0)
         rdmft_tprof(4) = rdmft_tprof(4)+abi_wtime()-tpk
         rdmft_tprof(7) = rdmft_tprof(7)+one
         if (wf) then
           ! F_pj += sum_G conjg(M_pq) Y_jq with Y_jq = -2 c_jq/Omega v M_jq. By Parseval (full FFT mesh):
           ! F_pj = 1/nfft sum_r u_p(r) W_j(r), W_j(r) = sum_q conjg(u_q(r)) b_jq(r), b_jq = inverse FFT of Y_jq
           ! (accumulated over k' and q, one ZGEMM per k-point at the end)
           tpk = abi_wtime()
           fcy = -two*req%ccoef(req%rowmap(ib,ik),req%colmap(iq))/orbs%ucvol
           yk(1,:) = fcy*vc(:)*pgrid(1,:)
           yk(2,:) = fcy*vc(:)*pgrid(2,:)
           call fourdp(2,yk,wrs,+1,mpi_enreg,nfft,1,ngfft,0)
           wacc(1,:,ib) = wacc(1,:,ib)+uj1(1,:)*wrs(1,:)+uj1(2,:)*wrs(2,:)
           wacc(2,:,ib) = wacc(2,:,ib)+uj1(1,:)*wrs(2,:)-uj1(2,:)*wrs(1,:)
           rdmft_tprof(6) = rdmft_tprof(6)+abi_wtime()-tpk
         end if
         if (.not.needk) cycle
         tpk = abi_wtime()
         kval = sum(vc(:)*(pgrid(1,:)**2+pgrid(2,:)**2))/orbs%ucvol
         if (usesym) then
           ! The same value for all the pairs of k-points of the class (and the transposed pairs)
           do imem=sym%mptr(cc),sym%mptr(cc+1)-1
             ma = sym%mem(1,imem)
             mb = sym%mem(2,imem)
             ipm = ints%offset(ma)+ib
             iqm = ints%offset(mb)+jb
             if (isff) then
               ints%kint(ints%free_of(ipm),ints%free_of(iqm)) = kval
               ints%kint(ints%free_of(iqm),ints%free_of(ipm)) = kval
               cycle
             end if
             ! sums over pairs: the pair (ma,ib),(mb,ib) with ma > mb is the same as the one with the members exchanged
             if (cselft .and. ib == jb .and. ma > mb) cycle
             if (fp .and. bq) then
               ints%bk(ints%free_of(ipm)) = ints%bk(ints%free_of(ipm))+kval
             else if (bp .and. fq) then
               ints%bk(ints%free_of(iqm)) = ints%bk(ints%free_of(iqm))+kval
             else if (ipm == iqm) then
               ints%skd = ints%skd+kval
             else
               ints%sko = ints%sko+two*kval
             end if
           end do
         else if (isff) then
           ints%kint(ints%free_of(ip),ints%free_of(iq)) = kval
           ints%kint(ints%free_of(iq),ints%free_of(ip)) = kval
           if (red) then
             ! image under time reversal: K_{-p,-q} = K_{pq}
             ipb = ints%offset(trk(ik))+ib
             iqb = ints%offset(trk(jk))+jb
             ints%kint(ints%free_of(ipb),ints%free_of(iqb)) = kval
             ints%kint(ints%free_of(iqb),ints%free_of(ipb)) = kval
           end if
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
         rdmft_tprof(5) = rdmft_tprof(5)+abi_wtime()-tpk
       end do
     end do
   end do
   if (wantf) then
     tpk = abi_wtime()
     call zgemm('T','N',nb,nb,nfft,cmplx(one/real(nfft,dp),zero,kind=dpc),ui,nfft,wacc,nfft,cone,req%fx(1,1,ik),&
                size(req%fx,1))
     rdmft_tprof(6) = rdmft_tprof(6)+abi_wtime()-tpk
     ABI_FREE(wacc)
   end if
   ABI_FREE(ui)
 end do
 rdmft_tint(4) = rdmft_tint(4)+abi_wtime()-tph
 call rdmft_mem_report('ints: after the K sweep')
 tph = abi_wtime()
 if (nproc > 1) call xmpi_barrier(comm_)
 rdmft_tint(5) = rdmft_tint(5)+abi_wtime()-tph
 tph = abi_wtime()
 if (nproc > 1) then
   call xmpi_sum(ints%kint,comm_,ierr)
   call xmpi_sum(ints%bk,comm_,ierr)
   call xmpi_sum(ints%sko,comm_,ierr)
   call xmpi_sum(ints%skd,comm_,ierr)
   if (fused) call xmpi_sum(req%fx,comm_,ierr)
 end if
 rdmft_tint(6) = rdmft_tint(6)+abi_wtime()-tph
 call rdmft_mem_report('ints: after the sums of the K results')
 ABI_FREE(uj1)
 ABI_FREE(trk)
 if (fused) then
   ABI_FREE(yk)
   ABI_FREE(wrs)
 end if

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
