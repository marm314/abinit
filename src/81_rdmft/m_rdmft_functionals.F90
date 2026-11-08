!!****m* ABINIT/m_rdmft_functionals
!! NAME
!!  m_rdmft_functionals
!!
!! FUNCTION
!!  Energy functionals of the Reduced Density Matrix Functional Theory (RDMFT) for fixed Bloch orbitals,
!!  written for closed-shell (spin-restricted) occupations n_p in [0,1] per spin, p=(n,k) in the full BZ.
!!  With the cell-normalized integrals J_pq, K_pq (see m_rdmft_integrals), the electronic energy per unit cell is
!!    JK-only functionals:  E = (1/N_k) sum_p 2 n_p h_pp + (1/N_k^2) sum_pq [ 2 f_H(n_p,n_q) J_pq - f_X(n_p,n_q) K_pq ]
!!    PNOF functionals:     E = (1/N_k) sum_a 2 n_a h_aa + (1/N_k^2) { sum_a n_a J_aa
!!                              + sum_{a/=b} [ H_ab J_ab - X_ab K_ab + Pi_ab K_ab ] }
!!  where a,b run over the pairs of orbitals (geminals) of the subspaces. Time reversal gives L_pq=K_pq.
!!  The constant terms (Ewald, pseudopotential core) are not included. The expressions are those of ReRDMFT
!!  (non-relativistic branch); see RDMFT/functionals_kpoints.tex.
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

module m_rdmft_functionals

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private
!!***

 ! Identifiers of the JK-only functionals
 integer,parameter,public :: RDMFT_SD        = 1   ! Single determinant (Hartree-Fock)
 integer,parameter,public :: RDMFT_MBB       = 2   ! Muller / Buijse-Baerends
 integer,parameter,public :: RDMFT_CA        = 3   ! Csanyi-Arias
 integer,parameter,public :: RDMFT_RHFB      = 4   ! Hartree-Fock-Bogoliubov (attractive interactions)
 integer,parameter,public :: RDMFT_CGA       = 5   ! Csanyi-Goedecker-Arias
 integer,parameter,public :: RDMFT_ML        = 6   ! Marques-Lathiotakis
 integer,parameter,public :: RDMFT_MLSIC     = 7   ! Marques-Lathiotakis, self-interaction corrected
 integer,parameter,public :: RDMFT_GU        = 8   ! Goedecker-Umrigar
 integer,parameter,public :: RDMFT_POWER     = 9   ! (n_p n_q)^alpha
 integer,parameter,public :: RDMFT_BBC2      = 10  ! Buijse-Baerends corrected (BBC2)
 integer,parameter,public :: RDMFT_MULLER_AS = 11  ! Antisymmetrized Muller (not HF at integer occupations)
 ! Identifiers of the PNOF functionals
 integer,parameter,public :: RDMFT_PNOF5     = 21
 integer,parameter,public :: RDMFT_PNOF7     = 22
 integer,parameter,public :: RDMFT_PNOF7S    = 23
 integer,parameter,public :: RDMFT_GNOF      = 24

 ! Fitted coefficients of the Marques-Lathiotakis functionals
 real(dp),parameter :: ml_a0 = 126.3101_dp, ml_a1 = 2213.33_dp, ml_b1 = 2338.64_dp
 real(dp),parameter :: mlsic_a0 = 1298.78_dp, mlsic_a1 = 35114.4_dp, mlsic_b1 = 36412.2_dp
 ! GNOF: 0.02*sqrt(2)
 real(dp),parameter :: gnof_c = 0.02828427124746190_dp

!!****t* m_rdmft_functionals/rdmft_pnof_table_t
!! NAME
!! rdmft_pnof_table_t
!!
!! FUNCTION
!!  Pairs of orbitals (geminals) taking part in a PNOF energy, and their coupling subspaces.
!!  A subspace has one principal (occupied) pair and coupled virtual pairs. The deep occupied pairs are
!!  trivial subspaces made of their own principal pair. The pairs may belong to different k-points.
!!
!! SOURCE

 type,public :: rdmft_pnof_table_t

   integer :: ng = 0
   ! Number of pairs.

   integer :: nsub = 0
   ! Number of subspaces (including the trivial ones).

   integer,allocatable :: gem_orb(:)
   ! gem_orb(ng)
   ! Index of the orbital p=(n,k) (composite index) of each pair.

   integer,allocatable :: gem_sub(:)
   ! gem_sub(ng)
   ! Subspace of each pair.

   logical,allocatable :: gem_principal(:)
   ! gem_principal(ng)
   ! True for the principal pair of its subspace.

 contains

   procedure :: free => rdmft_pnof_table_free
   procedure :: build => rdmft_pnof_table_build

 end type rdmft_pnof_table_t
!!***

 public :: rdmft_energy_jk
 public :: rdmft_energy_pnof
 public :: rdmft_functional_name
!!***

contains
!!***

!!****f* m_rdmft_functionals/rdmft_functional_name
!! NAME
!!  rdmft_functional_name
!!
!! FUNCTION
!!  Name of a functional from its identifier.
!!
!! SOURCE
function rdmft_functional_name(func) result(name)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 character(len=10) :: name

! *************************************************************************

 select case (func)
 case (RDMFT_SD);        name = 'SD'
 case (RDMFT_MBB);       name = 'MBB'
 case (RDMFT_CA);        name = 'CA'
 case (RDMFT_RHFB);      name = 'RHFB'
 case (RDMFT_CGA);       name = 'CGA'
 case (RDMFT_ML);        name = 'ML'
 case (RDMFT_MLSIC);     name = 'MLSIC'
 case (RDMFT_GU);        name = 'GU'
 case (RDMFT_POWER);     name = 'POWER'
 case (RDMFT_BBC2);      name = 'BBC2'
 case (RDMFT_MULLER_AS); name = 'MULLER_AS'
 case (RDMFT_PNOF5);     name = 'PNOF5'
 case (RDMFT_PNOF7);     name = 'PNOF7'
 case (RDMFT_PNOF7S);    name = 'PNOF7S'
 case (RDMFT_GNOF);      name = 'GNOF'
 case default;           name = 'UNKNOWN'
 end select

end function rdmft_functional_name
!!***

!!****f* m_rdmft_functionals/jk_exchange
!! NAME
!!  jk_exchange
!!
!! FUNCTION
!!  Exchange function f_X(n_i,n_j) of the JK-only functionals.
!!
!! INPUTS
!!  func = functional
!!  ni,nj = occupation numbers (per spin)
!!  same = true if i and j are the same orbital
!!  weak_i,weak_j = true if the orbitals are weakly occupied (only used by BBC2)
!!  alpha = exponent of the POWER functional
!!
!! SOURCE
function jk_exchange(func,ni,nj,same,weak_i,weak_j,alpha) result(fx)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 real(dp),intent(in) :: ni,nj,alpha
 logical,intent(in) :: same,weak_i,weak_j
 real(dp) :: fx

 !Local variables-------------------------------
 !scalars
 real(dp) :: x

! *************************************************************************

 x = ni*nj
 select case (func)
 case (RDMFT_SD)
   fx = x
 case (RDMFT_MBB)
   fx = sqrt(max(zero,x))
 case (RDMFT_CA)
   fx = sqrt(max(zero,ni*(one-ni)*nj*(one-nj)))+x
 case (RDMFT_RHFB)
   fx = -sqrt(max(zero,ni*(one-ni)*nj*(one-nj)))+x
 case (RDMFT_CGA)
   fx = half*(x+sqrt(max(zero,ni*(two-ni)*nj*(two-nj))))
 case (RDMFT_ML)
   fx = x*(ml_a0+ml_a1*x)/(one+ml_b1*x)
 case (RDMFT_MLSIC)
   if (same) then
     fx = x
   else
     fx = x*(mlsic_a0+mlsic_a1*x)/(one+mlsic_b1*x)
   end if
 case (RDMFT_GU)
   if (same) then
     fx = x
   else
     fx = sqrt(max(zero,x))
   end if
 case (RDMFT_POWER)
   fx = x**alpha
 case (RDMFT_BBC2)
   if (same) then
     fx = ni
   else if (weak_i .and. weak_j) then
     fx = -sqrt(max(zero,x))
   else if ((.not.weak_i) .and. (.not.weak_j)) then
     fx = x
   else
     fx = sqrt(max(zero,x))
   end if
 case (RDMFT_MULLER_AS)
   fx = half*(x-sqrt(max(zero,x)))
 case default
   ABI_ERROR('RDMFT: unknown JK-only functional.')
 end select

end function jk_exchange
!!***

!!****f* m_rdmft_functionals/rdmft_energy_jk
!! NAME
!!  rdmft_energy_jk
!!
!! FUNCTION
!!  Electronic energy per unit cell (without Ewald and core terms) of a JK-only functional.
!!
!! INPUTS
!!  func = functional (RDMFT_SD ... RDMFT_MULLER_AS)
!!  nk = number of k-points of the full BZ (all with weight 1/nk)
!!  n(nbk) = occupation numbers per spin of all the orbitals (composite index)
!!  h1(nbk) = diagonal elements of the core Hamiltonian
!!  jint(nbk,nbk), kint(nbk,nbk) = J and K integrals (cell-normalized)
!!  weak(nbk) = true for the weakly occupied orbitals (BBC2)
!!  alpha = exponent of the POWER functional
!!
!! SOURCE
function rdmft_energy_jk(func,nk,nbk,n,h1,jint,kint,weak,alpha) result(energy)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk
 real(dp),intent(in) :: alpha
 real(dp) :: energy
 !arrays
 real(dp),intent(in) :: n(nbk),h1(nbk),jint(nbk,nbk),kint(nbk,nbk)
 logical,intent(in) :: weak(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ip,iq
 real(dp) :: e1,e2,fh,fx

! *************************************************************************

 e1 = zero
 e2 = zero
 do ip=1,nbk
   e1 = e1+two*n(ip)*h1(ip)
   do iq=1,nbk
     fx = jk_exchange(func,n(ip),n(iq),ip==iq,weak(ip),weak(iq),alpha)
     fh = n(ip)*n(iq)
     if (func == RDMFT_MULLER_AS) fh = fx
     e2 = e2+two*fh*jint(ip,iq)-fx*kint(ip,iq)
   end do
 end do
 energy = e1/nk+e2/real(nk,dp)**2

end function rdmft_energy_jk
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_table_free
!! NAME
!!  rdmft_pnof_table_free
!!
!! FUNCTION
!!  Free the memory of a table of pairs.
!!
!! SOURCE
subroutine rdmft_pnof_table_free(table)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_pnof_table_t),intent(inout) :: table

! *************************************************************************

 ABI_SFREE(table%gem_orb)
 ABI_SFREE(table%gem_sub)
 ABI_SFREE(table%gem_principal)
 table%ng = 0
 table%nsub = 0

end subroutine rdmft_pnof_table_free
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_table_build
!! NAME
!!  rdmft_pnof_table_build
!!
!! FUNCTION
!!  Build the pairs and subspaces from the orbital energies, over all k-points:
!!  the nsub highest occupied orbitals are the principal pairs of nsub subspaces, each one coupled with the
!!  (ncoupl-1) lowest unused virtual orbitals (the first subspace takes the lowest virtuals, and so on). Pairs
!!  of different k-points can share a subspace. The remaining occupied orbitals are trivial subspaces.
!!  The closure of the subspaces under k -> -k is not enforced here.
!!
!! INPUTS
!!  nbk = number of orbitals
!!  eig(nbk) = orbital energies
!!  occupied(nbk) = true for the (strongly) occupied orbitals
!!  nsub = number of coupling subspaces
!!  ncoupl = number of pairs of each subspace (1 principal + ncoupl-1 virtual)
!!
!! OUTPUT
!!  ierr = 0 if ok, 1 if there are not enough occupied or virtual orbitals
!!
!! SOURCE
subroutine rdmft_pnof_table_build(table,nbk,eig,occupied,nsub,ncoupl,ierr)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_pnof_table_t),intent(inout) :: table
 integer,intent(in) :: nbk,nsub,ncoupl
 integer,intent(out) :: ierr
 !arrays
 real(dp),intent(in) :: eig(nbk)
 logical,intent(in) :: occupied(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ip,ig,is,iv,jp,nocc,nvirt
 !arrays
 logical,allocatable :: used(:)

! *************************************************************************

 call table%free()
 ierr = 0
 nocc = count(occupied)
 nvirt = nbk-nocc
 if (nsub < 1 .or. ncoupl < 2 .or. nsub > nocc .or. nsub*(ncoupl-1) > nvirt) then
   ierr = 1
   return
 end if

 table%nsub = nocc-nsub+nsub   ! deep occupied (own subspace) + coupling subspaces
 table%ng = nocc+nsub*(ncoupl-1)
 ABI_MALLOC(table%gem_orb,(table%ng))
 ABI_MALLOC(table%gem_sub,(table%ng))
 ABI_MALLOC(table%gem_principal,(table%ng))
 ABI_MALLOC(used,(nbk))
 used = .false.

 ig = 0
 ! Principal pairs: the nsub highest occupied orbitals; subspaces 1..nsub
 do is=1,nsub
   jp = 0
   do ip=1,nbk
     if (occupied(ip) .and. .not.used(ip)) then
       if (jp == 0) then
         jp = ip
       else if (eig(ip) > eig(jp)) then
         jp = ip
       end if
     end if
   end do
   used(jp) = .true.
   ig = ig+1
   table%gem_orb(ig) = jp
   table%gem_sub(ig) = is
   table%gem_principal(ig) = .true.
 end do
 ! Coupled virtual pairs: the lowest unused virtual orbitals, block after block
 do is=1,nsub
   do iv=1,ncoupl-1
     jp = 0
     do ip=1,nbk
       if (.not.occupied(ip) .and. .not.used(ip)) then
         if (jp == 0) then
           jp = ip
         else if (eig(ip) < eig(jp)) then
           jp = ip
         end if
       end if
     end do
     used(jp) = .true.
     ig = ig+1
     table%gem_orb(ig) = jp
     table%gem_sub(ig) = is
     table%gem_principal(ig) = .false.
   end do
 end do
 ! Deep occupied orbitals: trivial subspaces
 do ip=1,nbk
   if (occupied(ip) .and. .not.used(ip)) then
     ig = ig+1
     table%gem_orb(ig) = ip
     table%gem_sub(ig) = nsub+(ig-nsub*ncoupl)
     table%gem_principal(ig) = .true.
   end if
 end do
 table%nsub = nsub+(nocc-nsub)
 ABI_FREE(used)

end subroutine rdmft_pnof_table_build
!!***

!!****f* m_rdmft_functionals/rdmft_energy_pnof
!! NAME
!!  rdmft_energy_pnof
!!
!! FUNCTION
!!  Electronic energy per unit cell (without Ewald and core terms) of a PNOF functional.
!!
!! INPUTS
!!  func = RDMFT_PNOF5, RDMFT_PNOF7, RDMFT_PNOF7S or RDMFT_GNOF
!!  table <type(rdmft_pnof_table_t)> = pairs and subspaces
!!  nk = number of k-points of the full BZ
!!  n(nbk) = occupation numbers per spin of all the orbitals (only those in the table are used)
!!  h1(nbk), jint(nbk,nbk), kint(nbk,nbk) = core Hamiltonian diagonal, J and K integrals
!!  selfint = rdmft_selfint: 0 for the self-interaction free diagonal term n_a J_aa, 1 to add n_a^2 (J_aa-K_aa)
!!
!! SOURCE
function rdmft_energy_pnof(func,table,nk,nbk,n,h1,jint,kint,selfint) result(energy)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk,selfint
 type(rdmft_pnof_table_t),intent(in) :: table
 real(dp) :: energy
 !arrays
 real(dp),intent(in) :: n(nbk),h1(nbk),jint(nbk,nbk),kint(nbk,nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ia,ib,ig,ip,iq
 real(dp) :: cross,e1,e2,hab,hha,hhb,na,nb,nda,ndb,pi_ab,xab
 logical :: same_sub
 !arrays
 real(dp),allocatable :: hprin(:)

! *************************************************************************

 ! Holes of the principal pair of each subspace (needed by GNOF)
 ABI_MALLOC(hprin,(table%nsub))
 hprin = one
 do ig=1,table%ng
   if (table%gem_principal(ig)) hprin(table%gem_sub(ig)) = one-n(table%gem_orb(ig))
 end do

 e1 = zero
 e2 = zero
 do ia=1,table%ng
   ip = table%gem_orb(ia)
   na = n(ip)
   e1 = e1+two*na*h1(ip)
   e2 = e2+na*jint(ip,ip)
   ! rdmft_selfint=1: same-spin self term of the orbital, n^2 (J_pp-K_pp), nonzero when K_pp /= J_pp (periodic systems)
   if (selfint == 1) e2 = e2+na**2*(jint(ip,ip)-kint(ip,ip))
   do ib=1,table%ng
     if (ib == ia) cycle
     iq = table%gem_orb(ib)
     nb = n(iq)
     same_sub = (table%gem_sub(ia) == table%gem_sub(ib))
     xab = na*nb
     if (same_sub) xab = zero
     ! Hartree and exchange weights: H_ab=2 X_ab with X_ab = n_a n_b outside the subspace
     hab = two*xab
     ! Pairing term Pi_ab
     hha = one-na
     hhb = one-nb
     if (same_sub) then
       pi_ab = sqrt(max(zero,na*nb))
       if (table%gem_principal(ia) .or. table%gem_principal(ib)) pi_ab = -pi_ab
     else
       select case (func)
       case (RDMFT_PNOF5)
         pi_ab = zero
       case (RDMFT_PNOF7)
         pi_ab = -sqrt(max(zero,na*nb*hha*hhb))
       case (RDMFT_PNOF7S)
         pi_ab = -four*na*nb*hha*hhb
       case (RDMFT_GNOF)
         if (table%gem_principal(ia) .and. table%gem_principal(ib)) then
           pi_ab = zero
         else
           nda = na*exp(-(hprin(table%gem_sub(ia))/gnof_c)**2)
           ndb = nb*exp(-(hprin(table%gem_sub(ib))/gnof_c)**2)
           cross = sqrt(max(zero,nda*ndb))
           pi_ab = nda*ndb-sqrt(max(zero,na*nb*hha*hhb))
           if (table%gem_principal(ia) .neqv. table%gem_principal(ib)) then
             pi_ab = pi_ab-cross
           else
             pi_ab = pi_ab+cross
           end if
         end if
       case default
         ABI_ERROR('RDMFT: unknown PNOF functional.')
       end select
     end if
     e2 = e2+hab*jint(ip,iq)-xab*kint(ip,iq)+pi_ab*kint(ip,iq)
   end do
 end do
 energy = e1/nk+e2/real(nk,dp)**2

 ABI_FREE(hprin)

end function rdmft_energy_pnof
!!***

end module m_rdmft_functionals
!!***
