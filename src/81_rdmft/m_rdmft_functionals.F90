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
 integer,parameter,public :: RDMFT_MULLER_AS = 10  ! Antisymmetrized Muller (not HF at integer occupations)
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

!!****t* m_rdmft_functionals/rdmft_bg_t
!! NAME
!! rdmft_bg_t
!!
!! FUNCTION
!!  Constants of the background orbitals (frozen occupied, occupation 1): their one-electron energy and the sums of
!!  the Coulomb and exchange integrals over their pairs, which are needed for the total energy only.
!!
!! SOURCE

 type,public :: rdmft_bg_t

   real(dp) :: sh = zero
   ! Sum of h_pp over the background orbitals.

   real(dp) :: sjo = zero, sko = zero
   ! Sums of J_pq and K_pq over the ordered pairs p/=q of background orbitals.

   real(dp) :: sjd = zero, skd = zero
   ! Sums of J_pp and K_pp over the background orbitals.

 end type rdmft_bg_t
!!***

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

   integer :: nclass = 0
   ! Number of classes (principal star with its coupled virtual stars). A class has one occupation variable per
   ! star (the principal one and nvirt virtual ones) and nvirt angles gamma.

   integer :: nvar = 0
   ! Number of occupation variables (principal stars and virtual stars of all the classes).

   integer :: ngam = 0
   ! Number of angles gamma (the free parameters of the PNOF optimization).

   integer,allocatable :: class_var(:)
   ! class_var(nclass)
   ! Index of the principal variable of each class; its virtual variables follow.

   integer,allocatable :: class_nvirt(:)
   ! class_nvirt(nclass)
   ! Number of virtual stars of each class.

   integer,allocatable :: var_mult(:)
   ! var_mult(nvar)
   ! Number of orbitals of the same subspace sharing the variable (2 for a generic-k virtual star in the subspace of a
   ! time-reversal-invariant principal pair, else 1).

   integer,allocatable :: orb_var(:)
   ! orb_var(nbk)
   ! Occupation variable of each orbital (0 for the frozen orbitals); the time-reversal partners share the variable.

 contains

   procedure :: free => rdmft_pnof_table_free
   procedure :: build => rdmft_pnof_table_build
   procedure :: build_sym => rdmft_pnof_table_build_sym

 end type rdmft_pnof_table_t
!!***

 public :: rdmft_active_space
 public :: rdmft_coef_jk
 public :: rdmft_coef_pnof
 public :: rdmft_energy_bg
 public :: rdmft_fermi_dirac
 public :: rdmft_pnof_gamma
 public :: rdmft_pnof_table_remap
 public :: rdmft_grad_pnof
 public :: rdmft_pnof_occ
 public :: rdmft_jk_diag
 public :: rdmft_jk_pair
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
!!  alpha = exponent of the POWER functional
!!
!! SOURCE
function jk_exchange(func,ni,nj,same,alpha) result(fx)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 real(dp),intent(in) :: ni,nj,alpha
 logical,intent(in) :: same
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
 case (RDMFT_MULLER_AS)
   fx = half*(x-sqrt(max(zero,x)))
 case default
   ABI_ERROR('RDMFT: unknown JK-only functional.')
 end select

end function jk_exchange
!!***

!!****f* m_rdmft_functionals/fx_pair
!! NAME
!!  fx_pair
!!
!! FUNCTION
!!  Exchange function f_X(a,b) of two different orbitals and its derivatives with respect to the first argument
!!  (d1, d11) and the mixed second derivative (d12).
!!
!! SOURCE
subroutine fx_pair(func,a,b,alpha,f,d1,d11,d12)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 real(dp),intent(in) :: a,b,alpha
 real(dp),intent(out) :: f,d1,d11,d12

 !Local variables-------------------------------
 !scalars
 real(dp) :: dd,gp,gpp,mm,phi,phip,psi,psip,s,sp,sq,su,sv,x
 real(dp),parameter :: tiny_den = 1.0e-300_dp

! *************************************************************************

 x = a*b
 s = sqrt(max(zero,x))
 select case (func)
 case (RDMFT_SD)
   f = x; d1 = b; d11 = zero; d12 = one
 case (RDMFT_MBB, RDMFT_GU)
   f = s; d1 = half*sqrt(max(zero,b/a)); d11 = -0.25_dp*sqrt(max(zero,b))*a**(-1.5_dp)
   d12 = zero; if (s > tiny_den) d12 = 0.25_dp/s
 case (RDMFT_CA, RDMFT_RHFB)
   su = sqrt(max(zero,a*(one-a))); sv = sqrt(max(zero,b*(one-b)))
   phi = (one-two*a)/(two*su)
   phip = -one/su-0.25_dp*(one-two*a)**2/(a*(one-a)*su)
   if (func == RDMFT_CA) then
     f = su*sv+x; d1 = phi*sv+b; d11 = phip*sv; d12 = one
     if (sv > tiny_den) d12 = phi*(one-two*b)/(two*sv)+one
   else
     f = -su*sv+x; d1 = -phi*sv+b; d11 = -phip*sv; d12 = one
     if (sv > tiny_den) d12 = -phi*(one-two*b)/(two*sv)+one
   end if
 case (RDMFT_CGA)
   sp = sqrt(max(zero,a*(two-a))); sq = sqrt(max(zero,b*(two-b)))
   psi = half*(one-a)/sp
   psip = -half/sp-half*(one-a)**2/(a*(two-a)*sp)
   f = half*(x+sp*sq); d1 = half*b+psi*sq; d11 = psip*sq; d12 = half
   if (sq > tiny_den) d12 = half+psi*(one-b)/sq
 case (RDMFT_ML, RDMFT_MLSIC)
   if (func == RDMFT_ML) then
     dd = one+ml_b1*x; mm = ml_a0+two*ml_a1*x+ml_a1*ml_b1*x*x
     f = x*(ml_a0+ml_a1*x)/dd
     gp = mm/dd**2; gpp = two*(ml_a1*dd**2-ml_b1*mm)/dd**3
   else
     dd = one+mlsic_b1*x; mm = mlsic_a0+two*mlsic_a1*x+mlsic_a1*mlsic_b1*x*x
     f = x*(mlsic_a0+mlsic_a1*x)/dd
     gp = mm/dd**2; gpp = two*(mlsic_a1*dd**2-mlsic_b1*mm)/dd**3
   end if
   d1 = gp*b; d11 = gpp*b*b; d12 = gpp*x+gp
 case (RDMFT_POWER)
   f = x**alpha; d1 = alpha*a**(alpha-one)*b**alpha; d11 = alpha*(alpha-one)*a**(alpha-two)*b**alpha
   d12 = alpha*alpha*x**(alpha-one)
 case (RDMFT_MULLER_AS)
   f = half*(x-s); d1 = half*(b-half*sqrt(max(zero,b/a))); d11 = half*0.25_dp*sqrt(max(zero,b))*a**(-1.5_dp)
   d12 = half; if (s > tiny_den) d12 = half*(one-0.25_dp/s)
 case default
   ABI_ERROR('RDMFT: unknown JK-only functional.')
 end select

end subroutine fx_pair
!!***

!!****f* m_rdmft_functionals/rdmft_jk_pair
!! NAME
!!  rdmft_jk_pair
!!
!! FUNCTION
!!  For two different orbitals with occupations a and b and integrals jv, kv: the pair energy function
!!  F(a,b) = 2 f_H(a,b) J - f_X(a,b) K and its derivatives Fa=dF/da, Faa=d2F/da2, Fab=d2F/dadb.
!!
!! SOURCE
subroutine rdmft_jk_pair(func,a,b,jv,kv,alpha,ff,fa,faa,fab)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 real(dp),intent(in) :: a,b,jv,kv,alpha
 real(dp),intent(out) :: ff,fa,faa,fab

 !Local variables-------------------------------
 !scalars
 real(dp) :: fx,x1,x11,x12,fh,h1,h11,h12

! *************************************************************************

 call fx_pair(func,a,b,alpha,fx,x1,x11,x12)
 if (func == RDMFT_MULLER_AS) then
   fh = fx; h1 = x1; h11 = x11; h12 = x12
 else
   fh = a*b; h1 = b; h11 = zero; h12 = one
 end if
 ff = two*fh*jv-fx*kv
 fa = two*h1*jv-x1*kv
 faa = two*h11*jv-x11*kv
 fab = two*h12*jv-x12*kv

end subroutine rdmft_jk_pair
!!***

!!****f* m_rdmft_functionals/rdmft_jk_diag
!! NAME
!!  rdmft_jk_diag
!!
!! FUNCTION
!!  Same-orbital term F_d(n) = F(n,n) of the energy (with the special diagonal definitions of MLSIC and GU)
!!  and its first and second derivatives with respect to n, for an orbital with integrals jv=J_pp, kv=K_pp.
!!
!! SOURCE
subroutine rdmft_jk_diag(func,n,jv,kv,alpha,fd,fd1,fd2)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 real(dp),intent(in) :: n,jv,kv,alpha
 real(dp),intent(out) :: fd,fd1,fd2

 !Local variables-------------------------------
 !scalars
 real(dp) :: ff,fa,faa,fab

! *************************************************************************

 select case (func)
 case (RDMFT_MLSIC, RDMFT_GU)
   fd = n*n*(two*jv-kv); fd1 = two*n*(two*jv-kv); fd2 = two*(two*jv-kv)
 case default
   ! The general pair function at equal occupations 
   call rdmft_jk_pair(func,n,n,jv,kv,alpha,ff,fa,faa,fab)
   fd = ff; fd1 = two*fa; fd2 = two*faa+two*fab
 end select

end subroutine rdmft_jk_diag
!!***

!!****f* m_rdmft_functionals/rdmft_fermi_dirac
!! NAME
!!  rdmft_fermi_dirac
!!
!! FUNCTION
!!  Fermi-Dirac occupations n_p = 1/(exp((e_p-mu)/kt)+1) of orbitals with energies eig, with the chemical potential
!!  mu fixed by bisection so that the sum of the occupations is ntarget (occupations per spin, between 0 and 1).
!!
!! OUTPUT
!!  n(nf) = occupations
!!  ok = false if ntarget is not in (0,nf)
!!
!! SOURCE
subroutine rdmft_fermi_dirac(nf,eig,ntarget,kt,n,ok)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nf
 real(dp),intent(in) :: ntarget,kt
 logical,intent(out) :: ok
 !arrays
 real(dp),intent(in) :: eig(nf)
 real(dp),intent(out) :: n(nf)

 !Local variables-------------------------------
 !scalars
 integer :: it
 real(dp) :: mu,mu_hi,mu_lo,tot

! *************************************************************************

 ok = (ntarget > zero .and. ntarget < real(nf,dp) .and. kt > zero)
 n = zero
 if (.not.ok) return
 mu_lo = minval(eig)-50.0_dp*kt-one
 mu_hi = maxval(eig)+50.0_dp*kt+one
 do it=1,200
   mu = half*(mu_lo+mu_hi)
   n = one/(exp(max(-200.0_dp,min(200.0_dp,(eig-mu)/kt)))+one)
   tot = sum(n)
   if (abs(tot-ntarget) < 1.0e-13_dp*max(one,ntarget)) exit
   if (tot > ntarget) then
     mu_hi = mu
   else
     mu_lo = mu
   end if
 end do

end subroutine rdmft_fermi_dirac
!!***

!!****f* m_rdmft_functionals/rdmft_coef_jk
!! NAME
!!  rdmft_coef_jk
!!
!! FUNCTION
!!  Coefficients of the pair (j,q) in the energy written as
!!    E = sum_p nu_p h_pp + (1/N_k^2) sum_{pq} [ A_pq J_pq - C_pq K_pq ]
!!  for a JK-only functional: A = 2 f_H(n_j,n_q), C = f_X(n_j,n_q) (with the special definitions of the diagonal
!!  pair j=q). They are the input of the orbital gradient (generalized Fock matrix).
!!
!! SOURCE
subroutine rdmft_coef_jk(func,nj,nq,same,alpha,acoef,ccoef)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func
 real(dp),intent(in) :: nj,nq,alpha
 logical,intent(in) :: same
 real(dp),intent(out) :: acoef,ccoef

 !Local variables-------------------------------
 !scalars
 real(dp) :: fh,fx

! *************************************************************************

 fx = jk_exchange(func,nj,nq,same,alpha)
 fh = nj*nq
 if (func == RDMFT_MULLER_AS) fh = fx
 acoef = two*fh
 ccoef = fx

end subroutine rdmft_coef_jk
!!***

!!****f* m_rdmft_functionals/rdmft_coef_pnof
!! NAME
!!  rdmft_coef_pnof
!!
!! FUNCTION
!!  Coefficients A_ab, C_ab of the pair of orbitals (a,b) for a PNOF functional, as rdmft_coef_jk: the energy terms
!!  H_ab J_ab - X_ab K_ab + Pi_ab K_ab (a/=b) and n_a J_aa + selfint n_a^2 (J_aa-K_aa) (a=b) of rdmft_energy_pnof, so
!!  that A_ab = H_ab, C_ab = X_ab - Pi_ab, A_aa = n_a + selfint n_a^2, C_aa = selfint n_a^2.
!!  The background (frozen occupied) orbitals are trivial subspaces with n=1; an orbital that is not in a subspace
!!  (frozen virtual, n=0) has sub=0 and zero coefficients.
!!
!! INPUTS
!!  func = PNOF functional
!!  na, nb = occupations
!!  same_orb = true if a and b are the same orbital
!!  sub_a, sub_b = subspaces (any different integers for different subspaces)
!!  prin_a, prin_b = true for a principal pair
!!  hpa, hpb = 1 - occupation of the principal pair of the subspace of a, b (GNOF)
!!  selfint = rdmft_selfint
!!
!! SOURCE
subroutine rdmft_coef_pnof(func,na,nb,same_orb,sub_a,sub_b,prin_a,prin_b,hpa,hpb,selfint,acoef,ccoef)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,sub_a,sub_b,selfint
 real(dp),intent(in) :: na,nb,hpa,hpb
 logical,intent(in) :: same_orb,prin_a,prin_b
 real(dp),intent(out) :: acoef,ccoef

 !Local variables-------------------------------
 !scalars
 real(dp) :: cross,hha,hhb,nda,ndb,pi_ab,xab
 logical :: same_sub

! *************************************************************************

 if (same_orb) then
   acoef = na+real(selfint,dp)*na**2
   ccoef = real(selfint,dp)*na**2
   return
 end if
 same_sub = (sub_a == sub_b .and. sub_a /= 0)
 xab = na*nb
 if (same_sub) xab = zero
 hha = one-na
 hhb = one-nb
 pi_ab = zero
 if (same_sub) then
   pi_ab = sqrt(max(zero,na*nb))
   if (prin_a .or. prin_b) pi_ab = -pi_ab
 else
   select case (func)
   case (RDMFT_PNOF5)
     pi_ab = zero
   case (RDMFT_PNOF7)
     pi_ab = -sqrt(max(zero,na*nb*hha*hhb))
   case (RDMFT_PNOF7S)
     pi_ab = -four*na*nb*hha*hhb
   case (RDMFT_GNOF)
     if (.not.(prin_a .and. prin_b)) then
       nda = na*exp(-(hpa/gnof_c)**2)
       ndb = nb*exp(-(hpb/gnof_c)**2)
       cross = sqrt(max(zero,nda*ndb))
       pi_ab = nda*ndb-sqrt(max(zero,na*nb*hha*hhb))
       if (prin_a .neqv. prin_b) then
         pi_ab = pi_ab-cross
       else
         pi_ab = pi_ab+cross
       end if
     end if
   case default
     ABI_ERROR('RDMFT: unknown PNOF functional.')
   end select
 end if
 acoef = two*xab
 ccoef = xab-pi_ab

end subroutine rdmft_coef_pnof
!!***

!!****f* m_rdmft_functionals/rdmft_energy_bg
!! NAME
!!  rdmft_energy_bg
!!
!! FUNCTION
!!  Energy per unit cell of the background orbitals (frozen occupied, n=1): their one-electron energy and the
!!  Coulomb and exchange energy among them. For the PNOF functionals (func>=RDMFT_PNOF5) these orbitals are trivial
!!  subspaces; selfint is rdmft_selfint.
!!
!! SOURCE
function rdmft_energy_bg(func,nk,bg,alpha,selfint) result(energy)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,selfint
 real(dp),intent(in) :: alpha
 type(rdmft_bg_t),intent(in) :: bg
 real(dp) :: energy

 !Local variables-------------------------------
 !scalars
 real(dp) :: faa,fab,fa,fdj,fdk,fdj1,fdj2,fdk1,fdk2,ffj,ffk

! *************************************************************************

 if (func >= RDMFT_PNOF5) then
   energy = two*bg%sh/nk+(two*bg%sjo-bg%sko+(one+real(selfint,dp))*bg%sjd-real(selfint,dp)*bg%skd)/real(nk,dp)**2
 else
   call rdmft_jk_pair(func,one,one,one,zero,alpha,ffj,fa,faa,fab)
   call rdmft_jk_pair(func,one,one,zero,one,alpha,ffk,fa,faa,fab)
   call rdmft_jk_diag(func,one,one,zero,alpha,fdj,fdj1,fdj2)
   call rdmft_jk_diag(func,one,zero,one,alpha,fdk,fdk1,fdk2)
   energy = two*bg%sh/nk+(bg%sjo*ffj+bg%sko*ffk+bg%sjd*fdj+bg%skd*fdk)/real(nk,dp)**2
 end if

end function rdmft_energy_bg
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
!!  alpha = exponent of the POWER functional
!!  bj(nbk), bk(nbk) = sums of J_pq and K_pq over the background orbitals q
!!  bg <type(rdmft_bg_t)> = constants of the background orbitals
!!
!! SOURCE
function rdmft_energy_jk(func,nk,nbk,n,h1,jint,kint,alpha,bj,bk,bg) result(energy)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk
 real(dp),intent(in) :: alpha
 type(rdmft_bg_t),intent(in) :: bg
 real(dp) :: energy
 !arrays
 real(dp),intent(in) :: n(nbk),h1(nbk),jint(nbk,nbk),kint(nbk,nbk),bj(nbk),bk(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ip,iq
 real(dp) :: e1,e2,fa,faa,fab,fh,ff,fx

! *************************************************************************

 e1 = zero
 e2 = zero
 do ip=1,nbk
   e1 = e1+two*n(ip)*h1(ip)
   do iq=1,nbk
     fx = jk_exchange(func,n(ip),n(iq),ip==iq,alpha)
     fh = n(ip)*n(iq)
     if (func == RDMFT_MULLER_AS) fh = fx
     e2 = e2+two*fh*jint(ip,iq)-fx*kint(ip,iq)
   end do
   ! Interaction with the background orbitals (occupation 1), counted for the two orders of the pairs
   call rdmft_jk_pair(func,n(ip),one,bj(ip),bk(ip),alpha,ff,fa,faa,fab)
   e2 = e2+two*ff
 end do
 energy = e1/nk+e2/real(nk,dp)**2+rdmft_energy_bg(func,nk,bg,alpha,0)

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
 ABI_SFREE(table%class_var)
 ABI_SFREE(table%class_nvirt)
 ABI_SFREE(table%var_mult)
 ABI_SFREE(table%orb_var)
 table%nclass = 0
 table%nvar = 0
 table%ngam = 0
 table%ng = 0
 table%nsub = 0

end subroutine rdmft_pnof_table_free
!!***

!!****f* m_rdmft_functionals/rdmft_active_space
!! NAME
!!  rdmft_active_space
!!
!! FUNCTION
!!  Select the active orbitals, whose occupations are optimized. The other occupied orbitals are frozen at 1 and
!!  the other virtual orbitals at 0.
!!  Occupied orbitals: if nocc_act>0, the nocc_act highest occupied bands of each k-point; else if ecut_act>0,
!!  those with an energy larger than E_HOMO-ecut_act; else all of them.
!!  Virtual orbitals: if nvir_act>0, the nvir_act lowest virtual bands of each k-point; else if ecut_act>0,
!!  those with an energy smaller than E_LUMO+ecut_act; else all of them.
!!  Orbitals of a degenerate shell on the edge of an energy window are all kept (tolerance 1e-8 Ha).
!!
!! INPUTS
!!  nbk = number of orbitals (composite index over n and k)
!!  eig(nbk) = orbital energies
!!  band(nbk) = band index of each orbital
!!  nocc_band(nbk) = number of occupied bands of the k-point of each orbital
!!  occupied(nbk) = true for the occupied orbitals
!!  nocc_act, nvir_act, ecut_act = rdmft_nocc_act, rdmft_nvir_act, rdmft_ecut_act
!!
!! OUTPUT
!!  active(nbk) = true for the active orbitals
!!
!! SOURCE
subroutine rdmft_active_space(nbk,eig,band,nocc_band,occupied,nocc_act,nvir_act,ecut_act,active)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nbk,nocc_act,nvir_act
 real(dp),intent(in) :: ecut_act
 !arrays
 integer,intent(in) :: band(nbk),nocc_band(nbk)
 real(dp),intent(in) :: eig(nbk)
 logical,intent(in) :: occupied(nbk)
 logical,intent(out) :: active(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ip
 real(dp),parameter :: tol_edge = 1.0e-8_dp
 real(dp) :: ehomo,elumo

! *************************************************************************

 ehomo = -huge(one)
 elumo = huge(one)
 do ip=1,nbk
   if (occupied(ip)) then
     ehomo = max(ehomo,eig(ip))
   else
     elumo = min(elumo,eig(ip))
   end if
 end do

 do ip=1,nbk
   if (occupied(ip)) then
     if (nocc_act > 0) then
       active(ip) = (band(ip) > nocc_band(ip)-nocc_act)
     else if (ecut_act > zero) then
       active(ip) = (eig(ip) >= ehomo-ecut_act-tol_edge)
     else
       active(ip) = .true.
     end if
   else
     if (nvir_act > 0) then
       active(ip) = (band(ip) <= nocc_band(ip)+nvir_act)
     else if (ecut_act > zero) then
       active(ip) = (eig(ip) <= elumo+ecut_act+tol_edge)
     else
       active(ip) = .true.
     end if
   end if
 end do

end subroutine rdmft_active_space
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_table_build
!! NAME
!!  rdmft_pnof_table_build
!!
!! FUNCTION
!!  Build the pairs and coupling subspaces from the orbital energies, over all k-points, closed under k -> -k.
!!  The time-reversal partner of an orbital is partner(p) (p itself for the time-reversal-invariant k-points);
!!  a "star" is {p, partner(p)}. Every active occupied star is a principal star, from the highest energy to the lowest,
!!  coupled with the (ncoupl-1) lowest unused virtual stars (over all k-points):
!!   - principal orbital p with partner(p)=p: one subspace with p and all the orbitals of its virtual stars (a
!!     generic virtual star has two orbitals with the same occupation);
!!   - principal orbital p with partner(p)/=p: two mirror subspaces, {p, V_j} and {partner(p), partner(V_j)}, with the
!!     same occupations; only generic virtual stars (two orbitals) can be used since a time-reversal-invariant
!!     virtual orbital would belong to both.
!!  The other occupied orbitals are frozen at 1 (background orbitals, not in the table).
!!
!! INPUTS
!!  nbk = number of orbitals
!!  eig(nbk) = orbital energies
!!  occupied(nbk) = true for the occupied orbitals
!!  active(nbk) = true for the active orbitals (see rdmft_active_space), only the occupied ones are used here
!!  partner(nbk) = time-reversal partner of each orbital
!!  ncoupl = number of stars of each class (1 principal + ncoupl-1 virtual)
!!
!! OUTPUT
!!  ierr = 0 if ok, 1 if there are no principal pair or not enough virtual orbitals,
!!         2 if the active window is not closed under time reversal
!!
!! SOURCE
subroutine rdmft_pnof_table_build(table,nbk,eig,occupied,active,partner,ncoupl,ierr)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_pnof_table_t),intent(inout) :: table
 integer,intent(in) :: nbk,ncoupl
 integer,intent(out) :: ierr
 !arrays
 integer,intent(in) :: partner(nbk)
 real(dp),intent(in) :: eig(nbk)
 logical,intent(in) :: occupied(nbk),active(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ic,iv,ip,jp,m,nc,ng,nsub,nvar,sub1,sub2
 logical :: trim_p
 !arrays
 integer,allocatable :: cp(:),cv(:,:)
 logical,allocatable :: used(:)

! *************************************************************************

 call table%free()
 ierr = 0
 m = ncoupl-1
 if (m < 1 .or. count(occupied .and. active) < 1) then
   ierr = 1
   return
 end if
 ABI_MALLOC(used,(nbk))
 ABI_MALLOC(cp,(nbk))
 ABI_MALLOC(cv,(m,nbk))
 used = .false.
 nc = 0

 ! Principal stars (highest energy first), each with its virtual stars
 do
   jp = 0
   do ip=1,nbk
     if (occupied(ip) .and. active(ip) .and. .not.used(ip)) then
       if (jp == 0) then
         jp = ip
       else if (eig(ip) > eig(jp)) then
         jp = ip
       end if
     end if
   end do
   if (jp == 0) exit
   if (.not.active(partner(jp))) then
     ierr = 2
     exit
   end if
   trim_p = (partner(jp) == jp)
   used(jp) = .true.
   used(partner(jp)) = .true.
   nc = nc+1
   cp(nc) = jp
   do iv=1,m
     ip = 0
     do jp=1,nbk
       if (occupied(jp) .or. used(jp)) cycle
       if (.not.trim_p .and. partner(jp) == jp) cycle
       if (ip == 0) then
         ip = jp
       else if (eig(jp) < eig(ip)) then
         ip = jp
       end if
     end do
     if (ip == 0) then
       ierr = 1
       exit
     end if
     used(ip) = .true.
     used(partner(ip)) = .true.
     cv(iv,nc) = ip
   end do
   if (ierr /= 0) exit
 end do
 if (ierr /= 0) then
   ABI_FREE(used)
   ABI_FREE(cp)
   ABI_FREE(cv)
   return
 end if

 ! Tables of pairs and variables
 ABI_MALLOC(table%gem_orb,(nbk))
 ABI_MALLOC(table%gem_sub,(nbk))
 ABI_MALLOC(table%gem_principal,(nbk))
 ABI_MALLOC(table%orb_var,(nbk))
 ABI_MALLOC(table%class_var,(nc))
 ABI_MALLOC(table%class_nvirt,(nc))
 ABI_MALLOC(table%var_mult,(nc*(m+1)))
 table%orb_var = 0
 table%nclass = nc
 ng = 0
 nsub = 0
 nvar = 0
 do ic=1,nc
   jp = cp(ic)
   table%class_var(ic) = nvar+1
   table%class_nvirt(ic) = m
   ! Principal variable
   nvar = nvar+1
   table%var_mult(nvar) = 1
   table%orb_var(jp) = nvar
   table%orb_var(partner(jp)) = nvar
   nsub = nsub+1
   sub1 = nsub
   ng = ng+1
   table%gem_orb(ng) = jp; table%gem_sub(ng) = sub1; table%gem_principal(ng) = .true.
   sub2 = sub1
   if (partner(jp) /= jp) then
     nsub = nsub+1
     sub2 = nsub
     ng = ng+1
     table%gem_orb(ng) = partner(jp); table%gem_sub(ng) = sub2; table%gem_principal(ng) = .true.
   end if
   ! Virtual variables
   do iv=1,m
     ip = cv(iv,ic)
     nvar = nvar+1
     table%orb_var(ip) = nvar
     table%orb_var(partner(ip)) = nvar
     ng = ng+1
     table%gem_orb(ng) = ip; table%gem_sub(ng) = sub1; table%gem_principal(ng) = .false.
     if (partner(jp) == jp) then
       ! Time-reversal-invariant principal pair: the partner of a generic virtual orbital is in the same subspace
       table%var_mult(nvar) = 1
       if (partner(ip) /= ip) then
         table%var_mult(nvar) = 2
         ng = ng+1
         table%gem_orb(ng) = partner(ip); table%gem_sub(ng) = sub1; table%gem_principal(ng) = .false.
       end if
     else
       ! Generic principal pair: mirror subspace
       table%var_mult(nvar) = 1
       ng = ng+1
       table%gem_orb(ng) = partner(ip); table%gem_sub(ng) = sub2; table%gem_principal(ng) = .false.
     end if
   end do
 end do
 table%nvar = nvar
 table%ngam = nc*m
 ! The other occupied orbitals are frozen at 1: they are background orbitals and are not in the table
 table%ng = ng
 table%nsub = nsub

 ABI_FREE(used)
 ABI_FREE(cp)
 ABI_FREE(cv)

end subroutine rdmft_pnof_table_build
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_table_build_sym
!! NAME
!!  rdmft_pnof_table_build_sym
!!
!! FUNCTION
!!  Table of pairs and subspaces for the hybrid symmetry reduction. The units of occupation are the stars of the orbital
!!  optimization (star_rep): for the generic k-points (trivial little group) the orbit of the orbital (n,k) under the
!!  symmetry operations of the ground state, for the other k-points the time-reversal pair {k,-k} (or k alone if it is its
!!  own -k), as in rdmft_pnof_table_build. All the orbitals of a unit have the same occupation (one variable).
!!  Every active occupied unit is a principal unit, from the highest energy to the lowest, coupled with the (ncoupl-1)
!!  lowest unused virtual units for which an assignment exists:
!!   - generic principal unit: the subspaces form a family invariant under the operations, S Sigma(p) = Sigma(S p). The
!!     virtual unit V must be generic and have a member v whose stabilizer (operations, including antiunitary ones, that
!!     leave its k-point unchanged) is contained in the one of the principal p; the subspace of p contains the orbit of v
!!     under the stabilizer of p and the one of the image S p its image, so that the subspaces of the principal unit use
!!     every orbital of V exactly once;
!!   - principal unit {p} alone (time-reversal-invariant k-point): one subspace with all the orbitals of V;
!!   - principal pair {p,partner(p)}: two mirror subspaces, with one orbital of each pair {v,partner(v)} of V and its
!!     partner (V without time-reversal-invariant k-points).
!!  The operations are given by their action on the k-points, kimg(k,op) (0 for all k if the operation does not map the
!!  mesh onto itself). The degenerate shells at the non-generic k-points need no special treatment: their orbitals are
!!  independent and free to mix.
!!
!! INPUTS
!!  nbk, nkpt = numbers of orbitals and of k-points; offset(nkpt), nband(nkpt) = composite index of the orbitals
!!  eig, occupied, active = as rdmft_pnof_table_build
!!  ncoupl = number of units of each class (1 principal + ncoupl-1 virtual)
!!  star_rep(nkpt) = representative of the star of each k-point
!!  generic(nkpt) = generic k-points; trk(nkpt) = index of -k (0 if not in the mesh)
!!  kimg(nkpt,nops) = image of each k-point by each operation (unitary and antiunitary)
!!
!! OUTPUT
!!  ierr = 0 if ok, 1 no principal unit or not enough eligible virtual units, 3 the active window or the occupied set is
!!         not the same inside a unit
!!
!! SOURCE
subroutine rdmft_pnof_table_build_sym(table,nbk,nkpt,offset,nband,eig,occupied,active,ncoupl,star_rep,generic,trk,kimg,ierr)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_pnof_table_t),intent(inout) :: table
 integer,intent(in) :: nbk,nkpt,ncoupl
 integer,intent(out) :: ierr
 !arrays
 integer,intent(in) :: offset(nkpt),nband(nkpt),star_rep(nkpt),trk(nkpt),kimg(:,:)
 logical,intent(in) :: generic(nkpt)
 real(dp),intent(in) :: eig(nbk)
 logical,intent(in) :: occupied(nbk),active(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: i,ig,iop,iv,ip,jk,jp,kp,m,maxnc,n,nc,ng,nops,nsub,nvar,nw,s_p,sub0,vv,w
 logical :: found
 !arrays
 integer,allocatable :: orb_k(:),orb_n(:),members(:),w0(:),wk(:,:),cvar(:),cnv(:),vmult(:),gorb(:),gsub(:),ovar(:)
 logical,allocatable :: used(:),skipv(:),usable(:),fixk(:,:),gprin(:),tm(:)

! *************************************************************************

 call table%free()
 ierr = 0
 m = ncoupl-1
 maxnc = count(occupied .and. active)
 if (m < 1 .or. maxnc < 1) then
   ierr = 1
   return
 end if
 nops = size(kimg,2)
 ABI_MALLOC(usable,(nops))
 ABI_MALLOC(fixk,(nops,nkpt))
 do iop=1,nops
   usable(iop) = (kimg(1,iop) /= 0)
   do jk=1,nkpt
     fixk(iop,jk) = (usable(iop) .and. kimg(jk,iop) == jk)
   end do
 end do
 ABI_MALLOC(orb_k,(nbk))
 ABI_MALLOC(orb_n,(nbk))
 do jk=1,nkpt
   do n=1,nband(jk)
     orb_k(offset(jk)+n) = jk
     orb_n(offset(jk)+n) = n
   end do
 end do
 ABI_MALLOC(used,(nbk))
 ABI_MALLOC(skipv,(nbk))
 ABI_MALLOC(members,(nkpt))
 ABI_MALLOC(w0,(nkpt))
 ABI_MALLOC(wk,(nkpt,nkpt))
 ABI_MALLOC(tm,(nkpt))
 ABI_MALLOC(cvar,(maxnc))
 ABI_MALLOC(cnv,(maxnc))
 ABI_MALLOC(vmult,(maxnc*(m+1)))
 ABI_MALLOC(gorb,(nbk))
 ABI_MALLOC(gsub,(nbk))
 ABI_MALLOC(gprin,(nbk))
 ABI_MALLOC(ovar,(nbk))
 used = .false.
 ovar = 0
 nc = 0
 ng = 0
 nsub = 0
 nvar = 0

 ! The occupied and active flags are the same inside a unit
 do ip=1,nbk
   do jk=1,nkpt
     if (star_rep(jk) /= star_rep(orb_k(ip))) cycle
     if (orb_n(ip) > nband(jk)) then
       ierr = 3
       exit
     end if
     jp = offset(jk)+orb_n(ip)
     if ((occupied(jp) .neqv. occupied(ip)) .or. (active(jp) .neqv. active(ip))) ierr = 3
   end do
   if (ierr /= 0) exit
 end do

 ! Principal units (highest energy first), each with its virtual units
 do while (ierr == 0)
   jp = 0
   do ip=1,nbk
     if (occupied(ip) .and. active(ip) .and. .not.used(ip)) then
       if (jp == 0) then
         jp = ip
       else if (eig(ip) > eig(jp)) then
         jp = ip
       end if
     end if
   end do
   if (jp == 0) exit
   kp = orb_k(jp)
   n = orb_n(jp)
   nc = nc+1
   cvar(nc) = nvar+1
   cnv(nc) = m
   s_p = 1
   members(1) = kp
   do jk=1,nkpt
     if (jk /= kp .and. star_rep(jk) == star_rep(kp)) then
       s_p = s_p+1
       members(s_p) = jk
     end if
   end do
   ! Principal variable and subspaces: one subspace per member of the principal unit
   nvar = nvar+1
   vmult(nvar) = 1
   sub0 = nsub
   do i=1,s_p
     used(offset(members(i))+n) = .true.
     ovar(offset(members(i))+n) = nvar
     nsub = nsub+1
     ng = ng+1
     gorb(ng) = offset(members(i))+n
     gsub(ng) = nsub
     gprin(ng) = .true.
   end do
   skipv = .false.
   do iv=1,m
     ! lowest-energy unused virtual unit with an assignment
     found = .false.
     do
       ip = 0
       do vv=1,nbk
         if (occupied(vv) .or. used(vv) .or. skipv(vv)) cycle
         if (ip == 0) then
           ip = vv
         else if (eig(vv) < eig(ip)) then
           ip = vv
         end if
       end do
       if (ip == 0) exit
       call try_virtual(ip,found)
       if (found) exit
       do jk=1,nkpt
         if (star_rep(jk) == star_rep(orb_k(ip)) .and. orb_n(ip) <= nband(jk)) skipv(offset(jk)+orb_n(ip)) = .true.
       end do
     end do
     if (.not.found) then
       ierr = 1
       exit
     end if
     ! Accept: variable of the unit, orbitals of each subspace
     nvar = nvar+1
     vmult(nvar) = nw
     do i=1,s_p
       do w=1,nw
         ig = offset(wk(w,i))+orb_n(ip)
         used(ig) = .true.
         ovar(ig) = nvar
         ng = ng+1
         gorb(ng) = ig
         gsub(ng) = sub0+i
         gprin(ng) = .false.
       end do
     end do
   end do
   if (ierr /= 0) exit
 end do
 if (ierr == 0 .and. nc == 0) ierr = 1

 if (ierr == 0) then
   table%nclass = nc
   table%nvar = nvar
   table%ngam = nc*m
   table%ng = ng
   table%nsub = nsub
   ABI_MALLOC(table%gem_orb,(nbk))
   ABI_MALLOC(table%gem_sub,(nbk))
   ABI_MALLOC(table%gem_principal,(nbk))
   ABI_MALLOC(table%orb_var,(nbk))
   ABI_MALLOC(table%class_var,(nc))
   ABI_MALLOC(table%class_nvirt,(nc))
   ABI_MALLOC(table%var_mult,(nvar))
   table%gem_orb(1:ng) = gorb(1:ng)
   table%gem_sub(1:ng) = gsub(1:ng)
   table%gem_principal(1:ng) = gprin(1:ng)
   table%orb_var = ovar
   table%class_var = cvar(1:nc)
   table%class_nvirt = cnv(1:nc)
   table%var_mult = vmult(1:nvar)
 end if

 ABI_FREE(usable)
 ABI_FREE(fixk)
 ABI_FREE(orb_k)
 ABI_FREE(orb_n)
 ABI_FREE(used)
 ABI_FREE(skipv)
 ABI_FREE(members)
 ABI_FREE(w0)
 ABI_FREE(wk)
 ABI_FREE(tm)
 ABI_FREE(cvar)
 ABI_FREE(cnv)
 ABI_FREE(vmult)
 ABI_FREE(gorb)
 ABI_FREE(gsub)
 ABI_FREE(gprin)
 ABI_FREE(ovar)

contains

 subroutine try_virtual(ipv,ok_v)
   ! Assignment of the virtual unit of the orbital ipv to the subspaces of the principal unit (kp, members, s_p)
   integer,intent(in) :: ipv
   logical,intent(out) :: ok_v
   integer :: iq,kv,kvm,nv,s_v,mk,ii,iw,kk
   logical :: sub
   ok_v = .false.
   kv = orb_k(ipv)
   nv = orb_n(ipv)
   s_v = 0
   do mk=1,nkpt
     if (star_rep(mk) == star_rep(kv) .and. nv <= nband(mk)) s_v = s_v+1
   end do
   if (generic(kp)) then
     ! invariant family of subspaces: the virtual unit must be generic too
     if (.not.generic(kv)) return
     kvm = 0
     do mk=1,nkpt
       if (star_rep(mk) /= star_rep(kv) .or. nv > nband(mk)) cycle
       sub = .true.
       do iq=1,nops
         if (fixk(iq,mk) .and. .not.fixk(iq,kp)) sub = .false.
       end do
       if (sub) then
         kvm = mk
         exit
       end if
     end do
     if (kvm == 0) return
     ! orbit of the member under the stabilizer of the principal
     tm = .false.
     nw = 0
     do iq=1,nops
       if (.not.fixk(iq,kp)) cycle
       kk = kimg(kvm,iq)
       if (.not.tm(kk)) then
         tm(kk) = .true.
         nw = nw+1
         w0(nw) = kk
       end if
     end do
     if (s_p*nw /= s_v) return
     ! images under an operation that maps kp onto each member of the principal unit
     do ii=1,s_p
       if (ii == 1) then
         wk(1:nw,1) = w0(1:nw)
       else
         iq = 0
         do iw=1,nops
           if (usable(iw) .and. kimg(kp,iw) == members(ii)) then
             iq = iw
             exit
           end if
         end do
         if (iq == 0) return
         do iw=1,nw
           wk(iw,ii) = kimg(w0(iw),iq)
         end do
       end if
     end do
   else
     ! non-generic principal: time-reversal type table
     nw = 0
     do mk=1,nkpt
       if (star_rep(mk) /= star_rep(kv) .or. nv > nband(mk)) cycle
       if (s_p == 1) then
         nw = nw+1
         wk(nw,1) = mk
       else
         if (trk(mk) == mk .or. trk(mk) == 0) return
         if (mk < trk(mk)) then
           nw = nw+1
           wk(nw,1) = mk
           wk(nw,2) = trk(mk)
         end if
       end if
     end do
     if (s_p == 2 .and. 2*nw /= s_v) return
     if (s_p > 2) return
   end if
   ! the orbitals must be distinct and unused
   tm = .false.
   do ii=1,s_p
     do iw=1,nw
       kk = wk(iw,ii)
       if (tm(kk) .or. used(offset(kk)+nv)) return
       tm(kk) = .true.
     end do
   end do
   ok_v = .true.
 end subroutine try_virtual

end subroutine rdmft_pnof_table_build_sym
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_table_remap
!! NAME
!!  rdmft_pnof_table_remap
!!
!! FUNCTION
!!  Change the orbital indices of a table from the composite index (all the orbitals) to the index among the free
!!  orbitals: free_of(composite) = index (0 if not free).
!!
!! SOURCE
subroutine rdmft_pnof_table_remap(table,nbk,free_of,nfree)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_pnof_table_t),intent(inout) :: table
 integer,intent(in) :: nbk,nfree
 !arrays
 integer,intent(in) :: free_of(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ig,ip
 !arrays
 integer,allocatable :: ov(:)

! *************************************************************************

 do ig=1,table%ng
   if (free_of(table%gem_orb(ig)) == 0) then
     ABI_BUG('RDMFT: an orbital of the PNOF table is not in the free orbitals.')
   end if
   table%gem_orb(ig) = free_of(table%gem_orb(ig))
 end do
 ABI_MALLOC(ov,(nfree))
 ov = 0
 do ip=1,nbk
   if (free_of(ip) > 0) ov(free_of(ip)) = table%orb_var(ip)
 end do
 call move_alloc(ov,table%orb_var)

end subroutine rdmft_pnof_table_remap
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_gamma
!! NAME
!!  rdmft_pnof_gamma
!!
!! FUNCTION
!!  Angles gamma giving (as closely as the bounds allow) the occupations n(orbital) of the variables of the
!!  classes: inverse of rdmft_pnof_occ, gamma_0=acos(sqrt(2 n_0-1)), gamma_j=asin(sqrt(y_j/R_{j-1})).
!!
!! SOURCE
subroutine rdmft_pnof_gamma(table,nfree,n,gam)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_pnof_table_t),intent(in) :: table
 integer,intent(in) :: nfree
 !arrays
 real(dp),intent(in) :: n(nfree)
 real(dp),intent(out) :: gam(table%ngam)

 !Local variables-------------------------------
 !scalars
 integer :: ic,ip,j,m,v0,g0
 real(dp),parameter :: pi4=0.78539816339744830962_dp
 real(dp) :: rr
 !arrays
 real(dp),allocatable :: nv(:)

! *************************************************************************

 ABI_MALLOC(nv,(table%nvar))
 nv = zero
 do ip=1,nfree
   if (table%orb_var(ip) > 0) nv(table%orb_var(ip)) = n(ip)
 end do
 gam = pi4
 g0 = 0
 do ic=1,table%nclass
   m = table%class_nvirt(ic)
   v0 = table%class_var(ic)
   gam(g0+1) = acos(sqrt(min(one,max(zero,two*nv(v0)-one))))
   rr = one-nv(v0)
   do j=1,m-1
     if (rr > 1.0e-14_dp) gam(g0+1+j) = asin(sqrt(min(one,max(zero,nv(v0+j)*real(table%var_mult(v0+j),dp)/rr))))
     rr = rr-nv(v0+j)*real(table%var_mult(v0+j),dp)
   end do
   g0 = g0+m
 end do
 ABI_FREE(nv)

end subroutine rdmft_pnof_gamma
!!***

!!****f* m_rdmft_functionals/rdmft_pnof_occ
!! NAME
!!  rdmft_pnof_occ
!!
!! FUNCTION
!!  Occupation numbers of the variables of the classes from the angles gamma, and the Jacobian dn_var/dgamma.
!!  For a class with principal variable n_0 and m virtual variables with multiplicities mu_j:
!!    n_0 = 1/2 + cos^2(gamma_0)/2,  y_j = R_{j-1} sin^2(gamma_j) (j<m),  y_m = R_{m-1},  R_j = R_{j-1} cos^2(gamma_j),
!!    R_0 = 1 - n_0,  n_j = y_j / mu_j
!!  so that n_0 + sum_j mu_j n_j = 1 for any angles and all the occupations are in (0,1).
!!
!! INPUTS
!!  table <type(rdmft_pnof_table_t)> = classes and variables
!!  gam(ngam) = angles
!!
!! OUTPUT
!!  nvar_occ(nvar) = occupations of the variables
!!  jac(nvar,ngam) = derivatives of the occupations with respect to the angles
!!
!! SOURCE
subroutine rdmft_pnof_occ(table,gam,nvar_occ,jac)

 !Arguments ------------------------------------
 !scalars
 type(rdmft_pnof_table_t),intent(in) :: table
 !arrays
 real(dp),intent(in) :: gam(table%ngam)
 real(dp),intent(out) :: nvar_occ(table%nvar),jac(table%nvar,table%ngam)

 !Local variables-------------------------------
 !scalars
 integer :: ic,i,j,m,v0,g0
 real(dp) :: c,s
 !arrays
 real(dp),allocatable :: dr(:),dy(:),rr(:)

! *************************************************************************

 nvar_occ = zero
 jac = zero
 g0 = 0
 do ic=1,table%nclass
   m = table%class_nvirt(ic)
   v0 = table%class_var(ic)
   ABI_MALLOC(rr,(0:m))
   ABI_MALLOC(dr,(m+1))   ! dR_{j-1}/dgamma_i for the current j
   ABI_MALLOC(dy,(m+1))
   ! principal: n_0 = 1 - R_0, R_0 = sin^2(gamma_0)/2
   s = sin(gam(g0+1)); c = cos(gam(g0+1))
   rr(0) = half*s*s
   nvar_occ(v0) = one-rr(0)
   jac(v0,g0+1) = -s*c
   dr = zero
   dr(1) = s*c                         ! dR_0/dgamma_0
   do j=1,m
     if (j < m) then
       s = sin(gam(g0+1+j)); c = cos(gam(g0+1+j))
       ! y_j = R_{j-1} sin^2(gamma_j); index i refers to gamma_{i-1}
       dy = zero
       do i=1,j
         dy(i) = dr(i)*s*s
       end do
       dy(j+1) = rr(j-1)*two*s*c
       nvar_occ(v0+j) = rr(j-1)*s*s/real(table%var_mult(v0+j),dp)
       do i=1,j+1
         jac(v0+j,g0+i) = dy(i)/real(table%var_mult(v0+j),dp)
       end do
       ! R_j = R_{j-1} cos^2(gamma_j)
       rr(j) = rr(j-1)*c*c
       do i=1,j
         dr(i) = dr(i)*c*c
       end do
       dr(j+1) = -rr(j-1)*two*s*c
     else
       ! last virtual variable takes the rest: y_m = R_{m-1}
       nvar_occ(v0+j) = rr(j-1)/real(table%var_mult(v0+j),dp)
       do i=1,j
         jac(v0+j,g0+i) = dr(i)/real(table%var_mult(v0+j),dp)
       end do
     end if
   end do
   g0 = g0+m
   ABI_FREE(rr)
   ABI_FREE(dr)
   ABI_FREE(dy)
 end do

end subroutine rdmft_pnof_occ
!!***

!!****f* m_rdmft_functionals/rdmft_grad_pnof
!! NAME
!!  rdmft_grad_pnof
!!
!! FUNCTION
!!  Gradient of the PNOF electronic energy with respect to the occupation number of each orbital of the table,
!!  treating the occupations of all the orbitals as independent (the time-reversal partners are tied by the caller).
!!  For GNOF the dependence of the pairing term on the occupation of the principal pair of the subspace is included.
!!
!! INPUTS
!!  as rdmft_energy_pnof
!!
!! OUTPUT
!!  grad(nbk) = dE/dn_p (zero for the orbitals that are not in the table)
!!
!! SOURCE
subroutine rdmft_grad_pnof(func,table,nk,nbk,n,h1,jint,kint,selfint,bj,bk,grad)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk,selfint
 type(rdmft_pnof_table_t),intent(in) :: table
 !arrays
 real(dp),intent(in) :: n(nbk),h1(nbk),jint(nbk,nbk),kint(nbk,nbk),bj(nbk),bk(nbk)
 real(dp),intent(out) :: grad(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ia,ib,ig,ip,ng,sa
 real(dp),parameter :: tiny_den = 1.0e-300_dp
 real(dp) :: dpi,dpr,dnd,dsa,ea,f,hpa,na,nda,nk2,sigma,viap,sga
 logical :: pra,same_sub
 !arrays
 integer,allocatable :: iorb(:),prin_orb(:),sub(:)
 logical,allocatable :: pr(:)
 real(dp),allocatable :: gp(:),gt(:),hprin(:),jcol(:),kcol(:),nd(:),nt(:),sd(:),sn(:),st(:)

! *************************************************************************

 if (func /= RDMFT_PNOF5 .and. func /= RDMFT_PNOF7 .and. func /= RDMFT_PNOF7S .and. func /= RDMFT_GNOF) then
   ABI_ERROR('RDMFT: unknown PNOF functional.')
 end if
 grad = zero
 nk2 = real(nk,dp)**2
 ng = table%ng
 ABI_MALLOC(prin_orb,(table%nsub))
 ABI_MALLOC(hprin,(table%nsub))
 ABI_MALLOC(gp,(table%nsub))
 prin_orb = 0
 hprin = one
 gp = zero
 do ig=1,ng
   if (table%gem_principal(ig)) then
     prin_orb(table%gem_sub(ig)) = table%gem_orb(ig)
     hprin(table%gem_sub(ig)) = one-n(table%gem_orb(ig))
   end if
 end do

 ! Per-orbital quantities of the table (the square roots and exponentials are computed once per orbital, not per pair)
 ABI_MALLOC(iorb,(ng))
 ABI_MALLOC(sub,(ng))
 ABI_MALLOC(pr,(ng))
 ABI_MALLOC(nt,(ng))
 ABI_MALLOC(gt,(ng))
 ABI_MALLOC(st,(ng))
 ABI_MALLOC(sn,(ng))
 ABI_MALLOC(nd,(ng))
 ABI_MALLOC(sd,(ng))
 ABI_MALLOC(jcol,(ng))
 ABI_MALLOC(kcol,(ng))
 do ig=1,ng
   iorb(ig) = table%gem_orb(ig)
   sub(ig) = table%gem_sub(ig)
   pr(ig) = table%gem_principal(ig)
   nt(ig) = n(iorb(ig))
   gt(ig) = nt(ig)*(one-nt(ig))
   st(ig) = sqrt(max(zero,gt(ig)))
   sn(ig) = sqrt(max(zero,nt(ig)))
   nd(ig) = nt(ig)*exp(-(hprin(sub(ig))/gnof_c)**2)
   sd(ig) = sqrt(max(zero,nd(ig)))
 end do

 do ia=1,ng
   ip = iorb(ia)
   sa = sub(ia)
   na = nt(ia)
   f = gt(ia)
   pra = pr(ia)
   hpa = hprin(sa)
   ea = exp(-(hpa/gnof_c)**2)
   nda = nd(ia)
   ! The integrals are symmetric: the column ip is used (contiguous in memory)
   do ib=1,ng
     jcol(ib) = jint(iorb(ib),ip)
     kcol(ib) = kint(iorb(ib),ip)
   end do
   grad(ip) = grad(ip)+two*h1(ip)/nk+(jint(ip,ip)+two*real(selfint,dp)*na*(jint(ip,ip)-kint(ip,ip)))/nk2
   do ib=1,ng
     if (ib == ia) cycle
     same_sub = (sa == sub(ib))
     if (same_sub) then
       sigma = one
       if (pra .or. pr(ib)) sigma = -one
       dpi = zero
       if (na > tiny_den) dpi = sigma*half*sn(ib)/sn(ia)
       grad(ip) = grad(ip)+two*dpi*kcol(ib)/nk2
     else
       dpi = zero
       viap = zero
       select case (func)
       case (RDMFT_PNOF7)
         if (f > tiny_den) dpi = -half*(st(ib)/st(ia))*(one-two*na)
       case (RDMFT_PNOF7S)
         dpi = -four*(one-two*na)*gt(ib)
       case (RDMFT_GNOF)
         if (.not.(pra .and. pr(ib))) then
           sga = one
           if (pra .neqv. pr(ib)) sga = -one
           dnd = nd(ib)
           if (nda > tiny_den) dnd = nd(ib)+sga*half*sd(ib)/sd(ia)
           dsa = zero
           if (f > tiny_den) dsa = half*(st(ib)/st(ia))*(one-two*na)
           dpi = dnd*ea-dsa
           viap = dnd*nda*two*hpa/gnof_c**2
         end if
       end select
       dpr = nt(ib)*(two*jcol(ib)-kcol(ib))+dpi*kcol(ib)
       grad(ip) = grad(ip)+two*dpr/nk2
       if (viap /= zero) gp(sa) = gp(sa)+two*viap*kcol(ib)/nk2
     end if
   end do
   ! Interaction with the background orbitals (see rdmft_energy_pnof)
   dpi = zero
   if (func == RDMFT_GNOF .and. .not.pra) then
     dnd = one
     if (nda > tiny_den) dnd = one-half/sqrt(nda)
     dpi = dnd*ea
     gp(sa) = gp(sa)+two*dnd*nda*two*hpa/gnof_c**2*bk(ip)/nk2
   end if
   grad(ip) = grad(ip)+two*((two*bj(ip)-bk(ip))+dpi*bk(ip))/nk2
 end do
 ! Dependence through the occupation of the principal pair of the subspace (GNOF)
 do ig=1,table%nsub
   if (prin_orb(ig) > 0) grad(prin_orb(ig)) = grad(prin_orb(ig))+gp(ig)
 end do

 ABI_FREE(prin_orb)
 ABI_FREE(hprin)
 ABI_FREE(gp)
 ABI_FREE(iorb)
 ABI_FREE(sub)
 ABI_FREE(pr)
 ABI_FREE(nt)
 ABI_FREE(gt)
 ABI_FREE(st)
 ABI_FREE(sn)
 ABI_FREE(nd)
 ABI_FREE(sd)
 ABI_FREE(jcol)
 ABI_FREE(kcol)

end subroutine rdmft_grad_pnof
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
function rdmft_energy_pnof(func,table,nk,nbk,n,h1,jint,kint,selfint,bj,bk,bg) result(energy)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk,selfint
 type(rdmft_pnof_table_t),intent(in) :: table
 type(rdmft_bg_t),intent(in) :: bg
 real(dp) :: energy
 !arrays
 real(dp),intent(in) :: n(nbk),h1(nbk),jint(nbk,nbk),kint(nbk,nbk),bj(nbk),bk(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: ia,ib,ig,ip,iq,ng
 real(dp) :: e1,e2,na,pi_ab,pifz,xab,nda
 !arrays
 integer,allocatable :: iorb(:),sub(:)
 logical,allocatable :: pr(:)
 real(dp),allocatable :: gt(:),hprin(:),jcol(:),kcol(:),nd(:),nt(:),sd(:),sn(:),st(:)

! *************************************************************************

 if (func /= RDMFT_PNOF5 .and. func /= RDMFT_PNOF7 .and. func /= RDMFT_PNOF7S .and. func /= RDMFT_GNOF) then
   ABI_ERROR('RDMFT: unknown PNOF functional.')
 end if
 ng = table%ng
 ! Holes of the principal pair of each subspace (needed by GNOF)
 ABI_MALLOC(hprin,(table%nsub))
 hprin = one
 do ig=1,ng
   if (table%gem_principal(ig)) hprin(table%gem_sub(ig)) = one-n(table%gem_orb(ig))
 end do

 ! Per-orbital quantities of the table (the square roots and exponentials are computed once per orbital, not per pair)
 ABI_MALLOC(iorb,(ng))
 ABI_MALLOC(sub,(ng))
 ABI_MALLOC(pr,(ng))
 ABI_MALLOC(nt,(ng))
 ABI_MALLOC(gt,(ng))
 ABI_MALLOC(st,(ng))
 ABI_MALLOC(sn,(ng))
 ABI_MALLOC(nd,(ng))
 ABI_MALLOC(sd,(ng))
 ABI_MALLOC(jcol,(ng))
 ABI_MALLOC(kcol,(ng))
 do ig=1,ng
   iorb(ig) = table%gem_orb(ig)
   sub(ig) = table%gem_sub(ig)
   pr(ig) = table%gem_principal(ig)
   nt(ig) = n(iorb(ig))
   gt(ig) = nt(ig)*(one-nt(ig))
   st(ig) = sqrt(max(zero,gt(ig)))
   sn(ig) = sqrt(max(zero,nt(ig)))
   nd(ig) = nt(ig)*exp(-(hprin(sub(ig))/gnof_c)**2)
   sd(ig) = sqrt(max(zero,nd(ig)))
 end do

 e1 = zero
 e2 = zero
 ! The pair terms are symmetric in (a,b): the loop over a runs inside a column of the integrals (contiguous in memory)
 do ib=1,ng
   iq = iorb(ib)
   do ia=1,ng
     jcol(ia) = jint(iorb(ia),iq)
     kcol(ia) = kint(iorb(ia),iq)
   end do
   do ia=1,ng
     if (ia == ib) cycle
     if (sub(ia) == sub(ib)) then
       ! same subspace: no Hartree/exchange weight; pairing term
       pi_ab = sn(ia)*sn(ib)
       if (pr(ia) .or. pr(ib)) pi_ab = -pi_ab
       e2 = e2+pi_ab*kcol(ia)
     else
       xab = nt(ia)*nt(ib)
       select case (func)
       case (RDMFT_PNOF5)
         pi_ab = zero
       case (RDMFT_PNOF7)
         pi_ab = -st(ia)*st(ib)
       case (RDMFT_PNOF7S)
         pi_ab = -four*gt(ia)*gt(ib)
       case (RDMFT_GNOF)
         if (pr(ia) .and. pr(ib)) then
           pi_ab = zero
         else
           pi_ab = nd(ia)*nd(ib)-st(ia)*st(ib)
           if (pr(ia) .neqv. pr(ib)) then
             pi_ab = pi_ab-sd(ia)*sd(ib)
           else
             pi_ab = pi_ab+sd(ia)*sd(ib)
           end if
         end if
       end select
       e2 = e2+two*xab*jcol(ia)-xab*kcol(ia)+pi_ab*kcol(ia)
     end if
   end do
 end do
 do ia=1,ng
   ip = iorb(ia)
   na = nt(ia)
   e1 = e1+two*na*h1(ip)
   e2 = e2+na*jint(ip,ip)
   ! rdmft_selfint=1: same-spin self term of the orbital, n^2 (J_pp-K_pp), nonzero when K_pp /= J_pp (periodic systems)
   if (selfint == 1) e2 = e2+na**2*(jint(ip,ip)-kint(ip,ip))
   ! Interaction with the background orbitals (frozen occupied pairs, n=1, trivial subspaces), for the two orders of
   ! the pairs. Their Pi term only survives for GNOF and a non-principal pair a.
   pifz = zero
   if (func == RDMFT_GNOF .and. .not.pr(ia)) then
     nda = nd(ia)
     pifz = nda-sqrt(max(zero,nda))
   end if
   e2 = e2+two*(na*(two*bj(ip)-bk(ip))+pifz*bk(ip))
 end do
 energy = e1/nk+e2/real(nk,dp)**2+rdmft_energy_bg(func,nk,bg,zero,selfint)

 ABI_FREE(hprin)
 ABI_FREE(iorb)
 ABI_FREE(sub)
 ABI_FREE(pr)
 ABI_FREE(nt)
 ABI_FREE(gt)
 ABI_FREE(st)
 ABI_FREE(sn)
 ABI_FREE(nd)
 ABI_FREE(sd)
 ABI_FREE(jcol)
 ABI_FREE(kcol)

end function rdmft_energy_pnof
!!***

end module m_rdmft_functionals
!!***
