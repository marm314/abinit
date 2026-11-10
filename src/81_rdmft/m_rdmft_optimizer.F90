!!****m* ABINIT/m_rdmft_optimizer
!! NAME
!!  m_rdmft_optimizer
!!
!! FUNCTION
!!  Optimization of the occupation numbers of the JK-only functionals for fixed Bloch orbitals.
!!  Newton-SQP with the exact gradient and Hessian, ported from ReRDMFT (Utils/SQP): minimize E(n) with the
!!  linear equality sum_s mult_s x_s = N and the box eps <= x_s <= 1-eps. At each iteration the quadratic model
!!  is minimized with a primal active-set method (equality-constrained QP solved through its KKT system,
!!  Levenberg-Marquardt damping if singular), followed by a backtracking Armijo line search.
!!  The variables x_s are the occupations of the active orbitals, with the time-reversal partners (k and -k)
!!  sharing one variable.
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

module m_rdmft_optimizer

 use defs_basis
 use m_abicore
 use m_errors

 use m_rdmft_functionals, only : rdmft_bg_t, rdmft_energy_jk, rdmft_energy_pnof, rdmft_grad_pnof, rdmft_jk_diag, &
                                 rdmft_jk_pair, rdmft_pnof_gamma, rdmft_pnof_occ, rdmft_pnof_table_t
 use m_specialmsg,        only : wrtout
 use m_time,              only : abi_wtime

 implicit none

 private
!!***

 public :: rdmft_grad_jk
 public :: rdmft_tpnof

 !> Profile of the evaluation of the PNOF energy and gradient: 1 occupations and Jacobian, 2 energy, 3 orbital gradient,
 !> 4 contraction with the Jacobian, 5 number of evaluations
 real(dp),save :: rdmft_tpnof(5) = 0.0_dp
 public :: rdmft_optimize_jk
 public :: rdmft_optimize_pnof
!!***

 !> Inverse of the KKT matrix of the free variables of the QP (plus the equality constraint, index n+1), stored with the
 !> indices of all the variables (zero for the variables at a bound), kept up to date when a variable is fixed or freed
 type kkt_inc_t
   logical :: valid = .false.
   integer :: nupd = 0
   logical,allocatable :: fr(:)
   real(dp),allocatable :: minv(:,:)
 end type kkt_inc_t

contains
!!***

!!****f* m_rdmft_optimizer/solve_dense
!! NAME
!!  solve_dense
!!
!! FUNCTION
!!  Solve the dense linear system a x = b by an LU factorization with partial pivoting (LAPACK dgetrf).
!!
!! OUTPUT
!!  ok = false if the matrix is numerically singular
!!
!! SOURCE
subroutine solve_dense(n,a,b,x,ok)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: n
 logical,intent(out) :: ok
 !arrays
 real(dp),intent(in) :: a(n,n),b(n)
 real(dp),intent(out) :: x(n)

 !Local variables-------------------------------
 !scalars
 integer :: i,info,k
 real(dp) :: tmp
 !arrays
 integer,allocatable :: ipiv(:)
 real(dp),allocatable :: aa(:,:),colscale(:)

! *************************************************************************

 ABI_MALLOC(aa,(n,n))
 ABI_MALLOC(ipiv,(n))
 ABI_MALLOC(colscale,(n))
 aa = a
 x = b
 ok = .true.
 ! LU factorization with partial pivoting (LAPACK). Each pivot (diagonal of U) is compared with the scale of its own
 ! column (the entries of a KKT matrix differ by many orders)
 do k=1,n
   colscale(k) = max(maxval(abs(a(:,k))),tiny(one))
 end do
 call dgetrf(n,n,aa,n,ipiv,info)
 if (info /= 0) ok = .false.
 if (ok) then
   do k=1,n
     tmp = abs(aa(k,k))
     if (tmp <= 1.0e-14_dp*colscale(k) .or. tmp /= tmp) then
       ok = .false.
       exit
     end if
   end do
 end if
 if (ok) then
   ! P b, then L y = Pb (unit lower), then U x = y
   do i=1,n
     if (ipiv(i) /= i) then
       tmp = x(i); x(i) = x(ipiv(i)); x(ipiv(i)) = tmp
     end if
   end do
   call dtrsv('L','N','U',n,aa,n,x,1)
   call dtrsv('U','N','N',n,aa,n,x,1)
 else
   x = zero
 end if
 ABI_FREE(aa)
 ABI_FREE(ipiv)
 ABI_FREE(colscale)

end subroutine solve_dense
!!***

!!****f* m_rdmft_optimizer/qp_given_fixed
!! NAME
!!  qp_given_fixed
!!
!! FUNCTION
!!  Minimize g.d + d.B.d/2 with the single equality a.d = 0 for the free variables, the variables with
!!  status /= 0 being fixed at fixed_value.
!!
!! OUTPUT
!!  d = the step, rgrad = g + B d - a lambda
!!
!! SOURCE
subroutine qp_given_fixed(n,bmat,g,a,status,fixed_value,d,rgrad)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: n
 !arrays
 integer,intent(in) :: status(n)
 real(dp),intent(in) :: bmat(n,n),g(n),a(n),fixed_value(n)
 real(dp),intent(out) :: d(n),rgrad(n)

 !Local variables-------------------------------
 !scalars
 integer :: attempt,i,j,jj,kk,nf
 real(dp) :: damping,lambda,rhat
 logical :: ok
 !arrays
 integer,allocatable :: fidx(:)
 real(dp),allocatable :: kkt(:,:),rhs(:),sol(:)

! *************************************************************************

 ABI_MALLOC(fidx,(n))
 nf = 0
 do i=1,n
   if (status(i) == 0) then
     nf = nf+1; fidx(nf) = i
   end if
 end do
 d = zero
 do i=1,n
   if (status(i) /= 0) d(i) = fixed_value(i)
 end do
 lambda = zero

 if (nf > 0) then
   ABI_MALLOC(kkt,(nf+1,nf+1))
   ABI_MALLOC(rhs,(nf+1))
   ABI_MALLOC(sol,(nf+1))
   kkt = zero
   do jj=1,nf
     j = fidx(jj)
     rhs(jj) = -(g(j)+sum(bmat(j,:)*d(:),mask=(status /= 0)))
     do kk=1,nf
       kkt(jj,kk) = bmat(j,fidx(kk))
     end do
     kkt(jj,nf+1) = a(j)
     kkt(nf+1,jj) = a(j)
   end do
   rhat = -sum(a(:)*d(:),mask=(status /= 0))
   rhs(nf+1) = rhat
   damping = zero
   ok = .false.
   do attempt=1,20
     do jj=1,nf
       kkt(jj,jj) = bmat(fidx(jj),fidx(jj))+damping
     end do
     call solve_dense(nf+1,kkt,rhs,sol,ok)
     if (ok) exit
     if (damping <= zero) then
       damping = 1.0e-10_dp
     else
       damping = damping*ten
     end if
   end do
   if (.not.ok) then
     ABI_ERROR('RDMFT: KKT system of the QP is singular even with Levenberg-Marquardt damping.')
   end if
   do jj=1,nf
     d(fidx(jj)) = sol(jj)
   end do
   lambda = sol(nf+1)
   ABI_FREE(kkt)
   ABI_FREE(rhs)
   ABI_FREE(sol)
 end if

 do i=1,n
   rgrad(i) = g(i)+sum(bmat(i,:)*d(:))-a(i)*lambda
 end do
 ABI_FREE(fidx)

end subroutine qp_given_fixed
!!***

!!****f* m_rdmft_optimizer/kkt_build
!! NAME
!!  kkt_build
!!
!! FUNCTION
!!  Build the inverse of the KKT matrix [B_FF a_F; a_F^T 0] of the free variables F (status = 0) from scratch (LU).
!!
!! SOURCE

subroutine kkt_build(st,n,bmat,a,status,ok)

 !Arguments ------------------------------------
 integer,intent(in) :: n
 type(kkt_inc_t),intent(inout) :: st
 logical,intent(out) :: ok
 integer,intent(in) :: status(n)
 real(dp),intent(in) :: bmat(n,n),a(n)

 !Local variables-------------------------------
 integer :: i,info,j,jj,kk,mm,nf
 real(dp) :: tmp
 integer,allocatable :: fidx(:),idx(:),ipiv(:)
 real(dp),allocatable :: kkt(:,:),colscale(:),work(:)

! *************************************************************************

 ok = .false.
 st%valid = .false.
 st%nupd = 0
 if (.not.allocated(st%fr)) then
   ABI_MALLOC(st%fr,(n))
   ABI_MALLOC(st%minv,(n+1,n+1))
 end if
 ABI_MALLOC(fidx,(n))
 nf = 0
 do i=1,n
   if (status(i) == 0) then
     nf = nf+1; fidx(nf) = i
   end if
 end do
 st%fr = (status == 0)
 st%minv = zero
 if (nf == 0) then
   ABI_FREE(fidx)
   return
 end if
 mm = nf+1
 ABI_MALLOC(kkt,(mm,mm))
 ABI_MALLOC(colscale,(mm))
 ABI_MALLOC(ipiv,(mm))
 ABI_MALLOC(idx,(mm))
 ABI_MALLOC(work,(64*mm))
 idx(1:nf) = fidx(1:nf)
 idx(mm) = n+1
 kkt = zero
 do jj=1,nf
   j = fidx(jj)
   do kk=1,nf
     kkt(jj,kk) = bmat(j,fidx(kk))
   end do
   kkt(jj,mm) = a(j)
   kkt(mm,jj) = a(j)
 end do
 do kk=1,mm
   colscale(kk) = max(maxval(abs(kkt(:,kk))),tiny(one))
 end do
 call dgetrf(mm,mm,kkt,mm,ipiv,info)
 ok = (info == 0)
 if (ok) then
   do kk=1,mm
     tmp = abs(kkt(kk,kk))
     if (tmp <= 1.0e-14_dp*colscale(kk) .or. tmp /= tmp) ok = .false.
   end do
 end if
 if (ok) then
   call dgetri(mm,kkt,mm,ipiv,work,64*mm,info)
   ok = (info == 0)
 end if
 if (ok) then
   do kk=1,mm
     do jj=1,mm
       st%minv(idx(jj),idx(kk)) = kkt(jj,kk)
     end do
   end do
   st%valid = .true.
 end if
 ABI_FREE(fidx)
 ABI_FREE(kkt)
 ABI_FREE(colscale)
 ABI_FREE(ipiv)
 ABI_FREE(idx)
 ABI_FREE(work)

end subroutine kkt_build
!!***

!!****f* m_rdmft_optimizer/kkt_update
!! NAME
!!  kkt_update
!!
!! FUNCTION
!!  Update the inverse of the KKT matrix when the variable j is fixed (removed: Schur complement of the inverse) or
!!  freed (added: bordering formula). O(n^2) instead of the O(n^3) of a new factorization.
!!
!! SOURCE

subroutine kkt_update(st,n,bmat,a,j,ok)

 !Arguments ------------------------------------
 integer,intent(in) :: n,j
 type(kkt_inc_t),intent(inout) :: st
 logical,intent(out) :: ok
 real(dp),intent(in) :: bmat(n,n),a(n)

 !Local variables-------------------------------
 integer :: i,n1
 real(dp) :: piv,s,uw
 real(dp),allocatable :: u(:),w(:)

! *************************************************************************

 n1 = n+1
 ok = .false.
 ABI_MALLOC(u,(n1))
 ABI_MALLOC(w,(n1))
 if (st%fr(j)) then
   ! j leaves the free set: Minv <- Minv - Minv(:,j) Minv(j,:) / Minv(j,j) (symmetric matrix)
   piv = st%minv(j,j)
   if (abs(piv) > 1.0e-10_dp*max(maxval(abs(st%minv(:,j))),tiny(one)) .and. piv == piv) then
     w = st%minv(:,j)
     call dger(n1,n1,-one/piv,w,1,w,1,st%minv,n1)
     st%minv(:,j) = zero
     st%minv(j,:) = zero
     st%fr(j) = .false.
     ok = .true.
   end if
 else
   ! j enters the free set: u = [B_Fj; a_j], s = B_jj - u^T Minv u
   u = zero
   do i=1,n
     if (st%fr(i)) u(i) = bmat(i,j)
   end do
   u(n1) = a(j)
   call dgemv('N',n1,n1,one,st%minv,n1,u,1,zero,w,1)
   uw = sum(u*w)
   s = bmat(j,j)-uw
   if (abs(s) > 1.0e-10_dp*(abs(bmat(j,j))+abs(uw)) .and. s == s) then
     call dger(n1,n1,one/s,w,1,w,1,st%minv,n1)
     st%minv(:,j) = -w/s
     st%minv(j,:) = -w/s
     st%minv(j,j) = one/s
     st%fr(j) = .true.
     ok = .true.
   end if
 end if
 if (ok) st%nupd = st%nupd+1
 ABI_FREE(u)
 ABI_FREE(w)

end subroutine kkt_update
!!***

!!****f* m_rdmft_optimizer/qp_given_fixed_inc
!! NAME
!!  qp_given_fixed_inc
!!
!! FUNCTION
!!  Same as qp_given_fixed (without damping) with the inverse of the KKT matrix of the free variables kept up to date
!!  between the calls. ok = false if the matrix is (nearly) singular or the solution is not accurate: the caller then uses
!!  qp_given_fixed.
!!
!! SOURCE

subroutine qp_given_fixed_inc(st,n,bmat,g,a,status,fixed_value,d,rgrad,ok)

 !Arguments ------------------------------------
 integer,intent(in) :: n
 type(kkt_inc_t),intent(inout) :: st
 logical,intent(out) :: ok
 integer,intent(in) :: status(n)
 real(dp),intent(in) :: bmat(n,n),g(n),a(n),fixed_value(n)
 real(dp),intent(out) :: d(n),rgrad(n)

 !Local variables-------------------------------
 integer :: attempt,i,n1
 real(dp) :: lambda,res,scale
 logical :: ok1
 real(dp),allocatable :: dfix(:),rhsf(:),solf(:),tmp(:),tmp2(:)

! *************************************************************************

 n1 = n+1
 ok = .true.
 if (.not.st%valid) call kkt_build(st,n,bmat,a,status,ok)
 if (ok .and. st%valid) then
   do i=1,n
     if ((status(i) == 0) .neqv. st%fr(i)) then
       call kkt_update(st,n,bmat,a,i,ok1)
       if (.not.ok1) then
         call kkt_build(st,n,bmat,a,status,ok)
         exit
       end if
     end if
   end do
   if (ok .and. st%nupd > 100) call kkt_build(st,n,bmat,a,status,ok)
 end if
 if (.not.(ok .and. st%valid)) then
   ok = .false.
   st%valid = .false.
   return
 end if

 ABI_MALLOC(dfix,(n))
 ABI_MALLOC(rhsf,(n1))
 ABI_MALLOC(solf,(n1))
 ABI_MALLOC(tmp,(n))
 ABI_MALLOC(tmp2,(n))
 do attempt=1,2
   dfix = zero
   do i=1,n
     if (status(i) /= 0) dfix(i) = fixed_value(i)
   end do
   call dgemv('N',n,n,one,bmat,n,dfix,1,zero,tmp,1)
   rhsf = zero
   do i=1,n
     if (st%fr(i)) rhsf(i) = -(g(i)+tmp(i))
   end do
   rhsf(n1) = -sum(a*dfix)
   call dgemv('N',n1,n1,one,st%minv,n1,rhsf,1,zero,solf,1)
   ! residual of the KKT system on the free variables
   lambda = solf(n1)
   call dgemv('N',n,n,one,bmat,n,solf,1,zero,tmp2,1)
   res = abs(sum(a(1:n)*solf(1:n))-rhsf(n1))
   do i=1,n
     if (st%fr(i)) res = max(res,abs(tmp2(i)+a(i)*lambda-rhsf(i)))
   end do
   scale = max(one,maxval(abs(rhsf)))
   if (res <= 1.0e-9_dp*scale) exit
   ! inaccurate inverse: new factorization, once
   ok = .false.
   if (attempt == 1) call kkt_build(st,n,bmat,a,status,ok)
   if (.not.ok) exit
 end do
 if (ok) then
   d = dfix
   do i=1,n
     if (st%fr(i)) d(i) = solf(i)
   end do
   call dgemv('N',n,n,one,bmat,n,d,1,zero,tmp,1)
   do i=1,n
     rgrad(i) = g(i)+tmp(i)-a(i)*lambda
   end do
 else
   st%valid = .false.
 end if
 ABI_FREE(dfix)
 ABI_FREE(rhsf)
 ABI_FREE(solf)
 ABI_FREE(tmp)
 ABI_FREE(tmp2)

end subroutine qp_given_fixed_inc
!!***

!!****f* m_rdmft_optimizer/qp_box_eq
!! NAME
!!  qp_box_eq
!!
!! FUNCTION
!!  Minimize g.d + d.B.d/2 subject to a.d = 0 and lb <= d <= ub (d=0 is feasible), primal active-set method.
!!
!! OUTPUT
!!  d = solution, converged = true if the KKT conditions are satisfied
!!
!! SOURCE
subroutine qp_box_eq(n,bmat,g,a,lb,ub,d,converged)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: n
 logical,intent(out) :: converged
 !arrays
 real(dp),intent(in) :: bmat(n,n),g(n),a(n),lb(n),ub(n)
 real(dp),intent(out) :: d(n)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: max_iter=500
 integer :: blocking_bound,blocking_idx,i,iter,worst_idx
 real(dp),parameter :: tolz=1.0e-14_dp
 real(dp) :: alpha,cand,tolm,worst_violation
 !arrays
 integer,allocatable :: status(:)
 logical :: ok_inc
 type(kkt_inc_t) :: st
 real(dp),allocatable :: dsol(:),fixed_value(:),p(:),rgrad(:)

! *************************************************************************

 ABI_MALLOC(status,(n))
 ABI_MALLOC(dsol,(n))
 ABI_MALLOC(fixed_value,(n))
 ABI_MALLOC(p,(n))
 ABI_MALLOC(rgrad,(n))
 d = zero
 status = 0
 fixed_value = zero
 converged = .false.
 ! Tolerance on the multipliers of the bounds, relative to the scale of the gradient
 tolm = 1.0e-12_dp*max(one,maxval(abs(g)))

 do iter=1,max_iter
   do i=1,n
     if (status(i) == -1) fixed_value(i) = lb(i)
     if (status(i) == 1) fixed_value(i) = ub(i)
   end do
   call qp_given_fixed_inc(st,n,bmat,g,a,status,fixed_value,dsol,rgrad,ok_inc)
   if (.not.ok_inc) call qp_given_fixed(n,bmat,g,a,status,fixed_value,dsol,rgrad)
   p = dsol-d

   if (sqrt(sum(p**2)) <= tolz) then
     ! Optimal for the working set; check the signs of the multipliers of the active bounds
     worst_idx = 0
     worst_violation = tolm
     do i=1,n
       if (status(i) == -1 .and. -rgrad(i) > worst_violation) then
         worst_violation = -rgrad(i); worst_idx = i
       else if (status(i) == 1 .and. rgrad(i) > worst_violation) then
         worst_violation = rgrad(i); worst_idx = i
       end if
     end do
     if (worst_idx == 0) then
       converged = .true.
       exit
     end if
     status(worst_idx) = 0
     cycle
   end if

   ! Ratio test
   alpha = one
   blocking_idx = 0
   blocking_bound = 0
   do i=1,n
     if (status(i) /= 0) cycle
     if (p(i) > tolz) then
       cand = (ub(i)-d(i))/p(i)
       if (cand < alpha) then
         alpha = cand; blocking_idx = i; blocking_bound = 1
       end if
     else if (p(i) < -tolz) then
       cand = (lb(i)-d(i))/p(i)
       if (cand < alpha) then
         alpha = cand; blocking_idx = i; blocking_bound = -1
       end if
     end if
   end do
   if (alpha < zero) alpha = zero
   d = d+alpha*p
   if (blocking_idx > 0 .and. alpha < one) then
     status(blocking_idx) = blocking_bound
     if (blocking_bound == 1) then
       d(blocking_idx) = ub(blocking_idx)
     else
       d(blocking_idx) = lb(blocking_idx)
     end if
   end if
 end do

 ABI_FREE(status)
 ABI_SFREE(st%fr)
 ABI_SFREE(st%minv)
 ABI_FREE(dsol)
 ABI_FREE(fixed_value)
 ABI_FREE(p)
 ABI_FREE(rgrad)

end subroutine qp_box_eq
!!***

!!****f* m_rdmft_optimizer/rdmft_optimize_jk
!! NAME
!!  rdmft_optimize_jk
!!
!! FUNCTION
!!  Optimize the occupation numbers of the active orbitals for a JK-only functional.
!!
!! INPUTS
!!  func = JK-only functional
!!  nk = number of k-points of the full BZ
!!  nbk = number of free orbitals (all of them are variables)
!!  h1(nbk), jint(nbk,nbk), kint(nbk,nbk) = core Hamiltonian diagonal, J and K integrals
!!  alpha = exponent of the POWER functional
!!  partner(nbk) = time-reversal partner of each orbital (itself for the time-reversal-invariant k-points)
!!  n0(nbk) = starting occupations
!!  ntarget = sum of the occupations of the free orbitals
!!  bj(nbk), bk(nbk), bg = background orbitals (frozen occupied): sums of the integrals and constants
!!  eps = the occupations are kept in [eps,1-eps]
!!
!! OUTPUT
!!  nout(nbk) = optimized occupations (the frozen ones are nfix)
!!  energy = optimized electronic energy per unit cell
!!  converged, niter
!!
!! SOURCE
subroutine rdmft_optimize_jk(func,nk,nbk,h1,jint,kint,alpha,partner,n0,ntarget,eps,bj,bk,bg,&
                             nout,energy,converged,niter,verbose,e_shift)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk
 real(dp),intent(in) :: alpha,ntarget,eps
 logical,optional,intent(in) :: verbose
 real(dp),optional,intent(in) :: e_shift
 type(rdmft_bg_t),intent(in) :: bg
 real(dp),intent(out) :: energy
 logical,intent(out) :: converged
 integer,intent(out) :: niter
 !arrays
 integer,intent(in) :: partner(nbk)
 real(dp),intent(in) :: h1(nbk),jint(nbk,nbk),kint(nbk,nbk),n0(nbk),bj(nbk),bk(nbk)
 real(dp),intent(out) :: nout(nbk)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: max_iter=100,max_ls=50
 integer :: attempt,iter,ils,ip,is,nv
 real(dp),parameter :: armijo_c1=1.0e-4_dp,step_tol=1.0e-10_dp
 real(dp) :: alph,dirder,f0,lambda,ftrial,headroom_tot,resid,gr,fa,faa,fab,ff,fd,fd1,fd2
 logical :: qp_ok,talk
 real(dp) :: eshift
 character(len=500) :: msg
 !arrays
 integer,allocatable :: var_of(:)
 real(dp),allocatable :: aeq(:),bmat(:,:),bshift(:,:),d(:),g(:),headroom(:),lb(:),ub(:),x(:),xtrial(:),gsum(:),hdiag(:)

! *************************************************************************

 talk = .true.
 if (present(verbose)) talk = verbose
 eshift = zero
 if (present(e_shift)) eshift = e_shift
 if (talk) then
   write(msg,'(3a)') ch10,' Occupation numbers, Newton-SQP:   iter     Etotal(Ha)       |step|        g.d       alpha',&
     '    lambda  qp_ok'
   call wrtout([std_out],msg)
 end if

 ! Variables: one per orbital, shared with the time-reversal partner
 ABI_MALLOC(var_of,(nbk))
 var_of = 0
 nv = 0
 do ip=1,nbk
   if (var_of(ip) /= 0) cycle
   nv = nv+1
   var_of(ip) = nv
   var_of(partner(ip)) = nv
 end do
 ABI_MALLOC(aeq,(nv))
 aeq = zero
 do ip=1,nbk
   if (var_of(ip) > 0) aeq(var_of(ip)) = aeq(var_of(ip))+one
 end do

 ABI_MALLOC(x,(nv))
 ABI_MALLOC(xtrial,(nv))
 ABI_MALLOC(lb,(nv))
 ABI_MALLOC(ub,(nv))
 ABI_MALLOC(d,(nv))
 ABI_MALLOC(g,(nv))
 ABI_MALLOC(bmat,(nv,nv))
 ABI_MALLOC(bshift,(nv,nv))
 ABI_MALLOC(headroom,(nv))
 ABI_MALLOC(gsum,(nbk))
 ABI_MALLOC(hdiag,(nbk))
 lb = eps
 ub = one-eps

 ! Starting point: aufbau occupations moved into the interior of the box, keeping the sum
 do ip=1,nbk
   x(var_of(ip)) = n0(ip)
 end do
 x = min(one-eps,max(eps,x))
 resid = ntarget-sum(aeq*x)
 if (abs(resid) > zero) then
   do is=1,nv
     if (resid > zero) then
       headroom(is) = (one-eps)-x(is)
     else
       headroom(is) = x(is)-eps
     end if
   end do
   headroom_tot = sum(aeq*headroom)
   if (headroom_tot > zero) x = x+resid*headroom/headroom_tot
 end if

 converged = .false.
 niter = 0
 do iter=1,max_iter
   niter = iter
   call expand(x,nout)
   call gradient_hessian(nout,g,bmat)
   if (any(g /= g) .or. any(bmat /= bmat) .or. maxval(abs(bmat)) > huge(one)/ten .or. maxval(abs(g)) > huge(one)/ten) then
     write(msg,'(a,i0,a,es12.4,a,es12.4,a,es12.4)') 'RDMFT: non-finite gradient or Hessian at iteration ',iter,&
       ': max|g| = ',maxval(abs(g)),', max|B| = ',maxval(abs(bmat)),', x range ',minval(x)
     ABI_ERROR(msg)
   end if
   ! The Hessian is not positive definite in general: if the QP does not converge or gives an ascent direction,
   ! a growing multiple of the identity (Levenberg-Marquardt) is added until it gives a descent direction.
   lambda = zero
   do attempt=1,30
     bshift = bmat
     do is=1,nv
       bshift(is,is) = bshift(is,is)+lambda
     end do
     call qp_box_eq(nv,bshift,g,aeq,lb-x,ub-x,d,qp_ok)
     if (qp_ok .and. sqrt(sum(d**2)) <= step_tol) exit
     if (qp_ok .and. sum(g*d) < -1.0e-14_dp) exit
     if (lambda <= zero) then
       lambda = 1.0e-6_dp*max(maxval(abs(bmat)),one)
     else
       lambda = lambda*ten
     end if
   end do

   if (qp_ok .and. sqrt(sum(d**2)) <= step_tol) then
     converged = .true.
     exit
   end if

   ! Backtracking Armijo line search along d
   dirder = sum(g*d)
   call expand(x,nout)
   f0 = evalue(nout)
   alph = one
   xtrial = x
   do ils=1,max_ls
     xtrial = x+alph*d
     call expand(xtrial,nout)
     ftrial = evalue(nout)
     if (ftrial <= f0+armijo_c1*alph*dirder) exit
     alph = alph*half
   end do
   x = xtrial
   if (talk) then
     write(msg,'(a,i5,f20.12,es12.3,es11.3,2es10.2,l4)') '  OCC',iter,ftrial+eshift,sqrt(sum(d**2)),dirder,alph,lambda,qp_ok
     call wrtout([std_out],msg)
   end if
 end do

 call expand(x,nout)
 energy = evalue(nout)

 ABI_FREE(var_of)
 ABI_FREE(aeq)
 ABI_FREE(x)
 ABI_FREE(xtrial)
 ABI_FREE(lb)
 ABI_FREE(ub)
 ABI_FREE(d)
 ABI_FREE(g)
 ABI_FREE(bmat)
 ABI_FREE(bshift)
 ABI_FREE(headroom)
 ABI_FREE(gsum)
 ABI_FREE(hdiag)

contains

 subroutine expand(xv,nfull)
   real(dp),intent(in) :: xv(nv)
   real(dp),intent(out) :: nfull(nbk)
   integer :: jp
   do jp=1,nbk
     nfull(jp) = xv(var_of(jp))
   end do
 end subroutine expand

 function evalue(nfull) result(e)
   real(dp),intent(in) :: nfull(nbk)
   real(dp) :: e
   e = rdmft_energy_jk(func,nk,nbk,nfull,h1,jint,kint,alpha,bj,bk,bg)
 end function evalue

 subroutine gradient_hessian(nfull,gx,hx)
   ! Gradient and Hessian with respect to the variables x_s
   real(dp),intent(in) :: nfull(nbk)
   real(dp),intent(out) :: gx(nv),hx(nv,nv)
   integer :: jp,jq
   real(dp) :: nk2
   nk2 = real(nk,dp)**2
   gx = zero
   hx = zero
   gsum = zero
   hdiag = zero
   do jp=1,nbk
     do jq=1,nbk
       if (jq == jp) cycle
       call rdmft_jk_pair(func,nfull(jp),nfull(jq),jint(jp,jq),kint(jp,jq),alpha,ff,fa,faa,fab)
       gsum(jp) = gsum(jp)+two*fa
       hdiag(jp) = hdiag(jp)+two*faa
       hx(var_of(jp),var_of(jq)) = hx(var_of(jp),var_of(jq))+two*fab/nk2
     end do
     ! Interaction with the background orbitals (occupation 1)
     call rdmft_jk_pair(func,nfull(jp),one,bj(jp),bk(jp),alpha,ff,fa,faa,fab)
     gsum(jp) = gsum(jp)+two*fa
     hdiag(jp) = hdiag(jp)+two*faa
     call rdmft_jk_diag(func,nfull(jp),jint(jp,jp),kint(jp,jp),alpha,fd,fd1,fd2)
     gr = two*h1(jp)/nk+(fd1+gsum(jp))/nk2
     gx(var_of(jp)) = gx(var_of(jp))+gr
     hx(var_of(jp),var_of(jp)) = hx(var_of(jp),var_of(jp))+(fd2+hdiag(jp))/nk2
   end do
 end subroutine gradient_hessian

end subroutine rdmft_optimize_jk
!!***

!!****f* m_rdmft_optimizer/rdmft_grad_jk
!! NAME
!!  rdmft_grad_jk
!!
!! FUNCTION
!!  Exact gradient of the JK-only energy with respect to the occupation of each orbital (the other occupations fixed,
!!  without the ties between k and -k): g_p = dE/dn_p.
!!
!! INPUTS
!!  func, nk, nbk, h1, jint, kint, alpha, bj, bk = as rdmft_energy_jk
!!  n(nbk) = occupations
!!
!! OUTPUT
!!  gorb(nbk) = dE/dn_p
!!
!! SOURCE
subroutine rdmft_grad_jk(func,nk,nbk,h1,jint,kint,alpha,n,bj,bk,gorb)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk
 real(dp),intent(in) :: alpha
 !arrays
 real(dp),intent(in) :: n(nbk),h1(nbk),jint(nbk,nbk),kint(nbk,nbk),bj(nbk),bk(nbk)
 real(dp),intent(out) :: gorb(nbk)

 !Local variables-------------------------------
 !scalars
 integer :: jp,jq
 real(dp) :: ff,fa,faa,fab,fd,fd1,fd2,gsum,nk2

! *************************************************************************

 nk2 = real(nk,dp)**2
 do jp=1,nbk
   gsum = zero
   do jq=1,nbk
     if (jq == jp) cycle
     call rdmft_jk_pair(func,n(jp),n(jq),jint(jp,jq),kint(jp,jq),alpha,ff,fa,faa,fab)
     gsum = gsum+two*fa
   end do
   call rdmft_jk_pair(func,n(jp),one,bj(jp),bk(jp),alpha,ff,fa,faa,fab)
   gsum = gsum+two*fa
   call rdmft_jk_diag(func,n(jp),jint(jp,jp),kint(jp,jp),alpha,fd,fd1,fd2)
   gorb(jp) = two*h1(jp)/nk+(fd1+gsum)/nk2
 end do

end subroutine rdmft_grad_jk
!!***

!!****f* m_rdmft_optimizer/rdmft_optimize_pnof
!! NAME
!!  rdmft_optimize_pnof
!!
!! FUNCTION
!!  Optimize the occupation numbers of a PNOF functional through the angles gamma of its subspaces (see
!!  rdmft_pnof_occ), which make the constraints and the box automatic: unconstrained limited-memory BFGS
!!  (two-loop recursion, 10 pairs) with a backtracking Armijo line search, starting from all gamma = pi/4
!!  as in ReRDMFT.
!!
!! INPUTS
!!  func = PNOF functional
!!  table <type(rdmft_pnof_table_t)> = pairs, subspaces, classes and variables
!!  nk, nbk, h1, jint, kint, selfint = as rdmft_energy_pnof
!!  bj(nbk), bk(nbk), bg = background orbitals (frozen occupied): sums of the integrals and constants
!!  gam_start(ngam) = optional: starting angles (they take precedence over use_start/n0)
!!  use_start, n0(nbk) = if use_start, the angles are those that give the occupations n0 (as closely as possible),
!!    else all the angles start at pi/4
!!  gtol = convergence threshold on the largest component of the gradient
!!  maxit = maximum number of iterations
!!  check_grad = if true, compare the analytic gradient with finite differences at the start
!!
!! OUTPUT
!!  nout(nbk) = optimized occupations of all the orbitals
!!  energy = optimized electronic energy per unit cell
!!  converged, niter
!!  gam_out(ngam) = optional: final angles
!!
!! SOURCE
subroutine rdmft_optimize_pnof(func,table,nk,nbk,h1,jint,kint,selfint,bj,bk,bg,use_start,n0,gtol,maxit,check_grad,&
                               nout,energy,converged,niter,verbose,e_shift,gam_start,gam_out)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: func,nk,nbk,selfint,maxit
 logical,optional,intent(in) :: verbose
 real(dp),optional,intent(in) :: e_shift
 real(dp),intent(in) :: gtol
 logical,intent(in) :: check_grad,use_start
 type(rdmft_pnof_table_t),intent(in) :: table
 type(rdmft_bg_t),intent(in) :: bg
 real(dp),intent(out) :: energy
 logical,intent(out) :: converged
 integer,intent(out) :: niter
 !arrays
 real(dp),intent(in) :: h1(nbk),jint(nbk,nbk),kint(nbk,nbk),bj(nbk),bk(nbk),n0(nbk)
 real(dp),intent(out) :: nout(nbk)
 real(dp),optional,intent(in) :: gam_start(:)
 real(dp),optional,intent(out) :: gam_out(:)

 !Local variables-------------------------------
 !scalars
 integer,parameter :: nhist=10,max_ls=60,nstag_max=10
 integer :: i,iter,ils,ih,nstag,nstored,ngam
 real(dp),parameter :: pi4=0.78539816339744830962_dp,c1=1.0e-4_dp
 real(dp) :: alph,beta,dirder,eshift,f0,ftrial,gamma_scale,sy,yy,fp,fm,step
 logical :: talk
 character(len=500) :: msg
 !arrays
 real(dp),allocatable :: gam(:),gnew(:),g(:),dvec(:),gtrial(:),xtrial(:),smat(:,:),ymat(:,:),rho(:),alpha_l(:),q(:)

! *************************************************************************

 talk = .true.
 if (present(verbose)) talk = verbose
 eshift = zero
 if (present(e_shift)) eshift = e_shift
 if (talk) then
   write(msg,'(2a)') ch10,' Occupation numbers, L-BFGS:   iter     Etotal(Ha)       max|g|'
   call wrtout([std_out],msg)
 end if

 ngam = table%ngam
 ABI_MALLOC(gam,(ngam))
 ABI_MALLOC(g,(ngam))
 ABI_MALLOC(gnew,(ngam))
 ABI_MALLOC(dvec,(ngam))
 ABI_MALLOC(gtrial,(ngam))
 ABI_MALLOC(xtrial,(ngam))
 ABI_MALLOC(smat,(ngam,nhist))
 ABI_MALLOC(ymat,(ngam,nhist))
 ABI_MALLOC(rho,(nhist))
 ABI_MALLOC(alpha_l,(nhist))
 ABI_MALLOC(q,(ngam))

 gam = pi4
 if (use_start) call rdmft_pnof_gamma(table,nbk,n0,gam)
 if (present(gam_start)) then
   if (size(gam_start) == ngam) gam = gam_start
 end if
 call eval(gam,f0,g)

 if (check_grad) then
   do i=1,min(3,ngam)
     step = 1.0e-5_dp
     xtrial = gam; xtrial(i) = gam(i)+step; call eval(xtrial,fp,gtrial)
     xtrial = gam; xtrial(i) = gam(i)-step; call eval(xtrial,fm,gtrial)
     write(msg,'(a,i0,a,es16.8,a,es16.8)') ' RDMFT PNOF gradient check, angle ',i,': analytic ',g(i),&
       ' finite differences ',(fp-fm)/(two*step)
     call wrtout([ab_out,std_out],msg)
   end do
 end if

 nstored = 0
 ih = 0
 nstag = 0
 converged = .false.
 niter = 0
 do iter=1,maxit
   niter = iter
   if (maxval(abs(g)) < gtol) then
     converged = .true.
     exit
   end if

   ! Two-loop recursion: d = -H g
   q = g
   do i=nstored,1,-1
     alpha_l(i) = rho(i)*sum(smat(:,i)*q)
     q = q-alpha_l(i)*ymat(:,i)
   end do
   gamma_scale = one
   if (nstored > 0) then
     sy = sum(smat(:,nstored)*ymat(:,nstored)); yy = sum(ymat(:,nstored)**2)
     if (yy > zero) gamma_scale = sy/yy
   end if
   q = gamma_scale*q
   do i=1,nstored
     beta = rho(i)*sum(ymat(:,i)*q)
     q = q+smat(:,i)*(alpha_l(i)-beta)
   end do
   dvec = -q
   dirder = sum(g*dvec)
   if (dirder >= zero) then
     ! not a descent direction: restart with the steepest descent
     nstored = 0
     dvec = -g
     dirder = -sum(g**2)
   end if

   ! Backtracking Armijo line search
   alph = one
   if (nstored == 0) alph = min(one,one/max(maxval(abs(g)),tiny(one)))
   do ils=1,max_ls
     xtrial = gam+alph*dvec
     call eval(xtrial,ftrial,gtrial)
     if (ftrial <= f0+c1*alph*dirder) exit
     alph = alph*half
   end do

   ! Update of the history if the curvature is positive
   sy = sum((xtrial-gam)*(gtrial-g))
   if (sy > 1.0e-10_dp*sqrt(sum((xtrial-gam)**2)*sum((gtrial-g)**2))) then
     if (nstored == nhist) then
       smat(:,1:nhist-1) = smat(:,2:nhist)
       ymat(:,1:nhist-1) = ymat(:,2:nhist)
       rho(1:nhist-1) = rho(2:nhist)
       nstored = nhist-1
     end if
     nstored = nstored+1
     smat(:,nstored) = xtrial-gam
     ymat(:,nstored) = gtrial-g
     rho(nstored) = one/sy
   end if
   ! Noise floor: the energy no longer decreases (to the precision of its evaluation) although the gradient is above
   ! gtol, because of the round-off of the energy and of the gradient. The optimization is stopped after nstag_max such
   ! iterations if the gradient is already small (energy error ~ |g|^2/curvature, far below 1e-12 Ha).
   if (f0-ftrial <= 1.0e-13_dp*max(one,abs(ftrial))) then
     nstag = nstag+1
   else
     nstag = 0
   end if
   gam = xtrial
   g = gtrial
   f0 = ftrial
   if (nstag >= nstag_max .and. maxval(abs(g)) < 1.0e-5_dp) then
     converged = .true.
     if (talk) then
       write(msg,'(a,i5,f20.12,es12.3,a)') '  OCC',iter,f0+eshift,maxval(abs(g)),'  (energy stationary to round-off)'
       call wrtout([std_out],msg)
     end if
     exit
   end if
   if (talk .and. (mod(iter,10) == 1 .or. maxval(abs(g)) < gtol)) then
     write(msg,'(a,i5,f20.12,es12.3)') '  OCC',iter,f0+eshift,maxval(abs(g))
     call wrtout([std_out],msg)
   end if
 end do

 call occupations(gam,nout)
 energy = rdmft_energy_pnof(func,table,nk,nbk,nout,h1,jint,kint,selfint,bj,bk,bg)
 if (present(gam_out)) gam_out(1:ngam) = gam

 ABI_FREE(gam)
 ABI_FREE(g)
 ABI_FREE(gnew)
 ABI_FREE(dvec)
 ABI_FREE(gtrial)
 ABI_FREE(xtrial)
 ABI_FREE(smat)
 ABI_FREE(ymat)
 ABI_FREE(rho)
 ABI_FREE(alpha_l)
 ABI_FREE(q)

contains

 subroutine occupations(gv,nfull)
   ! Occupations of all the orbitals from the angles
   real(dp),intent(in) :: gv(ngam)
   real(dp),intent(out) :: nfull(nbk)
   integer :: jp
   real(dp),allocatable :: nv(:),jac(:,:)
   ABI_MALLOC(nv,(table%nvar))
   ABI_MALLOC(jac,(table%nvar,ngam))
   call rdmft_pnof_occ(table,gv,nv,jac)
   nfull = zero
   do jp=1,nbk
     if (table%orb_var(jp) > 0) nfull(jp) = nv(table%orb_var(jp))
   end do
   ABI_FREE(nv)
   ABI_FREE(jac)
 end subroutine occupations

 subroutine eval(gv,e,gg)
   ! Energy and gradient with respect to the angles
   real(dp),intent(in) :: gv(ngam)
   real(dp),intent(out) :: e,gg(ngam)
   integer :: jp,jv,ji
   real(dp) :: tev
   real(dp),allocatable :: nv(:),jac(:,:),gorb(:),gvar(:),nfull(:)
   tev = abi_wtime()
   ABI_MALLOC(nv,(table%nvar))
   ABI_MALLOC(jac,(table%nvar,ngam))
   ABI_MALLOC(gorb,(nbk))
   ABI_MALLOC(gvar,(table%nvar))
   ABI_MALLOC(nfull,(nbk))
   call rdmft_pnof_occ(table,gv,nv,jac)
   nfull = zero
   do jp=1,nbk
     if (table%orb_var(jp) > 0) nfull(jp) = nv(table%orb_var(jp))
   end do
   rdmft_tpnof(1) = rdmft_tpnof(1)+abi_wtime()-tev
   tev = abi_wtime()
   e = rdmft_energy_pnof(func,table,nk,nbk,nfull,h1,jint,kint,selfint,bj,bk,bg)
   rdmft_tpnof(2) = rdmft_tpnof(2)+abi_wtime()-tev
   tev = abi_wtime()
   call rdmft_grad_pnof(func,table,nk,nbk,nfull,h1,jint,kint,selfint,bj,bk,gorb)
   rdmft_tpnof(3) = rdmft_tpnof(3)+abi_wtime()-tev
   tev = abi_wtime()
   gvar = zero
   do jp=1,nbk
     if (table%orb_var(jp) > 0) gvar(table%orb_var(jp)) = gvar(table%orb_var(jp))+gorb(jp)
   end do
   gg = zero
   do ji=1,ngam
     do jv=1,table%nvar
       gg(ji) = gg(ji)+gvar(jv)*jac(jv,ji)
     end do
   end do
   ABI_FREE(nv)
   ABI_FREE(jac)
   ABI_FREE(gorb)
   rdmft_tpnof(4) = rdmft_tpnof(4)+abi_wtime()-tev
   rdmft_tpnof(5) = rdmft_tpnof(5)+one
   ABI_FREE(gvar)
   ABI_FREE(nfull)
 end subroutine eval

end subroutine rdmft_optimize_pnof
!!***

end module m_rdmft_optimizer
!!***
