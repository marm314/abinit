!!****m* ABINIT/m_rdmft_adam
!! NAME
!!  m_rdmft_adam
!!
!! FUNCTION
!!  ADAM optimizer for the rotation of the orbitals in RDMFT (complex version of the ADAM of ReRDMFT), and the
!!  unitary matrix exponential of the rotation generators.
!!  Each parameter is complex (kappa_pq = t + i y): the first moment m1 is complex, the second moment m2 is real
!!  and shared by the real and imaginary parts (it accumulates |g|^2), with the running maximum of the bias-corrected
!!  moment (AMSGrad). The step is  -lr * m1hat / (sqrt(m2max) + eps).
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

module m_rdmft_adam

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private
!!***

 !> State of the ADAM optimizer
 type,public :: rdmft_adam_t
   integer :: n = 0
   integer :: nstep = 0
   real(dp) :: beta1 = 0.7_dp
   real(dp) :: beta2 = 0.9_dp
   real(dp) :: eps = 1.0e-10_dp
   complex(dpc),allocatable :: m1(:)
   real(dp),allocatable :: m2(:)
   real(dp),allocatable :: m2max(:)
 contains
   procedure :: init => rdmft_adam_init
   procedure :: clean => rdmft_adam_clean
   procedure :: step => rdmft_adam_step
   procedure :: free => rdmft_adam_free
 end type rdmft_adam_t

 public :: rdmft_expm
!!***

contains
!!***

!!****f* m_rdmft_adam/rdmft_adam_init
!! NAME
!!  rdmft_adam_init
!!
!! FUNCTION
!!  Allocate the moments for n complex parameters.
!!
!! INPUTS
!!  n = number of parameters
!!  beta1, beta2 = decay rates of the first and second moments
!!  eps = small constant added to the denominator of the step: the components of the gradient below it (zero by symmetry,
!!        only round-off) are not amplified into steps of the size of the learning rate
!!
!! SOURCE
subroutine rdmft_adam_init(adam,n,beta1,beta2,eps)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_adam_t),intent(inout) :: adam
 integer,intent(in) :: n
 real(dp),intent(in) :: beta1,beta2,eps

! *************************************************************************

 call adam%free()
 adam%n = n
 adam%beta1 = beta1
 adam%beta2 = beta2
 adam%eps = eps
 ABI_MALLOC(adam%m1,(n))
 ABI_MALLOC(adam%m2,(n))
 ABI_MALLOC(adam%m2max,(n))
 call adam%clean()

end subroutine rdmft_adam_init
!!***

!!****f* m_rdmft_adam/rdmft_adam_clean
!! NAME
!!  rdmft_adam_clean
!!
!! FUNCTION
!!  Reset the moments and the step counter.
!!
!! SOURCE
subroutine rdmft_adam_clean(adam)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_adam_t),intent(inout) :: adam

! *************************************************************************

 adam%nstep = 0
 adam%m1 = czero
 adam%m2 = zero
 adam%m2max = zero

end subroutine rdmft_adam_clean
!!***

!!****f* m_rdmft_adam/rdmft_adam_step
!! NAME
!!  rdmft_adam_step
!!
!! FUNCTION
!!  One ADAM step for the complex gradient g (g = dE/dRe + i dE/dIm of each parameter).
!!
!! INPUTS
!!  grad(n) = complex gradient
!!  lr = learning rate
!!
!! OUTPUT
!!  dx(n) = change of the parameters
!!
!! SOURCE
subroutine rdmft_adam_step(adam,grad,lr,dx)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_adam_t),intent(inout) :: adam
 real(dp),intent(in) :: lr
 !arrays
 complex(dpc),intent(in) :: grad(:)
 complex(dpc),intent(out) :: dx(:)

 !Local variables-------------------------------
 !scalars
 integer :: i
 real(dp) :: bc1,bc2,v2

! *************************************************************************

 adam%nstep = adam%nstep+1
 bc1 = one-adam%beta1**adam%nstep
 bc2 = one-adam%beta2**adam%nstep
 do i=1,adam%n
   adam%m1(i) = adam%beta1*adam%m1(i)+(one-adam%beta1)*grad(i)
   adam%m2(i) = adam%beta2*adam%m2(i)+(one-adam%beta2)*(real(grad(i),dp)**2+aimag(grad(i))**2)
   v2 = adam%m2(i)/bc2
   adam%m2max(i) = max(adam%m2max(i),v2)
   dx(i) = -lr*(adam%m1(i)/bc1)/(sqrt(adam%m2max(i))+adam%eps)
 end do

end subroutine rdmft_adam_step
!!***

!!****f* m_rdmft_adam/rdmft_adam_free
!! NAME
!!  rdmft_adam_free
!!
!! FUNCTION
!!  Deallocate the optimizer.
!!
!! SOURCE
subroutine rdmft_adam_free(adam)

 !Arguments ------------------------------------
 !scalars
 class(rdmft_adam_t),intent(inout) :: adam

! *************************************************************************

 ABI_SFREE(adam%m1)
 ABI_SFREE(adam%m2)
 ABI_SFREE(adam%m2max)
 adam%n = 0

end subroutine rdmft_adam_free
!!***

!!****f* m_rdmft_adam/rdmft_expm
!! NAME
!!  rdmft_expm
!!
!! FUNCTION
!!  Exponential of an anti-Hermitian matrix (scaling and squaring with a Taylor series): a unitary matrix.
!!
!! INPUTS
!!  kappa(n,n) = anti-Hermitian matrix
!!
!! OUTPUT
!!  umat(n,n) = exp(kappa)
!!
!! SOURCE
function rdmft_expm(kappa) result(umat)

 !Arguments ------------------------------------
 !arrays
 complex(dpc),intent(in) :: kappa(:,:)
 complex(dpc) :: umat(size(kappa,1),size(kappa,1))

 !Local variables-------------------------------
 !scalars
 integer :: i,is,nb,nsq
 real(dp) :: anorm
 !arrays
 complex(dpc),allocatable :: a(:,:),term(:,:)

! *************************************************************************

 nb = size(kappa,1)
 anorm = maxval(abs(kappa))*nb
 nsq = 0
 do while (anorm > 0.25_dp)
   anorm = anorm*half
   nsq = nsq+1
 end do
 ABI_MALLOC(a,(nb,nb))
 ABI_MALLOC(term,(nb,nb))
 a = kappa/real(2**nsq,dp)
 umat = czero
 term = czero
 do i=1,nb
   umat(i,i) = cone
   term(i,i) = cone
 end do
 do i=1,16
   term = matmul(term,a)/real(i,dp)
   umat = umat+term
 end do
 do is=1,nsq
   umat = matmul(umat,umat)
 end do
 ABI_FREE(a)
 ABI_FREE(term)

end function rdmft_expm
!!***

end module m_rdmft_adam
!!***
