---
description: How to perform Reduced Density Matrix Functional Theory (RDMFT) calculations
authors: MRM
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives some information on how to perform Reduced Density Matrix Functional Theory (RDMFT) calculations
with ABINIT.

!!! warning

    RDMFT is under active development and should thus be used with caution!

## Introduction

In RDMFT the total energy is a functional of the one-body reduced density matrix. In the present implementation
the natural orbitals are the Bloch orbitals of a previous ground-state calculation (DFT or Hartree-Fock), which are kept
fixed, and the occupation numbers $n_{n\mathbf{k}}$ (per spin, between 0 and 1) are optimized under the constraint
$\sum_{n\mathbf{k}} w_\mathbf{k} n_{n\mathbf{k}} = N_e/2$. For each pair of orbitals $p=(n,\mathbf{k})$ and
$q=(m,\mathbf{k}')$ of the full Brillouin zone the Coulomb and exchange integrals

$$
J_{pq}=\langle pq|pq\rangle ,\qquad K_{pq}=\langle pq|qp\rangle
$$

are computed with the Coulomb kernel used by ABINIT (the Hartree term of the ground state for $J$, and
the exact exchange one, [[fock_icutcoul]], for $K$). With time-reversal symmetry (no spin-orbit coupling, collinear
spin) the time-inversion integrals of the PNOF functionals are equal to $K_{pq}$.
The electronic energy per unit cell is

$$
E = \frac{1}{N_k}\sum_p 2 n_p h_{pp} + \frac{1}{N_k^2}\sum_{pq}\Big[2 f_H(n_p,n_q) J_{pq}-f_X(n_p,n_q) K_{pq}\Big]
$$

for the JK-only functionals, with the one-electron part $h=T+V_\mathrm{loc}+V_\mathrm{nl}$ computed with the Hamiltonian
of ABINIT, to which the Ewald and pseudopotential core energies are added. The single-determinant functional
($f_H=f_X=n_pn_q$) is the Hartree-Fock energy: with the occupations of an Hartree-Fock ground state it reproduces the
total energy of the ground-state calculation.

## Functionals

The functional is selected with [[rdmft_functional]]:

* JK-only: SD (Hartree-Fock), MBB (Muller), CA, RHFB, CGA, ML, MLSIC, GU, POWER ([[rdmft_power]]) and MULLER_AS. They
  modify the exchange part of the energy, $f_X(n_p,n_q)$. The occupations are optimized with a Newton-SQP method
  with the exact gradient and Hessian (box $10^{-6}\le n\le 1-10^{-6}$ and the electron-count equality).
* PNOF5, PNOF7, PNOF7s and GNOF (Piris natural orbital functionals). The occupied orbitals are the principal pairs
  of coupling subspaces, each one coupled with ([[rdmft_ncoupl]] - 1) virtual orbitals. The occupations of a
  subspace are parametrized by angles that make its electron count and the bounds automatic, and are optimized with
  L-BFGS with the analytic gradient. Orbitals of different k-points can share a subspace.

With integer occupations all the functionals except MULLER_AS (and ML, MLSIC to $10^{-6}$ Ha) give the
Hartree-Fock energy. For periodic systems the exchange integral of an orbital with itself, $K_{pp}$, differs
from $J_{pp}$ by a finite-size term of the treatment of the divergence of the Coulomb interaction, and the PNOF
energies then differ from the Hartree-Fock one at integer occupations: with [[rdmft_selfint]] = 1 this
term is added and they agree.

## Time reversal and the k-point mesh

The PNOF subspaces must be closed under $\mathbf{k}\rightarrow-\mathbf{k}$. The orbital at $-\mathbf{k}$ is the complex conjugate
of the one at $\mathbf{k}$, $c_{n,-\mathbf{k}}(\mathbf{G})=c^*_{n\mathbf{k}}(-\mathbf{G})$. If the mesh of the ground-state calculation does
not contain $-\mathbf{k}$ for all its points (a Monkhorst-Pack mesh either contains all of them or none), the PNOF
optimization completes it with these orbitals, which doubles the number of k-points. A mesh with only some
of the partners is refused. The JK-only functionals use the mesh as it is: the energy only needs the
integrals among the k-points of the mesh, and the occupations of $\mathbf{k}$ and $-\mathbf{k}$ are
tied if both are present.

The energy functionals only use the integrals $K_{pq}$ (not the time-inversion integrals $L_{pq}$, which equal $K_{pq}$ only
when the orbital at $-\mathbf{k}$ is the conjugate of the one at $\mathbf{k}$ with a fixed phase), and $K_{pq}$ does not depend on the
phase of the orbitals. No alignment of the phases of the time-reversed partners is therefore needed. The ties of the
occupations of $\mathbf{k}$ and $-\mathbf{k}$ do assume that the orbital $n$ at $-\mathbf{k}$ is the conjugate of the orbital $n$ at
$\mathbf{k}$ up to a phase; this is checked at the start (largest $1-|\langle\psi_{n,-\mathbf{k}}|\Theta\psi_{n\mathbf{k}}\rangle|$ is
printed, with a warning if it is larger than $10^{-4}$, for example for a mixture of degenerate states).

## Using time reversal to reduce the k-points

The orbital at $-\mathbf{k}$ is the complex conjugate of the one at $\mathbf{k}$. A ground-state calculation can therefore use a
mesh reduced by time reversal ([[kptopt]] = 2, only one of each pair $\pm\mathbf{k}$), which halves its cost: for all the
functionals the missing $-\mathbf{k}$ states are built by conjugation and the completed mesh gets uniform weights (the points that
are their own $-\mathbf{k}$ are kept once). Then, with [[rdmft_tr_reduce]] = 1 (default), the Coulomb and exchange integrals are
computed only for the independent blocks of k-points and the rest of the matrices is filled by the symmetry relations
$J_{pq}=J_{p\bar q}=J_{\bar p\bar q}$ and $K_{\bar p\bar q}=K_{pq}$; the generalized Fock matrices already use
the representative k-points. [[rdmft_tr_reduce]] = 0 computes all the pairs explicitly, for debugging; the results must be the same.

## Active space and cost

The occupations that are optimized are those of the active orbitals: the occupied orbitals with an energy larger
than E$_\mathrm{HOMO}$ - [[rdmft_ecut_act]] (or the [[rdmft_nocc_act]] highest occupied bands of each k-point) and, for
the JK-only functionals, the virtual orbitals with an energy smaller than E$_\mathrm{LUMO}$ + [[rdmft_ecut_act]] (or
the [[rdmft_nvir_act]] lowest virtual bands). The deeper occupied orbitals are frozen at 1 and the other virtual
orbitals at 0. The frozen virtual orbitals do not contribute to the energy and are dropped. The frozen occupied orbitals
only interact with the active ones through sums of integrals, and the integrals are computed only for the pairs
with an active orbital. The exchange among the frozen occupied orbitals is a constant of the total
energy that needs all their pairs; it is computed if [[rdmft_core_ints]] = 1 (default), and can be skipped
(0) to save time when only energy differences matter.

## Workflow

1. A ground-state calculation (dataset 1) with norm-conserving pseudopotentials, [[nsppol]] = 1, [[nspinor]] = 1,
   [[istwfk]] = 1 at all the k-points, the full Brillouin zone ([[kptopt]] 3 with [[nsym]] 1) and enough bands
   ([[nband]]) for the virtual orbitals that are needed, writing the wavefunctions ([[prtwf]] 1).
   For Hartree-Fock orbitals use [[ixc]] = 40 (a first DFT dataset helps its convergence).
2. An RDMFT dataset with [[optdriver]] = 11 that reads these orbitals with [[getwfk]] (or [[getwfk_filepath]]),
   with the same cell, k-points, [[ecut]] and pseudopotentials, and the input variables of the RDMFT:
   [[rdmft_functional]], the active space and [[rdmft_ncoupl]], [[rdmft_selfint]]...
3. The optimized occupation numbers are printed for the first k-point and written in the file
   `*_RDMFT_OCC`. Copying it to `RDMFT_OCC.in` and setting [[rdmft_restart]] = 1 restarts a calculation from them
   ([[rdmft_start]] = 1 starts the JK-only functionals from Fermi-Dirac occupations instead).

An example (N$_2$ molecule in a cubic cell of 4 angstrom with 8 k-points; the Hartree-Fock orbitals, then the Muller
functional and GNOF with one active occupied band):

```
ndtset 4
acell 3*4 angstrom
ntypat 1  znucl 7  natom 2  typat 1 1
xcart -0.55 0 0  0.55 0 0 angstrom
ecut 8  nband 8  nsym 1  istwfk *1
kptopt 3  ngkpt 2 2 2  nshiftk 1  shiftk 0 0 0
pseudos "N.psp8"  tolvrs 1e-12  nstep 150  prtwf 1

# 1: PBE ground state
ixc1 11
# 2: Hartree-Fock ground state (plain Coulomb exchange), starting from the PBE orbitals
ixc2 40  fock_icutcoul2 3  getwfk2 -1  iscf2 7  diemix2 0.3  diemac2 1.0
# 3: RDMFT, Muller functional
optdriver3 11  ixc3 40  fock_icutcoul3 3  getwfk3 2  rdmft_functional3 2
# 4: RDMFT, GNOF, one active occupied band per k-point, perfect pairing
optdriver4 11  ixc4 40  fock_icutcoul4 3  getwfk4 2  rdmft_functional4 24
rdmft_nocc_act4 1  rdmft_ncoupl4 2  rdmft_selfint4 1
```

## Checking the orbital-rotation gradient numerically

With [[rdmft_check_rot]] the code rotates two states of the same k-point (the rotation of the orbitals can only mix
states of the same k-point, and the time-reversed orbitals at $-\mathbf{k}$ are rotated with the conjugate matrix),
recomputes the energy of the functional at the optimized (fixed) occupations, and prints the central-difference
gradient and curvature with respect to the rotation angle, in the real and the imaginary directions of the
anti-Hermitian generator. For Hartree-Fock orbitals and the SD functional the gradient vanishes (Brillouin theorem), and
it vanishes for any functional for two orbitals with the same occupation. The gradient is complex (one real and one
imaginary direction per pair); only its modulus does not depend on the phases of the orbitals. This is the numerical
reference to validate the analytic orbital gradient, which is printed next to it.

## Analytic orbital gradient

For the rotation of the orbitals the energy is written with coefficient matrices,
$E=\sum_p \nu_p h_{pp}+N_k^{-2}\sum_{pq}[A_{pq}J_{pq}-C_{pq}K_{pq}]$, which covers the JK-only functionals
($A=2f_H$, $C=f_X$) and the PNOF ones, whose self-interaction terms ([[rdmft_selfint]]) are diagonal entries of $A$ and $C$.
The generalized Fock matrix $F_{pj}$ of a k-point is then built with the plane-wave version of the resolution-of-identity
contraction: the auxiliary index is the reciprocal vector $\mathbf{G}$, the Hartree part uses the coupling-weighted density
$W_j(\mathbf{G})=\sum_q a_{jq}\rho_{qq}(\mathbf{G})$ and the pair densities of the states of the same k-point, and the exchange
part the pair densities of all the bands of the k-point with every orbital $q$ of any k-point. $F$ is block diagonal in k:
only states of the same k-point are rotated. The gradient of the pair $(p,q)$ is $2\,\mathrm{Re}(F_{pq}-F_{qp}^*)$ in the
real direction and $-2\,\mathrm{Im}(F_{pq}-F_{qp}^*)$ in the imaginary one.

## Optimization of the orbitals (ADAM)

With [[rdmft_orb_macro]] > 0 the orbitals are optimized too. Each macro-iteration rotates the orbitals at fixed occupations
with the ADAM optimizer (never a Hessian) and then re-optimizes the occupations for the new orbitals, until the energy
changes by less than [[rdmft_orb_etol]]. The orbitals are rotated with $\psi'_j=\sum_i\psi_i U_{ij}$, $U=\exp(\kappa)$,
where $\kappa$ is anti-Hermitian and **block diagonal in k**: only states of the same k-point are mixed. Time reversal is
preserved: the rotation at $-\mathbf{k}$ is the conjugate matrix, so only one of the two k-points is a set of independent
parameters (and the rotations at time-reversal-invariant k-points are real). The parameters are the complex
$\kappa_{pq}$ ($p<q$), and the gradient is the analytic one of the previous section.

The ADAM optimizer is the complex version of the one of ReRDMFT: complex first moment (decay [[rdmft_adam_beta1]]), a real
second moment of $|g|^2$ shared by the real and imaginary parts ([[rdmft_adam_beta2]]), bias corrections and the running
maximum of the second moment. A run has [[rdmft_adam_iter]] iterations (learning rate [[rdmft_adam_lr]]) and stops earlier
if the largest gradient is below [[rdmft_orb_gtol]] or the energy changes by less than [[rdmft_orb_etol]]. If a run does not
lower the energy, it restarts from the best point with the learning rate multiplied by [[rdmft_adam_rate]] and
[[rdmft_adam_incr]] more iterations (up to [[rdmft_adam_maxit]]). The rotations are limited to the bands of the WFK file:
the minimum can only be reached if these bands span the orbitals of the optimum, so a large [[nband]] is needed. A good test
is to start from PBE orbitals with the SD functional: the energy tends to the Hartree-Fock one as [[nband]] grows.

## Orbital space and the WFK file

The orbitals, and therefore the occupations and the rotations of the orbitals, are those of the bands stored in the WFK file
of the ground-state calculation: the [[nband]] of the RDMFT dataset is ignored. The cutoff [[ecut]] (together with the cell, the
k-points and the pseudopotentials) must be the one of the ground-state calculation: it fixes the plane-wave basis in which
the orbitals are expanded, the FFT mesh of the pair densities and the one-electron operator. The number of
states that [[ecut]] can produce, however, is the number of plane waves (hundreds per k-point even for a small cell), and
the WFK file usually contains only a few of them. All the Coulomb and exchange integrals are computed exactly for the
loaded orbitals, but the minimization is restricted to the space they span: the energies are variational upper
bounds of the ones in the complete basis. For example, the single-determinant functional started from Kohn-Sham orbitals
([[rdmft_functional]] = 1 with [[rdmft_orb_macro]] > 0) must reach the Hartree-Fock energy, and it only does if the loaded
bands span the Hartree-Fock orbitals. The output reports the bands loaded and the plane waves available.

With [[rdmft_complete]] = 1 the band set of every k-point is completed to all the plane waves (a Gram-Schmidt procedure on plane
waves). Then the rotations span the whole basis, and PBE orbitals optimized with the SD functional reproduce the Hartree-Fock
total energy of the ground-state calculation (to about 1e-9 Ha, with [[rdmft_occ_min]] very small). This is expensive
(the number of orbitals is the number of plane waves times the number of k-points) and is intended for small cells and
tests. The cost is greatly reduced by [[rdmft_occ_skip]], which drops the orbitals that are (almost) unoccupied from the
sums of pair densities.

## Degenerate states and the choice of the orbitals

The energy of a functional with fractional occupations (MBB, CGA, MLSIC, GU, ...) at **fixed** orbitals is not invariant under
unitary rotations of degenerate orbitals that have different occupations: the Coulomb and exchange integrals of the rotated
orbitals are different. At time-reversal-invariant k-points and along symmetry lines the Kohn-Sham or Hartree-Fock orbitals are
degenerate, and the basis returned by the diagonalization within a degenerate shell is arbitrary (it depends on the solver and on
the starting vectors). The fixed-orbital energies of these functionals then depend on that arbitrary choice: in N2 two ground
states that only differ by the eigensolver gave MBB and CGA energies that differ by 5e-5 Ha, while the single-determinant
functional (invariant under such rotations) and the PNOF with one active occupied band do not change. This also applies to the
bands that are read in a degenerate shell at the time-reversal-invariant points, which are replaced by real orbitals (needed
for the real rotations at those points). The dependence disappears when the orbitals are optimized ([[rdmft_orb_macro]] > 0), and it is
absent for non-degenerate states. When comparing fixed-orbital energies of such functionals, use the same ground-state run.

## Present limitations

Norm-conserving pseudopotentials, closed shells ([[nsppol]] = 1, [[nspinor]] = 1), no band or k-point parallelization
(each process computes everything), and the integrals of the active space are kept in memory as dense
matrices. The BBC2 functional of the literature is not implemented.
