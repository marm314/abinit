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

## Present limitations

Norm-conserving pseudopotentials, closed shells ([[nsppol]] = 1, [[nspinor]] = 1), no band or k-point parallelization
(each process computes everything), and the integrals of the active space are kept in memory as dense
matrices. The BBC2 functional of the literature is not implemented.
