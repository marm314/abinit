---
description: How to perform Reduced Density Matrix Functional Theory (RDMFT) calculations
authors: MRM
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives some information on the Reduced Density Matrix Functional Theory (RDMFT) calculations of ABINIT.

!!! warning

    RDMFT is under active development and should thus be used with caution!

## Introduction

In RDMFT the total energy is a functional of the one-body reduced density matrix. For a fixed set of orbitals
(the natural orbitals), the energy is a function of the occupation numbers, which are optimized under the constraints
$0 \le n_p \le 1$ and $\sum_p n_p = N_e/2$.

## Implementation in ABINIT

The RDMFT driver is selected with [[optdriver]] = 11. It starts from the orbitals of a ground-state DFT or Hartree-Fock
calculation (read with [[getwfk]], using [[istwfk]] = 1 at all k-points) and works in the full Brillouin zone.
For the moment it computes the Coulomb and exchange integrals between all the orbitals, the one-electron part of the
energy, and the energy of the JK-only and PNOF functionals for the ground-state occupations.
With integer occupations the single-determinant functional gives the Hartree-Fock energy.
