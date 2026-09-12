# Timoshenko Beam FEM

A small, general-purpose MATLAB library for static/dynamic finite element
analysis of planar (2-D) beams using **Timoshenko beam elements** (shear
deformation included; reduces to Euler-Bernoulli when shear stiffness is
large).

The library grew out of a university coursework project (Politecnico di
Torino, *Dynamic Design of Machines*, 2018–2019) that identified the
dynamic behavior of a cantilever beam fitted with a tuned dynamic
vibration absorber, using experimental modal testing (impact hammer +
accelerometer) to validate a hand-built FE model.

The **original** coursework scripts worked, but were written as three
near-duplicate files (one per test configuration) with hand-picked
matrix indices and a fixed element count baked into the code. This
repository is a rewrite of the same underlying element formulation and
assembly logic into a small, reusable, **data-driven** library: change
a segment's length, cross-section, material, mesh density, or boundary
conditions by changing *input arguments*, not by editing matrix index
arithmetic.

## What this is / is not

- **Is**: a from-scratch implementation of the 2-D Timoshenko beam
  element (consistent mass and stiffness matrices), a generic
  multi-segment mesh assembler, boundary condition handling, modal
  damping construction, and a state-space / FRF utility — all unit
  tested against a closed-form analytical solution.
- **Is not**: a general-purpose commercial-grade FE package. It handles
  planar bending of prismatic beam segments only (no torsion/axial
  coupling, no 3-D frames, no gyroscopic/rotordynamic terms, no
  nonlinearities). Scope is intentionally narrow and matches what the
  original coursework problem needed.

## Repository layout

```
timoshenko-beam-fem/
├── src/                          Core library (no hard-coded problem data)
│   ├── timoshenkoElement.m       Single-element mass & stiffness matrices
│   ├── assembleBeamFE.m          Generic multi-segment mesh + assembly
│   ├── addLumpedMass.m           Add point mass/inertia at any node
│   ├── applyBoundaryConditions.m Remove any subset of DOF (not just a range)
│   ├── modalDampingMatrix.m      Build C from per-mode damping ratios
│   └── beamStateSpace.m          Constrained {M,C,K} -> state-space / FRF
├── examples/
│   ├── example1_cantilever_absorber.m   Reproduces the original 3-configuration
│   │                                     coursework problem using the generic library
│   └── example2_generality_demo.m       Mesh refinement, different BCs, and a
│                                         stepped (multi-segment) shaft — all with
│                                         zero changes to src/
├── tests/
│   └── test_cantilever_analytical.m     Validates FE results against the
│                                         closed-form Euler-Bernoulli cantilever
│                                         solution (< 0.1% error)
└── docs/
    └── theory.md                        Element derivation notes and references
```

## Quick start

```matlab
addpath('src');

% Define one prismatic segment (SI units)
seg.L   = 0.28;                        % length [m]
seg.E   = 70e9;                        % Young's modulus [Pa]
seg.nu  = 0.33;
seg.G   = seg.E / (2*(1+seg.nu));      % shear modulus [Pa]
seg.rho = 2700;                        % density [kg/m^3]
seg.A   = 0.03 * 0.003;                % cross-section area [m^2]
seg.Iy  = 0.03 * 0.003^3 / 12;         % second moment of area [m^4]
seg.chi = (12+11*seg.nu) / (10*(1+seg.nu)); % rectangular-section shear factor

% Mesh: 28 elements along the segment
model = assembleBeamFE(seg, 28);

% Clamp node 1 (cantilever boundary condition)
[Mc, Kc] = applyBoundaryConditions(model.M, model.K, [1, 2]);

% Natural frequencies
wn = sort(sqrt(eig(Kc, Mc))) / (2*pi);
disp(wn(1:3))   % first three natural frequencies, in Hz
```

See `examples/` for a full worked problem (beam + tuned vibration
absorber, matched against experimental data) and a demonstration of
mesh refinement, alternate boundary conditions, and a stepped shaft —
all using the same unmodified library code.

## Validation

`tests/test_cantilever_analytical.m` checks the FE natural frequencies
of a slender cantilever against the closed-form Euler-Bernoulli
solution

```
f_i = (beta_i^2 / (2*pi*L^2)) * sqrt(E*I / (rho*A))
```

with `beta = [1.875, 4.694, 7.855]` for the first three clamped-free
modes. The FE result converges to the analytical one to within ~0.01%
even at a coarse mesh (10 elements), confirming the element formulation
and assembly are implemented correctly.

The library was also checked against the first natural frequency
reported in the original coursework's experimental configuration
(bare aluminum cantilever, tuned Young's modulus): the generic
assembler reproduces the ~30 Hz first resonance visible in that
report's measured/simulated FRF overlay.

## Background / provenance

This library reimplements, in a reusable form, the Timoshenko beam
element formulation (consistent mass matrix per Friedman & Kosmatka,
1993, and the classical Przemieniecki stiffness form) that was
originally hand-derived and coded (in a problem-specific, non-reusable
form) as part of a graduate course project:

> *Dynamic Design of Machines* — Politecnico di Torino, M.Sc. Mechanical
> Engineering, A.Y. 2018–2019. Project: "Identification of vibrating
> structures, effects of dynamic dampers" — flexural dynamic behavior of
> a cantilever beam equipped with a tuned dynamic vibration absorber,
> validated against impact-hammer experimental modal testing.

The original submission solved three specific configurations (bare
beam, beam + rigid end mass, beam + tuned absorber) with duplicated,
hard-coded scripts. This repository factors the same physics into a
tested, general-purpose library, and uses it to reproduce all three
original configurations as a worked example (`examples/example1_*.m`).

## Requirements

- MATLAB (Control System Toolbox required only for `beamStateSpace.m` /
  the FRF plots in `example1`; the core assembly, eigenvalue and
  boundary-condition functions have no toolbox dependencies).

## License

MIT — see `LICENSE`.
