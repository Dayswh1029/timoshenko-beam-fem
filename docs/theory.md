# Theory notes: the Timoshenko beam element used in this library

## 1. Why Timoshenko instead of Euler-Bernoulli

The Euler-Bernoulli beam theory assumes plane sections remain
perpendicular to the neutral axis after deformation (i.e. shear
deformation is neglected). This is a good approximation for slender
beams, but becomes inaccurate as the length-to-depth ratio decreases —
which matters for short, stubby beam segments (e.g. a stepped shaft
transition, or a stiff support spring modeled as an equivalent beam
element, as in the dynamic vibration absorber problem this library
originated from).

The Timoshenko beam theory adds an independent rotation degree of
freedom for the cross-section, decoupled from the slope of the
deflection curve, and introduces a shear correction factor to account
for the non-uniform shear stress distribution across the section.

## 2. Non-dimensional shear parameter

For a two-node beam element of length `L`, the ratio between shear and
bending flexibility is captured by the non-dimensional parameter

```
Phi = chi * 12 * E * Iy / (G * A * L^2)
```

where:
- `E` = Young's modulus
- `G` = shear modulus
- `A` = cross-sectional area
- `Iy` = second moment of area about the bending axis
- `chi` = shear correction factor (section-shape dependent)

Setting `Phi = 0` recovers the classical Euler-Bernoulli beam element
exactly — this library's Timoshenko formulation is a strict
generalization, not a separate code path.

Common shear correction factors:
- Rectangular section: `chi = (12 + 11*nu) / (10*(1+nu))`
- Circular section: `chi = 6*(1+nu) / (7 + 6*nu)`

## 3. Element degrees of freedom

Each node has 2 planar DOF:
- `u`     transverse displacement
- `theta` bending rotation

Element DOF order: `[u1, theta1, u2, theta2]`.

## 4. Stiffness matrix

```
Ke = E*Iy / (L^3*(1+Phi)) *
     [ 12        6L         -12        6L       ;
       6L        (4+Phi)L^2 -6L        (2-Phi)L^2;
      -12       -6L          12       -6L       ;
       6L        (2-Phi)L^2 -6L        (4+Phi)L^2]
```

## 5. Consistent mass matrix

The mass matrix is the sum of a translational-inertia contribution
(scales with `rho*A*L`) and a rotary-inertia contribution (scales with
`rho*Iy/L`), both functions of `Phi`:

```
Me = M_translational(Phi) + M_rotary(Phi)
```

See `src/timoshenkoElement.m` for the full closed-form expressions
(coefficients `m1..m10`), which follow the standard consistent-mass
Timoshenko beam formulation (Friedman & Kosmatka, "An improved two-node
Timoshenko beam finite element", Computers & Structures, 1993; see also
Przemieniecki, "Theory of Matrix Structural Analysis", 1968, for the
Phi=0 Euler-Bernoulli limit).

## 6. Assembly

Global matrices are built by looping over elements, mapping each
element's local 4x4 matrices into the appropriate 4x4 block of the
global matrix via a **DOF map** (node `i` -> global DOF `[2i-1, 2i]`).
This is the same idea as the "map matrix" used in the original
coursework report, generalized so that segment boundaries, element
counts, and node numbering are computed automatically instead of
hand-tabulated.

## 7. Boundary conditions

Fixed DOF are removed by deleting the corresponding rows/columns from
the global `M` and `K` matrices (static condensation of constrained
DOF). Unlike hand-picking a contiguous DOF range (`M(3:end, 3:end)`,
valid only when the fixed node happens to be node 1), this library
accepts an arbitrary list of DOF indices, so any combination of
clamped/pinned/free supports at any location can be modeled.

## 8. Modal damping

Because Timoshenko/Euler-Bernoulli FE models have no first-principles
damping mechanism, damping is introduced empirically as a per-mode
damping ratio (identified, in the original coursework, by tuning
against a measured frequency response function). Given a set of modal
damping ratios `zeta_i`, the physical damping matrix is reconstructed
by:

1. Solving the undamped eigenproblem to get mode shapes `phi`
2. Computing modal mass/stiffness `Mm = phi'*M*phi`, `Km = phi'*K*phi`
3. Building the diagonal modal damping matrix
   `Cm(i,i) = 2*zeta_i*sqrt(Km(i,i)*Mm(i,i))`
4. Transforming back to physical coordinates:
   `C = phi' \ Cm / phi`

## 9. State-space form and frequency response

For harmonic/transient response analysis, the constrained `{M, C, K}`
system is converted to first-order state-space form with state vector
`z = [qdot; q]`:

```
zdot = [ -M^-1*C   -M^-1*K ] z + [ M^-1*e_input ] f(t)
       [    I          0   ]     [      0       ]

y = e_output' * q
```

where `e_input`/`e_output` are unit vectors selecting the DOF where the
force is applied / the response is measured. This lets `bode()` /
`freqresp()` (MATLAB Control System Toolbox) compute the frequency
response function directly, without hand-deriving a transfer function
for each new problem.

## 10. Validation

See `tests/test_cantilever_analytical.m`, which compares the FE natural
frequencies of a slender clamped-free beam against the closed-form
Euler-Bernoulli solution

```
f_i = (beta_i^2 / (2*pi*L^2)) * sqrt(E*I / (rho*A)),  beta = [1.875, 4.694, 7.855, ...]
```

The FE result should converge to the analytical one as the shear
parameter `Phi -> 0` (slender beam) and the mesh is refined; the
included test checks agreement to within 1% even at a coarse (10
element) mesh.
