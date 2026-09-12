function [Mc, Kc, freeDofs] = applyBoundaryConditions(M, K, fixedDofs)
% APPLYBOUNDARYCONDITIONS  Remove fixed DOFs from global M and K.
%
%   [Mc, Kc, freeDofs] = applyBoundaryConditions(M, K, fixedDofs)
%
%   Inputs
%     M, K         Global (unconstrained) mass and stiffness matrices
%     fixedDofs    Vector of global DOF indices to constrain (e.g. for a
%                  cantilever fixed at node 1, fixedDofs = [1, 2] removes
%                  the translational and rotational DOF of node 1).
%
%   Outputs
%     Mc, Kc       Constrained (reduced) mass/stiffness matrices
%     freeDofs     Indices (into the original numbering) of the DOFs
%                  that remain -- useful for mapping a reduced-space
%                  result (e.g. mode shapes) back onto the full mesh.
%
%   This replaces the "Mc = M(3:nf,3:nf)" style of hand-picking a
%   contiguous DOF range: fixedDofs can be any subset (e.g. clamped at
%   one end, simply supported at another, or a support with only the
%   translational DOF constrained) and the function still works.

    nDof = size(M, 1);
    freeDofs = setdiff(1:nDof, fixedDofs);

    Mc = M(freeDofs, freeDofs);
    Kc = K(freeDofs, freeDofs);
end
