function sys = beamStateSpace(M, C, K, inputDof, outputDof)
% BEAMSTATESPACE  Convert a constrained {M, C, K} beam model to a
%   state-space system driven by a point force, with a chosen output
%   DOF, ready for frequency-response (FRF) computation.
%
%   sys = beamStateSpace(M, C, K, inputDof, outputDof)
%
%   Inputs
%     M, C, K      Constrained (reduced) mass, damping and stiffness
%                  matrices (same DOF numbering)
%     inputDof     Index (into the reduced DOF numbering) where a point
%                  force is applied
%     outputDof    Index (into the reduced DOF numbering) of the
%                  response to output (e.g. an accelerometer location)
%
%   Output
%     sys          MATLAB Control System Toolbox ss() object relating
%                  the input force to the output displacement.
%
%   State vector: z = [qdot; q], so that
%       zdot = [ -M\C   -M\K ; I  0 ] z + [ M\ei ; 0 ] f(t)
%       y    = ej' * q
%   where ei/ej are unit vectors selecting the input/output DOF.
%
%   This is a direct generalization of the ad hoc
%       AA = [-inv(Mc)*C_damp, -inv(Mc)*Kc; eye(size(Mc)), zeros(size(Mc))];
%       B  = [inv(Mc); zeros(size(Mc))];
%       C  = AA(2*n_acc-1-2, :);
%       D  = B(2*n_acc-1-2, :);
%   pattern: instead of hand-computing a row index into AA/B (fragile --
%   breaks the moment the mesh or DOF numbering changes), the caller
%   just names the input/output DOF directly and this function builds
%   the correct rows regardless of model size.

    n = size(M, 1);
    ei = zeros(n, 1); ei(inputDof) = 1;
    ej = zeros(n, 1); ej(outputDof) = 1;

    Minv = M \ eye(n);

    Adyn = [-Minv*C, -Minv*K; eye(n), zeros(n)];
    Bdyn = [Minv*ei; zeros(n,1)];
    Cdyn = [zeros(1,n), ej'];   % output = displacement, not velocity
    Ddyn = 0;

    sys = ss(Adyn, Bdyn, Cdyn, Ddyn);
end
