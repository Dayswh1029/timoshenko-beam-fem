function C = modalDampingMatrix(M, K, zeta)
% MODALDAMPINGMATRIX  Build a physical damping matrix from a set of
%   per-mode modal damping ratios.
%
%   C = modalDampingMatrix(M, K, zeta)
%
%   Inputs
%     M, K    Constrained (reduced) mass and stiffness matrices
%     zeta    Modal damping ratio, applied uniformly to every mode
%             (scalar) or specified individually as a vector
%             (one entry per mode; must be as long as size(M,1)).
%
%   Output
%     C       Physical-coordinate damping matrix, same size as M and K,
%             obtained by:
%               1. solving the undamped eigenproblem (K - lambda*M)phi=0
%               2. computing modal mass/stiffness: Mm = phi'*M*phi,
%                  Km = phi'*K*phi
%               3. building the modal damping matrix
%                  Cm(i,i) = 2*zeta(i)*sqrt(Km(i,i)*Mm(i,i))
%               4. transforming back to physical coordinates:
%                  C = phi' \ Cm / phi   (equivalent to
%                  C = inv(phi)' * Cm * inv(phi))
%
%   This generalizes the fixed damping-ratio "tuning" step used to match
%   an experimental FRF (e.g. zeta = 0.01 for an undamped cantilever, or
%   a different value per configuration) into a reusable utility that
%   accepts either one damping ratio for all modes or a full vector.

    n = size(M, 1);
    [phi, lambda] = eig(K, M);
    lambda = diag(lambda); %#ok<NASGU>  % (kept for clarity / future use)

    if isscalar(zeta)
        zeta = zeta * ones(n, 1);
    end
    assert(numel(zeta) == n, ...
        'zeta must be a scalar or have one entry per mode (%d).', n);

    Mm = phi' * M * phi;
    Km = phi' * K * phi;

    Cm = zeros(n);
    for i = 1:n
        Cm(i,i) = 2 * zeta(i) * sqrt(Km(i,i) * Mm(i,i));
    end

    C = phi' \ Cm / phi;
end
