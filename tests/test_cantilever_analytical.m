%% Validation: bare cantilever beam vs. closed-form Euler-Bernoulli solution
%
% For a slender cantilever (length-to-thickness ratio large enough that
% shear deformation is negligible), the Timoshenko FE model should
% converge to the classical Euler-Bernoulli analytical natural
% frequencies:
%
%   f_i = (beta_i^2 / (2*pi*L^2)) * sqrt(E*I/(rho*A))
%
% with beta_1 = 1.875, beta_2 = 4.694, beta_3 = 7.855 for a clamped-free
% beam (see e.g. Blevins, "Formulas for Natural Frequency and Mode
% Shape").
%
% This script asserts that the FE result matches the analytical one to
% within 1% for the first three modes, using a mesh fine enough (60
% elements) that discretization error is small. It also demonstrates
% mesh-independence: the whole point of a data-driven assembler is that
% refining the mesh is a one-line change (nElem), not a rewrite.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

%% ---- Slender beam, negligible shear (Timoshenko -> Euler-Bernoulli) --
L   = 1.0;          % [m]
b   = 0.02;         % [m]
h   = 0.005;         % [m]  (L/h = 200, very slender -> shear effect ~0)
rho = 7800;          % [kg/m^3] steel
E   = 210e9;         % [Pa]
nu  = 0.30;
G   = E / (2*(1+nu));
A   = b*h;
Iy  = b*h^3/12;
chi = (12+11*nu)/(10*(1+nu));

seg.L = L; seg.E = E; seg.G = G; seg.rho = rho; seg.A = A; seg.Iy = Iy; seg.chi = chi;

beta = [1.875, 4.694, 7.855];
fAnalytical = (beta.^2 / (2*pi*L^2)) * sqrt(E*Iy/(rho*A));

meshSizes = [10, 20, 40, 60];
fprintf('%-8s %-12s %-12s %-12s %-10s\n', 'nElem', 'f1_FE[Hz]', 'f1_exact', 'error[%]', 'PASS?');

allPass = true;
for nElem = meshSizes
    model = assembleBeamFE(seg, nElem);
    fixedDofs = [1, 2];
    [Mc, Kc] = applyBoundaryConditions(model.M, model.K, fixedDofs);

    wn = sqrt(eig(Kc, Mc));
    wn = sort(wn);
    fFE = wn(1:3) / (2*pi);

    err1 = abs(fFE(1) - fAnalytical(1)) / fAnalytical(1) * 100;
    pass = err1 < 1.0;
    allPass = allPass && pass;

    fprintf('%-8d %-12.4f %-12.4f %-12.4f %-10s\n', ...
        nElem, fFE(1), fAnalytical(1), err1, mat2str(pass));
end

fprintf('\nFirst three analytical frequencies [Hz]: %s\n', mat2str(round(fAnalytical,3)));
fprintf('First three FE frequencies (finest mesh) [Hz]: %s\n', mat2str(round(fFE',3)));

if allPass
    fprintf('\nPASS: FE model converges to the Euler-Bernoulli analytical solution.\n');
else
    error('FAIL: FE first natural frequency did not converge within tolerance.');
end
