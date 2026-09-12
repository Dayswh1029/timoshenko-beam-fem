%% Example 2: Demonstrating generality across geometry, BCs and mesh
%
% The point of a data-driven assembler is that none of the following
% require touching assembleBeamFE.m / timoshenkoElement.m:
%
%   (a) changing the number of elements (mesh refinement study)
%   (b) changing boundary conditions (cantilever vs. simply supported)
%   (c) changing the number of segments (a stepped shaft, e.g. two
%       diameters bonded together)
%
% This is the key difference from the original coursework scripts, where
% the element count (28), the node used for the accelerometer, and the
% constrained DOF range were hard-coded into three near-duplicate files.

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

E = 210e9; rho = 7800; nu = 0.30; G = E/(2*(1+nu));
b = 0.02; h = 0.01; A = b*h; Iy = b*h^3/12;
chi = (12+11*nu)/(10*(1+nu));

seg.L = 1.0; seg.E = E; seg.G = G; seg.rho = rho; seg.A = A; seg.Iy = Iy; seg.chi = chi;

%% (a) Mesh refinement: same physical beam, increasing element count
fprintf('--- (a) Mesh refinement study (cantilever) ---\n');
for nElem = [5, 10, 20, 40]
    model = assembleBeamFE(seg, nElem);
    [Mc, Kc] = applyBoundaryConditions(model.M, model.K, [1, 2]);
    wn = sort(sqrt(eig(Kc, Mc)));
    fprintf('  nElem=%-4d  f1=%.4f Hz  f2=%.4f Hz\n', nElem, wn(1)/2/pi, wn(2)/2/pi);
end

%% (b) Boundary conditions: cantilever vs. simply supported
fprintf('\n--- (b) Same mesh, different boundary conditions ---\n');
model = assembleBeamFE(seg, 30);

% Cantilever: clamp node 1 (both DOF)
[Mc, Kc] = applyBoundaryConditions(model.M, model.K, [1, 2]);
wn_cant = sort(sqrt(eig(Kc, Mc)));
fprintf('  Cantilever   f1=%.4f Hz\n', wn_cant(1)/2/pi);

% Simply supported: constrain translation only at both ends
uFirst = model.dofMap(1, 1);
uLast  = model.dofMap(model.nNodes, 1);
[Mc2, Kc2] = applyBoundaryConditions(model.M, model.K, [uFirst, uLast]);
wn_ss = sort(sqrt(eig(Kc2, Mc2)));
fprintf('  Simply-supp. f1=%.4f Hz\n', wn_ss(1)/2/pi);

%% (c) Stepped shaft: two segments with different diameters, one call
fprintf('\n--- (c) Stepped shaft (two segments, different sections) ---\n');
segA.L = 0.5; segA.E = E; segA.G = G; segA.rho = rho; segA.chi = chi;
segA.A = pi*0.02^2/4;      segA.Iy = pi*0.02^4/64;   % 20 mm diameter

segB.L = 0.5; segB.E = E; segB.G = G; segB.rho = rho; segB.chi = chi;
segB.A = pi*0.03^2/4;      segB.Iy = pi*0.03^4/64;   % 30 mm diameter

stepped = [segA, segB];
modelStep = assembleBeamFE(stepped, [15, 15]);  % 15 elements per segment
[McS, KcS] = applyBoundaryConditions(modelStep.M, modelStep.K, [1, 2]);
wnS = sort(sqrt(eig(KcS, McS)));
fprintf('  Stepped shaft (clamped-free)  f1=%.4f Hz  f2=%.4f Hz\n', ...
    wnS(1)/2/pi, wnS(2)/2/pi);

fprintf('\nAll three studies used the SAME assembleBeamFE.m / timoshenkoElement.m\n');
fprintf('with no changes to the library code -- only the input arguments changed.\n');
