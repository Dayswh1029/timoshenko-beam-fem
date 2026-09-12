%% Example 1: Cantilever beam with a tuned dynamic vibration absorber
%
% Reproduces (using the generic library in ../src) the three
% configurations of the "Identification of vibrating structures, effects
% of dynamic dampers" coursework problem:
%   Configuration 1: bare cantilever beam
%   Configuration 2: beam + non-suspended (rigidly attached) end mass
%   Configuration 3: beam + tuned dynamic vibration absorber (spring +
%                    suspended mass + eddy-current damping)
%
% Unlike the original coursework script (three near-duplicate scripts,
% one per configuration, with hand-picked matrix indices), this example
% builds all three configurations by calling the SAME generic functions
% with different segment/mesh/mass arguments -- change nElem, beam
% dimensions, or absorber properties in one place and everything
% downstream (mesh, matrices, state space, FRF) follows automatically.

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

%% ---- Beam geometry & material (AISI 2024 aluminum) -------------------
beam.L   = 280e-3;          % [m] free length
beam.b   = 30e-3;           % [m] width
beam.h   = 3e-3;            % [m] thickness
beam.rho = 2700;            % [kg/m^3]
beam.nu  = 0.33;

nElemBeam = 28;             % same discretization density as coursework

%% ---- Accelerometer & hammer (measurement setup) ----------------------
m_acc   = 18e-3;            % [kg]
L_acc   = 100e-3;           % [m] from clamp
L_ham   = 100e-3;           % [m] from clamp (same node in this setup)

%% ---- Dynamic absorber (Configuration 3 only) --------------------------
absorber.ms   = 15.3e-3;      % [kg] suspended mass (magnets + circuit)
absorber.Jgs  = 923e-9;       % [kg*m^2] its moment of inertia about its own CG
absorber.zs   = -1.67e-3;     % [m] offset of CG from beam free end
absorber.mns  = 19.9e-3;      % [kg] non-suspended (rigidly attached) mass
absorber.Jyns_cg = 2323.6e-9; % [kg*m^2] its moment of inertia about its own CG
absorber.d_ns = 18.01e-3;     % [m] offset of its CG from beam free end

absorber.Ls  = 28e-3;  absorber.bs = 6e-3;  absorber.hs = 0.2e-3;
absorber.rho_s = 7800; absorber.nu_s = 0.30;
absorber.c   = 0.5;         % [N*s/m]   translational eddy-current damping
absorber.ctheta = 1e-4;     % [N*s/rad] rotational eddy-current damping

fprintf('=== Configuration 1: bare cantilever beam ===\n');
res1 = solveConfiguration(1, beam, nElemBeam, m_acc, L_acc, L_ham, absorber, 1.0, 0.01);
report(res1);

fprintf('\n=== Configuration 2: beam + non-suspended end mass ===\n');
res2 = solveConfiguration(2, beam, nElemBeam, m_acc, L_acc, L_ham, absorber, 0.77, 0.006);
report(res2);

fprintf('\n=== Configuration 3: beam + dynamic vibration absorber ===\n');
res3 = solveConfiguration(3, beam, nElemBeam, m_acc, L_acc, L_ham, absorber, 0.77, 0.001);
report(res3);

%% ---- Overlay FRF plot --------------------------------------------------
figure; hold on; grid on;
freqHz = logspace(0.5, 2.2, 4000);
plotFRF(res1.sys, freqHz, 'k-',  'Config 1: bare beam');
plotFRF(res2.sys, freqHz, 'b--', 'Config 2: + end mass');
plotFRF(res3.sys, freqHz, 'r-.', 'Config 3: + absorber');
xlabel('Frequency [Hz]'); ylabel('|Displacement/Force| [m/N]');
title('Cantilever beam FRF at accelerometer location');
legend show; set(gca, 'XScale', 'log', 'YScale', 'log');

%% ======================================================================
function res = solveConfiguration(configId, beam, nElemBeam, m_acc, L_acc, L_ham, absorber, E_tuning, zeta)
% Build, assemble, constrain and (for config 3) add the absorber branch,
% then convert to state space at the accelerometer DOF.

    E0 = 70e9 * E_tuning;   % tuned Young's modulus (matches experimental FRF)
    G0 = E0 / (2*(1+beam.nu));
    A0 = beam.b * beam.h;
    Iy0 = beam.b * beam.h^3 / 12;
    chi0 = (12 + 11*beam.nu) / (10*(1+beam.nu));

    seg.L = beam.L; seg.E = E0; seg.G = G0; seg.rho = beam.rho;
    seg.A = A0; seg.Iy = Iy0; seg.chi = chi0;

    model = assembleBeamFE(seg, nElemBeam);

    % Node nearest to the accelerometer / hammer location
    accNode = round(L_acc / beam.L * nElemBeam) + 1;

    M = model.M;
    M = addLumpedMass(M, model.dofMap, accNode, m_acc, 0);

    if configId == 2
        % Rigidly attached non-suspended mass at the FREE end (last node)
        endNode = model.nNodes;
        M = addLumpedMass(M, model.dofMap, endNode, absorber.mns, ...
                           absorber.Jyns_cg + absorber.mns*absorber.d_ns^2);
    end

    K = model.K;
    Cdd = [];
    if configId == 3
        % Absorber modeled as an extra 2-node "spring" element appended
        % after the beam's free end: reuse the SAME generic element
        % routine with the absorber's own section properties.
        Gs = 78.5e9 * 0.8;                     % tuned shear modulus [Pa]
        Es = 2 * Gs * (1 + absorber.nu_s);      % consistent Young's modulus
        As = absorber.bs * absorber.hs;
        Iys = absorber.bs * absorber.hs^3 / 12;
        chis = (12+11*absorber.nu_s)/(10*(1+absorber.nu_s));

        segAbs.L = absorber.Ls; segAbs.E = Es; segAbs.G = Gs;
        segAbs.rho = absorber.rho_s; segAbs.A = As; segAbs.Iy = Iys;
        segAbs.chi = chis;

        [Me_abs, Ke_abs] = timoshenkoElement(segAbs.E, segAbs.G, segAbs.rho, ...
                                              segAbs.A, segAbs.Iy, segAbs.L, segAbs.chi);
        Me_abs = 2*Me_abs; Ke_abs = 2*Ke_abs;  % two parallel spring legs

        % Extend M, K by one extra node (the absorber's suspended mass)
        nOld = model.nDofFull;
        M(nOld+1, nOld+1) = 0; K(nOld+1, nOld+1) = 0; % pad
        M(nOld+2, nOld+2) = 0; K(nOld+2, nOld+2) = 0;

        endNodeDofs = model.dofMap(model.nNodes, :);
        newDofs     = [nOld+1, nOld+2];
        gdof = [endNodeDofs, newDofs];

        M(gdof, gdof) = M(gdof, gdof) + Me_abs;
        K(gdof, gdof) = K(gdof, gdof) + Ke_abs;

        M(newDofs(1), newDofs(1)) = M(newDofs(1), newDofs(1)) + absorber.ms;
        M(newDofs(2), newDofs(2)) = M(newDofs(2), newDofs(2)) + ...
            absorber.Jgs + absorber.ms*absorber.zs^2;

        Cdd = zeros(4);
        Cdd([1 3],[1 3]) = absorber.c     * [1 -1; -1 1];
        Cdd([2 4],[2 4]) = absorber.ctheta * [1 -1; -1 1];
        cddDofs = gdof; % [end u, end theta, new u, new theta]
    end

    fixedDofs = [1, 2]; % clamp at node 1 (first two global DOF)
    [Mc, Kc, freeDofs] = applyBoundaryConditions(M, K, fixedDofs);

    C = modalDampingMatrix(Mc, Kc, zeta);

    if configId == 3
        localIdx = arrayfun(@(d) find(freeDofs == d), cddDofs);
        C(localIdx, localIdx) = C(localIdx, localIdx) + Cdd;
    end

    accDofU = model.dofMap(accNode, 1);
    inputLocal  = find(freeDofs == accDofU);
    outputLocal = inputLocal; % collocated force/response, as in coursework

    sys = beamStateSpace(Mc, C, Kc, inputLocal, outputLocal);

    [V, D] = eig(Kc, Mc);
    wn = sqrt(diag(D));
    [wn, order] = sort(wn);
    V = V(:, order);

    res.configId = configId;
    res.model    = model;
    res.sys      = sys;
    res.naturalFreqHz = wn(1:min(3,end)) / (2*pi);
end

function report(res)
    fprintf('  First natural frequencies [Hz]: %s\n', ...
        mat2str(round(res.naturalFreqHz, 2)));
end

function plotFRF(sys, freqHz, style, labelStr)
    [mag, ~] = bode(sys, 2*pi*freqHz);
    mag = squeeze(mag);
    plot(freqHz, mag, style, 'LineWidth', 1.5, 'DisplayName', labelStr);
end
