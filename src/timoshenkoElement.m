function [Me, Ke] = timoshenkoElement(E, G, rho, A, Iy, L, shearCorrection)
% TIMOSHENKOELEMENT  Consistent mass and stiffness matrices for a 2-D
%   Timoshenko beam element (planar bending only).
%
%   [Me, Ke] = timoshenkoElement(E, G, rho, A, Iy, L, shearCorrection)
%
%   Each node has 2 degrees of freedom: transverse displacement (u) and
%   bending rotation (theta). Element DOF order is:
%       [u1, theta1, u2, theta2]
%
%   Inputs
%     E                 Young's modulus                [Pa]
%     G                 Shear modulus                  [Pa]
%     rho               Material density                [kg/m^3]
%     A                 Cross-section area              [m^2]
%     Iy                Second moment of area about the bending axis [m^4]
%     L                 Element length                  [m]
%     shearCorrection   Timoshenko shear correction factor (chi), e.g.
%                        chi = (12+11*nu)/(10*(1+nu)) for a rectangular
%                        section, or 6*(1+nu)/(7+6*nu) for a circular one.
%
%   Outputs
%     Me   4x4 consistent mass matrix (element local frame)
%     Ke   4x4 stiffness matrix (element local frame)
%
%   The formulation follows the standard consistent Timoshenko beam
%   element (e.g. Friedman & Kosmatka 1993; Przemieniecki 1968), where
%   the non-dimensional shear parameter is
%
%       Phi = shearCorrection * 12 * E * Iy / (G * A * L^2)
%
%   and Phi = 0 recovers the classical Euler-Bernoulli beam element.
%
%   This function is purely geometric/material -> element matrices. It
%   has no knowledge of global numbering, boundary conditions, or how
%   many elements a beam is discretized into: those are handled by
%   assembleBeamFE.m, which can build a mesh of ANY length and ANY
%   number of elements from a simple table of segment properties.
%
%   Reference formulation cross-checked against the derivation used in
%   Politecnico di Torino's "Dynamic Design of Machines" course project
%   (cantilever beam + dynamic vibration absorber identification).

    Phi = shearCorrection * 12 * E * Iy / (G * A * L^2);

    % --- Mass matrix -----------------------------------------------
    % Translational (bending) inertia contribution
    m1  = 156 + 294*Phi + 140*Phi^2;
    m2  = 22  + 38.5*Phi + 17.5*Phi^2;
    m3  = 54  + 126*Phi + 70*Phi^2;
    m4  = 13  + 31.5*Phi + 17.5*Phi^2;
    m5  = 4   + 7*Phi + 3.5*Phi^2;
    m6  = 3   + 7*Phi + 3.5*Phi^2;

    Mtrans = rho * A * L / (420 * (1 + Phi)^2) * [ ...
         m1,     L*m2,   m3,    -L*m4; ...
         L*m2,   L^2*m5, L*m4,  -L^2*m6; ...
         m3,     L*m4,   m1,    -L*m2; ...
        -L*m4,  -L^2*m6, -L*m2,  L^2*m5 ];

    % Rotary inertia contribution
    m7  = 36;
    m8  = 3 - 15*Phi;
    m9  = 4 + 5*Phi + 10*Phi^2;
    m10 = 1 + 5*Phi - 5*Phi^2;

    Mrot = rho * Iy / (30 * L * (1 + Phi)^2) * [ ...
         m7,    L*m8,   -m7,    L*m8; ...
         L*m8,  L^2*m9, -L*m8, -L^2*m10; ...
        -m7,   -L*m8,    m7,   -L*m8; ...
         L*m8, -L^2*m10, -L*m8, L^2*m9 ];

    Me = Mtrans + Mrot;

    % --- Stiffness matrix --------------------------------------------
    Ke = E * Iy / (L^3 * (1 + Phi)) * [ ...
         12,       6*L,           -12,       6*L; ...
         6*L,      (4+Phi)*L^2,   -6*L,      (2-Phi)*L^2; ...
        -12,      -6*L,            12,      -6*L; ...
         6*L,      (2-Phi)*L^2,   -6*L,      (4+Phi)*L^2 ];

end
