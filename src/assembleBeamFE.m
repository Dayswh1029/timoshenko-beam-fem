function model = assembleBeamFE(segments, nElemPerSegment)
% ASSEMBLEBEAMFE  Build a global mass/stiffness model for a beam made of
%   one or more prismatic segments, discretized with Timoshenko beam
%   elements. Fully data-driven: change the segment table and/or the
%   element count and the mesh, matrices and node coordinates are
%   rebuilt automatically -- no hand-editing of matrix indices required.
%
%   model = assembleBeamFE(segments, nElemPerSegment)
%
%   Inputs
%     segments          Struct array, one entry per prismatic segment,
%                        in order from the fixed/left end to the free/
%                        right end. Each entry must have fields:
%                          .L      segment length            [m]
%                          .E      Young's modulus            [Pa]
%                          .G      shear modulus              [Pa]
%                          .rho    density                    [kg/m^3]
%                          .A      cross-section area         [m^2]
%                          .Iy     second moment of area      [m^4]
%                          .chi    Timoshenko shear correction factor
%
%     nElemPerSegment   Scalar, or vector the same length as segments,
%                        giving how many elements to use to discretize
%                        each segment (e.g. 28, or [10 5 13]).
%
%   Output (struct)
%     model.nNodes      total number of nodes
%     model.nDofFull    2 * nNodes (before constraints)
%     model.x           nodal coordinates along the beam axis [m]
%     model.M           global (unconstrained) mass matrix
%     model.K           global (unconstrained) stiffness matrix
%     model.dofMap      nNodes x 2 matrix, dofMap(i,1:2) = global DOF
%                        indices [u_i, theta_i] for node i
%
%   Example
%     seg.L = 0.28; seg.E = 70e9; seg.G = seg.E/(2*(1+0.33));
%     seg.rho = 2700; seg.A = 30e-3*3e-3; seg.Iy = 30e-3*(3e-3)^3/12;
%     seg.chi = (12+11*0.33)/(10*(1+0.33));
%     model = assembleBeamFE(seg, 28);   % same as a 28-element cantilever
%
%   This is the generalization of the fixed 28-element / 4-DOF-window
%   assembly loop used in the original coursework script: any beam
%   length, any number of segments (e.g. a beam + a bonded damper
%   spring, as in cantilever + dynamic absorber problems), and any mesh
%   density can be built without touching the assembly logic itself.

    nSeg = numel(segments);
    if isscalar(nElemPerSegment)
        nElemPerSegment = repmat(nElemPerSegment, 1, nSeg);
    end
    assert(numel(nElemPerSegment) == nSeg, ...
        'nElemPerSegment must be a scalar or have one entry per segment.');

    % --- Build element list: element i has length Le(i) and belongs to
    %     segment segId(i), so we know which material/section to use.
    elemLength = [];
    segId      = [];
    for s = 1:nSeg
        nE = nElemPerSegment(s);
        elemLength = [elemLength, repmat(segments(s).L / nE, 1, nE)]; %#ok<AGROW>
        segId      = [segId, repmat(s, 1, nE)];                       %#ok<AGROW>
    end

    nElem  = numel(elemLength);
    nNodes = nElem + 1;
    nDof   = 2 * nNodes;

    % --- Nodal coordinates (cumulative element length)
    x = zeros(nNodes, 1);
    for e = 1:nElem
        x(e+1) = x(e) + elemLength(e);
    end

    % --- DOF map: node i -> global dof [2i-1 (u), 2i (theta)]
    dofMap = [(1:nNodes)'*2 - 1, (1:nNodes)'*2];

    % --- Assemble
    M = zeros(nDof);
    K = zeros(nDof);
    for e = 1:nElem
        s = segments(segId(e));
        [Me, Ke] = timoshenkoElement(s.E, s.G, s.rho, s.A, s.Iy, ...
                                      elemLength(e), s.chi);
        gdof = [dofMap(e,1), dofMap(e,2), dofMap(e+1,1), dofMap(e+1,2)];
        M(gdof, gdof) = M(gdof, gdof) + Me;
        K(gdof, gdof) = K(gdof, gdof) + Ke;
    end

    model.nNodes   = nNodes;
    model.nDofFull = nDof;
    model.x        = x;
    model.M        = M;
    model.K        = K;
    model.dofMap   = dofMap;
end
