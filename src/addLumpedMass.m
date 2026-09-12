function M = addLumpedMass(M, dofMap, nodeId, mass, inertia)
% ADDLUMPEDMASS  Add a lumped point mass and/or rotary inertia to a node.
%
%   M = addLumpedMass(M, dofMap, nodeId, mass, inertia)
%
%   Inputs
%     M         Global mass matrix to modify (returned, modified copy)
%     dofMap    nNodes x 2 DOF map, as returned by assembleBeamFE
%     nodeId    Node index (1-based) to which the lumped mass is added
%     mass      Lumped mass [kg] (added to the translational DOF)
%     inertia   Lumped rotary inertia [kg*m^2] (added to the rotational
%               DOF). Pass 0 if not applicable.
%
%   Typical uses: an accelerometer mass at a measurement point, a
%   sensor, or (with two calls, one per relevant node) the suspended /
%   non-suspended mass contributions of a dynamic vibration absorber.

    uDof     = dofMap(nodeId, 1);
    thetaDof = dofMap(nodeId, 2);

    M(uDof, uDof)         = M(uDof, uDof) + mass;
    M(thetaDof, thetaDof) = M(thetaDof, thetaDof) + inertia;
end
