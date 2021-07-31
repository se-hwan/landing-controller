function idxFunction = generateIndexFunction(time_breakPt)
    % generates casadi function to find index of spline

    queryPt = casadi.MX.sym('queryPt', 1, 1);                          % query point for spline index
    L = low(time_breakPt, queryPt) + 1;                         % index for spline, +1 bc of MATLAB indexing

    idxFunction = casadi.Function('splineIndex', {queryPt}, {L});
end