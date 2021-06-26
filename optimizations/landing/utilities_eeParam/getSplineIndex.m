% Description: returns spline index given time breakpoints and a time query point

function idx = getSplineIndex(time_breakPt, time_queryPt)
    
    % INPUT:    time_breakPt [1 x n] - vector of time breakpoints
    %           time_queryPt [1] - numeric value to find spline index of
    
    % OUPUT:    idx [1] - integer index of spline at query point

    f = generateIndexFunction(time_breakPt);
    
    idx = f(time_queryPt);
end

function idxFunction = generateIndexFunction(time_breakPt)
    % generates casadi function to find index of spline

    queryPt = casadi.MX.sym('queryPt', 1, 1);                          % query point for spline index
    L = low(time_breakPt, queryPt) + 1;                         % index for spline, +1 bc of MATLAB indexing

    idxFunction = casadi.Function('splineIndex', {queryPt}, {L});
end