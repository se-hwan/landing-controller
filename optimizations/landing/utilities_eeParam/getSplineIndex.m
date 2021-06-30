% Description: returns spline index given time breakpoints and a time query point

function [idx, f] = getSplineIndex(time_breakPt, time_queryPt)
    
    % INPUT:    time_breakPt [1 x n] - vector of time breakpoints
    %           time_queryPt [1] - numeric value to find spline index of
    
    % OUPUT:    idx [1] - integer index of spline at query point

    f = generateIndexFunction(time_breakPt);
    
    idx = f(time_queryPt);
end

