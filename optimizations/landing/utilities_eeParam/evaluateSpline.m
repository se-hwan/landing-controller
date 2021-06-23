% Description: evaluates spline with NUMERICAL coefficients

function [eval] = evaluateSpline(coeff, durations, time, polyType)

    % INPUT:    coeff [polyOrder x nSplines] - matrix of coefficients that parametrize polynomial
    %           durations [1 x nSplines] - durations of each piecewise spline
    %           time [1] - query time to evaluate spline
    %           polyType - 'Hermite' or 'Power', polynomial parametrization type
    
    % OUTPUT:   eval - evaluation of spline at queried time
    
    polyOrder = length(coeff(:,1)) - 1;         % polynomial order
    N = length(durations) + 1;
    breakpoints = [0];
    
    for i = 1:N-1
        breakpoints(i+1) = sum(durations(1:i));
    end
    
    idx = find(time >= breakpoints);
    idx = idx(end);
    if idx > length(durations)
        idx = length(durations);
    end
    currentCoeff = coeff(:, idx);
    
    switch polyType
        case 'Hermite'
            x0 = currentCoeff(1); x0_dot = currentCoeff(2);
            x1 = currentCoeff(3); x1_dot = currentCoeff(4);
            
            DeltaT = durations(idx);
            
            a0 = x0; a1 = x0_dot; a2 = (-DeltaT^(-2))*(3*(x0 - x1) + DeltaT*(2*x0_dot + x1_dot));
            a3 = (DeltaT^(-3))*(2*(x0 - x1) + DeltaT*(x0_dot + x1_dot));
            
            hermitianCoeff = [a3 a2 a1 a0]';
            splineEval = dot(hermitianCoeff, [time.^[polyOrder:-1:0]]');
            
        case 'Power'
            splineEval = dot(currentCoeff, [time.^[polyOrder:-1:0]]');
            
        otherwise
            error('Invalid selection. Choose ''Hermite'' or ''Power'' for evaluation');
    end
    
    eval = splineEval;
end