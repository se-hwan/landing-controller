% Description: generates casadi function that can finds the lower index for a piecewise function

function [evalFunction, idxFunction, L] = getSplineFunctions(phaseDurations, coeff, splineDurations, polyType)

    % INPUT:    phaseDurations [1 x n] - vector of phase durations of contact, can be symbolic
    %           coeff [polyOrder x nSplines] - matrix of coefficients that parametrize polynomial
    %           polyType - 'Hermite' or 'Power', polynomial parametrization type
    %           splineDurations [1 x nSplines] - vector of durations for each spline
    
    % OUTPUT:   evalFunction - returns casadi function that evaluates spline at query point
    %           idxFunction - returns casadi function that finds lower index of spline phase given some query point
    %           L [1] - index of spline at query point
    
    addpath(genpath('../../utilities')); 
    import casadi.*
    
    queryPt = MX.sym('queryPt', 1, 1);          % query point for spline index
    polyOrder = length(coeff(:,1)) - 1;         % polynomial order

    if (isempty(splineDurations))
        N = length(phaseDurations) + 1;         % length of time breakpoint vector
        breakpoints = casadi.MX.zeros(1, N);    % breakpoints in time for polynomial switching
        for i = 1:N-1
            breakpoints(i+1) = sum(phaseDurations(1:i));
        end
    else
        N = length(splineDurations) + 1;
        breakpoints = casadi.MX.zeros(1, N);    % breakpoints in time for polynomial switching
        for i = 1:N-1
            breakpoints(i+1) = sum(splineDurations(1:i));
        end
    end
    
    L = low(breakpoints, queryPt) + 1;      % index for spline, +1 bc of MATLAB indexing
    currentCoeff = coeff(:, L);             % coefficients at index point
    
    switch polyType
        case 'Hermite'     % Hermite parameterization of polynomial
            x0 = currentCoeff(1); x0_dot = currentCoeff(2);
            x1 = currentCoeff(3); x1_dot = currentCoeff(4);
            
            DeltaT = splineDurations(L);
            
            a0 = x0; a1 = x0_dot; a2 = (-DeltaT^(-2))*(3*(x0 - x1) + DeltaT*(2*x0_dot + x1_dot));
            a3 = (DeltaT^(-3))*(2*(x0 - x1) + DeltaT*(x0_dot + x1_dot));
            
            hermitianCoeff = [a3 a2 a1 a0]';
            hermitianCoeff_dot = [3*a3 2*a2 a1]';
            
            splineEval = dot(hermitianCoeff, [queryPt.^[polyOrder:-1:0]]');
            splineEval_dot = dot(hermitianCoeff_dot, [queryPt.^[polyOrder-1:-1:0]]');
            
        case 'Power'        % "Standard" power definition of polynomial
            currentCoeff_dot = [polyOrder:-1:1].*currentCoeff(1:end-1)';
            currentCoeff_ddot = [(length(currentCoeff_dot)-1):-1:1].*currentCoeff_dot(1:end-1);
            
            splineEval = dot(currentCoeff, [queryPt.^[polyOrder:-1:0]]');
            splineEval_dot = dot(currentCoeff_dot', [queryPt.^[polyOrder-1:-1:0]]');
            splineEval_ddot = dot(currentCoeff_ddot', [queryPt.^[polyOrder-2:-1:0]]');
            
            % not needed for Hermite splines from foot forces, foot positions
            evalFunction.x_ddot = casadi.Function('splineEval_ddot', {queryPt}, {splineEval_ddot});
            
        otherwise
            error('Polynomial parametrization type invalid, choose ''Power'' or ''Hermite.''');
    end

    idxFunction = casadi.Function('splineIndex', {queryPt}, {L});
    evalFunction.x = casadi.Function('splineEval', {queryPt}, {splineEval});
    evalFunction.x_dot = casadi.Function('splineEval_dot', {queryPt}, {splineEval_dot});
    
end