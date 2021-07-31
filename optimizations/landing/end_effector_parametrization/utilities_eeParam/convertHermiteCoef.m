% Description: returns power polynomial coefficients of Hermite parameterized polynomial

function p_coef = convertHermiteCoef(h_coef, duration)
    
    % INPUT:    h_coef [n x m] - n rows of order m Hermite polynomials
    %           duration [n x 1] - n rows of durations for the coefficients
    
    % OUPUT:    p_coef [n x m] - converted power polynomial coefficients

    n = length(h_coef(:, 1));
    m = length(h_coef(1, :));
    switch class(h_coef)
        case 'casadi.MX'
            coef = casadi.MX.zeros(n, m);
            for i = 1:n
                x0 = h_coef(i, 1); x0_dot = h_coef(i, 2);
                x1 = h_coef(i, 3); x1_dot = h_coef(i, 4);
                a0 = x0; a1 = x0_dot;
                a2 = (-duration(i)^(-2))*(3*(x0 - x1) + duration(i)*(2*x0_dot + x1_dot));
                a3 = (duration(i)^(-3))*(2*(x0 - x1) + duration(i)*(x0_dot + x1_dot));
                coef(i, :) = [a3 a2 a1 a0];
            end
        otherwise
            coef = zeros(n, m);
            for i = 1:n
                x0 = h_coef(i, 1); x0_dot = h_coef(i, 2);
                x1 = h_coef(i, 3); x1_dot = h_coef(i, 4);
                a0 = x0; a1 = x0_dot;
                a2 = (-duration(i)^(-2))*(3*(x0 - x1) + duration(i)*(2*x0_dot + x1_dot));
                a3 = (duration(i)^(-3))*(2*(x0 - x1) + duration(i)*(x0_dot + x1_dot));
                coef(i, :) = [a3 a2 a1 a0];
            end
    end
    
    p_coef = coef;
end
