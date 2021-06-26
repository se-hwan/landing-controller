% Description: given set of coefficients, returns coefficients of their derivatives

function deriv_coef = getDerivCoef(coef)

    % INPUT:    coef [n x order + 1] - matrix of polynomial coefficients of order 'order'
    %                                - ordered from highest polynomial to lowest
    
    % OUTPUT:   deriv_coef [n x order] - matrix of polynomial derivative coefficients
    
    n = length(coef(:, 1));
    m = length(coef(1, :));
    order = m - 1;
    switch class(coef)
        case 'casadi.MX'
            dCoef = casadi.MX.zeros(n, m - 1);
            for i = 1:n
                dCoef(i, :) = [order:-1:1].*coef(i, 1:end - 1);
            end
        otherwise
            dCoef = zeros(n, m - 1);
            for i = 1:n
                dCoef(i, :) = [order:-1:1].*coef(i, 1:end - 1);
            end
    end
    
    deriv_coef = dCoef;
end