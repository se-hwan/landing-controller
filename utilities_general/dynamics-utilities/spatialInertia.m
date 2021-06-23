function I = spatialInertia(a, b, c)
% 
% In the case of 3 arguments, a = mass, b = com, c = 3x3 inertia

if nargin == 1
    I = eye(6);
    I(1,1) = a(5);
    I(2,2) = a(6);
    I(3,3) = a(7);
    I(4:6,4:6) = a(1)*eye(3);
    cSkew = skew([a(2),a(3),a(4)]);
    
    % I(1:3,1:3) = [a(5) a(10) a(9);
    %     a(10) a(6) a(8);
    %     a(9) a(8) a(7)];
    I(1:3,4:6) = cSkew;
    I(4:6,1:3) = -cSkew';
    % I(4:6,4:6) = a(1)*eye(3);
    
else
    
    cSkew = skew(b);
    I = [c + a*(cSkew*cSkew'), a*cSkew;
        a*cSkew', a * eye(3)];
    
end