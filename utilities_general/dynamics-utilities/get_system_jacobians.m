function J = get_system_jacobians( model, q, is_symbolic)
% computes system Jacobians via steps outlined in Orin "Centroidal Momentum 
% Matrix of a Humanoid Robot: Structure and Properties"
%
% Xup (transforms parent frame to ith frame)
n = model.NB; % number of degrees of freedom
N = model.NB-5; % number of links in tree (we count floating base as single link)

if is_symbolic
    P = sym(zeros(6*N,6*N));
    Phi = sym(zeros(6*N,n));
else
    P = zeros(6*N,6*N);
    Phi = zeros(6*N,n);
end

% Floating base
P(1:6,1:6) = eye(6);
Phi(1:6,1:6) = eye(6);

for i = 7:model.NB
    P(6*(i-6)+1:6*(i-6)+6,6*(i-6)+1:6*(i-6)+6) = eye(6);
    
    [ XJ, S ] = jcalc( model.jtype{i}, q(i));        % joint xform, joint motion subspace
    Xup = XJ * model.Xtree{i};                       % xform from parent
    
    if model.parent(i) ~= 0
        j = model.parent(i);
        P(6*(i-6)+1:6*(i-6)+6,6*(j-6)+1:6*(j-6)+6) = -Xup;
    end

    Phi(6*(i-6)+1:6*(i-6)+6,i) = S;
     
end

J_ = P\Phi;

% Parse into individual Jacobians?
for i = 1:N
   J{i} = J_(6*(i-1)+1:6*(i-1)+6,:);
end