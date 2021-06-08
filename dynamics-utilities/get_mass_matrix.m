function  [H,Ic] = get_mass_matrix( model, q, is_opt)

% Computes mass matrix, H, for a floating base system

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
for i = 1:5
    Xup{i} = zeros(6,6);
end

Xup{6} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew(q(1:3)) R_world_to_body];

for i = 7:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end

%% Composite Inertia
IC = model.I;
for i = model.NB:-1:7
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
end

%% Mass Matrix
if (is_opt==1)
    H = casadi.MX(zeros(model.NB));
elseif (is_opt==0)
    H = zeros(model.NB);
else %is_opt == 2
    H = casadi.SX.sym('H',model.NB,model.NB);
end
H(1:6,1:6) = IC{6};

for i = 7:model.NB
    fh = IC{i} * S{i};
    H(i,i) = S{i}' * fh;
    
    fh = Xup{i}' * fh;
    j = model.parent(i);
    while j > 6
        H(i,j) = S{j}' * fh;
        H(j,i) = H(i,j);
        fh = Xup{j}' * fh;
        j = model.parent(j);
    end
    
    H(1:6,i) = fh;
    H(i,1:6) = fh';
end

% Composite Rigid Body Inertia - Body frame
U = [eye(6) zeros(6,model.NB-6)];
Ic = U*H*U';
