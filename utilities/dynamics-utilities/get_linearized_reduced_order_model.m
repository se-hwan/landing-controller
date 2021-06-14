function get_linearized_reduced_order_model(r)
% Creates symbolic functions for the linearized A and B matrices of the
% reduced order "potato" model for the desired robot

if strcmp(r,'mc')
    % Robot parameters
    m = sym('m','real');
    I = sym('I',[3 1]);
    Ib = diag(I);
    % Reference State
    p = sym('p',[3 1],'real'); % com position
    R = sym('R',[3 3],'real'); % orientation (via rotation matrix)
    v = sym('v',[3 1],'real'); % com velocity
    omega = sym('w',[3 1],'real'); % body angular velocity (in body frame)
    pf0 = sym('pf0',[3 1],'real'); % foot 0 position
    pf1 = sym('pf1',[3 1],'real'); % foot 1 position
    pf2 = sym('pf2',[3 1],'real'); % foot 2 position
    pf3 = sym('pf3',[3 1],'real'); % foot 3 position
    cs = sym('cs',[4 1],'real'); % contact schedule
    % Error State
    delta_p = sym('delta_p',[3 1],'real'); % com position
    delta_eta = sym('delta_eta',[3 1],'real'); % orientation (via rotation matrix)
    delta_v = sym('delta_v',[3 1],'real'); % com velocity
    delta_omega = sym('delta_w',[3 1],'real'); % body angular velocity (in body frame)
    delta_pf0 = sym('delta_pf0',[3 1],'real'); % foot 0 position
    delta_pf1 = sym('delta_pf1',[3 1],'real'); % foot 1 position
    delta_pf2 = sym('delta_pf2',[3 1],'real'); % foot 2 position
    delta_pf3 = sym('delta_pf3',[3 1],'real'); % foot 3 position
    % Control
    f0 = sym('f0',[3 1],'real'); % foot 0 position
    f1 = sym('f1',[3 1],'real'); % foot 1 position
    f2 = sym('f2',[3 1],'real'); % foot 2 position
    f3 = sym('f3',[3 1],'real'); % foot 3 position
    % Variation to Control
    delta_f0 = sym('delta_f0',[3 1],'real'); % foot 0 position
    delta_f1 = sym('delta_f1',[3 1],'real'); % foot 1 position
    delta_f2 = sym('delta_f2',[3 1],'real'); % foot 2 position
    delta_f3 = sym('delta_f3',[3 1],'real'); % foot 3 position
    
    % Variational Dynamics
    delta_u = [delta_f0;delta_f1;delta_f2;delta_f3];
    
    delta_xDot_footpos = sym(zeros(24,1));
    delta_x_footpos = [delta_p;delta_eta;delta_v;delta_omega;...
        delta_pf0;delta_pf1;delta_pf2;delta_pf3];
    
    delta_xDot = sym(zeros(12,1));
    delta_x = [delta_p;delta_eta;delta_v;delta_omega];
    
    delta_xDot(1:3,1) = delta_v;
    delta_xDot(4:6,1) = -skew(omega)*delta_eta+delta_omega;
    delta_xDot(7:9,1) = 1/m * (cs(1)*delta_f0 + ...
        cs(2)*delta_f1 + ...
        cs(3)*delta_f2 + ...
        cs(4)*delta_f3);
    t1 = cs(1)*R'*skew(pf0-p)*delta_f0 + ...
        cs(2)*R'*skew(pf1-p)*delta_f1 + ...
        cs(3)*R'*skew(pf2-p)*delta_f2 + ...
        cs(4)*R'*skew(pf3-p)*delta_f3;
    t2 = -cs(1)*R'*skew(f0)*delta_pf0 - ...
        cs(2)*R'*skew(f1)*delta_pf1 - ...
        cs(3)*R'*skew(f2)*delta_pf2 - ...
        cs(4)*R'*skew(f3)*delta_pf3;
    t3 = R'*skew(cs(1)*f0+cs(2)*f1+cs(3)*f2+cs(4)*f3)*delta_p; % done
    t4 = cs(1)*R'*skew(skew(pf0-p)*f0)*delta_eta + ...
        cs(2)*R'*skew(skew(pf1-p)*f1)*delta_eta + ...
        cs(3)*R'*skew(skew(pf2-p)*f2)*delta_eta + ...
        cs(4)*R'*skew(skew(pf3-p)*f3)*delta_eta; % done
    t5 = skew(Ib*omega)*delta_omega - skew(omega)*Ib*delta_omega;
    delta_xDot(10:12,1) = Ib\(t1 + t2 + t3 + t4 + t5);
    
    delta_xDot_footpos(1:12,1) = delta_xDot;
    delta_xDot_footpos(13:24,1) = zeros(12,1)-0.00001*delta_x_footpos(13:24,1); % small stabilizing term
    
    % Linear System Matrices
    A = jacobian(delta_xDot,delta_x);
    disp('Done with A matrix')
    Afp = jacobian(delta_xDot_footpos,delta_x_footpos);
    disp('Done with A matrix (foot pos included)')
    B = jacobian(delta_xDot,delta_u);
    disp('Done with B matrix')
    Bfp = jacobian(delta_xDot_footpos,delta_u);
    disp('Done with B matrix (foot pos included)')
    
    % Create Matlab functions
    matlabFunction(A,'File','Dynamics_VBL/A_sym','Vars',{m,I,p,R,v,omega,[pf0;pf1;pf2;pf3],[f0;f1;f2;f3],cs});
    disp('Done with A function')
    
    matlabFunction(Afp,'File','Dynamics_VBL/Afp_sym','Vars',{m,I,p,R,v,omega,[pf0;pf1;pf2;pf3],[f0;f1;f2;f3],cs});
    disp('Done with A function (foot pos included)')
    
    matlabFunction(B,'File','Dynamics_VBL/B_sym','Vars',{m,I,p,R,v,omega,[pf0;pf1;pf2;pf3],[f0;f1;f2;f3],cs});
    disp('Done with B function')
    
    matlabFunction(Bfp,'File','Dynamics_VBL/Bfp_sym','Vars',{m,I,p,R,v,omega,[pf0;pf1;pf2;pf3],[f0;f1;f2;f3],cs});
    disp('Done with B function (foot pos included)')
    
else
    error('Invalid robot. Use mc for mini cheetah and h for humanoid.... Humanoid not currently ready')
end