function [x, fval, exitflag] = inverse_kinematics(q_foot, model, q_body, q_foot_guess, cs)
    [x, fval, exitflag] = fsolve(@(q) foot_residual(q_foot, model, q_body, q, cs), q_foot_guess);
end

function residual = foot_residual(q_foot, model, q_body, q_foot_guess, cs)

% inverse kinematics solver for joint angles given foot positions in absolute global frame

% IN: q_foot        - global foot locations (12 x 1)
%     model         - robot model
%     q_body        - floating base state (6 x 1)
%     q_foot_guess  - symbolic variable solution, with initial guess (12 x 1)
%     cs            - contact schedule (4 x 1)

    pf0 = get_forward_kin_foot(model,[q_body;q_foot_guess]);
    footstep_guess = reshape(cell2mat(pf0),[12,1]);
    residuals = q_foot-footstep_guess;
    cs_array = [ones(3,1)*cs(1);ones(3,1)*cs(2);ones(3,1)*cs(3);ones(3,1)*cs(4)];
    cs_array = ones(12,1);
    residuals2 = residuals.*cs_array;
    residual = residuals2.^2;
end
