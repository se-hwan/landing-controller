function [x, fval, exitflag] = inverse_kinematics(q_foot, model, q_body, q_foot_guess)
    [x, fval, exitflag] = fsolve(@(q) foot_residual(q_foot, model, q_body, q), q_foot_guess, optimset('Display','off'));
end

function residual = foot_residual(q_foot, model, q_body, q_foot_guess)

% inverse kinematics solver for joint angles given foot positions in absolute global frame

% IN: q_foot        - global foot locations (12 x 1)
%     model         - robot model
%     q_body        - floating base state (6 x 1)
%     q_foot_guess  - symbolic variable solution, with initial guess (12 x 1)
%     cs            - contact schedule (4 x 1)

    pf0 = get_forward_kin_foot(model,[q_body;q_foot_guess]);
    footstep_guess = reshape(cell2mat(pf0),[12,1]);
    residuals2 = q_foot-footstep_guess;
    residual = residuals2.^2;
end
