function residual = fwd_kin_foot_resid(q_foot,model,q_body,q_foot_guess,cs)
pf0=get_forward_kin_foot(model,[q_body;q_foot_guess]); %this way q_body isnt varied to solve
footstep_guess=reshape(cell2mat(pf0),[12,1]);
residuals=q_foot-footstep_guess;
cs_array=[ones(3,1)*cs(1);ones(3,1)*cs(2);ones(3,1)*cs(3);ones(3,1)*cs(4)];%only care about the foot touching the ground
cs_array=ones(12,1);
residuals2=residuals.*cs_array;
%residual=sum(residuals2.^2);
residual=residuals2.^2;
end
