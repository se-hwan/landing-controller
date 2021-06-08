function [ params, model ] = get_3d_robot_model( )
% currently only does mini cheetah

params = get_robot_params(1);
[ pmw_model, jtype, Xtree, I, Xrot, Irot ] = mc_sim_model(  );

model = sim_model_to_opt_model(pmw_model,jtype,Xtree,I,Xrot,Irot,params);

end

