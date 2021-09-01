function showmotion_floatingBase( model, t_data, q_data )
% Use this when rendering any floating base robot with spatial v2. It
% corrects for the mismatch between RPY euler angles and the index ordering
% of ZYX needed for proper rendering

roll_temp_var = q_data(4,:);
q_data(4,:) = q_data(6,:);
q_data(6,:) = roll_temp_var;
showmotion( model, t_data, q_data ); 