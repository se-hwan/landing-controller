% checks for ground penetration from elbows and for resteps
function isValidSoln = check_solution(model, q_star, f_star)
%%
    N = length(q_star(1,:));
    
    td_idx = zeros(4,1);
    for leg = 1:4
        xyz_idx = 3*leg-2 : 3*leg;
        f_leg = f_star(xyz_idx, :);
        td_idx(leg) = find(f_leg(3,:) > 1, 1);
    end
    
    for i = 1:N-1
       pos_knee = get_forward_kin_knee(model, q_star(:, i));
       pos_foot = get_forward_kin_foot(model, q_star(:, i));
       for leg = 1:4
           xyz_idx = 3*leg-2 : 3*leg;
           if ((pos_knee{leg}(3) < -0.025) || (pos_foot{leg}(3) < -0.025))
               isValidSoln = false;
               return;
           elseif ((i > td_idx(leg)) && (f_star(xyz_idx(3), i) < 5))
               isValidSoln = false;
               return;
           else
               isValidSoln = true;
           end
       end
    end

end