function box_lim = kin_box_limits(v, dir)
    
    v_max = 2; % velocity limit of box adjustment
    switch dir
        case 'x'
            box_max = 0.2;  % maximum box adjustment in direction
            if (abs(v) < v_max)
                box_lim = abs(v*(box_max/v_max));
            else
                box_lim = box_max;
            end
        case 'y'
            box_max = 0.15; % maximum box adjustment in direction
            if (abs(v) < v_max)
                box_lim = abs(v*(box_max/v_max));
            else
                box_lim = box_max;
            end
        otherwise
            error('Invalid option selected for kinematic box limits.')
    end
end