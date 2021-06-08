function [ H ] = add_offboard_rotors( H_in, I_rot, gr, has_rot )

H = H_in + diag(has_rot) * gr * I_rot * gr;


end

