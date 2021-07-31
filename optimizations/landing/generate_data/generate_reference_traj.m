function [X_ref, U_ref] = generate_reference_traj(q_init, qd_init, q_term, qd_term, N)

    c_ref = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.2 0.1 -0.35],1,4)';
    f_ref = zeros(12,1);

    for i = 1:6
        Xref_val(i,:)   = linspace(q_init(i),q_term(i),N);
        Xref_val(6+i,:) = linspace(qd_init(i),qd_term(i),N);
    end
    for leg = 1:4
        for xyz = 1:3
            Uref_val(3*(leg-1)+xyz,:)    = Xref_val(xyz,1:end-1) + c_ref(3*(leg-1)+xyz);
            Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N-1);
        end
    end
    
    X_ref = Xref_val; U_ref = Uref_val;
end