function  write_3d_dynamics_to_file( H,C,p,pf,Jf,Jdqdf,vf )
q = sym('q',[11 1],'real');
qd = sym('qd',[11 1],'real');

Q = [q;qd];

    matlabFunction(H,'File','3d_dynamics_out/H_sym','Vars',{q});
    disp('Done with H');

    matlabFunction(C,'File','3d_dynamics_out/C_sym','Vars',{[Q]});
    disp('Done with C');

    write_3d_cell_to_file(p,'p_sym',{q}); 
    disp('Done with p');

    write_3d_cell_to_file(pf,'pf_sym',{q});
    disp('Done with pf');

    write_3d_cell_to_file(Jf,'Jf_sym',{q});
    disp('Done with Jf');
    
%     write_3d_cell_to_file(Jdqdf,'Jdqdf_sym',{[q, qd]});
%     disp('Done with Jdqdf');
% 
%     write_3d_cell_to_file(vf,'vf_sym',{[q,qd]});
%     disp('Done with vf');

end

