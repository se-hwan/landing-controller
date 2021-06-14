function write_dynamics_to_file( params,H,C,p,pf,Jf,vf)

switch params.model
    case 'quad'
        q = sym('q',[7 1],'real');
        qd = sym('qd',[7 1],'real');
        Q = [q;qd];
        folder_out = 'Dynamics_mc2D/';
    case 'humanoid'
        q = sym('q',[8 1],'real');
        qd = sym('qd',[8 1],'real');
        Q = [q;qd];
        folder_out = 'Dynamics_h2D/';
    case 'quad3D'
        q = sym('q',[18 1],'real');
        qd = sym('qd',[18 1],'real');
        Q = [q;qd];
        folder_out = 'Dynamics_mc3D/';
end

parfor i = 1:6
    if i == 1
        disp('Starting H');
        matlabFunction(H,'File',[folder_out,'H_sym'],'Vars',{q});
        disp('Done with H');
    end
    
    if i == 2
        disp('Starting C');
        matlabFunction(C,'File',[folder_out,'C_sym'],'Vars',{Q});
        
        disp('Done with C');
    end
    
    if i == 3
        write_cell_to_file(p,'p_sym',{q},folder_out);
        disp('Done with p');
    end
    
    if i == 4
        write_cell_to_file(pf,'pf_sym',{q},folder_out);
        disp('Done with pf');
    end
    
    if i == 5
        write_cell_to_file(Jf,'Jf_sym',{q},folder_out);
        disp('Done with Jf');
    end
    
    if i == 6
        write_cell_to_file(vf,'vf_sym',{[q,qd]},folder_out);
        disp('Done with vf');
    end
    
end
end

