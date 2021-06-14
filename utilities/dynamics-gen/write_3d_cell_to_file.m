function write_3d_cell_to_file( C, name, vars )

parfor i = 1:length(C)
file_name = ['3d_dynamics_out/' name num2str(i)];
matlabFunction(C{i},'File',file_name,'Vars',vars);
end


end

