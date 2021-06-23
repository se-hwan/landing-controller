function write_cell_to_file( C, name, vars , out )

parfor i = 1:length(C)
file_name = [out name num2str(i)];
matlabFunction(C{i},'File',file_name,'Vars',vars);
end


end

