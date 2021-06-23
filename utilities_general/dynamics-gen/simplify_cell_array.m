function [ C_out ] = simplify_cell_array( C, in_parallel )
C_out = cell(1,length(C));
if(in_parallel)
    parfor i = 1:length(C)
        C_out{i} = simplify(C{i},'Seconds',60);
    end
else
    for i = 1:length(C)
        C_out{i} = simplify(C{i},'Seconds',60);
    end
end

end

