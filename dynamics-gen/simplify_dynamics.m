function [ D_out ] = simplify_dynamics( D, in_parallel )
D_out = cell(1,length(D));
disp(['Simplify dynamics called with ' num2str(length(D)) ' cells']);
if(in_parallel)
    parfor i = 1:length(D)
        D_out{i} = simplify_cell_array(D{i},1);
        disp(['Done with cell ' num2str(i)]);
    end
else
    parfor i = 1:length(D)
        D_out{i} = simplify_cell_array(D{i},0);
    end
end

end

