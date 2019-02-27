function [output] = if_else( cond, A, B )
n_cond = length( cond );
output = B;
for ind_cond = 1:n_cond
    if cond(ind_cond)
        output(ind_cond) = A(ind_cond);
    end
end

end
