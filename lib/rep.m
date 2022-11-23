function x = rep(x, oldvalues, newvalues)
% Function to replace numbers in x by new numbers.
%
% Input
% x = cell or double array for whose numbers will be replaced.
%     If x is a cell, replace is applied to each element of the cell array.
% oldvalues = array representing numbers that will be replaced.
% newvalues = array representing the replacement of each number in
% oldvalues.
if iscell(x)
    for n = 1:numel(x)
        x{n} = rep(x{n}, oldvalues, newvalues);
    end
elseif ~isempty(x)
    [is_x_old, x_old_locs] = ismember(x, oldvalues);
    if any(is_x_old)
        x(is_x_old) = newvalues(x_old_locs(is_x_old));
    end
end
end