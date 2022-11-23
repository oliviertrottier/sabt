function [L_autocorr, Paths] = lengthcorr(Tree)
% Function to calculate the average branch length correlation across depth in a
% tree.
% 
% Input
% Tree = Tree structure
%
% Output
% L_autocorr = Length auto-correlation for a given depth difference
% Paths = structure containing all paths used to calculate the average auto-correlation.
%%
% Find the fieldname that gives the branch ID of the children branches.
if isfield(Tree, 'DaughtersID')
    ChildrenID_name = 'DaughtersID';
elseif isfield(Tree, 'ChildrenID')
    ChildrenID_name = 'ChildrenID';
else
    error('The children IDs fieldname cannot be determined');
end

% Find the branches ID that corresponds to tips.
Tips_ID = find(cellfun(@isempty, {Tree.(ChildrenID_name)}));

% Remove the tips that have a depth < 6.
Depth_min = 6;
Tips_ID([Tree(Tips_ID).Depth] < Depth_min) = [];
N_Tips = numel(Tips_ID);

% For all tips, record the lengths of the branches located on the path joining
% the said tip to the soma.
Paths = struct('Lautocorr', cell(N_Tips, 1), 'Lags', cell(N_Tips, 1));
Depth_max = max([Tree(Tips_ID).Depth]);
Branch_lengths = zeros(1, Depth_max);

for i = 1:N_Tips
    Tip_ID = Tips_ID(i);
    Tip_depth = Tree(Tip_ID).Depth;
    
    % Visit all branches on the path and record their length.
    Branch_ID = Tip_ID;
    Branch_depth = Tip_depth;
    while Branch_depth > 0
        Branch_lengths(Branch_depth) = Tree(Branch_ID).Length;
        
        Branch_ID = Tree(Branch_ID).ParentID;
        Branch_depth = Branch_depth-1;
    end
    
    % Calculate the auto correlation of the lengths with a maximum lag of
    % Tip_depth/2.
    [Paths(i).Lautocorr, Paths(i).Lags] = autocorr(Branch_lengths(1:Tip_depth), ceil(Tip_depth/2));
end

% Combine the auto-correlation function of all paths and calculate the
% average correlation at each lag.
L_autocorr_all = cell2mat({Paths.Lautocorr});
Lags_all = cell2mat({Paths.Lags});
L_autocorr = accumarray(Lags_all(:)+1, L_autocorr_all(:), [], @mean); % +1 because subs need to be > 0.
end