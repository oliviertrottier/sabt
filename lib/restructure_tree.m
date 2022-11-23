function Tree = restructure_tree(Tree,varargin)
% Function to restructure a tree structure. 
%
% Zero-length branches are removed and branches are reordered according to
% their 1) Tree depth, 2) Parent ID and 3) Branching angle.

% In addition, unnecessary fields are removed leaving only the following fields:
% 'PointsPos', 'PointsBirthTime', 'ParentID', SiblingID_name, ChildrenID_name,
% 'Length', 'Depth', 'BranchingAngle'.
%% Parse optional parameters
p = inputParser;
addParameter(p,'RemoveZeroLength',true); % Remove zero-length branches.
addParameter(p,'RemoveInactive',false); % Remove inactive branches using the "isActive" field in Tree.
parse(p,varargin{:});
options = p.Results;
%% Determine the field names for the children and sibling ID.
% Determine the field name of the children IDs.
if isfield(Tree, 'DaughtersID')
    ChildrenID_name = 'DaughtersID';
elseif isfield(Tree, 'ChildrenID')
    ChildrenID_name = 'ChildrenID';
else
    error('The children IDs fieldname cannot be determined');
end

% Determine the field name of the sibling IDs.
if isfield(Tree, 'SisterID')
    SiblingID_name = 'SisterID';
elseif isfield(Tree, 'SiblingID')
    SiblingID_name = 'SiblingID';
else
    error('The sibling ID fieldname cannot be determined');
end
%% Remove zero-length branches
% Determine which branches are removed.
if options.RemoveInactive
    if isfield(Tree,'isActive')
        is_branch_removed = ~[Tree.isActive];
    else
        error('The field isActive is absent from the tree structure.');
    end
else
    is_branch_removed = false(size(Tree));
end

if options.RemoveZeroLength
    is_branch_removed = is_branch_removed | [Tree.Length] == 0;
end

if nnz(is_branch_removed) > 0
    % Delete the removed branches.
    Deleted_branches_ID = find(is_branch_removed);
    Tree = delete_branches(Tree,Deleted_branches_ID,'Disconnect',true);
else
    % If no branches are removed, simply reorder the tree.
    Tree = order_tree(Tree);
end
%% Change the data type to reduce memory.
Tree = reduce_memory(Tree);
%% Remove unnecessary fields and order the remaining fieldnames.
Tree_fieldnames = fieldnames(Tree);
Fieldnames_final = {'PointsPos', 'PointsBirthTime', 'ParentID', SiblingID_name, ChildrenID_name, 'Length', 'Depth', 'BranchingAngle'};
is_field_kept = ismember(Tree_fieldnames,Fieldnames_final);
Fields_kept = Tree_fieldnames(is_field_kept);
Fields_removed = Tree_fieldnames(~is_field_kept);

Tree = rmfield(Tree, Fields_removed);
Tree = orderfields(Tree, Fields_kept);
end