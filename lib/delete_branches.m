function Tree = delete_branches(Tree, Deleted_branches_ID, varargin)
% Function to delete branches from a tree structure.
%
% Input
% Tree = structure containing the branches
% Deleted_branches_ID = array containing IDs of the branches that are deleted.
%% Parse optional parameters
p = inputParser;
addParameter(p, 'Disconnect', true); % Disconnect branches before deletion. See below for details.
parse(p, varargin{:});
options = p.Results;
%%
% Function that deletes branches with ID=branchIDs and connects their
% respective sibling branch to their parent.
Tree_fieldnames = fieldnames(Tree);

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

% First, disconnect branches before deleting them. Disconnecting a branch
% removes its reference in the sibling and children IDs but does not remove
% its entry in the tree structure.
if options.Disconnect
    [Tree, Removed_branches_ID] = disconnect_branches(Tree, Deleted_branches_ID);
else
    Removed_branches_ID = Deleted_branches_ID;
end
%% Assign branch ID
% Find an ID name to use for storing the Branch ID.
% Record the rng state to restore it after finding a fieldname.
Rng = rng();
ID_fieldname = sprintf('ID%d',randi(1e6));
while ismember(ID_fieldname, Tree_fieldnames)
    ID_fieldname = sprintf('ID%d',randi(1000));
end
rng(Rng);

% Assign ID to be used for further reference.
IDs = num2cell(1:numel(Tree));
[Tree.(ID_fieldname)] = deal(IDs{:});
IDs = cell2mat(IDs);
Retained_branches_ID = setdiff(IDs,Removed_branches_ID);
%% Perform sanity checks on deleted and non-deleted branches.
% Check if deleted branches are referenced by non-deleted branches.
Tree_new = Tree(Retained_branches_ID);
Retained_parents_ID = [Tree_new.ParentID];
if any(ismember(Retained_parents_ID, Removed_branches_ID))
    error('Some removed branches are referenced as Parent IDs.')
end
Retained_siblings_ID = [Tree_new.(SiblingID_name)];
if any(ismember(Retained_siblings_ID, Removed_branches_ID))
    error('Some removed branches are referenced as Sibling IDs.')
end
Retained_children_ID = [Tree_new.(ChildrenID_name)];
if any(ismember(Retained_children_ID, Removed_branches_ID))
    error('Some removed branches are referenced as Children IDs.')
end

% Check that the sibling ID of deleted branches is non-empty.
Siblings_ID = {Tree_new.(SiblingID_name)};
has_empty_siblings = cellfun(@isempty, Siblings_ID);

% Remove soma branches from the check.
is_soma_branch = [Tree_new(has_empty_siblings).ParentID] == 0;
has_empty_siblings(is_soma_branch) = false;
if any(has_empty_siblings)
    has_empty_siblings_ID = Retained_branches_ID(has_empty_siblings);
    
%     % Prune soma branches.
%     % After deleting branches, there may be only 1 branch that connects to
%     % the soma. When this is the case, delete the first branch and connect its
%     % children to the soma. The children will form the new soma branches.
%     Soma_branches_ID = has_empty_siblings_ID([Tree_new(has_empty_siblings_ID).ParentID] == 0);
%     if ~isempty(Soma_branches_ID)
%         Tree_new = delete_branches(Tree_new, Soma_branches_ID);
%         
%         % Recheck for empty siblings.
%         SiblingIDs = {Tree_new.(SiblingID_name)};
%         has_empty_siblings = cellfun(@isempty, SiblingIDs);
%         if any(has_empty_siblings)
%             error('Some branches have no Sibling ID.');
%         end
%     else
%         
%     end
    error('Some branches have no Sibling ID.');
end

% Check that the Parent ID of deleted branchs is non-empty.
ParentIDs = {Tree_new.ParentID};
has_empty_parent = cellfun(@isempty, ParentIDs);
if any(has_empty_parent)
    disp(Retained_branches_ID(has_empty_parent));
    error('Some branches have no Parent ID.')
end
%% Overwrite the old tree structure.
Tree = Tree_new;
%% Fix ID references before reordering the branches.
Old_IDs = [Tree.(ID_fieldname)];
N_branches = numel(Tree);
New_IDs = 1:N_branches;

ParentIDs = rep([Tree.ParentID], Old_IDs, New_IDs);
SiblingID = rep({Tree.(SiblingID_name)}, Old_IDs, New_IDs);
ChildrenIDs = rep({Tree.(ChildrenID_name)}, Old_IDs, New_IDs);
for i = 1:N_branches
    Tree(i).(ID_fieldname) = New_IDs(i);
    Tree(i).ParentID = ParentIDs(i);
    Tree(i).(SiblingID_name) = SiblingID{i};
    Tree(i).(ChildrenID_name) = ChildrenIDs{i};
end
%% Reorder the tree structure.
Tree = order_tree(Tree);

% Check if ID references go out of bounds.
All_IDs = cell2mat(cellfun(@double,{[Tree.ParentID], [Tree.(SiblingID_name)], [Tree.(ChildrenID_name)]},'Uni',0));
if any(All_IDs > N_branches)
    error('Some branch ID references are out-of-bound.')
end

% Remove ID fieldname.
Tree = rmfield(Tree, ID_fieldname);
end