function Tree_ordered = order_tree(Tree)
% Function to order a Tree branches in ascending depth, orientation and branching angle.
% Root branches are defined as branches whose Parent ID=0 and initialized
% at depth = 1.

% Input
% Tree = input tree structure.

% Output
% Tree_ordered = ordered tree structure.
%% Determine field names of sibling and children,
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
%% Define the depth of each branch.
% Find the number of branches at the soma and initialize their depth to 1.
Parent_IDs = reshape([Tree.ParentID],[],1);
Root_branches_ID = find(Parent_IDs==0);

Branches_ID = Root_branches_ID;
depth = 1;
while ~isempty(Branches_ID)
    % Assign depth.
    [Tree(Branches_ID).Depth] = deal(depth);
    
    % Move to the next depth level.
    Branches_ID = [Tree(Branches_ID).(ChildrenID_name)];
    depth = depth + 1;
end
Branches_depth = reshape([Tree(:).Depth],[],1);
%% Calculate the branching angle of each branch. 
% This will be used to order the branches.
N_branches = numel(Tree);
Branches_branching_angle = nan(N_branches, 1);
if isfield(Tree,'PointsPos')
    for i = 1:N_branches
        if size(Tree(i).PointsPos,1) > 1
            Branching_vec = diff(Tree(i).PointsPos(1:2, :));
            Branches_branching_angle(i) = mod(atan2(Branching_vec(2), Branching_vec(1)),2*pi);
        end
    end
end
%% Calculate the orientation of each branch's branchpoint with respect to the soma.
Branchpoints_theta = nan(N_branches, 1);
if isfield(Tree,'PointsPos')
    Soma_position = Tree(Root_branches_ID(1)).PointsPos(1, :);
    for i=1:N_branches
        Branchpoints_vec = Tree(i).PointsPos(1, :) - Soma_position;
        Branchpoints_theta(i) = atan2(Branchpoints_vec(2), Branchpoints_vec(1));
    end
end
%% Reorder the Tree in ascending depth, branchpoint orientation and branching angle.
% First, find the map that maps the old IDs to the new IDs such that
% i = ID_old2new(ID_new2old(i)).

% ID_new2old(i) gives the old branch ID of the branch whose new ID is i.
[~,ID_new2old] = sortrows([Branches_depth, Branchpoints_theta, Branches_branching_angle]);

% ID_old2new(i) gives the new branch ID of the branch whose old ID is i.
[~,ID_old2new] = sort(ID_new2old);

Tree_ordered = reshape(Tree(ID_new2old),size(Tree));

% Fix the parent, sibling and children ID references.
for i=1:N_branches
    ParID = Tree_ordered(i).ParentID;
    if ParID > 0
        Tree_ordered(i).ParentID = ID_old2new(ParID);
    end
    
    SiblingIDs = Tree_ordered(i).(SiblingID_name);
    if SiblingIDs > 0
        Tree_ordered(i).(SiblingID_name) = ID_old2new(SiblingIDs)';
    end
    
    ChildrenIDs = Tree_ordered(i).(ChildrenID_name);
    Tree_ordered(i).(ChildrenID_name) = ID_old2new(ChildrenIDs)';
end
end