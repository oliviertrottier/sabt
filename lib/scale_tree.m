function Tree = scale_tree(Tree,scale)
% Function to the scale the positions of a tree structure.
%
% Input
% Tree = structure of the tree to scale
% scale = scalar used to scale the positions and diameters of the tree.
N_branches = numel(Tree);
for i=1:N_branches
    Tree(i).PointsPos = Tree(i).PointsPos * scale;
end

% Scale the dimensionful length.
if isfield(Tree,'Length_dim')
    for i=1:N_branches
        Tree(i).Length_dim = Tree(i).Length_dim * scale;
    end
end

% Scale the nodes' diameter.
if isfield(Tree,'PointsDiameter')
    for i=1:N_branches
        Tree(i).PointsDiameter = Tree(i).PointsDiameter * scale;
    end
end
end