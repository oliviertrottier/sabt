function Tree = rotate_tree(Tree,Angle,varargin)
% Function to rotate a tree.

% Input
% Tree = tree structure
% Angle = angle of rotation in degrees.
%% Parse optional parameters
p = inputParser;
addParameter(p, 'Center', Tree(1).PointsPos(1,:)); % Center of rotation. The default center of rotation is set to the first node of the first branch.
parse(p, varargin{:});
options = p.Results;
%% Rotate the branches with respect to the center of rotation.
rotation_matrix = [cosd(Angle) sind(Angle); -sind(Angle) cosd(Angle)];
N_branches = numel(Tree);
for i=1:N_branches
    Tree(i).PointsPos = (Tree(i).PointsPos - options.Center)*rotation_matrix + options.Center;
end