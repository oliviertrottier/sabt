function [Rg,COM] = rg(varargin)
% Function to calculate the radius of gyration and the center of mass of a distribution of points.
%
% Call options
% rg(R), R = position vectors.
% rg(occ,x,y), occ = occupancy matrix, x,y = lattice points coordinate meshgrids.
% rg(Tree), Tree = tree structure where the nodes position of branch i are stored in Tree(i).PointsPos.
%
% Output
% Rg = radius of gyration
% COM = position of center of mass
%% Determine the format of the input.

if numel(varargin)==1
    if isa(varargin{1},'struct')
        Tree = varargin{1};
        Branches_nodes_cell = {Tree(:).PointsPos};
        
        % Collect all branch nodes, omitting branch points. Branchpoints 
        % correspond to the parent's last node.
        Branches_nodes_cell = cellfun(@(x) x(2:end,:),Branches_nodes_cell,'uni',0);
        R = cat(1,Branches_nodes_cell{:});
        
        % Add soma node position.
        R = [Tree(1).PointsPos(1,:); R];
    else
        R = varargin{1};
    end
elseif numel(varargin)==3
    occ = varargin{1};
    x = varargin{2};
    y = varargin{3};
    R = [x(occ) y(occ)];
end
%% Calculate the center of mass and radius of gyration
N_R = size(R,1);
COM = mean(R);
R_centered = (R-repmat(COM,N_R,1)).^2;
Rg = sqrt(sum(R_centered(:))/N_R);
end