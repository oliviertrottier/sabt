function [CH_area, CH_diam, CH_DT] = convexhull_analysis(varargin)
% Function to analyze the convex hull of an input tree structure.
%
% Call options
% convexhull_analysis(Nodes_pos,...)   = calculates the convex hull given a
%                                        set of nodes positions. Nodes_pos is
%                                        an nx2 array.
% convexhull_analysis(Tree_struct,...) = calculates the convex hull from a
%                                        given tree structure. The nodes
%                                        position are defined by the field
%                                        "PointsPos".
%
% Output
%
% CH_area = area of the convex hull
% CH_diam = diameter of the convex hull
% CH_DT = Delaunay triangulation structure associated with the convex hull
%% Parse parameters
if isa(varargin{1},'struct')
    Tree_struct = varargin{1}(:);
    Nodes_pos = cell2mat({Tree_struct.PointsPos}');
else
    Nodes_pos = varargin{1};
end
varargin(1) = [];

p = inputParser;
addParameter(p,'Lengthscale',1); % Scale all tree nodes positions.
parse(p,varargin{:});
options = p.Results;
%% Calculate the convex hull and construct Delaunay triangulation.
% Scale positions with the lengthscale.
Nodes_pos = Nodes_pos*options.Lengthscale;

% Calculate the Delaunay triangulation out of the convex hull points.
CH_ind = convhull(Nodes_pos);
CH_DT = delaunayTriangulation(Nodes_pos(CH_ind(1:end-1),:));

% Calculate the area and the diameter.
CH_area = delaunayArea(CH_DT);
CH_diam = max(pdist(CH_DT.Points));
end