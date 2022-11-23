function [C, df, df_err] = corrdim2(nodes_coordinates, rad, varargin)
% Function to calculate the correlation dimension of set of points.
%
% Input
% nodes_coordinates = coordinates of the set of nodes.
% rad = radius at which the correlation integral (C) is evaluated and
%       fitted to determine the correlation dimension.
%
% Output
% C = array giving the correlation integral at each of the input radius (rad)
% df = correlation dimension
% df_err = error of the correlation dimension inferred from non-linear fit.
%
% In the 'standard' method, all pairwise distances between the nodes are
% calculated and binned and the bin counts are cumulated to calculate the
% correlation integral (C). This method uses a MEX function corrdim2mex.c
% that must be compiled.
% 
% Two alternative methods are supported to calculate the fractal dimension.
% 'periodic_boundaries': Use periodic boundary conditions when calculating 
%                        the correlation integral. This method also
%                        requires the MEX function corrdim2mex.c.
% 'mask': Consider only points that are a distance R or larger from the
%         rectangular boundary of the shape when calculating the correlation
%         integral at radius R.
%% Parse inputs
rad_size = size(rad);
rad = rad(:);

% Calculate the radius of gyration and center of mass.
[Rg,R_COM] = rg(nodes_coordinates);
%% Parse optional parameters
p = inputParser;
addParameter(p, 'R_min', Rg/10); % Lower radius limit of the C vs rad fit.
addParameter(p, 'R_max', Rg); % Upper radius limit of the C vs rad fit.
addParameter(p, 'Method', 'periodic_boundaries',@(x) ismember(x,{'standard','periodic_boundaries','mask'}));
addParameter(p, 'Boundaries_extent', []); % (Deprecated) Extension of the boundary, in each dimension, expressed as a ratio of the range of the data.
addParameter(p, 'Boundary', 1, @(x) (isnumeric(x) && x>=0 && x<=1) || strcmp(x,'uni')); % Type of boundary used for the periodic boundaries method. By default the smallest boundary that encloses all points is used.
addParameter(p, 'Plot', nargout==0); % Plot the fit.
parse(p, varargin{:});
options = p.Results;

% Send warning if the Boundaries_extent option is used.
if ~isempty(options.Boundaries_extent)
    warning('Boundaries_extent is deprecated. Use "Boundary" instead.');
end
%% Center and rotate the nodes coordinates
% This step is not necessary for some of the methods since pairwise
% distances are not affected by translation or rotation.
N_points = size(nodes_coordinates, 1);
bin_rad = [0;rad];

% Center the points at the center of mass and calculate the PCA vectors of the positions.
nodes_coordinates = bsxfun(@minus,nodes_coordinates,R_COM);
PCA_vecs = pca(nodes_coordinates);

% Rotate the nodes positions to align the x,y axes with the PCA vectors.
% This rotation aligns the tree with its symmetry axes. This is useful when
% using periodic boundary conditions.
rotation_angle = pi/2 - atan2(PCA_vecs(2,1),PCA_vecs(1,1));
rotation_matrix  = [cos(rotation_angle) sin(rotation_angle); -sin(rotation_angle) cos(rotation_angle)];
nodes_coordinates = nodes_coordinates*rotation_matrix;

% Calculate the position of the boundaries.
%Boundaries_position = minmax(nodes_coordinates')';
if isnumeric(options.Boundary)
    [~,Boundaries_position] = calculate_tree_size(nodes_coordinates,'Method','percentile','Extent_ratio',options.Boundary,'Rotation',0);
elseif strcmp(options.Boundary,'uni')
    [~,Boundaries_position] = calculate_tree_size(nodes_coordinates,'Method','uniform','Rotation',0);
else
    error('Boundary option not supported');
end

% Remove points that are outside the boundary.
is_node_inside_boundary = nodes_coordinates(:,1) >= Boundaries_position(1,1) ...
    & nodes_coordinates(:,1) <= Boundaries_position(2,1) ...
    & nodes_coordinates(:,2) >= Boundaries_position(1,2) ...
    & nodes_coordinates(:,2) <= Boundaries_position(2,2);
if any(~is_node_inside_boundary)
    nodes_coordinates = nodes_coordinates(is_node_inside_boundary,:);
end
%% Calculate the corelation integral (C) at each input radius (R).
switch options.Method
    case 'standard'
        bincounts = corrdim2mex(nodes_coordinates', rad);
        assert(all(bincounts >= 0),'Bin counts are negative.')
        assert(sum(bincounts) <= N_points*(N_points-1)/2,'The total number of pairwise distances is higher than expected.')
        C = cumsum(bincounts)/(N_points*(N_points-1)/2);
    case 'mask'
        % Calculate the distance to the boundary for each node. 
        % In 2D, this is the shortest distance to any of the 4 boundaries.
        nodes_boundary_dist = min([bsxfun(@minus,nodes_coordinates,Boundaries_position(1,:)) bsxfun(@minus,Boundaries_position(2,:),nodes_coordinates)],[],2);
        boundary_dist_counts = histc(nodes_boundary_dist,bin_rad);
        
        % Calculate the number of points that at least a distance R from
        % any of the boundaries.
        N_points_R = N_points - cumsum(boundary_dist_counts);
        
        % For each point X, calculate all pairwise distances with other
        % points Y up to the shortest boundary distance. This essentially
        % masks the point X once the radius R is bigger than its distance to
        % the boundary, because all pairwise distances formed with other points Y 
        % are not counted if they are bigger than the distance to the boundary.
        bincounts2 = zeros(size(rad));
        
        % Use a KDtree to accelerate distance calculation and masking.
        KDTree = KDTreeSearcher(nodes_coordinates);
        for i=1:N_points
            %
            [~,dists] = rangesearch(KDTree,nodes_coordinates(i,:),nodes_boundary_dist(i),'SortIndices',false);
            Counts = histc(dists{1},bin_rad);
            
            % dists{1} will always contain a distance in its first entry 
            % due to the fact that the point i also appears in the KDtree.
            % To correct this, simply remove one count from the lowest bin.
            Counts(1) = Counts(1) - 1;
            bincounts2 = bincounts2 + Counts(1:end-1)';
        end
        C = cumsum(bincounts2)/N_points./N_points_R(1:end-1);
    case 'periodic_boundaries'
        bincounts = corrdim2mex(nodes_coordinates', rad, Boundaries_position);
        
        assert(all(bincounts >= 0),'Bin counts are negative.')
        assert(sum(bincounts) <= N_points*(N_points-1)/2,'The total number of pairwise distances is higher than expected.')
        C = cumsum(bincounts)/(N_points*(N_points-1)/2);
end
%% Plot the points and the position of the boundaries.
if options.Plot
    figure; hold on;
    title('Nodes position');
    a = gca; axis equal;
    Nodes_h = plot(nodes_coordinates(:,1),nodes_coordinates(:,2),'o');
    
    Boundary_endpoints = [
        Boundaries_position(1,1),Boundaries_position(1,2);
        Boundaries_position(1,1),Boundaries_position(2,2);
        Boundaries_position(2,1),Boundaries_position(2,2);
        Boundaries_position(2,1),Boundaries_position(1,2);
        Boundaries_position(1,1),Boundaries_position(1,2);
        ];
    boundary_h = plot(Boundary_endpoints(:,1),Boundary_endpoints(:,2),'-');
    xlim(a.XLim + 0.2*diff(a.XLim)*[-1,1]);
    ylim(a.YLim + 0.2*diff(a.YLim)*[-1,1]);
    legend([Nodes_h,boundary_h],{'Nodes','Boundary'})
end

% Reshape rad and C according to the original rad size.
C = reshape(C, rad_size);

% Fit the correlation function with a power law to find fractal dimension
[df, df_err] = corrdim_fit(rad, C, 'R_min', options.R_min, 'R_max', options.R_max, 'Plot', options.Plot);
end