function [Hitting_prob, Boxsizes, Meshsize, varargout] = hitprob(Tree_nodes_position, varargin)
% Function to calculate the hitting probability of a set of Tree nodes in 2D.
%
% occ is a occupied points on a lattice. 0=unoccupied, 1=occupied.
% latticex,latticey = mesh coordinates of the lattice.
% Boxsizes = sizes of the boxes
% hitting_prob = prob. of being greater or equal than the critical value (crit).

% Call options
% hitprob(Tree_nodes_position):
%   calculates hitting probability from a set of nodes position.
% hitprob(Tree_nodes_position, Boxsizes):
%   calculates hitting probability for a given set of box sizes.
% hitprob(Tree_nodes_position, N_boxsizes):
%   calculates hitting probability for a given number of boxes. The box
%   sizes are inferred from the radius of gyration.
%
% Input
% Tree_nodes_position = Nx2 array representing the 2D positions of the tree
%                       nodes.
%
% Output
% Hitting_prob = 1xm array representing the probability that a box of size
%                Boxsizes(i) hits the set of nodes.
% Boxsizes = size of boxes used to calculate the hitting probability.
% Meshsize = approximate size of box that hits the structure 50% of the time.
% varargout{1} = Box centers used (if options.Method=3).
% varargout{2} = Number of hits of each box (if options.Method=3).
%% Parse inputs
if nargin > 1 && isnumeric(varargin{1})
    if numel(varargin{1}) == 1
        N_boxsizes = varargin{1};
        Boxsizes = [];
    else
        Boxsizes = varargin{1};
        N_boxsizes = numel(Boxsizes);
    end
    varargin(1) = [];
else
    % Use 20 box sizes by default.
    N_boxsizes = 20;
    Boxsizes = [];
end
varargout{1} = [];
%% Parse optional parameters
p = inputParser;
addParameter(p, 'Plot', nargout == 0); % Plot the hitting probability curve.
addParameter(p, 'PlotIterations', false); % Plot each iteration of calculating the box hitting probability. Only applicable to certain methods.
addParameter(p, 'Method', 3); % Method used to calculate the hitting probability. See below.
addParameter(p, 'Mask', 'rg',@(x) ismember(x,{'rg','D_uni','convexhull','','none'})); % Remove boxes outside the mask (convex hull or radius of gyration).
addParameter(p, 'Rotate', true); % Rotate the tree to align its symmetry axes (calculated with PCA) with the cartesian plane.
parse(p, varargin{:});
options = p.Results;
%%
% Calculate the radius of gyration and center the tree positions at the center-of-mass.
[Rg,R_com] = rg(Tree_nodes_position);
Tree_nodes_position = bsxfun(@minus, Tree_nodes_position, R_com);

% Rotate the nodes positions to align the x,y axes with the PCA vectors.
% This rotation aligns the tree with its symmetry axes.
if options.Rotate
    PCA_vecs = pca(Tree_nodes_position);
    rotation_angle = pi/2 - atan2(PCA_vecs(2,1),PCA_vecs(1,1));
    rotation_matrix  = [cos(rotation_angle) sin(rotation_angle); -sin(rotation_angle) cos(rotation_angle)];
    Tree_nodes_position = Tree_nodes_position*rotation_matrix;
end

% Calculate the maximal box size.
X_minmax = minmax(Tree_nodes_position(:, 1)')';
Y_minmax = minmax(Tree_nodes_position(:, 2)')';
Tree_bounding_box = [X_minmax,Y_minmax];
Tree_range_max = max(abs(Tree_bounding_box(:)));

% Calculate the smallest distance between two points positions.
[~, Points_smallest_dist] = knnsearch(Tree_nodes_position, Tree_nodes_position, 'K', 2);
Points_smallest_dist_min = min(Points_smallest_dist(:, 2));
Boxsize_min = Points_smallest_dist_min;
%% First method: rectangular array of uniform square boxes.
% The advantage of this method is that the array of boxes preserve the
% aspect ratio of the shape. However, for bigger boxes, the array can
% sometimes be too big and the hitting probability is then underestimated.
% Also, the size of Boxsizes may be inappropriate for the size of the
% array.
if options.Method == 1
    % Initialize output
    N_boxes_total = nan(4, N_boxsizes);
    N_non_empty_boxes = nan(4, N_boxsizes);
    
    % Define the default box sizes
    if isempty(Boxsizes)
        Boxsizes = logspace(log10(Boxsize_min), log10(Tree_range_max), N_boxsizes);
    end
    
    for i = 1:N_boxsizes
        % Ceil the tree range to a multiple of the box size.
        Boxsize = Boxsizes(i);
        
        % Define multiple sets of boxes to average the hitting probability over all of them.
        % For each array, one of the 4 corners of the tree bounding box
        % perfectly fits with the respective bin corner.
        Bottomleftbin_pos(1,:) = round([X_minmax(1) Y_minmax(1)]/Boxsize)*Boxsize;
        Bottomleftbin_pos = round(Bottomleftbin_pos/Boxsize)*Boxsize;
        for j=1:size(Bottomleftbin_pos,1)
            % Determine the 2d bin indices of all points. The origin of the bin
            % index is set to the bottom left corner.
            Bottomleftbin_pos_temp = Bottomleftbin_pos(j,:);
            Points_boxes_ind = ceil(bsxfun(@minus, Tree_nodes_position, Bottomleftbin_pos_temp)/Boxsize);
            
            % Include points that are right on the left or bottom boundary.
            Points_boxes_ind(Points_boxes_ind(:,1)==0,1) = 1;
            Points_boxes_ind(Points_boxes_ind(:,2)==0,2) = 1;
            
            % Determine the index of the topright-most bin.
            % Do not count bins that are completely filled in the tree's
            % bounding box.
            Topright_bin_ind = round(([X_minmax(2) Y_minmax(2)] - Bottomleftbin_pos_temp)/Boxsize);
            
            % Remove points that are outside the bins.
            Points_boxes_ind = Points_boxes_ind(all(bsxfun(@le,Points_boxes_ind,Topright_bin_ind),2),:);
            Points_boxes_ind = Points_boxes_ind(all(bsxfun(@ge,Points_boxes_ind,[1 1]),2),:);
            
            % Determine the total number of boxes in each dimension.
            N_boxes_x = Topright_bin_ind(1);
            N_boxes_y = Topright_bin_ind(2);
            Box_counts = accumarray(Points_boxes_ind, 1, [N_boxes_x, N_boxes_y]);
            
            if options.PlotIterations && i >= 3/4*N_boxsizes
                Box_edges_x = (0:N_boxes_x)*Boxsize + Bottomleftbin_pos_temp(1);
                Box_edges_y = (0:N_boxes_y)*Boxsize + Bottomleftbin_pos_temp(2);
                %Box_edges_x = (Bottomleftbin_pos(1)/boxsize:ceil(X_minmax(2)/boxsize))*boxsize;
                %Box_edges_y = (Bottomleftbin_pos(2)/boxsize:ceil(Y_minmax(2)/boxsize))*boxsize;
                plot_boxes(Tree_nodes_position, Box_edges_x, Box_edges_y, Box_counts)
            end
            
            % Calculate the hitting probability.
            N_boxes_total(j,i) = numel(Box_counts);
            N_non_empty_boxes(j,i) = nnz(Box_counts);
        end
    end
    
    % Calculate hitting probability by counting the number of non-empty
    % boxes.
    Hitting_prob = sum(N_non_empty_boxes(1:size(Bottomleftbin_pos,1),:),1)./sum(N_boxes_total(1:size(Bottomleftbin_pos,1),:),1);
end
%% Second method: square array where each bin is divided into 4 at every iteration.
% Calculate minimal box sizes until the total N_octaves is greater than
% N_boxsizes.
if options.Method == 2
    % Calculate the number of octaves to span the range of sizes of the
    % boxes.
    N_octaves = floor(log(Tree_range_max./Boxsize_min)/log(2));
    N_points_per_octave = ceil(N_boxsizes/N_octaves);
    Boxsize_maxima = Tree_range_max.*2.^((0:N_points_per_octave-1)/N_points_per_octave);
    
    % Define the default box sizes.
    if isempty(Boxsizes)
        Boxsizes = bsxfun(@times, Boxsize_maxima, (1/2).^(N_octaves-1:-1:0)');
    end
    
    % Initialize output
    Hitting_prob = nan(size(Boxsizes));
    N_boxes_total = nan(size(Boxsizes));
    
    % Filter used for convolution. See below.
    filter = ones(2);
    
    for i = 1:N_points_per_octave
        for j = 1:N_octaves
            Boxsize = Boxsizes(j, i);
            if j == 1
                % Calculate the bin ind of the points.
                Bottomleftbin_pos = floor(-Boxsize_maxima(i)*ones(1, 2)/Boxsize)*Boxsize;
                Points_boxes_ind = ceil(bsxfun(@minus, Tree_nodes_position, Bottomleftbin_pos)/Boxsize);
                
                % Determine the total number of bins in each dimension.
                N_boxes_x = 2*Boxsize_maxima(i)/Boxsize;
                N_boxes_y = N_boxes_x;
                Box_counts = accumarray(Points_boxes_ind, 1, [N_boxes_x, N_boxes_y]);
            else
                % Calculate the boxcounts of bigger boxes by convolving the Box_counts array.
                Box_counts = conv2(Box_counts, filter, 'same');
                Box_counts = Box_counts(1:2:end, 1:2:end);
            end
            
            % Calculate the hitting probability.
            N_boxes_total(j, i) = numel(Box_counts);
            Hitting_prob(j, i) = nnz(Box_counts)/N_boxes_total(j, i);
            
            % Plot bins with the tree positions.
            if options.PlotIterations
                Box_edges_x = (Bottomleftbin_pos(1)/Boxsize:ceil(Boxsize_maxima(i)/Boxsize))*Boxsize;
                Box_edges_y = (Bottomleftbin_pos(2)/Boxsize:ceil(Boxsize_maxima(i)/Boxsize))*Boxsize;
                plot_boxes(Tree_nodes_position, Box_edges_x, Box_edges_y, Box_counts)
            end
        end
    end
end
%% Third method: overlapping boxes on a square grid.
if options.Method == 3
    % Define the range of the boxes in each dimension.
    switch options.Mask
        case 'D_uni'
            % Define the boxes inside the uniform boundary.
            [~,Boxes_limits] = calculate_tree_size(Tree_nodes_position,'Method','uniform','Rotation',0);
        case 'rg'
            % Define the boxes inside the radius of gyration for each
            % dimension.
            Boxes_limits = Rg*[-ones(1,2);ones(1,2)];
        case {'','none','convexhull'}
            % Use the full extension of the tree for the ranges of the box.
            Boxes_limits = Tree_bounding_box;
    end
    Boxes_range = diff(Boxes_limits);
    
    % Define the smallest boxsize such that there is at least 500 boxes in each dimension.
    N_centers = 500;
    Boxsize = min(Boxes_range'/N_centers);
    
    % Round the box limits in units of the boxsize. This will center the
    % lattice of boxes.
    Boxes_limits = [ceil(Boxes_limits(1,:)/Boxsize); floor(Boxes_limits(2,:)/Boxsize)]*Boxsize;
    
    % Define the box centers.
    Box_centers_x = (Boxes_limits(1,1)+Boxsize/2):Boxsize:(Boxes_limits(2,1)-Boxsize/2);
    Box_centers_y = (Boxes_limits(1,2)+Boxsize/2):Boxsize:(Boxes_limits(2,2)-Boxsize/2);
    [Box_centers_xmesh, Box_centers_ymesh] = meshgrid(Box_centers_x, Box_centers_y);
    Box_centers = [Box_centers_xmesh(:), Box_centers_ymesh(:)];
    Box_lattice_size = size(Box_centers_xmesh);
    
    % Remove box centers that are outside the convex hull.
    switch options.Mask
        case 'convexhull'
            Convexhull_indices = convhull(Tree_nodes_position(:,1),Tree_nodes_position(:,2));
            Convexhull_vertices = Tree_nodes_position(Convexhull_indices,:);
            is_box_inside_mask = inpolygon(Box_centers(:,1),Box_centers(:,2),Convexhull_vertices(:,1),Convexhull_vertices(:,2));
            Box_centers = Box_centers(is_box_inside_mask,:);
    end
    
    % For each boxcenter, find the smallest size of the box that hits a tree
    % node (hitting condition).
    Tree_KDTree = KDTreeSearcher(Tree_nodes_position,'Distance','chebychev');
    [Ind,Cheby_Dist] = knnsearch(Tree_KDTree,Box_centers,'k',1);
    %[Cheby_Dist, Ind] = pdist2(Tree_nodes_position, box_centers, 'chebychev', 'Smallest', 1);
    
    % Multiply the chebychev distance by 2 to obtain the smallest box size that hits the tree.
    Boxsize_hits = 2*Cheby_Dist;
    
    % Calculate the default box sizes.
    if isempty(Boxsizes)
        Boxsizes = logspace(floor(log10(0.5*min(Boxsize_hits))),ceil(log10(2*max(Boxsize_hits))),N_boxsizes);
    end
    
    % Calculate the cumulative distribution of the hitting box sizes. This corresponds to
    % the hitting probability curve.
    Boxsize_binedges = [10^(log10(Boxsizes(1)) - log10(Boxsizes(2)/Boxsizes(1))) Boxsizes];
    Hitting_prob = histcounts(Boxsize_hits,Boxsize_binedges,'Normalization','cdf');
    
    % Save secondary outputs.
    varargout{1} = reshape(Box_centers + R_com ,[Box_lattice_size,2]);
    varargout{2} = reshape(Boxsize_hits,Box_lattice_size);
end
%% Format output and calculate mesh size.
Boxsizes = reshape(Boxsizes,[],1);
Hitting_prob = reshape(Hitting_prob,[],1);
Meshsize = calculate_meshsize(Boxsizes,Hitting_prob);
%% Plot the hitting probability curve.
if options.Plot
    % Plot the hitting probability curve.
    f = figure;
    plot(Boxsizes, Hitting_prob, '--ko')
    hold on;
    a = gca;
    
    ylabel('Hitting Probability ($H$)', 'Interpreter', 'latex');
    xlabel('Box Size ($B$) [$\mu$m]', 'Interpreter', 'latex');
    ylim([0, 1]);
    %a.XScale = 'linear';
    a.XScale = 'log';
    a.YTick = 0:0.25:1;
    
    % Add minor ticks.
    a.XMinorTick = 'on';
    a.YMinorTick = 'on';
    a.XAxis.MinorTickValues = a.XTick(2:end)-diff(a.XTick(1:2))/2;
    a.YAxis.MinorTickValues = a.YTick(2:end)-diff(a.YTick(1:2))/2;
    
    % Plot an horizontal dashed red line at a hitting probability of 0.5.
    Meshsize_hitprob = 0.5;
    plot([a.XLim(1),Meshsize], Meshsize_hitprob*ones(1, 2), 'r-');
    
    % Plot a vertical line indicating the boxsize that achieves a hitting
    % probability of 0.5;
    plot(Meshsize*ones(1, 2), [0, Meshsize_hitprob], 'r-');
    
    % Add arrow and textbox indicating where the meshsize is.
    t = text(a, Meshsize+0.1*diff(a.XLim), Meshsize_hitprob/2, sprintf('Meshsize (%.2f \\mum)',Meshsize), 'Color', 'r','HorizontalAlignment','center');
    
    arrowpos = t.Extent(1:2);
    arrowpos(2, :) = [Meshsize, 0];
    drawnow;
    arrowpos = du2nfu(a, arrowpos);
    annotation(f, 'arrow', arrowpos(:, 1), arrowpos(:, 2), 'Color', 'r')
end
end
%% Function to plot the boxes with a set of 2D points.
function plot_boxes(Points_positions, Box_edges_x, Box_edges_y, Box_counts)
% Create the box meshgrids.
[box_x_mesh, box_y_mesh] = ndgrid(Box_edges_x, Box_edges_y);

% Plot the points positions with the mesh grids.
figure;
Box_counts = [Box_counts, Box_counts(:, end)];
Box_counts = [Box_counts; Box_counts(1, :)];
surf_h = surf(box_x_mesh, box_y_mesh, Box_counts, 'FaceAlpha', 0.5);
view(2);
xlabel('x')
ylabel('y')

% Color the boxes based on their box counts.
c=colorbar;
c.Label.String = 'Box counts';
cmap = jet;

% Add white color for empty boxes.
cmap = [ones(1, 3);cmap];
N_colors = size(cmap, 1);
Box_counts_min = 0;
Box_counts_range = max(surf_h.ZData(:))-Box_counts_min;

cmap_index = ceil((surf_h.ZData-Box_counts_min)/Box_counts_range*(N_colors-1))+1;
CData = reshape(cmap(cmap_index(:), :), [size(cmap_index), 3]);
surf_h.CData = CData;
colormap(cmap);
caxis([0, Box_counts_range]+Box_counts_min);

hold on
plot(Points_positions(:, 1), Points_positions(:, 2), 'kx')
axis equal
axis tight
end
%% Function to calculate the mesh size.
function Meshsize = calculate_meshsize(Radius,Hitting_prob)
% Interpolate r and h.
N = max(1000,10*numel(Radius));
Radius_interp = linspace(min(Radius),max(Radius),N);
Hitting_prob_interp = spline(Radius,Hitting_prob,Radius_interp);

% Find r such that h(r) = 0.5;
h_intersect = 0.5;
[~,intersectind] = min(abs(Hitting_prob_interp-h_intersect));
Meshsize = Radius_interp(intersectind);
end