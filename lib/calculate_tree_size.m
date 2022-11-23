function [Sizes,Boundary_positions,f] = calculate_tree_size(varargin)
% Function to calculate the size of a 2D tree.
% The size is calculated by fitting a rectangle to the set of nodes' position using PCA.

% Call options
% calculate_tree_size(Positions,...)
% Positions = n x 2 array corresponding to the x,y coordinates of a set of points.

% calculate_tree_size(struct,...)
% struct = tree structure that contains the branches.

% Output
% Boundary_positions = x,y positions of the 4 boundaries.
% Boundary_positions(:,1) corresponds to the lower and upper x limits of
% the vertical boundaries. If rotation is allowed, the boundary positions
% are given in the rotated coordinate system.
%
% History
% August 12, 2022
% Ensure that the intrinsic neuron coordinate system is right-handed.
% Depending on the input, the output orthogonal vectors of pca may sometimes by
% left-handed for 2D data. If the coordinate system is left-handed, 
% the neuron coordinate system is flipped along the long axis (Axes_vec(:,1) = -Axes_vec(:,1)).
%% Parse Inputs
if isnumeric(varargin{1})
    Positions = varargin{1};
elseif isa(varargin{1},'struct')
    Tree = varargin{1};
    Positions = {Tree.PointsPos};
    Positions = cell2mat(cellfun(@(x) x(2:end,:),Positions(:),'uni',0));
    Positions = [Tree(1).PointsPos(1,:); Positions];
else
    error('The input must be a struct or an array');
end
assert(size(Positions,2) == 2,'Positions must be a list of 2D row vectors.');
varargin(1) = [];
%% Parse optional parameters
p = inputParser;
addParameter(p, 'Method', 'percentile',@(x) ismember(x,{'percentile','uniform'})); % Method used to calculate the size in each dimension.
addParameter(p, 'Extent_ratio', 0.99); % Central proportion of the mass distribution contained in each dimension.
addParameter(p, 'AxesLabels', {'Axis 1','Axis 2'}); % Labels of axes. The first axis corresponds to the longest dimension when rotation=true.
addParameter(p, 'Plot', nargout == 0); % Plot the fitted rectangle.
addParameter(p, 'Plot_align_axes', false); % Align first PCA axis along the y-axis.
addParameter(p, 'Rotation', true); % Allow rotation of the rectangle.
addParameter(p, 'COM', []); % Center of mass used to calculate the sizes.
parse(p, varargin{:});
options = p.Results;
%%
% Center the positions to the center of mass.
if isempty(options.COM)
    COM = mean(Positions);
else
    COM = options.COM;
end
Positions_centered = bsxfun(@minus,Positions,COM);

% Find the symmetry axes of the shape.
% Each column of Axes_vec represents one axis vector.
if options.Rotation
    % Use PCA to find the symmetry axes using the centered positions.
    % 'Centered' is false in the pca routine to prevent recentering. This
    % may cause differences if the center of mass is not equal to
    % mean(Positions).
    [Axes_vec,Positions_axes_projection] = pca(Positions_centered,'Centered',false);
    
    % The PCA vectors are ordered in decreasing variance. Order the axis
    % vectors in increasing variance instead.
    Axes_vec = flip(Axes_vec,2);
    Positions_axes_projection = flip(Positions_axes_projection,2);
    
    % Ensure that the axes vector output by pca corresponds to a right-handed
    % coordinate system. If not, flip the coordinate system along the
    % long-axis.
    Axes_vec_cross = round(det(Axes_vec)); % will equal 1 or -1
    if Axes_vec_cross == -1
        Axes_vec(:,1) = -Axes_vec(:,1);
        Positions_axes_projection(:,1) = -Positions_axes_projection(:,1);
    end
else
    Axes_vec = [[1;0] [0;1]];
    Positions_axes_projection = Positions_centered;
end

switch options.Method
    case 'percentile'
        % Calculate the size of the rectangle by finding the position of the low and
        % high percentiles along each dimension. The low percentile is
        % options.Extent_ratio/2 and the high percentile is
        % 1-options.Extent_ratio/2. This way, a  maximum of
        % options.Extent_ratio*100 % of the points are enclosed by the rectangle in each dimension.
        Boundary_percentiles = [0 100] + (1-options.Extent_ratio)*100/2 * [1, -1];
        Boundary_positions = prctile(Positions_axes_projection,Boundary_percentiles);
        Sizes = sort(diff(Boundary_positions));
        
    case 'uniform'
        % Calculate the size  assuming that the node positions are
        % uniformly distributed. In such case, we expect the variance the of
        % positions to equal D.^2/12 where D is the size.
        Sizes = sqrt(12*var(Positions_axes_projection));
        Boundary_positions = [-1; 1].*Sizes/2;
end
%% Plot the rectangle and the PCA vectors.
if options.Plot
    f = figure;
    
    % Rotate the points and axes for plotting.
    if options.Plot_align_axes
        % Align the shortest dimension along the x axis.
        rotation_angle = 0 - atan2(Axes_vec(2,1),Axes_vec(1,1));
        rotation_matrix  = [cos(rotation_angle) sin(rotation_angle); -sin(rotation_angle) cos(rotation_angle)];
        
        Positions_rotated = Positions_centered * rotation_matrix;
        Positions = bsxfun(@plus,Positions_rotated,COM);
        Axes_vec = (Axes_vec'*rotation_matrix)';
    end
    
    % Plot the points and color code according to their projections on the
    % PCA vectors.
%     [27,158,119;
%     217,95,2;
%     117,112,179;
%     231,41,138];

    Colors = [Positions_axes_projection(:,1) > 0, zeros(size(Positions_axes_projection,1),1), Positions_axes_projection(:,2) > 0];
    scatter(Positions(:,1),Positions(:,2),1,Colors,'.');
    axis equal
    
    % Plot the axes vectors.
    hold on;
    axis_1_pos = bsxfun(@plus,bsxfun(@times,Boundary_positions(:,1),Axes_vec(:,1)'),COM);
    axis_2_pos = bsxfun(@plus,bsxfun(@times,Boundary_positions(:,2),Axes_vec(:,2)'),COM);
    axis1_h = plot(axis_1_pos(:,1),axis_1_pos(:,2),'r--','DisplayName',options.AxesLabels{1});
    axis2_h = plot(axis_2_pos(:,1),axis_2_pos(:,2),'b--','DisplayName',options.AxesLabels{2});
    
    % Plot the rectangle edges.
    Boundary_axes_projections = [repmat(Boundary_positions(:,1),[2,1]) repelem(Boundary_positions(:,2),2)];
    Rectangle_vertices = Boundary_axes_projections * Axes_vec';
    
    % Calculate the angle of the vertices and sort them by their angle.
    angles = atan2(Rectangle_vertices(:,2),Rectangle_vertices(:,1));
    [~,sorted_ind] = sort(angles);
    Rectangle_vertices = Rectangle_vertices(sorted_ind,:);
    
    % Add COM to the vertices and plot.
    Rectangle_vertices = [Rectangle_vertices; Rectangle_vertices(1,:)]; % Close the rectangle.
    Rectangle_vertices = bsxfun(@plus,Rectangle_vertices,COM);
    boundary_h = plot(Rectangle_vertices(:,1),Rectangle_vertices(:,2),'-','Color',[56,108,176]/255,'DisplayName','Boundary');
    
    xlabel('Anterior-posterior position (\mum)');
    ylabel('Left-right position (\mum)');
    legend([axis1_h,axis2_h,boundary_h],[options.AxesLabels,{'Boundary'}],'Location','Northwest');
else
    f = [];
end

% Recenter the boundary positions with respect to the center of mass, if
% no rotation were performed.
Boundary_positions = bsxfun(@plus,Boundary_positions,COM(:)');
% if ~options.Rotation
%     Boundary_positions = bsxfun(@plus,Boundary_positions,COM(:)');
% end
end