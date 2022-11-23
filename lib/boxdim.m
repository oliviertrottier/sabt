function [df, df_err, N_boxes, Boxsizes] = boxdim(Nodes_pos, varargin)
% Function to calculate the box fractal dimension also know as the
% Haussdorff dimension.
%
% Call options:
% boxdim(Nodes_pos,...) = calculates box dimension using a default set of box sizes (see below)
% boxdim(Nodes_pos,Boxsizes...) = calculates box dimension using a given set of box sizes.
%
% Input
% Nodes_pos = Nx2 array representing the 2D position of tree nodes.
% Boxsizes = Sizes of box used to cover the shape.
%
% Output
% df = box fractal dimension
% df_err = error of box fractal dimension calculated from non-linear fit
% N_boxes = number of boxes needed to cover the shape for each input Boxsizes.
% Boxsizes = Box sizes used to query the build the coverings.
%% Parse input
[Rg,R_com] = rg(Nodes_pos);
if isnumeric(varargin{1})
    Boxsizes = reshape(varargin{1},[],1);
else
    % Delay definition until R_min and R_max are defined.
    Boxsizes = [];
end
%% Parse optional parameters
p = inputParser;
addParameter(p, 'R_min', Rg/20); % Lower radius limit of the N vs Boxsize fit.
addParameter(p, 'R_max', Rg/2); % Upper radius limit of the N vs Boxsize fit.
addParameter(p, 'Periodic_boundaries', true); % Use periodic boundary conditions
addParameter(p, 'Plot', nargout==0); % Plot the fit.
addParameter(p, 'Translations', true); % Translate each covering of a given size to find a more optimal covering.
parse(p, varargin{:});
options = p.Results;
%%
% Use a KD tree to accelerate nearest neighbor searches.
KDTree_chebychev = KDTreeSearcher(Nodes_pos,'Distance','chebychev');

% Define default box sizes that will be used to cover the shape.
if isempty(Boxsizes)
    % Find the smallest distance between any two points in the dataset.
    KDTree_seuclidean = KDTreeSearcher(Nodes_pos,'Distance','euclidean');
    [~,smallestdists] = knnsearch(KDTree_seuclidean,Nodes_pos,'K',2);
    Min_dist = min(smallestdists(:,2));
    Boxsizes = logspace(log10(min(options.R_min,Min_dist)),log10(max(options.R_max,2*Rg)),50)';
end

N_boxsizes = numel(Boxsizes);
Nodes_pos_minmax = minmax(Nodes_pos')';

% Initialize the translations of the box centers. The translations are
% given in units of the box sizes.
N_translations = zeros(N_boxsizes,1);
if options.Translations
    Translations = cell(N_boxsizes,1);
    max_translation_size = Rg/10;
    min_translation_size = Rg/50;
    
    % Extend the boundary of the shape.
    Nodes_pos_minmax = Nodes_pos_minmax + max_translation_size*[-1 -1;1 1];
    
    for i = 1:N_boxsizes
        Translation_stepsize = min(max(1/10,min_translation_size/Boxsizes(i)),max_translation_size/Boxsizes(i));
        Translations_x = -.5:Translation_stepsize:.5;
        Translations_y = -.5:Translation_stepsize:.5;
        N_translations(i) = numel(Translations_x)*numel(Translations_y);
        
        [Translations_xgrid,Translations_ygrid] = ndgrid(Translations_x,Translations_y);
        Translations{i} = [Translations_xgrid(:) Translations_ygrid(:)];
    end
    
    if any(N_translations > 10000)
        warning('The number of translations is higher than 10000 for some box sizes. This may not be efficient.');
    end
end

% Extend the size of the box bin boundary slightly to make sure that all points
% are covered.
Nodes_pos_minmax = Nodes_pos_minmax + abs(Nodes_pos_minmax)*0.01.*[-1 -1;1 1];

% Count the number of boxes needed to cover the shape for each box size.
N_boxes = zeros(size(Boxsizes));
for i = 1:N_boxsizes
    % Define the box centers by centering them on the center of mass.
    boxcenters_x = (floor((Nodes_pos_minmax(1,1)-R_com(1))/Boxsizes(i)):ceil((Nodes_pos_minmax(2,1)-R_com(1))/Boxsizes(i)))*Boxsizes(i) + R_com(1);
    boxcenters_y = (floor((Nodes_pos_minmax(1,2)-R_com(2))/Boxsizes(i)):ceil((Nodes_pos_minmax(2,2)-R_com(1))/Boxsizes(i)))*Boxsizes(i) + R_com(2);
    [boxcenters_xgrid, boxcenters_ygrid] = ndgrid(boxcenters_x,boxcenters_y);
    box_centers = [boxcenters_xgrid(:) boxcenters_ygrid(:)];
    
    % For each box center, find the Cehbyshev distance of the closest point in the dataset.
    %box_closest_ChebyDist = pdist2(nodes_coordinates, box_centers, 'chebychev', 'Smallest', 1);
    [~,box_closest_ChebyDist] = knnsearch(KDTree_chebychev,box_centers,'K',1);
    
    % Initialize the number of boxes needed to cover.
    N_boxes(i) = nnz(box_closest_ChebyDist < Boxsizes(i)/2);
    
    % Attempt to find a better covering by translating the box centers.
    if N_translations(i) > 0
        N_boxes_temp = zeros(N_translations(i),1);
        Translations_temp = Translations{i};
        for j = 1:N_translations(i)
            box_centers_translated = bsxfun(@plus,box_centers,Translations_temp(j,:)*Boxsizes(i));
            [~,box_closest_ChebyDist] = knnsearch(KDTree_chebychev,box_centers_translated,'K',1);
            N_boxes_temp(j) = nnz(box_closest_ChebyDist < Boxsizes(i)/2);
        end
        
        % Find the minimal number of boxes needed to cover the points.
        N_boxes(i) = min([N_boxes(i); N_boxes_temp]);
    else
        N_boxes(i) = nnz(box_closest_ChebyDist < Boxsizes(i)/2);
    end
end
%% Fit the number of boxes with a power law to find the fractal dimension.
[df, df_err] = boxdim_fit(Boxsizes, N_boxes, 'R_min', options.R_min, 'R_max', options.R_max, 'Plot', options.Plot);
end