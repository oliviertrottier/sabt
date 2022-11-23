function metrics = analyze_tree(varargin)
% Function to analyze a tree.
%
% Call options
% analyze_tree(tree_filename,options):
%   analyze tree saved in tree_filename
% analyze_tree(treedata_structure,options):
%   analyze tree in treedata_structure. The tree structure must be stored
%   in treedata_structure.struct. Other metadata can also be stored in
%   treedata_structure.
if isa(varargin{1},'char')
    Treedata = load(varargin{1});
elseif isa(varargin{1},'struct')
    Treedata = varargin{1};
else
    error('The first input is not recognized.');
end
varargin(1) = [];
%% Parse the optional parameters
p = inputParser;
addParameter(p, 'misc', true); % Miscellaneous calculations like total branch length, number of branchpoints, etc.
addParameter(p, 'df', false); % Calculate the correlation fractal dimension.
addParameter(p, 'df_box', false); % Calculate the box fractal dimension.
addParameter(p, 'dftime', false); % Calculate the fractal dimension timeseries.
addParameter(p, 'lacunarity', false); % Calculate the lacunarity.
addParameter(p, 'hitprob', false); % Calculate the hitting probability.
addParameter(p, 'hitprobtime', false); % Calculate the hitting probability timeseries.
addParameter(p, 'interbranch_dist', false); % Calculate the interbranch distances by line scanning at different angles.
addParameter(p, 'collisions', false); % Analyze the statistics of the collisions in simulated trees.
addParameter(p, 'length', false); % Analyze the statistics of the branch length.
addParameter(p, 'persistencelength', false); % Calculate the persistence length.
addParameter(p, 'tipdynamics', false); % Analyze the tip dynamics in simulated trees.
addParameter(p, 'sizes', false); % Calculate various tree sizes.
addParameter(p, 'speciesdensity', false); % Calculate the species densities.
addParameter(p, 'tangential_angles', false); % Calculate the branch tangential angles over space.
addParameter(p, 'branching_angles', false); % Calculate the branching angles in the tree.
addParameter(p, 'Metrics', struct()); % Structure used to initialize the output metrics structure.
parse(p, varargin{:});
options = p.Results;
%% Calculate basic metrics
% Initialize the metric structure
metrics = options.Metrics;

% Return an empty structure if the tree is empty.
if numel(Treedata.struct) == 0
    warning('Tree structure has no branches.');
    return;
end

% Define the timescale of the tree (minutes per time frame).
if isfield(Treedata,'timescale')
    % Load the timescale saved in the input structure.
    Timescale = Treedata.timescale;
    Timescale_unit = Treedata.timescale_unit;
elseif isfield(Treedata, 'minpertimeframe')
    % Load the timescale saved in the input structure (old format).
    Timescale = Treedata.minpertimeframe;
    Timescale_unit = 'min';
else
    % Fix the timescale to 5 mins. Used for old tree structures that didn't save
    % the timescale.
    Timescale = 5;
    Timescale_unit = 'min';
end

% Define the lengthscale and timescale of the tree.
if isfield(Treedata, 'lengthscale')
    Lengthscale = Treedata.lengthscale;
    Lengthscale_unit = Treedata.lengthscale_unit;
else
    Lengthscale = 1;
    Lengthscale_unit = 'none';
end
metrics.lengthscale = Lengthscale;
metrics.lengthscale_unit = Lengthscale_unit;
if isfield(Treedata,'timescale')
    metrics.timescale = Timescale;
    metrics.timescale_unit = Timescale_unit;
end

% Save basic tree info.
if isfield(Treedata,'forestid')
    metrics.forestid = Treedata.forestid;
end

if isfield(Treedata,'image_filename')
    metrics.image_filename = Treedata.image_filename;
end

if isfield(Treedata, 'inputnames')
    N_inputs = numel(Treedata.inputnames);
    metrics.inputnames = Treedata.inputnames;
    for j = 1:N_inputs
        if isfield(Treedata,Treedata.inputnames{j})
            metrics.(Treedata.inputnames{j}) = Treedata.(Treedata.inputnames{j});
        end
    end
end

% Save the total number of timesteps in the simulation and the total number
% of branches.
if isfield(Treedata,'MaxT')
    metrics.MaxT = Treedata.MaxT;
end
metrics.NBranches = numel(Treedata.struct);

% Save the simulation parameters.
if isfield(Treedata,'MaxT') || isfield(Treedata,'TotalSimulationTime') || isfield(Treedata,'Start_time')
    Sim_params = rmfield(Treedata,{'struct','Tips','timeseries','Collisions','Boundary_size'});
    Sim_params_metadata = whos('Sim_params');
    if Sim_params_metadata.bytes > 1e6
        warning('The size of the simulation parameters(%.2f MB) are bigger than 1MB.',Sim_params_metadata.bytes/(1e6));
    end
    metrics.Simulation_params = Sim_params;
end

% Calculate the number of trifurcations. These are the number of
% branches that separate into 3 children.
Branches_NChildren = cellfun(@numel,{Treedata.struct.ChildrenID});
metrics.N_branchtips = nnz(Branches_NChildren == 0);
metrics.N_branchpoints = nnz(Branches_NChildren == 2);
metrics.N_trifurcations = nnz(Branches_NChildren == 3);

% Count the number of tip artifacts
if isfield(Treedata,'TipArtifacts_BranchID')
    metrics.N_tip_artifacts = numel(Treedata.TipArtifacts_BranchID);
end

% Add developmental time, if present.
if isfield(Treedata,'Hours_AEL')
    metrics.Hours_AEL = Treedata.Hours_AEL;
end

% Add the neuron position with respect to the dorsal line.
if isfield(Treedata,'NeuronPosition')
    metrics.NeuronPosition = Treedata.NeuronPosition;
end

% If a timeseries field exist, determine the time indices that will
% be used to perform the timseries analyses.
if isfield(Treedata,'timeseries')
    N_timeseries_samples = metrics.MaxT;
    Simulation_sampling_rate = 1/Timescale;
    Timeseries_sampling_rate = 1/60; %(samples/min) 1 sample per hour.
    Downsampling_factor = Simulation_sampling_rate/Timeseries_sampling_rate;
    Timeseries_sampled_indices = (1:floor(N_timeseries_samples/Downsampling_factor))*Downsampling_factor;
    
    % Add timeseries times when the whole tree structure was saved.
    if isfield(Treedata.timeseries,'Tree')
        metrics.timeseries.Tree_times = Treedata.timeseries.Tree_times;
        metrics.timeseries.Tree_times_units = Treedata.timeseries.Tree_times_units;
    end
else
    Downsampling_factor = [];
    Timeseries_sampled_indices = [];
end

% Calculate the position of the tree nodes.
if isfield(Treedata.struct, 'PointsInd')
    occ = logical(Treedata.occmatrix);
    occsize = size(occ);
    d = 1;
    dy = d/2;
    dx = sqrt(3)/2*d;
    [x, y] = meshgrid(1:occsize(2), 1:occsize(1));
    x = x*dx;
    y = y*dy;
    Tree_nodes_position = [x(occ), y(occ)];
else
    occ = [];
    Tree_nodes_position = {Treedata.struct.PointsPos};
    Tree_nodes_position = cell2mat(cellfun(@(x) x(2:end,:),Tree_nodes_position(:),'uni',0));
    Tree_nodes_position = [Treedata.struct(1).PointsPos(1,:); Tree_nodes_position];
end

% Rescale the tree nodes position and diameter with the lengthscale to create a dimensionful tree.
Tree_struct_dim = scale_tree(Treedata.struct,Lengthscale);

% If an AP axis is defined, rotate the tree to align the Anterior end with the negative x axis.
if isfield(Treedata,'AP_axis_angle') && ~isempty(Treedata.AP_axis_angle)
    Rotation_angle = 180 - Treedata.AP_axis_angle;
    %Tree_struct_dim_rotated = rotate_tree(Tree_struct_dim,Rotation_angle);
    Tree_struct_dim = rotate_tree(Tree_struct_dim,Rotation_angle);
else
    Rotation_angle = [];
end

% Calculate the nodes positions of the dimensionfull tree.
Tree_dim_nodes_position = {Tree_struct_dim.PointsPos};
Tree_dim_nodes_position = cell2mat(cellfun(@(x) x(2:end,:),Tree_dim_nodes_position(:),'uni',0));
Tree_dim_nodes_position = [Tree_struct_dim(1).PointsPos(1,:); Tree_dim_nodes_position];

% Fit a rectangle to the spatial distribution of tree nodes and find the tree aspect ratio.
[rect_size,Rect_boundary_positions] = calculate_tree_size(Tree_nodes_position,'Method','uniform','Rotation',0);
metrics.Rect_boundary_size = rect_size;
metrics.Rect_boundary_positions = Rect_boundary_positions;
metrics.Aspectratio = rect_size(2)/rect_size(1);

% Calculate the internode distance. This is used to scale the length
% measurements.
if ~isfield(Treedata,'Internode_distance')
    Internode_distances = cell2mat(cellfun(@(x) sqrt(sum(diff(x).^2,2)),{Treedata.struct.PointsPos}','Uni',0));
    Internode_distance = uniquetol(Internode_distances);
    assert(numel(Internode_distance) == 1,'The internode distance must be unique');
else
    Internode_distance = Treedata.Internode_distance;
end
metrics.Internode_distance = Internode_distance;

% Calculate the mean length.
meanL = mean([Treedata.struct.Length]) * Internode_distance;
metrics.meanL = meanL;

% Calculate the radius of gyration.
[Rg,COM] = rg(Tree_nodes_position);
metrics.COM = COM;
metrics.Rg_final = Rg;

% Define the polygon that defines the soma region.
if isfield(Treedata,'SomaMaskInd')
    % Define the position of the soma pixels in
    % dimensionful coordinates, centered at the soma.
    if isfield(Treedata,'Image_size')
        Image_size = Treedata.Image_size;
    else
        Image_size = [Treedata.Tiff_metadata.Height, Treedata.Tiff_metadata.Width];
    end
    
    % Change the Soma mask indices to subscripted indices.
    SomaMask_sub_ind = [mod(Treedata.SomaMaskInd - 1,Image_size(1))+1 ceil(Treedata.SomaMaskInd/Image_size(1))];
    
    % Center the indices to the soma and calculate the
    % soma pixel positions in cartesian coordinates.
    SomaMask_sub_ind = bsxfun(@minus,SomaMask_sub_ind,Treedata.SomaImageInd);
    SomaMask_pixel_pos = [SomaMask_sub_ind(:,2), -SomaMask_sub_ind(:,1)]*Lengthscale;
    
    % Rotate the soma pixel positions, if the tree was
    % rotated.
    if ~isempty(Rotation_angle)
        rot_matrix = [cosd(Rotation_angle) -sind(Rotation_angle); sind(Rotation_angle) cosd(Rotation_angle)]; % Rotation matrix for column vectors.
        SomaMask_pixel_pos = SomaMask_pixel_pos*rot_matrix'; % Rotate each row vector using the transpose.
    end
    
    % Calculate the soma polygon vertices.
    Vertices_ind = boundary(SomaMask_pixel_pos(:,1),SomaMask_pixel_pos(:,2));
    Soma_polygon_vertices = {SomaMask_pixel_pos(Vertices_ind,:)};
else
    Soma_polygon_vertices = {};
end

% Calculate various neurons sizes.
metrics.D_95 = calculate_tree_size(Tree_dim_nodes_position,'Method','percentile','Extent_ratio',0.95,'Rotation',0);
metrics.D_99 = calculate_tree_size(Tree_dim_nodes_position,'Method','percentile','Extent_ratio',0.99,'Rotation',0);
[D_uniform, Boundary_uniform] = calculate_tree_size(Tree_dim_nodes_position,'Method','uniform','Rotation',0);
metrics.D_uniform = D_uniform;
metrics.Boundary_uniform = Boundary_uniform;
metrics.Area_uniform = prod(metrics.D_uniform);
%% Calculate miscellaneous metrics.
if options.misc
    % Save the final total length.
    metrics.totL_final = sum([Treedata.struct.Length]*Internode_distance);
    
    % Calculate the mean path length and distance between
    % branchpoints.
    is_branch_terminal = cellfun(@isempty,{Treedata.struct.ChildrenID});
    metrics.meanL_BpToBp = mean([Treedata.struct(~is_branch_terminal).Length])*Internode_distance;
    
    % Calculate the average distance between branchpoints.
    for j=1:numel(Treedata.struct)
        Treedata.struct(j).Branch_dist = sqrt(sum(diff(Treedata.struct(j).PointsPos([1 end],:)).^2,2));
    end
    metrics.meandist_BpToBp = mean([Treedata.struct(~is_branch_terminal).Branch_dist]);
    
    if isfield(Treedata,'timeseries')
        % Define the window size for moving averages.
        Window_size = round(1/(Timescale*Timeseries_sampling_rate));
        
        % Save the a downsampled version of relevant
        % timeseries.
        metrics.timeseries.timescale = Downsampling_factor*Timescale;
        metrics.timeseries.timescale_unit = Treedata.timescale_unit;
        
        % Calculate the rectangle boundary size.
        if isfield(Treedata.timeseries,'Tree')
            for i=1:numel(Treedata.timeseries.Tree)
                Tree_nodes_position_temp = {Treedata.timeseries.Tree{i}.PointsPos};
                Tree_nodes_position_temp = unique(cell2mat(Tree_nodes_position_temp(:)), 'rows');
                
                metrics.timeseries.Rect_boundary_size(i,:) = calculate_tree_size(Tree_nodes_position_temp,'Extent_ratio',0.99);
                
                % Calculate the uniform size.
                Tree_dim_temp = scale_tree(Treedata.timeseries.Tree{i},Lengthscale);
                metrics.timeseries.D_uniform(i,:) = calculate_tree_size(Tree_dim_temp,'Method','uniform','Rotation',0);
            end
        end
        
        % Calculate the number of branches, branchpoints and branchtips.
        NBranches_timeseries = cellfun(@numel,Treedata.timeseries.Lengths);
        NBranchpoints_timeseries = double([Treedata.timeseries.NBranchpoints]);
        NBranchtips_timeseries = NBranches_timeseries - NBranchpoints_timeseries;
        
        metrics.timeseries.NBranches = NBranches_timeseries(Timeseries_sampled_indices);
        metrics.timeseries.NBranchpoints = NBranchpoints_timeseries(Timeseries_sampled_indices);
        metrics.timeseries.NBranchtips = NBranchtips_timeseries(Timeseries_sampled_indices);
        
        % Calculate the mean and total length timeseries.
        metrics.timeseries.totL = cellfun(@sum,Treedata.timeseries.Lengths(Timeseries_sampled_indices)) * Internode_distance;
        metrics.timeseries.meanL = cellfun(@mean,Treedata.timeseries.Lengths(Timeseries_sampled_indices)) * Internode_distance;
        
        % Calculate Rg.
        metrics.Rg = [Treedata.timeseries.Rg];
        
        % Fit the last 25 hours to find Vg, the velocity of Rg.
        N_pointstofit = 25*60/Timescale;
        metrics.Vg = mean(diff(metrics.Rg(end-N_pointstofit+1:end)));
        
        % Find the effective branching frequency.
        Time = 1:N_timeseries_samples;
        
        % Make sure xdata and ydata are column vectors.
        Time = Time(:);
        NBranchpoints_timeseries = NBranchpoints_timeseries(:);
        
        % Fit last 20% of the time points.
        N_fittedpoints = ceil(0.2*N_timeseries_samples);
        Time_fitted = Time(end-N_fittedpoints+1:end);
        NBranchpoints_fitted = NBranchpoints_timeseries(end-N_fittedpoints+1:end);
        
        % Fit the number of branch points to a line. The effective
        % branching frequency (omega_eff) will be given by the slope of
        % that line.
        Fit_params_value = polyfit(Time_fitted, NBranchpoints_fitted, 1);
        metrics.omega_eff = Fit_params_value(1)/Timescale;
        metrics.omega_eff_unit = [Timescale_unit, '^-1'];
        
        % Fit the number of branch points from 30 AEL as was done on
        % experimental images. Morphogenesis starts at 14 AEL.
        N_fittedpoints = N_timeseries_samples-(30-14)*60/Timescale;
        Time_fitted = Time(end-N_fittedpoints+1:end);
        NBranchpoints_fitted = NBranchpoints_timeseries(end-N_fittedpoints+1:end);
        
        % Fit the number of branch points to a line. The effective
        % branching frequency (omega_eff) will be given by the slope of
        % that line.
        Fit_params_value = polyfit(Time_fitted, NBranchpoints_fitted, 1);
        metrics.omega_eff_30AEL = Fit_params_value(1);
        
        % Calculate the phi parameter.
        metrics.phi = metrics.meanL*metrics.omega_eff/metrics.Vg;
        
        % Calculate the branching rate.
        metrics.timeseries.BranchingRate_abs = movmean(double(Treedata.timeseries.NNewBranchpoints),Window_size)/Timescale;
        metrics.timeseries.BranchingRate_abs = metrics.timeseries.BranchingRate_abs(Timeseries_sampled_indices);
        metrics.timeseries.BranchingRate_abs_units = 'min^{-1}';
        metrics.timeseries.Branching_rate_length_normalized = metrics.timeseries.BranchingRate_abs./double(metrics.timeseries.totL*Lengthscale);
        metrics.timeseries.Branching_rate_length_normalized_units = 'min^{-1} \mum^{-1}';
        
        % Calculate the average collision rate.
        metrics.timeseries.CollisionRate_abs = movmean(double(Treedata.timeseries.NCollisions),Window_size)/Timescale;
        metrics.timeseries.CollisionRate_abs = metrics.timeseries.CollisionRate_abs(Timeseries_sampled_indices);
        metrics.timeseries.CollisionRate_abs_units = 'min^{-1}';
        metrics.timeseries.CollisionRate_tip = metrics.timeseries.CollisionRate_abs./double(metrics.timeseries.NBranchtips);
        metrics.timeseries.CollisionRate_tip_units = 'min^{-1}';
        
        % Calculate the average death rate.
        if isfield(Treedata.timeseries,'NDeaths')
            metrics.timeseries.DeathRate_abs = movmean(double(Treedata.timeseries.NDeaths),Window_size)/Timescale;
        else
            metrics.timeseries.DeathRate_abs = movmean(double(Treedata.timeseries.NFullRetractions),Window_size)/Timescale;
        end
        metrics.timeseries.DeathRate_abs = metrics.timeseries.DeathRate_abs(Timeseries_sampled_indices);
        metrics.timeseries.DeathRate_abs_units = 'min^{-1}';
        metrics.timeseries.DeathRate_tip = metrics.timeseries.DeathRate_abs./double(metrics.timeseries.NBranchtips);
        metrics.timeseries.DeathRate_tip_units = 'min^{-1}';
        
        % Calculate the net branching rate.
        metrics.timeseries.BranchingRate_net = movmean(double(Treedata.timeseries.NNewBranchpoints) - double(Treedata.timeseries.NFullRetractions),Window_size)/Timescale;
        metrics.timeseries.BranchingRate_net = metrics.timeseries.BranchingRate_net(Timeseries_sampled_indices);
        metrics.timeseries.BranchingRate_net_abs_units = 'min^{-1}';
    end
    
    % Calculate the convex hull properties.
    [ch_area, ch_diam] = convexhull_analysis(Treedata.struct, 'Lengthscale', Lengthscale);
    metrics.ConvexHull_area = ch_area;
    metrics.ConvexHull_diameter = ch_diam;
end
%% Calculate the fractal dimension.
if options.df
    % Define the radius over which the 2D density correlation is evaluated.
    rad = logspace(0, log10(2*Rg), 50)';
    
    % Calculate the fractal dimension by restricting the
    % boundaries of the tree to the uniform boundary in each dimension.
    % The uniform boundary is defined by the expected range of the branch
    % nodes position assuming that they are uniformly distributed.
    %[C, df_95, df_95_err] = corrdim2(Tree_nodes_position,rad,'R_min',meanL/2,'R_max',Rg,'Method','periodic_boundaries','Boundaries_extent',0.95);
    [C, df, df_err] = corrdim2(Tree_nodes_position,rad,'R_min',meanL/2,'R_max',Rg,'Method','periodic_boundaries','Boundary','uni');
    metrics.rad = rad;
    metrics.C = C;
    metrics.df = df;
    metrics.df_err = df_err;
    
    % Record or calculate the fractal dimension timeseries.
    if isfield(Treedata,'timeseries') && isfield(Treedata.timeseries,'Tree')
        Timeseries_trees = Treedata.timeseries.Tree;
        N_trees = numel(Timeseries_trees);
        metrics.timeseries.df_times = Treedata.timeseries.Tree_times;
        metrics.timeseries.df_times_units = Treedata.timeseries.Tree_times_units;
        metrics.timeseries.df = zeros(N_trees,1);
        for i=1:N_trees
            metrics_temp = analyze_tree(struct('struct',Timeseries_trees{i}),'df',1);
            metrics.timeseries.df(i) = metrics_temp.df;
        end
    end
end
%% Calculate the box counting dimension (type of fractal dimension).
if options.df_box
    R_boxes = logspace(0,log10(2*Rg),50)';
    [df, df_err, N_boxes] = boxdim(Tree_nodes_position,R_boxes,'R_min',meanL,'R_max',Rg,'Translations',true);
    
    metrics.R_boxes = R_boxes;
    metrics.N_boxes = N_boxes;
    metrics.df_box = df;
    metrics.df_box_err = df_err;
end
%% Calculate the hitting probability
if options.hitprob
    % Calculate the hitting probability on the rotated Tree nodes.
    [h, rh, Rh] = hitprob(Tree_dim_nodes_position,50,'Method',3,'Mask','rg','Rotate',0);
    
    % Remove the scale in rh and Rh since these quantities were usually
    % stored dimensionless.
    rh = rh/Lengthscale;
    Rh = Rh/Lengthscale;
    
    metrics.rh = rh;
    metrics.h = h;
    metrics.Rh = Rh;
    
    % Record the meshsize timeseries.
    if isfield(Treedata,'timeseries') && isfield(Treedata.timeseries,'Tree')
        Timeseries_trees = Treedata.timeseries.Tree;
        N_trees = numel(Timeseries_trees);
        metrics.timeseries.Rh_times = Treedata.timeseries.Tree_times;
        metrics.timeseries.Rh_times_units = Treedata.timeseries.Tree_times_units;
        metrics.timeseries.Rh = zeros(N_trees,1);
        for i=1:N_trees
            metrics_temp = analyze_tree(struct('struct',Timeseries_trees{i}),'hitprob',1);
            metrics.timeseries.Rh(i) = metrics_temp.Rh;
        end
    end
end
%% Calculate the branch length distribution.
if options.length
    % Calculate the length counts distribution.
    Lengths = double([Treedata.struct.Length] * Internode_distance * Lengthscale);
    Binwidth = 1; % um
    Binedges = (0:ceil(max(Lengths)/Binwidth))*Binwidth;
    Bincenters = Binedges(1:end-1) + Binwidth/2;
    Bincounts = histc(Lengths, Binedges);
    
    metrics.Length_bincenters = Bincenters;
    metrics.Length_counts = reshape(Bincounts(1:end-1),1,[]);
    
    % Calculate the autocorrelation function of the branch length across depth.
    metrics.Lautocorr = lengthcorr(Treedata.struct);
end
%% Calculate the persistence length averaged over all tree branches.
if options.persistencelength
    Fit_Method = 'tangent_vec';
    
    switch Fit_Method
        case 'tangent_vec'
            % Calculate the persistence length by fitting the average of the
            % tangent vectors inner product
            [Persistence_length, Persistence_length_err] = calculate_PersistenceLength(Treedata.struct,'Method',Fit_Method,'Lengthscale',Lengthscale,'Max_fit_length','mean_length');
            metrics.PersistenceLength = Persistence_length;
            metrics.PersistenceLength_err = Persistence_length_err;
            metrics.PersistenceLength_units = Treedata.lengthscale_unit;
            
        case 'squared_dist'
            % Calculate the persistence length by fitting the mean squared
            % end-to-end distance.
            [Persistence_length, Persistence_length_err] = calculate_PersistenceLength(Treedata.struct,'Method',Fit_Method,'Lengthscale',Lengthscale);
            metrics.PersistenceLength = Persistence_length;
            metrics.PersistenceLength_err = Persistence_length_err;
            metrics.PersistenceLength_units = Treedata.lengthscale_unit;
            
        case 'curvature'
            % Calculate the persistence length by fitting the inverse decimated
            % curvature.
            [Persistence_length, Persistence_length_err] = calculate_PersistenceLength(Treedata.struct,'Method',Fit_Method,'Lengthscale',Lengthscale);
            metrics.PersistenceLength = Persistence_length;
            metrics.PersistenceLength_err = Persistence_length_err;
            metrics.PersistenceLength_units = Treedata.lengthscale_unit;
            
    end
    metrics.PersistenceLength_method = Fit_Method;
    
    % Calculate the persistence for each tree saved in the timeseries.
    if isfield(Treedata,'timeseries') && isfield(Treedata.timeseries,'Tree')
        metrics.timeseries.perL_times = Treedata.timeseries.Tree_times;
        metrics.timeseries.perL_times_units = Treedata.timeseries.Tree_times_units;
        N_trees = numel(Treedata.timeseries.Tree);
        perL = zeros(N_trees,1);
        for i=1:numel(Treedata.timeseries.Tree)
            %perL(i) = calculate_PersistenceLength(Treedata.timeseries.Tree{i},'Method','squared_dist','Lengthscale', Lengthscale);
            perL(i) = calculate_PersistenceLength(Treedata.timeseries.Tree{i},'Method',Fit_Method,'Lengthscale',Lengthscale,'Max_fit_length','mean_length');
        end
        metrics.timeseries.perL = perL;
        metrics.timeseries.perL_units = Treedata.lengthscale_unit;
    end
end