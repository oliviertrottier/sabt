function [Params,Optional_params] = define_params(Params_name)
%% Define common parameters
Base_directory = fileparts(fileparts(mfilename('fullpath')));
Base_params_filename = fullfile(Base_directory,'data','base_parameters.mat');
Params = load(Base_params_filename);
%% Define optional parameters
Optional_params.Branching_rule = 'intensive'; % Constant branching rate at all times.
Optional_params.MaxT = 20; % Simulate for 20 hours
Optional_params.Rng = 10; % Fix RNG for reproducibility
%% Non-default parameters
if nargin > 0
    switch Params_name
        case 'Constant_branching'
            % Add a square boundary that prevents growth. The size of the
            % boundary increases over time.
            Optional_params.Boundary_shape = 'square';
            Optional_params.Boundary_type = 0;
            Optional_params.Boundary_time = [0;Optional_params.MaxT]*60;
            Optional_params.Boundary_size = [150;250];
            
        case 'Extensive_branching_with_decay'
            % Make the total branching rate extensive by scaling with the total
            % branch length.
            Params.omega = 0.0375;
            Optional_params.Branching_rule = 'extensive';
            Optional_params.Branching_decay_type = 'exp+b_capped';
            Optional_params.Branching_decay_params = [0.0025, 0.00125, 0.01];
        
        case 'Extensive_branching_with_decay_and_boundary'
            % Make the total branching rate extensive by scaling with the total
            % branch length. Add square boundary.
            Optional_params.Boundary_shape = 'square';
            Optional_params.Boundary_type = 0;
            Optional_params.Boundary_time = [0;Optional_params.MaxT]*60;
            Optional_params.Boundary_size = [150;250];
            
            Params.omega = 0.0375;
            Optional_params.Branching_rule = 'extensive';
            Optional_params.Branching_decay_type = 'exp+b_capped';
            Optional_params.Branching_decay_params = [0.0025, 0.00125, 0.01];
            
        case 'Drift_diff_extensive_branching_with_decay_and_boundary'
            % Make the total branching rate extensive by scaling with the total
            % branch length. Add square boundary.
            Params.TipDynamicsModel = 'drift_diff';
            Params = rmfield(Params,{'N_states','state_labels','transition_rate_matrix','V_dist_fit'});
            Params.Drift_velocity = 0.05; % um/min
            Params.Diffusion_coefficient = 0.2; % um^2/min
            
            % Initial conditions
            Optional_params.Initial_growth_duration = 1; % min
            Optional_params.Initial_growth_velocity = 1; % um/min
            
            Optional_params.Boundary_shape = 'square';
            Optional_params.Boundary_type = 0;
            Optional_params.Boundary_time = [0;Optional_params.MaxT]*60;
            Optional_params.Boundary_size = [150;250];
            
            Params.omega = 0.0375;
            Optional_params.Branching_rule = 'extensive';
            Optional_params.Branching_decay_type = 'exp+b_capped';
            Optional_params.Branching_decay_params = [0.0025, 0.00125, 0.01];
    end
end
%% Convert optional parameters to a cell array.
Optional_params = reshape([fieldnames(Optional_params)';struct2cell(Optional_params)'],1,[]);
end