function [Treeout, Tips_dynamics, Collisions, Sim_params, Timeseries, Movieframes] = sabt(Tips_params_in, varargin)
% Function that simulates the growth of self-avoiding branching (SABT) tree.

% Branches grow off-lattice at an angle determined by the persistence
% length. The initial branching angle is uniformly distributed between theta_min
% pi-theta_min and pi + theta_min and 2*pi - theta_min where the angle is
% defined with respect to the axis determined by the parent branch.

% The persistence of each branch is coded through the angles that each step
% takes. At each step, the change in the angle is determined by sampling a
% gaussian distribution with mean=0 and std=sqrt(2*delta/beta).
% See details in "The persistence length of double stranded DNA determined
% using dark field tethered particle motion", Brinkers 2009.

% Branches retract and retraction continues through the parent branch if the
% sister branch has been previously deleted.

% Different models for adding new branches to the tree are implemented.
% The different branching models are controlled with the 'Branching_rule'
% option. See the 'Branching_rules' section for more details.

% Branches can regrow when they fully retract to their branchpoint. The
% regrowth probability is given in tips_params_input.

% Stop on error
dbstop if error;

% Clear global variables.
clear global
%% Input summary
% Tips_params_in = structure containing various branch tips parameters
% Fields description:
%
% TipDynamicsModel: 
% Dynamics model used to define the tip velocity. Two models are supported:
% '3state': define the velocity by sampling an input velocity distribution.
%           The distribution changes depending on the state of the tip when
%           sampling is requested.
% 'drift_diff': velocity is determined by a drift-diffusion model with
%               input drift velocity and diffusion constant.
%
% beta:
% Persistence length (um) of the branches.
%
% omega:
% Initial value of the branching rate (min^-1 or min^-1 um^-1).
% See 'Branching_rule' option to see how omega influences the branching rate. 
%
% rebranching_prob:
% Probability that a branch regrows after it shrinks to a zero-length.
% 
% rebranching_delta_theta_dist (optional):
% Sample of rebranching angles (with respect to parent branch) to use when
% rebranching occurs. The set of angles is sampled uniformly when drawing a
% branching angle.
%
% rebranching_delta_theta_dist_fit (optional):
% Fit of the rebranching angles sample. When a fit is defined, drawing a
% rebranching angle is done with this fit assuming it is a 
% prob.ProbabilityDistribution object.
%
%
% 
% Parameters for TipDynamicsModel = '3state'
%
% N_states:
% Number of states in the tip dynamical system. Only N_states=3 is
% currently implemented.
%
% transition_rate_matrix :
% transition rate matrix determining the rate constants (min^-1) of each 
% transition. The (i,j)^th entry corresponds to the rate of transitioning
% from state i to state j.
%
% transition_rate_matrix_time (optional) :
% times (min) defining the schedule for changing the transition rate matrix.
% This field is necessary when transition_rate_matrix is a 3D array.
%
% V_dist_fit: 
% cit object definining the 3-component velocity distribution. The object
% must define the parameters mum, sipm, mu0, sig0, mup, sigp.
% The shrinking state velocity distribution is log-normal with parameters mum, sigm
% The paused state velocity distribution is normal with parameter mu0, sig0
% The growing state velocity distribution is log-normal with parameters mup, sigp
%
% V_dist_time (optional):
% times (min) defining the schedule for changing the velocity distribution.
% This field is necessary when V_dist_fit is a cell array of cfit objects.
%
%
%
% Parameters for TipDynamicsModel = 'drift_diff':
% 
% Drift_velocity:
% Drift velocity (mum/min) of the tips.
%
% Drift_velocity_time (optional):
% times (min) defining the schedule for changing the drift velocity.
%
% Diffusion_coefficient:
% diffusion coefficient (mum^2/min) of the tips.
%
% Diffusion_coefficient_time (optional):
% times (min) defining the schedule for changing the diffusion coefficient.
%% Parse parameters
% In older versions, the 2nd and 3rd inputs correspond to the value of beta
% and omega, respectively. Parse the input values to ensure backward
% compatibility.
if all(cellfun(@isnumeric,varargin(1:2)))
    [~,funcname] = fileparts(mfilename);
    warning('Using old calling method:\n%s(Tips_params,beta,omega,...)\nbeta and omega are now defined inside Tips_params.',funcname);
    Tips_params_in.beta = varargin{1};
    Tips_params_in.omega = varargin{2};
    varargin(1:2) = [];
end
%% Parse optional parameters
% Options are defined by calling the function as follows:
% sabt(Tips_params_in,...,'Option_name',Option_value,...)
% where Option_name is the name of one of the options defined below and
% Option_value is its value.
p = inputParser;

% Simulation scales.
% The length scale (um) is used to define the distance between two branch
% nodes when sampling the 2D path of the branches.
% The time scale (min) is used to define the duration of one time step.
addParameter(p, 'lengthscale', 0.1); 
addParameter(p, 'timescale', 0.1);

% Contact distance
% When checking branch-branch collisions, a branch is in contact with
% another branch if any of its branch nodes are within a distance
% Contact_distance from the nodes of other neighboring branches.
addParameter(p, 'Contact_distance', 0.4); % um

% Option to add a velocity or duration bias (um/min) on the shrinking stage
% after collision.
addParameter(p, 'Collision_v_bias', 0);
addParameter(p, 'Collision_t_bias', 0);

% Change the state that the tip transits to upon collision with other dendrites.
addParameter(p, 'Collision_state', 'shrinking',@(x) ismember(x,{'shrinking','paused'}));

% Define the time that the tip stays idle after collision before switching to the collision
% state.
addParameter(p, 'Collision_duration', 0, @(x) isa(x,'prob.ToolboxFittableParametricDistribution') || x >= 0);

% Define the initial growth velocity and duration when a new branch tip is
% born.
addParameter(p, 'Initial_conditions_time',0,@(x) isnumeric(x));
addParameter(p, 'Initial_growth_velocity','Random',@(x) strcmp(x,'Random') || isnumeric(x));
addParameter(p, 'Initial_growth_duration','Random',@(x) strcmp(x,'Random') || isnumeric(x));

% Option for the branching mechanism.
addParameter(p, 'Branching_rule', 'uniform', @(x) ismember(x,{'uniform','age_dependent','exp_decay_boundary','extensive','extensive_branch','intensive'}));

% Define the branching angle distribution in radians. Use a uniform distribution by
% default (see below).
addParameter(p, 'Branching_angle_dist', []);

% Option to add a velocity bias (um/min) on the shrinking stage after collision.
addParameter(p, 'Branching_decay_params', []);
addParameter(p, 'Branching_decay_type', 'none');

% Option to add a physical boundary around the tree.
addParameter(p, 'Boundary_size', inf); % Size that define the shape of the boundary.
addParameter(p, 'Boundary_time', 0); % Time associated with the boundary size array.
addParameter(p, 'Boundary_type', -1, @(x) ismember(x,[0,-1,1])); % -1=repulsive, 0=pause, 1=periodic.
addParameter(p, 'Boundary_shape', 'square', @(x) ismember(x,{'circle','square','rectangle'})); % 'square' or 'circle'

% Add an exclusion zone representing the Soma. The zone is
% defined by a circle around the soma with a given radius (um).
addParameter(p, 'Soma_radius', 0);

% Position of the soma on the 2D plane.
addParameter(p, 'Soma_position', [0,0], @(x) numel(x) == 2);

% Initial number of branches at the soma.
addParameter(p, 'N_soma_branches', 3);

% Change the maximal simulation time. The default total growth time is 100 hrs.
addParameter(p, 'MaxT', 100);

% Use different tip dynamics parameters after collision.
% Definition of the fields follows the same rules as the input tips'
% parameters structure.
addParameter(p, 'TipsParamsPostCollision', []);

% Duration of the post collision dynamics (must be non zero when TipsParamsPostCollision is not empty).
addParameter(p, 'alpha', 0);

% Record movie. The recorded frames are changed with Movie_recorded_frames.
addParameter(p, 'Movie', false);
addParameter(p, 'Movie_recorded_frames', 'All', @(x) strcmp(x,'All') || (isnumeric(x) && all(x==round(x))));

% Record timeseries.
addParameter(p, 'Timeseries', false); % Activate timeseries recording.
addParameter(p, 'Timeseries_metric', {},@(x)iscell(x) && all(ismember(x,{'df','meshsize'}))); % Cell array of strings determined the types of metrics to calculate.
addParameter(p, 'Timeseries_metric_times', []); % Times (minutes) when the metrics are calculated.
addParameter(p, 'Timeseries_tree_times', []); % Save the structure of the tree at the given input times (minutes).
addParameter(p, 'Timeseries_Rg', true); % Calculate the radius of gyration at each simulation time step.
addParameter(p, 'Timeseries_tracks', []); % Save tip growth tracks in given time intervals during the simulation.

% Set the initial growth velocity (um/min) of the soma branches.
addParameter(p, 'SomaInitialVelocity', 1, @(x) x > 0);

% Record the tip dynamics (State and number of steps per timeframe).
addParameter(p, 'RecordTipDynamics', false);

% Record information about collisions.
addParameter(p, 'RecordCollisions', false);

% Calculate convex hull timeseries (considerably lengthen the simulation time).
addParameter(p, 'ConvexHull', false);

% Set the RNG at the beginning of the simulation.
addParameter(p, 'Rng', 0);

% Show waitbar.
addParameter(p, 'waitbar', false);

% Display warnings.
addParameter(p, 'warnings', false);

% Profile the code.
addParameter(p, 'Profile', false);

parse(p, varargin{:});
global options
options = p.Results;
options_isdefault = cell2struct(num2cell(ismember(p.Parameters, p.UsingDefaults)), p.Parameters, 2);
%% Miscellaneous parameters
% Turn on/off debugging.
global Debugmode
Debugmode = false;
%Debugmode = true;

% Initialize dummy variables for debugging the static workspace.
global a s

% Create the folders where the outputs will be saved.
savingfolder = cell2mat(regexp(mfilename('fullpath'), '/.*/', 'match'));
if ~isdir(savingfolder)
    mkdir(savingfolder);
end

% Unpack some parameters to get concise code.
beta = Tips_params_in.beta;
omega = Tips_params_in.omega;

% Initialize structure to save the simulation parameters
Sim_params = struct();
Sim_params.beta = beta;
Sim_params.beta_unit = 'um';
Sim_params.omega = omega;
Sim_params.omega_unit = ''; % Defined in the branching rule section.

% Set and save the state of the rng.
if ~options_isdefault.Rng
    rng(options.Rng);
end
Sim_params.Rng = rng;

% Upack the length and time scales.
lengthscale = options.lengthscale; % um
timescale = options.timescale; % min

Sim_params.lengthscale = lengthscale;
Sim_params.lengthscale_unit = 'um';
Sim_params.timescale = timescale;
Sim_params.timescale_unit = 'min';

% Non-dimensionalize the persistence length.
beta = beta/lengthscale;

% Start the growth at ~25 min (which corresponds to the start time in sabt14 (5th timestep * 5 min/timestep)).
Branching_start_time = round(25/timescale);
Sim_params.Branching_start_time = Branching_start_time;

% Define the initial number of branches at the soma.
N_soma_branches = options.N_soma_branches;
Sim_params.N_soma_branches = N_soma_branches;
Soma_pos = reshape(options.Soma_position,[1,2])/lengthscale;
Sim_params.Soma_position = options.Soma_position;
Soma_radius = options.Soma_radius/lengthscale;
Soma_radius_squared = Soma_radius^2;
Sim_params.Soma_radius = options.Soma_radius;

% Define the total number of temporal steps.
MaxT = round(options.MaxT/(timescale/60)); % options.MaxT is given in hours.
Sim_params.MaxT = MaxT;

% Define the size of each growth or shrinkage step as equal to the lengthscale.
delta = 1;
Sim_params.delta = delta;

% Set the contact distance.
global Contact_dist
Contact_dist = options.Contact_distance/lengthscale; % um;
Sim_params.Contact_dist = Contact_dist;

% Define the maximal number of nodes that each branch can have. This is
% used to define the size of the branch arrays. Simulation errors out if
% the number of nodes of a given branch exceeds this number.
N_nodes_per_branch = 10000;

% Define the maximal number of branches.
N_branches_max = 40000;

% Define the maximal number of tips
% Use a higher maximal number of tips when running on a remote server with
% more memory.
if RAM_total() > 32
    N_tips_max = 600000;
else
    N_tips_max = 200000;
end
%% Tips parameters
% Non-dimensionalize the tips parameters.
global Tips_params
Tips_params = struct();

% Define the tip dynamics model.
if isfield(Tips_params_in,'TipDynamicsModel') && ~isempty(Tips_params_in.TipDynamicsModel)
    Tips_params.TipDynamicsModel = Tips_params_in.TipDynamicsModel;
    assert(ismember(Tips_params.TipDynamicsModel,{'3state','drift_diff'}),'The tip dynamics model is undefined.');
else
    % The default is the 3-state model.
    Tips_params.TipDynamicsModel = '3state';
end
if strcmp(Tips_params.TipDynamicsModel,'3state')
    Tips_params.N_states = Tips_params_in.N_states;
    assert(Tips_params.N_states==3,'The transition dynamics is implemented for 3 states only.');
else
    Tips_params.N_states = 3; % Use 3 for compatibility with other functions.
end

% '3state' parameters.
% Define a map between the string labels and the integer labels of the
% states.
State_char2int = containers.Map('KeyType','char','ValueType','double');
State_char2int('shrinking') = 1; State_char2int('S') = 1;
State_char2int('paused') = 2; State_char2int('P') = 2;
State_char2int('growing') = 3; State_char2int('G') = 3;

global Tips_params_post_coll
Tips_params_post_coll = [];
alpha = [];
if strcmp(Tips_params.TipDynamicsModel,'3state')
    % Define the transition probability matrix. The first two dimensions of the
    % input transition rate matrix determine the transition rates.
    % Specifically, T(i,j) corresponds to the transition rate from state i to
    % state j where i,j = 1,...,3 and 1=S, 2=P, 3=G is the mapping of the
    % states.
    
    % If the input transition rate matrix array is 3-dimensional, the 3rd dimension
    % corresponds to time. In this case, the index of the 3rd dimension refers
    % to the specific time indicated by transition_rate_matrix_time.
    % Specifically, T(i,j,k) corresponds to the transition rate from the i
    % state to the j state at time transition_rate_matrix_time(k). The time is
    % asumed to be given in minutes (dimensionfull).
    is_transition_rate_matrix_constant = ndims(Tips_params_in.transition_rate_matrix) == 2;
    if is_transition_rate_matrix_constant
        Tips_params.transition_prob_matrix = Tips_params_in.transition_rate_matrix*timescale;
    else
        assert(ndims(Tips_params_in.transition_rate_matrix) == 3,'Transition matrix is not 3-dimensional.');
        assert(isfield(Tips_params_in,'transition_rate_matrix_time'),'The field transition_rate_matrix_time is needed to perform the transition rate updates.');
        assert(size(Tips_params_in.transition_rate_matrix,3) == numel(Tips_params_in.transition_rate_matrix_time),'The number of transition matrix update times must match the number of transition matrices.');
        
        % Determine the times when the transition matrix will be updated.
        Transition_matrix_update_times = [round(Tips_params_in.transition_rate_matrix_time/timescale), inf];
        Transition_matrix_update_ind = 1;
        Transition_matrix_update_time = Transition_matrix_update_times(Transition_matrix_update_ind + 1);
        
        % Initialize the transition matrix.
        Tips_params.transition_prob_matrix = Tips_params_in.transition_rate_matrix(:,:,Transition_matrix_update_ind)*timescale;
    end
    
    % Define useful routines to sample transitions.
    Tips_params.Calculate_exit_rates = @(x) sum(x.transition_prob_matrix,2);
    Tips_params.Exit_rates = Tips_params.Calculate_exit_rates(Tips_params);
    Tips_params.Transition_states = cell2mat(arrayfun(@(x)[1:x-1, x+1:Tips_params.N_states] ,(1:Tips_params.N_states)','Uni',0));
    Tips_params.Calculate_transition_dest_prob = @(x) x.transition_prob_matrix((x.Transition_states-1)*x.N_states + (1:x.N_states)')./x.Exit_rates;
    Tips_params.Transition_dest_prob = Tips_params.Calculate_transition_dest_prob(Tips_params);
    
    N_states = Tips_params.N_states;
    States_datatype = min_int_type(N_states,false);
    
    % Rescale parameters of the velocity distribution fit.
    V_dist_fits = Tips_params_in.V_dist_fit;
    N_dists = numel(V_dist_fits);
    if N_dists == 1 && ~iscell(V_dist_fits)
        V_dist_fits = {V_dist_fits};
    end
    
    % If there are more than one distribution, assume the fit objects are
    % concatenated in a cell array.
    % The sigma parameter of a log-normal distribution is unaffected by
    % rescaling.
    warning('off','curvefit:cfit:subsasgn:coeffsClearingConfBounds')
    for i=1:N_dists
        V_dist_fits{i}.mum = V_dist_fits{i}.mum - log(lengthscale/timescale);
        V_dist_fits{i}.mup = V_dist_fits{i}.mup - log(lengthscale/timescale);
        V_dist_fits{i}.mu0 = V_dist_fits{i}.mu0/lengthscale*timescale;
        V_dist_fits{i}.sig0 = V_dist_fits{i}.sig0/lengthscale*timescale;
    end
    warning('on','curvefit:cfit:subsasgn:coeffsClearingConfBounds')
    Tips_params.V_dist_fit = V_dist_fits{1};
    
    % Determine if the velocity distribution is constant in time.
    is_V_dist_constant = ~isfield(Tips_params_in,'V_dist_time');
    if ~is_V_dist_constant
        assert(numel(Tips_params_in.V_dist_time) == numel(Tips_params_in.V_dist_fit),'The number of velocity distributions must match the number of velocity distribution update times.');
        
        V_dist_update_times = [round(Tips_params_in.V_dist_time/timescale) Inf];
        V_dist_update_ind = 2;
        V_dist_update_time = V_dist_update_times(V_dist_update_ind);
    end
    
    % Do the same non-dimensionalization for the post collision tips
    % parameters.
    post_coll_dynamics = ~isempty(options.TipsParamsPostCollision);
    if post_coll_dynamics
        % Define the duration of the post collision dynamics.
        if options.alpha > 0
            % Assume the input alpha is dimensionful.
            alpha = options.alpha/timescale;
        else
            error('alpha must be non-zero when the post collision dynamics is different');
        end
        Sim_params.alpha = alpha;
        
        Tips_params_post_coll_input = options.TipsParamsPostCollision;
        Tips_params_post_coll = struct();
        Tips_params_post_coll.N_states = Tips_params_post_coll_input.N_states;
        
        % Determine if the post-collision transition rate matrix is constant.
        Is_transition_rate_matrix_post_coll_constant = ndims(Tips_params_post_coll_input.transition_rate_matrix) == 2;
        if Is_transition_rate_matrix_post_coll_constant
            Tips_params_post_coll.transition_prob_matrix = Tips_params_post_coll_input.transition_rate_matrix*timescale;
        else
            % Verify inputs on the transition matrix.
            assert(ndims(Tips_params_post_coll_input.transition_rate_matrix) == 3,'Transition matrix is not 3-dimensional.');
            assert(isfield(Tips_params_post_coll_input,'transition_rate_matrix_time'),'The field transition_rate_matrix_time is needed to perform the transition rate updates.');
            assert(size(Tips_params_post_coll_input.transition_rate_matrix,3) == numel(Tips_params_post_coll_input.transition_rate_matrix_time),'The number of transition matrix update times must match the number of transition matrices.');
            
            % Determine the times when the transition matrix will be updated.
            Transition_matrix_post_coll_update_times = [round(Tips_params_post_coll_input.transition_rate_matrix_time/timescale) inf];
            Transition_matrix_post_coll_update_ind = 1;
            Transition_matrix_post_coll_update_time = Transition_matrix_post_coll_update_times(Transition_matrix_post_coll_update_ind + 1);
            
            % Initialize the transition matrix.
            Tips_params_post_coll.transition_prob_matrix = Tips_params_post_coll_input.transition_rate_matrix(:,:,Transition_matrix_post_coll_update_ind)*timescale;
        end
        Tips_params_post_coll.Calculate_exit_rates = Tips_params.Calculate_exit_rates;
        Tips_params_post_coll.Exit_rates = Tips_params.Calculate_exit_rates(Tips_params_post_coll);
        Tips_params_post_coll.Transition_states = cell2mat(arrayfun(@(x)[1:x-1, x+1:Tips_params_post_coll.N_states] ,(1:Tips_params_post_coll.N_states)','Uni',0));
        Tips_params_post_coll.Calculate_transition_dest_prob = Tips_params.Calculate_transition_dest_prob;
        Tips_params_post_coll.Transition_dest_prob = Tips_params_post_coll.Calculate_transition_dest_prob(Tips_params_post_coll);
        
        % Rescale the parameters of the velocity distribution fit.
        V_dist_post_coll_fits = Tips_params_post_coll_input.V_dist_fit;
        N_dists = numel(V_dist_post_coll_fits);
        if N_dists == 1 && ~iscell(V_dist_post_coll_fits)
            V_dist_post_coll_fits = {V_dist_post_coll_fits};
        end
        
        warning('off','curvefit:cfit:subsasgn:coeffsClearingConfBounds')
        % If there are more than one distribution, assume the fit objects are
        % concatenated in a cell array.
        for i=1:N_dists
            V_dist_post_coll_fits{i}.mum = V_dist_post_coll_fits{i}.mum - log(lengthscale/timescale);
            V_dist_post_coll_fits{i}.mup = V_dist_post_coll_fits{i}.mup - log(lengthscale/timescale);
            V_dist_post_coll_fits{i}.mu0 = V_dist_post_coll_fits{i}.mu0/lengthscale*timescale;
            V_dist_post_coll_fits{i}.sig0 = V_dist_post_coll_fits{i}.sig0/lengthscale*timescale;
        end
        Tips_params_post_coll.V_dist_fit = V_dist_post_coll_fits{1};
        warning('on','curvefit:cfit:subsasgn:coeffsClearingConfBounds')
        
        % Determine if the velocity distribution post-collision is constant in time.
        constant_V_dist_post_coll = ~isfield(Tips_params_post_coll_input,'V_dist_time');
        if ~constant_V_dist_post_coll
            assert(numel(Tips_params_post_coll_input.V_dist_time) == numel(Tips_params_post_coll_input.V_dist_fit),'The number of velocity distributions must match the number of velocity distribution update times.');
            
            V_dist_post_coll_update_times = [round(Tips_params_post_coll_input.V_dist_time/timescale) Inf];
            V_dist_post_coll_update_ind = 2;
            V_dist_post_coll_update_time = V_dist_post_coll_update_times(V_dist_post_coll_update_ind);
        end
    end
end

% Non-dimensionalize the drift and diffusion coefficients.
if strcmp(Tips_params.TipDynamicsModel,'drift_diff')
    Tips_params.Drift_velocity = Tips_params_in.Drift_velocity(1)/lengthscale*timescale;
    if isfield(Tips_params_in,'Drift_velocity_time')
        V_drift_update_times = [round(Tips_params_in.Drift_velocity_time/timescale), Inf];
    else
        V_drift_update_times = [0,Inf];
    end
    V_drift_update_ind = 2;
    V_drift_update_time = V_drift_update_times(V_drift_update_ind);
    
    Tips_params.Diffusion_coefficient = Tips_params_in.Diffusion_coefficient(1)/lengthscale^2*timescale;
    Tips_params.Diffusion_displacement_sigma = sqrt(2*Tips_params.Diffusion_coefficient);
    if isfield(Tips_params_in,'Diffusion_coefficient_time')
        D_coeff_update_times = [round(Tips_params_in.Diffusion_coefficient_time/timescale), Inf];
    else
        D_coeff_update_times = [0,Inf];
    end
    D_coeff_next_update_ind = 2;
    D_coeff_next_update_time = D_coeff_update_times(D_coeff_next_update_ind);
    
    % Save some fields.
    Sim_params.Drift_velocity = Tips_params_in.Drift_velocity/lengthscale*timescale;
    Sim_params.Drift_velocity_time = V_drift_update_times(1:end-1);
    Sim_params.Diffusion_coefficient = Tips_params_in.Diffusion_coefficient/lengthscale^2*timescale;
    Sim_params.Diffusion_coefficient_time = D_coeff_update_times(1:end-1);
end

% Non-dimensionalize the collision bias parameters.
options.Collision_v_bias = options.Collision_v_bias/lengthscale*timescale;
Sim_params.Collision_v_bias = options.Collision_v_bias;
options.Collision_t_bias = round(options.Collision_t_bias/timescale);
Sim_params.Collision_t_bias = options.Collision_t_bias;
Sim_params.Collision_state = options.Collision_state;

% Redefine the collision state using the integer labels.
options.Collision_state = State_char2int(options.Collision_state);

% Define the collision duration sampler and non-dimensionalize the duration time.
switch class(options.Collision_duration)
    case 'double'
        Collision_duration_sampler = @() options.Collision_duration/timescale;
    case 'prob.ExponentialDistribution'
        Collision_duration_dist = options.Collision_duration;
        Collision_duration_dist.mu = Collision_duration_dist.mu/timescale;
        Collision_duration_sampler = @() random(Collision_duration_dist);
    otherwise
        error('%s is not a supported distribution for the contact duration',class(options.Collision_duration));
end
Sim_params.Collision_duration = options.Collision_duration;

% Determine the rebranching parameters.
if isfield(Tips_params_in, 'rebranching_prob') && Tips_params_in.rebranching_prob > 0
    Rebranching_prob = Tips_params_in.rebranching_prob;
    Rebranching_angle_dist = Tips_params_in.rebranching_delta_theta_dist;
    
    % Check that the angles are within a range of pi. This will create an
    % error if the distribution is in degrees and is wide enough to span
    % further than pi rad.
    assert(all(abs(Rebranching_angle_dist)<=pi),'The rebranching angle distribution is not in the range [-pi,pi].');
    
    % Load the distribution fit, if any.
    if isfield(Tips_params_in,'rebranching_delta_theta_dist_fit')
        Rebranching_angle_dist_fit = Tips_params_in.rebranching_delta_theta_dist_fit;
        
        % Check that the distribution fit is a probability distribution
        % object.
        assert(isa(Rebranching_angle_dist_fit,'prob.ProbabilityDistribution'),'The rebranching angle distribution ist not a prob. dist. object');
    else
        Rebranching_angle_dist_fit = [];
    end
else
    % By default, branches cannot rebranch upon reaching their branch
    % point.
    Rebranching_prob = 0;
    Rebranching_angle_dist = [];
end
% Initialize rebranching parameters used by local functions.
Sim_params.Rebranching_prob = Rebranching_prob;
Rebranching_angle_ecdf = [];
Rebranching_angle_ecdf_x = [];

% Define the initial growth velocity (um/min) and duration (min).
if isnumeric(options.Initial_conditions_time)
    Initial_conditions_time_timeseries = round([options.Initial_conditions_time(:)'/timescale,inf]);
    Initial_conditions_update_ind = 2;
else
    error('The initial conditions time is not a numerical value.');
end

if strcmp(options.Initial_growth_velocity,'Random')
    Initial_growth_velocity_timeseries = {options.Initial_growth_velocity};
elseif iscell(options.Initial_growth_velocity)
    Initial_growth_velocity_timeseries = options.Initial_growth_velocity;
elseif isnumeric(options.Initial_growth_velocity)
    Initial_growth_velocity_timeseries = options.Initial_growth_velocity*timescale/lengthscale;
else
    error('The initial growth velocity is neither "Random" nor a numerical value.');
end
Initial_growth_velocity_timeseries = Initial_growth_velocity_timeseries(:)';

if strcmp(options.Initial_growth_duration,'Random')
    Initial_growth_duration_timeseries = {options.Initial_growth_duration};
elseif iscell(options.Initial_growth_duration)
    Initial_growth_duration_timeseries = options.Initial_growth_duration;
elseif isnumeric(options.Initial_growth_duration)
    Initial_growth_duration_timeseries = options.Initial_growth_duration/timescale;
else
    error('The initial growth duration is neither "Random" nor a numerical value.');
end
Initial_growth_duration_timeseries = Initial_growth_duration_timeseries(:)';

N_initial_conditions_timeseries = [numel(Initial_conditions_time_timeseries)-1,numel(Initial_growth_velocity_timeseries),numel(Initial_growth_duration_timeseries)];
assert(all(N_initial_conditions_timeseries == N_initial_conditions_timeseries(1)),'The initial conditions options are not the same length.');

Sim_params.Initial_conditions_time = Initial_conditions_time_timeseries(1:end-1);
Sim_params.Initial_growth_velocity = Initial_growth_velocity_timeseries;
Sim_params.Initial_growth_duration = Initial_growth_duration_timeseries;
Initial_growth_velocity = Initial_growth_velocity_timeseries(1);
Initial_growth_duration = Initial_growth_duration_timeseries(1);

% Determine how to handle branches that retract back to their branchpoint.
N_retractions_max = Inf;
Rebranching_scenario = 'pause_tip';
switch Rebranching_scenario
    case {'pause_tip','pause_tip_in_dying_state','delete_branches'}
    otherwise
        error('The rebranching scenario (%s) is undefined.',Rebranching_scenario);
end
Sim_params.Rebranching_scenario = Rebranching_scenario;
Sim_params.N_retractions_max = N_retractions_max;
%% Boundary
% Determine if a boundary exists.
if options.Boundary_size(1) < inf
    % Vectorize the input boundary info.
    N_boundary_inputs = numel(options.Boundary_time);
    if ~any(size(options.Boundary_size)==N_boundary_inputs)
        error('The size of the boundary_time and boundary_size arrays do not match.')
    end
    
    % Make sure the temporal axis of the boundary arrays corresponds to the
    % 1st dimension.
    if size(options.Boundary_size,1) ~= N_boundary_inputs
        options.Boundary_size = options.Boundary_size';
    end
    Boundary_times = options.Boundary_time(:)/timescale;
    Boundary_size = options.Boundary_size/lengthscale;
    
    % Make sure that the boundary is defined until the end of the
    % simulation.
    assert(Boundary_times(end)<= MaxT,'The boundary must be defined until the end of the simulation');
    
    % Linearly interpolate the distance to the boundary using the input data.
    Boundary_size = interp1(Boundary_times,Boundary_size,(1:MaxT)');
    
    % Check that the boundary size is positive definite at all times.
    if any(Boundary_size <= 0,[1,2])
        error('The boundary size must be positive definite.');
    end
    
    % Check that the boundary sizes are monotonically increasing.
    if any(diff(Boundary_size) < 0,[1,2])
        error('The boundary sizes are not monotonically increasing.');
    end
    is_boundary_existent = true;
    
    % Define the distance to the boundary from the origin to check boundary collisions.
    Boundary_positions = Boundary_size/2;
    
    % Define the current boundary size at the initial time used in
    % the search_nodes() function.
    Boundary_size_temp = Boundary_size(1,:);
else
    Boundary_size = nan;
    is_boundary_existent = false;
    Sim_params.Boundary_pos = [];
    Boundary_size_temp = [];
end
Sim_params.Boundary_size = Boundary_size;
Sim_params.Boundary_type = options.Boundary_type;
Sim_params.Boundary_shape = options.Boundary_shape;
%% Branching rules.
% Non-dimensionalize the branching rate.
switch options.Branching_rule
    case 'intensive'
        omega = omega*timescale;
        Sim_params.omega_unit = 'min^{-1}';
    case {'uniform','extensive','extensive_branch','age_dependent','exp_decay_boundary','exp_decay_pathlength'}
        omega = omega*timescale*lengthscale;
        Sim_params.omega_unit = 'min^{-1} um^{-1}';
    otherwise
        error('Branching rule %s not supported',options.Branching_decay_type);
end

switch options.Branching_rule
    case 'intensive'
        % Branching is an intensive property of the tree. This implies that
        % it does not change with the size of the tree. In other words, the
        % total branching rate of the tree doesn't increase with the total branch length.
        switch options.Branching_decay_type
            case 'none'
                omega_timeseries = omega*ones(1,MaxT);
            case 'linear_interp'
                omega_timeseries = interp1(options.Branching_decay_params(:,1)/timescale,options.Branching_decay_params(:,2)*timescale,(0:MaxT-1),'linear','extrap');
            otherwise
                error('Decay type %s not supported',options.Branching_decay_type);
        end
    case {'uniform','extensive','extensive_branch'}
        % Branching is an extensive property of the tree. This implies that
        % the branching rate changes as a function of the size (or amount
        % of matter) in the tree. Practically, this implies that the
        % total branching rate increases with the total branch length.
        % 'uniform' was used previously as a synonym of 'extensive'.
        switch options.Branching_decay_type
            case 'none'
                omega_timeseries = omega*ones(1,MaxT);
            case 'linear'
                % The branching rate decays linearly from the initial
                % branching rate.
                omega_timeseries = omega - options.Branching_decay_params*timescale*(0:MaxT-1);
            case 'exp'
                % The branching rate decays exponentially from the initial
                % branching rate.
                omega_timeseries = omega*exp(-options.Branching_decay_params*timescale*(0:MaxT-1));
            case 'exp+b'
                % The branching rate decays exponentially from the initial
                % branching rate with an additive constant.
                omega_timeseries = omega*exp(-options.Branching_decay_params(1)*timescale*(0:MaxT-1)) + options.Branching_decay_params(2)*timescale*lengthscale;
            case 'exp+b_capped'
                % The branching rate decays exponentially from the initial
                % branching rate with an additive constant. The rate is
                % capped to a given maximum.
                omega_timeseries = omega*exp(-options.Branching_decay_params(1)*timescale*(0:MaxT-1)) + options.Branching_decay_params(2)*timescale*lengthscale;
                omega_timeseries = min(omega_timeseries,options.Branching_decay_params(3)*timescale*lengthscale);
            case 'zero_exp'
                % The branching rate does not decay until a certain time
                % and then starts decaying exponentially from the initial branching rate.
                decay_activation_time = round(options.Branching_decay_params(2)/timescale);
                omega_timeseries = omega*[ones(1,decay_activation_time), exp(-options.Branching_decay_params(1)*timescale*(0:MaxT-decay_activation_time-1))];
            case 'exp_stop'
                % The branching rate decays exponentially from the initial
                % branching rate and stops at a given time.
                omega_timeseries = omega*exp(-options.Branching_decay_params(1)*timescale*(0:MaxT-1));
                decay_stop_ind = round(options.Branching_decay_params(2)/timescale);
                omega_timeseries(decay_stop_ind:end) = omega_timeseries(decay_stop_ind);
            otherwise
                error('Decay type %s not supported',options.Branching_decay_type);
        end
    case 'age_dependent'
        % Branching is a function of the age of the branch. Each node of
        % each branch is given an age based on the first time they
        % appeared.
        switch options.Branching_decay_type
            case 'none'
                % The branching rate does not depend on time. Used for
                % debugging.
                is_node_branching_func = @(t) rand(size(t)) < (1 - exp(-omega*ones(size(t))));
            case 'exp'
                % The probability that a given node spawns a new branch
                % is an exponentiall decaying function of its age.
                % Define a function that will determine if a node of age t 
                % will spawn a new branchpoint. t is dimensionless in the
                % function defined below.
                is_node_branching_func = @(t) rand(size(t)) < (1 - exp(-omega*exp(-options.Branching_decay_params*timescale.*t)));
            otherwise
                error('Decay type %s not supported',options.Branching_decay_type);
        end
    case 'exp_decay_boundary'
        % The rate of branching is a function of its position to the
        % boundary. Branches closer to the boundary are more likely to
        % spawn new branches. The rate decays exponentially from the
        % boundary towards the center and it is azimuthally symmetric.
        % lambda defines the lengthscale of the decay.
        switch options.Branching_decay_type
            case 'exp'
                % Non-dimensionalize the decay length.
                lambda = options.Branching_decay_params(1)/lengthscale;
                branching_rate_func = @(r,R_tree) omega*exp((min(r,R_tree) - R_tree)/lambda);
            case 'spatial_decay+temporal_capped_decay'
                % The branching rate decays exponentially from the initial
                % branching rate with an additive constant. The rate is
                % capped to a given maximum.
                omega_timeseries = omega * exp(-options.Branching_decay_params(1)*timescale*(0:MaxT-1)) + options.Branching_decay_params(2)*timescale*lengthscale;
                omega_timeseries = min(omega_timeseries,options.Branching_decay_params(3)*timescale*lengthscale);
                
                % Non-dimensionalize the decay length.
                lambda = options.Branching_decay_params(4)/lengthscale;
                branching_rate_func = @(d,t) omega_timeseries(t)*exp(-d/lambda);
            otherwise
                error('Decay type %s not supported',options.Branching_decay_type);
        end
    case 'exp_decay_pathlength'
        % Branching is a function of the path length to the soma. Nodes
        % further in path length from the soma are less likely to spawn new
        % branches. The input omega corresponds to the branching rate at
        % the furthest point on the tree. lambda defines the lengthscale of
        % the decay.
        switch options.Branching_decay_type
            case 'exp'
                % Non-dimensionalize the decay length.
                lambda = options.Branching_decay_params(1)/lengthscale;
                branching_rate_func = @(L) omega*exp(-(L - max(L))/lambda);
            otherwise
                error('Decay type %s not supported',options.Branching_decay_type);
        end
end
Sim_params.Branching_rule = options.Branching_rule;
Sim_params.Branching_decay_type = options.Branching_decay_type;
Sim_params.Branching_decay_params = options.Branching_decay_params;
%% Memory saving
global N_nodes_max
N_nodes_max = 2^32-1;

% Define the most efficient data types to use for the NodeID and BranchID.
NodeID_datatype = min_int_type(N_nodes_max,false);
BranchID_datatype = min_int_type(N_branches_max,false);
Sim_params.NodeID_format = NodeID_datatype;
Sim_params.BranchID_format = BranchID_datatype;

% Check if the maximum number of branches exceeds
% the maximal number allowed for 32 bit integer.
if N_branches_max >= intmax(BranchID_datatype)
    error('Error! The maximum number of branches is bigger than the data type allows (%s).', BranchID_datatype)
end
%% Profiling
if options.Profile
    % Save the profile several times over the duration of the simulation.
    Profile_filename = fullfile(mfilename('fullpath'), 'Profile');
    N_profiles = 10;
    Profile_recording_period = round(MaxT/N_profiles);
    
    % Start profiling.
    profile clear;
    %profile -memory on; Significantly lengthen the simulation time.
    profile on;
end
%% Random numbers
% To sample efficiently, batches of random numbers are created and used
% when necessary. When a batch is fully exhausted, a new batch of number is
% created.

% Branching angle
Branching_angle_samples = [];
Branching_i = [];
Branching_imax = 10000;

% Rebranching angle
Rebranching_angle_samples = [];
Rebranching_i = [];
Rebranching_imax = 10000;
N_rebranching_trials_max = 10;

% Initialize batches of normally and uniformly distributed random numbers
% used to sample the tip's states and velocities.

% Clear persistent variables in random sampler and initialize it.
clear random_sampler
random_sampler('rand','Initialize');
random_sampler('randn','Initialize');
%% Branch angles
global Exclusion_radius delta_theta theta_min
% Define the default branching angle distribution.
if isempty(options.Branching_angle_dist)
    options.Branching_angle_dist = makedist('Uniform',0,pi);
end

% The branching angles are defined with respect to the parent branch.
% The minimum branching angle is defined as a function of the Exclusion_radius R.
% In general, R*sin(theta_min/2) = Contact_dist/2. Then, theta_min =
% 2*asin(Contact_dist/2R). See the Mathematica file Branchpoint_Geometry.m.

% We add +1 to make sure that branches with branching angle close to
% theta_min don't collide with their sister or parent branch right when
% they exit the excluded area.
%Exclusion_radius = Contact_dist;
Exclusion_radius = 2*Contact_dist;
Exclusion_radius15 = 1.5*Exclusion_radius;
Exclusion_radius_squared = Exclusion_radius^2;
theta_min = 2*asin(Contact_dist/(2*(Exclusion_radius)));
Exclusion_radius = ceil(Exclusion_radius/delta); % Make Exclusion_radius dimensionless.
Sim_params.Exclusion_radius = Exclusion_radius;

% Truncate the branching angle distribution with the minimal angle and save it.
options.Branching_angle_dist = truncate(options.Branching_angle_dist,theta_min,pi - theta_min);
Sim_params.Branching_angle_dist = options.Branching_angle_dist;

% In general, branches are curved 2D paths with a given persistence length.
% Roughly speaking, the persistence length beta defines the lengthscale over
% which the branch starts curving. Assuming that branches are worm-like chain 
% where each link has a length delta (in this case, this corresponds to the
% size of the step), the change in the branch orientation (delta_theta) is
% normally distributed with a std = 2*delta/beta.

% Define a function that will sample a batch of growth angles.
Growth_angle_i = 1;
Growth_angle_i_max = 10000;
sample_deltatheta = @() sqrt(delta/(beta/2))*randn(Growth_angle_i_max, 1);
delta_theta = sample_deltatheta();

% Test the persistence length of a sequence of angular changes.
%test_persistence_length(reshape(delta_theta,[],10),lengthscale);
%% Tips arrays
% Initialize the tips arrays.
global Tips_State Tips_Pos Tips_NodeID
Tips_State = nan(N_tips_max, 1); % State of the tip.
Tips_Pos = nan(N_tips_max, 2); % 2D position of the tip.
Tips_NodeID = nan(N_tips_max, 1); % Node ID of the tip.
Tips_BirthPos = nan(N_tips_max, 2); % 2D birth position of the tip.
Tips_CollisionPos = nan(N_tips_max, 2); % 2D position of the last tip collision.
Tips_DeathPos = nan(N_tips_max, 2); % 2D position where the tip died.
Tips_BirthTime = nan(N_tips_max, 1); % Time that the tip was born.
Tips_DeathTime = nan(N_tips_max, 1); % Time that the tip died.
Tips_CollisionTime = nan(N_tips_max, 1); % Last time that the tip collided.
Tips_CollisionDuration = zeros(N_tips_max, 1); % Duration of the last collision.
Tips_isColliding = false(N_tips_max,1); % Boolean to indicate if the tip is colliding. A collision can occur over many time steps.
Tips_Theta = nan(N_tips_max, 1); % Angle along which the tip is growing.
Tips_CollisionAngle = nan(N_tips_max, 1); % Growth angle along which the last tip collision occurred.
Tips_NCollisions = zeros(N_tips_max, 1); % Number of collisions that the tip experienced.
Tips_NRetractions = zeros(N_tips_max, 1); % Number of complete retractions that the tip experienced.
Tips_Velocity = nan(N_tips_max, 1); % Velocity of the tip.
Tips_TransitionTime = nan(N_tips_max, 1); % Future time where the tip's state will change.
Tips_FloatingLength = zeros(N_tips_max,1); % Floating length representing the distance between the tip and the last established branch node.

% The state of each tip is identified with integers where
% 1=shrinking, 2=paused, 3=growing.
%% Tips dynamics recording
NSteps_datatype = 'int8';
N_recorded_tips = 0;
Tips_isrecorded = false(N_tips_max,1);
if options.RecordTipDynamics
    N_recorded_tips_max = 30000;
    Tips_record_length_increment = 1000;
    if options.warnings && N_tips_max > N_recorded_tips_max
        warning('The maximal number of tips (%d) is too high to record tips''s dynamics. Only the first %d tips will be recorded.',N_tips_max,N_recorded_tips_max);
    end
    
    Tips_record_Ind = nan(N_tips_max,1);
    Recorded_tips_temporal_ind = ones(N_recorded_tips_max, 1);
    Recorded_tips_record_length = Tips_record_length_increment*ones(N_recorded_tips_max,1);
    Recorded_tips_states = num2cell(zeros(N_recorded_tips_max, Tips_record_length_increment, States_datatype),2);
    Recorded_tips_NSteps = num2cell(zeros(N_recorded_tips_max, Tips_record_length_increment, NSteps_datatype),2);
else
    Recorded_tips_states = [];
    Recorded_tips_NSteps = [];
end
%% Search bins
% To check collisions, the distance between a new branch node and existent
% branch nodes must be checked to ensure that it is greater than Contact_dist.
% This implies that each branch can be thought of as a sequence of disks 
% with radius Contact_dist/2.
%
% Checking collisions with all existent nodes is expensive. To reduce
% computation time, a set of 2D bins that bin the positions of the branch
% nodes is used to only check the neigborhood of a temptative new node.
% More efficient methods using kd-trees could be implemented to perform
% this nearest-neighbour check more efficiently.

% Construct the bins structure and initialize it.
% The first bin is indexed (1,1) and is located at the top left corner.
% The binwidth is such that each bin can be occupied by a number of nodes <
% Max_N_nodes_per_bin.
global binwidth N_search_bins_per_dim Last_bin_edge_pos Bins_NodesID Branchpoint_search_range Branchpoint_radius
binwidth = 10*delta;
N_nodes_per_bin = ceil(binwidth/delta+1)*ceil(binwidth/Contact_dist+1);
N_nodes_per_bin_max = 2^8-1;
if N_nodes_per_bin > N_nodes_per_bin_max
    error(['The number of nodes per bin exceed the memory limit(', num2str(N_nodes_per_bin_max), ')']);
end
Sim_params.binwidth = binwidth;

% If a boundary is defined, use the maximum size of the boundary to define
% the size of the search bins array. Otherwise, estimate the maximal size
% of the tree from the averate velocity of the growth state.
Boundary_size_max = max(Boundary_size(:));
if ~isnan(Boundary_size_max)
    if options.Boundary_type == 1
        % Double the size of the boundary if periodic boundaries are used
        % because the branches can extend beyond the apparent boundary.
        L = 2*Boundary_size_max;
    else
        L = Boundary_size_max;
    end
else
    % To estimate the size of the bin lattice using the end-to-end distance
    % prediction of a worm-like chain of size LongestBranchLength.
    lognormal_mean = @(mu,sigma) exp(mu + sigma.^2/2);
    lognormal_var = @(mu,sigma) (exp(sigma.^2)-1).*exp(2.*mu + sigma.^2);
    %Trans_mat_cells = num2cell(Tips_params.transition_prob_matrix,[1,2]);
    %Trans_mat_P_ss = cellfun(@calculate_P_ss,Trans_mat_cells,'Uni',0);
    State_V_means = cellfun(@(V_dist) [-lognormal_mean(V_dist.mum,V_dist.sigm),V_dist.mu0,lognormal_mean(V_dist.mup,V_dist.sigp)],V_dist_fits,'Uni',0);
    %State_V_vars = cellfun(@(V_dist) [lognormal_var(V_dist.mum,V_dist.sigm),V_dist.sig0^2,lognormal_var(V_dist.mup,V_dist.sigp)],V_dist_fits,'Uni',0);
    
    Longest_branch_length = max(cellfun(@(x) x(3),State_V_means))*MaxT;
    L = sqrt(2*beta*Longest_branch_length) + Contact_dist;
end
N_search_bins_per_dim = 2*ceil(L/binwidth-1/2) + 1; % +1 such that (0,0) is at the center of the center bin.
Search_bins_size = [N_search_bins_per_dim,N_search_bins_per_dim]; % Use a square array
N_search_bins = prod(Search_bins_size);
Last_bin_edge_pos = binwidth*(ceil(L/binwidth-1/2)+1/2);

% Calculate the position of the bins' center. The first bin(1,1) is at the
% top leftmost position in the cartesian coordinate system.
Bins_edges = -Last_bin_edge_pos:binwidth:Last_bin_edge_pos;
Bins_centers_xy = zeros([Search_bins_size,2]);
[Bins_centers_xy(:,:,1),Bins_centers_xy(:,:,2)] = meshgrid(Bins_edges(1:end-1) + binwidth/2,flip(Bins_edges(2:end)) - binwidth/2);
Bins_centers_r = sqrt(sum(Bins_centers_xy.^2,3));
Bins_centers_xy(:,:,1) = Bins_centers_xy(:,:,1) + Soma_pos(1);
Bins_centers_xy(:,:,2) = Bins_centers_xy(:,:,2) + Soma_pos(2);

% Define a 3D array that stores the nodes ID in each bin.
if N_search_bins > 1e7
    warning('The number of bins (%d) is large. Large amount of memory may be needed.',N_search_bins)
end
Bins_NodesID = zeros([Search_bins_size, N_nodes_per_bin], NodeID_datatype);

% Define search bin arrays to store various info.
Bins_total_counts = []; % Total number of nodes in each bin
Bins_branchpoint_counts = []; % Total number of branchpoint nodes in each bin
Bins_tip_counts = []; % Total number of branchtip nodes in each bin
Bins_dendrites_counts = []; % Total number of dendrite nodes (or core nodes) in each bin
Bins_isempty = []; % Boolean indicating if the bin is empty.

% Define the period at which the search bins are updated.
%Bincounts_update_period = 1/timescale; % 1 min;
%Bincounts_update_period = 0.5/timescale; % 0.5 min;
%Bincounts_update_period = floor(0.2/timescale); % 0.2 min; 
Bincounts_update_period = 2; % Every 2 time steps.
Sim_params.Bincounts_update_period = Bincounts_update_period;

% Send warning if the update period is too slow. A slow update period may
% cause collisions to be missed.
if Bincounts_update_period > 3
    warning('The update period (%d) of the search bins is too high. Collision checks may miss some collisions.',Bincounts_update_period);
end

% Create a branch point mask that will be used to store the bins contained
% within the branch point area (circle centered at branch point with radius
% Branchpoint_radius). Another branch point cannot appear in this area.
Branchpoint_radius = Contact_dist;
Branchpoint_radius_squared = Branchpoint_radius^2;
Sim_params.Branchpoint_radius = Branchpoint_radius;
Branchpoint_bin_search_radius = ceil(Branchpoint_radius/binwidth);
Branchpoint_search_range = -Branchpoint_bin_search_radius:Branchpoint_bin_search_radius;
%% Collisions recording
if options.RecordCollisions
    % Initialize the structure that records the collisions' data.
    N_recorded_collisions_max = N_tips_max;
    Collisions = struct();
    Collisions.TipID = nan(N_recorded_collisions_max, 1);
    Collisions.Type = nan(N_recorded_collisions_max, 1);
    Collisions.Pos = nan(N_recorded_collisions_max, 2);
    Collisions.Time = nan(N_recorded_collisions_max, 1);
    Collisions.Angle = nan(N_recorded_collisions_max, 1);
    Collisions.Length = nan(N_recorded_collisions_max, 1);
else
    Collisions = [];
end

% Count the number of collisions.
Collisions_counter = 0; 
Collision_counter_prev = 0;

% Count the number of branches that fully retracted.
Retractions_counter = 0;
Retractions_counter_prev = 0;
%% Branches arrays
% Define global arrays that store the branch properties.
global Branches_isDynamic Branches_isActive Branches_BranchingAngle Branches_TipID ...
    Branches_ID Branches_NodesID Branches_ParentID ...
    Branches_SiblingID Branches_ChildrenID Branches_Length

Branches_ID = 1:N_branches_max; % ID of each branch.
Branches_NodesID = zeros(N_branches_max, N_nodes_per_branch, NodeID_datatype); % Nodes ID belonging to each branch
Branches_ParentID = zeros(N_branches_max, 1, BranchID_datatype); % ID of the parent branch in the tree hierarchy.
Branches_SiblingID = zeros(N_branches_max, 1, BranchID_datatype); % ID of the sibling branch in the tree hierarchy.
Branches_ChildrenID = zeros(N_branches_max, 2, BranchID_datatype); % IDs of the children branches in the tree hierarchy.
Branches_Length = zeros(N_branches_max, 1); % Length of the branch.
Branches_isDynamic = false(N_branches_max, 1); % Boolean indicating if the branch is moving.
Branches_isActive = false(N_branches_max, 1); % Boolean indicating if the branch is active. A branch becomes inactive when it disappears permanently due to a complete retraction.
Branches_TipID = zeros(N_branches_max, 1); % ID of the tip associated to each branch.
Branches_BranchingAngle = zeros(N_branches_max, 1); % Angle at which each branch initially grew.

% Array to store branches that have fully retracted.
Retracted_branches_ID = zeros(1, floor(N_branches_max/10));
Retracted_branches_counter = 0;

% Array to store branches that have rebranched.
Rebranched_branches_ID = zeros(1, floor(N_branches_max/10));
Rebranched_branches_counter = 0;

% Make sure that all branch arrays are recycled to reduce the memory footprint.
% The variable name of branches arrays must match the pattern 'Branches*';
s = whos;
s = s(cellfun(@(x) ~isempty(x) && x==1,regexp({s.name},'Branches*')'));
Branches_arrays_name = {s.name};
Recyled_branches_arrays = {'Branches_ID','Branches_NodesID','Branches_ParentID','Branches_SiblingID',...
    'Branches_ChildrenID','Branches_Length','Branches_isDynamic','Branches_isActive',...
    'Branches_TipID','Branches_BranchingAngle'};
is_array_recycled = ismember(Branches_arrays_name,Recyled_branches_arrays);
if any(~is_array_recycled)
    err_msg = sprintf('The following branches arrays are not recycled:\n');
    err_msg = [err_msg sprintf('%s',Branches_arrays_name{~is_array_recycled})];
    error(err_msg);
end
%% Nodes arrays
% Define global arrays that store the branch node properties.
global Nodes_counter Nodes_freeindices Nodes_Pos Nodes_BirthTime Nodes_isActive ...
    Nodes_isBranchpoint Nodes_isNextToBP Nodes_BranchID Nodes_BinInd
N_nodes = min(ceil(MaxT*N_branches_max), 1000000);
Nodes_Pos = zeros(N_nodes, 2); % 2D position of each node.
Nodes_BirthTime = zeros(N_nodes, 1); % Time when each node appeared.
Nodes_isActive = false(N_nodes, 1); % Boolean indicating if each node is active. Inactive nodes denote nodes that have been deleted due to branch retractions.
Nodes_isBranchpoint = false(N_nodes, 1); % Boolean indicating if node is a branchpoint.
Nodes_isNextToBP = false(N_nodes, 1); % Boolean indicating if node is within a distance Branchpoint_radius to a branchpoint.
Nodes_ParentID = zeros(N_nodes, 1); % ID of the parent node in the tree hierarchy. 
Nodes_PathLength = zeros(N_nodes, 1); % Path length between each node and the soma.
Nodes_BranchID = zeros(N_nodes, 1, BranchID_datatype); % ID of the branch to which each node belongs.
Nodes_BinInd = zeros(N_nodes, 3); % 3D indices of the search bins where each node belongs.

Nodes_counter = 1;
Nodes_freeindices = 1:N_nodes;
%% Initialize tree
% The tree is initialized with N_soma_branches located at the soma (0,0).
% Soma branches are oriented at equal angular distance from one another
% such that they cover a range of 360 deg.
% Soma branches are in the growth state until branching starts. This 
% ensures that the tree will not disappear due to random fluctuations of
% the initial tips.

% In addition, ensure that the initial soma tips grow at a minimum positive velocity.
% This is necessary to avoid rare cases where all soma branches disappear.
Sim_params.SomaInitialVelocity = options.SomaInitialVelocity;
Sim_params.SomaInitialVelocity_unit = 'um/min';
Soma_initial_velocity = round(options.SomaInitialVelocity/lengthscale*timescale);

% Initialize the global time index.
global t Soma_branches_ID
t = 0;
Soma_branches_ID = 1:N_soma_branches;
for i = Soma_branches_ID
    Branches_ParentID(i) = 0;
    Branches_SiblingID(i) = 0;
    Branches_NodesID(i, 1) = 1;
    Branches_TipID(i) = i;
    Branches_isDynamic(i) = true;
    Branches_isActive(i) = true;

    % Initialize Tips arrays
    Tips_BirthTime(i) = 1;
    Tips_Pos(i, :) = Soma_pos;
    Tips_NodeID(i) = 1;
    Tips_Theta(i) = 0+(i-1)*2*pi/N_soma_branches;
    Branches_BranchingAngle(i) = Tips_Theta(i);
    
    % Assign the initial growth tracks for soma branches.
    Tips_State(i) = 3;
    Tips_Velocity(i) = Soma_initial_velocity;
    Tips_TransitionTime(i) = Branching_start_time;
end

% Initialize Nodes arrays with the soma node.
Soma_nodeID = add_node(Soma_pos, 0, 0);
Nodes_isBranchpoint(Soma_nodeID) = true;
Nodes_isNextToBP(Soma_nodeID) = true;

% Initialize the number of tips and the branch counter.
N_tips = N_soma_branches;
N_branches = N_soma_branches;
%% Initialize timeseries
% Define a structure to keep track of the timeseries data.
if options.Timeseries
    Timeseries_Lengths = cell(1, MaxT);
    Timeseries_Statecounts = zeros(3, MaxT, 'uint16');
    Timeseries_NBranchpoints = zeros(1, MaxT, 'uint16');
    Timeseries_NNewBranchpoints = zeros(1, MaxT, 'uint16');
    Timeseries_NCollisions = zeros(1, MaxT, 'uint16');
    Timeseries_NRetractions = zeros(1, MaxT, 'uint16');
    Timeseries_NDeaths = zeros(1, MaxT, 'uint16');
    Timeseries_Rg = zeros(1, MaxT);
    Timeseries_D_uniform = zeros(MaxT,2);
    if options.ConvexHull
        Timeseries_CH_Diameter = zeros(1, MaxT);
        Timeseries_CH_Area = zeros(1, MaxT);
    else
        Timeseries_CH_Diameter = [];
        Timeseries_CH_Area = [];
    end

    % Radial density measurements.
    Radial_hist_binsize = 20/lengthscale; % 20 um
    Radial_hist_Rmax = 1000/lengthscale; % 1000 um
    Radial_hist_binedges = 0:Radial_hist_binsize:Radial_hist_Rmax;
    Radial_hist_binedges2 = Radial_hist_binedges.^2;
    Radial_hist_update_period = 1/timescale; % 1 min
    
    N_radial_hist_updates = ceil(MaxT/Radial_hist_update_period);
    Timeseries_dendrites_counts_r = nan(numel(Radial_hist_binedges)-1, N_radial_hist_updates);
    Timeseries_branchpoints_counts_r = nan(numel(Radial_hist_binedges)-1, N_radial_hist_updates);
    Timeseries_branchtips_counts_r = nan(numel(Radial_hist_binedges)-1, N_radial_hist_updates);
end

% Timeseries of the requested metrics.
if ~isempty(options.Timeseries_metric)
    Timeseries_metric_sampled_time_ind = round(options.Timeseries_metric_times/timescale);
    Timeseries_metric_sample_ind = 1;
    Timeseries_metric_N_samples = numel(Timeseries_metric_sampled_time_ind);
    
    % Define the cell array of options that is given to the analyze_tree
    % function.
    Metrics_cell = [strrep(options.Timeseries_metric,'meshsize','hitprob'); num2cell(ones(1,numel(options.Timeseries_metric)))];
    
    Timeseries_df = zeros(Timeseries_metric_N_samples,1);
    Timeseries_df_rad = cell(Timeseries_metric_N_samples,1);
    Timeseries_df_C = cell(Timeseries_metric_N_samples,1);
    Timeseries_meshsize = zeros(Timeseries_metric_N_samples,1);
end

% Timeseries of the entire tree structure
if ~isempty(options.Timeseries_tree_times)
    Timeseries_Tree_sampled_time_ind = round(options.Timeseries_tree_times/timescale);
    Timeseries_Tree_sample_ind = 1;
    Timeseries_Tree_N_samples = numel(options.Timeseries_tree_times);
    if Timeseries_Tree_N_samples > 10
        warning('The number of saved tree structures is %d. This may increase memory demands',Timeseries_Tree_N_samples);
    end
    Timeseries_Tree = cell(Timeseries_Tree_N_samples,1);
end

% Initialize schedule to record tip dynamics tracks in given intervals.
is_timeseries_tracks_recorded = false;
Timeseries_tracks_interval_ind = 0;
if ~isempty(options.Timeseries_tracks)
    Timeseries_tracks_record_time_start = round(options.Timeseries_tracks.Record_times_start/timescale);
    Timeseries_tracks_record_time_end = round(options.Timeseries_tracks.Record_times_end/timescale);
    Timeseries_tracks_N_intervals = numel(Timeseries_tracks_record_time_start);
    
    % Ensure that the record end times are greater than the start times.
    assert(all(Timeseries_tracks_record_time_start < Timeseries_tracks_record_time_end),'All recording intervals must be non-empty.');
    
    % Initialize the structure that will store the timeseries
    Timeseries_tracks = struct();
    Timeseries_tracks.Tracks_L = {};
    Timeseries_tracks.Tracks_t = {};
    Timeseries_tracks.Tracks_length = [];
    Timeseries_tracks.Tracks_hascollided = {};
    
    % Define dictionary to map Tip ID to Track ID.
    Timeseries_tracks.Tracks_ID = containers.Map('KeyType','double','ValueType','double');
else
    Timeseries_tracks_N_intervals = 0;
end
%% Initialize movie.
if options.Movie
    % Find the time index of the frames that are recorded.
    if strcmp(options.Movie_recorded_frames,'All')
        % Record the tree at a given period to reduce memory.
        MovieRecordingPeriod = round(5/timescale); % 5 min
        Movie_RecordedFramesInd = 0:MovieRecordingPeriod:MaxT; % Initial structure is recorded.
    else
        Movie_RecordedFramesInd = options.Movie_recorded_frames;
    end
    
    N_movieframes = numel(Movie_RecordedFramesInd);
    Movieframes = struct('tree', cell(N_movieframes, 1), 'time', cell(N_movieframes, 1));
    Movie_frame_counter = 1;
    Movie_next_frame_ind = Movie_RecordedFramesInd(Movie_frame_counter);

    % Record the zeroth frame, the frame of the initial tree structure before the time loop.
    if Movie_next_frame_ind == 0
        Tree = build_tree();
        Movieframes(Movie_frame_counter).tree = restructure_tree(Tree, 'RemoveZeroLength', false,'RemoveInactive',true);
        Movieframes(Movie_frame_counter).time = t*timescale;
        
        % Add position of boundary (if any).
        if is_boundary_existent
            Movieframes(Movie_frame_counter).boundary_size = Boundary_size(1,:)*lengthscale;
            Movieframes(Movie_frame_counter).boundary_shape = options.Boundary_shape;
            Movieframes(Movie_frame_counter).boundary_type = options.Boundary_type;
        end
        
        % Move the counter to the next frame.
        Movie_frame_counter = Movie_frame_counter + 1;
        Movie_next_frame_ind = Movie_RecordedFramesInd(Movie_frame_counter);
    end
else
    Movieframes = [];
end
%% Time loop

% Initialize waitbar
if options.waitbar
    wb = waitbar(0, [mfilename(), ' simulation']);
end

% Start simulation timer.
tic;
for t = 1:MaxT
    %% Save profile.
    if options.Profile && mod(t,Profile_recording_period) == 0
        save_profile(Profile_filename);
        vars = whos;
        [~,ind] = sort([vars.bytes],'descend');
        vars = vars(ind);
        fprintf('Total Memory Usage: %.2f GB\n',sum([vars.bytes])/1e9);
    end
    %% Handle growth and retraction of each dynamic branch in this timeframe.
    % Randomize the growth order of the dynamic tips.
    Dynamic_branches_ID = Branches_ID(Branches_isDynamic(1:N_branches));
    N_dynamic_branches = numel(Dynamic_branches_ID);
    Growth_order = randperm(N_dynamic_branches);

    % Update the position of the boundary.
    if is_boundary_existent
        Boundary_pos_temp = Boundary_positions(t,:);
        Boundary_pos_squared_temp = Boundary_pos_temp.^2;
        Boundary_size_temp = Boundary_size(t,:);
    end
    
    % Check if growth tracks timeseries are recorded in this timeframe.
    if is_timeseries_tracks_recorded
        % If growth tracks were recorded in the previous frame, check if
        % they are still recorded in this frame.
        is_timeseries_tracks_recorded = t < Timeseries_tracks_record_time_end(Timeseries_tracks_interval_ind);
    else
        % If growth tracks were NOT recorded in the previous frame, check if
        % they are recorded in this frame.
        is_timeseries_tracks_recorded = Timeseries_tracks_interval_ind < Timeseries_tracks_N_intervals && t >= Timeseries_tracks_record_time_start(Timeseries_tracks_interval_ind + 1);
        if is_timeseries_tracks_recorded
            % If recording is turned on in this frame, increase the
            % interval index.
            Timeseries_tracks_interval_ind = Timeseries_tracks_interval_ind + 1;
        end
    end
    
    % Loop through each dynamic branch and simulate its dynamics.
    for ii = 1:N_dynamic_branches
        % Get info about the current moving branch.
        BranchID = Dynamic_branches_ID(Growth_order(ii));
        TipID = Branches_TipID(BranchID);
        Branch_length_initial = Branches_Length(BranchID);
        Branchpoint_NodeID = Branches_NodesID(BranchID, 1);
        Branchpoint_pos = Nodes_Pos(Branchpoint_NodeID, :);
        
        % Check if the tip is colliding. If it is, check if the collision
        % duration has ended. If collision has end, switch the state of the
        % tip to the defined post-collision state.
        if Tips_isColliding(TipID) && t > Tips_CollisionTime(TipID) + Tips_CollisionDuration(TipID)
            % Change the state of the tip to the user-defined collision state.
            change_tip_state(TipID,options.Collision_state);
            Tips_isColliding(TipID) = false;
            Tips_CollisionDuration(TipID) = 0;
            
            % Add the bias on the duration and velocity of the
            % post-collision state.
            if options.Collision_t_bias ~= 0
                Tips_TransitionTime(TipID) = Tips_TransitionTime(TipID) + options.Collision_t_bias;
            end
            if options.Collision_v_bias ~= 0
                Tips_Velocity(TipID) = Tips_Velocity(TipID) + options.Collision_v_bias;
            end
        end
        
        % Determine if the tip's state has changed in this timeframe.
        Tip_TransitionTime = Tips_TransitionTime(TipID);
        if Tip_TransitionTime <= t
            % If the tip state has changed, record the proportion of the
            % time spent in the previous state.
            State_duration_prev = Tip_TransitionTime - (t-1);
            State_velocity_prev = Tips_Velocity(TipID);
            
            % Perform a state transition following the transition matrix.
            change_tip_state(TipID);
        else
            % No transition has occurred.
            State_duration_prev = 0;
            State_velocity_prev = 0;
        end
        
        % Determine the change in branch length given the tip dynamics model.
        Tip_Velocity = Tips_Velocity(TipID);
        Tip_FloatingLength = Tips_FloatingLength(TipID);
        switch Tips_params.TipDynamicsModel
            case '3state'
                % Determine the length change of the branch based on the tip
                % velocity.
                Delta_L = Tip_FloatingLength + State_duration_prev*State_velocity_prev + (1-State_duration_prev)*Tip_Velocity;
            case 'drift_diff'
                % In this case, simply use the tips velocity to determine the
                % change in length. The effect of diffusion must be added
                % on top of the tip's drift velocity.
                
                % Add diffusion transport after the initial growth phase.
                if Tips_params.Diffusion_displacement_sigma > 0 && t >= Tips_BirthTime(TipID) + Initial_growth_duration
                    Randn_num = random_sampler('randn');
                    Delta_L = Tip_FloatingLength + Tip_Velocity + Tips_params.Diffusion_displacement_sigma * Randn_num;
                else
                    Delta_L = Tip_FloatingLength + Tip_Velocity;
                end
        end
        
        % Determine the number of steps taken. The remainder is
        % saved in Tips_FloatingLength and added in the next time step.
        N_steps = floor(Delta_L);
        Tips_FloatingLength(TipID) = Delta_L - N_steps;
        %% Find neighboring nodes that are excluded in the collision check.
        % If the branch length is < Exclusion_radius, find the nodes that
        % are at most a path length of Exclusion_radius15 away from the branch
        % point. Among these nodes, remove those that belong to either the 
        % sibling or parent branch. Omitting these nodes removes spurious
        % collisions with the sibling or parent due to the fact that the 
        % stepsize may be smaller than Contact_dist.
        
        % Calculate the distance between the tip and its branch point.
        if Branch_length_initial <= Exclusion_radius15
            Neigh_parent_nodes_ID = [];
            Neigh_sibling_nodes_ID = [];
            SiblingID = Branches_SiblingID(BranchID);
            ParentID = Branches_ParentID(BranchID);

            if BranchID > N_soma_branches
                % Find the neighboring parent's nodes.
                Parent_Length = Branches_Length(ParentID);
                if Parent_Length > 0
                    Neigh_parent_nodes_ID = Branches_NodesID(ParentID, max(Parent_Length+1-2*Exclusion_radius, 1):Parent_Length+1)';
                end

                % Remove the neighboring sister's nodes.
                Sibling_Length = Branches_Length(SiblingID);
                if Sibling_Length > 0
                    Neigh_sibling_nodes_ID = Branches_NodesID(SiblingID, 1:min(Sibling_Length, 2*Exclusion_radius)+1)';
                end
            else
                % Soma branches only have siblings.
                BranchIDs = [1:BranchID-1, BranchID+1:N_soma_branches];
                Parent_Length = 0;

                for kk = BranchIDs
                    Sibling_Length = Branches_Length(kk);
                    if Sibling_Length > 0
                        Neigh_sibling_nodes_ID = [Neigh_sibling_nodes_ID; Branches_NodesID(kk, 1:min(Sibling_Length+1, 2*Exclusion_radius))'];
                    end
                end
            end

            % Define the excluded nodes as the parent and sister nodes that
            % are within a distance Exclusion_radius from the common branchpoint.
            if ~isempty(Neigh_parent_nodes_ID) || ~isempty(Neigh_sibling_nodes_ID)
                Neigh_ParSis_NodesID = [Neigh_parent_nodes_ID; Neigh_sibling_nodes_ID];
                Neigh_ParSis_Nodes_Pos = Nodes_Pos(Neigh_ParSis_NodesID, :) - Branchpoint_pos;
                Neigh_ParSis_Nodes_SquaredDist = sum(Neigh_ParSis_Nodes_Pos.^2, 2);
                Neigh_ParSis_Nodes_isexcluded = Neigh_ParSis_Nodes_SquaredDist < Exclusion_radius_squared;
                Neigh_ParSis_NodesID_excluded = Neigh_ParSis_NodesID(Neigh_ParSis_Nodes_isexcluded);
            else
                Neigh_ParSis_NodesID_excluded = [];
            end
        else
            Neigh_ParSis_NodesID_excluded = [];
        end
        
        %% Grow each segment separately and perform collision check.
        % When N_steps is negative, the for loop will be skipped.
        Collision_branch = false;
        Collision_boundary = false;
        Collision_soma = false;
        for kk = 1:N_steps
            % Get the tip position and node ID of the tip.
            Tip_Pos = Tips_Pos(TipID, :);

            % Send warning if the future length surpasses the array
            % dimension of Branch_PointsBinInd. Surpassing this limit will
            % slow down execution since the array will be reshaped.
            Branch_length_future = Branch_length_initial + kk;
            if options.warnings && Branch_length_future > N_nodes_per_branch
                warning('Branch %d (L=%d) is greater than N_nodes_per_branch (%d)',BranchID,Branch_length_future,N_nodes_per_branch);
            end

            % Find the nodesID that should be omitted when considering collisions.
            startpos = max(Branch_length_future - 3*ceil(Contact_dist), 1);
            ExcludedNodesID = Branches_NodesID(BranchID, startpos:Branch_length_future)';
            
            % Add the the sister and parent segments if the branch length is smaller than the exclusion radius.
            if Branch_length_future <= Exclusion_radius15
                ExcludedNodesID = [ExcludedNodesID; Neigh_ParSis_NodesID_excluded];
            end

            % Determine the angle of growth.
            if Branch_length_future > 0
                if Growth_angle_i > Growth_angle_i_max
                    delta_theta = sample_deltatheta();
                    Growth_angle_i = 1;
                end
                Theta = Tips_Theta(TipID) + delta_theta(Growth_angle_i);
                Growth_angle_i = Growth_angle_i + 1;
            else
                Theta = Tips_Theta(TipID);
            end

            % Determine the new 2D position of the tip.
            Tip_Pos_new = Tip_Pos + delta * [cos(Theta), sin(Theta)];
            
            % Calculate the position of the new tip relative to the
            % soma.
            Tip_Pos_new_centered = Tip_Pos_new - Soma_pos;
            Tip_Pos_new_soma_radius_squared = sum(Tip_Pos_new_centered.^2,2);

            % Check collision with the boundary (if any).
            if is_boundary_existent && (options.Boundary_type == 0 || options.Boundary_type == -1)
                switch options.Boundary_shape
                    case {'square','rectangle'}
                        Collision_boundary = any(Tip_Pos_new_centered.^2 > Boundary_pos_squared_temp);
                    case 'circle'
                        Collision_boundary = Tip_Pos_new_soma_radius_squared > Boundary_pos_squared_temp;
                end

                % Stop the growth if the tip has hit the boundary.
                if Collision_boundary
                    if options.Boundary_type == -1
                        % If the boundary is repulsive, change the tip to
                        % the shrinking state (1).
                        change_tip_state(TipID, 1);
                    end

                    % Stop further growth after a collision with the boundary.
                    break;
                end
            end

            % Check if the growth of non-Soma branches enter the Soma exclusion zone.
            Collision_soma = BranchID > N_soma_branches && Tip_Pos_new_soma_radius_squared < Soma_radius_squared;
            if Collision_soma
                break;
            end

            %% Check collisions with neighboring nodes.
            Segment_new = [Tip_Pos; Tip_Pos_new];
            [Collision_branch, Collided_nodes_ID, Collisions_time] = check_nodes_collision(Segment_new, ExcludedNodesID);
            
            % Check that no collision has happened if the current path length
            % is less than Exclusion_radius.
            if Branch_length_future <= Exclusion_radius && Collision_branch && BranchID <= N_soma_branches
                warning('Soma branch %d has collided with branches:\n%s', BranchID, sprintf('%d\n',Nodes_BranchID(Collided_nodes_ID)));
            end

            % Continue growth if there is no collision.
            % Otherwise, stop the growth for this timeframe.
            if ~Collision_branch
                % Add the new node to the nodes and bins arrays.
                Parent_node_ID = Branches_NodesID(BranchID, Branch_length_future);
                New_node_ID = add_node(Tip_Pos_new, BranchID, Parent_node_ID);

                % Update the branches arrays.
                Branches_Length(BranchID) = Branch_length_future;
                Branches_NodesID(BranchID, Branch_length_future+1) = New_node_ID;

                % Update the tips arrays.
                Tips_NodeID(TipID) = New_node_ID;
                Tips_Pos(TipID, :) = Tip_Pos_new;
                Tips_Theta(TipID) = Theta;
            else
                if all(Collisions_time < 0)
                    % If all collision times are negative, this indicates that
                    % the tip position was already in a colliding position.
                    % This can be due to a bad branching angle or a neigbouring
                    % tip that grew in the vicinity, while the current tip was
                    % growing without collision checks. In these cases,
                    % set the collision time to 0.
                    Collision_time = 0;
                    Collision_type = 0;
                else
                    % Find the closest positive collision time.
                    Collision_time = min(Collisions_time(Collisions_time >= 0));
                    Collision_type = 1;
                end
                
                % Test if Collision_time is calculated properly.
                if isempty(Collision_time) || isnan(Collision_time) || ~isreal(Collision_time) || Collision_time > 1
                    PlotTree();
                    error('The tip collision time was not calculated properly.')
                end
                
                % Set the floating length of the tip to the amount of
                % length that was grown before collision.
                Tips_FloatingLength(TipID) = Collision_time*delta;
                
                % Increase the collision counter.
                Tips_NCollisions(TipID) = Tips_NCollisions(TipID)+1;
                
                % Calculate the true collision position and time.
                TrueCollisionPos = (1-Collision_time)*Tip_Pos + Collision_time*Tip_Pos_new;
                TrueCollisionTime = t + Collision_time;
                
                % Record the collision time in the tips array.
                Tips_CollisionTime(TipID) = TrueCollisionTime;
                
                % Record the position and time of the collision.
                Collisions_counter = Collisions_counter+1;
                if options.RecordCollisions && Collisions_counter < N_recorded_collisions_max
                    Collisions.TipID(Collisions_counter) = TipID;
                    Collisions.Type(Collisions_counter) = Collision_type;
                    Collisions.Pos(Collisions_counter, :) = TrueCollisionPos;
                    Collisions.Time(Collisions_counter) = TrueCollisionTime;
                    Collisions.Angle(Collisions_counter) = mod(Theta, 2*pi);
                    Collisions.Length(Collisions_counter) = Branch_length_future-1;
                end
                
                % Determine the duration of the collision.
                Collision_duration = Collision_duration_sampler();
                Tips_CollisionDuration(TipID) = Collision_duration;
                
                % Change the state of the tip to the user-defined collision 
                % state if the collision duration is less than 1 frame.
                if Collision_time + Collision_duration < 1
                    change_tip_state(TipID, options.Collision_state);
                else
                    % In this case, the collision lasts more than 1
                    % timeframe. Delay the post-collision dynamics until 
                    % the collision has ended.
                    Tips_isColliding(TipID) = true;
                    
                    % Pause the growth of the tip during collision.
                    Tips_State(TipID) = 2;
                    Tips_Velocity(TipID) = 0;
                    Tips_TransitionTime(TipID) = inf; % While the tip is colliding, no spontaneous transitions can occur.
                end
                
                % Stop the growth of the remaining segments after collision.
                break;
            end
        end
        %% Determine the final length.
        % There are 4 cases to handle:
        % 1) a negative N_steps was chosen
        % 2) a collision happened when the growth wasn't finished
        % 3) the growth completed without collisions.
        % 4) the growth stopped because the tip entered an exclusion zone.
        has_tip_collided = Collision_branch || Collision_boundary || Collision_soma;
        if has_tip_collided
            Branch_length_final = Branch_length_initial + kk - 1;
        else
            Branch_length_final = Branch_length_initial + N_steps;
        end
        
        % Impose lower bound on the final length.
        if BranchID <= N_soma_branches
            % Ensure that the soma branches survive at all times.
            Branch_length_final = max(Branch_length_final,1);
        else
            Branch_length_final = max(Branch_length_final,0);
        end
        
        % Update the tips arrays if the branch length changed.
        if Branch_length_final ~= Branch_length_initial
            Tip_NodeID_new = Branches_NodesID(BranchID, Branch_length_final+1);
            Tip_Pos_new = Nodes_Pos(Tip_NodeID_new, :);
            Tips_NodeID(TipID) = Tip_NodeID_new;
            Tips_Pos(TipID, :) = Tip_Pos_new;
        end
        
        % In cases 1) and 2), delete branch nodes of deleted segments and
        % determine the new angle of growth.
        if N_steps < 0
            % Check if the branch has fully retracted.
            if Branch_length_final == 0
                % If the branch has reached back to the branch point, handle
                % the various full-retraction scenarios.
                
                % Reset the tip floating length to 0 since it has come back to its original branchpoint.
                Tips_FloatingLength(TipID) = 0;
                
                % If the tip has shrunk back to the branch point, handle the
                % complete retraction of the branch.
                Retractions_counter = Retractions_counter + 1;
                
                % Check if the tip hasn't reached the maximum number of
                % retractions.
                Tips_NRetractions(TipID) = Tips_NRetractions(TipID)+1;
                if Tips_NRetractions(TipID) > N_retractions_max
                    error(['Tip ', num2str(TipID), ' is still active and has fully retracted beyond the maximum number of full retractions.'])
                end

                % Determine the full retraction scenario.
                % Based on the maximum allowed number of full retractions
                % and the rebranching probability, the tip will either pause
                % or get deleted.
                is_rebranching = Rebranching_prob > 0 && ...
                                random_sampler('rand') < Rebranching_prob && ...
                                Tips_NRetractions(TipID) < N_retractions_max;
                if is_rebranching
                    Retraction_scenario_temp = Rebranching_scenario;
                else
                    Retraction_scenario_temp = 'delete_branches';
                end

                switch Retraction_scenario_temp
                    case 'pause_tip'
                        % Pause the tip.
                        change_tip_state(TipID, 2);
                        Tips_Velocity(TipID) = 0;
                        Rebranched_branches_counter = Rebranched_branches_counter+1;
                        Rebranched_branches_ID(Rebranched_branches_counter) = BranchID;
                    case 'pause_tip_in_dying_state'
                        % Pause the tip in a dying state. In such a state,
                        % the branch tip does not reach back to the branchpoint. 
                        % The branch length is minimal.
                        change_tip_state(TipID, 2);
                        Tips_Velocity(TipID) = 0;
                        Branch_length_final = 1;
                    case 'delete_branches'
                        % Delete branch.
                        Tips_DeathTime(TipID) = t;
                        if options.Timeseries
                            Timeseries_NDeaths(t) = Timeseries_NDeaths(t) + 1;
                        end
                        Retracted_branches_counter = Retracted_branches_counter+1;
                        Retracted_branches_ID(Retracted_branches_counter) = BranchID;
                    otherwise
                        error('The retraction scenario (%s) is undefined.',Retraction_scenario_temp);
                end
            end
            
            % Calculate the new growing angle Theta of the branch if it has
            % not reached back to the branchpoint.
            if Branch_length_final > 0
                Tip_vec_NodeIDs = Branches_NodesID(BranchID, Branch_length_final:Branch_length_final+1);
                Tip_vec = diff(Nodes_Pos(Tip_vec_NodeIDs, :));
                Tips_Theta(TipID) = atan2(Tip_vec(2), Tip_vec(1));
            end
            
            % Remove the branch nodes from the bin arrays.
            Removed_nodes_ID = Branches_NodesID(BranchID, Branch_length_final+2:Branch_length_initial+1);
            delete_nodes(Removed_nodes_ID);
        end
        
        % Update the branch length in the branches array.
        Branches_Length(BranchID) = Branch_length_final;
        %% Record the dynamics of the tip.
        if Tips_isrecorded(TipID)
            Record_ind = Tips_record_Ind(TipID);
            Temporal_ind = Recorded_tips_temporal_ind(Record_ind);
            
            % If the record index has reached the end of the record array,
            % increase the size of the arrays.
            Record_length = Recorded_tips_record_length(Record_ind);
            if Temporal_ind > Record_length
                Recorded_tips_NSteps{Record_ind} = [Recorded_tips_NSteps{Record_ind} zeros(1,Tips_record_length_increment,States_datatype)];
                Recorded_tips_states{Record_ind} = [Recorded_tips_states{Record_ind} zeros(1,Tips_record_length_increment,NSteps_datatype)];
                Recorded_tips_record_length(Record_ind) = numel(Recorded_tips_NSteps{Record_ind});
            end
            
            % Record the state and number of steps taken.
            Recorded_tips_states{Record_ind}(Temporal_ind) = Tips_State(TipID);
            Recorded_tips_NSteps{Record_ind}(Temporal_ind) = Branch_length_final - Branch_length_initial;
            Recorded_tips_temporal_ind(Record_ind) = Temporal_ind + 1;
        end       
        %% Add tips growth tracks, if growth tracks timeseries are recorded.
        if is_timeseries_tracks_recorded
            % Find the track index of the current tip.
            Delta_length_final = Branch_length_final - Branch_length_initial;
            if ~isKey(Timeseries_tracks.Tracks_ID,TipID)
                % Initialize a new track.
                Track_ID = double(Timeseries_tracks.Tracks_ID.Count) + 1;
                Timeseries_tracks.Tracks_ID(TipID) = Track_ID;
                Timeseries_tracks.Tracks_t{Track_ID} = [t-1; t; nan(100,1)];
                Timeseries_tracks.Tracks_L{Track_ID} = [Branch_length_initial; Delta_length_final; nan(100,1)];
                Timeseries_tracks.Tracks_hascollided{Track_ID} = [false; has_tip_collided; false(100,1)];
                Timeseries_tracks.Tracks_length(Track_ID) = 2;
            else
                % Add new record in the existing track.
                Track_ID = Timeseries_tracks.Tracks_ID(TipID);
                Track_Length = Timeseries_tracks.Tracks_length(Track_ID);
                
                % Increase size of track's array if all space is used.
                if Track_Length >= numel(Timeseries_tracks.Tracks_t{Track_ID})
                    Timeseries_tracks.Tracks_t{Track_ID} = [Timeseries_tracks.Tracks_t{Track_ID}; nan(100,1)];
                    Timeseries_tracks.Tracks_L{Track_ID} = [Timeseries_tracks.Tracks_L{Track_ID}; nan(100,1)];
                    Timeseries_tracks.Tracks_hascollided{Track_ID} = [Timeseries_tracks.Tracks_hascollided{Track_ID}; false(100,1)];
                end
                Timeseries_tracks.Tracks_t{Track_ID}(Track_Length+1) = t;
                Timeseries_tracks.Tracks_L{Track_ID}(Track_Length+1) = Timeseries_tracks.Tracks_L{Track_ID}(Track_Length) + Delta_length_final;
                Timeseries_tracks.Tracks_hascollided{Track_ID}(Track_Length+1) = has_tip_collided;
                Timeseries_tracks.Tracks_length(Track_ID) = Track_Length + 1;
            end
        end
        %% Debugging tests.
        if Debugmode
            Test('Growth');
            Test('180Reversal');
        end
    end
    %% Post-growth handling of rebranched and deleted branches.
    % Count the number of tips in each possible state before they
    % get deactivated by @deletebranch.
    if options.Timeseries
        Dynamic_tips_ID = Branches_TipID(Dynamic_branches_ID);
        Tip_state_counts = histc(Tips_State(Dynamic_tips_ID), 1:N_states);
        Timeseries_Statecounts(:, t) = Tip_state_counts;
    end

    % Calculate the rebranching angles of branches that rebranched.
    if Rebranched_branches_counter > 0
        calculate_rebranching_angle(Rebranched_branches_ID(1:Rebranched_branches_counter));
    end

    % Delete branches that have retracted back to a branch point.
    if Retracted_branches_counter > 0
        deletebranch(Retracted_branches_ID(1:Retracted_branches_counter));
    end
    
    if Debugmode
        % Check that all active tips have full retractions below the maximum.
        Active_branches_ID = Branches_ID(Branches_isDynamic(1:N_branches));
        Active_tips_ID = Branches_TipID(Active_branches_ID);
        if any(Tips_NRetractions(Active_tips_ID) >= N_retractions_max)
            error('Some tips have exceeded the maximum number of full retractions.')
        end
    end
    %% Update tip dynamics parameters
    % Update the transition matrix, if it is not constant in time.
    is_transition_matrix_updated = strcmp(Tips_params.TipDynamicsModel,'3state') ...
                                   && ~is_transition_rate_matrix_constant ...
                                   && t == Transition_matrix_update_time;
    if is_transition_matrix_updated
        Transition_matrix_update_ind = Transition_matrix_update_ind + 1;
        Tips_params.transition_prob_matrix = Tips_params_in.transition_rate_matrix(:,:,Transition_matrix_update_ind) * timescale;
        Transition_matrix_update_time = Transition_matrix_update_times(Transition_matrix_update_ind);
        
        % Recalculate the exit rates and destination probability based on the new rates.
        Tips_params.Exit_rates = Tips_params.Calculate_exit_rates(Tips_params);
        Tips_params.Transition_dest_prob = Tips_params.Calculate_transition_dest_prob(Tips_params);
    end
    
    % Update the post-collision transition matrix, if it is not constant in time.
    is_transition_matrix_post_coll_updated = strcmp(Tips_params.TipDynamicsModel,'3state') ...
                                            && post_coll_dynamics ...
                                            && ~Is_transition_rate_matrix_post_coll_constant ...
                                            && t == Transition_matrix_post_coll_update_time;
    if is_transition_matrix_post_coll_updated
        Transition_matrix_post_coll_update_ind = Transition_matrix_post_coll_update_ind + 1;
        Tips_params_post_coll.transition_prob_matrix = Tips_params_post_coll_input.transition_rate_matrix(:,:,Transition_matrix_post_coll_update_ind) * timescale;
        Transition_matrix_post_coll_update_time = Transition_matrix_post_coll_update_times(Transition_matrix_post_coll_update_ind);
        
        % Recalculate the exit rates and destination probability based on the new rates.
        Tips_params_post_coll.Exit_rates = Tips_params_post_coll.Calculate_exit_rates(Tips_params_post_coll);
        Tips_params_post_coll.Transition_dest_prob = Tips_params_post_coll.Calculate_transition_dest_prob(Tips_params_post_coll);
    end
    
    % Update the velocity distribution.
    is_V_dist_updated = strcmp(Tips_params.TipDynamicsModel,'3state') ...
                        && ~is_V_dist_constant && t == V_dist_update_time;
    if is_V_dist_updated
        Tips_params.V_dist_fit = V_dist_fits{V_dist_update_ind};
        V_dist_update_ind = V_dist_update_ind + 1;
        V_dist_update_time = V_dist_update_times(V_dist_update_ind);
    end
    
    % Update the post-collision velocity distribution.
    is_V_dist_post_coll_updated = strcmp(Tips_params.TipDynamicsModel,'3state')...
                                  && post_coll_dynamics ...
                                  && ~constant_V_dist_post_coll ...
                                  && t == V_dist_post_coll_update_time;
    if is_V_dist_post_coll_updated
        Tips_params_post_coll.V_dist_fit = V_dist_post_coll_fits{V_dist_post_coll_update_ind};
        V_dist_post_coll_update_ind = V_dist_post_coll_update_ind + 1;
        V_dist_post_coll_update_time = V_dist_post_coll_update_times(V_dist_post_coll_update_ind);
    end
    
    % Update the drift velocity.
    is_V_drift_updated = strcmp(Tips_params.TipDynamicsModel,'drift_diff') && t == V_drift_update_time;
    if is_V_drift_updated
        Tips_params.Drift_velocity = Tips_params_in.Drift_velocity(V_drift_update_ind)/lengthscale*timescale;
        V_drift_update_ind = V_drift_update_ind + 1;
        V_drift_update_time = V_drift_update_times(V_drift_update_ind);
    end
    
    % Update the diffusion coefficient.
    is_D_coeff_updated = strcmp(Tips_params.TipDynamicsModel,'drift_diff') && t == D_coeff_next_update_time;
    if is_D_coeff_updated
        Tips_params.Diffusion_coefficient = Tips_params_in.Diffusion_coefficient(D_coeff_next_update_ind)/lengthscale^2*timescale;
        Tips_params.Diffusion_displacement_sigma = sqrt(2*Tips_params.Diffusion_coefficient);
        D_coeff_next_update_ind = D_coeff_next_update_ind + 1;
        D_coeff_next_update_time = D_coeff_update_times(D_coeff_next_update_ind);
    end
    
    % Update the tips' initial conditions.
    if t == Initial_conditions_time_timeseries(Initial_conditions_update_ind)
        Initial_growth_velocity = Initial_growth_velocity_timeseries(Initial_conditions_update_ind);
        Initial_growth_duration = Initial_growth_duration_timeseries(Initial_conditions_update_ind);
        Initial_conditions_update_ind = Initial_conditions_update_ind + 1;
    end
    %% Build tree structure if it is needed for metric calculations or movies
    % Determine if any events that necessitates the tree structure occurs in this
    % timeframe.
    is_tree_saved = ~isempty(options.Timeseries_tree_times) && Timeseries_Tree_sample_ind <= Timeseries_Tree_N_samples && t == Timeseries_Tree_sampled_time_ind(Timeseries_Tree_sample_ind);
    is_metric_calculated = ~isempty(options.Timeseries_metric) && Timeseries_metric_sample_ind <= Timeseries_metric_N_samples && t == Timeseries_metric_sampled_time_ind(Timeseries_metric_sample_ind);
    is_movieframe_saved = options.Movie && t == Movie_next_frame_ind;
    
    % Build the tree structure.
    is_tree_structure_needed = is_tree_saved || is_metric_calculated || is_movieframe_saved;
    if is_tree_structure_needed
        Tree = restructure_tree(build_tree(), 'RemoveZeroLength', true);
    end
    %% Record timeseries and movie
    % Calculate the radius of gyration (if needed)
    Active_nodes_pos = Nodes_Pos(Nodes_isActive, :);
    if strcmp(options.Branching_rule,'exp_decay_boundary') || options.Timeseries_Rg
        Rg = rg(Active_nodes_pos);
    end
    
    % Record timeseries
    if options.Timeseries
        % Find the lengths of the branches and the number of branch points.
        Lengths = nonzeros(Branches_Length(1:N_branches));
        N_branches = numel(Lengths);
        N_branchpoints = (N_branches - nnz(Branches_Length(1:N_soma_branches) ~= 0))/2;
        Timeseries_Lengths{t} = uint16(Lengths);
        Timeseries_NBranchpoints(t) = uint32(N_branchpoints);

        % Record the radius of gyration.
        if options.Timeseries_Rg
            Timeseries_Rg(t) = Rg;
        end
        Timeseries_D_uniform(t,:) = calculate_tree_size(Active_nodes_pos,'Method','uniform','Rotation',0);
        
        % Count the number of full retractions that occured in this time step.
        if t > 1
            Timeseries_NRetractions(t) = Retractions_counter - Retractions_counter_prev;
        end
        
        % Count the number of collisions.
        if t > 1
            Timeseries_NCollisions(t) = Collisions_counter - Collision_counter_prev;
        end

        % Calculate the convex hull properties.
        if options.ConvexHull
            [CH_area, CH_diam] = convexhull_analysis(Active_nodes_pos, 'Lengthscale', lengthscale);
            Timeseries_CH_Diameter(t) = CH_diam;
            Timeseries_CH_Area(t) = CH_area;
        end
        
        % Save the Tree structure.
        if is_tree_saved
            Timeseries_Tree{Timeseries_Tree_sample_ind} = Tree;
            Timeseries_Tree_sample_ind = Timeseries_Tree_sample_ind + 1;
        end
    end

    % Record movie data.
    if is_movieframe_saved
        % Save the tree structure for movie.
        Movieframes(Movie_frame_counter).tree = Tree;
        Movieframes(Movie_frame_counter).time = t*timescale;

        % Add position of boundary (if any).
        if is_boundary_existent
            Movieframes(Movie_frame_counter).boundary_size = Boundary_size(t,:)*lengthscale;
        end
        
        % Increment the frame counter and determine the next recorded frames.
        Movie_frame_counter = Movie_frame_counter + 1;
        if Movie_frame_counter <= N_movieframes
            Movie_next_frame_ind = Movie_RecordedFramesInd(Movie_frame_counter);
        else
            Movie_next_frame_ind = inf;
        end
    end
    %% Record metrics
    % Calculate the requested metrics at the sampled times.
    if is_metric_calculated
        % Calculate the metrics.
        Metrics = analyze_tree(struct('struct',Tree),Metrics_cell{:});
        
        % Save the metrics.
        for kk = 1:numel(options.Timeseries_metric)
            switch options.Timeseries_metric{kk}
                case 'df'
                    % Save the fractal dimension.
                    Timeseries_df(Timeseries_metric_sample_ind) = Metrics.df;
                    Timeseries_df_rad{Timeseries_metric_sample_ind} = Metrics.rad(:)';
                    Timeseries_df_C{Timeseries_metric_sample_ind} = Metrics.C(:)';
                case 'meshsize'
                    % Save the meshsize.
                    Timeseries_meshsize(Timeseries_metric_sample_ind) = Metrics.Rh;
            end
        end
        Timeseries_metric_sample_ind = Timeseries_metric_sample_ind + 1;
    end
    %% Calculate bin counts and radial histograms
    % Determine if bincounts or radial histograms need updates.
    update_bincounts = mod(t, Bincounts_update_period) == 1;
    update_radial_hist = options.Timeseries && mod(t, Radial_hist_update_period) == 1;

    % Fetch positions of neurons, branchpoints and branchtips.
    if update_bincounts || update_radial_hist
        Neurons_Pos = Nodes_Pos(Nodes_isActive, :);
        Branchpoints_Pos = Nodes_Pos(Nodes_isBranchpoint, :);
        Branchtips_Pos = Tips_Pos(Tips_State > 0, :);
    end

    % Update bincounts
    if update_bincounts
        calculate_bin_counts();
    end

    % Update radial histograms.
    if update_radial_hist
        [Dendrites_counts, Branchpoints_counts, Branchtips_counts] = calculate_radial_hist();
        t_ind = ceil(t/Radial_hist_update_period);
        Timeseries_dendrites_counts_r(:, t_ind) = Dendrites_counts;
        Timeseries_branchpoints_counts_r(:, t_ind) = Branchpoints_counts;
        Timeseries_branchtips_counts_r(:, t_ind) = Branchtips_counts;
    end
    %% Stop time loop
    % Stop growth if maximal time is reached
    if t >= MaxT
        break;
    end
    %% Branching

    % When branching occurs, the branch that contains the chosen branchpoint is split in two
    % parts and two new branches are created.
    % The existing points downstream (towards the tip) of the branch point
    % form the first new branch.
    % The second new branch is initiated at the branch point in a growing state
    % with a branch length of value 0.
    % The points upstream of the branch point remain associated with the
    % separated branch.
    % Only nodes that are not a branchpoint and not within a distance
    % Branchpoint_radius of another branchpoint are considered as 
    % candidates to become branchpoints.
    
    % Track the total number of branchpoints added in this timeframe.
    N_branchpoints_added = 0; 
    
    switch options.Branching_rule
        case 'intensive'
            %% Intensive Branching
            % The number of new branchpoints is determined by sampling a
            % Poisson distribution with a rate defined by the current value
            % of the branching rate (omega).
            if t+1 >= Branching_start_time
                N_new_branchpoints = poissrnd(omega_timeseries(t));
            else
                N_new_branchpoints = 0;
            end
            
            % Recycle branches array if the number of new branchpoints exceeds the
            % available space in the arrays.
            if 2*N_new_branchpoints >= N_branches_max - N_branches
                recycle_branch_arrays();
            end
            
            if N_new_branchpoints > 0
                % Find free nodes.
                Free_nodes_ID = find(Nodes_isActive & ~Nodes_isNextToBP);
                
                for ii = 1:N_new_branchpoints
                    % Stop the growth if there are no more available nodes.
                    N_available_nodes = numel(Free_nodes_ID);
                    if N_available_nodes == 0
                        break;
                    end
                    
                    % Among the free nodes, randomly choose the position
                    % of the next branch point.
                    New_branchpoint_selector = randi(N_available_nodes, 1);
                    New_branchpoint_node_ID = Free_nodes_ID(New_branchpoint_selector);
                    
                    % Add the branchpoint
                    Nodes_ID_to_exclude = add_branchpoint(New_branchpoint_node_ID);
                    
                    % Remove the branchpoint nodes from the free nodes for the next
                    % iteration.
                    Free_nodes_ID = Free_nodes_ID(~builtin('_ismemberhelper', Free_nodes_ID, sort(Nodes_ID_to_exclude)));
                end
            end
        case {'uniform','extensive'}
            %% Extensive Branching
            if t+1 >= Branching_start_time
                % Determine the number of new branchpoints for each non_empty bins.
                Nonempty_bin_linind = find(~Bins_isempty);
                Nonempty_bin_dendrites_counts = Bins_total_counts(Nonempty_bin_linind);
                Bins_N_new_branchpoints = poissrnd(Nonempty_bin_dendrites_counts * omega_timeseries(t));
                N_new_branchpoints = sum(Bins_N_new_branchpoints);
            else
                N_new_branchpoints = 0;
            end
            
            % Recycle branches array if the number of new branchpoints exceeds the
            % available space in the arrays.
            if 2*N_new_branchpoints >= N_branches_max - N_branches
                recycle_branch_arrays();
            end
            
            if N_new_branchpoints > 0
                % Find the bins where branching happens.
                is_bin_branching = Bins_N_new_branchpoints > 0;
                Branching_bins_N_new_branchpoints = Bins_N_new_branchpoints(is_bin_branching);
                N_branching_bins = numel(Branching_bins_N_new_branchpoints);
                Branching_bins_linind = Nonempty_bin_linind(is_bin_branching);
                Branching_bins_subind = [mod(Branching_bins_linind-1, Search_bins_size(1))+1, ceil(Branching_bins_linind/Search_bins_size(2))];
                
                % Handle branching events for each bin.
                for jj = 1:N_branching_bins
                    Bin_subind = Branching_bins_subind(jj, :);
                    
                    % Find all the free nodes in the current bin.
                    % Free nodes are nodes that are not in proximity of branchpoints.
                    Bin_nodes_ID = double(nonzeros(Bins_NodesID(Bin_subind(1), Bin_subind(2), :)));
                    Free_nodes_ID = Bin_nodes_ID(~Nodes_isNextToBP(Bin_nodes_ID));
                    
                    for ii = 1:Branching_bins_N_new_branchpoints(jj)
                        % Stop the growth if there are no more available nodes.
                        N_available_nodes = numel(Free_nodes_ID);
                        if N_available_nodes == 0
                            break;
                        end
                        
                        % Among the free nodes, randomly choose the position
                        % of the next branch point.
                        New_branchpoint_selector = randi(N_available_nodes, 1);
                        New_branchpoint_node_ID = Free_nodes_ID(New_branchpoint_selector);
                        
                        % Add the branchpoint.
                        Nodes_ID_to_exclude = add_branchpoint(New_branchpoint_node_ID);
                        
                        % Remove the branchpoint nodes from the free nodes for the next
                        % iteration.
                        Free_nodes_ID = Free_nodes_ID(~builtin('_ismemberhelper', Free_nodes_ID, sort(Nodes_ID_to_exclude)));
                    end
                end
            end
        case 'extensive_branch'
            %% Extensive branching per branch.
            % For this branching mechanism, a branch has a certain
            % probability of spawning a new branchpoint which is proportional to
            % the extensive branching rate and its length.
            
            % Determine the number of new branchpoints for each branch.
            if t+1 >= Branching_start_time
                Active_branches_ID = find(Branches_isActive);
                Active_branches_Length = Branches_Length(Active_branches_ID);
                Active_branches_N_new_branchpoints = poissrnd(Active_branches_Length * omega_timeseries(t));
                N_new_branchpoints = sum(Active_branches_N_new_branchpoints);
            else
                N_new_branchpoints = 0;
            end
            
            % Recycle branches array if the number of new branchpoints exceeds the
            % available space in the arrays.
            if 2*N_new_branchpoints >= N_branches_max - N_branches
                recycle_branch_arrays();
                
                % Re-evaluate the branches ID, since the recycling changed
                % them.
                Active_branches_ID = find(Branches_isActive);
            end
            
            if N_new_branchpoints > 0
                % Find which branches are spawning new branchpoints.
                is_branch_spawning = Active_branches_N_new_branchpoints > 0;
                Spawning_branches_ID = Active_branches_ID(is_branch_spawning);
                Spawning_branches_Length = Active_branches_Length(is_branch_spawning);
                Spawning_branches_N_new_branchpoints = Active_branches_N_new_branchpoints(is_branch_spawning);
                N_spawning_branches = numel(Spawning_branches_ID);
                
                % Randomize the branching order.
                Branching_order = randperm(N_spawning_branches);
                
                % Handle branching events for each spawning branch.
                for jj = Branching_order
                    % Determine the node index where each new branch will be
                    % spawn. The location of the new branchpoint is uniformly
                    % distributed along the spawning branch.
                    
                    % Find all free nodes in the current branch.
                    % Free nodes are nodes that are not in proximity of branchpoints.
                    Branch_nodes_ID = Branches_NodesID(Spawning_branches_ID(jj), 1:Spawning_branches_Length(jj));
                    Free_nodes_ID = Branch_nodes_ID(~Nodes_isNextToBP(Branch_nodes_ID));
                    
                    for ii = 1:Spawning_branches_N_new_branchpoints(jj)
                        % Stop the growth if there are no more available nodes.
                        N_available_nodes = numel(Free_nodes_ID);
                        if N_available_nodes == 0
                            break;
                        end
                        
                        % Among the free nodes, randomly choose the position
                        % of the next branch point.
                        New_branchpoint_selector = randi(N_available_nodes, 1);
                        New_branchpoint_node_ID = Free_nodes_ID(New_branchpoint_selector);
                        
                        % Add the branchpoint.
                        Nodes_ID_to_exclude = add_branchpoint(New_branchpoint_node_ID);
                        
                        % Remove the branchpoint nodes from the free nodes for the next
                        % iteration.
                        Free_nodes_ID = Free_nodes_ID(~builtin('_ismemberhelper', Free_nodes_ID, sort(Nodes_ID_to_exclude)));
                    end
                end
            end
        case 'age_dependent'
            %% Branching decays exponentially with the age of the branch.
            % Find the nodes that will spawn new branchpoints.
            if t+1 >= Branching_start_time
                Branching_nodes_ID = find(Nodes_isActive & ~Nodes_isNextToBP);
                Branching_nodes_ID = Branching_nodes_ID(is_node_branching_func(t - Nodes_BirthTime(Branching_nodes_ID)));
                N_new_branchpoints = numel(Branching_nodes_ID);
                
                % Randomize the order of the branching.
                Branching_nodes_ID = Branching_nodes_ID(randperm(N_new_branchpoints));
            else
                N_new_branchpoints = 0;
            end
            
            % Recyle branches array if the number of new branchpoints exceeds the
            % available space in the arrays.
            if 2*N_new_branchpoints >= N_branches_max - N_branches
                recycle_branch_arrays();
            end
            
            % Loop through each node and create a new branchpoint.
            for ii = 1:N_new_branchpoints
                New_branchpoint_node_ID = Branching_nodes_ID(ii);
                % Skip branching if the node has become close to a
                % branchpoint, caused by the branching of a previous node.
                if Nodes_isNextToBP(New_branchpoint_node_ID)
                    continue
                end
                
                % Add the branchpoint.
                Nodes_ID_to_exclude = add_branchpoint(New_branchpoint_node_ID);
            end
        case 'exp_decay_boundary'
            %% Branching decays exponentially from the position of the boundary.
            % Determine the number of new branchpoints for each non_empty
            % bins. This number will depend on the dendrites content of the
            % bin and the position of the bin with respect to the boundary of the tree.
            
            if t+1 >= Branching_start_time
                Nonempty_bin_linind = find(~Bins_isempty);
                Nonempty_bin_dendrites_counts = Bins_total_counts(Nonempty_bin_linind);
                Nonempty_bin_radius = Bins_centers_r(Nonempty_bin_linind);
                
                if strcmp(options.Branching_decay_type,'spatial_decay+temporal_capped_decay')
                    % Calculate the total branching rate across the whole tree.
                    omega_tot = sum(omega_timeseries(t)*Nonempty_bin_dendrites_counts(:));
                    
                    % For each bin, calculate the distance to the radius of
                    % gyration.
                    Nonempty_bin_Rg_dist = Rg - min(Nonempty_bin_radius,Rg);
                    
                    % Calculate the fraction of the total branching rate that goes into each bin.
                    Branching_rate_weight = Nonempty_bin_dendrites_counts.*branching_rate_func(Nonempty_bin_Rg_dist,t);
                    Branching_rate_weight = Branching_rate_weight/sum(Branching_rate_weight(:));
                    Nonempty_bin_branching_rate = omega_tot * Branching_rate_weight;
                else
                    Nonempty_bin_branching_rate = Nonempty_bin_dendrites_counts.*branching_rate_func(Nonempty_bin_radius,Rg);
                end
                
                % Sample the number of new branchpoints in each bin.
                Bins_N_new_branchpoints = poissrnd(Nonempty_bin_branching_rate);
                N_new_branchpoints = sum(Bins_N_new_branchpoints);
            else
                N_new_branchpoints = 0;
            end
            
            % Recyle branches array if the number of new branchpoints exceeds the
            % available space in the arrays.
            if 2*N_new_branchpoints >= N_branches_max - N_branches
                recycle_branch_arrays();
            end
            
            if N_new_branchpoints > 0
                % Find the bins where branching happens.
                is_bin_branching = Bins_N_new_branchpoints > 0;
                Branching_bins_N_new_branchpoints = Bins_N_new_branchpoints(is_bin_branching);
                N_branching_bins = numel(Branching_bins_N_new_branchpoints);
                Branching_bins_linind = Nonempty_bin_linind(is_bin_branching);
                Branching_bins_subind = [mod(Branching_bins_linind-1, Search_bins_size(1))+1, ceil(Branching_bins_linind/Search_bins_size(2))];
                
                % Handle branching events for each bin.
                for jj = 1:N_branching_bins
                    Bin_subind = Branching_bins_subind(jj, :);
                    
                    % Find all the free nodes in the current bin.
                    % Free nodes are nodes that are not in proximity of branchpoints.
                    Bin_nodes_ID = double(nonzeros(Bins_NodesID(Bin_subind(1), Bin_subind(2), :)));
                    Free_nodes_ID = Bin_nodes_ID(~Nodes_isNextToBP(Bin_nodes_ID));
                    N_new_branchpoints = Branching_bins_N_new_branchpoints(jj);
                    
                    for ii = 1:N_new_branchpoints
                        % Stop the growth if there are no more available nodes.
                        N_available_nodes = numel(Free_nodes_ID);
                        if N_available_nodes == 0
                            break;
                        end
                        
                        % Among the free nodes, randomly choose the position
                        % of the next branch point.
                        New_branchpoint_selector = randi(N_available_nodes, 1);
                        New_branchpoint_node_ID = Free_nodes_ID(New_branchpoint_selector);
                        
                        % Add the branchpoint
                        Nodes_ID_to_exclude = add_branchpoint(New_branchpoint_node_ID);
                        
                        % Remove the branchpoint nodes from the free nodes for the next
                        % iteration.
                        Free_nodes_ID = Free_nodes_ID(~builtin('_ismemberhelper', Free_nodes_ID, sort(Nodes_ID_to_exclude)));
                    end
                end
            end
        case 'exp_decay_pathlength'
            %% Extensive branching per node, which decays as the node is farther from the soma.
            % Find free nodes.
            if t+1 >= Branching_start_time
                Free_nodes_ID = find(Nodes_isActive & ~Nodes_isNextToBP);
                Free_nodes_path_length = Nodes_PathLength(Free_nodes_ID);
                
                % Determine if each node is spawning a new branchpoint.
                Free_nodes_branching_rate = branching_rate_func(Free_nodes_path_length);
                Free_nodes_branching_prob = 1 - exp(-Free_nodes_branching_rate);
                Free_nodes_isbranching = rand(size(Free_nodes_branching_prob)) < Free_nodes_branching_prob;
                Branching_nodes_ind = find(Free_nodes_ID);
                N_new_branchpoints = numel(Free_nodes_isbranching);
            else
                N_new_branchpoints = 0;
            end
            
            % Recycle branches array if the number of new branchpoints exceeds the
            % available space in the arrays.
            if 2*N_new_branchpoints >= N_branches_max - N_branches
                recycle_branch_arrays();
            end
            
            if N_new_branchpoints > 0
                % Randomize the branching order.
                Branching_nodes_ind = Branching_nodes_ind(randperm(N_new_branchpoints));
                
                for ii = 1:N_new_branchpoints
                    % Skip branching if the node is not free anymore.
                    New_branchpoint_node_ID = Free_nodes_ID(Branching_nodes_ind(ii));
                    if Nodes_isNextToBP(New_branchpoint_node_ID)
                        continue;
                    end
                    
                    % Add the branchpoint
                    Nodes_ID_to_exclude = add_branchpoint(New_branchpoint_node_ID);
                end
            end
    end
    %% Miscellaneous
    % Record the number of new branchpoints added.
    if options.Timeseries && t+1 >= Branching_start_time
        Timeseries_NNewBranchpoints(t+1) = N_branchpoints_added;
    end
    
    % Update counters of the previous timestep.
    Collision_counter_prev = Collisions_counter;
    Retractions_counter_prev = Retractions_counter;
    
    % Update waitbar.
    if options.waitbar
        waitbar(t/MaxT, wb);
    end
end
%% Format output and exit
% Record the total simulation time.
Sim_params.TotalSimulationTime = toc;

% Close the waitbar.
if options.waitbar
    close(wb);
end

% Remove branches that have been initialized, but that haven't grown yet.
Dynamic_branches_not_grown = Branches_ID(Branches_isDynamic & Branches_Length == 0);
deletebranch(Dynamic_branches_not_grown);

% Build the tree structure from the branches arrays.
Tree = build_tree();

% Clean the tree structure by removing all non-necessary fields. Also,
% change the format of the recorded data to reduce memory needs.
if nnz([Tree.Length])==0
    % If the entire tree has disappeared, return single branch with the
    % initial point to preserve fieldnames.
    warning('The tree is empty.');
    Treeout = Tree(1); 
    Treeout.Length = 0;
    Treeout.PointsPos = Soma_pos;
else
    Treeout = restructure_tree(Tree,'RemoveZeroLength',true);
end

% Build the timeseries structure.
Timeseries = struct();
if options.Timeseries
    Timeseries.Lengths = Timeseries_Lengths(:);
    Timeseries.Statecounts = Timeseries_Statecounts';
    Timeseries.NCollisions = Timeseries_NCollisions(:);
    Timeseries.NRetractions = Timeseries_NRetractions(:);
    Timeseries.NDeaths = Timeseries_NDeaths(:);
    Timeseries.NBranchpoints = Timeseries_NBranchpoints(:);
    Timeseries.NNewBranchpoints = Timeseries_NNewBranchpoints(:);
    Timeseries.Rg = Timeseries_Rg(:);
    Timeseries.D_uniform = Timeseries_D_uniform;
    Timeseries.CH_Diameter = Timeseries_CH_Diameter(:);
    Timeseries.CH_Area = Timeseries_CH_Area(:);
    Timeseries.dendrites_counts = Timeseries_dendrites_counts_r';
    Timeseries.branchpoints_counts = Timeseries_branchpoints_counts_r';
    Timeseries.branchtips_counts = Timeseries_branchtips_counts_r';
    Timeseries.Radial_hist_binsize = Radial_hist_binsize;
    Timeseries.Radial_hist_update_period = Radial_hist_update_period;
end

if ~isempty(options.Timeseries_metric)
    Timeseries.Metric_times = Timeseries_metric_sampled_time_ind(:)*timescale;
    Timeseries.Metric_times_units = 'min';
    if ismember('df',options.Timeseries_metric)
        Timeseries.Metric_df = Timeseries_df;
        Timeseries.Metric_df_rad = cell2mat(Timeseries_df_rad);
        Timeseries.Metric_df_C = cell2mat(Timeseries_df_C);
    end
    if ismember('meshsize',options.Timeseries_metric)
        Timeseries.Metric_meshsize = Timeseries_meshsize(:);
    end
end

% Save the tree structure timeseries
if ~isempty(options.Timeseries_tree_times)
    Timeseries.Tree_times = Timeseries_Tree_sampled_time_ind(:)*timescale;
    Timeseries.Tree_times_units = 'min';
    Timeseries.Tree = Timeseries_Tree(:);
end

% Save the dynamics tracks' timeseries
if ~isempty(options.Timeseries_tracks)
    N_tracks = double(Timeseries_tracks.Tracks_ID.Count);
    Timeseries_tracks.N_tracks = N_tracks;
    Tracks_arrays_name = {'Tracks_t','Tracks_L','Tracks_length','Tracks_hascollided'};
    for j=1:numel(Tracks_arrays_name)
        Array_fieldname = Tracks_arrays_name{j};
        
        % Define the scale to use to rescale cell arrays into dimensionfull
        % values.
        switch Array_fieldname
            case 'Tracks_t'
                Field_scale = timescale;
            case 'Tracks_L'
                Field_scale = lengthscale;
            case {'Tracks_length','Tracks_hascollided'}
                Field_scale = 1;
            otherwise
                error('Fieldname %s not supported.',Array_fieldname);
        end
        
        % Truncate nans in tracks' cell arrays.
        if iscell(Timeseries_tracks.(Array_fieldname))
            for i=1:N_tracks
                % Truncate the tracks to remove trailing nans and scale the
                % values to reinsert dimensions.
                Track_length = Timeseries_tracks.Tracks_length(i);
                Timeseries_tracks.(Array_fieldname){i} = Field_scale*Timeseries_tracks.(Array_fieldname){i}(1:Track_length,:);
            end
        end
        
        % Columnize tracks' array.
        Timeseries_tracks.(Array_fieldname) = Timeseries_tracks.(Array_fieldname)(:);
    end
    Timeseries.Tracks = Timeseries_tracks;
end

% Assign an empty array to timeseries if no fields were defined.
if isempty(fieldnames(Timeseries))
    Timeseries = [];
end

% Shorten collision arrays.
if options.RecordCollisions
    NRecordedCollisions = min(N_recorded_collisions_max,Collisions_counter);
    Collisions.TipID = Collisions.TipID(1:NRecordedCollisions);
    Collisions.Type = Collisions.Type(1:NRecordedCollisions);
    Collisions.Pos = Collisions.Pos(1:NRecordedCollisions, :);
    Collisions.Time = Collisions.Time(1:NRecordedCollisions);
    Collisions.Angle = Collisions.Angle(1:NRecordedCollisions);
    Collisions.Length = Collisions.Length(1:NRecordedCollisions);
end

% Truncate the tip dynamics information
if options.RecordTipDynamics
    for i = 1:N_recorded_tips
        Recorded_tips_states{i} = Recorded_tips_states{i}(1:Recorded_tips_temporal_ind(i)-1);
        Recorded_tips_NSteps{i} = Recorded_tips_NSteps{i}(1:Recorded_tips_temporal_ind(i)-1);
    end
    Recorded_tips_states = Recorded_tips_states(1:N_recorded_tips);
    Recorded_tips_NSteps = Recorded_tips_NSteps(1:N_recorded_tips);
end

Tips_dynamics = struct();
Tips_dynamics.States = Recorded_tips_states;
Tips_dynamics.NSteps = Recorded_tips_NSteps;
Tips_dynamics.BirthTime = Tips_BirthTime(1:N_tips);
Tips_dynamics.CollisionTime = Tips_CollisionTime(1:N_tips);
Tips_dynamics.DeathTime = Tips_DeathTime(1:N_tips);
Tips_dynamics.BirthPos = Tips_BirthPos(1:N_tips, :);
Tips_dynamics.CollisionPos = Tips_CollisionPos(1:N_tips, :);
Tips_dynamics.DeathPos = Tips_DeathPos(1:N_tips, :);
Tips_dynamics.CollisionAngle = Tips_CollisionAngle(1:N_tips);
Tips_dynamics.NRetractions = Tips_NRetractions(1:N_tips);
Tips_dynamics.NCollisions = Tips_NCollisions(1:N_tips);

% Save last profile.
if options.Profile
    save_profile(Profile_filename);
    profile off;
end

% Clear global variables
clear global
%% Function to find the bin index of position x.
    function Bin_ind = pos2binind(x)
        % Input
        % x = Nx2 array representing positions whose bin index will be
        %     calcualted.
        %
        % Output
        % Bin_ind = Nx2 array representing the 2D bin index of each input
        %           positions
        
        % Calculate positions relative to the soma.
        x = x - Soma_pos;
        Bin_ind = ceil(([-x(:, 2), x(:, 1)] + Last_bin_edge_pos)/binwidth);
    end
%% Function that changes the state of a tip.
    function change_tip_state(Tip_ID, New_state)
        % Input
        % Tip_ID = ID of the tip whose state will changed.
        % New_state (optional) = new state of the tip. The new state can be
        %                        either an integer between 1 and N_states or
        %                        a string. The mapping between string and
        %                        integer is defined by State_char2int().
        
        % Map NewState string to integer if a char is given.
        if nargin == 2 && ischar(New_state)
            New_state = State_char2int(New_state);
        end
        
        switch Tips_params.TipDynamicsModel
            case '3state'
                % Determine which tips parameters (transition matrix and velocity dist) will be used.
                % If the last tip's collision occurred within a span of
                % alpha + Tips_CollisionDuration(TipID), use the
                % post-collision tip parameters.
                if post_coll_dynamics && t < Tips_CollisionTime(Tip_ID) + alpha + Tips_CollisionDuration(Tip_ID)
                    Tips_params_temp = Tips_params_post_coll;
                else
                    Tips_params_temp = Tips_params;
                end
                
                % Determine the next state if not given as input.
                Tip_current_state = Tips_State(TipID);
                if nargin < 2
                    if random_sampler('rand') < Tips_params_temp.Transition_dest_prob(Tip_current_state,1)
                        New_state = Tips_params_temp.Transition_states(Tip_current_state,1);
                    else
                        New_state = Tips_params_temp.Transition_states(Tip_current_state,2);
                    end
                end
                
                
                % Determine the duration of the state by sampling an exponential distribution with the state exit rate.
                % Use the duration to determine the time of the next transitions.
                State_duration = -log(random_sampler('rand'))/Tips_params_temp.Exit_rates(New_state);
                New_transition_time = t + State_duration;
                
                % Determine the velocity in the new state.
                Randn_num = random_sampler('randn');
                switch New_state
                    case 1
                        % Prepare the velocity samples in each state.
                        % Shrinking (State = 1)
                        New_velocity = -exp(Randn_num*Tips_params_temp.V_dist_fit.sigm + Tips_params_temp.V_dist_fit.mum);
                        
                    case 2
                        % Paused (State = 2)
                        New_velocity = Randn_num*Tips_params_temp.V_dist_fit.sig0 + Tips_params_temp.V_dist_fit.mu0;
                        
                    case 3
                        % Growing (State = 3)
                        New_velocity = exp(Randn_num*Tips_params_temp.V_dist_fit.sigp + Tips_params_temp.V_dist_fit.mup);
                end
                
            case 'drift_diff'
                % For this model, the velocity is the same in each of the
                % three state.
                if nargin == 1
                    New_state = 3;
                end
                New_velocity = Tips_params.Drift_velocity;
                New_transition_time = inf;
        end
        
        % Assign the new velocity and state.
        Tips_TransitionTime(Tip_ID) = New_transition_time;
        Tips_State(Tip_ID) = New_state;
        Tips_Velocity(Tip_ID) = New_velocity;
    end
%% Function to search for nodes surrounding a given position.
    function Neigh_nodes_ID = search_nodes(Pos, Search_radius)
        % Input
        % Pos = 1x2 array representing the position around which the search
        %       is performed.
        % Search_radius = radius of the search.
        %
        % Output
        % Neigh_nodes_ID = ID of the nodes that are within a distance of
        %                  Search_radius around Pos.
        
        % Find the bin index of position x.
        Bin_sub_ind = pos2binind(Pos);
        
        % Find the neigboring bins.
        Bin_search_radius = ceil(Search_radius/binwidth);
        Bin_search_range = -Bin_search_radius:Bin_search_radius;
        Row_ind_range = Bin_sub_ind(1) + Bin_search_range;
        Col_ind_range = Bin_sub_ind(2) + Bin_search_range;
        
        % If periodic boundary conditions are used, add bins from the other side,
        % if the position is close to the boundary.
        if options.Boundary_type == 1
            Bin_search_positions_original = [Pos(1) + Bin_search_range'*binwidth, Pos(2) + Bin_search_range'*binwidth];
            Bin_search_positions_translated_up = Bin_search_positions_original + Boundary_size_temp;
            Bin_search_positions_translated_down = Bin_search_positions_original - Boundary_size_temp;
            Bin_ind_translated_up = pos2binind(Bin_search_positions_translated_up);
            Bin_ind_translated_down = pos2binind(Bin_search_positions_translated_down);
            Row_ind_range = [Row_ind_range,Bin_ind_translated_up(:,1)',Bin_ind_translated_down(:,1)'];
            Col_ind_range = [Col_ind_range,Bin_ind_translated_up(:,2)',Bin_ind_translated_down(:,2)'];
        end
        
        % Find the neigboring nodes IDs.
        Neigh_nodes_ID = Bins_NodesID(Row_ind_range, Col_ind_range, :);
        Neigh_nodes_ID = Neigh_nodes_ID(logical(Neigh_nodes_ID));
    end
%% Function to check if a node at a given position collides with neighboring nodes.
    function [Collision, Collided_nodes_ID, Collided_nodes_time] = check_nodes_collision(Trajectory, Excluded_nodes_ID)
        % Input
        % Trajectory = 2x2 array representing the trajector that needs to
        %              be checked for collision. Trajectory(1,:) represents
        %              the starting position.
        % Excluded_nodes_ID = ID of nodes that will be excluded for the
        %                     collision check.
        %
        % Output
        % Collision = bool representing if a collision occurred.
        % Collided_nodes_ID = ID of the nodes that collided with the
        %                     trajectory.
        % Collided_nodes_time = time when the collision of each node
        %                       occurred. The time is given between 0 and 1.
        
        % Define the start and end position.
        Pos_start = Trajectory(1, :);
        Pos_end = Trajectory(2, :);

        % Find neigboring nodes IDs.
        Neigh_nodes_ID = search_nodes(Pos_end, Contact_dist);

        % Remove excluded nodes.
        if ~isempty(Excluded_nodes_ID)
            if ~issorted(Excluded_nodes_ID)
                Excluded_nodes_ID = sort(Excluded_nodes_ID);
            end
            Neigh_nodes_isexcluded = builtin('_ismemberhelper', Neigh_nodes_ID, Excluded_nodes_ID);
            Neigh_nodes_ID = Neigh_nodes_ID(~Neigh_nodes_isexcluded);
        end

        % Initialize output.
        Collision = false;
        Collided_nodes_ID = [];
        Collided_nodes_time = [];

        % If there are no neigboring nodes, terminate early.
        if isempty(Neigh_nodes_ID)
            return;
        end

        % Check collisions. If there is a collision, recalculate the collision
        % time.
        Neigh_nodes_pos = Nodes_Pos(Neigh_nodes_ID, :);
        Neigh_nodes_iscolliding = check_distance(Pos_end, Neigh_nodes_pos, Contact_dist);
        Collision = any(Neigh_nodes_iscolliding);
        if Collision
            % If collision, calculate the collision time by assuming a linear
            % trajectory between the start and end position.

            % Center the end position and the neighboring nodes position
            % with respect to the start position.
            Pos_end_centered = Pos_end - Pos_start;
            Neigh_nodes_pos_centered = Neigh_nodes_pos(Neigh_nodes_iscolliding, :) - Pos_start;
            
            % If periodic boundary conditions are used, use the shortest
            % distance between the clockwise and anti-clockwise direction.
            if options.Boundary_type == 1
                Neigh_nodes_pos_centered = mod(Neigh_nodes_pos_centered + Boundary_size_temp/2,Boundary_size_temp) - Boundary_size_temp/2;
            end
            
            % Calculate the projection of the nodes position vectors
            % along the line that interpolates between the start and end positions.
            Pos_end_norm_squared = sum(Pos_end_centered.^2);
            Neigh_nodes_proj_squared = (Neigh_nodes_pos_centered*Pos_end_centered').^2/(Pos_end_norm_squared);

            % Calculate the perpendicular distance of the nodes with
            % respect to the interpolating line.
            Neigh_nodes_dist_squared = sum(Neigh_nodes_pos_centered.^2, 2);
            Neigh_nodes_dist_perp_squared = Neigh_nodes_dist_squared - Neigh_nodes_proj_squared;

            % Find the nodes that intersect with the interpolating
            % line.
            Neigh_nodes_is_intersecting = Neigh_nodes_dist_perp_squared < Contact_dist^2;

            % Find the intersection time along the interpolating line.
            Intersect_pos = Neigh_nodes_pos_centered(Neigh_nodes_is_intersecting, :);

            Neigh_nodes_ID = Neigh_nodes_ID(Neigh_nodes_iscolliding);
            Collided_nodes_ID = Neigh_nodes_ID(Neigh_nodes_is_intersecting);
            Collided_nodes_time = (Intersect_pos*Pos_end_centered'- ...
                sqrt(Contact_dist^2*Pos_end_norm_squared-(Intersect_pos(:, 1)*Pos_end_centered(2)-Pos_end_centered(1)*Intersect_pos(:, 2)).^2))/Pos_end_norm_squared;
            
            % Set the small collision times below 1e-10 to zero.
            Collided_nodes_time(abs(Collided_nodes_time) <= 1e-10) = 0;
            
%             % Check that at least some collision times are positive.
%             if all(Collided_nodes_time < 0)
%                 figure; hold on; axis equal;
%                 viscircles(Neigh_nodes_pos(Neigh_nodes_iscolliding,:),Contact_dist/2*ones(nnz(Neigh_nodes_iscolliding),1),'Color','Red');
%                 viscircles(Neigh_nodes_pos(~Neigh_nodes_iscolliding,:),Contact_dist/2*ones(nnz(~Neigh_nodes_iscolliding),1),'Color','Green');
%                 viscircles(Pos_end,Contact_dist/2,'Color','Blue');
%                 plot(Pos_start(1),Pos_start(2),'o','DisplayName','Start Position');
%                 plot(Pos_end(1),Pos_end(2),'bx','DisplayName','End Position');
%                 legend;
%                 error('The collision time was not calculated accurately.')
%             end
        end
    end
%% Function to add a node to the tree.
    function Node_ID = add_node(Position, BranchID, Parent_node_ID)
        % Input
        % Position = 1x2 array representing position of the new node.
        % BranchID = ID of the branch to which the node belongs.
        % Parent_node_ID = ID of the parent node to which the new node is
        %                  attached.
        % 
        % Output
        % Node_ID = ID of the new node that is added.
        
        % Check that there is enough space in the nodes arrays. Otherwise, expand
        % them.
        N_freenodes = numel(Nodes_freeindices);
        if Nodes_counter > N_freenodes
            % Find rows in the rows arrays that can be reused. These rows
            % corresponds to Nodes that are not active anymore.
            Nodes_freeindices = find(~Nodes_isActive);

            % Extend the Nodes arrays if there are no more free indices.
            if isempty(Nodes_freeindices)
                N_nodes = numel(Nodes_isActive);
                N_new_nodes = max(1000, round(0.1*N_nodes));
                N_nodes_new = N_nodes + N_new_nodes;
                
                % Check that the number of nodes is not higher than the maximum.
                if N_nodes_new > N_nodes_max
                    error('The number of nodes (%d) has exceeded the maximum (%d).', N_nodes_new, N_nodes_max)
                end

                Nodes_Pos = [Nodes_Pos; zeros(N_new_nodes, 2, 'like', Nodes_Pos)];
                Nodes_BirthTime = [Nodes_BirthTime; zeros(N_new_nodes, 1, 'like', Nodes_BirthTime)];
                Nodes_isActive = [Nodes_isActive; zeros(N_new_nodes, 1, 'like', Nodes_isActive)];
                Nodes_isBranchpoint = [Nodes_isBranchpoint; zeros(N_new_nodes, 1, 'like', Nodes_isBranchpoint)];
                Nodes_isNextToBP = [Nodes_isNextToBP; zeros(N_new_nodes, 1, 'like', Nodes_isNextToBP)];
                Nodes_ParentID = [Nodes_ParentID; zeros(N_new_nodes, 1, 'like', Nodes_ParentID)];
                Nodes_PathLength = [Nodes_PathLength; zeros(N_new_nodes, 1, 'like', Nodes_PathLength)];

                Nodes_freeindices = (N_nodes_new - N_new_nodes + 1):N_nodes_new;
            end
            
            % Reset the counter for the free indices.
            Nodes_counter = 1;
        end
        
        % Find a new node ID.
        Node_ID = Nodes_freeindices(Nodes_counter);
        Nodes_counter = Nodes_counter + 1;

        % Change properties of the new node.
        Nodes_Pos(Node_ID, :) = Position;
        Nodes_BirthTime(Node_ID) = t;
        Nodes_isActive(Node_ID) = true;
        Nodes_isBranchpoint(Node_ID) = false;
        Nodes_BranchID(Node_ID) = BranchID;
        Nodes_ParentID(Node_ID) = Parent_node_ID;
        
        % Calculate the path length to the soma.
        if Parent_node_ID > 0
            Nodes_PathLength(Node_ID) = Nodes_PathLength(Parent_node_ID) + 1;
        end
        
        % Get the bin ID of the new node position.
        Bin_sub_ind = pos2binind(Position);
        
        % Update the bins arrays with the new node.
        Bin_bucket_ind = find(~Bins_NodesID(Bin_sub_ind(1), Bin_sub_ind(2), :), 1);
        if isempty(Bin_bucket_ind)
            error('Could not find an empty entry to store the new node');
        end
        Bins_NodesID(Bin_sub_ind(1), Bin_sub_ind(2), Bin_bucket_ind) = Node_ID;
        %Bins_isempty(Bin_sub_ind(1),Bin_sub_ind(2),Bin_bucket_ind) = false;
        Nodes_BinInd(Node_ID, :) = [Bin_sub_ind, Bin_bucket_ind];
        
        % Check if the new node is next to the soma.
        %is_next_to_BP = check_distance(Soma_pos,Position,Soma_radius);
        is_next_to_BP = sum((Position - Soma_pos).^2,2) < Soma_radius_squared;
        
        % Determine if the new node is next to a branch point.
        if Branchpoint_radius > 0
            % Get neigboring nodes IDs.
            Neigh_nodes_ID = search_nodes(Position, Branchpoint_radius);

            % Retain nodes that are branchpoints.
            Neigh_nodes_isbranchpoint = Nodes_isBranchpoint(Neigh_nodes_ID);
            Branchpoint_nodes_ID = Neigh_nodes_ID(Neigh_nodes_isbranchpoint);

            % Check if the new node is within a distance Branchpoint_radius from a branchpoint.
            if ~isempty(Branchpoint_nodes_ID)
                Branchpoint_nodes_Pos = Nodes_Pos(Branchpoint_nodes_ID, :);
                Branchpoints_isclose = check_distance(Position, Branchpoint_nodes_Pos, Branchpoint_radius);
                is_next_to_BP = any(Branchpoints_isclose);
            end
        end
        Nodes_isNextToBP(Node_ID) = is_next_to_BP;
    end
%% Function to delete nodes from the tree.
    function delete_nodes(Nodes_ID)
        % Input
        % Nodes_ID = ID of nodes that are deleted.
        for Node_ID = Nodes_ID
            % Send error if a branchpoint will be removed. This cannot
            % happen.
            if Nodes_isBranchpoint(Node_ID)
                error('A branchpoint cannot be removed');
            end

            % Deactivate the node.
            Nodes_isActive(Node_ID) = false;

            % Remove the node from the bins.
            Bin_ind = Nodes_BinInd(Node_ID, :);
            Bins_NodesID(Bin_ind(1), Bin_ind(2), Bin_ind(3)) = 0;
        end
    end
%% Function to add a branchpoint to the tree.
    function Nodes_next_to_BP_ID = add_branchpoint(Node_ID)
        % Input
        % Node_ID = ID of the node that becomes a branchpoint.
        
        % Increment the branches and tip counter
        N_branches = N_branches + 2;
        N_tips = N_tips + 1;
        if N_branches > N_branches_max
            error('The maximal number (%d) of branches has been reached.', N_branches_max);
        end
        if N_tips > N_tips_max
            error('The maximal number of tips (%d) has been reached.', N_tips_max);
        end
        
        % Define the ID of the new branches and tip.
        New_branch_1_ID = N_branches-1;
        New_branch_2_ID = N_branches;
        New_tip_ID = N_tips;
        
        Branchpoint_pos = Nodes_Pos(Node_ID, :);
        Separated_branch_ID = Nodes_BranchID(Node_ID);

        % Find the position of the new branch point along the separated
        % branch.
        Separated_branch_Length = Branches_Length(Separated_branch_ID);
        Separated_branch_NodesID = Branches_NodesID(Separated_branch_ID, 1:Separated_branch_Length+1);
        Separated_branch_BP_location = find(Node_ID == Separated_branch_NodesID, 1);
        if isempty(Separated_branch_BP_location)
            error('error. Branch Point was not found on the separated branch.')
        end
        
        % Calculate the branching angle of the new tip.
        if Separated_branch_BP_location < Separated_branch_Length+1
            Adjacent_nodes_index = [Separated_branch_BP_location-1, Separated_branch_BP_location+1];
            Adjacent_nodes_ID = Separated_branch_NodesID(Adjacent_nodes_index);
            Adjacent_nodes_pos = Nodes_Pos(Adjacent_nodes_ID, :);
        else
            % In this case, the branchpoint is at the separated branch tip.
            % To make sure that no intersections happen with future
            % growth of the separated branch, an extra node is
            % positioned along the projection of the parent segment and used
            % to calculate the branching angle.
            Parent_node_ind = Separated_branch_BP_location-1;
            Parent_node_ID = Separated_branch_NodesID(Parent_node_ind);
            Parent_node_pos = Nodes_Pos(Parent_node_ID, :);
            Extra_node_pos = Branchpoint_pos+(Branchpoint_pos-Parent_node_pos);
            Adjacent_nodes_pos = [Parent_node_pos; Extra_node_pos];
        end
        Branching_angle = sample_branching_angle(Branchpoint_pos, Adjacent_nodes_pos, theta_min);
        
        % Update the Tips arrays with the new tip.
        Tips_Theta(New_tip_ID) = Branching_angle;
        Tips_Pos(New_tip_ID, :) = Branchpoint_pos;
        Tips_BirthPos(New_tip_ID, :) = Branchpoint_pos;
        Tips_NodeID(New_tip_ID) = Node_ID;
        Tips_BirthTime(New_tip_ID) = t + 1;
        
        % Initialize the tip in the growing state.
        change_tip_state(New_tip_ID,3);
        
        % Define the initial tip velocity.
        % By default, the initial tip velocity is randomly sampled from the
        % growth state velocity distribution. If a constant
        % Initial_growth_velocity is used, overwrite this sample.
        if ~strcmp(Initial_growth_velocity,'Random') && isnumeric(Initial_growth_velocity)
            Tips_Velocity(New_tip_ID) = Initial_growth_velocity;
        end
        
        % Define the duration of the initial condition.
        % By default, the duration of the initial condition is randomly 
        % sampled from the growth state duration distribution, which is an
        % exponential distribution with parameter defined by the exit rate.
        % If a constant Initial_growth_duration is used, overwrite this sample.
        if ~strcmp(Initial_growth_duration,'Random') && isnumeric(Initial_growth_duration)
            Tips_TransitionTime(New_tip_ID) = Tips_BirthTime(New_tip_ID) + Initial_growth_duration;
        end
        
        % Determine if the new tip's dynamics will be recorded.
        if options.RecordTipDynamics && N_recorded_tips <= t/MaxT*N_recorded_tips_max
            N_recorded_tips = N_recorded_tips + 1;
            Tips_isrecorded(New_tip_ID) = true;
            Tips_record_Ind(New_tip_ID) = N_recorded_tips;
        end
        
        % Update the branchpoint proximity status of the nodes due the new
        % branch point.
        if Branchpoint_radius > 0
            % Find the nodes that are a distance Branchpoint_radius from the branchpoint.
            Neigh_nodes_ID = search_nodes(Branchpoint_pos, Branchpoint_radius);
            Neigh_nodes_pos = Nodes_Pos(Neigh_nodes_ID, :);
            Neigh_nodes_is_within_BP = check_distance(Branchpoint_pos, Neigh_nodes_pos, Branchpoint_radius);

            Nodes_next_to_BP_ID = Neigh_nodes_ID(Neigh_nodes_is_within_BP);
        else
            % When the branchpoint radius is zero, only the branchpoint itself
            % is considered "near" a branchpoint. Since branching candidates 
            % are defined by Nodes_isNextToBP, this is necessary to
            % choosing branchpoints as candidates for new branchpoints.
            Nodes_next_to_BP_ID = Node_ID;
        end
        Nodes_isNextToBP(Nodes_next_to_BP_ID) = true;
        
        % Change the branchID of the nodes associated with the separated branch.
        New_branch_1_nodes_ID = Separated_branch_NodesID(Separated_branch_BP_location:Separated_branch_Length+1);
        Nodes_BranchID(New_branch_1_nodes_ID) = New_branch_1_ID;

        % Change the node to a branchpoint node.
        Nodes_isBranchpoint(Node_ID) = true;

        % Find the length of the first new branch.
        New_branch_1_Length = Separated_branch_Length + 1 - Separated_branch_BP_location;

        % Update the branch structure of the 1st new
        % branch. This branch corresponds to the downstream
        % segment of the separated branch.
        Branches_ParentID(New_branch_1_ID) = Separated_branch_ID;
        Branches_SiblingID(New_branch_1_ID) = New_branch_2_ID;
        Branches_ChildrenID(New_branch_1_ID, :) = Branches_ChildrenID(Separated_branch_ID, :);
        Branches_Length(New_branch_1_ID) = New_branch_1_Length;
        Branches_NodesID(New_branch_1_ID, 1:New_branch_1_Length+1) = Branches_NodesID(Separated_branch_ID, Separated_branch_BP_location:Separated_branch_Length+1);
        Branches_TipID(New_branch_1_ID) = Branches_TipID(Separated_branch_ID);
        Branches_isDynamic(New_branch_1_ID) = Branches_isDynamic(Separated_branch_ID);
        Branches_isActive(New_branch_1_ID) = true;

        % Update the branch structure of the 2nd new
        % branch. This branch corresponds to the additional
        % segment with zero length.
        Branches_ParentID(New_branch_2_ID) = Separated_branch_ID;
        Branches_SiblingID(New_branch_2_ID) = New_branch_1_ID;
        Branches_Length(New_branch_2_ID) = 0;
        Branches_NodesID(New_branch_2_ID, 1) = Node_ID;
        Branches_TipID(New_branch_2_ID) = New_tip_ID;
        Branches_BranchingAngle(New_branch_2_ID) = Branching_angle;
        Branches_isDynamic(New_branch_2_ID) = true;
        Branches_isActive(New_branch_2_ID) = true;

        % Update the data of the separated branch's children (if any).
        Children_branch_ID = Branches_ChildrenID(Separated_branch_ID, :);
        if all(Children_branch_ID > 0)
            for ID = Children_branch_ID
                Branches_ParentID(ID) = New_branch_1_ID;
            end
        end

        % Update the data of the separated branch.
        Branches_Length(Separated_branch_ID) = Separated_branch_BP_location-1;
        Branches_ChildrenID(Separated_branch_ID, :) = [New_branch_1_ID, New_branch_2_ID];
        Branches_TipID(Separated_branch_ID) = nan;
        Branches_isDynamic(Separated_branch_ID) = false;
        
        % Increment the number branchpoints counter.
        N_branchpoints_added = N_branchpoints_added + 1;
        
        % Debugging test.
        if Debugmode
            Test('Branching');
        end
    end
%% Function to delete fully-retracted branches.
    function deletebranch(Deleted_branches_ID)
        % Input
        % Deleted_branches_ID = ID of branches that will be deleted.
        
        N_deleted_branches = numel(Deleted_branches_ID);
        for m = 1:N_deleted_branches
            Branch_ID = Deleted_branches_ID(m);
            Parent_ID = Branches_ParentID(Branch_ID);

            % If the branch is at the soma and there are more than 1 branch at the soma
            % simply stop the growth of Branch_ID. Since the branchpoint at the soma is
            % permanent, the bin counts are not changed due to the disappearance of
            % a soma branch.
            if Parent_ID == 0
                Branches_isDynamic(Branch_ID) = false;
                Branches_isActive(Branch_ID) = false;
                continue;
            end

            % Find the node ID of the branch point that will disappear.
            Disappearing_BP_node_ID = Branches_NodesID(Branch_ID, 1);
            Disappearing_BP_node_pos = Nodes_Pos(Disappearing_BP_node_ID, :);

            % Update the nodes arrays.
            
            % Change status of the removed branchpoint to reflect that it
            % has become a normal node.
            Nodes_isBranchpoint(Disappearing_BP_node_ID) = false;

            % Change the proximity-to-branchpoint status of nearby nodes.
            if Branchpoint_radius > 0
                % Find nodes next to the disappearing branchpoint.
                Bin_search_radius = 2*Branchpoint_radius;
                Neigh_node_ID = search_nodes(Disappearing_BP_node_pos, Bin_search_radius);
                Neigh_pos = Nodes_Pos(Neigh_node_ID, :);
                
                % Among the neighboring nodes, find those that are branch points.
                Neigh_isBP = Nodes_isBranchpoint(Neigh_node_ID);
                Neigh_branchpoints_pos = Neigh_pos(Neigh_isBP, :);
                
                % Among the neighbouring nodes, find nodes that are a
                % distance Branchpoint_radius from the disappearing branchpoint.
                Neigh_iswithinBP = check_distance(Disappearing_BP_node_pos, Neigh_pos, Branchpoint_radius);
                
                % Among the neighbouring nodes, find nodes that are close
                % to the soma.
                if Soma_radius > 0
                    Neigh_is_close_to_soma = check_distance(Soma_pos,Neigh_pos,Soma_radius);
                else
                    Neigh_is_close_to_soma = false;
                end
                
                % Find nodes whose proximity to a branchpoint needs to be
                % rechecked.
                Neigh_nodes_needs_check = ~Neigh_isBP & Neigh_iswithinBP & ~Neigh_is_close_to_soma;
                Neigh_nodes_checked_pos_and_ID = [Neigh_pos(Neigh_nodes_needs_check,:), double(Neigh_node_ID(Neigh_nodes_needs_check))];

                % For each node in the neighbourhood, check if they are
                % close to another branchpoint. Keep nodes that are not in
                % proximity to any other branchpoints in the neighbourhood.
                for k = 1:size(Neigh_branchpoints_pos, 1)
                    dist_squared = sum(bsxfun(@minus, Neigh_nodes_checked_pos_and_ID(:, 1:2), Neigh_branchpoints_pos(k, :)).^2, 2);
                    Neigh_nodes_checked_pos_and_ID = Neigh_nodes_checked_pos_and_ID(dist_squared > Branchpoint_radius_squared, :);
                end
                
                % Update the BP proximity status of the nodes that survived
                % all distance checks.
                for k = 1:size(Neigh_nodes_checked_pos_and_ID)
                    Nodes_isNextToBP(Neigh_nodes_checked_pos_and_ID(k, 3)) = false;
                end
            else
                % When the branchpoint radius is zero, only the removed
                % branchpoint proximity status needs to be updated.
                Nodes_isNextToBP(Disappearing_BP_node_ID) = false;
            end

            % Update the branches arrays.

            % Get info on the deleted branch.
            Sibling_ID = Branches_SiblingID(Branch_ID);
            Sibling_Length = Branches_Length(Sibling_ID);
            Sibling_isDynamic = Branches_isDynamic(Sibling_ID);
            Sibling_processed_ID_ind = find(Deleted_branches_ID == Sibling_ID, 1);
            is_Sibling_processed = ~isempty(Sibling_processed_ID_ind);
            Branch_and_sibling_ID = [Branch_ID, Sibling_ID];
            
            if Sibling_Length == 0 && ~Sibling_isDynamic && ~is_Sibling_processed
                error('The Sibling has zero length, but it is not accounted for in the zero-length IDs');
            end

            if ~is_Sibling_processed
                % Add the length and branch points of the sibling to the parent branch.
                ParBranchLength = Branches_Length(Parent_ID);
                Branches_Length(Parent_ID) = ParBranchLength+Sibling_Length;
                Branches_NodesID(Parent_ID, ParBranchLength+1:ParBranchLength+Sibling_Length+1) = Branches_NodesID(Sibling_ID, 1:Sibling_Length+1);
                Branches_ChildrenID(Parent_ID, :) = Branches_ChildrenID(Sibling_ID, :);
                Branches_TipID(Parent_ID) = Branches_TipID(Sibling_ID);
                Branches_isDynamic(Parent_ID) = Branches_isDynamic(Sibling_ID);

                % Update the parent ID of the daughters.
                ChildrenID = Branches_ChildrenID(Sibling_ID, :);
                for k = ChildrenID
                    if k > 0
                        Branches_ParentID(k) = Parent_ID;
                    end
                end

                % Update the branch ID of the sibling nodes.
                Sibling_nodes_ID = Branches_NodesID(Sibling_ID, 1:Sibling_Length+1);
                Nodes_BranchID(Sibling_nodes_ID) = Parent_ID;

                % Set the Length of the sibling branch to zero.
                Branches_Length(Sibling_ID) = 0;

                % Since the branchpoint associated to the sibling branch is disappearing,
                % reset the number of full retractions of the tip connected to the
                % sibling branch if it is active.
                if Sibling_isDynamic
                    Sibling_TipID = Branches_TipID(Sibling_ID);
                    Tips_NRetractions(Sibling_TipID) = 0;
                end

                % Deactivate the state of the current tip.
                TipID = Branches_TipID(Branch_ID);
                Tips_State(TipID) = nan;

            elseif Sibling_processed_ID_ind > m
                % If the sister is also being deleted, check that it hasn't
                % been processed yet. This will be true if it is a member
                % of the upcoming IDs in Deleted_branches_ID.

                % In this case, both branches retracted at the same time.
                % The parent branch starts moving and its tip ID is
                % chosen randomly among its daughter's tip ID.

                % Make the parent branch dynamic.
                Branches_isDynamic(Parent_ID) = true;

                % Among the branch and its sibling, select the candidate tip that
                % will continue moving.
                TipID_candidates = Branches_TipID(Branch_and_sibling_ID);
                Selected_ind = randi(2);
                Selected_TipID = TipID_candidates(Selected_ind);
                Branches_TipID(Parent_ID) = Selected_TipID;
                Tips_NRetractions(Selected_TipID) = 0;
                Tips_DeathTime(Selected_TipID) = nan;

                % Remove the children of the parent.
                Branches_ChildrenID(Parent_ID, :) = [0, 0];

                % Change the nodes branch ID of the branchpoint.
                Nodes_BranchID(Disappearing_BP_node_ID) = Parent_ID;

                % Deactivate the state of the tip that was not selected.
                NotSelected_ind = mod(Selected_ind, 2)+1;
                NotSelected_TipID = TipID_candidates(NotSelected_ind);
                Tips_State(NotSelected_TipID) = nan;
            end

            % Deactivate the branch and its sister.
            Branches_isDynamic(Branch_and_sibling_ID) = false;
            Branches_isActive(Branch_and_sibling_ID) = false;
            Branches_TipID(Branch_and_sibling_ID) = nan;
        end

        % Reset the zero-length branch counter.
        Retracted_branches_counter = 0;
    end
%% Function to calculate the rebranching angle of a set of branches.
    function calculate_rebranching_angle(Branches_ID)
        % Input
        % Branches_ID = ID of branches whose rebranchng angle is recalculated.
        for BranchID = Branches_ID
            % Sample the regrowth angle.
            TipID = Branches_TipID(BranchID);
            First_segment_NodesID = Branches_NodesID(BranchID, 1:2);

            % Calculate the previous branching angle.
            if First_segment_NodesID(2) > 0
                First_segment_vec = diff(Nodes_Pos(First_segment_NodesID, :));
                Branching_angle = atan2(First_segment_vec(2), First_segment_vec(1));
            else
                % In this case, the branch hasn't grown by
                % at least 1 segment, use the attempted angle
                % as the branching angle.
                Branching_angle = Tips_Theta(TipID);
            end

            % Determine the position of the nodes adjacent to
            % the branch point. These nodes will be used to determine
            % rebranching directions that are free of collisions.
            ParentID = Branches_ParentID(BranchID);
            if ParentID == 0
                % Handle the special case of soma branches.
                Siblings_ID = setdiff(Soma_branches_ID,BranchID);
                Siblings_NodeID = Branches_NodesID(Siblings_ID, 2);
                Adjacent_nodes_pos = Nodes_Pos(Siblings_NodeID, :);
            else
                % Find the node of the parent that is adjacent to the
                % branch point.
                Parent_Length = Branches_Length(ParentID);
                Parent_nexttolast_NodeID = Branches_NodesID(ParentID, Parent_Length);
                Parent_nexttolast_NodePos = Nodes_Pos(Parent_nexttolast_NodeID, :);

                Sibling_ID = Branches_SiblingID(BranchID);
                Sibling_Length = Branches_Length(Sibling_ID);
                if Sibling_Length > 0
                    Sibling_NodeID = Branches_NodesID(Sibling_ID, 2);
                    Sibling_NodePos = Nodes_Pos(Sibling_NodeID, :);
                else
                    % In this case, the sibling branch hasn't grown or is about
                    % to grow. Either way, it provides no obstacle to the growth
                    % of the current branch.
                    Sibling_NodePos = [];
                end

                Adjacent_nodes_pos = [Sibling_NodePos; Parent_nexttolast_NodePos];
            end

            % Sample the regrowth angle and save it as the new growth
            % angle of the branch tip.
            Tips_Theta(TipID) = sample_rebranching_angle(Branching_angle, Branchpoint_pos, Adjacent_nodes_pos);
        end

        % Reset the zero-length branch counter.
        Rebranched_branches_counter = 0;
    end
%% Function to calculate the counts of dendrites, branchpoints and tips in each bin.
    function calculate_bin_counts()
        % Count the total amount of neuron particles.
        Neurons_bin_ind = pos2binind(Neurons_Pos);
        Bins_total_counts = accumarray(Neurons_bin_ind, 1, Search_bins_size);

        % Count the number of branchpoints in each bin.
        Branchpoints_bin_ind = pos2binind(Branchpoints_Pos);
        Bins_branchpoint_counts = accumarray(Branchpoints_bin_ind, 1, Search_bins_size);

        % Count the number of tips in each bin.
        Tips_bin_ind = pos2binind(Branchtips_Pos);
        Bins_tip_counts = accumarray(Tips_bin_ind, 1, Search_bins_size);
        
        % Calculate the number of dendrites in each bin.
        Bins_dendrites_counts = Bins_total_counts - Bins_tip_counts - Bins_branchpoint_counts;
        
        % Determine which bin is empty.
        Bins_isempty = Bins_total_counts == 0;
    end
%% Function to calculate the radial histograms of the 3 different species
    function [Dendrites_radial_counts, Branchpoints_radial_counts, Branchtips_radial_counts] = calculate_radial_hist()
        % Each branch node is either a branchtip, a branchpoint or a core
        % node (dendrite node). The squared distance of all nodes
        % from the soma is used to bin the radial density.
        
        % Neurons
        Neurons_dist_squared = sum(Neurons_Pos.^2, 2);
        Neurons_radial_counts = histc(Neurons_dist_squared, Radial_hist_binedges2, 1)';
        Neurons_radial_counts = Neurons_radial_counts(1:end-1);

        % Branchpoints
        Branchpoints_dist_squared = sum(Branchpoints_Pos.^2, 2);
        Branchpoints_radial_counts = histc(Branchpoints_dist_squared, Radial_hist_binedges2, 1)';
        Branchpoints_radial_counts = Branchpoints_radial_counts(1:end-1);

        % Tips
        Branchtips_dist_squared = sum(Branchtips_Pos.^2, 2);
        Branchtips_radial_counts = histc(Branchtips_dist_squared, Radial_hist_binedges2, 1)';
        Branchtips_radial_counts = Branchtips_radial_counts(1:end-1);

        % Dendrites
        Dendrites_radial_counts = Neurons_radial_counts - Branchpoints_radial_counts - Branchtips_radial_counts;
    end
%% Function to sample the branching angle.
    function Branching_angle = sample_branching_angle(Branchpoint_pos, Adjacent_nodes_pos, Theta_min)
        % This function samples a random branching angle of a branch originating
        % from branchpoint_pos. The branch does not take a direction that is 
        % within an angle Theta_min of the two line segments defined by 
        % joining the two adjacent nodes position (adjacent_nodes_pos) with
        % the branchpoint position.
        %
        % Input
        % Branchpoint_pos = 1x2 array reprensenting the position of the branchpoint
        % Adjacent_nodes_pos = Nx2 array representing the positions of the
        %                      adjacent nodes.
        % Theta_min = Minimum angle (rad) between the obstructing line
        %             segments and the branching angle.
        %
        % Output
        % Branching_angle = Branching angle (rad).
        
        % Determine the free angle intervals around the branchpoint.
        Free_angular_intervals = calculate_free_angular_intervals(Branchpoint_pos, Adjacent_nodes_pos, Theta_min);
        
        % Calculate the parent branch angle by performing PCA on the branchpoint
        % and its adjacent nodes position to find the best fitted vector.
        Nodes_position = [Adjacent_nodes_pos(1,:);Branchpoint_pos;Adjacent_nodes_pos(2,:)];
        COM = mean(Nodes_position);
        [U,S,V] = svd(Nodes_position - COM);
        Nodes_projections = U*S;
        Symmetry_vecs = V;
        
        % Find the best-fit vector. The best-fit vector will be the vector with
        % the highest average projection.
        Total_projections = mean(abs(Nodes_projections));
        [~,Branch_vector_ind] = max(Total_projections);
        Bestfit_vector = Symmetry_vecs(:,Branch_vector_ind);
        
        % Flip the vector if the last point has a negative
        % projection (vector is pointing away from the last point).
        if Nodes_projections(end,Branch_vector_ind) < 0
            Bestfit_vector = -Bestfit_vector;
        end
        Parent_angle = atan2(Bestfit_vector(2),Bestfit_vector(1));
        
        % Select a random branching angle from the branching angle
        % distribution until it is in a free interval.
        is_angle_free = false;
        N_iterations_max = 20;
        iteration = 0;
        while ~is_angle_free
            if isempty(Branching_i) || Branching_i > Branching_imax
                Branching_angle_samples = random(options.Branching_angle_dist,Branching_imax,1);
                Branching_angle_samples = (2*double(rand(size(Branching_angle_samples)) < 0.5)-1).*Branching_angle_samples;
                Branching_i = 1;
            end
            Branching_angle = Branching_angle_samples(Branching_i);
            Branching_i = Branching_i+1;
            
            % Define the absolute branching angle by adding the parent angle.
            Branching_angle = mod(Branching_angle + Parent_angle,2*pi);
            
            % Check if the angle is free.
            is_angle_free = any(Branching_angle > Free_angular_intervals(:,1) & Branching_angle < Free_angular_intervals(:,2));
            
            % Stop loop if maximum is reached.
            iteration = iteration + 1;
            if iteration > N_iterations_max
                warning('Maximal number of iterations reached for branching angle sampling.');
                break;
            end
        end
        
        % Plot the line segments for debugging.
        if Debugmode
            figure; axis equal; hold on
            Legend_h = []; Leg_str = [];
            Line_vectors = bsxfun(@minus, Adjacent_nodes_pos, Branchpoint_pos);
            N_lines = size(Line_vectors, 1);
            for k = 1:N_lines
                Line_enpoints = bsxfun(@plus, [0, 0; Line_vectors(k, :)], Branchpoint_pos);
                h_line = plot(Line_enpoints(:, 1), Line_enpoints(:, 2));
                Legend_h = [Legend_h, h_line];
                Leg_str = [Leg_str, {sprintf('Segment %d', k)}];
                
                Excluded_boundary_angle = atan2(Line_vectors(k, 2), Line_vectors(k, 1))+[-1, 1]*Theta_min;
                Boundary_endpoints = bsxfun(@plus, [0, 0; cos(Excluded_boundary_angle(:)), sin(Excluded_boundary_angle(:))], Branchpoint_pos);
                Excluded_h = fill(Boundary_endpoints(:, 1), Boundary_endpoints(:, 2), h_line.Color, 'FaceAlpha', 0.3, 'EdgeColor', 'None');
                Legend_h = [Legend_h, Excluded_h(1)];
                Leg_str = [Leg_str, {sprintf('Segment %d Excluded Region', k)}];
            end
            
            % Plot the new line with the chosen angle.
            New_branch_endpoints = bsxfun(@plus, [0, 0; cos(Branching_angle), sin(Branching_angle)], Branchpoint_pos);
            Branch_h = plot(New_branch_endpoints(:, 1), New_branch_endpoints(:, 2), 'k');
            legend(Legend_h, Leg_str)
        end
    end
%% Function to sample the rebranching angle at a given branchpoint.
    function Branching_angle = sample_rebranching_angle(Branching_angle_prev, Branchpoint_pos, Neighbour_nodes_pos)
        % Input
        % Branching_angle_prev = Angle of previous branching angle. The previous
        %                        angle is used if a free angle cannot be found.
        % Branchpoint_pos = 1x2 array representing position of the branchpoint.
        % Neighbour_nodes_pos = Nx2 array representing position of neighbouring
        %                       nodes that may obstruct the growth of the branch.
        %
        % Output
        % Branching_angle = Angle (rad) of rebranching.
        
        % Define the empirical cdf of the regrowth angle distribution. This
        % is evaluated only once.
        if isempty(Rebranching_i)
            if isempty(Rebranching_angle_dist_fit)
                % If there is no distribution fit, create an empirical cdf
                % that will be used to sample rebranching angles.
                [Rebranching_angle_ecdf, Rebranching_angle_ecdf_x] = ecdf(Rebranching_angle_dist);
            end
            Rebranching_i = Rebranching_imax + 1;
        end

        % Sample a random regrowth angle until it falls within the free
        % angular intervals of the branch point.
        Free_angular_intervals = calculate_free_angular_intervals(Branchpoint_pos, Neighbour_nodes_pos, theta_min);
        is_theta_free = false;
        for iii = 1:N_rebranching_trials_max
            % Create a sample of rebranching angles.
            if Rebranching_i > Rebranching_imax
                if isempty(Rebranching_angle_dist_fit)
                    [~, bin_ind] = histc(rand(Rebranching_imax, 1), Rebranching_angle_ecdf);
                    Rebranching_angle_samples = Rebranching_angle_ecdf_x(bin_ind);
                else
                    Rebranching_angle_samples = random(Rebranching_angle_dist_fit,[Rebranching_imax, 1]);
                end
                Rebranching_i = 1;
            end

            Branching_angle = mod(Rebranching_angle_samples(Rebranching_i) + Branching_angle_prev, 2*pi);
            Rebranching_i = Rebranching_i+1;

            % Test if the rebranching angle is in a free angular region. If
            % so, end the sampling trials.
            is_theta_free = any(Branching_angle >= Free_angular_intervals(:, 1) & Branching_angle < Free_angular_intervals(:, 2));
            if is_theta_free
                break;
            end
        end

        % If a free value was not found, use the initial branching angle.
        if ~is_theta_free
            if options.warnings
                warning('A free branching angle cannot be found. The search stopped and the same branching angle was returned.')
            end
            Branching_angle = Branching_angle_prev;
        end
    end
%% Function to recycle memory by rearranging branch arrays.
    function recycle_branch_arrays()
        % Note that Branches_ID doesn't need to be recycled as it indicates
        % the index of each branch.
        is_branch_active = Branches_isActive;
        Real_branches_ID = Branches_ID(is_branch_active);
        N_real_branches = numel(Real_branches_ID);
        New_ParentID = rep(Branches_ParentID(Real_branches_ID),Real_branches_ID,1:N_real_branches);
        New_SiblingID = rep(Branches_SiblingID(Real_branches_ID),Real_branches_ID,1:N_real_branches);
        New_ChildrenID = rep(Branches_ChildrenID(Real_branches_ID,:),Real_branches_ID,1:N_real_branches);
        
        % Rearrange branches array by moving the real branches info at the beginning of each respective array.
        Branches_NodesID(1:N_real_branches,:) = Branches_NodesID(Real_branches_ID,:);
        Branches_ParentID(1:N_real_branches) = New_ParentID;
        Branches_SiblingID(1:N_real_branches) = New_SiblingID;
        Branches_ChildrenID(1:N_real_branches,:) = New_ChildrenID;
        Branches_Length(1:N_real_branches) = Branches_Length(Real_branches_ID);
        Branches_isDynamic(1:N_real_branches) = Branches_isDynamic(Real_branches_ID);
        Branches_isActive(1:N_real_branches) = Branches_isActive(Real_branches_ID);
        Branches_TipID(1:N_real_branches) = Branches_TipID(Real_branches_ID);
        Branches_BranchingAngle(1:N_real_branches) = Branches_BranchingAngle(Real_branches_ID);

        % Erase the remaining entries using the last entry as a fill value.
        if Branches_Length(end) > 0 || Branches_isDynamic(end)
            error('The fill value is non-zero');
        end
        Branches_NodesID(N_real_branches+1:end,:) = Branches_NodesID(end);
        Branches_ParentID(N_real_branches+1:end) = Branches_ParentID(end);
        Branches_SiblingID(N_real_branches+1:end) = Branches_SiblingID(end);
        Branches_ChildrenID(N_real_branches+1:end,:) = Branches_ChildrenID(end);
        Branches_Length(N_real_branches+1:end) = Branches_Length(end);
        Branches_isDynamic(N_real_branches+1:end) = Branches_isDynamic(end);
        Branches_isActive(N_real_branches+1:end) = Branches_isActive(end);
        Branches_TipID(N_real_branches+1:end) = Branches_TipID(end);
        Branches_BranchingAngle(N_real_branches+1:end) = Branches_BranchingAngle(end);

        % Rearrange nodes array that stores the BranchID.
        Nodes_BranchID = rep(Nodes_BranchID,Real_branches_ID,1:N_real_branches);
        N_branches = N_real_branches;
    end
%% Function that checks if locations y are within a distance d of x.
    function [y_is_within_distance, Closest_y_ind] = check_distance(x, y, d)
        % Input
        % x = 1x2 array representing the 2D position from which distances
        %     are calculated.
        % y = Nx2 array representing locations whose distance from x need
        %     to be checked.
        % d = 1x1 array representing distance used for check.
        %
        % Output
        % y_is_within_distance = Nx1 bool array indicated if y is within a
        %                        distance d from x.
        % Closest_y_ind = Row index of y array indicating which point is
        %                 closest to x.
        
        % Center the targets position at location x.
        y_centered = y;
        y_centered(:, 1) = y(:, 1) - x(1);
        y_centered(:, 2) = y(:, 2) - x(2);
        
        % If periodic boundary conditions are used, use the shortest
        % distance between the clockwise and anti-clockwise direction.
        if options.Boundary_type == 1 && any(abs(y_centered) > Boundary_size_temp/2,[1,2])
            y_centered = mod(y_centered + Boundary_size_temp/2,Boundary_size_temp) - Boundary_size_temp/2;
        end
        
        % Remove targets outside of a rectangle of size d.
        y_is_within_distance = all(abs(y_centered) <= d, 2);
        
        % If any targets remain in the rectangle centered on x, measure exact squared distances
        % to neighboring targets to check for collisions.
        if any(y_is_within_distance)
            % Calculate the distance between x and y.
            y_centered = y_centered(y_is_within_distance, :);
            y_dist_squared = sum(y_centered.^2, 2);
            
            % Determine which targets y is within a distance d of x.
            y_is_inside = y_dist_squared < ((1-100*eps)*d)^2;
            
            % If requested, find the target that is closest to x.
            if nargout > 1
                ind = 1:numel(y_is_within_distance);
                ind = ind(y_is_within_distance);
                [~, Closest_y_ind] = min(y_dist_squared(y_is_inside));
                Closest_y_ind = ind(Closest_y_ind);
            end
            
            % For targets that are within distance, correct isinside logical based
            % on the exact distance measurement.
            y_is_within_distance(y_is_within_distance) = y_is_inside;
        else
            Closest_y_ind = [];
        end
    end
%% Function to build and plot the current tree.
    function PlotTree()
        % Build the tree structure.
        Tree = build_tree();
        
        % Define useful cursor data for debugging.
        Active_nodesID = find(Nodes_isActive);
        Cursor_data = struct();
        for k = 1:numel(Active_nodesID)
            NodeID = Active_nodesID(k);
            Node_Pos = Nodes_Pos(NodeID, :);
            Cursor_data(k).X = Node_Pos(1);
            Cursor_data(k).Y = Node_Pos(2);
            Cursor_data(k).NodeID = NodeID;
            Cursor_data(k).Node_BirthTime = Nodes_BirthTime(NodeID);
            Cursor_data(k).Node_isActive = Nodes_isActive(NodeID);
            Cursor_data(k).Node_isBranchPoint = Nodes_isBranchpoint(NodeID);
            Cursor_data(k).Node_isnexttoBP = Nodes_isNextToBP(NodeID);
            Cursor_data(k).Node_BranchID = Nodes_BranchID(NodeID);
        end
        
        Title = ['t=', num2str(t)];
        f = plottree(Tree, 'BranchID', 0, 'Lengthscale', lengthscale,...
            'CollisionSpheres', Contact_dist/2, 'ExclusionRadius', Exclusion_radius,...
            'BranchPointRadius', Branchpoint_radius, 'Bins', 1, 'TipPos', 1, ...
            'CursorData', Cursor_data, 'Title', Title, 'Boundary_position',[-Boundary_size_temp/2;Boundary_size_temp/2]);
    end
end
%% Function to build the tree structure from the branch arrays.
function Tree = build_tree()
% The tree structure is build out of the branches arrays defined globally.
%
% Output
% Tree = 1xN structure representing a tree with N branches.
global Soma_branches_ID Branches_NodesID Nodes_Pos Nodes_BirthTime ...
    Branches_ParentID Branches_SiblingID Branches_ChildrenID ...
    Branches_Length Branches_isDynamic Branches_isActive ...
    Branches_TipID Branches_BranchingAngle
N_branches = size(Branches_ParentID, 1);

Tree = struct('PointsPos', cell(N_branches, 1), ...
    'PointsBirthTime', cell(N_branches, 1), ...
    'isActive', num2cell(Branches_isActive), ...
    'isDynamic', num2cell(Branches_isDynamic), ...
    'ParentID', num2cell(Branches_ParentID), ...
    'SiblingID', num2cell(Branches_SiblingID), ...
    'ChildrenID', squeeze(num2cell(Branches_ChildrenID, 2)), ...
    'Length', num2cell(Branches_Length), ...
    'TipID', num2cell(Branches_TipID), ...
    'BranchingAngle', num2cell(Branches_BranchingAngle));

for i = 1:N_branches
    Branch_Length = Branches_Length(i);
    % Assign PointsPos to:
    % 1) Active branches that are growing.
    % 2) Active branches that have yet to grow (zero length).
    % 3) Inactive branches that have non-zero length.
    if Branches_isDynamic(i) || Branch_Length > 0
        NodesID_temp = Branches_NodesID(i, 1:Branch_Length+1);
        Tree(i).PointsPos = Nodes_Pos(NodesID_temp, :);
        Tree(i).PointsBirthTime = Nodes_BirthTime(NodesID_temp, :);
    end
end

% Fix the sibling ID for branches at the soma.
% Keep only Soma branches that have non-zero length.
Soma_branches_ID_temp = Soma_branches_ID(Branches_Length(Soma_branches_ID)>0);
Soma_branches_ID_temp = Soma_branches_ID_temp(:)';
if numel(Soma_branches_ID_temp) > 1
    for i=Soma_branches_ID_temp
        Tree(i).SiblingID = setdiff(Soma_branches_ID_temp,i);
    end
elseif numel(Soma_branches_ID_temp) > 0
    Tree(Soma_branches_ID_temp).SiblingID = Soma_branches_ID_temp;
end

% Make sure the output tree is a row structure.
Tree = reshape(Tree,1,[]);

% Set the children IDs to an empty matrix if they are < 1.
for i = 1:numel(Tree)
    ChildrenIDs = Tree(i).ChildrenID;
    if any(ChildrenIDs < 1)
        Tree(i).ChildrenID = [];
    end
end
end
%% Function that checks the angles between branches at a branchpoint.
function Theta_deg = check_BPangles(Branchpoint_nodeID)
% Input
% Branchpoint_nodeID = ID of the node representing the branchpoint.
%
% Output
% Theta_deg = 1x3 array representing angles (deg) between the branch, sibling
%             and parent. Specifically,
%             Theta_deg(1) = smallest angle between branch and sibling
%             Theta_deg(2) = smallest angle between sibling and parent
%             Theta_deg(3) = smallest angle between branch and parent
global Nodes_Pos Nodes_BranchID Branches_NodesID Branches_SiblingID Branches_ParentID Branches_Length
BranchID = Nodes_BranchID(Branchpoint_nodeID);
SibID = Branches_SiblingID(BranchID);
ParID = Branches_ParentID(BranchID);
Vecs = zeros(3, 2);

% The first vector is defined by the first 2 nodes of the branch.
Vecs(1, :) = diff(Nodes_Pos(Branches_NodesID(BranchID, 1:2), :));

% The second vector is defined by the first 2 nodes of the sibling branch.
Vecs(2, :) = diff(Nodes_Pos(Branches_NodesID(SibID, 1:2), :));

% The third vector is defined by the last 2 nodes of the parent branch.
Parent_Length = Branches_Length(ParID);
Vecs(3, :) = -diff(Nodes_Pos(Branches_NodesID(ParID, Parent_Length:Parent_Length+1), :));

% Calculate the angle between each vector.
Theta_deg = zeros(1, 3);
for i = 1:3
    Ind = mod([i, i+1]'-1, 3)+1;
    Vec_norms = sqrt(sum(Vecs(Ind, :).^2, 2));
    Theta = acos(Vecs(Ind(1), :)*Vecs(Ind(2), :)'/(Vec_norms(1)*Vec_norms(2)));
    Theta_deg(i) = Theta/pi*180;
end
disp(['The two siblings are:', num2str(BranchID), ',', num2str(SibID)]);
disp(['The angle between ', num2str(BranchID), ' and ', num2str(SibID), ' is: ', num2str(Theta_deg(1))]);
disp(['The angle between ', num2str(SibID), ' and ', num2str(ParID), ' is: ', num2str(Theta_deg(2))]);
disp(['The angle between ', num2str(ParID), ' and ', num2str(BranchID), ' is: ', num2str(Theta_deg(3))]);
disp(['Angle sum: ', num2str(sum(Theta_deg))]);
end
%% Function to run Test on the branch structures.
function Test(Problem)
global Bin_BranchID Bin_Counts Tips_Times Branches_counter
global Branch_PointsPos Branch_PointsBinInd Branch_TipPos Branch_PointsBirthTime ...
    Branches_ParentID Branches_SiblingID Branches_ChildrenID Branches_isDynamic ...
    Branches_Length Branch_Theta Branches_TipID Branch_BranchingAngle Tips_DeathPos

% Caller variables.
Caller_vars_name = evalin('caller', 'who');

% Test if branches that are sent to the retraction handler have a non-nan
% death time.
Dynamic_branches_ID = find(Branches_isDynamic(1:Branches_counter))';

switch Problem
    case 'Growth'
        % Check if the full retracted times have nan death times.
        Retracted_branches_counter = evalin('caller', 'Retracted_branches_counter');
        Retracted_branches_ID = evalin('caller', 'Retracted_branches_ID');
        Errors = isnan(Tips_Times(Branches_TipID(Retracted_branches_ID(1:Retracted_branches_counter)), 3));
        Errors_ind = find(Errors);
        if ~isempty(Errors_ind)
            Branches_ind = Retracted_branches_ID(Errors_ind);
            error('Growth error');
        end

    case {'180Reversal','180Reversal_CurrentBranch'}
        switch Problem
            case '180Reversal'
                Scanned_branches_ID = Dynamic_branches_ID;
            case '180Reversal_CurrentBranch'
                Scanned_branches_ID = evalin('caller', 'BranchID');
        end
        
        for ID = Scanned_branches_ID
            Length = Branches_Length(ID);
            if Length >= 2
                Branch_PointsPos = Branch_PointsPos(Length-1:Length+1, :, ID);
                Branch_Vecs = diff(Branch_PointsPos);

                Delta_theta = acos(Branch_Vecs(1, :)*Branch_Vecs(2, :)');
                if any(abs(Delta_theta) > 7*pi/8)
                    error('The tip position of branch %d has jumped by %.2f deg.',ID, Delta_theta/pi*180);
                end
            end
        end
end


% Test if deleted branches have non-nan tips.
if 0
    Branches_tip_nonnan = ~isnan(Branches_TipID);
    Branches_tip_nonnan(Branches_tip_nonnan) = isnan(Tips_Times(Branches_TipID(Branches_tip_nonnan), 3));
    
    if any(Branches_tip_nonnan)
        Branches_ind = find(Branches_tip_nonnan);
        Tree = build_tree();
        Tree(Branches_ind(1));
        error('Active tips are associated with deleted branches.');
    end
end

% Test if the tips have stopped when they are expected to.
if 0
    Tips_hascollided = ~isnan(Tips_Times(:, 2));
    Tips_hasstopped = ~isnan(Tips_Times(:, 3));
    Tips_deathtime_max = 2*Tips_Times(:, 2); % Retraction time can be greater than growth time.
    
    t = evalin('caller', 't');
    Erroneous_tips = Tips_hascollided & ~Tips_hasstopped & ceil(Tips_deathtime_max) < t;
    if any(Erroneous_tips)
        Tip_ind = find(Erroneous_tips)
        Branches_ind = find(ismember(Branches_TipID, Tip_ind))
        Tree = build_tree();
        Tree(Branches_ind(1));
        Tips_Times(Tip_ind(1), :);
        error('error');
    end
end

% Test if all dynamic branches have a well-defined tip ID.
if 0
    Dynamic_branches_has_nan_tip = isnan(Branches_TipID(Dynamic_branches_ID));
    if any(Dynamic_branches_has_nan_tip)
        disp(Dynamic_branches_ID(Dynamic_branches_has_nan_tip));
        error('Some dynamic branches have nan tip.')
    end
end

% Test if all the branch IDs in the bins are existent (a branch with a -2
% state is considered removed).
if 0
    Bin_BranchIDs = Bin_BranchID(Bin_Counts ~= 0);
    Bin_BranchIDs(Bin_BranchIDs == 0) = [];
    if any(Branches_State(Bin_BranchIDs) == -2)
        disp(Bin_BranchIDs(Branches_State(Bin_BranchIDs) == -2));
        error('Some branches in bins have been removed.');
    end
end

% Test if there is any active branches whose sibling branch has a length of
% zero.
if 0
    Siblings_ID = Branches_SiblingID(Dynamic_branches_ID);
    if ismember('Retracted_branches_counter', Caller_vars_name)
        Retracted_branches_counter = evalin('caller', 'Retracted_branches_counter');
        Retracted_branches_ID = evalin('caller', 'Retracted_branches_ID');
        
        Error_ind = Branches_Length(Siblings_ID) == 0 ...
            & Branches_Length(Dynamic_branches_ID) ~= 0 ...
            & ~ismember(Siblings_ID, Retracted_branches_ID(1:Retracted_branches_counter));
        if any(Error_ind)
            disp(Dynamic_branches_ID(Error_ind));
            disp(Siblings_ID(Error_ind));
            error('Branches have active sibling with zero-length.')
        end
    end
end
end
%% Function to test the persistence length of the generated paths.
function test_persistence_length(delta_theta, stepsize)
% Input
% delta_theta = 2D array representing angular changes of each step (row)
%               for each branch (col).
% stepsize = size of each step used to reconstruct the 2D paths.
N_branches = size(delta_theta, 2);
Tree_test = struct();

% Construct a fictitious tree whose branches are formed by worm-like chains
% defined by delta_theta and stepsize.
Soma_pos = [0,0];
for i = 1:N_branches
    % Define the tangent vectors that form the path of each branch.
    Theta = cumsum([0; delta_theta(:, i)]);
    Tangent_vecs = [cos(Theta), sin(Theta)];
    
    % Cumulate the tangent vectors to find the position of each node of
    % each chain.
    Nodes_Pos = cumsum([Soma_pos;Tangent_vecs]);
    Tree_test(i).PointsPos = Nodes_Pos;
end

% Calculate the persistence length of the synthetic tree.
calculate_PersistenceLength(test_tree, 'Lengthscale', stepsize,'Method','squared_dist');
%calculate_PersistenceLength(Tree_test, 'Lengthscale', stepsize,'Method','curvature');
end
%% Function to calculate the free angular ranges of a branchpoint.
function Free_intervals_boundaries = calculate_free_angular_intervals(Ori_pos, Enpoints_pos, Theta)
% This function calculates the allowed angles of line segments originating from Ori_pos,
% without intersecting the line segments represented by joining each
% enpoint position with the origin position.
% The excluded regions correspond to circular arcs that subtend an angle
% 2*Theta. They are centered at Ori_pos and oriented along each line 
% segment formed by joining Ori_pos to Enpoints_pos.
% 
% Input
%
% Ori_pos = 1x2 array
% Enpoints_pos = nx2 array where the columns represent the x,y position of
%                each line segment endpoint.
% Theta = half-angle (rad) of the prohibited angular intervals.
%
% Output
%
% Free_intervals_boundaries = Nx2 array representing the angular intervals
%                             that do not intersect with the prohibited
%                             regions. The 1st column denotes the lower
%                             bound and the 2nd column denotes the upper
%                             bound.
N_lines = size(Enpoints_pos, 1);

% Calculate the angle of the line segments.
Segments_vec = bsxfun(@minus, Enpoints_pos, Ori_pos);
Segments_angle = mod(atan2(Segments_vec(:, 2), Segments_vec(:, 1)), 2*pi);

% Define the boundaries of the excluded angular intervals.
Excluded_intervals_boundaries = bsxfun(@plus, Segments_angle, [-1, 1]*Theta);

% Define the type of each boundary. Opening = 1, Closing = -1;
Boundaries_type = ones(size(Excluded_intervals_boundaries));
Boundaries_type(:, 2) = -1;

% Find an opening boundary that is not inside an excluded region. The
% search for allowed regions will begin there.
for i = 1:N_lines
    Angle_origin_cand = Excluded_intervals_boundaries(i, 1);
    Other_boundaries_recentered = mod(Excluded_intervals_boundaries([1:i-1, i+1:end], :) - Angle_origin_cand + pi,2*pi) - pi;
    is_origin_intersecting =  Other_boundaries_recentered(:,1) < 0 & Other_boundaries_recentered(:,2) > 0;
    if all(~is_origin_intersecting)
        Angle_origin = Angle_origin_cand;
        break;
    end
end

% Recenter angles to the new angle origin.
Excluded_intervals_boundaries = mod(Excluded_intervals_boundaries-Angle_origin, 2*pi);

% Sort the boundary locations.
[Excluded_intervals_boundaries, sorted_ind] = sort(Excluded_intervals_boundaries(:));
Boundaries_type = Boundaries_type(sorted_ind);

% Add an opening boundary at 2*pi.
Excluded_intervals_boundaries = [Excluded_intervals_boundaries(:); 2*pi];
Boundaries_type = [Boundaries_type(:); 1];

% Find allowed angular intervals by locating regions where all excluded
% regions have been closed.
N_open_boundaries = cumsum(Boundaries_type);
assert(all(N_open_boundaries >= 0),'The start of the search for the calculated free range is incorrect.');
Free_region_start_ind = find(N_open_boundaries == 0);
if isempty(Free_region_start_ind)
    error('A free angular region could not be found');
end
Free_intervals_boundaries = [Excluded_intervals_boundaries(Free_region_start_ind(:)), Excluded_intervals_boundaries(Free_region_start_ind(:)+1)];

% Recenter angles to the usual origin and reorder the intervals in
% ascending order.
Free_intervals_boundaries = mod(Free_intervals_boundaries+Angle_origin, 2*pi);
[~, Begin_ind] = min(Free_intervals_boundaries(:, 1));
Free_intervals_boundaries = Free_intervals_boundaries([Begin_ind:end, 1:Begin_ind-1], :);

% Check if the last allowed interval crosses the origin. It it does, the
% opening boundary will be higher than the closing boundary (since we take 2*pi modulo).
% In this case, break the interval into two intervals separated at the origin.
if Free_intervals_boundaries(end, 1) > Free_intervals_boundaries(end, 2)
    Free_intervals_boundaries = [0, Free_intervals_boundaries(end, 2); Free_intervals_boundaries];
    Free_intervals_boundaries(end, 2) = 2*pi;
end
end
%% Function to save profile.
function save_profile(Filename)
% Input
% Filename = Full filename of the profile.

% Append current date to the filename
Filename = [Filename,'_',datestr(datetime('now'))];
fprintf('Saving profile:%s\n', Filename)
Profile_info = profile('info');
save(Filename, 'Profile_info');
profsave(Profile_info, Filename);
profile resume;
end
%% Function to determine the most memory efficient integer data type of an array.
function Data_type = min_int_type(Max_value,is_signed)
% Input
% Max_value = Maximum value that the data can take.
% is_signed = Bool indicating if the data is signed.
%
% Output
% Data_type = string corresponding to the most efficient data type to store
%             the data.
Integer_classes = {'int8','int16','int32','int64'};
if ~is_signed
    Integer_classes = cellfun(@(s) ['u',s],Integer_classes,'Uni',0);
end
%Integer_classes_min = cellfun(@(x) double(intmin(x)),Integer_classes);
Integer_classes_max = cellfun(@(x) double(intmax(x)),Integer_classes);
Data_type_ind = find(Max_value <= Integer_classes_max,1);
if ~isempty(Data_type_ind)
    Data_type = Integer_classes{Data_type_ind};
else
    warning('Could not find an efficient integer type that handles max_value=%d. Will us "double" type.',Max_value);
    Data_type = 'double';
end
end
%% Function to sample random numbers in batches.
function Rand_num = random_sampler(Type,varargin)
% Generate random numbers in batches and exhaust them when needed since 
% it is more efficient to sample random numbers in chunks in Matlab.
%
% Call options
% random_sampler(Type,'Initialize') = initialize the batches of type 'Type'.
% random_sampler(Type,Size) = sample a random number array of type 'Type' and size 'Size'.
%                             The default size is 1x1.
%
% Input
% Type = char reprensenting the type of random number to generate. Only
%        'randn' and 'rand' are supported.
%
% Output
% Rand_num = random number array
persistent rand_batch randn_batch rand_batch_counter randn_batch_counter rand_num_batch_size

% Initialize the arrays
if nargin > 1 && strcmp(varargin{1},'Initialize')
    rand_num_batch_size = 500000;
    switch Type
        case 'rand'
            rand_batch = rand(rand_num_batch_size,1);
            rand_batch_counter = 0;
        case 'randn'
            randn_batch = randn(rand_num_batch_size,1);
            randn_batch_counter = 0;
        otherwise
            error('%s random numbers not supported.')
    end
    return;
end

% Sample the numbers.
if nargin > 1 && isnumeric(varargin{1})
    % In this case, the second input corresponds to the size of the random
    % number array that is requested.
    Size = varargin{1};
    N = prod(Size);
else
    Size = [1,1];
    N = 1;
end

% Checks
if N > rand_num_batch_size
    error('The requested number of random numbers (%d) is larger than the batch size (%d)',N,rand_num_batch_size);
end

% Get the random numbrs
switch Type
    case 'rand'
        if rand_batch_counter + N > rand_num_batch_size
            % Regenerate random numbers
            random_sampler('rand','Initialize');
        end
        Rand_num = rand_batch(rand_batch_counter + (1:N));
        rand_batch_counter = rand_batch_counter + N;
    
    case 'randn'
        if randn_batch_counter + N > rand_num_batch_size
            random_sampler('randn','Initialize');
        end
        Rand_num = randn_batch(randn_batch_counter + (1:N));
        randn_batch_counter = randn_batch_counter + N;
    otherwise
        error('%s random numbers not supported.')
end

% Reshape the random number array.
if N > 1
    Rand_num = reshape(Rand_num,Size);
end
end