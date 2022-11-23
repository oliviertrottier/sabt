function createmovie(Movie_data, Movie_filename, varargin)
% Function to create a movie of tree growth from a sequence of tree
% structures.
%
% Input
% Movie_data = scalar structure containing data for the movie.
%              The structure must contain the "timeframes" field that
%              describe the tree structure at each frame. Specifically,
%              Movie_data.timeframes(i).tree denotes the tree structure of
%              the i^th frame. An optional Movie_data.timeframes(i).time
%              field may also be used to denote the time of the i^th frame.
%% Check inputs
% Extract the movie timeframe structures.
assert(isfield(Movie_data,'timeframes'),'The input movie data require the "timeframes" field');
Timeframes = Movie_data.timeframes;
%% Parse the optional parameters
p = inputParser;
addParameter(p, 'AxisVisibility', 'on'); % Control axis visiblity.
addParameter(p, 'Boundary', isfield(Timeframes,'boundary_size')); % Option to plot the boundary.
addParameter(p, 'Framerate', 30); % Movie framerate (N_frames/sec).
addParameter(p, 'FontSize', 30); % Font size of text boxes.
addParameter(p, 'FrameDimension', [1280, 1280]); % Dimension of the frame in pixels.
addParameter(p, 'InvertColor', false); % Invert the color of the movie.
addParameter(p, 'Lengthscale', 0); % Set the lengthscale for plotting the trees.
addParameter(p, 'Lmax', 0); % Set the spatial limits.
addParameter(p, 'MovieLength', 10); % Length of movie in seconds.
addParameter(p, 'SaveLastFrame', true); % Save the last frame of the movie in a separate image.
addParameter(p, 'ScaleBarSize', 'Auto'); % Set the size of the scale bar in \mum. An automatic size is determined based on the size of the frames.
addParameter(p, 'TitleFrame', 1); % Add a title frame at the beginning of the movie.
addParameter(p, 'TimeOffset', duration(14,0,0),@isduration); % Temporal offset added to the displayed timestamp. 14 hours correspond to the approximate start of morphogenesis after egg-lay (AEL) in class IV neurons.
addParameter(p, 'Timescale', 5); % Default duration (minutes) between two frames.
parse(p, varargin{:});
options = p.Results;
options_isdefault = cell2struct(num2cell(ismember(p.Parameters, p.UsingDefaults)), p.Parameters, 2);
%% Initialization

% FontSize
Linewidth = 2;

% Determine if the trees were created off-lattice.
Tree_fieldnames = fieldnames(Timeframes(1).tree);
Off_lattice = ismember('PointsPos', Tree_fieldnames);

% Lengthscale
if options.Lengthscale
    Lengthscale = options.Lengthscale;
else
    if Off_lattice
        Lengthscale = 1;
    else
        Lengthscale = 0.4;
    end
end

% Timeframes definitions
Timescale = options.Timescale;
N_timeframes = numel(Timeframes);
if isfield(Timeframes,'time')
    % If a "time" field is defined in Timeframes, assume the time is given in 
    % minutes.
    Timeframes_time = duration(0,[Timeframes.time],0);
else
    Timeframes_time = duration(0,N_timeframes*Timescale,0); % hours
end

% Shift the time of the frames by the given offset.
Time_initial = options.TimeOffset;
Timeframes_time = Timeframes_time + options.TimeOffset;
Timeframes_time.Format = 'hh:mm';

% Lattice definitions.
if ~Off_lattice
    L = N_timeframes+100;
    N = round(4*L);
    d = Lengthscale;
    dy = d/2;
    dx = sqrt(3)/2*d;
    [x, y] = meshgrid(-N/2:1:N/2, N/2:-1:-N/2);
    x = x.*dx;
    y = y.*dy;
    occ = x.*0;
end
%% Initialize figure
f = figure('Position', [0, 0, options.FrameDimension], 'Visible', 'off');
colormap('gray');
axis_h = gca;
hold on;

% Determine the maximum size of the tree. Lmax denotes the maximum length
% between the center and the farthest point of the tree.
if options.Lmax > 0
    Lmax = options.Lmax;
else
    All_trees = [Timeframes.tree];
    if Off_lattice
        % Find the maximum and minimum position of all nodes (points) of
        % all trees. The input tree nodes position are assumed to be
        % dimensionless.
        PointsPos = cell2mat({All_trees.PointsPos}');
        PointsPos_max = max(PointsPos)*Lengthscale;
        PointsPos_min = min(PointsPos)*Lengthscale;
        Lmax = double(round(max(abs([PointsPos_max, PointsPos_min]))));
    else
        % Determine the maximum extent of the tree from the farthest
        % occupied site on the lattice.
        Occ_ind = cell2mat({All_trees.PointsInd}');
        Occ_rows = Occ_ind(:, 1);
        Occ_cols = Occ_ind(:, 2);
        Occ_minrow = min(Occ_rows);
        Occ_maxrow = max(Occ_rows);
        Occ_mincol = min(Occ_cols);
        Occ_maxcol = max(Occ_cols);
        Occ_xmin = x(1, Occ_mincol);
        Occ_xmax = x(1, Occ_maxcol);
        Occ_ymin = y(Occ_maxrow, 1);
        Occ_ymax = y(Occ_minrow, 1);
        Lmax = max(abs([Occ_xmin, Occ_xmax, Occ_ymin, Occ_ymax]));
    end
    
    % Make sure that Lmax is as large as the boundary (if any).
    if options.Boundary
        Boundary_sizes = cell2mat({Timeframes.boundary_size}');
        Lmax = max([Lmax;Boundary_sizes(:)/2]);
    end
end

% Format axes.
xlabel('x (\mum)')
ylabel('y (\mum)')
TickMax = round(Lmax, 1, 'Significant');
axis_h.XTick = -TickMax:TickMax/5:TickMax;
axis_h.YTick = axis_h.XTick;
axis_h.FontSize = options.FontSize;
axis_h.Visible = options.AxisVisibility;
axis(axis_h,'equal'); axis(axis_h,'square');
xlim(Lmax*1.05*[-1, 1]);
ylim(Lmax*1.05*[-1, 1]);

% Make sure the axis uses the entire frame when the axis is invisible.
if strcmp(axis_h.Visible, 'off')
    axis_h.Position = [0, 0, 1, 1];
end
%% Initialize movie
% Determine the full filename of the movie
is_filename_full = strcmp(Movie_filename(1),filesep);
if is_filename_full
    Movie_filename_full = Movie_filename;
    [~,Movie_filename] = fileparts(Movie_filename_full);
else
    % Define the default folder where movies are saved.
    Moviesfolder_default = fileparts(mfilename('fullpath'));
    Movie_filename_full = fullfile(Moviesfolder_default, Movie_filename);
end

% Determine extension.
[Movie_folder,Movie_name,Movie_ext] = fileparts(Movie_filename_full);

% Create the folder where the movie is written if it doesn't exist.
if ~isdir(Movie_folder)
    mkdir(Movie_folder);
end

% Use .mp4 for the default movie format.
Available_videoprofiles = {VideoWriter.getProfiles().Name};
switch Movie_ext
    case {'','.mp4'}
        videoprofile = 'MPEG-4';
        % Use .avi if .mp4 is unavailable.
        if ~ismember(videoprofile, Available_videoprofiles)
            warning('%s is unavailable. Switching to Motion JPEG AVI');
            videoprofile = 'Motion JPEG AVI';
        end
    case '.avi'
        videoprofile = 'Motion JPEG AVI';
    otherwise
        error('Movie extension %s is not supported.',Movie_ext);
end



% Initialize movie writer.
vid = VideoWriter(Movie_filename_full, videoprofile);
vid.FrameRate = options.Framerate;
open(vid);

% Determine indices of timeframes that will be used in the video.
N_recorded_frames = options.Framerate*options.MovieLength - 1; % -1 to account for initial frame.
Recorded_frames_ind = round(linspace(1,N_timeframes,N_recorded_frames));

% Initialize the structure where frames are saved.
Recorded_frames(1:N_recorded_frames) = struct('cdata', [], 'colormap', []);
%% Title frame
if options.TitleFrame
    % Remove the axes.
    if strcmp(options.AxisVisibility,'on')
        axis_h.Visible = 'off';
    end
    
    % Write the input parameters in the title.
    Title_str = {};
    if isfield(Movie_data, 'alpha')
        Title_str = {['\alpha = ', num2str(Movie_data.alpha*Lengthscale), ' \mum  (', num2str(Movie_data.alpha), ')']};
    end
    Title_str = [Title_str, ['\beta = ', num2str(Movie_data.beta*Lengthscale), ' \mum  (', num2str(Movie_data.beta), ')']];
    if isfield(Movie_data, 'gamma')
        Title_str = [Title_str, ['\omega = ', num2str(Movie_data.gamma/Timescale), ' min^-^1 (', num2str(Movie_data.gamma), ')']];
    elseif isfield(Movie_data, 'omega')
        Title_str = [Title_str, ['\omega = ', num2str(Movie_data.omega/Timescale), ' min^-^1 (', num2str(Movie_data.omega), ')']];
    end
    
    
    % Write the legend.
    Legend_str = {'Input parameters are given in real units and (dimensionless units).'};
    Legend_str = [Legend_str, '\alpha = Retraction scale'];
    Legend_str = [Legend_str, '\beta = Persistence Length'];
    Legend_str = [Legend_str, '\omega = Branching Frequency'];
    
    if Off_lattice
        Title_str = [Title_str, ['\delta = ', num2str(Movie_data.delta*Lengthscale), ' \mum (', num2str(Movie_data.delta), ')']];
        Legend_str = [Legend_str, '\delta = Step Size'];
    end
    
    % Insert the title and legend.
    inputinfo = text(-0.9*Lmax, 0.2*Lmax, Title_str, 'FontSize', 70);
    legend = text(-0.9*Lmax, -0.7*Lmax, Legend_str, 'FontSize', 20);
    
    % Write the title frame to the video.
    writeVideo(vid, getframe(f));
    
    % Delete the title and legend and turn the axes back on.
    delete(inputinfo);
    delete(legend);
    axis_h.Visible = options.AxisVisibility;
end
%% Time stamp and scale bar
% Iniatilize time stamp.
Timestamp_str = [char(Time_initial), ' AEL'];
Timestamp_text = text(Lmax, 0.85*Lmax, Timestamp_str, 'HorizontalAlignment', 'right');
Timestamp_text.FontSize = options.FontSize;

% Add scalebar.
if strcmp(options.ScaleBarSize,'Auto') || options.ScaleBarSize > 0
    xrange = diff(axis_h.XLim);
    yrange = diff(axis_h.YLim);
    
    % If the automatic scalebar size is chosen, choose a length that fits
    % within the frame.
    if strcmp(options.ScaleBarSize,'Auto')
        Frame_size_x_cand = [1,5,10,50,100];
        Frame_size_x_cand_ind = find(Frame_size_x_cand < 0.8*xrange, 1,'last');
        Scale_bar_size = Frame_size_x_cand(Frame_size_x_cand_ind);
    else
        Scale_bar_size = options.ScaleBarSize;
    end
    Scale_bar_length_axisu = Scale_bar_size/xrange;
    Scale_bar_length_figu = Scale_bar_length_axisu*axis_h.Position(3);
    Scale_bar_pos_axisu = [1 - Scale_bar_length_axisu - 0.1, 0.05];
    Scale_bar_pos_datau = Scale_bar_pos_axisu.*[xrange, yrange]+[axis_h.XLim(1), axis_h.YLim(1)];
    Scale_bar_pos_figu = du2nfu(axis_h, Scale_bar_pos_datau);
    
    Scale_bar_legend_pos = [Scale_bar_pos_datau(1) + (Scale_bar_length_axisu/2)*xrange, Scale_bar_pos_datau(2) + (0.50*Scale_bar_pos_axisu(2))*yrange];
    Scale_bar_text_h = text(Scale_bar_legend_pos(1), Scale_bar_legend_pos(2), [num2str(Scale_bar_size), ' \mum'], 'HorizontalAlignment', 'center', 'FontSize', options.FontSize);
    annotation('arrow', [Scale_bar_pos_figu(1), Scale_bar_pos_figu(1) + Scale_bar_length_figu], [Scale_bar_pos_figu(2), Scale_bar_pos_figu(2)], 'HeadStyle', 'none', 'LineWidth', 10);
end

%% Add boundary
if options.Boundary
    Boundary_sizes = cell2mat({Timeframes.boundary_size}');
    if size(Boundary_sizes,2)==1
        Boundary_sizes = repmat(Boundary_sizes,1,2);
    end
    Boundary_positions = [-Boundary_sizes/2 Boundary_sizes];
    
    % Define the color of the boundary based on its type.
    switch Movie_data.Boundary_type
        case -1
            % Red color for repulsing boundary.
            Boundary_color = [1,0,0];
        case 0
            % Yellow color for absorbing boundary.
            Boundary_color = [234,219,42]/255;
        case 1
            % Blue color for periodic boundary.
            Boundary_color = [43,140,190]/255;
        otherwise
            error('Boundary type (%d) is not supported',Movie_data.Boundary_type);
    end
    
    % If colors are inverted, invert the color of the boundary such that
    % its color is preserved when the color is inverted again when the
    % frame is recorded.
    if options.InvertColor
        Boundary_color = 1 - Boundary_color;
    end
    
    boundary_h = rectangle(axis_h,'Position',Boundary_positions(1,:),'EdgeColor',Boundary_color,'Linewidth',Linewidth);
end
%% Write frames
% Write initial frame.
initial_frame = getframe(f);
if options.InvertColor
    initial_frame.cdata = uint8(255) - initial_frame.cdata;
end
writeVideo(vid, initial_frame);

% Record the rest of the timeframes and save them in the video file.
wb = parwaitbar(N_recorded_frames, 'WaitMessage',['Generating movie: ', Movie_filename]);

for n = 1:N_recorded_frames
    % Plot an invisible figure with the tree.
    Frame_ind = Recorded_frames_ind(n);
    if Off_lattice
        f2 = plottree(Timeframes(Frame_ind).tree, 'Invisible', 1, 'Lengthscale', Lengthscale);
    else
        f2 = plottree(Timeframes(Frame_ind).tree, occ, 'Invisible', 1, 'Lengthscale', Lengthscale);
    end
    
    % Transfer branches handle to the main figure.
    Lines = f2.CurrentAxes.Children;
    [Lines.Parent] = deal(axis_h);
    
    % Increase linewidth of lines.
    [Lines.LineWidth] = deal(Linewidth);
    
    % Close the invisible figure.
    close(f2);
    
    % Update the size of the boundary.
    if options.Boundary
        boundary_h.Position = Boundary_positions(Frame_ind,:);
    end
    
    % Print the timestamp.
    Timestamp_text.String = [char(Timeframes_time(Frame_ind)) ' AEL'];
    
    % Record frame and delete the lines.
    Recorded_frames(n) = getframe(f);
    if options.InvertColor
        Recorded_frames(n).cdata = uint8(255)-Recorded_frames(n).cdata;
    end
    delete(Lines);
    
    wb.progress;
end

% Write the frames into the video and close it.
writeVideo(vid, Recorded_frames);
close(vid);
close(f);

% Save the last frame if needed.
if options.SaveLastFrame
    imwrite(Recorded_frames(end).cdata, fullfile(Movie_folder,[Movie_name, '_lastframe.jpg']));
end
end