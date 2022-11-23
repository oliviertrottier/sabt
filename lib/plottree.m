function f = plottree(Tree,varargin)
% Function to plot a 2D tree represented by a tree structure.
% Calling options
% plotree(tree_struct,options)
%% Parse Inputs
% Make Tree a row structure.
Tree = Tree(:)';
%% Parse Optional Inputs
p = inputParser;
addParameter(p, 'Invisible', false); % Make the figure invisible. Useful for producing movies.
addParameter(p, 'BranchingAngles', false); % Plot the branching angle of each branch near the branchpoint.
addParameter(p, 'BranchID', false); % Plot the ID of each branch.
addParameter(p, 'BranchColor', [0,0,0]); % Change color of branches.
addParameter(p, 'BranchWidth', get(0,'DefaultLineLineWidth')); % Change color of branches.
addParameter(p, 'Coordinates', 'cartesian', @(x) ismember(x,{'cartesian','pixel'})); % Determine the coordinates of the branch nodes. cartesian=usual cartesian coordinates with positive x,y pointing towards right and top. pixel=pixel coordinates.
addParameter(p, 'NodesMarker', 'None'); % Change the marker used for the branch nodes.
addParameter(p, 'Lengthscale', 1); % Change the lengthscale of the tree. By default, no changes are made (scale = 1).
addParameter(p, 'AxisLimits', []); % 1x4 array defining the axis limits as [xmin,xmax,ymin,ymax] or single entry x defining limits as:[-x,x,-x,x].
addParameter(p, 'SquareAxisLimits', true); % Force axis limits to be square. (0,0) is at the center of the axis.
addParameter(p, 'Center', 'soma',@(x) ismember(x,{'soma','com'})); % Change the origin of the axes to the soma or center-of-mass.
addParameter(p, 'NoAxis', false); % Makes the axis invisible.
addParameter(p, 'CursorData', false); % Add cursor data in plot.
addParameter(p, 'Title', ''); % Add title.
addParameter(p, 'FigureSize', 600); % Set figure size.
addParameter(p, 'ScaleBar', []); % Scalebar size in data units or 3x1 double array defining the size (1), x position (2) and y position(3) of the scale bar in data units.
addParameter(p, 'Timestamp', []); % Timestamp in hours or 3x1 cell array defining the timestamp (1), x position (2) and y position(3) in data units.
addParameter(p, 'Boundary_position', []); % Plot position of a rectangular boundary. Boundary_position(1,:)=lower limits, Boundary_position(2,:) = upper limits
parse(p, varargin{:});
options = p.Results;

% Redefine some default options.
% Do not use square axis limits with pixel coordinates.
if strcmp(options.Coordinates,'pixel')
    options.SquareAxisLimits = false;
end
%% Initialization
% Find the branches with a non-zero length.
BranchIDs = find([Tree.Length] > 0);
N_branches = numel(BranchIDs);
Tree_fieldnames = fieldnames(Tree);
textfontsize = 18;

% Initialize the figure.
if options.Invisible
    Fig_visibility = 'off';
else
    Fig_visibility = 'on';
end

% Determine the title.
if ~isempty(options.Title)
    % Use the input title
    Title = options.Title;
else
    % Use the name of the input.
    Title = inputname(1);
end

% Fix the figure size.
f = figure('Visible', Fig_visibility);
f.Position(3:4) = options.FigureSize;

% Add the title and labels.
title(Title,'interpreter', 'none');
xlabel(['X [\mu', 'm]']);
ylabel(['Y [\mu', 'm]']);

% Get the current axis handle.
a = f.CurrentAxes;
axis equal;
hold on;

if options.NoAxis
    a.Visible = 'off';
    a.Position = [0, 0, 1, 1];
end
%% Scale the Tree with the lengthscale.
Lengthscale = options.Lengthscale;
Tree = scale_tree(Tree,Lengthscale);
%% Soma and center position
Soma_position = Tree(1).PointsPos(1,:);
Tree_nodes_position = unique(cell2mat({Tree.PointsPos}'),'rows');
switch options.Center
    case 'soma'
        Center_position = Soma_position;
    case 'com'
        [~,Center_position] = rg(Tree_nodes_position);
end
%% Tree dimensions
% Find the maximal extension of the tree in each dimension.
Tree_nodes_position = unique(cell2mat({Tree.PointsPos}'),'rows');
Tree_XY_minmax = minmax(Tree_nodes_position')';
Tree_XY_range = diff(Tree_XY_minmax);
%% Plot the branches
for i = 1:N_branches
    BranchID = BranchIDs(i);
    Length = Tree(BranchID).Length;
    switch options.Coordinates
        case 'cartesian'
            PointsPos = Tree(BranchID).PointsPos(1:Length+1, :);
        case 'pixel'
            % In pixel coordinates, the first entry correspond to the row
            % dimension (vertical dimension).
            PointsPos = Tree(BranchID).PointsPos(1:Length+1, [2 1]);
    end
    Branch_line_h = line('Parent', a, ...
        'XData', PointsPos(:, 1), 'YData', PointsPos(:, 2), ...
        'Marker', options.NodesMarker, 'MarkerSize', 8, ...
        'Color',options.BranchColor,'LineWidth',options.BranchWidth);
    
    % Insert a textbox near the branchpoint giving the branching angle.
    if options.BranchingAngles
        theta = Tree(BranchID).BranchingAngle;
        BranchingAngle_text_pos = PointsPos(1, :);
        text(BranchingAngle_text_pos(1), BranchingAngle_text_pos(2), num2str(round(theta/pi*180, 1)),...
            'FontSize', 18, 'Color', 'red', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
    % Insert a textbox in the middle of the branch giving the branch
    % ID.
    if options.BranchID
        if Length==1
            BranchID_text_pos =  mean(PointsPos(1:2, :),1);
        else
            BranchID_text_pos =  PointsPos(ceil((Length+1)/2), :);
        end
        
        text(BranchID_text_pos(1), BranchID_text_pos(2), num2str(BranchID),...
            'FontSize', 25, 'Color', [0, 0.5, 0], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% If no branches are plotted (empty tree), plot the soma of the tree.
if N_branches == 0
    PointsPos = Tree(1).PointsPos(1,:);
    Branch_line_h = line('Parent', a, ...
        'XData', PointsPos(:, 1), 'YData', PointsPos(:, 2), ...
        'Marker', options.NodesMarker, 'MarkerSize', 8, ...
        'Color',options.BranchColor);
end

% Change the direction of the Y axis according to the coordinates option.
switch options.Coordinates
    case 'cartesian'
        a.YAxis.Direction = 'normal';
    case 'pixel'
        a.YAxis.Direction = 'reverse';
end
%% Plot the boundary.
if ~isempty(options.Boundary_position)
    Boundary_vertices = [repelem(options.Boundary_position(:,1),2), options.Boundary_position([1,2,2,1],2)];
    Boundary_vertices = [Boundary_vertices; Boundary_vertices(1,:)];
    plot(Boundary_vertices(:,1),Boundary_vertices(:,2),'r-');
end
%% Insert time stamp
if ~isempty(options.Timestamp)
    if iscell(options.Timestamp)
        Timestamp_string = sprintf('%.1f hrs AEL',options.Timestamp{1});
        Timestamp_pos = [options.Timestamp{2:3}];
    else
        Timestamp_string = sprintf('%.1f hrs AEL',options.Timestamp);
        Timestamp_pos = [Tree_XY_minmax(2,1),Tree_XY_minmax(2,2)];
    end
    
    Timestamp_text = text(Timestamp_pos(1), Timestamp_pos(2), Timestamp_string, ...
        'HorizontalAlignment', 'right','VerticalAlignment','bottom');
    Timestamp_text.FontSize = textfontsize;
else
    % When there is no timestamp, initialize the timestamp position at the
    % center since it is needed for calculating scale bar position later.
    Timestamp_pos = Center_position;
end
%% Insert scale bar
if options.ScaleBar
    % Define the color of the bar (black or white).
    Scale_bar_color = [0,0,0];
    
    % Insert the legend of the scale bar first.
    Scale_bar_legend_pos = [Tree_XY_minmax(2,1) - options.ScaleBar(1)/2 Tree_XY_minmax(1,2)];
    Scale_bar_text = text(Scale_bar_legend_pos(1), Scale_bar_legend_pos(2), [num2str(options.ScaleBar(1)), ' \mum'],...
        'HorizontalAlignment', 'center','VerticalAlignment','top', 'FontSize', textfontsize,'Color',Scale_bar_color);
    
    % Raise the text box to tighten it to the structure.
    is_node_in_vicinity = Tree_nodes_position(:,1) > Scale_bar_text.Extent(1) & Tree_nodes_position(:,1) < Scale_bar_text.Extent(1) + Scale_bar_text.Extent(3);
    Y_vicinity_min = min(Tree_nodes_position(is_node_in_vicinity,2));
    Scale_bar_text.Position(2) = Y_vicinity_min;
    
    % Insert the scale bar right below the text box.
    if numel(options.ScaleBar)<=1
        Scale_bar_pos_datau = [Scale_bar_legend_pos(1) - options.ScaleBar(1)/2, Scale_bar_text.Extent(2) - .1*Scale_bar_text.Extent(4)];
    else
        Scale_bar_pos_datau = options.ScaleBar(2:3);
    end
    
    plot(a,[0,options.ScaleBar(1)] + Scale_bar_pos_datau(1),Scale_bar_pos_datau(2)*ones(1,2),'-','LineWidth', 5,'Color',Scale_bar_color);
    Scale_bar_bounding_box = [Scale_bar_pos_datau, options.ScaleBar(1), 0.25*Scale_bar_text.Extent(4)];
else
    Scale_bar_bounding_box = [Center_position,0,0];
end
%% Fix axis limits.
drawnow;

% Adjust the axis limits.
Axis_limits_division = 5;
if isempty(options.AxisLimits)
    min_axis_limits = Axis_limits_division;
    X_min = floor(min([Scale_bar_bounding_box(1),Timestamp_pos(1),1.025*Tree_XY_minmax(1,1) - 1])/Axis_limits_division)*Axis_limits_division;
    Y_min = floor(min([Scale_bar_bounding_box(2) - Scale_bar_bounding_box(4),Timestamp_pos(2),1.025*Tree_XY_minmax(1,2) - 1])/Axis_limits_division)*Axis_limits_division;
    X_max = ceil(max([(Scale_bar_bounding_box(1) + Scale_bar_bounding_box(2)),Timestamp_pos(1),1.025*Tree_XY_minmax(2,1) + 1])/Axis_limits_division)*Axis_limits_division;
    Y_max = ceil(max([Scale_bar_bounding_box(2),Timestamp_pos(2),1.025*Tree_XY_minmax(2,2) + 1])/Axis_limits_division)*Axis_limits_division;
    
    XLimits = [X_min, X_max];
    YLimits = [Y_min, Y_max];
    if options.SquareAxisLimits
        XLimits = max(abs(XLimits - Center_position(1)))*[-1 1];
        YLimits = max(abs(YLimits - Center_position(2)))*[-1 1];
    end
    
    % Translate the limits using the positon of the center.
    XLimits = XLimits + Center_position(1);
    YLimits = YLimits + Center_position(2);
    
    xlim(a,XLimits);
    ylim(a,YLimits);
else
    % Set the axis limits to the values given as optional parameters.
    if numel(options.AxisLimits) == 1
        options.AxisLimits = options.AxisLimits*[-1,1,-1,1];
    end
    xlim(a,options.AxisLimits(1:2));
    ylim(a,options.AxisLimits(3:4));
end
%% Configure the figure data tooltip.
% Save the cursor data in the userdata field of the figure.
Cursor_data_fieldnames = {'Depth','PointsDiameter','ParentID','ChildrenID','Tangential_angle','Length_dim','Length'};
Cursor_data_fieldnames = Cursor_data_fieldnames(isfield(Tree,Cursor_data_fieldnames));
N_Cursor_data_fieldnames = numel(Cursor_data_fieldnames);

if isa(options.CursorData, 'struct')
    CursorData = options.CursorData;
else
    CursorData = struct();
    z = 0;
    for i = 1:N_branches
        BranchID = BranchIDs(i);
        BranchLength = Tree(BranchID).Length;
        for j=1:BranchLength+1
            z = z+1;
            switch options.Coordinates
                case 'cartesian'
                    CursorData(z).X = Tree(BranchID).PointsPos(j,1);
                    CursorData(z).Y = Tree(BranchID).PointsPos(j,2);
                case 'pixel'
                    CursorData(z).X = Tree(BranchID).PointsPos(j,2);
                    CursorData(z).Y = Tree(BranchID).PointsPos(j,1);
            end
            
            CursorData(z).ID = BranchID;
            CursorData(z).BranchLength = BranchLength;
            
            % Add optional fieldnames.
            for k=1:N_Cursor_data_fieldnames
                switch Cursor_data_fieldnames{k}
                    case {'PointsDiameter','Tangential_angle'}
                        CursorData(z).(Cursor_data_fieldnames{k}) = Tree(BranchID).(Cursor_data_fieldnames{k})(j);
                    otherwise
                        CursorData(z).(Cursor_data_fieldnames{k}) = Tree(BranchID).(Cursor_data_fieldnames{k});
                end
            end
        end
    end
end
if isempty(f.UserData)
    f.UserData = CursorData;
else
    f.UserData = structcat(f.UserData,CursorData);
end


% Set the datacursor mode function of the figure.
dcm_obj = datacursormode(f);
set(dcm_obj,'Interpreter','none');
set(dcm_obj, 'UpdateFcn', @cursor_fun);
%% Remove hold.
hold off
end
%% Custom cursor function
function Output_txt = cursor_fun(obj, event_obj)
% Display the position of the data cursor
%
% Input
%
% obj = Currently not used (empty)
% event_obj = Handle to event object
% output_txt = Data cursor text string (string or cell array of strings).

% Get the cursor position.
pos = event_obj.Position;

% Get the plot's data.
Plot = event_obj.Target;
fig_h = Plot.Parent.Parent;
UserData = fig_h.UserData;

% Assign the X and Y values in the output_txt.
Output_txt = cell(1, 2);
Output_txt{1} = ['X: ', num2str(pos(1), 4)];
Output_txt{2} = ['Y: ', num2str(pos(2), 4)];

% Find all the fields to be displayed;
if ~isempty(UserData)
    % Find the index of the point.
    Data_ind = find(pos(1) == [UserData.X] & pos(2) == [UserData.Y],1);
    
    % Get the relevant user data.
    if ~isempty(Data_ind)
        UserData = UserData(Data_ind);
        
        % Find the information to show in the tooltip.
        Fieldvalues = cellfun(@num2str, struct2cell(UserData), 'Uni', 0);
        Fieldnames = fieldnames(UserData);
        Output_txt2 = strcat(Fieldnames, {': '}, Fieldvalues);
        
        % Append it to the default tooltip.
        Output_txt = [Output_txt2(:)'];
    else
        Output_txt = [Output_txt, {'Tree data not found'}];
    end
end
end