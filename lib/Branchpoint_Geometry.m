%% Script to demonstrate geometry at a branch point.
%%
epsilon = 0.4; % Branch Thickness

% Define the radius of the branch point.
% 1/sqrt(2)*epsilon corresponds to the case where branches branch at 90 deg.
Branchpoint_radius = epsilon;
Branchpoint_radius = 1.5*epsilon;
Branchpoint_radius = 2*epsilon;
Branchpoint_radius = 1.9319*epsilon;
%Branchpoint_radius = (7*rand+1/sqrt(2))*epsilon; 
%Branchpoint_radius = 1/sqrt(2)*epsilon;
Branch_length = 3*max(epsilon,Branchpoint_radius);

% Legend structure.
H = struct('handle',{},'Legend',{});
%% Plot the skeleton of the first branch.
skeleton1x = [0 0];
skeleton1y = [0 1]*Branch_length;

Branch1LeftEdgex = skeleton1x-epsilon/2;
Branch1LeftEdgey = skeleton1y;
Branch1RightEdgex = skeleton1x+epsilon/2;
Branch1RightEdgey = skeleton1y;

figure
plot(skeleton1x,skeleton1y,'k')
hold on
plot(Branch1LeftEdgex,Branch1LeftEdgey,'--k')
plot(Branch1RightEdgex,Branch1RightEdgey,'--k')
%% Determine the minimal branching angle.
alpha = asin(epsilon/2/Branchpoint_radius);
Branching_Angle_min = 2*alpha;
Branching_Angle_min_abs = pi/2 - Branching_Angle_min;
%% Plot the skeleton of the second branch.
skeleton2x = [0,cos(Branching_Angle_min_abs)]*Branch_length;
skeleton2y = [0,sin(Branching_Angle_min_abs)]*Branch_length;

perpendicular_direction = [cos(Branching_Angle_min_abs+pi/2),sin(Branching_Angle_min_abs+pi/2)];
Branch2LeftEdgex = skeleton2x + epsilon/2*perpendicular_direction(1);
Branch2LeftEdgey = skeleton2y + epsilon/2*perpendicular_direction(2);
Branch2RightEdgex = skeleton2x - + epsilon/2*perpendicular_direction(1);
Branch2RightEdgey = skeleton2y - epsilon/2*perpendicular_direction(2);

plot(skeleton2x,skeleton2y,'k')
plot(Branch2LeftEdgex,Branch2LeftEdgey,'--k')
plot(Branch2RightEdgex,Branch2RightEdgey,'--k')
axis square equal
%% Plot an epsilon circle around the branch point.
xpos = skeleton1x(1);
ypos = skeleton1y(1);
clear h
h.handle = viscircles([xpos,ypos],epsilon/2,'EdgeColor','r');
h.Legend = ['Collision circle (r = ',num2str(epsilon/2,'%.2f'),')'];
H = [H,h];

% Plot the axis of the absorbing parallelogram.
%plot(skeleton1x,[0 epsilon/(2*sin(Branching_Angle_min))],'b')
%plot([0 epsilon/2],tan(pi/2-Branching_Angle_min)*[0 epsilon/2],'b')

% Plot the branch point circle. In this region, epsilon circles can
% overlap.
h.handle = viscircles([0,0],Branchpoint_radius,'EdgeColor','g');
h.Legend = ['Branch Point Circle (r = ',num2str(Branchpoint_radius,'%.2f'),')'];
H = [H,h];
plot([0,epsilon/2],tan(pi/2-Branching_Angle_min/2)*[0,epsilon/2],'g')

% Plot the arc of the minimal branching angle.
r_arc = Branchpoint_radius/2;
theta_arc = linspace(pi/2-Branching_Angle_min,pi/2,100)';
arc_xy = r_arc.*[cos(theta_arc) sin(theta_arc)];
h.handle = plot(arc_xy(:,1),arc_xy(:,2),'m');
h.Legend = ['Minimal branching angle = ',num2str(Branching_Angle_min/pi*180,'%.2f'),'\circ'];
H = [H,h];

% Add legend
Legend([H.handle],H.Legend,'Location','Northwest')
hold off