%% Preliminary definitions
Script_folder = fileparts(which('Examples.m'));
Movies_folder = fullfile(Script_folder,'movies');
%% Example 1:
% Summary of example objective
[Params,Optional_params] = define_params();
[Tree1,~,~,Sim1_params] = sabt(Params,Optional_params{:});
plottree(Tree1,'lengthscale',Sim1_params.lengthscale);
%% Create movie of example 1
[Tree1,~,~,Sim1_params,~,Movieframes] = sabt(Params,Optional_params{:},'waitbar',1,'Movie',1);
Moviedata = Sim1_params;
Moviedata.timeframes = Movieframes;
Moviedata.Movie_full_filename = fullfile(Movies_folder,'Tree1.mp4');
createmovie(Moviedata,Moviedata.Movie_full_filename,...
    'TitleFrame',0,'FrameDimension',[800 800],...
    'AxisVisibility','off','InvertColor',1,...
    'Lengthscale',Moviedata.lengthscale);
%% Example 2: tree growth with extensive branching
[Params,Optional_params] = define_params('Extensive_branching_with_decay');
[Tree2,~,~,Sim2_params] = sabt(Params,Optional_params{:},'waitbar',1);
plottree(Tree2,'lengthscale',Sim2_params.lengthscale);
%% Create movie of example 2
[Tree2,~,~,Sim2_params,~,Movieframes] = sabt(Params,Optional_params{:},'waitbar',1,'Movie',1);
Moviedata = Sim2_params;
Moviedata.timeframes = Movieframes;
Moviedata.Movie_full_filename = fullfile(Movies_folder,'Tree2.mp4');
createmovie(Moviedata,Moviedata.Movie_full_filename,...
    'TitleFrame',0,'FrameDimension',[800 800],...
    'AxisVisibility','off','InvertColor',1,...
    'Lengthscale',Moviedata.lengthscale);
%% Example 3: tree growth with extensive branching and boundary
[Params,Optional_params] = define_params('Extensive_branching_with_decay_and_boundary');
[Tree3,~,~,Sim3_params] = sabt(Params,Optional_params{:},'waitbar',1);
plottree(Tree3,'lengthscale',Sim3_params.lengthscale);
%% Create movie of example 3
[Tree3,~,~,Sim3_params,~,Movieframes] = sabt(Params,Optional_params{:},'waitbar',1,'Movie',1);
Moviedata = Sim3_params;
Moviedata.timeframes = Movieframes;
Moviedata.Movie_full_filename = fullfile(Movies_folder,'Tree3.mp4');
createmovie(Moviedata,Moviedata.Movie_full_filename,...
    'TitleFrame',0,'FrameDimension',[800 800],...
    'AxisVisibility','off','InvertColor',1,...
    'Lengthscale',Moviedata.lengthscale);
%% Example 4: Drift & diffusion tree growth with extensive branching and boundary.
[Params4,Optional_params4] = define_params('Drift_diff_extensive_branching_with_decay_and_boundary');
%%
[Tree4,~,~,Sim4_params] = sabt(Params4,Optional_params4{:},'waitbar',1);
plottree(Tree4,'lengthscale',Sim4_params.lengthscale);
%% Create movie of example 4
[Tree4,~,~,Sim4_params,~,Movieframes] = sabt(Params4,Optional_params4{:},'waitbar',1,'Movie',1);
%%
Moviedata = Sim4_params;
Moviedata.timeframes = Movieframes;
Moviedata.Movie_full_filename = fullfile(Movies_folder,'Tree4.mp4');
createmovie(Moviedata,Moviedata.Movie_full_filename,...
    'TitleFrame',0,'FrameDimension',[800 800],...
    'AxisVisibility','off','InvertColor',1,...
    'Lengthscale',Moviedata.lengthscale);