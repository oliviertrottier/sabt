% Add library and movies folder to the search path
main_folder = fileparts(mfilename('fullpath'));
lib_folder = fullfile(main_folder,'lib');
data_folder = fullfile(main_folder,'data');
movies_folder = fullfile(main_folder,'movies');
addpath(lib_folder);
addpath(data_folder);
addpath(movies_folder);