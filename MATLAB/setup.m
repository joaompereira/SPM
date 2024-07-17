%% To install just add folders to PATH
disp("Adding folders to path...")
fullfile = mfilename('fullpath');
filepath = fileparts(fullfile);
addpath(filepath + "/");
addpath(filepath + "/helper_functions/");