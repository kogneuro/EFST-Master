function [p_data] = exp_init_visual(exp_dir, thr_dir)
% exp_init_visual
% ----------------------------------------------------------
% Sets up the experiment environment for visual tasks.
% This includes:
%   • Adding all required folders to MATLAB path
%   • Initializing participant data (folder, ID, timestamp)
%   • Starting a diary logfile for debugging and traceability
%
% Inputs:
%   exp_dir  – name of experiment directory (e.g., 'respirationCA')
%   thr_dir  – name of threshold/aux directory (e.g., 'thr1F')
%
% Output:
%   p_data   – struct with participant metadata:
%                 • p_data.ID
%                 • p_data.dir   (participant folder)
%                 • p_data.ts    (timestamp)
%
% Author: Martin Grund
% Updated & commented by Alp (2025)
% ----------------------------------------------------------

%% ---- 1) Add required paths to MATLAB search path ----

% Add *experiment directory* and all its subfolders
addpath(genpath([pwd, '/', exp_dir]));

% Add *threshold directory* (if used) and all subfolders
addpath(genpath([pwd, '/', thr_dir]));

% Add "assets" folder (toolboxes, helper functions, images, etc.)
addpath(genpath([pwd, '/assets']));


%% ---- 2) Initialize participant data ----
% This function:
%    • asks for participant ID
%    • creates a data folder for that ID
%    • creates a timestamp
%    • returns all info as a struct (p_data)
%
% Example result:
%   p_data.ID  = 'S01'
%   p_data.dir = 'data/respirationCA/ID_S01/'
%   p_data.ts  = '2025_02_02_14_33'
%
p_data = participant_data(['data/', exp_dir, '/ID']);


%% ---- 3) Start diary logfile ----
% Every MATLAB command and error during the experiment is saved.
% Extremely important for debugging timing issues or crashes.
%
% Logfile example:
%   data/respirationCA/ID_S01/exp_S01_log.txt
%
diary([p_data.dir 'exp_' p_data.ID '_log.txt']);
end % ------------------------------------------------------
