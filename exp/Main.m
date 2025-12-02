%% ========================================================================
%  MAIN SCRIPT – Master Stroop / EFST Launcher
%
%  Description
%  -----------
%  High-level launcher script for running one of three experimental tasks:
%
%    1) Color Stroop (3 ink colors)
%    2) Emotional Face Stroop – 2 emotions (happy / sad)
%    3) Emotional Face Stroop – 3 emotions (happy / neutral / sad)
%
%  The script:
%    - Asks for session parameters via a GUI (task type, screen, EEG, eye tracking,
%      number of blocks, trials per block),
%    - Initializes participant info and logging,
%    - Calls the appropriate setup + sequence generator for the chosen task,
%    - Saves the settings and full trial sequence for reproducibility,
%    - Runs each block with the matching run-function (Color Stroop / EFST 2 / EFST 3),
%    - Computes timing intervals and saves per-block data.
%
%  Usage
%  -----
%    - Adjust the `cd(...)` and `addpath(...)` lines to point to your EFST
%      experiment folder.
%    - Make sure the following functions are on the MATLAB path:
%         exp_init_visual.m
%         setup_colorStroop_3c.m,   seq_colorStroop_3c.m,   run_exp_colorStroop_3c.m
%         setup_EFST_2s.m,          seq_EFST_2s.m,          run_exp_EFST_2s.m
%         setup_EFST_3s.m,          seq_EFST_3s.m,          run_exp_EFST_3s.m
%         intervals.m,              save_exp.m
%    - Run this script from MATLAB (F5 or “Run”), then choose task / options
%      in the popup window.
%
%  Outputs
%  -------
%    - One settings + sequence file per session:
%         [file_prefix 'settings_seq.mat']
%    - One MAT file per block with behavioral + timing data:
%         [file_prefix 'block_XX_ID.mat']
%
%  Notes
%  -----
%    - EEG and eyetracking flags are passed into the task-specific setup
%      and run functions via the struct `s`.
%    - The same block loop is used for all three tasks; only the setup,
%      sequence, and run functions change depending on `taskType`.
%
%  Author
%  ------
%    Alp, 2025
%% ========================================================================

clearvars;                     % clear all variables from workspace
close all;                     % close all open figure windows
clc;                           % clear command window

cd('C:\Users\xyz\Desktop\...\....');  % change working directory to experiment folder

%% 1) SESSION / RUN INFO GUI
prompt  = { ...                                     % cell array of prompt strings for input dialog
    'Task type: 1=Color Stroop, 2=EFST (2 emotions), 3=EFST (3 emotions):', ...
    'Screen number (usually 2):', ...
    'Record EEG? (1=yes,0=no):', ...
    'Eyetracking? (1=yes,0=no):', ...
    'Number of blocks:', ...
    'Trials per block:'};
title   = 'Stroop / EFST Session Setup';            % title of the input dialog

%        task, screen, EEG, eye, blocks, trials
defans  = {'2','1','0','0','1','40'};              % default answers for each input field (as strings)

answer  = inputdlg(prompt, title, 1, defans);      % show dialog and collect user input (cell array of strings)
if isempty(answer), disp('User cancelled.'); return; end   % if user cancels/escapes, print message and exit

taskType         = str2double(answer{1});          % convert first answer to numeric: 1, 2, or 3
s.screenNumber   = str2double(answer{2});          % screen index for Psychtoolbox
s.EEG            = logical(str2double(answer{3})); % convert '1'/'0' to logical true/false for EEG recording
s.eyetracking    = logical(str2double(answer{4})); % convert '1'/'0' to logical for eyetracking
s.blocks         = str2double(answer{5});          % number of blocks in the experiment
s.trials_per_block = str2double(answer{6});        % number of trials per block

% runtype tag just for logging if you want
switch taskType                                   % choose descriptive runtype string based on selected task
    case 1
        s.runtype = 'COLOR3';                     % label for 3-color Stroop
    case 2
        s.runtype = 'EFST_2';                     % label for 2-emotion EFST (happy/sad)
    case 3
        s.runtype = 'EFST_3';                     % label for 3-emotion EFST (happy/neutral/sad)
    otherwise
        error('Unknown task type %d. Use 1, 2, or 3.', taskType); % error if user entered something invalid
end

if s.eyetracking                                   % if eyetracking is enabled
    s.whicheyetracker = 2;                         % choose eyetracker type (2 = Pupil Labs in your setup)
    s.checkfix        = 0;                         % disable fixation break checking (0 = off)
end

try
    %% 2) PARTICIPANT INIT, PATHS
    [p_data] = exp_init_visual('pilot','exp'); % custom function: initialize participant info, logging, diary etc.
    addpath('C:\Users\xyz\Desktop\...\....'); % ensure experiment folder is on MATLAB path 
    %% 3) BRANCH BY TASK TYPE
    switch taskType

        case 1  % ---------- 3-Color Stroop ----------
            s = setup_colorStroop_3c(s);              % call setup function for 3-color Stroop (fills timing, colors, etc.)
            [exp_seq, s] = seq_colorStroop_3c(s);     % generate trial sequence (words, ink colors, ISIs…) for color Stroop
            runfun       = @run_exp_colorStroop_3c;   % function handle to the runtime experiment function
            s.file_prefix = 'colorStroop_';           % filename prefix for saving color Stroop data/settings

        case 2  % ---------- Emotional Face Stroop – 2 emotions ----------
            s = setup_EFST_2s(s);                        % setup for EFST with 2 emotions (happy/sad faces)
            [exp_seq, s] = seq_EFST_2s(s);            % build sequence: which image, congruent/incongruent, ISI per trial
            runfun       = @run_exp_EFST_2s;             % runtime function for 2-emotion EFST
            s.file_prefix = 'faceStroop_2em_';        % filename prefix for 2-emotion EFST data/settings

        case 3  % ---------- Emotional Face Stroop – 3 emotions ----------
            s = setup_EFST_3s(s);                     % setup for EFST with 3 emotions (happy/neutral/sad)
            [exp_seq, s] = seq_EFST_3s(s);         % build 3-emotion sequence (including 3-state congruency)
            runfun       = @run_exp_EFST_3s;          % runtime function for 3-emotion EFST
            s.file_prefix = 'faceStroop_3em_';        % filename prefix for 3-emotion EFST data/settings

    end

    %% 4) SAVE SETTINGS + SEQUENCE
    save(fullfile(p_data.dir, [s.file_prefix 'settings_seq.mat']), ... % save configuration + full sequence for reproducibility
         'p_data','exp_seq','s');

    %% 5) BLOCK LOOP (COMMON TO ALL TASKS)
    for block = 1:s.blocks                              % loop over all blocks

        this_block_seq = exp_seq([exp_seq.block] == block); % extract only the trials that belong to the current block

        % run correct version for this task
        exp_data = runfun(s, p_data, this_block_seq);      % call the appropriate run function with settings + block sequence

        % compute intervals (works for all, since onsets + seq.isi exist)
        exp_data = intervals(s, exp_data);            % post-process timing intervals (fix→cue, cue→stim, stim→resp, etc.)

        % save block data
        save_exp(p_data, exp_data, s, sprintf('%02d', block)); % custom function to save per-block behavioral + timing data

        % ask if we continue, except after last block
        if block < s.blocks                               % only ask if there is another block coming
            choice = questdlg( ...                        % popup dialog asking whether to continue
                sprintf('Block %d finished. Continue with block %d?', block, block+1), ...
                'Continue?', ...
                'Yes','No','Yes');                        % 'Yes' is default button
            if ~strcmp(choice,'Yes')                      % if participant/experimenter chooses 'No'
                break;                                    % exit the block loop early
            end
        end
    end

    diary off                                            % stop any active diary logging (exp_init_visual may have enabled it)
    clear all                                            % clear all variables from workspace (hard reset for next run)
    Screen('CloseAll')                                   % close all Psychtoolbox windows
    return;                                              % exit script cleanly

catch ERR                                                % if any error occurs in the try block
    ShowCursor;                                          % make sure mouse cursor is visible again
    try KbQueueRelease; end                             % attempt to release KB queue (if it exists) without crashing on failure
    ListenChar(0);                                       % re-enable keyboard input to MATLAB command window
    Screen('CloseAll')                                   % close any open Psychtoolbox screens to avoid stuck fullscreen
    rethrow(ERR)                                         % rethrow the original error so you see the stack in MATLAB
end
