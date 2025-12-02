function s = setup_EFST_3s(s)
% setup_EFST_3s – settings for Emotional Face Stroop (3 emotions: H / N / S)
%
% This function prepares a configuration struct "s" that controls:
%   • general experiment parameters (blocks, trials, stim folder)
%   • timing (fixation, cue, stimulus, ISI)
%   • screen / geometry settings (visual angle → pixels)
%   • stimulus logic (feelings, congruency states, transition patterns)
%   • response mapping (3 keys, 3 labels)
%   • EEG trigger / parallel port settings
%   • instruction image folder
%
% It mirrors your 2-state EFST setup, but now supports:
%   • 3 feelings in filenames: h (happy), n (neutral), s (sad)
%   • 3 congruency codes: C, I, IS
%   • 9 possible transitions between consecutive trials.
%
% Input:
%   s – (optional) struct; any fields not present will be filled by defaults.
%
% Output:
%   s – fully initialized setup struct, ready for seq_EFST_3s and run_exp_EFST_3s.
%
% Alp — 2025
% -------------------------------------------------------------------------

%% ---------- GENERAL ----------
% file name prefix for saved data files
if ~isfield(s, 'file_prefix')
    s.file_prefix = 'faceStroop_';
end

% number of blocks in the experiment
if ~isfield(s, 'blocks')
    s.blocks = 1;
end

% number of trials per block
if ~isfield(s, 'trials_per_block')
    s.trials_per_block = 40;
end

% folder that contains the face images (e.g. 001_y_f_h.jpg etc.)
if ~isfield(s, 'stim_dir')
    s.stim_dir = fullfile(pwd, 'stimuli');  % default: ./stimuli relative to current folder
end

%% ---------- TIMING ----------
% how many flips GetFlipInterval uses to estimate refresh rate
s.get_flip_i = 200;

% core timing (in seconds)
s.fix_t  = 0.040;  % fixation duration (you kept 40 ms as in your other code)
s.cue_t  = 1.000;  % cue duration
s.stim_t = 1.000;  % stimulus (face + word) duration

% ISI range and target mean (in seconds)
% total ISI is drawn in [isi_min, isi_max] and adjusted to have mean ≈ isi_mean
s.isi_min  = 2.000;
s.isi_max  = 4.000;
s.isi_mean = 3.000;

% Response timing:
%   resp_window = 'variable' → first keypress ends response interval
s.resp_window  = 'variable';
s.resp1_max_t  = 2.000;   % max time allowed to respond to first choice
s.resp2_max_t  = 2.000;   % kept for compatibility, not used here
s.resp_p_min_t = 0.300;   % minimal time between presses (if you ever use 2 responses)

%% ---------- SCREEN / VISUAL ----------
% background color of the PTB window (mid-grey)
s.window_color = [128 128 128];

% some convenience colors
s.white = [255 255 255];
s.black = [0   0   0];

% color and font for Stroop words drawn on the face
s.txt_color = [255 0 0];   % red text
s.txt_font  = 'Arial';

% fixation and cue symbols (same character, different color)
s.fix_symbol = '+';
s.cue        = '+';
s.cue_size   = 60;         % large font size for fixation/cue
s.cue_color  = [170 210 210];

% simple label vector (not used for drawing, but kept for compatibility)
% H / N / S = Happy / Neutral / Sad
s.resp1_txt = ['H'; 'N'; 'S'];   % columns of characters (legacy style)

%% ---------- STIMULUS SIZE ----------
% Physical screen + viewing geometry (can be overwritten in main script if needed)
s.screenWidth_cm  = 60;     % physical screen width in cm
s.viewDist_cm     = 70;     % viewing distance in cm
s.horizontalRes   = 1920;   % horizontal resolution in pixels
s.verticalRes     = 1080;   % vertical resolution in pixels
s.diodeSize       = 20;     % photodiode rectangle size in pixels

% Desired stimulus size in visual degrees (face image width)
s.stimSize_deg = 10;        % ~10° visual angle is a reasonable default

% Convert desired visual angle → cm → pixels
cm_per_pixel = s.screenWidth_cm / s.horizontalRes;                     % cm for one pixel
stim_size_cm = 2 * s.viewDist_cm * tan((s.stimSize_deg/2) * pi/180);   % chord length in cm
s.stimSize_px = round(stim_size_cm / cm_per_pixel);                    % final width in px

% Base text size for on-face word, relative to stimulus width
s.txt_size    = round(0.25 * s.stimSize_px);   % ≈ 25% of image width
s.resp_offset = 40;                             % spacing used in some layouts

%% ---------- STIMULUS META / LOGIC ----------
% The parser in seq_EFST_3s expects:
%   feeling: h / n / s
%   age:     y / m / o
%   gender:  f / m
s.feelings    = {'h','n','s'};   % happy, neutral, sad
s.age_groups  = {'y','m','o'};   % young, middle, old
s.genders     = {'f','m'};       % female, male

% Congruency labels to be stored in the sequence struct
% (these correspond to codes: 'C', 'I', 'IS')
s.cong_labels = {'congruent','incongruent','slightlyIncongruent'};

% All 9 possible congruency transitions (prevCong → currCong)
% C  = congruent
% I  = incongruent
% IS = slightly incongruent
s.patterns = {
    'C','C';
    'C','I';
    'C','IS';
    'I','C';
    'I','I';
    'I','IS';
    'IS','C';
    'IS','I';
    'IS','IS';
    };

% In the 2-state version we could fix e.g. 40 trials / 4 patterns = 10 each.
% Here we have 9 patterns; 40 is not divisible by 9.
% So we do NOT enforce exact pattern_repeats here – seq_EFST_3s will
% build a "best possible" 3-state sequence internally.
s.pattern_repeats = [];  % left empty on purpose – handled in seq_EFST_3s

%% ---------- RESPONSES (3 choices: H / N / S) ----------
% Labels that will be drawn above/below the response boxes:
s.resp_labels = {'GLÜCKLICH','NEUTRAL','TRAURIG'};

% Keyboard mapping:
%   RightArrow = Happy
%   UpArrow    = Neutral
%   LeftArrow  = Sad
% (run_exp_EFST_3s uses these to decode KbName into response codes 1/2/3)
s.responseKeys = {'LeftArrow','UpArrow','RightArrow'};

% Button codes (for parallel button box or logging)
% 1, 2, 3 correspond to the three response options;
% 4 is reserved as ESC / abort code.
s.debounceDelay = 0.025;       % small delay to avoid bouncing
s.btn_resp1     = [1 2 3];
s.btn_resp2     = [1 2 3];
s.btn_esc       = 4;

%% ---------- EEG / TRIGGERS ----------
% Keep same structure as in your other EFST / respiration experiments
s.use_eeg_trigger = 1;

% Parallel port initialisation (io64 from IOPort / IO32 driver)
s.ioObj    = io64;
s.address  = hex2dec('DFB8');  % base LPT address – must be adjusted per PC
s.lpt_adr1 = '4FF8';
s.lpt_adr2 = '4FF4';
s.lpt_dir  = 'bi';             % bidirectional mode
s.TTL_V    = 5;                % nominal TTL voltage

%% ---------- INSTRUCTIONS ----------
% Optional: you can show block-specific instruction images
s.instr_dir             = '/instr/';          % base folder (relative to this .m)
s.instr_subdir_wildcard = 'condition_*';      % which subfolder to use
s.instr_img_wildcard    = 'step*.jpg';        % instruction image pattern

end
