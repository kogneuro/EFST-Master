function s = setup_EFST_2s(s)
% setup_EFST_2s
% ----------------------------------------------------------
% Settings for the 2-emotion Emotional Face Stroop task.
%
% Task logic (your specs)
%  - Trial timing: fixation → cue → face+word → ISI
%      * s.fix_t  (fixation duration)
%      * s.cue_t  (cue "+")
%      * s.stim_t (face + word)
%      * ISI between s.isi_min and s.isi_max with mean ≈ s.isi_mean
%
%  - Stimuli are face images named like:
%        subjectNo_age_gender_feeling_a(or b).jpg
%        e.g. 001_y_f_h_a.jpg
%        age:    y / m / o  (young / middle / old)
%        gender: f / m
%        feeling: h / s     (happy / sad)
%
%  - Per block:
%        * equal gender (F/M)
%        * almost equal age groups (y/m/o)
%        * almost equal feelings (h/s)
%        * equal number of congruent vs incongruent trials
%        * transitions C→C, C→I, I→C, I→I roughly balanced
%
%  - Congruent: word matches face emotion (GLÜCKLICH / TRAURIG)
%  - Incongruent: word is opposite emotion.
%  - Participant responses: “GLÜCKLICH” vs “TRAURIG”.
%
% The fields here are used by:
%   - seq_EFST_2s     (sequence builder)
%   - run_exp_EFST    (runtime)
%   - intervals_EFST  (timing analysis)
%   - save_exp_EFST   (saving)
%
% Alp / 2025
% ----------------------------------------------------------

%% ---------- GENERAL ----------
% file name prefix for saving data
if ~isfield(s, 'file_prefix')
    s.file_prefix = 'faceStroop_';
end

% number of blocks in the experiment
if ~isfield(s, 'blocks')
    s.blocks = 1;                            % can be overridden from main
end

% stimulus folder (where the face jpgs are)
if ~isfield(s, 'stim_dir')
    % default: ./stimuli relative to current working directory
    s.stim_dir = fullfile(pwd, 'stimuli');
end

%% ---------- TIMING ----------
% how many flips PTB should sample to estimate refresh rate
s.get_flip_i = 200;

% NOTE: These comments were originally “1 second” but the value is 0.040.
% I leave the values as they are; if you want 1 s fixation, change to 1.000.
s.fix_t  = 0.040;     % fixation duration (in seconds)
s.cue_t  = 1.000;     % cue duration ("+" on screen)
s.stim_t = 1.000;     % face + word duration

% ISI (inter-stimulus interval) range and target mean
s.isi_min  = 2.000;   % lower bound for ISI  (seconds)
s.isi_max  = 4.000;   % upper bound for ISI  (seconds)
s.isi_mean = 3.000;   % target mean ISI over the block

% EFST-style response timing (kept for compatibility)
s.resp_window  = 'variable';   % stop checking after first press
s.resp1_max_t  = 2.000;        % max time for first response (s)
s.resp2_max_t  = 2.000;        % not really used but kept
s.resp_p_min_t = 0.300;        % minimal interval between presses (debounce logic)

%% ---------- SCREEN / VISUAL DESIGN ----------
% background color (mid-grey)
s.window_color = [128 128 128];

% convenience colors
s.white = [255 255 255];
s.black = [0   0   0];

% TEXT THAT APPEARS ON THE FACE (the Stroop word)
s.txt_color = [255 0 0];       % red word on top of face
s.txt_font  = 'Arial';

% fixation / cue
s.fix_symbol = '+';            % symbol used for fixation
s.cue        = '+';            % symbol used for cue
s.cue_size   = 60;             % font size for fixation/cue
s.cue_color  = [170 210 210];  % cue color (light salmon)

% response hint text (if used somewhere else)
s.resp1_txt  = ['H'; 'S'];     % “H” / “S” for happy/sad (legacy marker)

%% ---------- STIMULUS DISPLAY SIZE ----------
% We define face stimulus size using viewing geometry.

% physical monitor parameters (can be overwritten by main)
s.screenWidth_cm  = 60;        % monitor width in cm
s.viewDist_cm     = 70;        % viewing distance in cm
s.horizontalRes   = 1920;      % screen width in pixels
s.verticalRes     = 1080;      % screen height in pixels
s.diodeSize       = 20;        % photodiode patch size (px)

% desired stimulus size in visual degrees (face width)
s.stimSize_deg = 10;           % 10° visual angle for face width

% convert degrees → cm → pixels
cm_per_pixel = s.screenWidth_cm / s.horizontalRes;           % cm per pixel
stim_size_cm = 2 * s.viewDist_cm * tan((s.stimSize_deg/2) * pi/180);
s.stimSize_px = round(stim_size_cm / cm_per_pixel);          % face width in px

% main word size (relative to stimulus width)
s.txt_size = round(0.25 * s.stimSize_px);   % word ≈ 25% of face width
% (alternative: set a fixed value, e.g. s.txt_size = 90;)

% horizontal offset (in px) for response boxes relative to screen center
s.resp_offset = 40;

%% ---------- STIMULUS META / LOGIC ----------
% levels used by the sequence builder

s.feelings    = {'h','s'};           % happy / sad
s.age_groups  = {'y','m','o'};       % young / middle / old
s.genders     = {'f','m'};           % female / male

% congruency labels
s.cong_labels = {'congruent','incongruent'};

% transition patterns you care about (prev → current)
%   C → C
%   C → I
%   I → C
%   I → I
s.patterns = {
    'C','C';
    'C','I';
    'I','C';
    'I','I';
};

% number of times each pattern should appear per block:
%    pattern_repeats = trials_per_block / 4
% (trials_per_block is defined in main, not here)
s.pattern_repeats = s.trials_per_block / 4;

%% ---------- RESPONSES (HAPPY / SAD) ----------
% labels used on screen and for logging
s.resp_labels = {'GLÜCKLICH','TRAURIG'};    % happy / sad (German)

% Button box / logical button indices (EFST legacy)
s.lpt_adr1       = '4FF8';     % parallel port address (string form)
s.lpt_adr2       = '4FF4';     % typically adr1 + 4 in hex
s.lpt_dir        = 'bi';       % bidirectional
s.debounceDelay  = 0.025;      % software debounce for buttons (s)

% logical button assignments
s.btn_resp1 = [1 2];           % valid buttons for response 1
s.btn_resp2 = [1 2];           % valid buttons for response 2
s.btn_esc   = 3;               % button index for ESC / quit

% physical keyboard keys (for debug or if no button box)
% mapping: first = “GLÜCKLICH”, second = “TRAURIG”
s.responseKeys = {'LeftArrow', 'RightArrow'};

%% ---------- EEG / TRIGGERS ----------
% keep trigger structure compatible with previous EFST / respirationCA toolbox

s.use_eeg_trigger = 1;         % 1 = send triggers, 0 = no EEG triggers

% parallel port / LPT config (bit-banging triggers)
s.ioObj    = io64;             % I/O object handle (from io64 library)
s.address  = hex2dec('DFB8');  % base address for port (PC-specific)
s.lpt_adr1 = '4FF8';           % data port
s.lpt_adr2 = '4FF4';           % status/control port
s.lpt_dir  = 'bi';             % bidirectional
s.TTL_V    = 5;                % TTL voltage (for documentation)

%% ---------- INSTRUCTIONS (optional) ----------
% If you use image-based instructions, these fields tell the
% instruction function where to look for them.

s.instr_dir             = '/instr/';          % relative to current script path
s.instr_subdir_wildcard = 'condition_*';      % which subfolder to use
s.instr_img_wildcard    = 'step*.jpg';        % which images inside that folder

end
