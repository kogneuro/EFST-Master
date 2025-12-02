function s = setup_colorStroop_3c(s)
% setup_colorStroop_3c
% ----------------------------------------------------------
% Configure all parameters for the 3-color Stroop task:
%   • 3 color words: RED, GREEN, BLACK
%   • 2 congruency states: congruent (C), incongruent (I)
%   • Response as color naming (not reading)
%
% This function ONLY:
%   • Fills / updates fields in the parameter struct "s"
%   • Defines timing, visual properties, color words,
%     response mapping, and trigger settings
%
% This function does NOT:
%   • Initialize participant data (handled by exp_init_visual)
%   • Generate trial sequences (done in seq_colorStroop_3c)
%   • Run the experiment (done in run_exp_colorStroop_3c)
%   • Add paths or open PTB windows
%
% Input:
%   s  – struct (can be empty or partially filled)
%
% Output:
%   s  – struct with all required fields for the 3-color Stroop task
%
% Author: Alp, 2025 (based on EFST framework)
% ----------------------------------------------------------

%% ---------- GENERAL PARAMETERS ----------

% File prefix for saving settings and data
if ~isfield(s, 'file_prefix')
    s.file_prefix = 'colorStroop_';
end

% Number of blocks, if not already set from GUI / main script
if ~isfield(s, 'blocks')
    s.blocks = 1;
end

% Number of trials per block (default: 60)
if ~isfield(s, 'trials_per_block')
    s.trials_per_block = 60;
end


%% ---------- TIMING PARAMETERS ----------

% Number of samples for GetFlipInterval averaging
s.get_flip_i = 200;

% Durations (in seconds)
s.fix_t  = 0.040;   % fixation cross
s.cue_t  = 1.000;   % cue duration
s.stim_t = 1.000;   % Stroop stimulus duration

% Inter-stimulus interval range (for seq function to fill)
s.isi_min  = 2.000;
s.isi_max  = 4.000;
s.isi_mean = 3.000;


%% ---------- SCREEN / VISUAL APPEARANCE ----------

% Background and basic colors (RGB)
s.window_color = [128 128 128];   % mid-gray background
s.white        = [255 255 255];
s.black        = [  0   0   0];

% Default text color for fixation and general text
s.txt_color = [255 255 255];      % white by default
s.txt_font  = 'Arial';

% Fixation & cue symbols
s.fix_symbol = '+';
s.cue        = '+';
s.cue_size   = 60;
s.cue_color  = [170 210 210];     % light cyan/blue for the cue


%% ---------- STIMULUS SIZE / TEXT SCALING ----------

% Physical display parameters (for approximate visual angle)
s.screenWidth_cm  = 60;    % monitor width in cm
s.viewDist_cm     = 70;    % viewing distance in cm
s.horizontalRes   = 1920;  % resolution width in pixels
s.verticalRes     = 1080;  % resolution height in pixels
s.diodeSize       = 20;    % photodiode square size in pixels

% Desired stimulus size in degrees of visual angle
s.stimSize_deg = 10;

% Convert degrees → cm → pixels
cm_per_pixel   = s.screenWidth_cm / s.horizontalRes;
stim_size_cm   = 2 * s.viewDist_cm * tan((s.stimSize_deg/2) * pi/180);
s.stimSize_px  = round(stim_size_cm / cm_per_pixel);

% Font size of Stroop word = 20% of “stimulus size”
s.txt_size     = round(0.20 * s.stimSize_px);

% Horizontal spacing used for response boxes (left / middle / right)
s.resp_offset  = 40;


%% ---------- STIMULUS META / LOGIC (COLOR STROOP) ----------

% 3 word identities used as Stroop stimuli
s.color_words = {'RED','GREEN','BLACK'};

% Define ink colors (RGB values)
s.colors.RED   = [255   0   0];   % red ink
s.colors.GREEN = [  0 255   0];   % green ink
s.colors.BLACK = [  0   0   0];   % black ink

% Labels shown for the response options (same words)
s.resp_labels = {'RED','GREEN','BLACK'};

% Response keys:
%   • LeftArrow  → RED
%   • UpArrow    → GREEN
%   • RightArrow → BLACK
%
% (You can change these if you prefer other keys, but then make sure
%  run_exp_colorStroop_3c uses the same mapping for evaluation.)
s.responseKeys = {'LeftArrow','UpArrow','RightArrow'};

% Only two congruency types here (no IS condition in classic Stroop):
s.cong_labels = {'congruent','incongruent'};

% Transition patterns (previous → current congruency)
% These are used in seq_colorStroop_3c to build a balanced C/I path.
s.patterns = {
    'C','C';   % congruent → congruent
    'C','I';   % congruent → incongruent
    'I','C';   % incongruent → congruent
    'I','I';   % incongruent → incongruent
};
% NOTE:
%   The exact balancing of these transitions is implemented
%   in build_cong_sequence_Color inside seq_colorStroop_3c.


%% ---------- RESPONSE PARAMETERS ----------

% Response window style (we keep it like EFST)
s.resp_window  = 'variable';

% Maximum allowed response times (in seconds)
s.resp1_max_t  = 3.000;
s.resp2_max_t  = 3.000;

% Minimal latency (for debouncing, used in some variants)
s.resp_p_min_t = 0.300;

% Button indices for external response box compatibility
s.debounceDelay = 0.025;
s.btn_resp1     = [1 2 3];
s.btn_resp2     = [1 2 3];
s.btn_esc       = 4;


%% ---------- EEG / TRIGGER SETTINGS ----------

% Whether we intend to send EEG triggers (actual sending is done elsewhere)
s.use_eeg_trigger = 1;

% IO64 parallel port object and address
% (adapt these to your actual PC / LPT configuration if needed)
s.ioObj    = io64;               % handle from IO64 library
s.address  = hex2dec('DFB8');    % base address (example)
s.lpt_adr1 = '4FF8';
s.lpt_adr2 = '4FF4';
s.lpt_dir  = 'bi';
s.TTL_V    = 5;                  % TTL voltage level (5 V)


%% ---------- INSTRUCTION SLIDES (OPTIONAL) ----------

% If you use picture-based instructions, these fields help locate them.
s.instr_dir             = '/instr/';          % base instruction folder (relative)
s.instr_subdir_wildcard = 'condition_*';      % subfolder pattern
s.instr_img_wildcard    = 'step*.jpg';        % instruction image pattern

end  % ---- end of setup_colorStroop_3c ----
