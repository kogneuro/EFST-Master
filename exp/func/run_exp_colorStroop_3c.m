function data = run_exp_colorStroop_3c(s, p_data, seq_block)
% run_exp_colorStroop_3c
% ----------------------------------------------------------
% Executes a full block of the 3-color Stroop task 
% (RED / GREEN / BLACK) within the respirationCA experimental
% framework.
%
% This version follows the same structural conventions as the EFST
% run-functions:
%   - exp_init_visual.m for participant + diary setup
%   - intervals_EFST.m for timing diagnostics
%   - seq_colorStroop_3c.m for trial sequence generation
%
% Functionality:
%   • 3-response Stroop (RED / GREEN / BLACK)
%   • Congruent vs. incongruent mapping with balanced transitions
%   • Adjustable ISIs (3–5 sec, mean ~4)
%   • Dynamic word-ink assignment per trial
%   • Response screen embedded within ISI duration
%   • Optional EEG TTL pulses (BrainVision)
%   • Optional Pupil Labs eyetracking annotations
%
% Inputs:
%   s         – settings struct from setup_colorStroop_3c
%   p_data    – participant struct (ID, folder, diary)
%   seq_block – struct array describing each trial’s:
%                   .word
%                   .wordIdx
%                   .inkColor (RGB)
%                   .inkColorIdx
%                   .inkName
%                   .prevCong
%                   .currCong
%                   .congruency
%                   .isi
%
% Output:
%   data – struct containing:
%              .rng_state
%              .seq
%              .flip_t
%              .stimFrames
%              .onset_fix / cue / stim / resp1
%              .resp1, .resp1_btn, .resp1_t
%              .btn_instr / onset_instr (if used)
%
% Notes:
%   • Timing uses the same psychophysics-safe logic as EFST.
%   • Response collection uses KbQueue for millisecond accuracy.
%   • Response screen is shown immediately after stimulus offset
%     and remains visible for min(ISI, s.resp1_max_t).
%   • Remaining ISI is filled with blank screen (baseline).
%
% Alp, 2025
% ----------------------------------------------------------

try
    %% 1) EYETRACKING (optional, Pupil Labs)
    use_pupil = false;                                    % default: no Pupil recording

    % Only if eyetracking is enabled in settings
    if isfield(s,'eyetracking') && s.eyetracking == 1

        % We only set up Pupil Labs if whicheyetracker == 2
        if isfield(s,'whicheyetracker') && s.whicheyetracker == 2

            endpoint =  'tcp://000.000.00.00:00000';      % Pupil remote endpoint (REQ socket)
            ctx = zmq.core.ctx_new();                     % create new ZeroMQ context
            socket = zmq.core.socket(ctx, 'ZMQ_REQ');     % request socket for commands
            zmq.core.setsockopt(socket, 'ZMQ_RCVTIMEO', 1000); % 1s timeout for replies
            fprintf('Connecting to %s\n', endpoint);
            zmq.core.connect(socket, endpoint);           % connect to Pupil

            % Ask for SUB_PORT to subscribe to data
            zmq.core.send(socket, uint8('SUB_PORT'));
            sub_port = char(zmq.core.recv(socket));
            fprintf('Received sub port: %s\n', sub_port);

            % If we get false, something went wrong
            if isequal(sub_port, false)
                warning('No valid sub port received');
                return;
            end

            % Build SUB endpoint and connect subscriber socket
            ip_address   = '000.000.00.00';
            sub_endpoint = sprintf('tcp://%s:%s', ip_address, sub_port);
            sub_socket   = zmq.core.socket(ctx, 'ZMQ_SUB');
            zmq.core.connect(sub_socket, sub_endpoint);
            zmq.core.setsockopt(sub_socket, 'ZMQ_RCVTIMEO', 1000);
            fprintf('Connecting to SUB: %s\n', sub_endpoint);

            % Time sync: measure round-trip delay
            tic;
            zmq.core.send(socket, uint8('t'));            % 't' query
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));
            fprintf('Round trip command delay: %s\n', toc);

            % Reset Pupil time to 0.0
            zmq.core.send(socket, uint8('T 0.0'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));

            % Start recording
            pause(1.0);
            zmq.core.send(socket, uint8('R'));
            result = zmq.core.recv(socket);
            fprintf('Recording should start: %s\n', char(result));

            use_pupil = true;                             % flag that we have Pupil running
        else
            use_pupil = false;                            % other ET systems not handled here
        end
    else
        use_pupil = false;                                % eyetracking disabled
    end

    %% 2) RNG + sequence
    data.rng_state = set_rng_state(s);                    % store RNG state for reproducibility
    data.seq       = seq_block;                           % store the block sequence
    nTrials        = numel(seq_block);                    % number of trials in this block

    %% 3) KEYBOARD (3 response keys + ESC)
    KbName('UnifyKeyNames');                              % unify key names across OSes
    ListenChar(-1);                                       % suppress keypresses in MATLAB command window

    keylist = zeros(1,256);                              % 256-length mask for KbQueue
    keylist(KbName(s.responseKeys{1})) = 1;              % enable key for RED
    keylist(KbName(s.responseKeys{2})) = 1;              % enable key for GREEN
    keylist(KbName(s.responseKeys{3})) = 1;              % enable key for BLACK
    keylist(KbName('ESCAPE'))          = 1;              % ESC for abort

    KbQueueCreate(-1, keylist);                          % create queue on default device
    KbQueueStart;                                        % start listening

    %% 4) SCREEN / Psychtoolbox window
    AssertOpenGL;                                        % make sure PTB is working
    PsychDefaultSetup(2);                                % typical PTB defaults

    % Open onscreen window with background color s.window_color
    [window, windowRect] = Screen('OpenWindow', s.screenNumber, s.window_color);
    Screen('ColorRange', window, 255);                   % full 0–255 range
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;                                          % hide mouse cursor

    % Screen center coordinates
    [xCenter, yCenter] = RectCenter(windowRect);

    % Photodiode patch rectangle (e.g. top center)
    cueRect = [0 0 s.diodeSize s.diodeSize];             % base rect size
    cueRect = CenterRectOnPointd(cueRect, xCenter, s.diodeSize/2 + 10);

    % Set text properties
    Screen('TextFont', window, s.txt_font);
    Screen('TextSize', window, s.txt_size);

    % Measure flip interval (refresh rate) using PTB's routine
    Priority(1);
    flip_t = Screen('GetFlipInterval', window, s.get_flip_i);
    Priority(0);
    data.flip_t = flip_t;                                % store flip interval

    % Number of frames for stimulus based on desired duration
    data.stimFrames = round(s.stim_t / data.flip_t);

    % Get screen size in pixels
    [screenX, screenY] = Screen('WindowSize', window);

    % y-position row for the 3 response boxes (middle of screen vertically)
    yPos = screenY * 0.50;

    % Define width/height of the response rectangles
    w = 300;
    h = 80;

    % Center three rectangles at 1/6, 3/6, 5/6 of screen width
    rect1 = CenterRectOnPointd([0 0 w h], screenX*(1/6), yPos);  % left box
    rect2 = CenterRectOnPointd([0 0 w h], screenX*(3/6), yPos);  % middle box
    rect3 = CenterRectOnPointd([0 0 w h], screenX*(5/6), yPos);  % right box

    % Initial screen flip (clears backbuffer)
    Screen('Flip', window);

    %% 5) TIMING ARRAYS (pre-allocate)
    [data.onset_fix, ...
     data.onset_cue, ...
     data.onset_stim, ...
     data.onset_delay, ...
     data.onset_resp1, ...
     data.onset_resp1_p] = deal(cell(nTrials,5));

    % Round fixation duration to integer multiples of flip_t
    data.fix_t     = s.fix_t - mod(s.fix_t, flip_t);
    data.fixFrames = round(data.fix_t / flip_t);

    % Round cue duration similarly
    data.cue_t     = s.cue_t - mod(s.cue_t, flip_t);
    data.cueFrames = round(data.cue_t / flip_t);

    % Placeholders for instruction button + onsets (if used)
    [data.btn_instr, data.onset_instr] = deal(cell(1,5));

    % (Optional) 6) INSTRUCTIONS
    if seq_block(1).block == 1
        instr_dir   = [fileparts(mfilename('fullpath')) s.instr_dir];
        instr_subdir = dir([instr_dir s.instr_subdir_wildcard '1*']);
        if ~isempty(instr_subdir)
            instr_images = load_images([instr_dir instr_subdir.name '\'], s.instr_img_wildcard);
            [data.btn_instr, data.onset_instr] = show_instr_img_EFST(instr_images, window);
            clear instr_images
        end
    end

    %% 7) TRIAL LOOP
    for i = 1:nTrials
        tr = seq_block(i);                               % current trial struct

        % ---------- FIXATION ----------
        Screen('TextSize', window, s.cue_size);          % use cue size for fixation
        DrawFormattedText(window, s.fix_symbol, 'center','center', s.txt_color);
        fix_on = Screen('Flip', window);                 % flip to show fixation
        data.onset_fix{i,1} = fix_on;                    % log onset

        % EEG trigger: trial start
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 30);          % code 30 = trial start
            triggerout  = 1;                             %#ok<NASGU> (kept if you use it elsewhere)
            triggertime = GetSecs;                       %#ok<NASGU>
        end

        % Pupil Labs annotation: trial start
        if use_pupil
            zmq.core.send(socket, uint8('t'));           % request current Pupil time
            currentTime_bytes = zmq.core.recv(socket,20);
            currentTimeNum  = str2double(char(currentTime_bytes));

            keys_start = {'topic','label','timestamp','duration'};
            vals_start = {'annotation.TrialStart', ...
                          sprintf('Trial #%d started',i), ...
                          currentTimeNum, ...
                          1.0};

            start_annotation = containers.Map(keys_start, vals_start);
            send_annotation(socket, start_annotation);
            result = zmq.core.recv(socket);              %#ok<NASGU>
        end

        % ---------- CUE ----------
        Screen('TextSize', window, s.cue_size);          % set font size for cue
        DrawFormattedText(window, s.cue, 'center','center', s.cue_color);

        % schedule cue onset after fixation duration (minus half a frame)
        cue_on = Screen('Flip', window, fix_on + (data.fixFrames - 0.5)*flip_t);
        data.onset_cue{i,1} = cue_on;                    % log cue onset

        % EEG trigger: cue onset
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 11);          % code 11 = cue onset
            triggerout  = 1;                             %#ok<NASGU>
            triggertime = GetSecs;                       %#ok<NASGU>
        end

        % ---------- STIMULUS (COLORED WORD) ----------
        Screen('TextSize', window, s.txt_size);          % set font size for Stroop word

        % Draw the Stroop word in its ink color at screen center
        DrawFormattedText(window, tr.word, 'center', 'center', tr.inkColor);

        % Draw photodiode patch (white square)
        Screen('FillRect', window, s.white, cueRect);

        % Flip: show stimulus after cue duration
        stim_on = Screen('Flip', window, cue_on + (data.cueFrames - 0.5)*flip_t);
        data.onset_stim{i,1} = stim_on;                  % log stimulus onset

        % EEG trigger: stimulus onset
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 12);          % code 12 = stimulus
            triggerout  = 1;                             %#ok<NASGU>
            triggertime = GetSecs;                       %#ok<NASGU>
        end

        % Keep stimulus on screen for s.stim_t seconds
        WaitSecs(s.stim_t);
        KbQueueFlush;                                    % clear any keypresses during stim

        % ---------- RESPONSE SCREEN (within ISI) ----------
        isi_total = tr.isi;                              % planned ISI for this trial

        Screen('TextSize', window, s.txt_size);          % text size for response labels

        % Draw three response labels centered in rect1/rect2/rect3
        DrawFormattedText(window, s.resp_labels{1}, ...
                          'center','center', s.white, [], [], [], [], [], rect1);
        DrawFormattedText(window, s.resp_labels{2}, ...
                          'center','center', s.white, [], [], [], [], [], rect2);
        DrawFormattedText(window, s.resp_labels{3}, ...
                          'center','center', s.white, [], [], [], [], [], rect3);

        % Flip to show response screen
        resp_on = Screen('Flip', window);
        data.onset_resp1{i,1} = resp_on;                 % log response screen onset

        KbQueueFlush;                                    % clear buffer at start of response window

        % EEG trigger: response screen onset
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 13);          % code 13 = response screen
            triggerout  = 1;                             %#ok<NASGU>
            triggertime = GetSecs;                       %#ok<NASGU>
        end

        % Response screen is visible for min(resp1_max_t, ISI)
        resp_screen_window = min(s.resp1_max_t, isi_total);

        resp_btn = 0;                                    % 1=RED, 2=GREEN, 3=BLACK, 0=none
        rt       = NaN;                                  % reaction time
        start_t  = GetSecs;                              % start of response window

        % Collect first key within response window
        while GetSecs - start_t < resp_screen_window
            [pressed, firstPress] = KbQueueCheck;
            if pressed
                k     = find(firstPress,1);              % first key index
                kname = KbName(k);                       % key name string
                rt    = GetSecs - start_t;               % RT relative to response onset

                % Map key to response button index
                if strcmp(kname, s.responseKeys{1})
                    resp_btn = 1;                        % RED
                elseif strcmp(kname, s.responseKeys{2})
                    resp_btn = 2;                        % GREEN
                elseif strcmp(kname, s.responseKeys{3})
                    resp_btn = 3;                        % BLACK
                elseif strcmpi(kname, 'ESCAPE')
                    resp_btn = s.btn_esc;                % ESC pressed
                end
                break;
            end
        end

        % Store raw response button + RT
        data.resp1_btn(i,1) = resp_btn;
        data.resp1_t(i,1)   = rt;

        % Convert resp_btn into categorical code (1/2/3 or -1)
        if resp_btn == 1
            data.resp1(i,1) = 1;                         % RED
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 41);      % EEG code for RED response
                triggerout  = 1;                         %#ok<NASGU>
                triggertime = GetSecs;                   %#ok<NASGU>
            end

        elseif resp_btn == 2
            data.resp1(i,1) = 2;                         % GREEN
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 42);      % EEG code for GREEN response
                triggerout  = 1;                         %#ok<NASGU>
                triggertime = GetSecs;                   %#ok<NASGU>
            end

        elseif resp_btn == 3
            data.resp1(i,1) = 3;                         % BLACK
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 43);      % EEG code for BLACK response
                triggerout  = 1;                         %#ok<NASGU>
                triggertime = GetSecs;                   %#ok<NASGU>
            end

        elseif resp_btn == s.btn_esc
            % ESC: break out of trial loop (stop block)
            break;

        else
            % No valid response within window
            data.resp1(i,1) = -1;
        end

        % ---------- Fill remaining ISI (blank screen) ----------
        elapsed_resp  = GetSecs - resp_on;               % time spent on response screen
        remaining_ISI = isi_total - elapsed_resp;        % leftover ISI

        if remaining_ISI > 0
            Screen('FillRect', window, s.window_color);  % blank background
            Screen('Flip', window);                      % flip to blank
            WaitSecs(remaining_ISI);                     % wait until ISI is complete
        end
    end % trial loop

    %% END / SAVE
    sca;                                                 % close PTB window (Screen('CloseAll'))
    ListenChar(0);                                       % restore keyboard input to MATLAB
    KbQueueRelease;                                      % free keyboard queue

    % Save block data to participant folder, if dir exists
    if isfield(p_data,'dir')
        save(fullfile(p_data.dir, ...
            [s.file_prefix sprintf('block_%02d_%s.mat', seq_block(1).block, p_data.ID)]), ...
            'data','s','seq_block');
    end

catch ERR
    % On any error: cleanly close PTB & keyboard and rethrow
    sca;
    ListenChar(0);
    KbQueueRelease;
    rethrow(ERR);
end

end % function run_exp_colorStroop_3c
