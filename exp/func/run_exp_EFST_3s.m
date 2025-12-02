function data = run_exp_EFST_3s(s, p_data, seq_block)
% ============================================================
% run_exp_EFST_3s
% Emotional Face Stroop task – 3-response version
%
% Trial structure:
%   1) fixation cross:           0.04 s
%   2) cue (e.g. +):             1 s
%   3) stimulus (face + word):   1 s
%   4) response screen (3 labels: Happy / Neutral / Sad)
%        - appears immediately after stimulus
%        - stays on screen for min(ISI, s.resp1_max_t)
%   5) if the total ISI for this trial is longer than the response period,
%      the remaining time is filled with a blank screen.
%
% Requirements in s:
%   s.responseKeys = {'LeftArrow','DownArrow','RightArrow'};
%   s.resp_labels  = {'GLÜCKLICH','NEUTRAL','TRAURIG'};
%
% Input:
%   s         – setup struct (from setup_EFST_3s)
%   p_data    – participant data struct (fields: ID, dir, etc.)
%   seq_block – one block of trials (from seq_EFST_3s), struct array with
%               fields like img, word, isi, prevCong, currCong, ...
%
% Output:
%   data – struct with behavioral responses and timing info
%
% Alp — 2025
% ============================================================

try
    %% --------------------------------------------------------
    % 1) EYETRACKING INITIALIZATION (Pupil Labs, optional)
    %% --------------------------------------------------------
    if isfield(s,'eyetracking') && s.eyetracking == 1
        % check which eyetracker to use
        if isfield(s,'whicheyetracker') && s.whicheyetracker == 2
            % ---- Pupil Labs via ZMQ ----
            endpoint =  'tcp://000.000.00.00';   % Pupil Remote address

            % create ZMQ context and request socket
            ctx    = zmq.core.ctx_new();
            socket = zmq.core.socket(ctx, 'ZMQ_REQ');

            % set receive timeout to avoid blocking if server is not reachable
            zmq.core.setsockopt(socket, 'ZMQ_RCVTIMEO', 1000);

            fprintf('Connecting to %s\n', endpoint);
            zmq.core.connect(socket, endpoint);

            % request SUB port used for streaming data
            zmq.core.send(socket, uint8('SUB_PORT'));
            sub_port = char(zmq.core.recv(socket));
            fprintf('Received sub port: %s\n', sub_port);

            if isequal(sub_port, false)
                warning('No valid sub port received from Pupil.');
                return;  % exit experiment if eyetracking fails
            end

            % set up SUB socket for data stream (if needed)
            ip_address  = '000.000.00.00';
            sub_endpoint = sprintf('tcp://%s:%s', ip_address, sub_port);
            sub_socket   = zmq.core.socket(ctx, 'ZMQ_SUB');
            zmq.core.connect(sub_socket, sub_endpoint);
            zmq.core.setsockopt(sub_socket, 'ZMQ_RCVTIMEO', 1000);
            fprintf('Connecting to SUB: %s\n', sub_endpoint);

            % synchronize time with Pupil
            tic;
            zmq.core.send(socket, uint8('t'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));
            fprintf('Round trip command delay: %s\n', toc);

            % reset Pupil time to 0
            zmq.core.send(socket, uint8('T 0.0'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));

            % start Pupil recording
            pause(1.0);
            zmq.core.send(socket, uint8('R'));
            result = zmq.core.recv(socket);
            fprintf('Recording should start: %s\n', char(result));

            use_pupil = true;
        else
            % eyetracking flag is ON but not Pupil Labs mode
            use_pupil = false;
        end
    else
        % eyetracking disabled
        use_pupil = false;
    end

    %% --------------------------------------------------------
    % 2) RNG + ASSIGN SEQUENCE
    %% --------------------------------------------------------
    % store current RNG state to allow reproducing the same randomization
    data.rng_state = set_rng_state(s);

    % store the sequence for this block
    data.seq = seq_block;
    nTrials  = numel(seq_block);

    %% --------------------------------------------------------
    % 3) KEYBOARD INITIALIZATION (3 response keys)
    %% --------------------------------------------------------
    KbName('UnifyKeyNames');   % standardize key names across OS
    ListenChar(-1);            % prevent keypresses from leaking to MATLAB

    % build a key mask: 1 for keys we care about, 0 otherwise
    keylist = zeros(1,256);
    keylist(KbName(s.responseKeys{1})) = 1;   % Happy
    keylist(KbName(s.responseKeys{2})) = 1;   % Neutral
    keylist(KbName(s.responseKeys{3})) = 1;   % Sad
    keylist(KbName('ESCAPE'))          = 1;   % emergency quit

    KbQueueCreate(-1, keylist);   % use default keyboard (-1)
    KbQueueStart;                 % start key queue

    %% --------------------------------------------------------
    % 4) SCREEN SETUP
    %% --------------------------------------------------------
    AssertOpenGL;
    PsychDefaultSetup(2);

    % open window on specified screen with given background color
    [window, windowRect] = Screen('OpenWindow', s.screenNumber, s.window_color);
    Screen('ColorRange', window, 255);
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    HideCursor;

    % center coordinates of the window
    [xCenter, yCenter] = RectCenter(windowRect);

    % photodiode rectangle (small patch near top center)
    cueRect = [0 0 s.diodeSize s.diodeSize];
    cueRect = CenterRectOnPointd(cueRect, xCenter, s.diodeSize/2 + 10);

    % set base font properties
    Screen('TextFont', window, s.txt_font);
    Screen('TextSize', window, s.txt_size);

    % measure flip interval (frame duration)
    Priority(1);
    flip_t = Screen('GetFlipInterval', window, s.get_flip_i);
    Priority(0);
    data.flip_t = flip_t;

    % number of frames corresponding to stimulus duration
    data.stimFrames = round(s.stim_t / data.flip_t);

    % ---------------------------------------------------------------------
    % Define positions of the 3 response boxes
    % ---------------------------------------------------------------------
    [screenX, screenY] = Screen('WindowSize', window);

    % vertical position (middle of the screen)
    yPos = screenY * 0.50;

    % box size
    w = 300;
    h = 80;

    % three horizontally spaced rectangles at 1/6, 3/6, 5/6 of width
    rect1 = CenterRectOnPointd([0 0 w h], screenX*(1/6), yPos);  % left option
    rect2 = CenterRectOnPointd([0 0 w h], screenX*(3/6), yPos);  % middle option
    rect3 = CenterRectOnPointd([0 0 w h], screenX*(5/6), yPos);  % right option

    % initial buffer flip (clean gray screen)
    Screen('Flip', window);

    %% --------------------------------------------------------
    % 5) TIMING / ONSET ARRAYS
    %% --------------------------------------------------------
    % onsets are stored as cell arrays to stay compatible with older code
    [data.onset_fix, ...
     data.onset_cue, ...
     data.onset_stim, ...
     data.onset_delay, ...
     data.onset_resp1, ...
     data.onset_resp1_p] = deal(cell(nTrials,5));

    % make fixation duration an exact multiple of refresh interval
    data.fix_t     = s.fix_t - mod(s.fix_t, flip_t);
    data.fixFrames = round(data.fix_t / flip_t);

    % make cue duration an exact multiple of refresh interval
    data.cue_t     = s.cue_t - mod(s.cue_t, flip_t);
    data.cueFrames = round(data.cue_t / flip_t);

    % placeholders for instruction button / onset (if you later re-enable)
    [data.btn_instr, data.onset_instr] = deal(cell(1,5));

    % ---------------------------------------------------------
    % 6) INSTRUCTIONS instruction screens – currently commented out
    % ---------------------------------------------------------
    if seq_block(1).block == 1
        instr_dir   = [fileparts(mfilename('fullpath')) s.instr_dir];
        instr_subdir = dir([instr_dir s.instr_subdir_wildcard '3*']);
        if ~isempty(instr_subdir)
            instr_images = load_images([instr_dir instr_subdir.name '\'], s.instr_img_wildcard);
            [data.btn_instr, data.onset_instr] = show_instr_img_EFST(instr_images, window);
            clear instr_images
        end
    end

    %% =========================================================
    % 7) MAIN TRIAL LOOP
    %% =========================================================
    for i = 1:nTrials
        tr = seq_block(i);   % trial struct for this trial

        % ---------------- FIXATION ----------------
        Screen('TextSize', window, s.cue_size);
        DrawFormattedText(window, s.fix_symbol, 'center','center', s.txt_color);
        fix_on = Screen('Flip', window);           % flip to show fixation
        data.onset_fix{i,1} = fix_on;

        % EEG: trial start trigger
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 30);    % 30 = trial start
            triggerout = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % Pupil Labs: send trial-start annotation
        if use_pupil
            zmq.core.send(socket, uint8('t'));
            currentTime_bytes = zmq.core.recv(socket,20);
            currentTimeNum  = str2double(char(currentTime_bytes));

            keys_start = {'topic','label','timestamp','duration'};
            vals_start = {'annotation.TrialStart', ...
                          sprintf('Trial #%d started',i), ...
                          currentTimeNum, 1.0};
            start_annotation = containers.Map(keys_start, vals_start);
            send_annotation(socket, start_annotation);
            result = zmq.core.recv(socket); %#ok<NASGU>
        end

        % ---------------- CUE (1 s) ----------------
        Screen('TextSize', window, s.cue_size);
        DrawFormattedText(window, s.cue, 'center','center', s.cue_color);
        % schedule cue onset at fixation onset + (fixFrames - 0.5)*flip_t
        cue_on = Screen('Flip', window, fix_on + (data.fixFrames - 0.5)*flip_t);
        data.onset_cue{i,1} = cue_on;

        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 11);   % 11 = cue onset
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % ---------- STIMULUS (face + word, 1 s) ----------
        % read image from file
        img = imread(fullfile(s.stim_dir, tr.img));
        % make PTB texture
        tex = Screen('MakeTexture', window, img);

        % get original texture rectangle [0 0 width height]
        texRect = Screen('Rect', tex);

        % scale factor for stimulus size (here: quarter of full size)
        scaleFactor = 0.25;
        scaledRect  = CenterRectOnPointd(texRect * scaleFactor, xCenter, yCenter);

        % draw face texture
        Screen('DrawTexture', window, tex, [], scaledRect);

        % set text size for Stroop word relative to image height
        imgH     = texRect(3);                % width in px (assuming square-ish)
        wordSize = round(imgH * 0.045);
        Screen('TextSize', window, wordSize);

        % draw Stroop word over the face, centered
        DrawFormattedText(window, tr.word, 'center', 'center', s.txt_color);

        % draw photodiode patch
        Screen('FillRect', window, s.white, cueRect);

        % Flip at the correct time relative to cue
        stim_on = Screen('Flip', window, cue_on + (data.cueFrames - 0.5)*flip_t);
        data.onset_stim{i,1} = stim_on;

        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 12);   % 12 = stimulus onset
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % keep stimulus for s.stim_t seconds
        WaitSecs(s.stim_t);

        % free texture memory
        Screen('Close', tex);

        % clear any keypresses during the stimulus
        KbQueueFlush;

        % =====================================================
        % RESPONSE PERIOD (within ISI)
        % =====================================================
        isi_total = tr.isi;   % trial-specific ISI duration

        % ----- draw 3 response labels in their respective boxes -----
        Screen('TextSize', window, s.txt_size);

        DrawFormattedText(window, s.resp_labels{1}, ...
                          'center', 'center', s.txt_color, [], [], [], [], [], rect1);
        DrawFormattedText(window, s.resp_labels{2}, ...
                          'center', 'center', s.txt_color, [], [], [], [], [], rect2);
        DrawFormattedText(window, s.resp_labels{3}, ...
                          'center', 'center', s.txt_color, [], [], [], [], [], rect3);

        resp_on = Screen('Flip', window);
        data.onset_resp1{i,1} = resp_on;

        % flush any residual key presses before we start response timing
        KbQueueFlush;

        % EEG: response-screen onset
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 13);   % 13 = response screen onset
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % maximum time for response screen (capped by ISI)
        resp_screen_window = min(s.resp1_max_t, isi_total);

        resp_btn = 0;       % which button was pressed (1/2/3 or ESC)
        rt       = NaN;     % reaction time
        start_t  = GetSecs; % response period start time

        % ---------------- RESPONSE LOOP ----------------
        while GetSecs - start_t < resp_screen_window
            [pressed, firstPress] = KbQueueCheck;
            if pressed
                k     = find(firstPress,1);
                kname = KbName(k);
                rt    = GetSecs - start_t;

                if strcmp(kname, s.responseKeys{1})
                    resp_btn = 1;   % Happy
                elseif strcmp(kname, s.responseKeys{2})
                    resp_btn = 2;   % Neutral
                elseif strcmp(kname, s.responseKeys{3})
                    resp_btn = 3;   % Sad
                elseif strcmpi(kname, 'ESCAPE')
                    resp_btn = s.btn_esc;  % emergency quit
                end
                break;
            end
        end

        % store raw button and RT
        data.resp1_btn(i,1) = resp_btn;
        data.resp1_t(i,1)   = rt;

        % Evaluate response category and send EEG trigger
        if resp_btn == 1
            data.resp1(i,1) = 1;   % Happy
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 41);  % e.g. 41 = "happy" response
                triggerout  = 1; %#ok<NASGU>
                triggertime = GetSecs; %#ok<NASGU>
            end

        elseif resp_btn == 2
            data.resp1(i,1) = 2;   % Neutral
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 40);  % e.g. 40 = "neutral" response
                triggerout  = 1; %#ok<NASGU>
                triggertime = GetSecs; %#ok<NASGU>
            end

        elseif resp_btn == 3
            data.resp1(i,1) = 3;   % Sad
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 42);  % e.g. 42 = "sad" response
                triggerout  = 1; %#ok<NASGU>
                triggertime = GetSecs; %#ok<NASGU>
            end

        elseif resp_btn == s.btn_esc
            % ESC: stop experiment early
            break;

        else
            % no button pressed within response window
            data.resp1(i,1) = -1;
        end

        % ---------------- FILL REMAINING ISI ----------------
        elapsed_resp  = GetSecs - resp_on;
        remaining_ISI = isi_total - elapsed_resp;

        if remaining_ISI > 0
            % show blank (background) for the remaining ISI
            Screen('FillRect', window, s.window_color);
            Screen('Flip', window);
            WaitSecs(remaining_ISI);
        end
    end  % end trial loop

    %% --------------------------------------------------------
    % 8) END OF BLOCK / CLEANUP / SAVE
    %% --------------------------------------------------------
    sca;                 % close PTB window
    ListenChar(0);       % re-enable keyboard input
    KbQueueRelease;      % release keyboard queue

    % save block data if participant directory is defined
    if isfield(p_data,'dir')
        save(fullfile(p_data.dir, ...
            [s.file_prefix sprintf('block_%02d_%s.mat', seq_block(1).block, p_data.ID)]), ...
            'data','s','seq_block');
    end

catch ERR
    % In case of any error, make sure PTB and keyboard are cleaned up
    sca;
    ListenChar(0);
    KbQueueRelease;
    rethrow(ERR);
end

end
