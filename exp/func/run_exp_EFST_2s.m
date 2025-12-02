function data = run_exp_EFST_2s(s, p_data, seq_block)
% ============================================================
% run_exp_EFST_2s
% Emotional Face Stroop task – 2 emotions (Happy / Sad)
%
% Trial structure:
%   1) Fixation cross:      s.fix_t   (e.g. 0.04 s)
%   2) Cue:                 s.cue_t   (e.g. 1 s)
%   3) Stimulus:            s.stim_t  (e.g. 1 s; face + word)
%   4) Response screen:     up to 2 s (but limited by ISI)
%   5) ISI:                 2–4 s total, including response window
%
% Inputs:
%   s         – setup struct from setup_EFST_2s (timing, screen, keys, etc.)
%   p_data    – participant_data struct (ID, folder, etc.)
%   seq_block – 1 block of trials (struct array from seq_EFST_2s)
%
% Output:
%   data      – struct with:
%                  • rng_state
%                  • seq (copy of seq_block)
%                  • onsets: onset_fix, onset_cue, onset_stim, onset_resp1
%                  • resp1, resp1_t, resp1_btn
%                  • flip_t, stimFrames, window size, etc.
%
% All potential crashes are caught by try/catch so that the screen,
% keyboard, and queues are cleaned up properly.
%
% Alp — 2025
% ============================================================

try
    %% --------------------------------------------------------
    % 1) EYETRACKING INITIALIZATION (Pupil Labs, optional)
    %% --------------------------------------------------------
    % Only run this block if eyetracking is enabled and you are using
    % Pupil Labs (whicheyetracker == 2). The code sets up a ZMQ
    % connection, synchronizes time, and starts recording.
    use_pupil = false;
    if isfield(s,'eyetracking') && s.eyetracking == 1
        if isfield(s,'whicheyetracker') && s.whicheyetracker == 2
            % Pupil Remote address – configured in Pupil Capture
            endpoint = 'tcp://000.000.00.00:00000';

            % Create ZMQ context and REQ socket
            ctx    = zmq.core.ctx_new();
            socket = zmq.core.socket(ctx, 'ZMQ_REQ');

            % Set 1000 ms timeout so Matlab does not block forever if
            % the server is not reachable
            zmq.core.setsockopt(socket, 'ZMQ_RCVTIMEO', 1000);

            fprintf('Connecting to %s\n', endpoint);
            zmq.core.connect(socket, endpoint);

            % Request SUB port from Pupil Labs
            zmq.core.send(socket, uint8('SUB_PORT'));
            sub_port = char(zmq.core.recv(socket));
            fprintf('Received sub port: %s\n', sub_port);

            if isequal(sub_port, false)
                warning('No valid sub port received from PupilLabs');
                return;  % exit script if no valid port
            end

            % Create and connect SUB socket for streaming data
            ip_address  = '000.000.00.00';
            sub_endpoint = sprintf('tcp://%s:%s', ip_address, sub_port);
            sub_socket   = zmq.core.socket(ctx, 'ZMQ_SUB');

            zmq.core.connect(sub_socket, sub_endpoint);
            zmq.core.setsockopt(sub_socket, 'ZMQ_RCVTIMEO', 1000);
            fprintf('Connecting to SUB: %s\n', sub_endpoint);

            % Measure round-trip delay for a 't' command
            tic;
            zmq.core.send(socket, uint8('t'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));
            fprintf('Round trip command delay: %s\n', toc);

            % Set current Pupil time to 0.0
            zmq.core.send(socket, uint8('T 0.0'));
            result = zmq.core.recv(socket);
            fprintf('%s\n', char(result));

            % Start recording
            pause(1.0);
            zmq.core.send(socket, uint8('R'));
            use_pupil = true;
            fprintf('Pupil recording should start.\n');
        end
    end

    %% --------------------------------------------------------
    % 2) RANDOM SEED + SEQUENCE SETUP
    %% --------------------------------------------------------
    % Store the RNG state so you can reconstruct exactly the same
    % random behavior if you need to debug timing or sequences later.
    data.rng_state = set_rng_state(s);
    data.seq       = seq_block;
    nTrials        = numel(seq_block);

    %% --------------------------------------------------------
    % 3) KEYBOARD INITIALIZATION
    %% --------------------------------------------------------
    % Use KbQueue for responsive, non-blocking key collection.
    KbName('UnifyKeyNames');    % consistent key codes across OS
    ListenChar(-1);             % suppress keyboard input to MATLAB window

    % Create a key mask (only our keys of interest are listened to)
    keylist = zeros(1,256);
    keylist(KbName(s.responseKeys{1})) = 1;  % first response key (Happy)
    keylist(KbName(s.responseKeys{2})) = 1;  % second response key (Sad)
    keylist(KbName('ESCAPE'))          = 1;  % emergency exit

    KbQueueCreate(-1, keylist);  % -1 = default keyboard
    KbQueueStart;                % start queue

    %% --------------------------------------------------------
    % 4) SCREEN SETUP (Psychtoolbox)
    %% --------------------------------------------------------
    AssertOpenGL;                 % PTB requires OpenGL
    PsychDefaultSetup(2);         % default settings, normalized color range

    % Open main experiment window on s.screenNumber, fill with background
    [window, windowRect] = Screen('OpenWindow', s.screenNumber, s.window_color);

    % Ensure 0–255 color range
    Screen('ColorRange', window, 255);

    % Enable alpha blending (needed for transparency, anti-aliasing, etc.)
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    HideCursor;                   % no mouse cursor during the experiment

    % Center of the window (for positioning stimuli)
    [xCenter, yCenter] = RectCenter(windowRect);

    % Photodiode rectangle at the top-center (optional)
    cueRect = [0 0 s.diodeSize s.diodeSize];
    cueRect = CenterRectOnPointd(cueRect, xCenter, s.diodeSize/2 + 10);

    % Set font and size for main Stroop text (will be adjusted for word)
    Screen('TextFont', window, s.txt_font);
    Screen('TextSize', window, s.txt_size);

    % Measure refresh interval (frame duration) using PTB’s GetFlipInterval
    Priority(1);
    flip_t = Screen('GetFlipInterval', window, s.get_flip_i);
    Priority(0);
    data.flip_t = flip_t;

    % Number of frames to show the stimulus for s.stim_t seconds
    data.stimFrames = round(s.stim_t / data.flip_t);

    % Store window size [width, height]
    [data.window(1), data.window(2)] = Screen('WindowSize', window);

    % For positioning response labels, we also grab screen size separately
    [screenX, screenY] = Screen('WindowSize', window);

    % Row (y-position) for the response boxes (middle of the screen)
    yPos = screenY * 0.50;

    % Width and height of each response rectangle
    w = 300;
    h = 80;

    % Two rectangles – left and right – centered at 1/4 and 3/4 width
    rect1 = CenterRectOnPointd([0 0 w h], screenX*(1/4), yPos); % left box
    rect2 = CenterRectOnPointd([0 0 w h], screenX*(3/4), yPos); % right box

    % Initial paint/flip to clear screen
    Screen('Flip', window);

    %% --------------------------------------------------------
    % 5) TIMING PREALLOCATION
    %% --------------------------------------------------------
    % We store onset times (in seconds since GetSecs) for:
    %   fixation, cue, stim, (optional delay), response screen, etc.
    [data.onset_fix, ...
     data.onset_cue, ...
     data.onset_stim, ...
     data.onset_delay, ...
     data.onset_resp1, ...
     data.onset_resp1_p] = deal(cell(nTrials,5));

    % Make fixation duration divisible by flip_t so we can align flips
    data.fix_t     = s.fix_t - mod(s.fix_t, flip_t);
    data.fixFrames = round(data.fix_t / flip_t);

    % Same for the cue
    data.cue_t     = s.cue_t - mod(s.cue_t, data.flip_t);
    data.cueFrames = round(data.cue_t / flip_t);

    % Reserved fields for instruction button presses and onsets
    [data.btn_instr, ...
     data.onset_instr] = deal(cell(1,5));

    % --------------------------------------------------------
    % 6) INSTRUCTIONS
    % --------------------------------------------------------
    if seq_block(1).block == 1
        instr_dir   = [fileparts(mfilename('fullpath')) s.instr_dir];
        instr_subdir = dir([instr_dir s.instr_subdir_wildcard '2*']);
        if ~isempty(instr_subdir)
            instr_images = load_images([instr_dir instr_subdir.name '\'], ...
                                       s.instr_img_wildcard);
            [data.btn_instr, data.onset_instr] = ...
                show_instr_img(instr_images, window);
            clear instr_images
        end
    end

    %% =========================================================
    % 7) MAIN TRIAL LOOP
    %% =========================================================
    for i = 1:nTrials

        % Current trial structure (from seq_EFST_2s)
        tr = seq_block(i);

        % ------------------ FIXATION (s.fix_t) ------------------
        Screen('TextSize', window, s.cue_size);
        DrawFormattedText(window, s.fix_symbol, 'center', 'center', s.txt_color);

        % Flip immediately to show fixation
        fix_on = Screen('Flip', window);
        data.onset_fix{i,1} = fix_on;

        % Send EEG trigger for "trial start" if EEG is enabled
        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 30);  % 30 = trial start
            % (triggerout/triggertime kept for compatibility,
            % though not used below)
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % (Optional) Send annotation to Pupil Labs: TrialStart
        if use_pupil
            zmq.core.send(socket, uint8('t'));
            currentTime_bytes = zmq.core.recv(socket,20);
            currentTimeChar   = char(currentTime_bytes);
            currentTimeNum    = str2double(currentTimeChar);
            keys_start  = {'topic','label','timestamp','duration'};
            values_start = {'annotation.TrialStart', ...
                            sprintf('Trial #%d started', i), ...
                            currentTimeNum, 1.0};
            start_annotation = containers.Map(keys_start, values_start);
            send_annotation(socket, start_annotation);
            result = zmq.core.recv(socket); %#ok<NASGU>
        end

        % ------------------ CUE (s.cue_t) ------------------------
        Screen('TextSize', window, s.cue_size);
        DrawFormattedText(window, s.cue, 'center', 'center', s.cue_color);

        % Flip after the fixation period has elapsed
        cue_on = Screen('Flip', window, fix_on + (data.fixFrames - 0.5)*flip_t);
        data.onset_cue{i,1} = cue_on;

        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 11); % 11 = warning cue onset
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % --------------- STIMULUS (face + word) -----------------
        % Load face image from disk and make PTB texture
        img = imread(fullfile(s.stim_dir, tr.img));
        tex = Screen('MakeTexture', window, img);

        % Original size of the texture [0 0 width height]
        texRect = Screen('Rect', tex);

        % Scale factor – currently 0.25 (25% of native size)
        scaleFactor = 0.25;

        % Scaled rectangle centered in the screen
        scaledRect = CenterRectOnPointd(texRect * scaleFactor, xCenter, yCenter);

        % Draw the face
        Screen('DrawTexture', window, tex, [], scaledRect);

        % Adjust word size relative to texture size
        % NOTE: texRect(3) is the width, so imgW would be better name.
        imgW     = texRect(3);
        wordSize = round(imgW * 0.045);
        Screen('TextSize', window, wordSize);

        % Draw the Stroop word (GLÜCKLICH / TRAURIG) in front
        DrawFormattedText(window, tr.word, 'center', 'center', s.txt_color);

        % Photodiode patch set to white
        Screen('FillRect', window, s.white, cueRect);

        % Flip to show the face + word exactly after the cue period
        stim_on = Screen('Flip', window, cue_on + (data.cueFrames - 0.5)*flip_t);
        data.onset_stim{i,1} = stim_on;

        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 12); % 12 = stimulus onset
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % Keep the stimulus on for s.stim_t seconds
        WaitSecs(s.stim_t);

        % Free texture memory (important in longer experiments)
        Screen('Close', tex);

        % Clear any keypresses that happened during the stimulus
        KbQueueFlush;

        % =========================================================
        %  RESPONSE SCREEN (within ISI)
        % =========================================================
        % Total ISI for this trial (in seconds) was pre-assigned in seq
        isi_total = tr.isi;

        % Show response options immediately after stimulus
        Screen('TextSize', window, s.txt_size);

        % Left box: GLÜCKLICH, Right box: TRAURIG (from s.resp_labels)
        DrawFormattedText(window, s.resp_labels{1}, ...
            'center','center', s.txt_color, [], [], [], [], [], rect1);
        DrawFormattedText(window, s.resp_labels{2}, ...
            'center','center', s.txt_color, [], [], [], [], [], rect2);

        % Flip to response screen
        resp_on = Screen('Flip', window);
        data.onset_resp1{i,1} = resp_on;

        % Flush queue again to avoid pre-existing keypresses
        KbQueueFlush;

        if isfield(s,'EEG') && s.EEG == 1
            SendSignal(s.ioObj, s.address, 13); % 13 = response screen onset
            triggerout  = 1; %#ok<NASGU>
            triggertime = GetSecs; %#ok<NASGU>
        end

        % Response window = min(3 s, total ISI)
        resp_screen_window = min(s.resp1_max_t, isi_total);

        resp_btn = 0;      % which key was pressed (1=Happy, 2=Sad, 3=ESC)
        rt       = NaN;    % reaction time (s)
        start_t  = GetSecs;

        % Poll KbQueue until time runs out or a valid key is pressed
        while GetSecs - start_t < resp_screen_window
            [pressed, firstPress] = KbQueueCheck;
            if pressed
                k     = find(firstPress, 1);
                kname = KbName(k);
                rt    = GetSecs - start_t;

                % Map pressed key to response code
                if strcmp(kname, s.responseKeys{1})
                    resp_btn = 1;   % Happy
                elseif strcmp(kname, s.responseKeys{2})
                    resp_btn = 2;   % Sad
                elseif strcmpi(kname, 'ESCAPE')
                    resp_btn = s.btn_esc;
                end
                break;
            end
        end

        % Log response button and reaction time
        data.resp1_btn(i,1) = resp_btn;
        data.resp1_t(i,1)   = rt;

        % Interpret response + send EEG triggers
        if resp_btn == 1
            data.resp1(i,1) = 1;  % Happy
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 41);  % example code for "Happy"
            end

        elseif resp_btn == 2
            data.resp1(i,1) = 2;  % Sad
            if isfield(s,'EEG') && s.EEG == 1
                SendSignal(s.ioObj, s.address, 40);  % example code for "Sad"
            end

        elseif resp_btn == s.btn_esc
            % Emergency escape: break out of the trial loop
            break;

        else
            % No response
            data.resp1(i,1) = -1;
        end

        % Time spent from response screen onset until we exit the loop
        elapsed_resp  = GetSecs - resp_on;

        % Remaining ISI (could be 0 or negative if response took longer)
        remaining_ISI = isi_total - elapsed_resp;

        % If there is remaining ISI > 0, we show a blank (background) screen
        % to fill up the planned ISI before next trial.
        if remaining_ISI > 0
            Screen('FillRect', window, s.window_color);
            Screen('Flip', window);
            WaitSecs(remaining_ISI);
        end
    end  % end of trial loop

    %% --------------------------------------------------------
    % 8) END OF BLOCK / CLEANUP / SAVE
    %% --------------------------------------------------------
    sca;               % close Psychtoolbox window
    ListenChar(0);     % re-enable keyboard in MATLAB
    KbQueueRelease;    % release KbQueue

    % Save data for this block if participant directory exists
    if isfield(p_data,'dir')
        save(fullfile(p_data.dir, ...
            [s.file_prefix sprintf('block_%02d_%s.mat', ...
             seq_block(1).block, p_data.ID)]), ...
             'data','s','seq_block');
    end

catch ERR
    % --------- EMERGENCY CLEANUP ON CRASH --------------------
    sca;
    ListenChar(0);
    KbQueueRelease;
    rethrow(ERR);  % bubble the error up so you still see the message
end

end
