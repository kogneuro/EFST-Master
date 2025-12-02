function [exp_seq, s] = seq_colorStroop_3c(s)
% seq_colorStroop_3c
% ----------------------------------------------------------
% Build the full trial sequence for the 3-color Stroop task.
%
% DESIGN:
%   • 3 color words:      RED, GREEN, BLACK
%   • 3 ink colors:       RED, GREEN, BLACK
%   • 2 congruency states:
%         'C' – congruent  (word color == ink color)
%         'I' – incongruent (word color ~= ink color)
%
% BALANCING CONSTRAINTS (PER BLOCK):
%   • 50% congruent (C), 50% incongruent (I)
%   • color WORDS are roughly equally frequent
%   • INK colors are roughly equally frequent
%   • transitions (CC, CI, IC, II) are roughly balanced
%   • ISIs are in [s.isi_min, s.isi_max] with mean ≈ s.isi_mean
%
% INPUT:
%   s       – settings struct (comes from setup_colorStroop_3c)
%
% OUTPUT:
%   exp_seq – struct array with one element per trial:
%             .block, .trial
%             .word, .wordIdx
%             .inkColor (RGB), .inkColorIdx, .inkName
%             .prevCong, .currCong, .congruency
%             .isi
%   s       – same settings struct, with rng_state recorded
%
% Alp, 2025
% ----------------------------------------------------------

%% 1) RNG STATE
% Store RNG state for reproducibility and debugging
s.rng_state = set_rng_state(s);

% Initialize final sequence container
exp_seq = [];

% Convenience variables
nTrialsPerBlock = s.trials_per_block;
nBlocks         = s.blocks;

%% 2) BASIC CHECKS
% We want an exact 50/50 split of C / I per block,
% so we require an even number of trials per block.
if mod(nTrialsPerBlock,2) ~= 0
    error('ColorStroop:trials_per_block', ...
          'trials_per_block must be even (for 50/50 C/I).');
end

%% 3) COLOR WORD DEFINITIONS & TARGETS
% Color word labels (e.g., {'RED','GREEN','BLACK'})
colorWords = s.color_words;
nColors    = numel(colorWords);

% Target count per block for WORDS and for INK
% Start with equal distribution
baseColor = floor(nTrialsPerBlock / nColors);
remColor  = mod(nTrialsPerBlock, nColors);

% wordTargets(b, c) = how many times word c should appear in block b
% inkTargets(b, c)  = how many times ink color c should appear in block b
wordTargets = zeros(nBlocks, nColors);
inkTargets  = zeros(nBlocks, nColors);

for b = 1:nBlocks
    % base count for each color
    wordTargets(b,:) = baseColor;
    inkTargets(b,:)  = baseColor;

    % distribute remainders (if nTrialsPerBlock is not a multiple of 3)
    if remColor > 0
        idxExtra = mod(b-1, nColors) + 1;  % rotate the extra across blocks
        wordTargets(b, idxExtra) = wordTargets(b, idxExtra) + remColor;
        inkTargets(b, idxExtra)  = inkTargets(b, idxExtra)  + remColor;
    end
end

%% 4) BUILD SEQUENCE BLOCK-BY-BLOCK
for b = 1:nBlocks

    % 4.1) Build congruency path C/I with balanced transitions
    %      This returns a cell array of length nTrialsPerBlock
    %      with entries 'C' or 'I'
    cong_seq = build_cong_sequence_Color(nTrialsPerBlock);

    % 4.2) Initialize counters for word and ink usage (per color)
    cntWord = zeros(1,nColors);          % how many times each word used so far
    cntInk  = zeros(1,nColors);          % how many times each ink used so far
    tgtWord = wordTargets(b,:);          % target counts this block
    tgtInk  = inkTargets(b,:);           % target counts this block

    % Anti-run: avoid more than maxRunInk same ink color in a row
    maxRunInk = 3;

    % Container for this block’s trials
    block_trials = struct([]);

    % 4.3) TRIAL LOOP WITHIN BLOCK
    for t = 1:nTrialsPerBlock

        % Current congruency type for this trial
        currCong = cong_seq{t};

        % Previous congruency type (for logging)
        if t == 1
            prevCong = currCong;
        else
            prevCong = cong_seq{t-1};
        end

        % ------------------------------------------------------
        % 4.3.1) CHOOSE WORD INDEX (identity) NEAR TARGET COUNTS
        % ------------------------------------------------------
        % underWord(c) = true if we still need more of word c
        underWord = cntWord < tgtWord;
        idxUW = find(underWord);

        if isempty(idxUW)
            % all words at or above target → just pick randomly
            wordIdx = randi(nColors);
        elseif numel(idxUW) == 1
            % only one color word still under target
            wordIdx = idxUW;
        else
            % multiple words still under target → pick among them
            wordIdx = idxUW(randi(numel(idxUW)));
        end

        % ------------------------------------------------------
        % 4.3.2) CHOOSE INK INDEX GIVEN CONGRUENCY
        % ------------------------------------------------------
        if currCong == 'C'
            % Congruent: ink color must match the word
            inkCandidates = wordIdx;
        else
            % Incongruent: any ink color DIFFERENT from the word
            inkCandidates = setdiff(1:nColors, wordIdx);
        end

        % Among the candidate inks, prefer those under their ink target
        good = inkCandidates(cntInk(inkCandidates) < tgtInk(inkCandidates));
        if isempty(good)
            % all candidates at/above target → choose randomly among them
            inkIdx = inkCandidates(randi(numel(inkCandidates)));
        elseif numel(good) == 1
            inkIdx = good;
        else
            inkIdx = good(randi(numel(good)));
        end

        % ------------------------------------------------------
        % 4.3.3) ANTI-RUN: AVOID LONG RUNS OF SAME INK COLOR
        % ------------------------------------------------------
        if last_n_are(block_trials, 'inkColorIdx', inkIdx, maxRunInk)
            % Try any other ink color to break the run
            alt = setdiff(1:nColors, inkIdx);

            % Prefer alternatives that are still under their ink target
            altValid = alt(cntInk(alt) < tgtInk(alt));
            if isempty(altValid)
                % If none under target, just pick among all others
                altValid = alt;
            end

            % Replace inkIdx by a random valid alternative
            inkIdx = altValid(randi(numel(altValid)));
        end

        % ------------------------------------------------------
        % 4.3.4) FINALIZE WORD AND INK COLOR (RGB)
        % ------------------------------------------------------
        wordStr = colorWords{wordIdx};
        switch inkIdx
            case 1
                rgb     = s.colors.RED;
                inkName = 'RED';
            case 2
                rgb     = s.colors.GREEN;
                inkName = 'GREEN';
            case 3
                rgb     = s.colors.BLACK;
                inkName = 'BLACK';
            otherwise
                error('Invalid inkIdx %d', inkIdx);
        end

        % Label congruency as text
        if currCong == 'C'
            congLabel = 'congruent';
        else
            congLabel = 'incongruent';
        end

        % ------------------------------------------------------
        % 4.3.5) UPDATE WORD AND INK COUNTS
        % ------------------------------------------------------
        cntWord(wordIdx) = cntWord(wordIdx) + 1;
        cntInk(inkIdx)   = cntInk(inkIdx)   + 1;

        % ------------------------------------------------------
        % 4.3.6) STORE TRIAL INFO IN STRUCT
        % ------------------------------------------------------
        block_trials(t).block       = b;
        block_trials(t).trial       = t;

        block_trials(t).word        = wordStr;  % e.g., 'RED'
        block_trials(t).wordIdx     = wordIdx;  % 1/2/3 for RED/GREEN/BLACK
        block_trials(t).inkColor    = rgb;      % [R G B]
        block_trials(t).inkColorIdx = inkIdx;   % 1/2/3
        block_trials(t).inkName     = inkName;  % 'RED', 'GREEN', 'BLACK'

        block_trials(t).prevCong    = prevCong; % previous C/I
        block_trials(t).currCong    = currCong; % current  C/I
        block_trials(t).congruency  = congLabel; % 'congruent'/'incongruent'
    end

    % ----------------------------------------------------------
    % 4.4) ASSIGN ISI VALUES FOR THIS BLOCK
    % ----------------------------------------------------------
    % ISIs are random in [s.isi_min, s.isi_max] but shifted to match s.isi_mean
    isi_vec = makeISI_block(nTrialsPerBlock, s.isi_min, s.isi_max, s.isi_mean);
    for t = 1:nTrialsPerBlock
        block_trials(t).isi = isi_vec(t);
    end

    % Append this block’s trials to the global sequence
    exp_seq = [exp_seq; block_trials(:)]; %#ok<AGROW>
end

end  % ---- main function end ----


%% ========================================================================
function cong_seq = build_cong_sequence_Color(nTrials)
% build_cong_sequence_Color
% ----------------------------------------------------------
% Build a sequence of 'C' / 'I' (congruent / incongruent) with:
%   • 50% C, 50% I
%   • transitions CC, CI, IC, II approximately balanced
%
% INPUT:
%   nTrials – length of sequence (must be even)
%
% OUTPUT:
%   cong_seq – 1 x nTrials cell array of 'C' / 'I'
% ----------------------------------------------------------

% We want exactly half C and half I
if mod(nTrials,2) ~= 0
    error('build_cong_sequence_Color: nTrials must be even.');
end

% Total counts of states
nC = nTrials / 2;
nI = nTrials / 2;

% Number of transitions (between trials)
nTrans = nTrials - 1;

% Distribute transitions CC, CI, IC, II as evenly as possible
base = floor(nTrans / 4);
rem  = mod(nTrans, 4);

targ.CC = base;
targ.CI = base;
targ.IC = base;
targ.II = base;

% Assign the leftover transitions in order CC, CI, IC, II
order = {'CC','CI','IC','II'};
for k = 1:rem
    targ.(order{k}) = targ.(order{k}) + 1;
end

maxAttempts = 1000;

for attempt = 1:maxAttempts

    % Remaining counts for states
    remC  = nC;
    remI  = nI;

    % Remaining counts for transitions
    remCC = targ.CC;
    remCI = targ.CI;
    remIC = targ.IC;
    remII = targ.II;

    % Preallocate output
    cong_seq = cell(1, nTrials);

    % Randomly choose first state: C or I
    if rand < 0.5
        cong_seq{1} = 'C';
        remC = remC - 1;
    else
        cong_seq{1} = 'I';
        remI = remI - 1;
    end

    ok = true;

    % Fill rest of sequence
    for t = 2:nTrials

        prev = cong_seq{t-1};
        candidates = {};

        % Build possible next states, given remaining transition & state counts
        if strcmp(prev,'C')
            % Options: C->C (CC) or C->I (CI)
            if remCC > 0 && remC > 0
                candidates{end+1} = 'C';
            end
            if remCI > 0 && remI > 0
                candidates{end+1} = 'I';
            end
        else
            % prev == 'I' → options: I->C (IC) or I->I (II)
            if remIC > 0 && remC > 0
                candidates{end+1} = 'C';
            end
            if remII > 0 && remI > 0
                candidates{end+1} = 'I';
            end
        end

        % If we have no valid next state, abort this attempt
        if isempty(candidates)
            ok = false;
            break
        end

        % Choose next state from candidates
        if numel(candidates) == 1
            nxt = candidates{1};
        else
            nxt = candidates{randi(numel(candidates))};
        end

        % Place the chosen state
        cong_seq{t} = nxt;

        % Decrease remaining state count
        if strcmp(nxt,'C')
            remC = remC - 1;
        else
            remI = remI - 1;
        end

        % Decrease remaining transition count
        if strcmp(prev,'C') && strcmp(nxt,'C')
            remCC = remCC - 1;
        elseif strcmp(prev,'C') && strcmp(nxt,'I')
            remCI = remCI - 1;
        elseif strcmp(prev,'I') && strcmp(nxt,'C')
            remIC = remIC - 1;
        else % prev='I' & nxt='I'
            remII = remII - 1;
        end
    end

    % Check if we used all states and transitions exactly
    if ok && remC == 0 && remI == 0 && ...
             remCC == 0 && remCI == 0 && remIC == 0 && remII == 0
        % Success – return this sequence
        return
    end
end

% If we get here, all attempts failed
error('build_cong_sequence_Color:CouldNotBuild', ...
      'Could not build a valid 2-state C/I sequence in %d attempts.', maxAttempts);
end


%% ========================================================================
% Helper: inline if
function out = iff(cond,a,b)
% iff – simple inline if-else
if cond, out = a; else, out = b; end
end

%% Helper: ISI block
function isi = makeISI_block(n, lo, hi, target)
% makeISI_block
%   Draw n random ISIs in [lo, hi] and shift them so that
%   the resulting mean is approximately "target".

% Uniform random in [lo, hi]
isi = lo + (hi - lo) * rand(n,1);

% Current mean
curMean = mean(isi);

% Shift all ISIs by (target - curMean)
isi = isi + (target - curMean);

% Clamp to [lo, hi]
isi(isi < lo) = lo;
isi(isi > hi) = hi;
end

%% Helper: anti-run checker
function tf = last_n_are(trials, fieldname, value, n)
% last_n_are
%   Return true if the last n entries in "trials" have
%   the given field "fieldname" equal to "value".

if numel(trials) < n
    tf = false;
    return;
end

tf = true;
for k = numel(trials)-n+1 : numel(trials)
    % field must exist and match the desired value
    if ~isfield(trials(k), fieldname) || ...
       ~isequal(trials(k).(fieldname), value)
        tf = false;
        return;
    end
end
end
