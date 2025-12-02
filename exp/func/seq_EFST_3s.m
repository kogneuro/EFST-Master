function [exp_seq, s] = seq_EFST_3s(s)
% seq_EFST_3s – build full sequence for Emotional Face Stroop (3 feelings)
%
% Faces:       h / n / s      → happy / neutral / sad
% Congruency:  C / I / IS     → congruent / incongruent / slightly incongruent
%
% Design goals:
%   • No image reuse across the whole experiment (each jpg used at most once)
%   • Balanced age (y/m/o) per block
%   • Balanced gender (F / M) per block
%   • Approx-balanced feeling × congruency cells (3 feelings × 3 cong = 9 cells)
%   • 3-state congruency path with balanced counts and transitions (via
%     build_cong_sequence_EFST3)
%   • ISI per trial in [isi_min, isi_max] with block-wise mean ≈ isi_mean
%
% Input:
%   s       – setup struct from setup_EFST_3s (fields: stim_dir, blocks, trials_per_block, etc.)
%
% Output:
%   exp_seq – column vector of trial structs (one per trial, across all blocks)
%   s       – setup struct with stored RNG state (s.rng_state)
%
% Alp / 2025
% -------------------------------------------------------------------------


%% 1) RNG
% Store the current RNG state in s (for reproducibility) and create output container
s.rng_state = set_rng_state(s);
exp_seq = [];

%% 2) Load all jpg stimuli from stimulus folder
files_all = dir(fullfile(s.stim_dir, '*.jpg'));
keep_idx  = endsWith({files_all.name}, '.jpg');  % here it's redundant but explicit
files     = files_all(keep_idx);
clear files_all keep_idx

if isempty(files)
    error('seq_EFST_3s:NoStim', 'No *.jpg files found in %s', s.stim_dir);
end

% -------------------------------------------------------------------------
% Parse filenames into a "pool" array of structs.
% Expecting scheme: subj_age_gender_feeling*.jpg
%   subj   – subject ID
%   age    – y/m/o
%   gender – f/m
%   feeling– h/n/s (possibly + suffix)
% -------------------------------------------------------------------------
pool = [];
for iFile = 1:numel(files)
    fname = files(iFile).name;
    parts = split(fname, '_');
    if numel(parts) < 4
        % if filename does not have at least 4 underscore parts, skip it
        continue
    end

    subj    = parts{1};
    age     = lower(parts{2});     % y / m / o
    gender  = lower(parts{3});     % f / m
    feeling = lower(parts{4});     % h / n / s + maybe extra characters

    % strip ".jpg" from the last token (in case it’s directly appended)
    feeling = erase(feeling, '.jpg');

    % keep only valid categories
    if ~ismember(age, {'y','m','o'}), continue, end
    if ~ismember(gender, {'f','m'}), continue, end
    if ~ismember(feeling, {'h','n','s'}), continue, end

    % add to pool
    pool(end+1).fname   = fname;   %#ok<AGROW>
    pool(end).subj      = subj;
    pool(end).age       = age;
    pool(end).gender    = gender;
    pool(end).feeling   = feeling;
end

if isempty(pool)
    error('seq_EFST_3s:NoValidStim', ...
        'Could not parse any *.jpg into (age,gender,feeling). Check filenames.');
end


%% 3) Trial counts and basic checks
nTrialsPerBlock = s.trials_per_block;
nBlocks         = s.blocks;
nTotalTrials    = nTrialsPerBlock * nBlocks;

% Make sure we have enough unique images for all trials
if numel(pool) < nTotalTrials
    error('seq_EFST_3s:TooFewImages', ...
        'Need %d images (no repeats) but found %d.', nTotalTrials, numel(pool));
end


%% 4) Per-block target distributions
% Define category labels
ageCats      = {'y','m','o'};   nAge  = numel(ageCats);
feelingCats  = {'h','n','s'};   nFeel = numel(feelingCats);
congCats     = {'C','I','IS'};  nCong = numel(congCats);

% -------- Age targets per block (as evenly as possible) --------
baseAge = floor(nTrialsPerBlock / nAge);   % minimal count per age per block
remAge  = mod(nTrialsPerBlock, nAge);      % leftover trials to distribute

ageTargets = zeros(nBlocks, nAge);
for b = 1:nBlocks
    ageTargets(b,:) = baseAge;
    if remAge > 0
        % give the leftover trials to age categories, rotating across blocks
        idxExtra = mod(b-1, nAge) + 1;
        ageTargets(b, idxExtra) = ageTargets(b, idxExtra) + remAge;
    end
end

% -------- Feeling × congruency targets per block --------
% 3 feelings × 3 congruency states = 9 cells per block
baseCell = floor(nTrialsPerBlock / (nFeel * nCong));  % minimal count per cell
remCell  = mod(nTrialsPerBlock,  nFeel * nCong);      % leftover trials

% We use the SAME cell targets in each block (simpler, usually fine).
cellTarget = struct();
for f = 1:nFeel
    for c = 1:nCong
        cellTarget.(feelingCats{f}).(congCats{c}) = baseCell;
    end
end

% Distribute leftovers over some cells in a deterministic but arbitrary way
% This ensures that each cell differs at most by 1 trial.
if remCell > 0
    idx = 1;
    for r = 1:remCell
        f = mod(idx-1, nFeel) + 1;             % feeling index
        c = floor((idx-1) / nFeel) + 1;        % congruency index
        cellTarget.(feelingCats{f}).(congCats{c}) = ...
            cellTarget.(feelingCats{f}).(congCats{c}) + 1;
        idx = idx + 1;
    end
end
% If nTrialsPerBlock is a multiple of 9, all cells have exactly baseCell.
% Otherwise some cells have baseCell+1, but never more than +1.


%% 5) Build each block
for b = 1:nBlocks

    % 5a) Create a 3-state congruency path (balanced C/I/IS)
    cong_seq = build_cong_sequence_EFST3(nTrialsPerBlock);

    % 5b) Gender target: enforce 50% F, 50% M within each block
    if mod(nTrialsPerBlock,2) ~= 0
        error('seq_EFST_3s:EvenTrialsNeeded', ...
            'trials_per_block must be even to balance gender.');
    end
    target_gender = nTrialsPerBlock / 2;   % target F and M counts

    % gender counters
    cntF = 0; 
    cntM = 0;

    % age counters for current block
    cntAge = struct('y',0,'m',0,'o',0);
    tgtAge = struct('y', ageTargets(b,1), ...
                    'm', ageTargets(b,2), ...
                    'o', ageTargets(b,3));

    % feeling × congruency counters for this block
    cellCount = struct();
    for f = 1:nFeel
        ff = feelingCats{f};          % e.g. 'h', 'n', or 's'
        cellCount.(ff).C  = 0;
        cellCount.(ff).I  = 0;
        cellCount.(ff).IS = 0;
    end

    block_trials = struct([]);        % container for this block’s trials

    % Anti-run limits for stimulus properties
    maxRunGender  = 3;  % no more than 3 same gender in a row
    maxRunAge     = 3;  % no more than 3 same age group in a row
    maxRunFeeling = 3;  % no more than 3 same feeling in a row

    % ----------------- TRIAL LOOP FOR THIS BLOCK -----------------
    for t = 1:nTrialsPerBlock

        currCong = cong_seq{t};   % current congruency code ('C','I','IS')
        if t == 1
            prevCong = currCong;  % define prevCong for first trial (self)
        else
            prevCong = cong_seq{t-1};
        end

        %% -------- Choose gender to move towards 50/50 --------
        if cntF < target_gender && cntM < target_gender
            % both under target → random split
            wantGender = iff(rand < 0.5, 'f', 'm');
        elseif cntF < target_gender
            wantGender = 'f';
        else
            wantGender = 'm';
        end

        %% -------- Choose feeling to balance feeling × congruency cell --------
        % Compute which feelings are still under their target in the current
        % congruency column (C/I/IS).
        underCell = false(1,nFeel);
        for f = 1:nFeel
            ff = feelingCats{f};
            if cellCount.(ff).(currCong) < cellTarget.(ff).(currCong)
                underCell(f) = true;
            end
        end

        candFeelIdx = find(underCell);
        if isempty(candFeelIdx)
            % If all feelings are at or above their targets for this
            % congruency column, allow any feeling (rare, late in block).
            candFeelIdx = 1:nFeel;
        end

        % Randomly pick one of the allowed feelings
        wantFeeling = feelingCats{ candFeelIdx(randi(numel(candFeelIdx))) };

        %% -------- Choose age to approach target distribution --------
        underAge = [cntAge.y < tgtAge.y, ...
                    cntAge.m < tgtAge.m, ...
                    cntAge.o < tgtAge.o];
        idxUA = find(underAge);
        if isempty(idxUA)
            % all ages at or above target → choose any
            wantAge = ageCats{randi(nAge)};
        elseif numel(idxUA) == 1
            % only one age still under target → use that
            wantAge = ageCats{idxUA};
        else
            % several under target → randomly choose among them
            wantAge = ageCats{ idxUA(randi(numel(idxUA))) };
        end

        %% -------- First pass: find matching stimulus in pool --------
        cand_idx = find(strcmp({pool.gender},  wantGender) & ...
                        strcmp({pool.feeling}, wantFeeling) & ...
                        strcmp({pool.age},     wantAge));
        if isempty(cand_idx)
            error('seq_EFST_3s:MissingStratum', ...
                'No remaining image for block %d trial %d with gender=%s, feeling=%s, age=%s', ...
                b, t, wantGender, wantFeeling, wantAge);
        end

        % Pick one candidate at random
        pick = cand_idx(randi(numel(cand_idx)));
        stim = pool(pick);

        %% =====================================================
        % ANTI-RUN SECTION (same logic as 2-state EFST)
        %% =====================================================

        % 1) Avoid long runs of same gender
        if last_n_are(block_trials, 'gender', stim.gender, maxRunGender)
            % try other gender (but keep feeling and age)
            otherGender = iff(strcmp(stim.gender,'f'), 'm', 'f');
            alt_idx = find_stim_3s(pool, otherGender, wantFeeling, wantAge);
            if ~isempty(alt_idx)
                pick = alt_idx(randi(numel(alt_idx)));
                stim = pool(pick);
            end
        end

        % 2) Avoid long runs of same age
        if last_n_are(block_trials, 'age', stim.age, maxRunAge)
            altAges = {'y','m','o'};
            altAges(strcmp(altAges, stim.age)) = [];   % remove current age
            for aa = 1:numel(altAges)
                alt_idx = find_stim_3s(pool, stim.gender, wantFeeling, altAges{aa});
                if ~isempty(alt_idx)
                    pick = alt_idx(randi(numel(alt_idx)));
                    stim = pool(pick);
                    break;
                end
            end
        end

        % 3) Avoid long runs of same feeling
        if last_n_are(block_trials, 'face_feeling', stim.feeling, maxRunFeeling)
            altFeels = {'h','n','s'};
            altFeels(strcmp(altFeels, stim.feeling)) = [];  % remove current feeling
            for af = 1:numel(altFeels)
                alt_idx = find_stim_3s(pool, stim.gender, altFeels{af}, stim.age);
                if ~isempty(alt_idx)
                    pick = alt_idx(randi(numel(alt_idx)));
                    stim = pool(pick);
                    break;
                end
            end
        end

        % ---- Finalize chosen stimulus for this trial and update counters ----
        pool(pick) = [];   % remove used image so it is never reused

        % Gender counters
        if strcmp(stim.gender,'f')
            cntF = cntF + 1;
        else
            cntM = cntM + 1;
        end

        % Age counter
        cntAge.(stim.age) = cntAge.(stim.age) + 1;

        % feeling × congruency cell counter
        cellCount.(stim.feeling).(currCong) = ...
            cellCount.(stim.feeling).(currCong) + 1;

        % Word logic (which German word to show on the face)
        [word, congLabel] = choose_word_for_face(stim.feeling, currCong);

        % Store trial in block_trials
        block_trials(t).block        = b;
        block_trials(t).trial        = t;
        block_trials(t).img          = stim.fname;
        block_trials(t).age          = stim.age;
        block_trials(t).gender       = stim.gender;
        block_trials(t).face_feeling = stim.feeling;

        block_trials(t).prevCong     = prevCong;
        block_trials(t).currCong     = currCong;
        block_trials(t).congruency   = congLabel;   % 'congruent' / 'incongruent' / 'slightlyIncongruent'
        block_trials(t).word         = word;        % actual German label shown
    end

    % ---------- Assign ISIs for this block ----------
    isi_vec = makeISI_block(nTrialsPerBlock, s.isi_min, s.isi_max, s.isi_mean);
    for t = 1:nTrialsPerBlock
        block_trials(t).isi = isi_vec(t);
    end

    % Append this block’s trials to the global sequence
    exp_seq = [exp_seq; block_trials(:)]; %#ok<AGROW>
end

end  % end of main seq_EFST_3s function



%% NEW, MORE ROBUST 3-STATE CONGRUENCY BUILDER
function cong_seq = build_cong_sequence_EFST3(nTrials)
% build_cong_sequence_EFST3
%
% Build an approximately balanced sequence of length nTrials using 3 states:
%   'C'  – congruent
%   'I'  – incongruent
%   'IS' – slightly incongruent
%
% Properties:
%   • Counts of C / I / IS differ at most by 1 (balanced frequencies)
%   • Transitions between states are "as even as possible" (3×3 matrix)
%   • Long runs of the same state are penalized (limited by maxRunCongruency)

if nTrials < 3
    error('Need at least 3 trials.');
end

% ---- 1) target counts for C / I / IS (balanced ±1) ----
base = floor(nTrials / 3);
rem  = mod(nTrials, 3);

nC  = base;
nI  = base;
nIS = base;

if rem >= 1
    nC = nC + 1;
end
if rem == 2
    nI = nI + 1;
end
% nIS implicitly gets the remaining trials because nC + nI + nIS = nTrials

% Create initial label multiset
labels = [repmat({'C'}, 1, nC), ...
          repmat({'I'}, 1, nI), ...
          repmat({'IS'},1, nIS)];

maxAttempts      = 5000;
maxRunCongruency = 6;   % maximum allowed length of same state run
bestSeq   = [];
bestScore = Inf;

for attempt = 1:maxAttempts

    % random permutation of the label pool
    perm = labels(randperm(nTrials));

    % ---- check max run length of congruency states ----
    runLen = 1;
    badRun = false;
    for t = 2:nTrials
        if strcmp(perm{t}, perm{t-1})
            runLen = runLen + 1;
            if runLen > maxRunCongruency
                badRun = true;
                break;
            end
        else
            runLen = 1;
        end
    end
    if badRun
        % too long run → discard this permutation
        continue;
    end

    % ---- compute transition counts (3×3 matrix) ----
    % index mapping: C=1, I=2, IS=3
    transCounts = zeros(3,3);
    for t = 2:nTrials
        a = stateIdx(perm{t-1});
        b = stateIdx(perm{t});
        transCounts(a,b) = transCounts(a,b) + 1;
    end

    % we want transitions as equal as possible across the 9 entries
    vals  = transCounts(:);
    score = max(vals) - min(vals);  % the smaller the better

    % keep the best sequence seen so far
    if score < bestScore
        bestScore = score;
        bestSeq   = perm;
        % If transitions are almost perfectly balanced, we can stop early
        if score <= 2
            cong_seq = bestSeq;
            return;
        end
    end
end

% If the loop finishes, we fall back to the best sequence we found
if ~isempty(bestSeq)
    warning('build_cong_sequence_EFST3:UsingBestApprox', ...
        'Could not get perfectly balanced transitions; using best approx (range = %.1f).', bestScore);
    cong_seq = bestSeq;
else
    error('build_cong_sequence_EFST3:Failed', ...
        'Could not build a 3-state sequence after %d attempts.', maxAttempts);
end

end


function idx = stateIdx(label)
% stateIdx – map state labels 'C' / 'I' / 'IS' to indices 1 / 2 / 3
switch label
    case 'C'
        idx = 1;
    case 'I'
        idx = 2;
    case 'IS'
        idx = 3;
    otherwise
        error('Unknown state label %s', label);
end
end


%% ======================================================================
% WORD SELECTION
% ======================================================================

function [word, congLabel] = choose_word_for_face(faceFeeling, congCode)
% choose_word_for_face
%
% Inputs:
%   faceFeeling – 'h' | 'n' | 's'  (happy | neutral | sad face)
%   congCode    – 'C' | 'I' | 'IS' (congruent, incongruent, slightly incongruent)
%
% Outputs:
%   word      – German word to draw on the face ("GLÜCKLICH", "NEUTRAL", "TRAURIG")
%   congLabel – textual label for the trial ("congruent", "incongruent", "slightlyIncongruent")

happyWord   = 'GLÜCKLICH';
neutralWord = 'NEUTRAL';
sadWord     = 'TRAURIG';

switch faceFeeling

    case 'h'  % happy face
        switch congCode
            case 'C'
                word = happyWord;
                congLabel = 'congruent';
            case 'I'
                word = sadWord;
                congLabel = 'incongruent';
            case 'IS'
                word = neutralWord;
                congLabel = 'slightlyIncongruent';
        end

    case 's'  % sad face
        switch congCode
            case 'C'
                word = sadWord;
                congLabel = 'congruent';
            case 'I'
                word = happyWord;
                congLabel = 'incongruent';
            case 'IS'
                word = neutralWord;
                congLabel = 'slightlyIncongruent';
        end

    case 'n'  % neutral face
        switch congCode
            case 'C'
                word = neutralWord;
                congLabel = 'congruent';
            case 'I'
                word = happyWord;
                congLabel = 'incongruent';
            case 'IS'
                word = sadWord;
                congLabel = 'slightlyIncongruent';
        end
end
end


%% ======================================================================
% ISI HELPER
% ======================================================================

function isi = makeISI_block(n, lo, hi, target)
% makeISI_block – generate n ISIs in [lo, hi] with mean ≈ target
%
% Simple approach:
%   1) draw n uniform random values in [lo, hi]
%   2) shift them all so the mean matches "target"
%   3) clip any that fell outside [lo, hi] after shifting

isi = lo + (hi - lo) * rand(n,1);   % initial uniform distribution
curMean = mean(isi);                % current mean
isi = isi + (target - curMean);     % shift to match target mean
isi(isi < lo) = lo;                 % clip to lower bound
isi(isi > hi) = hi;                 % clip to upper bound
end


%% ======================================================================
% SMALL HELPERS
% ======================================================================

function out = iff(cond,a,b)
% iff – inline if (ternary) helper: out = a if cond, otherwise b
if cond, out = a; else, out = b; end
end


function tf = last_n_are(trials, fieldname, value, n)
% last_n_are – check if the last n entries of "trials" have field == value
%
% Used for anti-run constraints (gender, age, feeling).

if numel(trials) < n
    tf = false;
    return;
end

tf = true;
for k = numel(trials)-n+1 : numel(trials)
    if ~strcmp(trials(k).(fieldname), value)
        tf = false;
        return;
    end
end
end


function idx = find_stim_3s(pool, gender, feeling, age)
% find_stim_3s – find indices in pool matching given gender, feeling, age
idx = find( strcmp({pool.gender},  gender) & ...
            strcmp({pool.feeling}, feeling) & ...
            strcmp({pool.age},     age) );
end
