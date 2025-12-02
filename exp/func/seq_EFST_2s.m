function [exp_seq, s] = seq_EFST_2s(s)
% seq_EFST_2s  –  Sequence builder for 2-emotion EFST (happy / sad)
% -------------------------------------------------------------------------
% faces:       feeling = h / s   (happy / sad)
% congruency:  C / I             (congruent / incongruent)
%
% Design goals:
%   • No image reuse across the whole experiment.
%   • Per block, balance:
%       - gender (F / M)
%       - age (y / m / o)
%       - feeling × congruency cells: (h,C), (h,I), (s,C), (s,I)
%         ==> each cell gets ≈ same number of trials.
%   • Congruency sequence C/I is approx-balanced:
%       - roughly 50/50 C vs I
%       - transitions (C→C, C→I, I→C, I→I) approximately balanced
%       - no excessively long runs of the same state.
%   • ISI per trial is in [s.isi_min, s.isi_max] with mean ≈ s.isi_mean.
%
% Output:
%   exp_seq : 1D struct array, length = blocks * trials_per_block
%             Fields per trial (block_trials(t)):
%               .block        – block index
%               .trial        – trial index within block
%               .img          – filename of chosen face
%               .age          – 'y','m','o'
%               .gender       – 'f','m'
%               .face_feeling – 'h' or 's'
%               .prevCong     – previous congruency ('C' or 'I')
%               .currCong     – current congruency ('C' or 'I')
%               .congruency   – 'congruent' or 'incongruent'
%               .word         – 'GLÜCKLICH' or 'TRAURIG'
%               .isi          – ISI in seconds for this trial
%
% Alp, 2025
% -------------------------------------------------------------------------

%% 1) RNG – store RNG state for reproducibility
% This lets you later reconstruct the exact same sequence if needed.
s.rng_state = set_rng_state(s);
exp_seq = [];  % will collect all block_trials here

%% 2) Load *.jpg stimuli and build GLOBAL pool
% We load all images once and consume them across all blocks,
% so no image is ever re-used.

files_all = dir(fullfile(s.stim_dir, '*.jpg'));   % all .jpg in the folder
keep_idx  = endsWith({files_all.name}, '.jpg');   % trivial filter here
files     = files_all(keep_idx);
clear files_all keep_idx

if isempty(files)
    error('seq_EFST_2s:NoStim', ...
        'No *.jpg files found in %s', s.stim_dir);
end

% The pool holds all usable stimuli with parsed metadata.
pool = [];
for iFile = 1:numel(files)
    fname = files(iFile).name;
    parts = split(fname, '_');          % expect subj_age_gender_feeling.jpg

    % We need at least 4 parts to interpret subj/age/gender/feeling
    if numel(parts) < 4
        continue
    end

    subj    = parts{1};                 % subject ID (string)
    age     = lower(parts{2});         % 'y','m','o'
    gender  = lower(parts{3});         % 'f','m'
    feeling = lower(parts{4});         % 'h','s' + possible suffix

    % Remove endings like ".jpg" if they appear
    feeling = erase(feeling, '.jpg');

    % Only accept legal categories
    if ~ismember(age, {'y','m','o'}), continue, end
    if ~ismember(gender, {'f','m'}),  continue, end
    if ~ismember(feeling, {'h','s'}), continue, end

    % Store in pool
    pool(end+1).fname   = fname;   %#ok<AGROW>
    pool(end).subj      = subj;
    pool(end).age       = age;
    pool(end).gender    = gender;
    pool(end).feeling   = feeling;
end

if isempty(pool)
    error('seq_EFST_2s:NoValidStim', ...
        'Found *.jpg but none parsed as (age,gender,feeling h/s).');
end

%% 3) Basic counts (blocks, trials, total needed images)
nTrialsPerBlock = s.trials_per_block;
nBlocks         = s.blocks;
nTotalTrials    = nTrialsPerBlock * nBlocks;

% Ensure we have enough *distinct* images for all trials
if numel(pool) < nTotalTrials
    error('seq_EFST_2s:TooFewImages', ...
        'Need %d images (no repeats) but found %d.', nTotalTrials, numel(pool));
end

%% 4) Define target distributions: age + feeling×congruency
% -------------------------------------------------------------------------
% Age categories (for age balancing)
ageCats     = {'y','m','o'};     nAge   = numel(ageCats);

% Feelings (h/s) and congruency (C/I) – we will balance per cell
feelingCats = {'h','s'};         nFeel  = numel(feelingCats);
congCats    = {'C','I'};         nCong  = numel(congCats);

% ---- 4a) Age targets per BLOCK -----------------------------------------
% We split trials per block as evenly as possible across y/m/o.
baseAge = floor(nTrialsPerBlock / nAge);   % minimal number per age
remAge  = mod(nTrialsPerBlock, nAge);      % leftover trials

ageTargets = zeros(nBlocks, nAge);         % block × age
for b = 1:nBlocks
    ageTargets(b,:) = baseAge;
    if remAge > 0
        % Distribute leftovers across blocks in a rotating fashion
        idxExtra = mod(b-1, nAge) + 1;
        ageTargets(b, idxExtra) = ageTargets(b, idxExtra) + remAge;
    end
end

% ---- 4b) Feeling×Congruency cell targets per BLOCK ---------------------
% Cells are: (h,C), (h,I), (s,C), (s,I).
% We want them as equal as possible, but nTrialsPerBlock may not be
% divisible by 4 → we allow ±1 difference.

baseCell = floor(nTrialsPerBlock / (nFeel * nCong));  % minimal per cell
remCell  = mod(nTrialsPerBlock,  nFeel * nCong);      % leftover trials

% cellTarget.(feeling).(cong) = target count per cell
cellTarget = struct();
for f = 1:nFeel
    ff = feelingCats{f};           % 'h' or 's'
    for c = 1:nCong
        cc = congCats{c};          % 'C' or 'I'
        cellTarget.(ff).(cc) = baseCell;
    end
end

% Distribute the leftover trials across cells in a fixed order:
% (h,C) → (h,I) → (s,C) → (s,I) → then repeat if remCell>4.
if remCell > 0
    idx = 1;
    for r = 1:remCell
        f = mod(idx-1, nFeel) + 1;
        c = floor((idx-1) / nFeel) + 1;
        ff = feelingCats{f};
        cc = congCats{c};
        cellTarget.(ff).(cc) = cellTarget.(ff).(cc) + 1;
        idx = idx + 1;
    end
end
% Now sum(cellTarget.h.C, h.I, s.C, s.I) = nTrialsPerBlock.

%% 5) Build each BLOCK independently
for b = 1:nBlocks

    % ---- 5a) Congruency sequence (C/I) for this block -------------------
    % 2-state builder that:
    %   • approx balances number of C vs I
    %   • approx balances transitions between them
    cong_seq = build_cong_sequence_EFST2(nTrialsPerBlock);

    % ---- 5b) Gender target (F/M) ---------------------------------------
    if mod(nTrialsPerBlock,2) ~= 0
        error('seq_EFST_2s:EvenTrialsNeeded', ...
            'trials_per_block must be even to balance gender.');
    end
    target_gender = nTrialsPerBlock / 2;   % half F, half M

    % ---- 5c) Counters for gender, age, cells ---------------------------
    cntF = 0; cntM = 0;                   % #female, #male used so far

    % Age counters per block:
    cntAge = struct('y',0,'m',0,'o',0);   % how many y/m/o used in this block
    tgtAge = struct('y', ageTargets(b,1), ...
                    'm', ageTargets(b,2), ...
                    'o', ageTargets(b,3));

    % Feeling×congruency counters:
    cellCount = struct();
    for f = 1:nFeel
        ff = feelingCats{f};     % 'h' or 's'
        cellCount.(ff).C = 0;    % count of (ff,C)
        cellCount.(ff).I = 0;    % count of (ff,I)
    end

    % Placeholder for all trials in this block
    block_trials = struct([]);

    % ---- 5d) Anti-run settings -----------------------------------------
    maxRunGender  = 3;   % avoid >3 same gender in a row
    maxRunAge     = 3;   % optional: avoid long age runs
    maxRunFeeling = 3;   % avoid long 'h' or 's' runs

    %% 5e) Trial loop within this block
    for t = 1:nTrialsPerBlock

        currCong = cong_seq{t};    % 'C' or 'I' for current trial

        % previous congruency: used for logging pattern (prev→curr)
        if t == 1
            prevCong = currCong;   % for the first trial, define prev = curr
        else
            prevCong = cong_seq{t-1};
        end

        % -------- choose GENDER (F/M) to balance block-wise -------------
        if cntF < target_gender && cntM < target_gender
            % both under target → free to pick random F/M
            wantGender = iff(rand < 0.5, 'f', 'm');
        elseif cntF < target_gender
            wantGender = 'f';
        else
            wantGender = 'm';
        end

        % -------- choose FEELING (h/s) based on cell balancing ----------
        % We look at the column for currCong (C or I) and see which
        % feelings are still under their cellTarget.(feeling).(currCong).
        underCell = false(1,nFeel);
        for f = 1:nFeel
            ff = feelingCats{f};   % 'h' or 's'
            if cellCount.(ff).(currCong) < cellTarget.(ff).(currCong)
                underCell(f) = true;
            end
        end

        candFeelIdx = find(underCell);    % index/indices of underfilled cells
        if isempty(candFeelIdx)
            % In practice, this should only happen near the very end of
            % the block when the column is already exactly at target.
            % Then we allow any feeling (we’ve already hit the target).
            candFeelIdx = 1:nFeel;
        end
        % Random pick among the allowed feelings
        wantFeeling = feelingCats{ candFeelIdx(randi(numel(candFeelIdx))) };

        % -------- choose AGE (y/m/o) to approach age targets -------------
        underAge = [cntAge.y < tgtAge.y, ...
                    cntAge.m < tgtAge.m, ...
                    cntAge.o < tgtAge.o];
        idxUA = find(underAge);
        if isempty(idxUA)
            % all ages at or above target (should only occur at the tail)
            wantAge = ageCats{randi(nAge)};
        elseif numel(idxUA) == 1
            % only one age category still under target
            wantAge = ageCats{idxUA};
        else
            % multiple age categories under target → random among them
            wantAge = ageCats{ idxUA(randi(numel(idxUA))) };
        end

        % -------- first pass: find matching stimulus in the pool ---------
        cand_idx = find(strcmp({pool.gender},  wantGender) & ...
                        strcmp({pool.feeling}, wantFeeling) & ...
                        strcmp({pool.age},     wantAge));
        if isempty(cand_idx)
            error('seq_EFST_2s:MissingStratum', ...
                'No image left for block %d trial %d (gender=%s, feeling=%s, age=%s)', ...
                b, t, wantGender, wantFeeling, wantAge);
        end

        % Choose one candidate at random
        pick = cand_idx(randi(numel(cand_idx)));
        stim = pool(pick);

        %% =========================
        % ANTI-RUN SECTION
        % =========================

        % ---- 1) Gender anti-run ----------------------------------------
        if last_n_are(block_trials, 'gender', stim.gender, maxRunGender)
            % If we already have 3 same-gender in a row, try to switch to
            % the other gender (but keep feeling + age if possible).
            otherGender = iff(strcmp(stim.gender,'f'),'m','f');
            alt_idx = find_stim_2s(pool, otherGender, wantFeeling, wantAge);
            if ~isempty(alt_idx)
                pick = alt_idx(randi(numel(alt_idx)));
                stim = pool(pick);
            end
        end

        % ---- 2) Age anti-run -------------------------------------------
        if last_n_are(block_trials, 'age', stim.age, maxRunAge)
            % If same age repeats too often, try other ages while keeping
            % gender + feeling fixed.
            altAges = {'y','m','o'};
            altAges(strcmp(altAges, stim.age)) = [];   % remove current age
            for aa = 1:numel(altAges)
                alt_idx = find_stim_2s(pool, stim.gender, wantFeeling, altAges{aa});
                if ~isempty(alt_idx)
                    pick = alt_idx(randi(numel(alt_idx)));
                    stim = pool(pick);
                    break;
                end
            end
        end

        % ---- 3) Feeling anti-run ---------------------------------------
        if last_n_are(block_trials, 'face_feeling', stim.feeling, maxRunFeeling)
            % If feeling repeats too long, try switching to the other
            % feeling while keeping gender + age fixed as much as possible.
            altFeels = {'h','s'};
            altFeels(strcmp(altFeels, stim.feeling)) = []; % remove current
            for af = 1:numel(altFeels)
                alt_idx = find_stim_2s(pool, stim.gender, altFeels{af}, stim.age);
                if ~isempty(alt_idx)
                    pick = alt_idx(randi(numel(alt_idx)));
                    stim = pool(pick);
                    break;
                end
            end
        end

        %% ---- finalize chosen stimulus & update counters ----------------
        % Remove the chosen stimulus from the pool so it can't be reused.
        pool(pick) = [];

        % Update gender counters
        if strcmp(stim.gender,'f')
            cntF = cntF + 1;
        else
            cntM = cntM + 1;
        end

        % Update age counts
        cntAge.(stim.age) = cntAge.(stim.age) + 1;

        % Update feeling×congruency cell counts
        cellCount.(stim.feeling).(currCong) = ...
            cellCount.(stim.feeling).(currCong) + 1;

        % ---- Map face feeling + congruency → word ----------------------
        % For happy faces:  C → "GLÜCKLICH", I → "TRAURIG"
        % For sad   faces:  C → "TRAURIG",   I → "GLÜCKLICH"

        if strcmp(stim.feeling,'h')
            faceWord   = 'GLÜCKLICH';
            incongWord = 'TRAURIG';
        else % 's'
            faceWord   = 'TRAURIG';
            incongWord = 'GLÜCKLICH';
        end

        if strcmp(currCong,'C')
            word      = faceWord;
            congLabel = 'congruent';
        else
            word      = incongWord;
            congLabel = 'incongruent';
        end

        % ---- Store trial in block_trials --------------------------------
        block_trials(t).block        = b;
        block_trials(t).trial        = t;
        block_trials(t).img          = stim.fname;
        block_trials(t).age          = stim.age;
        block_trials(t).gender       = stim.gender;
        block_trials(t).face_feeling = stim.feeling;

        block_trials(t).prevCong     = prevCong;
        block_trials(t).currCong     = currCong;
        block_trials(t).congruency   = congLabel;
        block_trials(t).word         = word;
    end  % end trial loop within block

    %% 5f) Assign ISIs per trial in this block
    % We generate random ISIs within [isi_min, isi_max] and then shift them
    % so that the average ≈ s.isi_mean, then clip to bounds.
    isi_vec = makeISI_block(nTrialsPerBlock, s.isi_min, s.isi_max, s.isi_mean);
    for t = 1:nTrialsPerBlock
        block_trials(t).isi = isi_vec(t);
    end

    % Append this block to global sequence
    exp_seq = [exp_seq; block_trials(:)]; %#ok<AGROW>
end  % end block loop

end  % ================== end of seq_EFST_2s ==============================


%% ======================================================================
% 2-STATE CONGRUENCY BUILDER (C / I) – approx-balanced
% =======================================================================
function cong_seq = build_cong_sequence_EFST2(nTrials)
% build_cong_sequence_EFST2
% -------------------------------------------------------------------------
% Builds a sequence of length nTrials with labels 'C' / 'I':
%   • approx 50/50 C and I (±1 if nTrials is odd)
%   • transitions CC, CI, IC, II approximately balanced
%   • runs of the same label limited (maxRun)
%   • uses a best-of-many-random-permutations approach

if nTrials < 2
    error('build_cong_sequence_EFST2:NeedAtLeast2', ...
        'Need at least 2 trials.');
end

% target counts for C and I (±1 if odd)
base = floor(nTrials / 2);
rem  = mod(nTrials, 2);

nC = base;
nI = base;
if rem == 1
    % Simple rule: if odd, put the extra trial in 'C'
    % (could randomize if you want)
    nC = nC + 1;
end

% Build label list, e.g. {'C','C','C','I','I','I',...}
labels = [repmat({'C'},1,nC), repmat({'I'},1,nI)];

maxAttempts      = 5000;  % how many random permutations to try
maxRun           = 6;     % max allowed length of same-state run
bestSeq   = [];
bestScore = Inf;          % we minimize "score" = max(trans) - min(trans)

for attempt = 1:maxAttempts

    % Randomly permute labels
    perm = labels(randperm(nTrials));

    % ---- 1) Check run length constraint --------------------------------
    runLen = 1;
    badRun = false;
    for t = 2:nTrials
        if strcmp(perm{t}, perm{t-1})
            runLen = runLen + 1;
            if runLen > maxRun
                badRun = true;
                break;
            end
        else
            runLen = 1;
        end
    end
    if badRun
        % This permutation has too long a run, skip it
        continue;
    end

    % ---- 2) Compute transition counts (2×2 matrix) ---------------------
    % Rows = previous state, columns = current state:
    %   row 1 = C, row 2 = I
    %   col 1 = C, col 2 = I
    %
    % So:
    %   transCounts(1,1) = # of C→C,
    %   transCounts(1,2) = # of C→I,
    %   transCounts(2,1) = # of I→C,
    %   transCounts(2,2) = # of I→I.

    transCounts = zeros(2,2);
    for t = 2:nTrials
        % prev index: C→1, I→2
        a = strcmp(perm{t-1},'C') + 1;
        % current index: C→1, I→2
        b = strcmp(perm{t},  'C') + 1;
        transCounts(a,b) = transCounts(a,b) + 1;
    end

    vals  = transCounts(:);          % flatten to vector of 4 elements
    score = max(vals) - min(vals);   % smaller range = more balanced

    % ---- 3) Keep best sequence so far ----------------------------------
    if score < bestScore
        bestScore = score;
        bestSeq   = perm;
        % Early stop: if transitions already very balanced (range ≤ 1)
        if score <= 1
            cong_seq = bestSeq;
            return;
        end
    end
end

% If we exit the loop, we didn't hit a perfect/perfect solution,
% but we still have some "best" candidate.
if ~isempty(bestSeq)
    warning('build_cong_sequence_EFST2:UsingBestApprox', ...
        'Could not get perfectly balanced transitions; using best approx (range = %.1f).', bestScore);
    cong_seq = bestSeq;
else
    error('build_cong_sequence_EFST2:Failed', ...
        'Could not build a 2-state sequence after %d attempts.', maxAttempts);
end

end  % end build_cong_sequence_EFST2


%% ======================================================================
% Small HELPER FUNCTIONS
% ======================================================================

function out = iff(cond,a,b)
% inline-if replacement: out = a if cond is true, otherwise b
if cond
    out = a;
else
    out = b;
end
end

function tf = last_n_are(trials, fieldname, value, n)
% last_n_are
% -------------------------------------------------------------------------
% Returns true if the last n elements of trials all have
%   trials(k).(fieldname) == value.
% If there are fewer than n trials so far, returns false.
%
% Used to detect long runs (e.g., gender/age/feeling) and then avoid them.

if numel(trials) < n
    tf = false;
    return;
end
tf = true;
for k = numel(trials)-n+1 : numel(trials)
    if ~isfield(trials(k),fieldname) || ~strcmp(trials(k).(fieldname), value)
        tf = false;
        return;
    end
end
end

function idx = find_stim_2s(pool, gender, feeling, age)
% find_stim_2s
% -------------------------------------------------------------------------
% Find indices in the pool that match EXACTLY:
%   gender, feeling, and age.
%
% Used in anti-run logic to search for alternative stimuli.

idx = find( strcmp({pool.gender},  gender) & ...
            strcmp({pool.feeling}, feeling) & ...
            strcmp({pool.age},     age) );
end

function isi = makeISI_block(n, lo, hi, target)
% makeISI_block
% -------------------------------------------------------------------------
% Create a vector of n ISIs:
%   - initially uniform in [lo,hi]
%   - then shifted so that the mean ≈ target
%   - then clipped to stay within [lo,hi]
%
% Used to approximate your desired ISI distribution per block.

% 1) random values in [lo, hi]
isi = lo + (hi - lo) * rand(n,1);

% 2) shift the whole vector so its mean becomes "target"
curMean = mean(isi);
isi = isi + (target - curMean);

% 3) clip to [lo, hi]
isi(isi < lo) = lo;
isi(isi > hi) = hi;
end
