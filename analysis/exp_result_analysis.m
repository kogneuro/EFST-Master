function [result, report] = exp_result_analysis(exp_data)
% EXP_RESULT_ANALYSIS
% Build trial-wise result struct and detailed summary report.
%
% Assumes exp_data has:
%   exp_data.seq       -> 1xN struct with fields incl. face_feeling, prevCong, currCong, congruency
%   exp_data.resp1_btn -> Nx1 numeric (1 = sad, 2 = neutral, 3 = happy)
%   exp_data.resp1_t   -> Nx1 numeric response times

    %% --- Build RESULT struct ---
    N = numel(exp_data.seq);

    % Map button responses to 's' / 'n' / 'h'
    resp_map    = {'s','n','h'};          % 1 -> 's', 2 -> 'n', 3 -> 'h'
    resp_btn    = exp_data.resp1_btn(:);  % ensure column vector
    resp_labels = cell(N,1);

    for k = 1:N
        idx = resp_btn(k);
        if ~isnan(idx) && idx >= 1 && idx <= 3
            resp_labels{k} = resp_map{idx};
        else
            resp_labels{k} = '';  % fallback
        end
    end

    % Create unified struct array: one element per stimulus/trial
    result = exp_data.seq;   % start with all fields from seq

    for k = 1:N
        result(k).resp1_btn   = resp_btn(k);          % numeric button (1/2/3)
        result(k).resp1_label = resp_labels{k};       % 's','n','h'
        result(k).resp1_t     = exp_data.resp1_t(k);  % response time
        result(k).correct     = strcmp(result(k).resp1_label, result(k).face_feeling);
    end

    %% --- REPORT struct ---
    report = struct();

    %% A) Overall correctness
    correct_vec = [result.correct]';
    report.total_correct    = sum(correct_vec);
    report.total_incorrect  = sum(~correct_vec);
    report.total_trials     = N;
    report.accuracy_percent = (report.total_correct / N) * 100;

    %% B) Overall RT metrics (all trials)
    RT = exp_data.resp1_t(:);

    report.rt_min  = min(RT);
    report.rt_max  = max(RT);
    report.rt_mean = mean(RT);
    report.rt_std  = std(RT);
    report.rt_p25  = prctile(RT, 25);
    report.rt_p75  = prctile(RT, 75);

    %% C) prevCong × currCong combinations (C/I/IS → 3×3 = 9)
    levels = {'C','I','IS'};

    report.combinations = struct();

    comb_list = struct( ...
        'prevCong',    {}, ...
        'currCong',    {}, ...
        'n_correct',   [], ...
        'n_incorrect', [], ...
        'rt_min',      [], ...
        'rt_max',      [], ...
        'rt_mean',     [], ...
        'rt_std',      [], ...
        'rt_p25',      [], ...
        'rt_p75',      [] );

    idx_comb = 0;

    for a = 1:3
        for b = 1:3

            prev = levels{a};
            curr = levels{b};

            % Trials belonging to this combination
            mask   = strcmp({result.prevCong}', prev) & strcmp({result.currCong}', curr);
            trials = result(mask);

            if isempty(trials)
                stats.n_correct   = 0;
                stats.n_incorrect = 0;
                stats.rt_min  = NaN;
                stats.rt_max  = NaN;
                stats.rt_mean = NaN;
                stats.rt_std  = NaN;
                stats.rt_p25  = NaN;
                stats.rt_p75  = NaN;
            else
                corr_vec = [trials.correct]';
                rt_vec   = [trials.resp1_t]';

                stats.n_correct   = sum(corr_vec);
                stats.n_incorrect = sum(~corr_vec);

                stats.rt_min  = min(rt_vec);
                stats.rt_max  = max(rt_vec);
                stats.rt_mean = mean(rt_vec);
                stats.rt_std  = std(rt_vec);
                stats.rt_p25  = prctile(rt_vec, 25);
                stats.rt_p75  = prctile(rt_vec, 75);
            end

            % store under e.g. report.combinations.C_I
            keyname = sprintf('%s_%s', prev, curr);
            report.combinations.(keyname) = stats;

            % append to combined list
            idx_comb = idx_comb + 1;
            comb_list(idx_comb).prevCong    = prev;
            comb_list(idx_comb).currCong    = curr;
            comb_list(idx_comb).n_correct   = stats.n_correct;
            comb_list(idx_comb).n_incorrect = stats.n_incorrect;
            comb_list(idx_comb).rt_min      = stats.rt_min;
            comb_list(idx_comb).rt_max      = stats.rt_max;
            comb_list(idx_comb).rt_mean     = stats.rt_mean;
            comb_list(idx_comb).rt_std      = stats.rt_std;
            comb_list(idx_comb).rt_p25      = stats.rt_p25;
            comb_list(idx_comb).rt_p75      = stats.rt_p75;

        end
    end

    report.combination_list  = comb_list;
    report.combination_table = struct2table(comb_list);

    %% D) Feeling (s/n/h) × Congruency (congruent/incongruent)
    % Treat 'slightlyIncongruent' as 'incongruent'

    % Normalize congruency into: 'congruent' or 'incongruent'
    normCong = cell(N,1);
    for k = 1:N
        c = result(k).congruency;
        if strcmpi(c,'incongruent') || strcmpi(c,'slightlyIncongruent')
            normCong{k} = 'incongruent';
        else
            normCong{k} = 'congruent';
        end
        result(k).normCongruency = normCong{k}; % store inside result struct
    end

    feeling_codes  = {'s','n','h'};
    feeling_labels = {'sad','neutral','happy'};   % just for readability
    cong_levels    = {'congruent','incongruent'};

    feel_list = struct( ...
        'feeling_code',  {}, ...
        'feeling_label', {}, ...
        'congruency',    {}, ...
        'n_correct',     [], ...
        'n_incorrect',   [], ...
        'rt_min',        [], ...
        'rt_max',        [], ...
        'rt_mean',       [], ...
        'rt_std',        [], ...
        'rt_p25',        [], ...
        'rt_p75',        [] );

    idx_f = 0;

    for f = 1:numel(feeling_codes)
        for c = 1:numel(cong_levels)

            fc = feeling_codes{f};
            fl = feeling_labels{f};
            cg = cong_levels{c};

            % Select trials with this feeling + (normalized) congruency
            mask   = strcmp({result.face_feeling}', fc) & strcmp({result.normCongruency}', cg);
            trials = result(mask);

            if isempty(trials)
                stats2.n_correct   = 0;
                stats2.n_incorrect = 0;
                stats2.rt_min  = NaN;
                stats2.rt_max  = NaN;
                stats2.rt_mean = NaN;
                stats2.rt_std  = NaN;
                stats2.rt_p25  = NaN;
                stats2.rt_p75  = NaN;
            else
                corr_vec = [trials.correct]';
                rt_vec   = [trials.resp1_t]';

                stats2.n_correct   = sum(corr_vec);
                stats2.n_incorrect = sum(~corr_vec);

                stats2.rt_min  = min(rt_vec);
                stats2.rt_max  = max(rt_vec);
                stats2.rt_mean = mean(rt_vec);
                stats2.rt_std  = std(rt_vec);
                stats2.rt_p25  = prctile(rt_vec, 25);
                stats2.rt_p75  = prctile(rt_vec, 75);
            end

            idx_f = idx_f + 1;
            feel_list(idx_f).feeling_code  = fc;
            feel_list(idx_f).feeling_label = fl;
            feel_list(idx_f).congruency    = cg;
            feel_list(idx_f).n_correct     = stats2.n_correct;
            feel_list(idx_f).n_incorrect   = stats2.n_incorrect;
            feel_list(idx_f).rt_min        = stats2.rt_min;
            feel_list(idx_f).rt_max        = stats2.rt_max;
            feel_list(idx_f).rt_mean       = stats2.rt_mean;
            feel_list(idx_f).rt_std        = stats2.rt_std;
            feel_list(idx_f).rt_p25        = stats2.rt_p25;
            feel_list(idx_f).rt_p75        = stats2.rt_p75;

        end
    end

    report.feeling_list  = feel_list;
    report.feeling_table = struct2table(feel_list);

end
