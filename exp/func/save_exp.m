function save_exp(p_data, exp_data, s, file_name_end)
% save_exp_EFST
% ----------------------------------------------------------
% Generic saver for EFST-style experiments:
%   - Emotional Face Stroop (with img / face_feeling / age / gender_face)
%   - Color Stroop (no img/face fields; writes 'NA' in those columns)
%
% It saves:
%   1) .mat file with p_data, exp_data, s
%   2) a trial-wise .txt file with timing + behavioral info
%
% Alp, 2025
% ----------------------------------------------------------

%% ---------- base filename ----------
file_name = [s.file_prefix p_data.ID];

%% ---------- make sure participant dir exists ----------
if ~exist(p_data.dir, 'dir')
    mkdir('.', p_data.dir);
end

%% ---------- compute intervals (if not yet done) ----------
% (Safe to call again; if already computed, this overwrites with same values)
exp_data = intervals_EFST(s, exp_data);

d = exp_data;                    % shorter handle
date_str = datestr(now,'yyyy/mm/dd');

%% ---------- save .mat ----------
mat_file_tmp = fullfile(p_data.dir, [file_name '_data_' file_name_end '.mat']);

if exist(mat_file_tmp, 'file')
    % if file exists, don’t overwrite — add timestamp
    disp('MAT-file exists. Creating file with time-stamp to avoid overwrite.');
    save(fullfile(p_data.dir, ...
        [file_name '_data_' file_name_end '_' num2str(round(sum(100*clock))) '.mat']), ...
        'p_data', 'exp_data', 's');
else
    save(mat_file_tmp, 'p_data', 'exp_data', 's');
end

%% ---------- open TXT trial file ----------
txt_path = fullfile(p_data.dir, [file_name '_trials_' file_name_end '.txt']);
data_file = fopen(txt_path, 'a');

% ---------- header ----------
fprintf(data_file, ['ID\tage\tgender\tdate\tblock\ttrial\t' ...
                    'img\tface_feeling\tgender_face\tage_group\t' ...
                    'congruency\tword\t' ...
                    'isi_planned_ms\t' ...
                    'onset_fix\tonset_cue\tonset_stim\tonset_resp1\t' ...
                    'resp1\tresp1_t\tresp1_btn\t' ...
                    't_fix_cue_ms\t' ...
                    't_cue_stim_ms\t' ...
                    't_stim_resp_ms\t' ...
                    't_resp_fix_ms\t' ...
                    't_trial_fix_ms\n']);

%% ---------- loop over trials ----------
nTrials = numel(d.seq);
for i = 1:nTrials

    tr = d.seq(i);   % one trial in the sequence

    % ====== 1) FACE / IMAGE FIELDS (OPTIONAL) ======
    % EFST: img, face_feeling, gender, age exist
    % Color Stroop: they don't; we write 'NA'
    if isfield(tr, 'img')
        img_str = tr.img;
    else
        img_str = 'NA';
    end

    if isfield(tr, 'face_feeling')
        face_feel_str = tr.face_feeling;
    else
        face_feel_str = 'NA';
    end

    % careful: p_data also has gender; here we mean *face* gender
    if isfield(tr, 'gender')
        gender_face_str = tr.gender;
    else
        gender_face_str = 'NA';
    end

    if isfield(tr, 'age')
        age_group_str = tr.age;
    else
        age_group_str = 'NA';
    end

    % ====== 2) ONSETS (CELLS) ======
    onset_fix   = d.onset_fix{i,1};
    onset_cue   = d.onset_cue{i,1};
    onset_stim  = d.onset_stim{i,1};
    onset_resp1 = d.onset_resp1{i,1};

    % ====== 3) INTERVALS (MS) ======
    if isfield(d, 't_fix_cue')
        t_fix_cue   = d.t_fix_cue(i,1);
    else
        t_fix_cue   = NaN;
    end

    if isfield(d, 't_cue_stim')
        t_cue_stim  = d.t_cue_stim(i,1);
    else
        t_cue_stim  = NaN;
    end

    if isfield(d, 't_stim_resp')
        t_stim_resp = d.t_stim_resp(i,1);
    else
        t_stim_resp = NaN;
    end

    if isfield(d, 't_resp_fix')
        t_resp_fix  = d.t_resp_fix(i,1);
    else
        t_resp_fix  = NaN;
    end

    if isfield(d, 't_trial_fix')
        t_trial_fix = d.t_trial_fix(i,1);
    else
        t_trial_fix = NaN;
    end

    % ====== 4) PLANNED ISI (MS) ======
    if isfield(d, 'isi_planned')
        isi_planned_ms = d.isi_planned(i,1);
    elseif isfield(tr, 'isi')
        isi_planned_ms = tr.isi * 1000;
    else
        isi_planned_ms = NaN;
    end

    % ====== 5) RESPONSE FIELDS ======
    if isfield(d, 'resp1')
        resp1_val = d.resp1(i,1);
    else
        resp1_val = NaN;
    end

    if isfield(d, 'resp1_t')
        resp1_t = d.resp1_t(i,1);
    else
        resp1_t = NaN;
    end

    if isfield(d, 'resp1_btn')
        resp1_btn = d.resp1_btn(i,1);
    else
        resp1_btn = NaN;
    end

    % ====== 6) SAFE CONGRUENCY + WORD ======
    if isfield(tr, 'congruency')
        cong_str = tr.congruency;
    else
        cong_str = 'NA';
    end

    if isfield(tr, 'word')
        word_str = tr.word;
    else
        word_str = 'NA';
    end

    % ====== 7) WRITE ONE LINE ======
    fprintf(data_file, ...
        ['%s\t%s\t%s\t%s\t' ...                           % ID, age, gender, date
         '%.0f\t%.0f\t' ...                               % block, trial
         '%s\t%s\t%s\t%s\t' ...                           % img, face_feeling, gender_face, age_group
         '%s\t%s\t' ...                                   % congruency, word
         '%.3f\t' ...                                     % isi_planned_ms
         '%.6f\t%.6f\t%.6f\t%.6f\t' ...                   % onsets
         '%.0f\t%.4f\t%.0f\t' ...                         % resp1, resp1_t, resp1_btn
         '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n'], ...           % intervals
         p_data.ID, p_data.age, p_data.gender, date_str, ...
         tr.block, i, ...
         img_str, face_feel_str, gender_face_str, age_group_str, ...
         cong_str, word_str, ...
         isi_planned_ms, ...
         onset_fix, onset_cue, onset_stim, onset_resp1, ...
         resp1_val, resp1_t, resp1_btn, ...
         t_fix_cue, t_cue_stim, t_stim_resp, t_resp_fix, t_trial_fix);
end

fclose(data_file);

end
