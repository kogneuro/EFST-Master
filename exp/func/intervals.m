function data = intervals(s, data)
% intervals_EFST  " compute timing diagnostics for the EFST version
%
% This is the EFST-compatible version of Martin Grund's intervals.m.
% It assumes your run_exp_EFST saved, per trial:
%   data.onset_fix{i,1}
%   data.onset_cue{i,1}
%   data.onset_stim{i,1}
%   data.onset_resp1{i,1}
% and that the sequence is a struct array with a field .isi:
%   data.seq(i).isi      % in seconds, from seq_EFST
%
% It does NOT try to compute AO / DAQ intervals, because EFST doesn™t use them.
%
% Output fields (all in ms):
%   data.t_fix_cue(i,1)    " fixation " cue
%   data.t_cue_stim(i,1)   " cue " stimulus
%   data.t_stim_resp(i,1)  " stimulus " response screen
%   data.isi_planned(i,1)  " ISI you planned in the seq (3-5 s)
%   data.t_trial_fix(i,1)  " fixation (this trial) " fixation (next trial), if possible
%   data.t_resp_fix(i,1)   " response screen onset " next fixation, if possible
%
% Alp / 2025

% how many trials do we have?
nTrials = numel(data.seq);

for i = 1:nTrials

    % ---------------------------
    % 1) fixation " cue
    % ---------------------------
    % both exist for every trial
    data.t_fix_cue(i,1) = (data.onset_cue{i,1}  - data.onset_fix{i,1})  * 1000;

    % ---------------------------
    % 2) cue " stimulus
    % ---------------------------
    data.t_cue_stim(i,1) = (data.onset_stim{i,1} - data.onset_cue{i,1}) * 1000;

    % ---------------------------
    % 3) stimulus " response screen
    %    (in your run_exp_EFST this is basically 0"a few ms,
    %     because you flip to response right after stim)
    % ---------------------------
    data.t_stim_resp(i,1) = (data.onset_resp1{i,1} - data.onset_stim{i,1}) * 1000;

    % ---------------------------
    % 4) planned ISI (from seq)
    %    this is what YOU asked for in the sequence: 3-5 s, mean 4 s
    % ---------------------------
    if isfield(data.seq(i), 'isi')
        data.isi_planned(i,1) = data.seq(i).isi * 1000;   % convert sec " ms
    else
        data.isi_planned(i,1) = NaN;
    end

    % ---------------------------
    % 5) trial-to-trial intervals
    %    we can only compute this if there IS a next trial,
    %    because we need the next trial's fixation time
    % ---------------------------
    if i < nTrials
        % fixation (this trial) " fixation (next trial)
        data.t_trial_fix(i,1) = (data.onset_fix{i+1,1} - data.onset_fix{i,1}) * 1000;

        % response onset (this trial) " fixation (next trial)
        % this tells you how long the tail of the ISI actually was
        data.t_resp_fix(i,1)  = (data.onset_fix{i+1,1} - data.onset_resp1{i,1}) * 1000;
    else
        % for the last trial in block we don't have œnext fixation?
        data.t_trial_fix(i,1) = NaN;
        data.t_resp_fix(i,1)  = NaN;
    end
end
end
