function [EEG,badtrials,EEGbu] = av101_findbaddata_postica(EEG,epochlims,eventcodes,resting)
% av101_findbaddata_preica
% function for identifying bad channels and trials to be removed prior to
% ica.
% inputs:
% EEG = eeglab structure

% epochlims = (in seconds) upper and lower limits for epoching.
% eventcodes = cell containing the string codes for each condition.
%%%% leave limits and codes blank if not planning to prune the data of
%%%% trials before ica (eg wavelet ica).
% resting = 0/1, not-resting/resting state
%
% outputs:
% EEG = eeglab structure with indices from processing added and data
% rejected where appropriate
% badtrials = indices of bad trials
% EEGbu = backup of the original, un-altered EEGlab structure

%% create backup of eeg data
EEGbu = EEG;

%% get properties of data matrix
[ch,tp,tr] = size(EEG.data);

%% find bad epochs
if ~isempty(epochlims) && tr == 1
    disp('Running bad trial detection routines')
    disp('EEG output will be chans x frames x trials')
    [EEGepoch, indices] = pop_epoch(EEG, eventcodes, [epochlims(1) epochlims(2)]);
    EEGepoch = eeg_checkset(EEGepoch);
else
    EEGepoch = EEG;
end
[badtrials1] = identifyarttrial_AV101(EEGepoch,resting);
badtrials1 = badtrials1';
[badtrials2] = bad_trials_av101(EEGepoch,[]);
if ~isempty(epochlims)
[~, badtrials3, ~, ~] = eegthresh( EEGepoch.data, EEGepoch.pnts, ...
    1:length(EEGepoch.chanlocs), -500, 500, [epochlims(1) epochlims(2)], epochlims(1),epochlims(2));
else
    [~, badtrials3, ~, ~] = eegthresh( EEGepoch.data, EEGepoch.pnts, ...
    1:length(EEGepoch.chanlocs), -500, 500, [EEGepoch.xmin,EEGepoch.xmax], EEGepoch.xmin,EEGepoch.xmax);
end
[OUTEEG, ~, ~, ~] = ...
    pop_jointprob( EEGepoch, 1, 1:length(EEGepoch.chanlocs), ...
    6, 2, 0, 0, 0);
badtrials4 = find(OUTEEG.reject.rejjp == 1);
badtrials = sort(unique([badtrials1,badtrials2,badtrials3,badtrials4]));
clear OUTEEG
EEG = EEGepoch;
EEG = pop_select(EEG, 'notrial', badtrials);

%% Re-reference the data to average
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:,:) = zeros(1, EEG.pnts,EEG.trials);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, []);
EEG = pop_select(EEG,'nochannel',{'initialReference'});

end