function [EEG,badchans,badtrials,EEGbu,asr_used] = av101_findbaddata_preica(EEG,fmin,fmax,epochlims,eventcodes,chanrej)
% av101_findbaddata_preica
% function for identifying bad channels and trials to be removed prior to
% ica.
% inputs:
% EEG = eeglab structure
% fmin/fmax = min and max frequencies for assessing bad spectral
% contributions
% epochlims = (in seconds) upper and lower limits for epoching.
% eventcodes = cell containing the string codes for each condition.
%%%% leave limits and codes blank if not planning to prune the data of
%%%% trials before ica (eg wavelet ica).
% chanrej = 0/1, keep/reject bad channels (will be interpolated later)
%
% outputs:
% EEG = eeglab structure with indices from processing added and data
% rejected where appropriate
% badchans = indices of bad channels
% badtrials = indices of bad trials
% EEGbu = backup of the original, un-altered EEGlab structure


%% create a back up of the OG eeglab structure

EEGbu = EEG;
[nchans,tp,tr] = size(EEG.data);
%% run faster and happe bad channel routines
if chanrej == 1
    disp('Running bad channel detection routines')
    [badinds1] = bad_channels(EEG);
    [~, badinds2,~] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],'threshold',[-3 3],...
        'norm','on','measure','spec','freqrange',[fmin fmax]);
    tempchans = 1:EEG.nbchan;tempchans(badinds2) = [];
    if ~isempty(badinds2)
        [~, badinds3,~] = pop_rejchan(EEG, 'elec',tempchans,'threshold',[-3 3],...
            'norm','on','measure','spec','freqrange',[fmin fmax]);
        badinds2 = [badinds2,tempchans(badinds3)];
    end
    badchans = sort(unique([badinds1,badinds2]));
else
    badchans = [];
end

%% find bad epochs
if numel(badchans)>=floor(nchans/2)
    % if we have an extreme number of bad channels (50% or more) then clean
    % using ASR - this will likely mean that the event related data might
    % not be useable.
    EEG = clean_artifacts(EEG);
    %     vis_artifacts(EEG,EEGbu);
    if size(EEG.data,1)<size(EEGbu.data,1)
        badchans = find(EEG.etc.clean_channel_mask == 0);
        asr_used = 1;
        badtrials = [];
    else 
        badchans = [];
        asr_used = 0;
        badtrials = [];
    end
else
    if ~isempty(epochlims)
        disp('Running bad trial detection routines')
        disp('EEG output will be chans x frames x trials')
        [EEGepoch, indices] = pop_epoch(EEG, eventcodes, [epochlims(1) epochlims(2)]);
        EEGepoch = eeg_checkset(EEGepoch);
        [badtrials1] = identifyarttrial(EEGepoch, []);
        badtrials1 = badtrials1';
        [badtrials2] = bad_trials_av101(EEGepoch,badchans);
        
        [~, badtrials3, ~, ~] = eegthresh( EEGepoch.data, EEGepoch.pnts, ...
            1:length(EEGepoch.chanlocs), -500, 500, [epochlims(1) epochlims(2)], epochlims(1),epochlims(2));
        
        [OUTEEG, ~, ~, ~] = ...
            pop_jointprob( EEGepoch, 1, 1:length(EEGepoch.chanlocs), ...
            6, 2, 0, 0, 0);
        badtrials4 = find(OUTEEG.reject.rejjp == 1);
        badtrials = sort(unique([badtrials1,badtrials2,badtrials3,badtrials4]));
        clear OUTEEG
        EEG = EEGepoch;
        EEG = pop_select(EEG, 'notrial', badtrials);
    else
        disp('EEG output will be chans x frames')
        badtrials = [];
    end
    asr_used = 0;
end


%% Re-reference the data to average
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:,:) = zeros(1, EEG.pnts,EEG.trials);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, [],'exclude',badchans);
EEG = pop_select(EEG,'nochannel',{'initialReference'});

end