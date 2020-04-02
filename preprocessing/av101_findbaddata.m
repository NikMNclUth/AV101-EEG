function [EEG,badchans,badtrials,EEGbu] = av101_findbaddata(EEG,fmin,fmax,epochlims,eventcodes,chanrej)

EEGbu = EEG;

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




%% Re-reference the data to average
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:,:) = zeros(1, EEG.pnts,EEG.trials);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, [],'exclude',badchans);
EEG = pop_select(EEG,'nochannel',{'initialReference'});

end