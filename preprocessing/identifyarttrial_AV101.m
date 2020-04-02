function [arttrial] = identifyarttrial_AV101(EEG, resting)
% identifyarttrial_AV101
% statistical estimation bad data using baseline and active periods of
% epoched data.
% EEG = eeglab structure
% resting = 0/1, not-resting/resting.
% outputs
% arttrial = indices of bad trials
%
% adapted from the artist pipeline (yu et al, 2018).


% FIND BAD TRIALS
if resting == 1
    time = 1:size(EEG.data,2);
else
    time = (EEG.times >= 100)|(EEG.times <= 0);% Only look at baseline and beyond 300 ms
end
trialpow = squeeze(mean(EEG.data(:,time,:).^2,2));
meanpow = repmat(median(trialpow, 2), 1, size(EEG.data, 3));
stdpow = repmat(iqr(trialpow, 2), 1, size(EEG.data, 3));
datTmp = squeeze(nanmean(abs(EEG.data),2));
ChanWithLargeSTD = find(nanstd(datTmp')>30);
trialthr = 5;
artchannum = repmat(sum((trialpow-meanpow)./stdpow > trialthr, 1), size(EEG.data, 1), 1);
chanthr = round(0.2*EEG.nbchan);
[chan, trial] = find(((trialpow-meanpow)./stdpow > trialthr) & (artchannum < chanthr)); % only interplate trials with less than 10 noisy electrodes. Otherwise simply discard
arttrial = unique(trial);
for kk = 1:length(arttrial)
    artind = (trial == arttrial(kk));
    EEGart = pop_select(EEG, 'trial', arttrial(kk));
    EEGart = eeg_interp(EEGart,chan(artind));
    EEG.data(:, :, arttrial(kk)) = EEGart.data;
end
[~, trial] = find(((trialpow-meanpow)./stdpow > trialthr) & (artchannum >= chanthr));
arttrial = unique(trial);

