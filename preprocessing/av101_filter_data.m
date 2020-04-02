function [EEG] = av101_filter_data(EEG,lcut,hcut,newsr)
% av101_filter_data
% function for applying filtering and downsampling to the AV101 data
%
% inputs:
% EEG = eeglab structure
% lcut = lower band of the desired frequency range (high pass filter)
% hcut = upper band of the desiered frequency range (low pass filter)
% newsr = new sampling rate for downsampling.
%
% outputs:
% EEG = filtered eeglab structure. No modifications made to the rest of the
% data besides filtering.




%% step 1 - filter and downsample

% lcut = 1;hcut = 100;
if EEG.srate == hcut*2
    hcut = hcut-1;
end
EEG = pop_eegfiltnew(EEG, 0,hcut);
if ~isempty(newsr)
    EEG = pop_resample(EEG, newsr);
end
EEG = pop_eegfiltnew(EEG, lcut, 0);


%% step 2 - remove line noise


EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan],'computepower',0,'linefreqs',...
    [60 120 180 240] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
    'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
EEG.setname='rawEEG_f_cs_ln';
EEG = eeg_checkset(EEG);

end