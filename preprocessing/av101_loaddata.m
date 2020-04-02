function [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,code,visit,task,timepoint,chan_cut)
% av101_loaddata
% function for loading AV101 data. Requires that filenames for the study
% have been generated and organised.
% requires readlocs to be in the path
% inputs:
% - filenames = all task files -- see [filenames] = av101_filenames('/data/smathew/AV-101/data/');
% - code = numeric code for subject (as string), eg. '04'
% - visit = 'V2','V3', or 'V4'
% - task = 'mmn','resting','assr40','assr30','assr20'
% - timepoint = 'b','60','120','180','240','300'
% - chan_cut = 1/0 prune or do not prune to easycap 10-20 layout
% outputs:
% EEGloaded = original raw eeg (backup)
% EEG = active version of the eeg, including any edits made to layout
% full_selected_channels = list of active channels 


%% get subject name
fnames = fieldnames(filenames);
subjname = fnames(find(~cellfun(@isempty,strfind(lower(fieldnames(filenames)),lower(code)))));
loadtext = ['subjdata = filenames.',subjname{1},'.',visit,'.',task,';'];
%% get timepoint name
eval(loadtext)
subjdata = subjdata(find(~cellfun(@isempty,strfind(lower(subjdata),lower(timepoint)))));
%% get correct cap information
eval(['chanlocs = filenames.',subjname{1},'.',visit,'.chanlocs;']);
if strmatch('easycap',chanlocs)
    eloc = readlocs('/data/smathew/AV-101/data/easycap channel configuration.ced');
else
    load('/data/smathew/AV-101/data/black_cap_adj2.mat');
    eloc = elocB2;
    % re-label the mastoid electrodes to match brainvision 10-20 labels
    for ii = 1:length(eloc)
        if strmatch('M1',eloc(ii).labels)
            eloc(ii).labels = 'TP9';
        end
        if strmatch('M2',eloc(ii).labels)
            eloc(ii).labels = 'TP10';
        end
        
    end
end
chan_IDs={'FP1', 'FP2', 'F7','F3','FZ','F4','F8','FC5','FC1','FC2',...
    'FC6','T7','C3','CZ','C4','T8','TP9','CP5','CP1','CP2','CP6','TP10',...
    'P7','P3','PZ','P4','P8','O1','O2'};
%% load data
EEGloaded = loadcurry(subjdata{1},'CurryLocations','True');
EEGloaded = eeg_checkset(EEGloaded);
EEGloaded.chanlocs = eloc;
events=EEGloaded.event;
complete_event_info=EEGloaded.urevent;
srate=double(EEGloaded.srate);
EEGloaded.setname='rawEEG';
EEG = eeg_checkset(EEGloaded);
%% prune channels
if chan_cut == 1
    EEG = pop_select( EEG,'channel', chan_IDs);
    EEG.setname='rawEEG_f_cs';
    EEG = eeg_checkset( EEG );
end
full_selected_channels = EEG.chanlocs;

end