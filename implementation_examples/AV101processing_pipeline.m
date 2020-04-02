%% PROCESSING PIPELINE FOR AV101 DATA
% Dr. Nicholas Murphy, PhD, October 2019, Baylor College of Medicine
% Version 2.1
% last edited 10/09/2019
% Notes on usage:
% - pipeline is run through separately for each visit.
% - timepoints are identified using dynamic filename creations
% - make sure to set epoch and baseline limits for each event related data
%   set.
% --- MMN epochs: -0.1 to 0.4; baseline = -100 to 0
% --- ASSR epochs: -1.5 to 1.5; baseline = -100 to 0
% --- resting data is split into 4 second segments
% Preprocessing specifics:
% - final input for av101_handle_bad_data should be 0 for erp, and 1 for
%   resting.
% - ICA should be run using Wavelet ICA with continous data for best
%   performance. icmethod = 2
% - if >= 50% of the channelsare bad ASR is run and the data are flagged.
%   see summary.asr_used (1 = yes, 0 = no)
%
%
%% Set up loading paths
delete(gcp('nocreate'))
clear; clc;
%%%% add all necessary paths here
eeglab;
clear
close all
clc
[filenames] = av101_filenames('path with files to use');% generate nested structure with complete file paths 
% filenames is only really necessary with a high volume of data where naming conventions might share similarities but contain errors.
% this can be generated manually if necessary.

rmpath('/data/rcho/TOOLS/eeglab14_1_1b/plugins/AAR131130/');
% loading variables
visits = {'V2','V3','V4'};
codes = {'04','05','06','07','08','09','10','13','14','15','17','18'};
tasks = {'resting','assr40','assr30','assr20','mmn'};
timepoints = {'b','60','120','180','240','300'};
% pre-processing variables
chan_cut = 0; % use full montage
fmin = 1; % minimum frq of 1 Hz for artifact rejection
fmax = 100; % maximum frq of 100 Hz for artifact rejection
fnames = fieldnames(filenames); % subject codes
%% MMN
mmn_evcodes = [1,2,99];
mmn_altcodes = [5,6,103];
% Visit 2
savepath = 'path to save to';
for ii = 1:length(fnames)
    % not all data have a recording at each visit. it is faster to
    % initially identify if the data can be loaded or not, if it does not
    % exist the error will be caught and we move on to the next subject.
    % get existing time points
    % not all subjects have a recording at each stage, this section
    % identifies the existing ones and uses them to feedforwards into the
    % loading function.
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V2')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V2.mmn'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V2.mmn']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{1},tasks{5},tpp{timer},chan_cut);
            % get event codes
            alleventtypes = {EEG.event.type}';
            evs=double(cell2mat(alleventtypes));
            evs2 = unique(evs);
            % edit alternative codes
            % in this section we search for event codes that are equal to the
            % alternative event codes and then replace them with their
            % equivalent standard index.
            if sum(ismember(evs2,mmn_altcodes))>0
                f1 = find(evs==mmn_altcodes(1));
                f2 = find(evs==mmn_altcodes(2));
                f3 = find(evs==mmn_altcodes(3));
                if ~isempty(f1)
                    for ac1 = 1:length(f1);EEG.event(f1(ac1)).type = mmn_evcodes(1);end
                end
                if ~isempty(f2)
                    for ac1 = 1:length(f2);EEG.event(f2(ac1)).type = mmn_evcodes(2);end
                end
                if ~isempty(f3)
                    for ac1 = 1:length(f3);EEG.event(f3(ac1)).type = mmn_evcodes(3);end
                end
            end
            clear f1 f2 f3 ac1 evs evs2 alleventtypes iter f fff
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,[]);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-0.1,0.4],{1,2,99},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{5},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
    
end
% Visit 3
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V3')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V3.mmn'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V3.mmn']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{2},tasks{5},tpp{timer},chan_cut);
            % get event codes
            alleventtypes = {EEG.event.type}';
            evs=double(cell2mat(alleventtypes));
            evs2 = unique(evs);
            if sum(ismember(evs2,mmn_altcodes))>0
                f1 = find(evs==mmn_altcodes(1));
                f2 = find(evs==mmn_altcodes(2));
                f3 = find(evs==mmn_altcodes(3));
                if ~isempty(f1)
                    for ac1 = 1:length(f1);EEG.event(f1(ac1)).type = mmn_evcodes(1);end
                end
                if ~isempty(f2)
                    for ac1 = 1:length(f2);EEG.event(f2(ac1)).type = mmn_evcodes(2);end
                end
                if ~isempty(f3)
                    for ac1 = 1:length(f3);EEG.event(f3(ac1)).type = mmn_evcodes(3);end
                end
            end
            clear f1 f2 f3 ac1 evs evs2 alleventtypes iter f fff
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-0.1,0.4],{1,2,99},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{5},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 4
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V4')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V4.mmn'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V4.mmn']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{3},tasks{5},timepoints{timer},chan_cut);
            % get event codes
            alleventtypes = {EEG.event.type}';
            evs=double(cell2mat(alleventtypes));
            evs2 = unique(evs);
            if sum(ismember(evs2,mmn_altcodes))>0
                f1 = find(evs==mmn_altcodes(1));
                f2 = find(evs==mmn_altcodes(2));
                f3 = find(evs==mmn_altcodes(3));
                if ~isempty(f1)
                    for ac1 = 1:length(f1);EEG.event(f1(ac1)).type = mmn_evcodes(1);end
                end
                if ~isempty(f2)
                    for ac1 = 1:length(f2);EEG.event(f2(ac1)).type = mmn_evcodes(2);end
                end
                if ~isempty(f3)
                    for ac1 = 1:length(f3);EEG.event(f3(ac1)).type = mmn_evcodes(3);end
                end
            end
            clear f1 f2 f3 ac1 evs evs2 alleventtypes iter f fff
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-0.1,0.4],{1,2,99},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{5},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end

%% ASSR 40
% Visit 2
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V2')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V2.assr40'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V2.assr40']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{1},tasks{2},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{2},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 3
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V3')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V3.assr40'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V3.assr40']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{2},tasks{2},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{2},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 4
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V4')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V4.assr40'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V4.assr40']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{3},tasks{2},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{2},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end

%% ASSR 30
% Visit 2
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V2')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V2.assr30'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V2.assr30']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{1},tasks{3},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{3},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 3
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V3')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V3.assr30'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V3.assr30']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{2},tasks{3},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{3},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 4
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V4')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V4.assr30'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V4.assr30']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{3},tasks{3},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{3},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
        
    end
end
%% ASSR 20
% Visit 2
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V2')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V2.assr20'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V2.assr20']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{1},tasks{4},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{4},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 3
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V3')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V3.assr20'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V3.assr20']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{2},tasks{4},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{4},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end
% Visit 4
savepath = '';
for ii = 1:length(fnames)
    
    eval(['fn=fieldnames(filenames.',fnames{ii},');']);
    fn = find(~cellfun(@isempty,strfind(fn,'V4')));
    if ~isempty(fn)
        fnoms = lower(eval(['filenames.',fnames{ii},'.V4.assr20'])); % create lower case file names
        fff = [];
        for iter = 1:length(timepoints)
            f=find(~cellfun(@isempty,strfind(fnoms,timepoints{iter})));
            if ~isempty(f)
                fff(iter) = f;
            end
        end
        tpp = timepoints(sort(fff)); % timepoints available
        for timer = 1:length(eval(['filenames.',fnames{ii},'.V4.assr20']))
            % load appropriate data
            [EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes{ii},...
                visits{3},tasks{4},tpp{timer},chan_cut);
            alleventtypes = {EEG.event.type}';
            assrtype = mode(double(cell2mat(alleventtypes)));
            % apply filters
            [EEG] = av101_filter_data(EEG,1,125,250);
            % wavelet method (runica), without epoching
            [EEGfinal,~,~,summary] = av101_handle_bad_data(EEG,[-1.5,1.5],{assrtype},fmin,fmax,2,[],[-100,0],0);
            % save
            savename = [savepath,fnames{ii},'_',tpp{timer},'_',tasks{4},'.mat'];
            save(savename,'EEGfinal','summary','-v7.3');
            clear EEGloaded EEG full_selected_channels EEGfinal summary
        end
    end
end

