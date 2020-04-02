function [EEGfinal,EEGprerej,badcomps,summary] = av101_handle_bad_data(EEG,epochlims,eventcodes,fmin,fmax,icmethod,resting,bl,ep)
% av101_handle_bad_data
% function for managing the pre-processing of AV101 data.
% inputs:
% EEG = eeglab structure
% epochlims = (in seconds) upper and lower limits for epoching. leave blank
% if resting state.
% eventcodes = cell containing the string codes for each condition.
% fmin/fmax = min and max frequencies for assessing bad spectral
% contributions
% icamthethod = method for decomposing data with ICA
%   1 ) standard ica method, epoching, pruning, pca initialised
%   2 ) standardica method, no epoching or pruning, pca initialised
%   3 ) Wavelet ICA using RUNICA
%   4 ) Wavelet ICA using Radical ICA
% resting = 0/1, non-resting/resting state
%
%outputs
% EEGfinal  = EEGLab structure containing cleaned data
% EEGprerej = backup of EEGlab structure before component rejection
% badcomps = components marked for rejection
% summary = summary of the cleaning process


%% determine if the data will need to be pruned, and what ica method to use
if icmethod == 1 % data is pruned for the standard ica approach to help improve snr
    disp('preparing data for standard method')
    if ep == 1
        disp('data will be epoched and pruned before ica');
        epochflag = 1;
    else
        epochflag = 0;
    end
elseif icmethod == 2 % data is not pruned for wICA as this is antithetical to it's purpose
    disp('preparing data for wavelet ICA - RUNICA')
    epochflag = 0;
    icflag = 'runica';
    
elseif icmethod == 3
    disp('preparing data for wavelet ICA - RADICAL')
    epochflag = 0;
    icflag = 'radical';
end



%% find bad channels (& prune data if icm1)
full_selected_channels = EEG.chanlocs;
if epochflag == 1
    [EEG,badchans,badtrials1,EEGbu,asr_used] = av101_findbaddata_preica(EEG,fmin,fmax,epochlims,eventcodes,1);
else
    [EEG,badchans,badtrials1,EEGbu,asr_used] = av101_findbaddata_preica(EEG,fmin,fmax,[],[],1);
end
if asr_used == 0
    icchans = 1:length(EEG.chanlocs);
    icchans(badchans) = [];
    EEG = pop_select(EEG,'channel', icchans);
end
% if asr_used == 1
%     icmethod = 1;
% end

%% run ica
if asr_used == 0
    if icmethod == 1
        % run ica with pca and fastica initialisation for runica
        XX = EEG.data;
        [W,L,w,s,time]=decompose2(XX,[],[],[],[],[],[]);
        icaact=data2act(XX,W);
        %     IC = reshape(icaact,size(icaact,1),[]);
        %     EEG.icaact = reshape(IC,size(IC,1),size(EEG.data,2),size(EEG.data,3));
        EEG.icaact = icaact;
        EEG.icaweights = w;
        EEG.icawinv = W;
        EEG.icasphere = s;
        EEG = eeg_checkset( EEG );
    elseif icmethod == 2 || icmethod == 3
        [wIC, A, W, IC] = wICA(EEG,icflag, 1, 0, [], 5);
        artifacts = A*wIC;
        EEG2D=reshape(EEG.data, size(EEG.data,1), []);
        wavcleanEEG=EEG2D-artifacts;
        EEG.data = wavcleanEEG;
        EEG = pop_runica(EEG, 'extended',1,'interupt','off');
        EEG = eeg_checkset( EEG );
        
    end
    
    %% flag bad components
    [bad_comps] = av101_badcomponents(EEG);
    [~,temp,~]=processMARA ( EEG,EEG,EEG, [0, 0, 0, 0 , 0] );
    % test version with plotting
    % [~,temp,~]=processMARA ( EEG,EEG,EEG, [0, 0, 1,1 , 0] );
    mbc = zeros(size(temp.reject.gcompreject));
    mbc(temp.reject.MARAinfo.posterior_artefactprob > 0.5) = 1;
    
    %     if epochflag == 0
    %         disp('running additional blink component identification')
    %         % find blinks specifically --- only possible in continuous data
    %         EEG.icaquant = icablinkmetrics(EEG);
    %         [~,index] = sortrows([EEG.icaquant.metrics.convolution].');
    %         EEG.icaquant.metrics = EEG.icaquant.metrics(index(end:-1:1)); clear index
    %         if sum(EEG.icaquant.identifiedcomponents)>0
    %             ibc = EEG.icaquant.identifiedcomponents;
    %         else
    %             ibc = [];
    %         end
    %     else
    ibc = [];
    %     end
    
    
    % %%%% test notes %%%%%
    % % plot components
    % figure;viscomp(w,mean(EEG.data,3),EEG.chanlocs,[],[],'im',8,[],[],[],[]);
    % figure;viscomp(w,EEG.data(:,:),EEG.chanlocs,[],[],'im',8,[],[],[],[]);
    badcomps = sort(unique([bad_comps',find(mbc),ibc]));
    EEGprerej = EEG;
    EEG = pop_subcomp( EEG, badcomps, 0);
else
    EEGprerej = EEG;
    badcomps = [];
    %     ibc = [];
    
end
%% if resting data do segmenting and interpolate within segments as well as channels
if resting == 1
    EEG=eeg_regepochs(EEG,'recurrence',4,'limits',[0 4], 'rmbase', [NaN]);
    epochlims = [];
    EEG.times2 = EEG.times;
    EEG.times = linspace(0,4,size(EEG.data,2));
    EEG = eeg_checkset( EEG );
end

%% if data are not already epoched do epoching
[~,~,tr] = size(EEG.data);
% if ~isempty(epochlims) && tr == 1
%     [EEG,~,badtrials2,~] = av101_findbaddata(EEG,fmin,fmax,epochlims,eventcodes,0);
% else
%     [EEG,~,badtrials2,~] = av101_findbaddata(EEG,fmin,fmax,[],[],0);
% end
[EEG,badtrials2,EEGbu] = av101_findbaddata_postica(EEG,epochlims,eventcodes,resting);
%% interpolate bad channels
EEG = pop_interp(EEG, full_selected_channels, 'spherical');
EEG = eeg_checkset(EEG );
EEGfinal = EEG;
EEGfinal.nbchan = EEGfinal.nbchan+1;
EEGfinal.data(end+1,:,:) = zeros(1, EEGfinal.pnts,size(EEGfinal.data,3));
EEGfinal.chanlocs(1,EEGfinal.nbchan).labels = 'initialReference';
EEGfinal = pop_reref(EEGfinal, []);
EEGfinal = pop_select(EEGfinal,'nochannel',{'initialReference'});
if resting == 0
    EEGfinal = pop_rmbase(EEGfinal, [bl(1), 0]);
end

summary.badcomps = badcomps;
summary.badtrials = length([badtrials1,badtrials2]);
summary.badchans = badchans;
summary.chansused = full_selected_channels;
summary.bl = bl;
summary.evcodes = eventcodes;
summary.eplims = epochlims;
icmethods = {'standard','WICA_runica','WICA_radical'};
summary.icmethod = icmethods{icmethod};
summary.asr_used = asr_used;


end