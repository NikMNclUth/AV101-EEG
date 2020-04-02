function [EEGfinal] = AV101_RS_clean_loop(filenames,codes,visits,tasks,tpp,chan_cut)
% menninger_clean_loop
% Processing Pipeline for Menninger Project Data.
% Dependencies: EEGLab V13 or higher; MARA; FASTER; In-house tools ( see
% below ).
% setup_analysispath(1);
% setup_analysispath(5);
% setup_analysispath(6);
% setup_analysispath(7);
%
% Inputs:
%
% Outputs:
% EEGfinal = Pre-processed EEG structure.
% Processing Log (subjectlog) = table containing:
% 'Subject_Code','First_EEG_File_Code','Date_Processed','Total_Bad_Chans','MARA_Run','Total_Bad_Comps',...
%     'Total_Bad_Trials','capsize'
%
%
% The pipeline consists of the following steps:
% 1) Filter - 1 to 49 hz, cleaningline for removal of line noise.
% 2) ASR - remove bad channels and bad segments (repairs within reason)
% using threshold of 30 STD.
% 3) ICA - extended infomax, pca reduced to match rank.
% 4) MARA/FASTER - routines to identify bad components.
% 5) Reject bad data - segment resting data. regress out remaining ocular
% activity (blinks and advanced properties). Find bad trials using trial
% statistics.
% 6) Interpolate bad channels.
% 7) Re-reference to average.

%% load appropriate data
[EEGloaded,EEG,full_selected_channels] = av101_loaddata(filenames,codes,...
    visits,tasks,tpp,chan_cut);

%% Filtering
[EEG] = av101_filter_data(EEG,1,125,250);
EEG = eeg_checkset(EEG);
%% Create Backup
% EEGbu = EEG;
%%% add parpool options
%% ASR
EEGtemp = EEG;
EEGtemp = clean_rawdata(EEGtemp, 5, [0.25 0.75], 0.8, 4, 30, 0.25);
numbc = size(EEG.data,1)-size(EEGtemp.data,1);
EEGtemp = eeg_checkset( EEGtemp );
% vis_artifacts(EEGtemp,EEG);
%% ICA
EEGtemp = pop_runica(EEGtemp, 'extended',1,'interupt','off','reset_randomseed','off','pca',rank(EEGtemp.data(:,:)));
% Bad Components - FASTER
[bad_comps] = menninger_badcomponents(EEGtemp,3);
% Bad Components - MARA

% if floor(max(EEGtemp.times)/1000)>=15 % need minimum 15 sec for MARA
if ~isempty(1:15:length(EEGtemp.icaact(1,:))/EEG.srate-15) % check whether data will meet MARA criteria for time freq feature analysis
    %         how to remove dialogue box?
    [~,temp,~]=processMARA ( EEGtemp,EEGtemp,EEGtemp, [0, 0, 0, 0 , 0] );
    mbc = zeros(size(temp.reject.gcompreject));
    mbc(temp.reject.MARAinfo.posterior_artefactprob > 0.5) = 1;
    badcomps = sort(unique([bad_comps',find(mbc)]));
    mararun = 'yes';
else
    disp('issue running mara --- logging');
    badcomps = bad_comps;
    mararun = 'no';
end

% Remove Bad Components
EEGtemp = pop_subcomp( EEGtemp, badcomps, 0);
EEGtemp = eeg_checkset( EEGtemp );
%% Reject Bad Data
% Cut Data Into Pieces (not my last resort)
% EEGtemp=eeg_regepochs(EEGtemp,'recurrence',1,'limits',[0 2], 'rmbase', [NaN]); % some data length is lost as a result of this
EEGtemp = eeg_regepochs(EEGtemp,'recurrence',1); % use epochs of length 1 second for identifying artifacts
EEGtemp = eeg_checkset( EEGtemp );
bch = logical(zeros(size(EEGtemp.data,1),1));
% Regress Out Residual Ocular Artifacts

try
    EYE=findocularsignals(EEGtemp.data(:,:),EEGtemp.srate,EEGtemp.chanlocs,find(~bch),[],[]);
    EYE=resh(EYE,size(EEGtemp.data,2));
    X=EEGtemp.data(~bch,:,:);
    X=regressout(X,EYE(1,:,:));
catch
    disp('data not appropriate for advanced occular cleaning');
    X=EEGtemp.data(~bch,:,:);
end
try
    [X,bTR,~,~]=finaldatacleanup(X,EEGtemp.srate,EEGtemp.chanlocs,bch,3,3);
    badtrials0 = find(bTR==1);
    EEGtemp = pop_select(EEGtemp, 'notrial', badtrials0);
    EEGtemp.data = X;
catch
    disp('data quality not sufficient for rigorous cleanup');
end
EEGtemp = eeg_checkset(EEGtemp);

% Bad Trial Identification
[badtrials1] = identifyarttrial_AV101(EEGtemp, 1);
badtrials1 = badtrials1';
[badtrials2] = bad_trials_av101(EEGtemp,[]);
[~, badtrials3, ~, ~] = eegthresh( EEGtemp.data, EEGtemp.pnts, ...
    1:length(EEGtemp.chanlocs), -500, 500, [EEGtemp.xmin,EEGtemp.xmax], EEGtemp.xmin,EEGtemp.xmax);
[OUTEEG, ~, ~, ~] = ...
    pop_jointprob( EEGtemp, 1, 1:length(EEGtemp.chanlocs), ...
    6, 2, 0, 0, 0);
badtrials4 = find(OUTEEG.reject.rejjp == 1);
badtrials = sort(unique([badtrials1,badtrials2,badtrials3,badtrials4]));
EEGtemp = pop_select(EEGtemp, 'notrial', badtrials);



%% Finalise Data
% Interpolate Removed Channels
EEGtemp = pop_interp(EEGtemp, EEG.chanlocs, 'spherical');
EEGtemp = eeg_checkset(EEGtemp );
EEGfinal = EEGtemp;
% Average Re-Reference
EEGfinal.nbchan = EEGfinal.nbchan+1;
EEGfinal.data(end+1,:,:) = zeros(1, EEGfinal.pnts,size(EEGfinal.data,3));
EEGfinal.chanlocs(1,EEGfinal.nbchan).labels = 'initialReference';
EEGfinal = pop_reref(EEGfinal, []);
EEGfinal = pop_select(EEGfinal,'nochannel',{'initialReference'});
% figure;for ii = 1:10;subplot(5,2,ii);plot(EEGfinal.data(ii,1:2000));hold on; plot(EEG.data(ii,1:2000));end

end


function [bad_comps] = menninger_badcomponents(EEG,threshold)
% av101_badcomponents
% statistical and calculus derived estimation of bad independent components
% EEG = eeglab structure
%
% adapted from FASTER toolbox


for u = 1:size(EEG.icaact,1)
    [spectra(u,:) freqs] = pwelch(EEG.icaact(u,:),[],[],(EEG.srate),EEG.srate);
end
list_properties = zeros(size(EEG.icaact,1),5);
for u=1:size(EEG.icaact,1)
    measure = 1;
    % 1 Median gradient value, for high frequency stuff
    list_properties(u,measure) = median(diff(EEG.icaact(u,:)));
    measure = measure + 1;
    % 2 Mean slope around the LPF band (spectral)
    list_properties(u,measure) = mean(diff(10*log10(spectra(u,find(freqs>=5,1):find(freqs<=40,1,'last')))));
    measure = measure + 1;
    % 3 Kurtosis of spatial map (if v peaky, i.e. one or two points high
    % and everywhere else low, then it's probably noise on a single
    % channel)
    list_properties(u,measure) = kurt(EEG.icawinv(:,u));
    measure = measure + 1;
    % Hurst Exponent
    list_properties(u,measure) = hurst_exponent(EEG.icaact(u,:));
    measure = measure + 1;
    % Timeseries Kurtosis
    list_properties(u,measure) = kurt(EEG.icaact(u,:));
end
for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
    list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end
[lengths] = min_z(list_properties,threshold);
bad_comps=find(lengths);
end



function [lengths] = min_z(list_properties,threshold)

rejection_options.measure=ones(1,size(list_properties,2));
rejection_options.z=threshold*ones(1,size(list_properties,2));


rejection_options.measure=logical(rejection_options.measure);
zs=list_properties-repmat(mean(list_properties,1),size(list_properties,1),1);
zs=zs./repmat(std(zs,[],1),size(list_properties,1),1);
zs(isnan(zs))=0;
all_l = abs(zs) > repmat(rejection_options.z,size(list_properties,1),1);
lengths = any(all_l(:,rejection_options.measure),2);
end