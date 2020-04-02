function [bad_comps] = av101_badcomponents(EEG)
% av101_badcomponents
% statistical and calculus derived estimation of bad independent components
% EEG = eeglab structure
%
% adapted from FASTER toolbox 
% Nolan, H., Whelan, R.,&  Reilly, R.B. (2010). FASTER: 
% Fully Automated Statistical Thresholding for EEG artifact Rejection./Journal of Neuroscience Methods, 192/, 152-162.


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
[lengths] = min_z(list_properties);
bad_comps=find(lengths);
end



function [lengths] = min_z(list_properties,rejection_options)
if (~exist('rejection_options','var'))
    rejection_options.measure=ones(1,size(list_properties,2));
    rejection_options.z=3*ones(1,size(list_properties,2));
end

rejection_options.measure=logical(rejection_options.measure);
zs=list_properties-repmat(mean(list_properties,1),size(list_properties,1),1);
zs=zs./repmat(std(zs,[],1),size(list_properties,1),1);
zs(isnan(zs))=0;
all_l = abs(zs) > repmat(rejection_options.z,size(list_properties,1),1);
lengths = any(all_l(:,rejection_options.measure),2);
end