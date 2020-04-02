
function [orig_spectra,loginterp_pow,rp,ap,lp,rotatedpow] = compute_spectra_av101(EEGfinal,roi)

  % uses functions adapted from Colombo et al, (2019)
% use this function to estimate the spectra and it's properties.
% psd is put into a function for estimating the 1/f properties of the data,
% which can be converted into plotting variables also. Use the output from
% this for 1/f estimate and for rotated power measurement (psd without
% 1/f).
% Use the raw spectra (10*log10) for standard power estimates and
% comparison.



%% detrending & demeaning
dataf = zeros(size(EEGfinal.data(roi,:,:)));
for ii = 1:size(EEGfinal.data(roi,:),3)
    dataf(:,:,ii) = EEGfinal.data(roi,:,ii)-nanmean(EEGfinal.data(roi,:,ii),2);
    temp = detrend(dataf(:,:,ii)');
    dataf(:,:,ii) = temp';
    clear temp
end
%% compute spectra
for u = 1:size(dataf,1)
    if u == 1
        [s,freqs] = pwelch(dataf(u,:),hanning(EEGfinal.srate),EEGfinal.srate/2,(EEGfinal.srate),EEGfinal.srate);
        spectra = zeros(size(dataf,1),length(s));
        spectra(u,:) = s;
    else
        [spectra(u,:),~] = pwelch(dataf(u,:),hanning(EEGfinal.srate),EEGfinal.srate/2,(EEGfinal.srate),EEGfinal.srate);
    end
end
%% Estimate slope
spectra = nanmean(spectra,1);
% doPlot= 1;
% thisCol= [0 0 1];
% spectra = spectra./sum(spectra);
frBand=[3 45]; % frequencies to use
frBins= dsearchn(freqs, frBand(1) ):dsearchn(freqs,frBand(2));
XX= freqs(frBins); % cutdown frequencies% [global_power_vars_log_orig] = measure_power(global_orig_spectra.log_spectra,global_orig_spectra.freqs);

YY= spectra(frBins); % spectrum cutdown by bins
robRegMeth= 'ols';
[intSlo, ~, Pows, ~,  ~, ~] = fitPowerLaw3steps(XX,YY, robRegMeth,  0, []);
% 1/f removed from power
loginterp_pow.slope =intSlo(2);
loginterp_pow.intercept = intSlo(1);
loginterp_pow.observed_power = Pows.obs;
loginterp_pow.oof_fit = Pows.pred;
loginterp_pow.rotated_spectrum = Pows.res; % spectrum with 1/f removed (use for power measurements).
loginterp_pow.used_freqs = Pows.frex; % log spaced

%% estimate power from data including 1/f
orig_spectra.absolute_spectra = spectra(:,find(freqs==3):find(freqs==30));
orig_spectra.freqs = freqs(find(freqs==3):find(freqs==30));
orig_spectra.log_spectra = 10*log10(orig_spectra.absolute_spectra);
orig_spectra.relative_spectra = orig_spectra.absolute_spectra./sum(orig_spectra.absolute_spectra);

%% estimate band power

rp = measure_power(orig_spectra.relative_spectra,orig_spectra.freqs);
ap = measure_power(orig_spectra.absolute_spectra,orig_spectra.freqs);
lp = measure_power(orig_spectra.log_spectra,orig_spectra.freqs);
rotatedpow = measure_power(loginterp_pow.rotated_spectrum,loginterp_pow.used_freqs);

end

function [power_vars] = measure_power(spectral_data,frequencies)
power_vars.theta_power = nanmean(spectral_data(:,logical(frequencies>=5 & frequencies<=8)));
power_vars.alpha_power = nanmean(spectral_data(:,logical(frequencies>=8 & frequencies<=14)));
power_vars.beta_power = nanmean(spectral_data(:,logical(frequencies>=15 & frequencies<=25)));
power_vars.gamma_power = nanmean(spectral_data(:,logical(frequencies>=25 & frequencies<=40)));
end
