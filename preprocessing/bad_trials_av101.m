function [badtrials] = bad_trials_av101(EEGepoch,bad_chans)

data = double(EEGepoch.data);
data(bad_chans,:) = NaN;
[badindsA] = bt_amp(data);
[badindsV] = bt_var(data);
[badindsE] = bt_emd(data);
badtrials = sort(unique([badindsA,badindsV,badindsE]));
end

function [badindsA] = bt_amp(data)
for t = 1:size(data,1)
    for u = 1:size(data,3)
        ampdiffs(t,u) = max(data(t,:,u)) - min(data(t,:,u));
    end
end
ad = nanmean(ampdiffs,1);
adz = (ad-nanmean(ad))./nanstd(ad);
badindsA = find(adz>=3);
end

function [badindsV] = bt_var(data)
ev = nanmean(squeeze(nanvar(data,0,2)));
evz = (ev-nanmean(ev))./nanstd(ev);
badindsV = find(evz>=3);
end

function [badindsE] = bt_emd(data)
emd = zeros(1,size(data,3));
means = squeeze(nanmean(data(:,:),2));
for u = 1:size(data,3)
	emd(u) = nanmean(abs(squeeze(nanmean(data(:,:,u),2)) - means));
end
emdz = (emd-nanmean(emd))./nanstd(emd);
badindsE = find(emdz>=3);
end