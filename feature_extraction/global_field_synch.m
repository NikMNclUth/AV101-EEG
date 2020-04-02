function [smooth_gfs] = global_field_synch(data,evoked,sr,order,winlength)
rmpath('/data/rcho/TOOLS/eeglab14_1_1b/plugins/AAR131130/');

% if epoch == 0
%     data = data(:,:);
%     FT=fft(hanning(length(data))'.*data,[],2);
%     % get one sided spectrum
%     f = sr*(0:(length(data)/2))/length(data);
%     FT = FT(:,1:length(f));
%     gfs = zeros(1,length(f));
%     for ii = 1:length(FT)
%         temp = zeros(size(data,1),2);
%         temp(:,1) = imag(FT(:,ii));
%         temp(:,2) = real(FT(:,ii));
%         [~,~,latent,~,~,~] = pca(temp);
%         gfs(ii) = (abs(latent(1)-latent(2)))/(latent(1)+latent(2));
%     end
%     [coeff,score,latent,tsquared,explained,mu] = pca(temp);
%     gfs = (abs(latent(1)-latent(2)))/(latent(1)+latent(2));
if evoked == 1
    %     FT=fft(hanning(size(data,2))'.*data,[],2);
    %     f= sr*(0:(size(data,2)/2))/size(data,2);
    %     FT = FT(:,1:length(f),:);
    %     gfs = zeros(size(data,3),length(f));
    %     for iii = 1:size(data,3)
    %         for ii = 1:length(FT)
    %             temp = zeros(64,2);
    %             temp(:,1) = imag(FT(:,ii,iii));
    %             temp(:,2) = real(FT(:,ii,iii));
    %             [~,~,latent,~,~,~] = pca(temp);
    %             gfs(iii,ii) = (abs(latent(1)-latent(2)))/(latent(1)+latent(2));
    %         end
    %     end
    % end
    
    FT =fft(hanning(size(data,2))'.*mean(data,3),[],2);
    gfs = zeros(1,length(FT));
    for ii = 1:length(FT)
        temp = zeros(size(data,1),2);
        temp(:,1) = imag(FT(:,ii));
        temp(:,2) = real(FT(:,ii));
        [~,~,latent,~,~,~] = pca(temp);
        gfs(:,ii) = (abs(latent(1)-latent(2)))/(latent(1)+latent(2));
    end
    f= sr*(0:(size(data,2)/2))/size(data,2);
    gfs = gfs(1:length(f));
else
    disp('trialwise')
    FT=fft(hanning(size(data,2))'.*data,[],2);
    f= sr*(0:(size(data,2)/2))/size(data,2);
    FT = FT(:,1:length(f),:);
    gfs = zeros(size(data,3),length(f));
    for iii = 1:size(data,3)
        for ii = 1:length(FT)
            temp = zeros(64,2);
            temp(:,1) = imag(FT(:,ii,iii));
            temp(:,2) = real(FT(:,ii,iii));
            [~,~,latent,~,~,~] = pca(temp);
            gfs(iii,ii) = (abs(latent(1)-latent(2)))/(latent(1)+latent(2));
        end
    end
end

% apply smoothing
% smooth_gfs = sgolayfilt(gfs, 10, 301);
smooth_gfs = sgolayfilt(mean(gfs,1), order, winlength);


end