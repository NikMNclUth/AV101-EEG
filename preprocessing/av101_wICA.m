function [W,w,s,A,icaact,wavcleanEEG,artifacts] = av101_wICA(EEG,badchans,init)




%% first round of ica
X = EEG.data;
X(unique(badchans),:) = [];
X = double(X);
if init == 1
    [W,w,s,A,IC,icaact] = initialiseICA(EEG,badchans);
elseif init == 2
    [~,w0]=fastica(X,'approach','symm','verbose', 'off','g','tanh');
    [w,s]=runica(X,'weights',w0,'extended',1,'verbose','off');
    W=w*s;
    A = inv(W);
    icaact=data2act(X,W);
    IC = reshape(icaact,size(icaact,1),[]);
elseif init == 3
    [w,s]=runica(X,'extended',1,'verbose','off');
    W=w*s;
    A = inv(W);
    icaact=data2act(X,W);
    IC = reshape(icaact,size(icaact,1),[]);
end

%% pad data
% padding data for proper wavelet transform...data must be divisible by
% 2^L, where L = level set for the stationary wavelet transform

mult=1;
L=5;
wavename='coif5';
modulus = mod(size(IC,2),2^L); %2^level (level for wavelet)
if modulus ~=0
    extra = zeros(1,(2^L)-modulus);
else
    extra = [];
end

%% wavelet thresholding
% loop through ICs and perform wavelet thresholding
disp('Performing wavelet thresholding');
for s = 1:size(IC,1)
    if ~isempty(extra)
        sig = [IC(s,:),extra]; % pad with zeros
    else
        sig = IC(s,:);
    end
    [thresh,sorh,~] = ddencmp('den','wv',sig); % get automatic threshold value
    thresh = thresh*mult; % multiply threshold by scalar
    swc = swt(sig,L,wavename); % use stationary wavelet transform (SWT) to wavelet transform the ICs
    Y = wthresh(swc,sorh,thresh); % threshold the wavelet to remove small values
    wIC(s,:) = iswt(Y,wavename); % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
    clear y sig thresh sorh swc
end

%% remove padding
% remove extra padding
if ~isempty(extra)
    wIC = wIC(:,1:end-numel(extra));
end
% % plot the ICs vs. wICs
% if plotting>0
%     disp('Plotting');
%     subplot(3,1,1);
%     multisignalplot(IC,Fs,'r');
%     title('ICs');
%     subplot(3,1,2);
%     multisignalplot(wIC,Fs,'r');
%     title('wICs')
%     subplot(3,1,3);
%     multisignalplot(IC-wIC,Fs,'r');
%     title('Difference (IC - wIC)');
% end

%% create artifact reduced data
%reconstruct artifact signal as channelsxsamples format from the wavelet coefficients
artifacts = A*wIC;
%reshape EEG signal from EEGlab format to channelsxsamples format
EEG2D=reshape(X, size(X,1), []);

%subtract out wavelet artifact signal from EEG signal
wavcleanEEG=EEG2D-artifacts;

%% run ica on wavelet reduced data
[w,s]=runica(wavcleanEEG,'extended',1,'verbose','off');
if init == 1
    W = w*s*pinv(L);
else
    W=w*s;
end
A = inv(W);
icaact=data2act(X,W,s);

end

%% subfunctions

function [W,w,s,A,IC,icaact,L] = initialiseICA(EEG,badchans)


X = EEG.data;
X = EEG.data;
X(badchans,:) = [];
X = double(X);
nnans=~any(isnan(X));


% PCA
fprintf('Going to principal components [PCA]... ')
[L,PC,E]=princomp_BU(X(:,nnans)');
E=E/sum(E);
% see how many components fastica would retain
e=fastica(X(:,nnans),'verbose','off','only','pca');
ncomp=size(e,2);
fprintf('\nRetaining %s/%s components ',num2str(ncomp),num2str(size(L,2)))
PC=PC(:,1:ncomp)';
L=L(:,1:ncomp);
E=E(1:ncomp);
% run initial fastica to get weights
% - symm is faster, tanh handles better quasi gaussian components
[~,w0]=fastica(PC,'approach','symm','verbose', 'off','g','tanh');
if w0
    fprintf('\nRunning initialized extended ICA [RUNICA]... ');
    [w,s]=runica(PC,'weights',w0,'extended',1,'verbose','off');
else
    fprintf('\nRunning NON initialized extended ICA [RUNICA]... ');
    [w,s]=runica(PC,'extended',1,'verbose','off');
end
% algebra to get component activations
W=w*s*pinv(L);
A = inv(W);
icaact=data2act(X,W);
IC = reshape(icaact,size(icaact,1),[]);

end