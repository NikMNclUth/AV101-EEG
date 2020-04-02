function [output] = process_av101(EEGdata)


% %% GF ASSR
Y = baseline(abs(wavtransform(nanmean(EEGdata.data,3),1:60,250,10)).^2,1:375,2);
% [GFP,~]=gfp_new(Y,[],[],[],[]);
gfpwav = zeros(60,750);
for iter = 1:60
    [gfpwav(iter,:),~]=gfp_new(squeeze(Y(iter,:,:))',[],[],[],[]);
end
%% natural frequency methods
gpower = nanmean(Y,3);
gnat = nansum(gpower(15:60,401:626),2);% globalz
[gmaxsum,gmaxfrq] = max(gnat);
gmaxfrq = gmaxfrq+15;
%% synchronisation
[smooth_gfs_ev] = global_field_synch(EEGdata.data,1,250,10,301);
[smooth_gfs] = global_field_synch(EEGdata.data,1,250,10,301);

%% PLF
[~,~,PLF]=timefq_speedyv2(EEGdata.data,1:60,250,12,10,2,1);
gfpplf = nanmean(PLF,3);
%% output
output.gmaxsum =gmaxsum;
output.gmaxfrq = gmaxfrq;
output.gfpwav = gfpwav;
output.smooth_gfs_ev = smooth_gfs_ev;
output.smooth_gfs = smooth_gfs;
output.gfpplf = gfpplf;

end

function [GFP,LFP]=gfp_new(X,dim,bl,trialwise,chans)
% [GFP,LFP]=gpf(X,dim,bl,trialwise,indices);
%       Compute global (GFP) and local (LFP) field power for signal X.
%       X -- input data, assumes channel X time X trial
%	dim -- trial dimension
%	bl -- baseline time points if that needs to be done
%   trialwise -- keep trial dimension (true|false)
%	indices -- channels for LFP calculation
if nargin==0
    help gfp
else
    % initialize variables
    if ~exist('chans','var')||isempty(chans)
        chans=[];
        LFP=[];
    end
    if ~exist('bl','var')||isempty(bl)
        bl=-1;
    end
    if ~exist('dim','var')||isempty(dim)
        dim=3;
    end
    if ~exist('trialwise','var')||isempty(trialwise)
        trialwise=false;
    end
    %handle datasets larger than 3 dims
    origsize=size(X);
    X=X(:,:,:);
    % baseline data or remove mean
    if bl==-1
        %         X=bsxfun(@minus,X,nanmean(X,2));
    else
        X=bsxfun(@minus,X,nanmean(X(:,bl,:),2));
    end
    % calculate channel means mean and GFP
    if ~trialwise
        M=nanmean(X,dim);
    else
        M=X;
    end
    GFP=nanstd(M,1,1);
    % calculate LFP if channel indices are given
    if ~isempty(chans)
        LFP=nanstd(M(chans,:,:),1,1);
    end
    % if indices are passed and only one output is requested, return LFP
    % instead of GFP
    if nargout==1 && ~isempty(chans)
        GFP=LFP;
    end
    %for convenience, make sure the outputs have as few dimensions as
    %possible.
    if trialwise
        GFP=reshape(GFP,[1 origsize(2:end)]);
        GFP=squeeze(GFP);
        if nargout>1
            LFP=reshape(LFP,[1 origsize(2:end)]);
            LFP=squeeze(LFP);
        end
    else
        GFP=squeeze(GFP);
        LFP=squeeze(LFP);
    end
end
end