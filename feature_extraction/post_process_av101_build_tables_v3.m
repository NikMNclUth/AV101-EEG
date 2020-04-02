setup_analysispath(6)
clear;close all;clc
directory = '/data/smathew/AV-101/processed_data/assr40/V2/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist2 = cell(length(files),1);
for ii = 1:length(files)
    dirlist2{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist2)
    dirlist2{ii,2}='V2';
end

directory = '/data/smathew/AV-101/processed_data/assr40/V3/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist3 = cell(length(files),1);
for ii = 1:length(files)
    dirlist3{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist3)
    dirlist3{ii,2}='V3';
end

directory = '/data/smathew/AV-101/processed_data/assr40/V4/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist4 = cell(length(files),1);
for ii = 1:length(files)
    dirlist4{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist4)
    dirlist4{ii,2}='V4';
end

files = [dirlist2;dirlist3;dirlist4];
for ii = 1:length(files)
    load(files{ii,1});
    ROI={'F1';'FZ';'F2';'FC3';'FC1';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4';'CP1';'CP2';'CPZ'};
    chans = {EEGfinal.chanlocs.labels};
    candidate=find(ismember(chans,ROI));
    %% DB
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'D');
    Yydb = nanmean(Y(:,:,candidate),3);
    [~,~,PLF]=timefq_speedyv2(EEGfinal.data(candidate,:,:),38:42,250,12,10,2,1);
    PLFp(ii,1) = nanmean(nanmean(nanmean(PLF(:,401:626,:),3),2),1);
    gammadb(ii,1) = nanmean(nanmean(Yydb(38:42,401:626,:),2),1);
    
    %% Z
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'zscore');
    Yyz = nanmean(Y(:,:,candidate),3);
    gammaz(ii,1) = nanmean(nanmean(Yyz(38:42,401:626,:),2),1);
    %% Perc
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'percentage');
    Yyp = nanmean(Y(:,:,candidate),3);
    gammap(ii,1) = nanmean(nanmean(Yyp(38:42,401:626,:),2),1);
    
end
statsout_roi_40 = table(files(:,1),files(:,2),gammadb,gammaz,gammap,PLFp);
statsout_roi_40.Properties.VariableNames = {'contents','visit','roi_assr_db','roi_assr_z','roi_assr_prc','roi_plf'};
%%
gammap = [];gammadb=[];gammaz=[];PLF=[];PLFp=[];files=[];filesv2=[];filesv3=[];filesv4=[];
directory = '/data/smathew/AV-101/processed_data/assr30/V2/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist2 = cell(length(files),1);
for ii = 1:length(files)
    dirlist2{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist2)
    dirlist2{ii,2}='V2';
end

directory = '/data/smathew/AV-101/processed_data/assr30/V3/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist3 = cell(length(files),1);
for ii = 1:length(files)
    dirlist3{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist3)
    dirlist3{ii,2}='V3';
end

directory = '/data/smathew/AV-101/processed_data/assr30/V4/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist4 = cell(length(files),1);
for ii = 1:length(files)
    dirlist4{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist4)
    dirlist4{ii,2}='V4';
end

files = [dirlist2;dirlist3;dirlist4];
for ii = 1:length(files)
    load(files{ii});
    ROI={'F1';'FZ';'F2';'FC3';'FC1';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4';'CP1';'CP2';'CPZ'};
    chans = {EEGfinal.chanlocs.labels};
    candidate=find(ismember(chans,ROI));
    %% DB
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'D');
    Yydb = nanmean(Y(:,:,candidate),3);
    [~,~,PLF]=timefq_speedyv2(EEGfinal.data(candidate,:,:),38:42,250,12,10,2,1);
    PLFp(ii,1) = nanmean(nanmean(nanmean(PLF(:,401:626,:),3),2),1);
    gammadb(ii,1) = nanmean(nanmean(Yydb(38:42,401:626,:),2),1);
    
    %% Z
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'zscore');
    Yyz = nanmean(Y(:,:,candidate),3);
    gammaz(ii,1) = nanmean(nanmean(Yyz(38:42,401:626,:),2),1);
    %% Perc
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'percentage');
    Yyp = nanmean(Y(:,:,candidate),3);
    gammap(ii,1) = nanmean(nanmean(Yyp(38:42,401:626,:),2),1);
    
end
statsout_roi_30 = table(files(:,1),files(:,2),gammadb,gammaz,gammap,PLFp);
statsout_roi_30.Properties.VariableNames = {'contents','visit','roi_assr_db','roi_assr_z','roi_assr_prc','roi_plf'};

%%
gammap = [];gammadb=[];gammaz=[];PLF=[];PLFp=[];files=[];filesv2=[];filesv3=[];filesv4=[];
directory = '/data/smathew/AV-101/processed_data/assr20/V2/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist2 = cell(length(files),1);
for ii = 1:length(files)
    dirlist2{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist2)
    dirlist2{ii,2}='V2';
end

directory = '/data/smathew/AV-101/processed_data/assr20/V3/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist3 = cell(length(files),1);
for ii = 1:length(files)
    dirlist3{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist3)
    dirlist3{ii,2}='V3';
end

directory = '/data/smathew/AV-101/processed_data/assr20/V4/';
files = dir(directory);
files = {files(3:end).name}';
files = natsortfiles(files);
dirlist4 = cell(length(files),1);
for ii = 1:length(files)
    dirlist4{ii} = [directory,files{ii}];
end
for ii = 1:length(dirlist4)
    dirlist4{ii,2}='V4';
end

files = [dirlist2;dirlist3;dirlist4];

for ii = 1:length(files)
    load(files{ii});
    ROI={'F1';'FZ';'F2';'FC3';'FC1';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4';'CP1';'CP2';'CPZ'};
    chans = {EEGfinal.chanlocs.labels};
    candidate=find(ismember(chans,ROI));
    %% DB
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'D');
    Yydb = nanmean(Y(:,:,candidate),3);
    [~,~,PLF]=timefq_speedyv2(EEGfinal.data(candidate,:,:),38:42,250,12,10,2,1);
    PLFp(ii,1) = nanmean(nanmean(nanmean(PLF(:,401:626,:),3),2),1);
    gammadb(ii,1) = nanmean(nanmean(Yydb(38:42,401:626,:),2),1);
    
    %% Z
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'zscore');
    Yyz = nanmean(Y(:,:,candidate),3);
    gammaz(ii,1) = nanmean(nanmean(Yyz(38:42,401:626,:),2),1);
    %% Perc
    Y = baseline(abs(wavtransform(nanmean(EEGfinal.data,3),1:60,250,10)).^2,1:375,2,'percentage');
    Yyp = nanmean(Y(:,:,candidate),3);
    gammap(ii,1) = nanmean(nanmean(Yyp(38:42,401:626,:),2),1);
    
end
statsout_roi_20 = table(files(:,1),files(:,2),gammadb,gammaz,gammap,PLFp);
statsout_roi_20.Properties.VariableNames = {'contents','visit','roi_assr_db','roi_assr_z','roi_assr_prc','roi_plf
    

savepath = '/data/smathew/AV-101/processed_data/';
savename = [savepath,'stats_out_40hz_assr_roi_fixed.xls'];
writetable(statsout_roi_40,savename,'Sheet',1);
writetable(statsout_roi_30,savename,'Sheet',2);
writetable(statsout_roi_20,savename,'Sheet',3);



