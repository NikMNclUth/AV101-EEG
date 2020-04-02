%% Set up loading paths
delete(gcp('nocreate'))
clear; clc;
setup_analysispath(5); % files from Ketamine folder and misc folder
setup_analysispath(6); % files from AV101 folder
eeglab;
clear
close all
% [output] = process_av101(EEGdata)

%% assr 40
directory = '/data/smathew/AV-101/processed_data/assr40/V2/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
directory = '/data/smathew/AV-101/processed_data/assr40/V3/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
directory = '/data/smathew/AV-101/processed_data/assr40/V4/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii

%% assr 30

directory = '/data/smathew/AV-101/processed_data/assr30/V2/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
directory = '/data/smathew/AV-101/processed_data/assr30/V3/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
directory = '/data/smathew/AV-101/processed_data/assr30/V4/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii

%% assr 20

directory = '/data/smathew/AV-101/processed_data/assr20/V2/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
directory = '/data/smathew/AV-101/processed_data/assr20/V3/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
directory = '/data/smathew/AV-101/processed_data/assr20/V4/';
files = dir(directory);
files = {files(3:end).name}';
dirlist = cell(length(files),1);
for ii = 1:length(files)
    dirlist{ii} = [directory,files{ii}];
end
clear ii
for ii = 1:length(dirlist)
    load(dirlist{ii});
    [output] = process_av101(EEGfinal);
    save(dirlist{ii},'output','-append');
end
clear EEGfinal output summary ii
