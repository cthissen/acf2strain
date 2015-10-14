function [folder] = loadDirectoryInfo()
% script to load directory information

% conventions:
% - all folder names end with file sep (using system filesep command).



fs = filesep;
folder.currentDirectory = pwd;
folder.baseDirectory = [folder.currentDirectory,fs];

% setup useful folder structure
folder.images   = [folder.baseDirectory,'Images',fs];
folder.data     = [folder.baseDirectory,'Data',fs];
folder.results  = [folder.baseDirectory,'Results',fs];



% add subfunction folders
addpath([folder.baseDirectory,'private',fs,'fitting',fs]);
addpath([folder.baseDirectory,'private',fs,'misc',fs]);
addpath([folder.baseDirectory,'private',fs,'imageprocessing',fs]);
addpath([folder.baseDirectory,'private',fs,'output',fs]);
addpath([folder.baseDirectory,'private',fs,'output',fs,'export_fig',fs]);
addpath([folder.baseDirectory,'private',fs,'fitting',fs,'splinefit',fs]);

% make directories if they don't already exist
if ~exist([folder.baseDirectory,'Data'],'dir')
    mkdir(folder.baseDirectory,'Data');
end
if ~exist([folder.baseDirectory,'Results'],'dir')
    mkdir(folder.baseDirectory,'Results');
end
if ~exist([folder.baseDirectory,'Images'],'dir')
    mkdir(folder.baseDirectory,'Images');
end

end
