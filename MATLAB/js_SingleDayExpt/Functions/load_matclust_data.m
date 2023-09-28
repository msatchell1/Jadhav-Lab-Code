function [C_matclust] = load_matclust_data(dataDir,loadRats)
%LOAD_MATCLUST_DATA Loads matclust data structures from dataDir
%
% This function loads matclust data from rats specified in loadRats,
% grabbing the files under filetypes located at dataDir.
%
% Outputs: C_matclust - Cell array to hold data for all rats. Rows are
% tetrodes on which cells were clustered, columns are rats.
%   

filetype = "matclust_param";
C_matclust = {}; 

disp("Loading matclust data... ")
for r = 1:length(loadRats)
    fprintf("Loading animal: %s \n",loadRats{r})
   
   % Gets the directories to all matclust_param files, which are the files
   % for tetrodes that have clustered cells on them.
    FileDir = dir(dataDir+"/"+loadRats(r)+'.matclust'+"/"+filetype+"*");

    if isempty(FileDir)
        error("%s files do not exist for animal: %s \n",filetype,loadRats{r})
    end
    
    for f = 1:size(FileDir,1) % loops through tetrodes
        structFile = load(string(fullfile(FileDir(f).folder, FileDir(f).name))); % load data
        clustattrib = structFile.clustattrib;

        cutIdx = strfind(clustattrib.datafile,'waves_'); % find waves filename.
        wavesFile = load(string(fullfile(FileDir(f).folder, clustattrib.datafile(cutIdx:end))));
        
        C_matclust{f,r} = struct('clustattrib',clustattrib, 'waves',wavesFile.waves);
        fprintf("       Loaded file: %s  \n",  FileDir(f).name)

    end
end





fprintf("Finished loading matclust data.")














end

