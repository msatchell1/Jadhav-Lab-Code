function [matclustTet] = load_matclust_data(dataDir,loadRat,tetnum)
%LOAD_MATCLUST_DATA Loads matclust data structures from dataDir
%
% This function loads matclust data from a single rat,
% grabbing one tetrode's files located at dataDir. Unfortunately trying to
% load all of a rat's files at once crashes matlab, so this function
% returns the clustered data of one tetrode at a time.
%
% Inputs:
%   loadRat - string of a single rat name to load.
%
%   filenum - the index to be used for selecting a file out of the
%   directory. Can be between 1 and the number of matclust_param
%   files, otherwise an empty struct is returned.
%
%   tetnum - tetrode number to load.
%
% Outputs: 
%   C_matclust - Cell array to hold data for all rats. Rows are
%   tetrodes on which cells were clustered, columns are rats. In this
%   version of the function only one tetrode (tetnum) is loaded.
%   

filetype = "matclust_param_nt"+string(tetnum)+".mat";
matclustTet = struct();



% fprintf("Loading matclust animal: %s \n",loadRat)

% Identify underscores in any rat names and only load matclust file of the
% tetrode tetnum.
tgtIdx = strfind(loadRat,"_");
if ~isempty(tgtIdx)
    % Gets the directories to all matclust_param files, which are the files
    % for tetrodes that have clustered cells on them.
    FileDir = dir(dataDir+"/"+loadRat(1:tgtIdx-1)+'.matclust'+"/"+filetype);
else % If no '_' is present in rat name
    FileDir = dir(dataDir+"/"+loadRat+'.matclust'+"/"+filetype);
end

% C_matclust = cell(size(FileDir,1),size(loadRats,2));

if isempty(FileDir)
    fprintf("       File does not exist: %s \n",filetype)
elseif size(FileDir,1) > 1 % Only one tetrode file should be loaded per rat
    error("Cannot load multiple tetrodes: %s \n",FileDir.name)
else
    try
        structFile = load(string(fullfile(FileDir.folder, FileDir.name))); % load data
        clustattrib = structFile.clustattrib;
    
        cutIdx = strfind(clustattrib.datafile,'waves_'); % find waves filename.
        wavesFile = load(string(fullfile(FileDir.folder, clustattrib.datafile(cutIdx:end))));
        
        % Assigns data to struct for output.
        [matclustTet.clustattrib] = deal(clustattrib);
        [matclustTet.waves] = deal(wavesFile.waves);
        
        fprintf("       Loaded file: %s  \n",  FileDir.name)

    catch
        fprintf("       Error: Could not load: %s \n",FileDir.name)
    end


end








end

