function [matclustTet] = load_matclust_data(dataDir,loadRat,filenum)
%LOAD_MATCLUST_DATA Loads matclust data structures from dataDir
%
% This function loads matclust data from a single rat,
% grabbing one tetrode's files located at dataDir.
%
% Inputs:
%   loadRat - string of a single rat name to load.
%
%   filenum - the index to be used for selecting a file out of the
%   directory. Can be between 1 and the number of matclust_param
%   files, otherwise an empty struct is returned.
%
% Outputs: 
%   C_matclust - Cell array to hold data for all rats. Rows are
%   tetrodes on which cells were clustered, columns are rats.
%   

filetype = "matclust_param";
matclustTet = struct();



fprintf("Loading matclust animal: %s \n",loadRat)

% Identify underscores in any rat names and only load matclust files based on the 
% base rat name.
tgtIdx = strfind(loadRat,"_");
if ~isempty(tgtIdx)
    % Gets the directories to all matclust_param files, which are the files
    % for tetrodes that have clustered cells on them.
    FileDir = dir(dataDir+"/"+loadRat(1:tgtIdx-1)+'.matclust'+"/"+filetype+"*");
else % If no '_' is present in rat name
    FileDir = dir(dataDir+"/"+loadRat+'.matclust'+"/"+filetype+"*");
end

% C_matclust = cell(size(FileDir,1),size(loadRats,2));

if isempty(FileDir)
    error("%s files do not exist for animal: %s \n",filetype,loadRat)
end

% for f = 1:size(FileDir,1) % loops through tetrodes
structFile = load(string(fullfile(FileDir(filenum).folder, FileDir(filenum).name))); % load data
clustattrib = structFile.clustattrib;

cutIdx = strfind(clustattrib.datafile,'waves_'); % find waves filename.
wavesFile = load(string(fullfile(FileDir(filenum).folder, clustattrib.datafile(cutIdx:end))));

% Assigns data to struct for output.
[matclustTet.clustattrib] = deal(clustattrib);
[matclustTet.waves] = deal(wavesFile.waves);

fprintf("       Loaded file: %s  \n",  FileDir(filenum).name)

% end
















end

