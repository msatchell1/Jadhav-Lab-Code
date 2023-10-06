function [] = create_nrninfo_MES(dataDir,loadRats)
% Creates cellinfo_MES.mat structs that contain information pertaining to
% each neuron across all epochs.
%
% This function saves a file to each rat's folder called nrninfo_MES.mat.
% It is a cell array with dimensions 1 x (num rats), within which is a cell
% array (tot num neurons) x 1. Each cell contains a struct with that neuron's
% information.

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; 

loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

filetypes = {'spikes01','cellinfo'};

C_alldata = load_data(dataDir,loadRats,filetypes);

spikes_idx = find(contains(filetypes,'spikes01')); 
if isempty(spikes_idx)
    error("spikes01 data must be loaded to run this analysis.")
end
cellinfo_idx = find(contains(filetypes,'cellinfo')); 
if isempty(cellinfo_idx)
    error("spikes01 data must be loaded to run this analysis.")
end

C_allspikes = C_alldata(spikes_idx,:);
C_allcellinfo = C_alldata(cellinfo_idx,:);



%%

% The spikes file is the best to use for determining the total number of
% neurons per rat because it is constant across epochs. Cellinfo varies in
% the number of neurons per epoch.
nrninfo_MES = cell(1,size(loadRats,2));

nrnCounts = zeros(17,size(loadRats,2));
for r = 1:size(loadRats,2)
    
    spikeData = cell(1,17);

    for e = 1:size(C_allspikes{1,r},2)
        nrnsAllTets = [C_allspikes{1,r}{1,e}{:}];
        nrnCounts(e,r) = length(nrnsAllTets);

        for n = 1:size(nrnsAllTets,2)




    end

    if any(~(nrnCounts(:,r) == nrnCounts(1,r))) % Compares all epochs with first epoch to 
        % find any inconsistencies in the number of neurons within a rat.
        error("The number of neurons does not match across epochs in spikes01 file for %s",loadRats{1,r})
    end

    
end

% Now the number of neurons in the first epoch can simply be used.


% I should plot spike width across epochs for cells to make sure that it
% doesn't vary too much from epoch to epoch... otherwise I might not be
% recording the same cell.








end