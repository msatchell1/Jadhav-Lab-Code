% Script to identify trajectory-choice selective cells and analyze them
% Michael Satchell 11/02/23


%% Load data

clearvars;

dataDir = '/media/msatchell/10TBSpinDisk/js_SingleDayExpt'; % Location of data for all rats
% For this analysis I want to consider all rats, and all nrns on all
% tetrodes for each rat.

% All good rats: ZT2 ER1_NEW KL8 BG1 JS14 JS15 JS17 JS21 JS34 
loadRats = {'ZT2','ER1_NEW','KL8','BG1','JS14','JS15','JS17','JS21','JS34'};

% Common file types: 'cellinfo','sleep01','waking01','sws01','rem01','ripples01',
% 'spikes01','tetinfo','linfields01','rippletime01','pos01','nrninfo_MES','combstates_MES'
filetypes = {'nrninfo_MES','pos01','linfields01','combstates_MES'};

C_alldata = load_data(dataDir,loadRats,filetypes); % Load data to struct


C_linf = C_alldata(find(contains(filetypes,'linfields01')),:);
C_nrninfo = C_alldata(find(contains(filetypes,'nrninfo_MES')),:);
C_pos = C_alldata(find(contains(filetypes,'pos01')),:);
C_combstates = C_alldata(find(contains(filetypes,'combstates_MES')),:);

behEpochs = 2:2:17;
restEpochs = 1:2:17;
brainAreas = {'CA1','PFC'};


%% Identify choice-selective cells

% As a reminder: There are four columns here representing the four W-track trajectories.
% Col 1: out right, col 2: in right, col 3: out left, col 4: in left.

% For each cell, identify the mean FR along the center stem of the track
% when on outbound trajectories (0-80cm), and take the difference between
% this value for left and right trajectories, normalized by the sum of the FRs.
% This will give each cell an index from -1 to 1 on choice-selectivity for
% that epoch.


cellType = "CA1";
trajFR = cell(size(C_combstates{1,1})); % Firing rates on all trajectories for all epochs and cells
trajSC = cell(size(C_combstates{1,1})); % Spatial coverage
trajPk = cell(size(C_combstates{1,1})); % Peak firing rate
trajisPC = cell(size(C_combstates{1,1})); % Place cell (1) or not (0)
chcSltv = cell(size(C_combstates{1,1})); % Choice selectivity index.
for r = 1:size(C_nrninfo,2)
    
    for n = 1:size(C_nrninfo{1,r},1)

        nrn = C_nrninfo{1,r}{n,1};
        
        if strcmp(nrn.area,"CA1") && strcmp(nrn.type,"Pyr") %&& strcmp(nrn.LMRVtype,cellType)
    
            for e = 1:size(nrn.eTrajFR,2) % Loops through all epochs
                if ~isempty(nrn.eTrajFR{1,e}) % assume that if this nrn has FR data 
                    % for this epoch, then all other traj measures will
                    % exist.
                    trajFR{1,e} = [trajFR{1,e}; nrn.eTrajFR{1,e}];
                    trajSC{1,e} = [trajSC{1,e}; nrn.eTrajCoverage{1,e}];
                    trajPk{1,e} = [trajPk{1,e}; nrn.eTrajFRPeak{1,e}];
                    trajisPC{1,e} = [trajisPC{1,e}; nrn.eTrajisPC{1,e}];
                end
            end

        end


    end


end